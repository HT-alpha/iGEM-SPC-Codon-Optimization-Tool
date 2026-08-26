//! Composition: resolve terms, evaluate them, assemble a response.
//!
//! Ticket P1-22.
//!
//! ## The rule this module exists to enforce
//!
//! ADR-0002 fixes the landscape as a black-box oracle with an *optional* exact
//! delta. Optional is the load-bearing word. A composition layer has two easy
//! ways to get it wrong:
//!
//! 1. Require the delta, and no L0 predictor can ever be plugged in.
//! 2. Notice one term cannot answer a delta and rescore *everything*, which is
//!    correct and throws away the measured 640x-at-2 kb, 4800x-at-15 kb
//!    advantage that justified the whole capability model.
//!
//! [`Composer::evaluate`] does neither. Each term takes its own path: L1 terms
//! get their exact delta, L0 terms get a rescore against a mutated sequence
//! that is materialized once and shared, and `diagnostics.full_rescores`
//! reports how many terms took the expensive route. That counter is what makes
//! the rule testable — a composition that degraded globally would still return
//! correct numbers, and only the counter distinguishes it.
//!
//! No path here calls for a gradient, because `TermProvider` has none.
//!
//! ## `UNKNOWN` is propagated, never imputed
//!
//! A term that could not be computed within budget is absent from `terms`, is
//! listed in `diagnostics.terms_deferred`, and appears in the objective vector
//! as `None`. The NSGA-III constraint-domination rule defers such a candidate
//! rather than comparing it. Substituting a mean or a worst case would bias the
//! front in a way that leaves no trace in the output.

use std::collections::HashSet;
use std::time::Instant;

use cf_landscape_api::{
    Capability, Deltas, Diagnostics, Edit, EvalContext, LandscapeQuery, LandscapeResponse,
    ObjectiveAxis, ObjectiveVector, Registration, Status, Target, TermRegistry, TermRole,
    TermValue,
};
use cf_types::{code, CfError, CodonSequence, Digest, Result, TermId};

use crate::cache::EvaluationCache;

/// Terms a host has the reference data to support.
///
/// Supplied by the orchestrator from the host profile's `supported_terms`.
/// `None` means "not declared", which is treated as "no gate" rather than as
/// "nothing supported" — a host profile predating the field must not silently
/// disable every objective.
#[derive(Clone, Debug, Default)]
pub struct HostCapability(Option<HashSet<String>>);

impl HostCapability {
    pub fn unrestricted() -> Self {
        Self(None)
    }

    pub fn supporting(terms: impl IntoIterator<Item = String>) -> Self {
        Self(Some(terms.into_iter().collect()))
    }

    fn permits(&self, term: &TermId) -> bool {
        match &self.0 {
            None => true,
            Some(set) => set.contains(term.as_str()),
        }
    }
}

/// The composition layer.
pub struct Composer {
    registry: TermRegistry,
    cache: EvaluationCache,
    host: HostCapability,
    /// Base sequences seen in this session, by document digest.
    ///
    /// A delta query names its base by digest and carries no sequence. If the
    /// base is gone the query cannot be answered, and the contract says so:
    /// `CF-LAND-010`, resend as a full query.
    bases: std::collections::HashMap<String, CodonSequence>,
    base_capacity: usize,
    base_order: std::collections::VecDeque<String>,
}

impl std::fmt::Debug for Composer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Composer")
            .field("terms", &self.registry.term_ids())
            .field("namespace", &self.cache.namespace())
            .field("bases", &self.bases.len())
            .finish()
    }
}

/// Per-alternative term deltas, grouped by codon index.
///
/// Built up across terms — each L1 provider contributes its own column — so the
/// inner map has to be keyed by codon before it is keyed by term, which is the
/// opposite of how the response document reads.
type AlternativesByPosition =
    std::collections::BTreeMap<usize, Vec<(cf_types::Codon, Vec<(TermId, f64)>)>>;

/// How one term's value was obtained. Drives the diagnostics counters.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Route {
    Cached,
    Scored,
    ExactDelta,
    /// The term could not answer a delta query, so it was rescored in full.
    /// Recorded per term, which is the whole point.
    FullRescore,
}

struct Evaluated {
    value: f64,
    route: Route,
}

impl Composer {
    /// Default number of base sequences retained for delta queries.
    pub const DEFAULT_BASE_CAPACITY: usize = 256;

    /// Build a composer over a registry, namespaced by the landscape lockfile
    /// digest.
    pub fn new(registry: TermRegistry, lockfile_digest: impl Into<String>) -> Self {
        let namespace = lockfile_digest.into();
        Self {
            registry,
            cache: EvaluationCache::new(namespace),
            host: HostCapability::unrestricted(),
            bases: std::collections::HashMap::new(),
            base_capacity: Self::DEFAULT_BASE_CAPACITY,
            base_order: std::collections::VecDeque::new(),
        }
    }

    pub fn with_host_capability(mut self, host: HostCapability) -> Self {
        self.host = host;
        self
    }

    pub fn with_cache(mut self, cache: EvaluationCache) -> Self {
        self.cache = cache;
        self
    }

    pub fn registry(&self) -> &TermRegistry {
        &self.registry
    }

    pub fn cache(&self) -> &EvaluationCache {
        &self.cache
    }

    pub fn cache_mut(&mut self) -> &mut EvaluationCache {
        &mut self.cache
    }

    /// Remember a base sequence so later delta queries can name it.
    pub fn remember_base(&mut self, sequence: &CodonSequence) -> Result<Digest> {
        let digest = sequence.digest()?;
        let key = digest.as_str().to_owned();
        if !self.bases.contains_key(&key) {
            while self.bases.len() >= self.base_capacity {
                match self.base_order.pop_front() {
                    Some(oldest) => {
                        self.bases.remove(&oldest);
                    }
                    None => break,
                }
            }
            self.base_order.push_back(key.clone());
            self.bases.insert(key, sequence.clone());
        }
        Ok(digest)
    }

    pub fn forget_bases(&mut self) {
        self.bases.clear();
        self.base_order.clear();
    }

    /// Evaluate a query.
    pub fn evaluate(&mut self, query: &LandscapeQuery) -> Result<LandscapeResponse> {
        query.validate()?;
        self.check_host_support(query)?;

        let started = Instant::now();
        let hits_before = self.cache.hits();
        let misses_before = self.cache.misses();

        let (base, base_digest, edits) = self.resolve_target(query)?;
        let seed = self.effective_seed(query)?;

        // Materialized lazily and at most once, however many L0 terms need it.
        let mutated = if edits.is_empty() {
            None
        } else {
            let mut seq = base.clone();
            for edit in &edits {
                seq.substitute(edit.codon_index, edit.codon)?;
            }
            Some(seq)
        };

        let ctx_base = EvalContext::new(&base, &query.host_context, query.budget, seed);
        let ctx_target = match &mutated {
            Some(seq) => EvalContext::new(seq, &query.host_context, query.budget, seed),
            None => ctx_base,
        };

        let mut terms: Vec<(TermId, TermValue)> = Vec::new();
        let mut deferred: Vec<TermId> = Vec::new();
        let mut notes: Vec<String> = Vec::new();
        let mut evaluations = 0u64;
        let mut full_rescores = 0u64;

        // Resolved in the caller's order; the registry never imposes its own.
        let resolved: Vec<&Registration> = self.registry.resolve(&query.requested_terms)?;

        for (term, registration) in query.requested_terms.iter().zip(&resolved) {
            if over_deadline(query, started) {
                deferred.push(term.clone());
                notes.push(format!("{term} deferred: latency budget exhausted"));
                continue;
            }

            let key = self.cache.key(&base_digest, &edits, term);
            let outcome = match (query.cache_policy.read, self.cache.get(&key)) {
                (true, Some(value)) => Ok(Evaluated {
                    value,
                    route: Route::Cached,
                }),
                _ => Self::compute(registration, &ctx_base, &ctx_target, &edits),
            };

            match outcome {
                Ok(Evaluated { value, route }) => {
                    if route != Route::Cached {
                        evaluations += 1;
                        if query.cache_policy.write {
                            self.cache.put(key, value);
                        }
                    }
                    if route == Route::FullRescore {
                        full_rescores += 1;
                        notes.push(format!(
                            "{term} is {} for this input, so its delta was obtained by rescoring",
                            registration.provider.capability()
                        ));
                    }
                    let provider = registration.provider.as_ref();
                    let mut value = TermValue::new(
                        value,
                        provider.unit(),
                        provider.direction(),
                        provider.capability(),
                        provider.provenance(),
                    );
                    if query.outputs.uncertainty {
                        let scalar = value.value;
                        value = value.with_uncertainty(provider.uncertainty(&ctx_target, scalar));
                    }
                    terms.push((term.clone(), value));
                }
                Err(err) if err.code == code::BUDGET_EXCEEDED => {
                    deferred.push(term.clone());
                    notes.push(format!("{term} deferred: {}", err.detail));
                }
                Err(err) => return Err(err),
            }
        }

        let objective_vector = self.build_objective_vector(query, &terms)?;
        let deltas = self.build_deltas(query, &ctx_base, &base_digest)?;

        let status = if deferred.is_empty() {
            Status::Ok
        } else {
            Status::Partial
        };

        Ok(LandscapeResponse {
            query_id: query.query_id.clone(),
            status,
            terms,
            objective_vector,
            deltas,
            diagnostics: Diagnostics {
                cache_hits: self.cache.hits() - hits_before,
                cache_misses: self.cache.misses() - misses_before,
                evaluations,
                full_rescores,
                terms_deferred: deferred,
                engine_notes: notes,
            },
            problem: None,
        })
    }

    /// Evaluate one term, choosing its own path.
    ///
    /// The whole per-term-degradation rule is these six lines: an L1 term adds
    /// its exact delta to the base value, an L0 term is rescored against the
    /// mutated sequence, and neither decision is visible to the other terms.
    fn compute(
        registration: &Registration,
        ctx_base: &EvalContext<'_>,
        ctx_target: &EvalContext<'_>,
        edits: &[Edit],
    ) -> Result<Evaluated> {
        let provider = registration.provider.as_ref();
        if edits.is_empty() {
            return Ok(Evaluated {
                value: provider.score(ctx_base)?,
                route: Route::Scored,
            });
        }
        if provider.capability().offers_exact_delta() {
            match provider.delta(ctx_base, edits) {
                Ok(delta) => {
                    return Ok(Evaluated {
                        value: provider.score(ctx_base)? + delta,
                        route: Route::ExactDelta,
                    })
                }
                // A provider may be L1 in general and L0 for this input — the
                // violation counter is, when the constraint set contains a
                // k-mer repeat rule. Falling through is correct; the fallback
                // is still confined to this term.
                Err(err) if err.code == code::LAND_CAPABILITY => {}
                Err(err) => return Err(err),
            }
        }
        Ok(Evaluated {
            value: provider.score(ctx_target)?,
            route: Route::FullRescore,
        })
    }

    /// The objective vector, built from Pareto axes and nothing else.
    ///
    /// `TST-LAND-005`. A diagnostic term appears in `terms` and never here, so
    /// enabling it cannot move a ranking. The gate is
    /// [`ObjectiveAxis::new`], which refuses any other role, so a future edit
    /// that pushed a term id straight into the order would not compile.
    fn build_objective_vector(
        &self,
        query: &LandscapeQuery,
        terms: &[(TermId, TermValue)],
    ) -> Result<Option<ObjectiveVector>> {
        let mut axes = Vec::new();
        let mut sense = Vec::new();
        let mut values = Vec::new();

        for term in &query.requested_terms {
            let Some(role) = self.registry.role_of(term) else {
                continue;
            };
            if role != TermRole::ParetoAxis {
                continue;
            }
            let registration = self.registry.get(term).expect("role_of found it");
            axes.push(ObjectiveAxis::new(term.clone(), role)?);
            sense.push(registration.provider.direction());
            values.push(terms.iter().find(|(t, _)| t == term).map(|(_, v)| v.value));
        }

        if axes.is_empty() {
            return Ok(None);
        }
        ObjectiveVector::new(axes, sense, values)
            .map(Some)
            .map_err(|why| CfError::new(code::LAND_DOCUMENT_INVALID, why))
    }

    fn build_deltas(
        &self,
        query: &LandscapeQuery,
        ctx_base: &EvalContext<'_>,
        base_digest: &Digest,
    ) -> Result<Option<Deltas>> {
        let Target::DeltaTable { positions, .. } = &query.target else {
            return Ok(None);
        };
        if !query.outputs.deltas {
            return Ok(None);
        }

        let resolved = self.registry.resolve(&query.requested_terms)?;
        let mut per_position: AlternativesByPosition = Default::default();
        let mut exact = true;

        for (term, registration) in query.requested_terms.iter().zip(&resolved) {
            let provider = registration.provider.as_ref();
            if !provider.capability().offers_exact_delta() {
                // Absent from every entry, per the response schema: "missing
                // terms in an entry mean the provider is L0 for that term".
                // Not zero, which would read as "this move changes nothing".
                exact = false;
                continue;
            }
            let table = match provider.delta_table(ctx_base, positions) {
                Ok(table) => table,
                Err(err) if err.code == code::LAND_CAPABILITY => {
                    exact = false;
                    continue;
                }
                Err(err) => return Err(err),
            };
            for entry in table {
                let slot = per_position.entry(entry.codon_index).or_default();
                for (codon, delta) in entry.alternatives {
                    match slot.iter_mut().find(|(c, _)| *c == codon) {
                        Some((_, deltas)) => deltas.push((term.clone(), delta)),
                        None => slot.push((codon, vec![(term.clone(), delta)])),
                    }
                }
            }
        }

        Ok(Some(Deltas {
            base_digest: base_digest.clone(),
            exact,
            entries: per_position
                .into_iter()
                .map(|(codon_index, alternatives)| cf_landscape_api::DeltaEntry {
                    codon_index,
                    alternatives: alternatives
                        .into_iter()
                        .map(|(codon, delta)| cf_landscape_api::AlternativeDelta {
                            codon: Some(codon),
                            delta,
                        })
                        .collect(),
                })
                .collect(),
        }))
    }

    /// Refuse a term the host has no data for, at the composition boundary.
    ///
    /// `CF-HOST-003`. ADR-0004 makes this a real case rather than a defensive
    /// one: a host with only a genome assembly and no tRNA data cannot support
    /// `burden.growth_cost`, and the ingest layer records that in the profile's
    /// `supported_terms`. Without a gate here, the record would be advisory and
    /// the job would run with the objective silently missing.
    fn check_host_support(&self, query: &LandscapeQuery) -> Result<()> {
        for term in &query.requested_terms {
            if !self.host.permits(term) {
                return Err(CfError::new(
                    code::HOST_TERM_UNSUPPORTED,
                    format!(
                        "host {} has no reference data supporting {term}; \
                         it is absent from the profile's supported_terms",
                        query.host_context.host_id
                    ),
                ));
            }
        }
        Ok(())
    }

    fn effective_seed(&self, query: &LandscapeQuery) -> Result<u64> {
        let stochastic = query.requested_terms.iter().any(|term| {
            self.registry
                .get(term)
                .is_some_and(|r| r.provider.uncertainty_method().is_sampled())
        });
        query.effective_seed(stochastic && query.outputs.uncertainty)
    }

    /// The sequence to score, its digest, and the edits relative to it.
    fn resolve_target(
        &mut self,
        query: &LandscapeQuery,
    ) -> Result<(CodonSequence, Digest, Vec<Edit>)> {
        match &query.target {
            Target::Full { sequence, .. } => {
                let digest = self.remember_base(sequence)?;
                Ok(((**sequence).clone(), digest, Vec::new()))
            }
            Target::Delta {
                base_digest, edits, ..
            } => {
                let base = self.base(base_digest)?;
                Ok((base, base_digest.clone(), edits.clone()))
            }
            Target::DeltaTable { base_digest, .. } => {
                let base = self.base(base_digest)?;
                Ok((base, base_digest.clone(), Vec::new()))
            }
        }
    }

    fn base(&self, digest: &Digest) -> Result<CodonSequence> {
        self.bases.get(digest.as_str()).cloned().ok_or_else(|| {
            CfError::new(
                code::LAND_BASE_EVICTED,
                format!("base {digest} is not in this session's cache; resend as a full query"),
            )
        })
    }
}

fn over_deadline(query: &LandscapeQuery, started: Instant) -> bool {
    match query.budget.max_latency_ms {
        None => false,
        Some(ms) => started.elapsed().as_millis() > u128::from(ms),
    }
}

/// Capability a term will actually exercise for a given query shape.
///
/// Exposed so a caller can plan a tier without evaluating: ADR-0009's
/// compatibility filter refuses a tier-2 job whose only provider for a term is
/// L0 rather than discovering it mid-run.
pub fn planned_capability(registry: &TermRegistry, terms: &[TermId]) -> Vec<(TermId, Capability)> {
    terms
        .iter()
        .filter_map(|term| {
            registry
                .get(term)
                .map(|r| (term.clone(), r.provider.capability()))
        })
        .collect()
}
