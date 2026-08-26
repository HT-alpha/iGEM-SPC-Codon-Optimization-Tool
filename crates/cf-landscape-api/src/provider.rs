//! The provider seam and the term registry.
//!
//! Ticket P1-20.
//!
//! This is the trait every predictor implements and the only trait the
//! composition layer knows about. Three things about its shape are load-bearing
//! rather than incidental.
//!
//! **There is no gradient method.** Not an unimplemented one, not one behind a
//! feature flag — none. ADR-0002 says the engine never requires L2, and the way
//! to make that true permanently is to leave the composition layer no way to
//! ask. A provider that wants to expose gradients to some other consumer may;
//! it cannot do so through here.
//!
//! **Delta methods default to refusing.** [`Capability::L0`] is a legal,
//! supported answer, so the defaults return `CF-LAND-002` rather than
//! panicking or returning zero. Composition reads that refusal as "rescore this
//! term in full" and carries on — for that term only.
//!
//! **A capability claim is checkable by calling.** [`verify_capability`]
//! probes the surface a provider says it has and checks that the delta it
//! returns is the change actually observed by rescoring. An L1 claim is a claim
//! about *correctness*, not about the presence of a method, so verifying the
//! method exists would verify the wrong thing.

use std::collections::HashMap;

use cf_types::{code, CfError, Codon, CodonSequence, Direction, Result, TermId, Unit};

use crate::capability::Capability;
use crate::query::{Budget, Edit, HostContext};
use crate::response::{Provenance, Uncertainty, UncertaintyMethod};
use crate::role::TermRole;

/// Everything a provider needs to score one sequence.
#[derive(Clone, Copy, Debug)]
pub struct EvalContext<'a> {
    pub sequence: &'a CodonSequence,
    pub host: &'a HostContext,
    pub budget: Budget,
    /// The RNG stream seed. Resolved by
    /// [`crate::query::LandscapeQuery::effective_seed`], which refuses to
    /// invent one for a sampled provider.
    pub seed: u64,
}

impl<'a> EvalContext<'a> {
    pub fn new(
        sequence: &'a CodonSequence,
        host: &'a HostContext,
        budget: Budget,
        seed: u64,
    ) -> Self {
        Self {
            sequence,
            host,
            budget,
            seed,
        }
    }

    /// The same context over a different sequence.
    pub fn with_sequence<'b>(&self, sequence: &'b CodonSequence) -> EvalContext<'b>
    where
        'a: 'b,
    {
        EvalContext {
            sequence,
            host: self.host,
            budget: self.budget,
            seed: self.seed,
        }
    }
}

/// Exact changes for every synonymous alternative at one position.
#[derive(Clone, PartialEq, Debug)]
pub struct PositionDeltas {
    pub codon_index: usize,
    /// `(alternative codon, change in this term)`. Excludes the codon already
    /// present, whose delta is zero by definition.
    pub alternatives: Vec<(Codon, f64)>,
}

/// One term's predictor.
pub trait TermProvider: Send + Sync {
    fn term_id(&self) -> &TermId;
    fn unit(&self) -> Unit;
    fn direction(&self) -> Direction;
    fn capability(&self) -> Capability;
    fn provenance(&self) -> Provenance;

    /// The strongest role this term is *permitted* to hold, whatever a job
    /// requests.
    ///
    /// `AMD-002` is the reason this exists. Folding kinetics is diagnostic-only
    /// until `EXP-06`, and a policy that lives in prose is one careless job
    /// configuration away from being violated. A provider that caps itself at
    /// [`TermRole::Diagnostic`] cannot be registered as an axis or as an
    /// epsilon constraint, so the cap holds against callers that have not been
    /// written yet.
    fn max_role(&self) -> TermRole {
        TermRole::ParetoAxis
    }

    fn uncertainty_method(&self) -> UncertaintyMethod {
        UncertaintyMethod::None
    }

    /// L0. Mandatory.
    fn score(&self, ctx: &EvalContext<'_>) -> Result<f64>;

    /// L1. The exact change in this term when `edits` are applied to
    /// `ctx.sequence`.
    fn delta(&self, ctx: &EvalContext<'_>, edits: &[Edit]) -> Result<f64> {
        let _ = (ctx, edits);
        Err(CfError::new(
            code::LAND_CAPABILITY,
            format!("{} is L0: no closed-form delta", self.term_id()),
        ))
    }

    /// L1. Every alternative at each listed position.
    fn delta_table(
        &self,
        ctx: &EvalContext<'_>,
        positions: &[usize],
    ) -> Result<Vec<PositionDeltas>> {
        let _ = (ctx, positions);
        Err(CfError::new(
            code::LAND_CAPABILITY,
            format!("{} is L0: no delta table", self.term_id()),
        ))
    }

    /// Predictive uncertainty for a computed value.
    fn uncertainty(&self, ctx: &EvalContext<'_>, value: f64) -> Uncertainty {
        let _ = (ctx, value);
        Uncertainty::deterministic()
    }
}

/// Tolerance for "the delta equals the finite change".
///
/// Relative, with an absolute floor for terms whose values pass through zero.
/// Scaled to leave room for the reassociation an incremental update performs —
/// it sums a different set of floats than a full rescore does, so bit-equality
/// is the wrong bar — while being far tighter than any real disagreement.
const DELTA_TOLERANCE_REL: f64 = 1e-9;
const DELTA_TOLERANCE_ABS: f64 = 1e-12;

fn close_enough(a: f64, b: f64) -> bool {
    let diff = (a - b).abs();
    diff <= DELTA_TOLERANCE_ABS + DELTA_TOLERANCE_REL * a.abs().max(b.abs())
}

/// Check a declared capability against the surface the provider actually has.
///
/// `TST-LAPI-002`. Probes with `ctx.sequence`, which must have at least one
/// editable position with a synonym; a provider cannot be verified against a
/// sequence that admits no moves, and saying so is better than passing
/// vacuously.
pub fn verify_capability(provider: &dyn TermProvider, ctx: &EvalContext<'_>) -> Result<()> {
    let base = provider.score(ctx)?;

    if !provider.capability().offers_exact_delta() {
        // An L0 provider must actually refuse. One that quietly answers delta
        // queries is understating itself, which is harmless for correctness and
        // expensive for the budget, so it is still reported.
        if provider.delta(ctx, &[]).is_ok() {
            return Err(CfError::new(
                code::LAND_CAPABILITY,
                format!(
                    "{} declares L0 but answers delta queries; declare L1",
                    provider.term_id()
                ),
            ));
        }
        return Ok(());
    }

    let Some((index, alternative)) = first_move(ctx.sequence) else {
        return Err(CfError::new(
            code::LAND_CAPABILITY,
            format!(
                "cannot verify {}: the probe sequence has no editable position with a synonym",
                provider.term_id()
            ),
        ));
    };

    let mutated = ctx.sequence.substituted(index, alternative)?;
    let truth = provider.score(&ctx.with_sequence(&mutated))? - base;

    let claimed = provider.delta(ctx, &[Edit::new(index, alternative)])?;
    if !close_enough(claimed, truth) {
        return Err(CfError::new(
            code::LAND_DELTA_INEXACT,
            format!(
                "{} declares {} but its delta at codon {index} is {claimed}, \
                 while rescoring gives {truth}",
                provider.term_id(),
                provider.capability()
            ),
        ));
    }

    let table = provider.delta_table(ctx, &[index])?;
    let entry = table
        .iter()
        .find(|e| e.codon_index == index)
        .ok_or_else(|| {
            CfError::new(
                code::LAND_DELTA_INEXACT,
                format!(
                    "{} returned a delta table with no entry for the requested position {index}",
                    provider.term_id()
                ),
            )
        })?;
    let tabled = entry
        .alternatives
        .iter()
        .find(|(c, _)| *c == alternative)
        .map(|(_, d)| *d)
        .ok_or_else(|| {
            CfError::new(
                code::LAND_DELTA_INEXACT,
                format!(
                    "{} omitted synonym {} from the delta table at codon {index}",
                    provider.term_id(),
                    alternative.to_string_upper()
                ),
            )
        })?;
    if !close_enough(tabled, truth) {
        return Err(CfError::new(
            code::LAND_DELTA_INEXACT,
            format!(
                "{} disagrees with itself at codon {index}: table says {tabled}, \
                 rescoring says {truth}",
                provider.term_id()
            ),
        ));
    }

    Ok(())
}

/// The first editable position that has a synonym, and one such synonym.
fn first_move(sequence: &CodonSequence) -> Option<(usize, Codon)> {
    sequence.editable_positions().find_map(|i| {
        let current = sequence.codons()[i];
        sequence
            .synonyms_at(i)
            .iter()
            .copied()
            .find(|c| *c != current)
            .map(|c| (i, c))
    })
}

/// A registered term: its provider and the role it holds in this job.
pub struct Registration {
    pub provider: Box<dyn TermProvider>,
    pub role: TermRole,
}

impl std::fmt::Debug for Registration {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Registration")
            .field("term", &self.provider.term_id().as_str())
            .field("role", &self.role)
            .field("capability", &self.provider.capability())
            .finish()
    }
}

/// Resolves term ids to providers.
///
/// Holds a map for lookup but never iterates it to produce output. Every
/// method that yields more than one term takes the caller's order and returns
/// that order, because map iteration order is an implementation detail and an
/// objective vector ordered by one is not reproducible.
#[derive(Default)]
pub struct TermRegistry {
    by_term: HashMap<String, Registration>,
}

impl std::fmt::Debug for TermRegistry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TermRegistry")
            .field("terms", &self.term_ids())
            .finish()
    }
}

impl TermRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a provider in a role.
    ///
    /// Refuses a role stronger than the provider permits (`CF-LAND-004`) and a
    /// duplicate term id (`CF-LAND-005`).
    pub fn register(&mut self, provider: Box<dyn TermProvider>, role: TermRole) -> Result<()> {
        let term = provider.term_id().clone();
        if self.by_term.contains_key(term.as_str()) {
            return Err(CfError::new(
                code::LAND_TERM_DUPLICATE,
                format!("term {term} is already registered"),
            ));
        }
        let cap = provider.max_role();
        if influence_rank(role) > influence_rank(cap) {
            return Err(CfError::new(
                code::LAND_ROLE_NOT_AXIS,
                format!(
                    "term {term} may hold at most the role {cap}, but was registered as {role}"
                ),
            ));
        }
        self.by_term
            .insert(term.as_str().to_owned(), Registration { provider, role });
        Ok(())
    }

    /// Register in the provider's own strongest permitted role.
    pub fn register_default(&mut self, provider: Box<dyn TermProvider>) -> Result<()> {
        let role = provider.max_role();
        self.register(provider, role)
    }

    pub fn len(&self) -> usize {
        self.by_term.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_term.is_empty()
    }

    pub fn contains(&self, term: &TermId) -> bool {
        self.by_term.contains_key(term.as_str())
    }

    pub fn get(&self, term: &TermId) -> Option<&Registration> {
        self.by_term.get(term.as_str())
    }

    pub fn role_of(&self, term: &TermId) -> Option<TermRole> {
        self.by_term.get(term.as_str()).map(|r| r.role)
    }

    /// Resolve a term list **in the caller's order**.
    ///
    /// Fails with `CF-LAND-001` naming the first missing term. It does not
    /// return the terms it could resolve: a job that optimizes four objectives
    /// when the user asked for five, without saying so, is worse than a job
    /// that fails.
    pub fn resolve<'a>(&'a self, terms: &[TermId]) -> Result<Vec<&'a Registration>> {
        terms
            .iter()
            .map(|term| {
                self.by_term.get(term.as_str()).ok_or_else(|| {
                    CfError::new(
                        code::LAND_TERM_UNKNOWN,
                        format!("no provider is registered for term {term}"),
                    )
                })
            })
            .collect()
    }

    /// Registered terms, sorted by id.
    ///
    /// Sorted rather than in insertion order so that a listing is stable across
    /// runs. Nothing ranking-relevant may use this; it is for reporting.
    pub fn term_ids(&self) -> Vec<TermId> {
        let mut ids: Vec<TermId> = self
            .by_term
            .values()
            .map(|r| r.provider.term_id().clone())
            .collect();
        ids.sort_by(|a, b| a.as_str().cmp(b.as_str()));
        ids
    }
}

/// How much a role lets a term affect selection.
///
/// An epsilon constraint sits above a diagnostic because a bound changes which
/// candidates are feasible, and that changes the ranking just as surely as an
/// axis does. A term capped at `diagnostic` must therefore be barred from both.
fn influence_rank(role: TermRole) -> u8 {
    match role {
        TermRole::Diagnostic => 0,
        TermRole::EpsilonConstraint => 1,
        TermRole::ParetoAxis => 2,
    }
}
