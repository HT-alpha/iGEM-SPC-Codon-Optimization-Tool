//! `LandscapeQuery` — the contract between the search engine and the landscape.
//!
//! Ticket P1-20. Contract: `schemas/v1/landscape_query.schema.json`.
//!
//! Three target shapes, because one shape cannot serve both the population loop
//! and the coordinate-ascent loop:
//!
//! | shape | cost | caller |
//! | --- | --- | --- |
//! | [`Target::Full`] | `O(L)` | tier 3, final re-ranking |
//! | [`Target::Delta`] | `O(k·W)` for a decomposable term | tier 2 accept/reject |
//! | [`Target::DeltaTable`] | `O(P·|Syn|·W)` | coordinate ascent, repair |
//!
//! The reason the two delta shapes are distinct rather than one general shape
//! is that they have different callers with different access patterns. Tier 2
//! knows the move it wants to evaluate; coordinate ascent wants every move at a
//! position so it can pick the best one, and issuing `|Syn|` separate delta
//! queries to learn that would re-walk the same window `|Syn|` times.

use cf_types::{code, CfError, Codon, Digest, HostId, Result, TermId};
use serde::{Deserialize, Serialize};

/// Named compute budget. Changes which terms are evaluated and with which
/// engine; see ADR-0008 for the degradation ladder.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum BudgetTier {
    Fast,
    #[default]
    Standard,
    Exact,
}

impl BudgetTier {
    pub fn as_str(self) -> &'static str {
        match self {
            BudgetTier::Fast => "fast",
            BudgetTier::Standard => "standard",
            BudgetTier::Exact => "exact",
        }
    }

    pub fn parse(text: &str) -> Option<Self> {
        match text {
            "fast" => Some(BudgetTier::Fast),
            "standard" => Some(BudgetTier::Standard),
            "exact" => Some(BudgetTier::Exact),
            _ => None,
        }
    }
}

impl std::fmt::Display for BudgetTier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// One synonymous substitution.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub struct Edit {
    pub codon_index: usize,
    pub codon: Codon,
}

impl Edit {
    pub fn new(codon_index: usize, codon: Codon) -> Self {
        Self { codon_index, codon }
    }
}

/// What to score.
#[derive(Clone, Debug)]
pub enum Target {
    /// A complete sequence.
    Full {
        sequence: Box<cf_types::CodonSequence>,
        cistron_id: Option<String>,
    },
    /// A base the provider has already seen, plus substitutions.
    Delta {
        base_digest: Digest,
        edits: Vec<Edit>,
        cistron_id: Option<String>,
    },
    /// Every synonymous alternative at each listed position.
    ///
    /// Omitting a position is how the caller excludes a locked codon; there is
    /// deliberately no "all positions" spelling, because a caller that meant
    /// "all editable positions" and got "all positions" would silently propose
    /// moves at locked sites.
    DeltaTable {
        base_digest: Digest,
        positions: Vec<usize>,
        cistron_id: Option<String>,
    },
}

impl Target {
    pub fn kind(&self) -> &'static str {
        match self {
            Target::Full { .. } => "full",
            Target::Delta { .. } => "delta",
            Target::DeltaTable { .. } => "delta_table",
        }
    }

    pub fn cistron_id(&self) -> Option<&str> {
        match self {
            Target::Full { cistron_id, .. }
            | Target::Delta { cistron_id, .. }
            | Target::DeltaTable { cistron_id, .. } => cistron_id.as_deref(),
        }
    }

    /// The base this target is expressed relative to, if any.
    pub fn base_digest(&self) -> Option<&Digest> {
        match self {
            Target::Full { .. } => None,
            Target::Delta { base_digest, .. } | Target::DeltaTable { base_digest, .. } => {
                Some(base_digest)
            }
        }
    }
}

/// Which host, in which physiological state.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct HostContext {
    pub host_id: HostId,
    /// Pins the exact `HostProfile` document. tRNA tables get corrected over
    /// time, so a result that did not record which table it used is not
    /// reproducible.
    pub host_profile_digest: Digest,
    pub condition_id: Option<String>,
}

impl HostContext {
    pub fn new(host_id: HostId, host_profile_digest: Digest) -> Self {
        Self {
            host_id,
            host_profile_digest,
            condition_id: None,
        }
    }

    pub fn with_condition(mut self, condition_id: impl Into<String>) -> Self {
        self.condition_id = Some(condition_id.into());
        self
    }
}

/// Which output blocks to populate. Requesting less is how the cheap tiers stay
/// cheap.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Outputs {
    pub scalars: bool,
    pub per_position: bool,
    pub deltas: bool,
    pub uncertainty: bool,
    pub violations: bool,
}

impl Default for Outputs {
    fn default() -> Self {
        Self {
            scalars: true,
            per_position: false,
            deltas: false,
            uncertainty: false,
            violations: false,
        }
    }
}

impl Outputs {
    /// Everything the composition layer can produce.
    ///
    /// `gradients` is absent from this type entirely. The response schema
    /// admits a gradient block for a provider that offers one to some other
    /// consumer, but nothing reachable from a `LandscapeQuery` can ask for it,
    /// so no query can put a term on a gradient-only path.
    pub fn all() -> Self {
        Self {
            scalars: true,
            per_position: true,
            deltas: true,
            uncertainty: true,
            violations: true,
        }
    }

    pub fn with_deltas(mut self) -> Self {
        self.deltas = true;
        self
    }

    pub fn with_uncertainty(mut self) -> Self {
        self.uncertainty = true;
        self
    }
}

/// Compute budget for one query.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct Budget {
    pub tier: BudgetTier,
    /// Soft deadline. A provider that cannot meet it yields `partial` with the
    /// terms it did compute; it does not block and it does not guess.
    pub max_latency_ms: Option<u64>,
    pub max_memory_mb: Option<u64>,
}

impl Budget {
    pub fn new(tier: BudgetTier) -> Self {
        Self {
            tier,
            max_latency_ms: None,
            max_memory_mb: None,
        }
    }

    pub fn with_latency_ms(mut self, ms: u64) -> Self {
        self.max_latency_ms = Some(ms);
        self
    }
}

/// Cache partition and read/write intent.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct CachePolicy {
    pub read: bool,
    pub write: bool,
    /// Defaults to the landscape lockfile digest, which is what makes
    /// cross-version contamination structurally impossible rather than
    /// merely unlikely.
    pub namespace: Option<String>,
}

impl Default for CachePolicy {
    fn default() -> Self {
        Self {
            read: true,
            write: true,
            namespace: None,
        }
    }
}

impl CachePolicy {
    pub fn namespaced(namespace: impl Into<String>) -> Self {
        Self {
            read: true,
            write: true,
            namespace: Some(namespace.into()),
        }
    }

    pub fn bypass() -> Self {
        Self {
            read: false,
            write: false,
            namespace: None,
        }
    }
}

/// A per-term artifact pin — the top rung of the ADR-0009 ladder.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct ModelPin {
    pub artifact_id: cf_types::ArtifactId,
    pub digest: Digest,
}

/// A landscape query.
#[derive(Clone, Debug)]
pub struct LandscapeQuery {
    pub query_id: String,
    pub target: Target,
    pub host_context: HostContext,
    /// The construct the sequence sits in. Carried, not interpreted — no
    /// Phase 1 term reads it. See [`crate::json::ConstructRef`].
    pub construct_context: Option<crate::json::ConstructRef>,
    /// Explicit term list. A wildcard is deliberately unsupported: a job whose
    /// objective set depends on which predictors happen to be installed is not
    /// reproducible.
    pub requested_terms: Vec<TermId>,
    pub outputs: Outputs,
    pub budget: Budget,
    pub cache_policy: CachePolicy,
    pub model_pins: Vec<(TermId, ModelPin)>,
    /// RNG stream seed for stochastic providers.
    pub seed: Option<u64>,
}

impl LandscapeQuery {
    pub fn new(
        query_id: impl Into<String>,
        target: Target,
        host_context: HostContext,
        requested_terms: Vec<TermId>,
        budget: Budget,
    ) -> Self {
        Self {
            query_id: query_id.into(),
            target,
            host_context,
            construct_context: None,
            requested_terms,
            outputs: Outputs::default(),
            budget,
            cache_policy: CachePolicy::default(),
            model_pins: Vec::new(),
            seed: None,
        }
    }

    pub fn with_outputs(mut self, outputs: Outputs) -> Self {
        self.outputs = outputs;
        self
    }

    pub fn with_cache_policy(mut self, policy: CachePolicy) -> Self {
        self.cache_policy = policy;
        self
    }

    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = Some(seed);
        self
    }

    /// Structural validation, independent of any registry.
    ///
    /// The seed rule is the one worth reading twice. A sampled uncertainty
    /// request without a seed produces a number that is correct on average and
    /// irreproducible in particular, and irreproducibility surfaces much later
    /// than the query — usually as a replay that will not replay. Refusing at
    /// the query is the only place the cost of the refusal is small.
    pub fn validate(&self) -> Result<()> {
        if self.requested_terms.is_empty() {
            return Err(CfError::new(
                code::LAND_TERM_UNKNOWN,
                "a query must request at least one term",
            ));
        }
        let mut seen = std::collections::HashSet::new();
        for term in &self.requested_terms {
            if !seen.insert(term.as_str()) {
                return Err(CfError::new(
                    code::LAND_TERM_DUPLICATE,
                    format!("term {term} appears twice in requested_terms"),
                ));
            }
        }
        if self.outputs.deltas && !matches!(self.target, Target::DeltaTable { .. }) {
            return Err(CfError::new(
                code::LAND_CAPABILITY,
                format!(
                    "outputs.deltas is only meaningful with a delta_table target, got {}",
                    self.target.kind()
                ),
            ));
        }
        if let Target::DeltaTable { positions, .. } = &self.target {
            if positions.is_empty() {
                return Err(CfError::new(
                    code::LAND_CAPABILITY,
                    "a delta_table target must list at least one position",
                ));
            }
        }
        if let Target::Delta { edits, .. } = &self.target {
            if edits.is_empty() {
                return Err(CfError::new(
                    code::LAND_CAPABILITY,
                    "a delta target must carry at least one edit",
                ));
            }
        }
        Ok(())
    }

    /// The seed a provider must use.
    ///
    /// Deliberately fallible rather than defaulting to zero. A defaulted seed
    /// is reproducible, which is exactly what makes it dangerous: the run would
    /// replay, so nothing would ever reveal that the caller's intent to control
    /// the stream had been silently discarded.
    pub fn effective_seed(&self, provider_is_stochastic: bool) -> Result<u64> {
        match (self.seed, provider_is_stochastic) {
            (Some(seed), _) => Ok(seed),
            (None, false) => Ok(0),
            (None, true) => Err(CfError::new(
                code::LAND_SEED_REQUIRED,
                "a stochastic provider was asked to evaluate a query with no seed",
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cf_types::CodonSequence;

    fn terms(names: &[&str]) -> Vec<TermId> {
        names.iter().map(|n| TermId::new(*n).unwrap()).collect()
    }

    fn host() -> HostContext {
        HostContext::new(
            HostId::new("escherichia_coli:511145/BL21").unwrap(),
            Digest::of_bytes(b"host-profile"),
        )
    }

    fn full_query(term_names: &[&str]) -> LandscapeQuery {
        let seq = CodonSequence::from_cds("ATGGCTTAA", 1).unwrap();
        LandscapeQuery::new(
            "q_01ARZ3NDEKTSV4RRFFQ69G5FAV",
            Target::Full {
                sequence: Box::new(seq),
                cistron_id: None,
            },
            host(),
            terms(term_names),
            Budget::new(BudgetTier::Standard),
        )
    }

    #[test]
    fn unit_query_validation_rejects_the_four_malformed_shapes() {
        assert!(full_query(&["usage.cai"]).validate().is_ok());

        let mut q = full_query(&["usage.cai"]);
        q.requested_terms.clear();
        assert_eq!(q.validate().unwrap_err().code, "CF-LAND-001");

        let q = full_query(&["usage.cai", "usage.cai"]);
        assert_eq!(q.validate().unwrap_err().code, "CF-LAND-005");

        let q = full_query(&["usage.cai"]).with_outputs(Outputs::default().with_deltas());
        assert_eq!(q.validate().unwrap_err().code, "CF-LAND-002");

        let q = LandscapeQuery::new(
            "q_01ARZ3NDEKTSV4RRFFQ69G5FAV",
            Target::DeltaTable {
                base_digest: Digest::of_bytes(b"base"),
                positions: vec![],
                cistron_id: None,
            },
            host(),
            terms(&["usage.cai"]),
            Budget::new(BudgetTier::Fast),
        );
        assert_eq!(q.validate().unwrap_err().code, "CF-LAND-002");
    }

    /// A sampled provider with no seed is refused rather than defaulted.
    #[test]
    fn unit_stochastic_provider_requires_an_explicit_seed() {
        let q = full_query(&["usage.cai"]);
        assert_eq!(q.effective_seed(false).unwrap(), 0);
        assert_eq!(
            q.effective_seed(true).unwrap_err().code,
            "CF-LAND-006",
            "defaulting the seed would make the omission invisible forever"
        );
        assert_eq!(q.with_seed(7).effective_seed(true).unwrap(), 7);
    }

    #[test]
    fn unit_outputs_cannot_request_a_gradient() {
        // Not an assertion about a value: `Outputs` has no gradient field, so
        // this is a compile-time property restated for the reader. `all()`
        // enables every block that exists.
        let all = Outputs::all();
        assert!(all.scalars && all.per_position && all.deltas && all.uncertainty && all.violations);
    }
}
