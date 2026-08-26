//! `LandscapeResponse` — what a landscape query returns.
//!
//! Ticket P1-20. Contract: `schemas/v1/landscape_response.schema.json`.
//!
//! Two properties of this type carry most of its weight.
//!
//! **`UNKNOWN` is not a number.** A term absent from [`LandscapeResponse::terms`]
//! under `status: partial` means "not computed", which is different from any
//! value it might plausibly have had. Imputing a mean or a worst case would
//! bias the front invisibly, so the objective vector holds `Option<f64>` and
//! nothing in this crate fills one in.
//!
//! **Ordering is carried, never inferred.** `objective_vector.order` is
//! explicit. Two runs whose objective order differed would produce fronts that
//! are not comparable while looking identical in the output — a failure with no
//! symptom until someone compares two runs and gets nonsense.

use cf_types::{ArtifactId, Digest, Direction, TermId, Unit};
use serde::{Deserialize, Serialize};

use crate::capability::Capability;
use crate::role::ObjectiveAxis;

/// Outcome of a query.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Status {
    Ok,
    /// A first-class outcome, not an error: a provider that missed its latency
    /// budget returned what it had.
    Partial,
    Error,
}

impl Status {
    pub fn as_str(self) -> &'static str {
        match self {
            Status::Ok => "ok",
            Status::Partial => "partial",
            Status::Error => "error",
        }
    }
}

/// Which rung of the ADR-0009 ladder selected an artifact.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ResolvedBy {
    Pin,
    TenantChannel,
    HostChannel,
    NearestHost,
    Base,
    /// Not an artifact at all — CAI, GC, k-mer counts, the native structure
    /// kernel. Still recorded, with the engine version, because "computed by
    /// code we shipped" is provenance too.
    Builtin,
}

impl ResolvedBy {
    pub fn as_str(self) -> &'static str {
        match self {
            ResolvedBy::Pin => "pin",
            ResolvedBy::TenantChannel => "tenant_channel",
            ResolvedBy::HostChannel => "host_channel",
            ResolvedBy::NearestHost => "nearest_host",
            ResolvedBy::Base => "base",
            ResolvedBy::Builtin => "builtin",
        }
    }

    pub fn parse(text: &str) -> Option<Self> {
        match text {
            "pin" => Some(ResolvedBy::Pin),
            "tenant_channel" => Some(ResolvedBy::TenantChannel),
            "host_channel" => Some(ResolvedBy::HostChannel),
            "nearest_host" => Some(ResolvedBy::NearestHost),
            "base" => Some(ResolvedBy::Base),
            "builtin" => Some(ResolvedBy::Builtin),
            _ => None,
        }
    }
}

/// Which artifact produced a value. Per term, not per response, because one
/// response mixes terms from many artifacts.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Provenance {
    pub resolved_by: ResolvedBy,
    pub artifact_id: Option<ArtifactId>,
    pub digest: Option<Digest>,
    pub engine: Option<String>,
    pub engine_version: Option<String>,
    pub fallback: bool,
}

impl Provenance {
    /// Provenance for a term computed by shipped code rather than an artifact.
    pub fn builtin(engine: impl Into<String>, engine_version: impl Into<String>) -> Self {
        Self {
            resolved_by: ResolvedBy::Builtin,
            artifact_id: None,
            digest: None,
            engine: Some(engine.into()),
            engine_version: Some(engine_version.into()),
            fallback: false,
        }
    }
}

/// How a term's uncertainty was produced.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum UncertaintyMethod {
    DeepEnsemble,
    McDropout,
    BayesianLastLayer,
    GlmAsymptotic,
    ParametricBootstrap,
    SplitConformal,
    MondrianConformal,
    /// The term is deterministic. Legal, and distinct from "unknown".
    None,
}

impl UncertaintyMethod {
    pub fn as_str(self) -> &'static str {
        match self {
            UncertaintyMethod::DeepEnsemble => "deep_ensemble",
            UncertaintyMethod::McDropout => "mc_dropout",
            UncertaintyMethod::BayesianLastLayer => "bayesian_last_layer",
            UncertaintyMethod::GlmAsymptotic => "glm_asymptotic",
            UncertaintyMethod::ParametricBootstrap => "parametric_bootstrap",
            UncertaintyMethod::SplitConformal => "split_conformal",
            UncertaintyMethod::MondrianConformal => "mondrian_conformal",
            UncertaintyMethod::None => "none",
        }
    }

    /// Whether evaluating this method consumes the RNG stream.
    ///
    /// Drives the seed requirement: a sampled method with no seed produces a
    /// result that cannot be replayed.
    pub fn is_sampled(self) -> bool {
        matches!(
            self,
            UncertaintyMethod::DeepEnsemble
                | UncertaintyMethod::McDropout
                | UncertaintyMethod::ParametricBootstrap
        )
    }
}

/// Predictive uncertainty for one term.
#[derive(Clone, PartialEq, Debug)]
pub struct Uncertainty {
    pub method: UncertaintyMethod,
    pub calibrated: bool,
    pub sigma: Option<f64>,
    pub interval: Option<(f64, f64)>,
    pub coverage: Option<f64>,
    pub n_members: Option<u32>,
    /// Fraction of predictive **variance** shared across candidates scored by
    /// the same artifact. See `common.schema.json`; `None` means the
    /// decomposition was never measured and must not be read as zero.
    pub shared_fraction: Option<f64>,
}

impl Uncertainty {
    /// The uncertainty of a term that has none.
    ///
    /// A point mass, so the shared/idiosyncratic decomposition is degenerate
    /// rather than unmeasured: there is no error to apportion. Recording `1.0`
    /// versus `0.0` here would be a distinction without a difference at
    /// `sigma == 0`, and `0.0` is chosen only because it is the value the
    /// sampler will multiply by zero anyway.
    pub fn deterministic() -> Self {
        Self {
            method: UncertaintyMethod::None,
            calibrated: true,
            sigma: Some(0.0),
            interval: None,
            coverage: None,
            n_members: None,
            shared_fraction: Some(0.0),
        }
    }

    /// Whether this term may enter ADR-0005 probabilistic dominance.
    ///
    /// Requires a measured decomposition *and* a calibrated interval. Both
    /// conditions fail in the flattering direction if waived: an unmeasured
    /// decomposition read as independence, or an uncalibrated sigma taken at
    /// face value, each make the front look better resolved than the evidence
    /// supports.
    pub fn is_dominance_ready(&self) -> bool {
        if self.method == UncertaintyMethod::None && self.sigma == Some(0.0) {
            return true;
        }
        self.calibrated && self.shared_fraction.is_some() && self.sigma.is_some()
    }
}

/// One term's scalar value and everything needed to interpret it.
#[derive(Clone, PartialEq, Debug)]
pub struct TermValue {
    pub value: f64,
    pub unit: Unit,
    pub direction: Direction,
    pub capability: Capability,
    pub provenance: Provenance,
    pub uncertainty: Option<Uncertainty>,
    pub reference_range: Option<(f64, f64)>,
}

impl TermValue {
    pub fn new(
        value: f64,
        unit: Unit,
        direction: Direction,
        capability: Capability,
        provenance: Provenance,
    ) -> Self {
        Self {
            value,
            unit,
            direction,
            capability,
            provenance,
            uncertainty: None,
            reference_range: None,
        }
    }

    pub fn with_uncertainty(mut self, uncertainty: Uncertainty) -> Self {
        self.uncertainty = Some(uncertainty);
        self
    }
}

/// The exact change in each term for one synonymous alternative.
#[derive(Clone, PartialEq, Debug, Default)]
pub struct AlternativeDelta {
    pub codon: Option<cf_types::Codon>,
    /// Per-term change. A term missing from this map means the provider is L0
    /// for it — not that the change is zero.
    pub delta: Vec<(TermId, f64)>,
}

impl AlternativeDelta {
    pub fn get(&self, term: &TermId) -> Option<f64> {
        self.delta.iter().find(|(t, _)| t == term).map(|(_, d)| *d)
    }
}

/// Every alternative at one position.
#[derive(Clone, PartialEq, Debug)]
pub struct DeltaEntry {
    pub codon_index: usize,
    pub alternatives: Vec<AlternativeDelta>,
}

/// The `deltas` block of a response.
#[derive(Clone, PartialEq, Debug)]
pub struct Deltas {
    pub base_digest: Digest,
    /// True when every delta in this block is closed-form. The repair operator
    /// requires it before it may claim minimality.
    pub exact: bool,
    pub entries: Vec<DeltaEntry>,
}

/// Counters and notes. Not decorative: `terms_deferred` is how a caller learns
/// which terms it must not compare on.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct Diagnostics {
    pub cache_hits: u64,
    pub cache_misses: u64,
    pub evaluations: u64,
    /// Full rescores performed because a provider could not answer a delta
    /// query in closed form. Reported because it is the quantity that
    /// distinguishes per-term degradation from per-query degradation, and a
    /// budget regression here is invisible in wall time on small inputs.
    pub full_rescores: u64,
    pub terms_deferred: Vec<TermId>,
    pub engine_notes: Vec<String>,
}

/// The objective vector, assembled only from axes.
///
/// Wraps `cf_types::ObjectiveVector` and adds the role gate: the constructor
/// takes [`ObjectiveAxis`] values, which cannot be built from a diagnostic
/// term. A caller holding a bare `TermId` cannot reach this constructor.
#[derive(Clone, PartialEq, Debug)]
pub struct ObjectiveVector(cf_types::ObjectiveVector);

impl ObjectiveVector {
    pub fn new(
        axes: Vec<ObjectiveAxis>,
        sense: Vec<Direction>,
        values: Vec<Option<f64>>,
    ) -> std::result::Result<Self, &'static str> {
        let order: Vec<TermId> = axes.into_iter().map(ObjectiveAxis::into_term).collect();
        cf_types::ObjectiveVector::new(order, sense, values).map(Self)
    }

    pub fn inner(&self) -> &cf_types::ObjectiveVector {
        &self.0
    }

    pub fn order(&self) -> &[TermId] {
        self.0.order()
    }

    pub fn values(&self) -> &[Option<f64>] {
        self.0.values()
    }

    pub fn sense(&self) -> &[Direction] {
        self.0.sense()
    }

    pub fn has_unknown(&self) -> bool {
        self.0.has_unknown()
    }
}

/// A landscape response.
#[derive(Clone, PartialEq, Debug)]
pub struct LandscapeResponse {
    pub query_id: String,
    pub status: Status,
    pub terms: Vec<(TermId, TermValue)>,
    pub objective_vector: Option<ObjectiveVector>,
    pub deltas: Option<Deltas>,
    pub diagnostics: Diagnostics,
    pub problem: Option<serde_json::Value>,
}

impl LandscapeResponse {
    pub fn get(&self, term: &TermId) -> Option<&TermValue> {
        self.terms.iter().find(|(t, _)| t == term).map(|(_, v)| v)
    }

    /// Terms that were requested but not computed.
    pub fn deferred(&self) -> &[TermId] {
        &self.diagnostics.terms_deferred
    }

    /// Terms whose uncertainty cannot enter probabilistic dominance.
    ///
    /// Reported rather than silently sampled under an assumption. Empty is the
    /// normal answer while every term is deterministic; it stops being empty
    /// the moment a learned predictor arrives without a measured
    /// shared/idiosyncratic split, and that is exactly when someone needs to
    /// hear about it.
    pub fn dominance_blocked_terms(&self) -> Vec<TermId> {
        self.terms
            .iter()
            .filter(|(_, v)| match &v.uncertainty {
                None => true,
                Some(u) => !u.is_dominance_ready(),
            })
            .map(|(t, _)| t.clone())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cf_types::TermId;

    fn term(name: &str) -> TermId {
        TermId::new(name).unwrap()
    }

    #[test]
    fn unit_status_and_resolved_by_round_trip() {
        assert_eq!(Status::Partial.as_str(), "partial");
        for r in [
            ResolvedBy::Pin,
            ResolvedBy::TenantChannel,
            ResolvedBy::HostChannel,
            ResolvedBy::NearestHost,
            ResolvedBy::Base,
            ResolvedBy::Builtin,
        ] {
            assert_eq!(ResolvedBy::parse(r.as_str()), Some(r));
        }
    }

    /// Absence of a measured split blocks dominance; it is never read as zero.
    #[test]
    fn unit_unmeasured_shared_fraction_blocks_probabilistic_dominance() {
        assert!(Uncertainty::deterministic().is_dominance_ready());

        let unmeasured = Uncertainty {
            method: UncertaintyMethod::DeepEnsemble,
            calibrated: true,
            sigma: Some(0.2),
            interval: None,
            coverage: Some(0.9),
            n_members: Some(5),
            shared_fraction: None,
        };
        assert!(
            !unmeasured.is_dominance_ready(),
            "an unmeasured decomposition must not be read as independence"
        );

        let measured = Uncertainty {
            shared_fraction: Some(0.6),
            ..unmeasured.clone()
        };
        assert!(measured.is_dominance_ready());

        let uncalibrated = Uncertainty {
            calibrated: false,
            ..measured.clone()
        };
        assert!(!uncalibrated.is_dominance_ready());
    }

    #[test]
    fn unit_objective_vector_is_built_from_axes_only() {
        use crate::role::TermRole;
        let axes = vec![
            ObjectiveAxis::new(term("usage.cai"), TermRole::ParetoAxis).unwrap(),
            ObjectiveAxis::new(term("usage.gc_content"), TermRole::ParetoAxis).unwrap(),
        ];
        let v = ObjectiveVector::new(
            axes,
            vec![Direction::Maximize, Direction::Minimize],
            vec![Some(0.8), None],
        )
        .unwrap();
        assert_eq!(v.order(), &[term("usage.cai"), term("usage.gc_content")]);
        assert!(v.has_unknown());
    }
}
