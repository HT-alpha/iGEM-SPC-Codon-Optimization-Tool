//! `cf-metrics` — support layer.
//!
//! Pure multi-objective indicator math shared by the search engine and the
//! decision-support layer: non-dominated sorting, crowding distance,
//! Das-Dennis reference directions, WFG hypervolume, IGD+, epsilon indicator.
//! Deliberately a leaf so that `cf-search` and `cf-decision` never need to
//! depend on each other.
//!
//! Declared dependencies: none.
//! Contract: `planning/03-algorithms/03-nsga-loop-and-operators.md`.

#![forbid(unsafe_code)]
#![warn(missing_debug_implementations)]

pub mod dominance;
pub mod indicators;

pub use dominance::{
    constraint_dominance, crowding_distance, fast_nondominated_sort, Dominance, Ranked,
};
pub use indicators::{das_dennis, das_dennis_two_layer, epsilon_additive, hypervolume, igd_plus};

/// This crate's error type.
///
/// Local, because `cf-metrics` is rank 0 in `modules.yaml` with no declared
/// dependencies — deliberately a leaf so that `cf-search` and `cf-decision`
/// never need to depend on each other through it. Reaching for
/// `cf_types::CfError` would acquire a dependency for the sake of an error
/// string and the module-layering check would reject it.
///
/// The `CF-METR-*` codes are duplicated in `cf-types`' registry so the rest of
/// the system can name them. That is a duplicated *string*, not duplicated
/// logic, and the same arrangement `jcs.py` has with
/// `VOLATILE_DIGEST_FIELDS`: two independent implementations whose values must
/// match, checked by `tools/check_schema_sync.py`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MetricsError {
    pub code: &'static str,
    pub detail: String,
}

impl MetricsError {
    /// `CF-METR-001`: objective counts that do not line up.
    pub fn dimension(detail: impl Into<String>) -> Self {
        Self {
            code: "CF-METR-001",
            detail: detail.into(),
        }
    }
}

impl std::fmt::Display for MetricsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.code, self.detail)
    }
}

impl std::error::Error for MetricsError {}

pub type Result<T> = std::result::Result<T, MetricsError>;
