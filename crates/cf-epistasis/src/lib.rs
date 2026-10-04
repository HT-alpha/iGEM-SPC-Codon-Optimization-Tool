//! `cf-epistasis` — epistatic-surface layer.
//!
//! Synonymous-restricted Potts context model. Scoring and the
//! pseudo-likelihood fit live here. The dwell head is not in this crate:
//! it needs A-site counts, and a function that returned zeros for those
//! counts would look like a fitted head.
//!
//! Declared dependencies: cf-types, cf-landscape-api.
//! Contract: `planning/03-algorithms/01-epistatic-surface-fit.md`.

#![forbid(unsafe_code)]

mod potts;

pub use potts::{
    fit_potts, score_codon_usage, score_naturalness, FitReport, PottsModel, UsageBaseline,
};
