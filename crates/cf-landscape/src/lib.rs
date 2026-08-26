//! `cf-landscape` — composition layer.
//!
//! Tickets P1-21 and P1-22.
//!
//! Owns the term registry resolution, the window-keyed evaluation cache, and
//! the three built-in terms. Does **not** own any predictor: it depends on
//! `cf-landscape-api` and never on `cf-epistasis`, `cf-biophys-structure`,
//! `cf-biophys-kinetics`, `cf-biophys-lm` or `cf-burden`. That is the
//! cycle-breaker that keeps predictors swappable, and `TST-LAND-003` checks it
//! against the real dependency tree rather than trusting this comment.
//!
//! It does depend on `cf-constraints`, because `constraint.violation_count` is
//! a built-in term and counting violations is that crate's job. `cf-constraints`
//! is a first-party core layer, not a swappable predictor.
//!
//! Declared dependencies: cf-types, cf-landscape-api, cf-constraints.
//! Contract: `schemas/v1/landscape_response.schema.json`.

#![forbid(unsafe_code)]
#![warn(missing_debug_implementations)]

pub mod cache;
pub mod compose;
pub mod construct;
pub mod terms;

pub use cache::{CacheKey, EvaluationCache};
pub use compose::{planned_capability, Composer, HostCapability};
pub use construct::ConstructLayout;
pub use terms::{CaiTerm, CodonUsageTable, GcContentTerm, UsageSubset, ViolationCountTerm};

/// Version reported in the provenance of every built-in term.
///
/// Built-ins have no artifact to pin, so the engine version is the only thing
/// that identifies what computed a value. A change to a built-in's definition
/// must bump it, or a lockfile replay would reproduce the coordinates of the
/// old computation and the numbers of the new one.
pub const ENGINE_VERSION: &str = env!("CARGO_PKG_VERSION");
