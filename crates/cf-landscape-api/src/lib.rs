//! `cf-landscape-api` — the landscape contract.
//!
//! Ticket P1-20. Contracts: `schemas/v1/landscape_query.schema.json` and
//! `schemas/v1/landscape_response.schema.json`.
//!
//! Types and one trait. No predictor, no composition, no cache. That emptiness
//! is the point: this crate is the cycle-breaker that lets `cf-landscape`
//! compose terms it has never heard of, and it stops being one the moment it
//! knows about a concrete term.
//!
//! ## What the contract guarantees
//!
//! ADR-0002 fixes the landscape as a **black-box oracle with an optional exact
//! single-substitution delta**. Concretely:
//!
//! * every provider can score a full sequence (L0);
//! * a decomposable provider can also return the exact change for a move (L1),
//!   which is what makes coordinate ascent affordable — measured at ~640x
//!   cheaper than rescoring at 2 kb and ~4800x at 15 kb;
//! * nothing requires a gradient, and [`TermProvider`] has no method that could
//!   return one.
//!
//! ## Three invariants that are not shape properties
//!
//! `TST-LAPI-003` exists because the first two contract tests check that
//! documents round-trip and that capability claims match method surfaces, and
//! composition depends on three further properties that neither one touches:
//!
//! 1. **Ordering.** [`TermRegistry::resolve`] returns the caller's order. The
//!    objective vector is built from that order, never from map iteration.
//! 2. **Seed determinism.** Two evaluations of one query under one seed agree
//!    exactly, which is what makes the evaluation cache sound.
//! 3. **Term independence.** Evaluating a set of terms in any permutation gives
//!    identical per-term values.
//!
//! Declared dependencies: cf-types.

#![forbid(unsafe_code)]
#![warn(missing_debug_implementations)]

pub mod capability;
pub mod json;
pub mod provider;
pub mod query;
pub mod response;
pub mod role;

pub use capability::Capability;
pub use json::ConstructRef;
pub use provider::{
    verify_capability, EvalContext, PositionDeltas, Registration, TermProvider, TermRegistry,
};
pub use query::{
    Budget, BudgetTier, CachePolicy, Edit, HostContext, LandscapeQuery, ModelPin, Outputs, Target,
};
pub use response::{
    AlternativeDelta, DeltaEntry, Deltas, Diagnostics, LandscapeResponse, ObjectiveVector,
    Provenance, ResolvedBy, Status, TermValue, Uncertainty, UncertaintyMethod,
};
pub use role::{ObjectiveAxis, TermRole};

/// Major version of the landscape contract.
///
/// A job may not mix majors across terms (`CF-RESOLVE-002`): two terms that
/// disagree about what a query means cannot be composed into one vector.
pub const CONTRACT_MAJOR: u32 = 1;

/// Schema version of the query and response documents.
pub const SCHEMA_VERSION: &str = "1.0.0";
