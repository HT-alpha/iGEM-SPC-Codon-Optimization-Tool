//! `cf-burden` — burden-model layer.
//!
//! Isoacceptor resource ledger, wobble decoding weights, supply/demand
//! fixed-point solver, shadow-price vector, and the saturating burden
//! penalty. L1 delta-capable.
//!
//! Declared dependencies: cf-types, cf-landscape-api.
//! Contract: `planning/03-algorithms/05-burden-penalty-and-capacity.md`.

#![forbid(unsafe_code)]
#![warn(missing_debug_implementations)]

pub mod capacity;
pub mod ledger;
pub mod penalty;
pub mod wobble;

pub use capacity::{
    fit_capacity_tier1, fit_capacity_tier2, fit_capacity_tier3, CapacityEstimate, DataCompleteness,
    ThetaUnit,
};
pub use ledger::DemandLedger;
pub use penalty::{growth_cost, power_mean, GrowthCostTerm, AGGREGATION_P, SATURATION_N};
pub use wobble::{decoders, primary_decoder, Decoder, DecodingWeights, Pairing, CRICK_WOBBLE_V1};
