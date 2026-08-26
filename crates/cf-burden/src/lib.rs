//! `cf-burden` — burden-model layer.
//!
//! Isoacceptor resource ledger, wobble decoding weights, supply/demand
//! fixed-point solver, shadow-price vector, and the saturating burden
//! penalty. L1 delta-capable.
//!
//! Declared dependencies: cf-types, cf-landscape-api.
//! Contract: `planning/03-algorithms/05-burden-penalty-and-capacity.md`.

#![forbid(unsafe_code)]
