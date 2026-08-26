//! `cf-search` — search-engine layer.
//!
//! Three-tier search engine - tier 1 greedy/DP seeding, tier 2 red-black
//! block-local coordinate ascent, tier 3 NSGA-II/III population loop - plus
//! the systems-level joint search driver. Talks only to `cf-landscape-api`;
//! it never names a concrete predictor.
//!
//! Declared dependencies: cf-types, cf-landscape-api, cf-constraints, cf-metrics.
//! Contract: `planning/03-algorithms/02-block-local-coordinate-ascent.md`.

#![forbid(unsafe_code)]
