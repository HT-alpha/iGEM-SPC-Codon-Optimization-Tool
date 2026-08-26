//! `cf-landscape` — support layer.
//!
//! Composition layer. Owns the term registry, the window-keyed evaluation
//! cache, uncertainty fusion, and the batch-yield protocol driver. Resolves
//! concrete predictors through a runtime registry, so it depends on
//! `cf-landscape-api` and NOT on any predictor crate - this is the
//! cycle-breaker that keeps predictors swappable.
//!
//! Declared dependencies: cf-types, cf-landscape-api, cf-modelreg.
//! Contract: `planning/02-schemas/landscape_response.schema.json`.

#![forbid(unsafe_code)]
