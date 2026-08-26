//! `cf-landscape-api` — support layer.
//!
//! The contract between the search engine and every model that scores a
//! sequence. Defines LandscapeQuery / LandscapeResponse, the predictor
//! capability levels (L0 scalar / L1 delta / L2 gradient), the batch-yield
//! evaluation protocol, and the TermId registry. Contains trait/Protocol
//! definitions and serialization only - never an implementation.
//!
//! Declared dependencies: cf-types.
//! Contract: `planning/02-schemas/landscape_query.schema.json`.

#![forbid(unsafe_code)]
