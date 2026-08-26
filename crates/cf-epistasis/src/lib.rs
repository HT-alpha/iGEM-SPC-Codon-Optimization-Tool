//! `cf-epistasis` — epistatic-surface layer.
//!
//! Synonymous-restricted Potts context model plus the negative-binomial dwell
//! head. Scoring (L1 delta-capable) lives in Rust; the pseudo-likelihood
//! fitter lives in `codonforge.epistasis.train` and is a separate build
//! target that may import `cf-ingest`.
//!
//! Declared dependencies: cf-types, cf-landscape-api.
//! Contract: `planning/03-algorithms/01-epistatic-surface-fit.md`.

#![forbid(unsafe_code)]
