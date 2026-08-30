//! `cf-biophys-structure` — biophysical-predictors layer.
//!
//! mRNA secondary-structure objective. Ships a permissive-licensed
//! windowed-accessibility kernel as the default engine and exposes a
//! `StructureEngine` plugin seam for ViennaRNA / LinearFold back ends that
//! the user installs separately (see RISK-DEP-001).
//!
//! `AMD-013` re-scoped the kernel: rather than implementing folding from
//! scratch, the default engine binds to a user-installed ViennaRNA, whose
//! licence permits use and restricts only charging a fee for redistribution.
//! What ships today is the parameter *reader*, which survives the correction
//! because recording which parameter set produced a run is a reproducibility
//! requirement independent of who may distribute it.
//!
//! Declared dependencies: cf-types, cf-landscape-api.
//! Contract: `planning/02-schemas/interfaces/predictors.pyi`.

#![forbid(unsafe_code)]
#![warn(missing_debug_implementations)]

pub mod parameters;

pub use parameters::{ParameterSource, Provenance, StructureParameters};
