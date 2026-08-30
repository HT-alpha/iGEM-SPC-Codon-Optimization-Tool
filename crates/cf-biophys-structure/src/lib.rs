//! `cf-biophys-structure` — biophysical-predictors layer.
//!
//! mRNA secondary-structure objective. Ships a permissive-licensed
//! windowed-accessibility kernel as the default engine and exposes a
//! `StructureEngine` plugin seam for ViennaRNA / LinearFold back ends that
//! the user installs separately (see RISK-DEP-001).
//!
//! The kernel itself is blocked on parameter data (`DEP-STRU-001`): what ships
//! today is the parameter *reader*, so that the data question is answered
//! before a kernel exists to answer it wrongly.
//!
//! Declared dependencies: cf-types, cf-landscape-api.
//! Contract: `planning/02-schemas/interfaces/predictors.pyi`.

#![forbid(unsafe_code)]
#![warn(missing_debug_implementations)]

pub mod parameters;

pub use parameters::{ParameterSource, Provenance, StructureParameters};
