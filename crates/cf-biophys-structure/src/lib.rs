//! `cf-biophys-structure` — biophysical-predictors layer.
//!
//! mRNA secondary-structure objective. Ships a permissive-licensed
//! windowed-accessibility kernel as the default engine and exposes a
//! `StructureEngine` plugin seam for ViennaRNA / LinearFold back ends that
//! the user installs separately (see RISK-DEP-001).
//!
//! Declared dependencies: cf-types, cf-landscape-api.
//! Contract: `planning/02-schemas/interfaces/predictors.pyi`.

#![forbid(unsafe_code)]
