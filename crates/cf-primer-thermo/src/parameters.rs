//! Load a `primer_parameters` document and refuse a digest that does not match.
//!
//! The numbers live in the JSON artifact. This module does not contain a
//! second copy of Table 2.

use std::collections::BTreeMap;

use cf_types::canonical::Digest;
use serde::Deserialize;
use serde_json::Value;

const BUNDLED: &str = include_str!("../data/santalucia1998.json");

#[derive(Debug)]
pub struct ParameterError(pub String);

impl std::fmt::Display for ParameterError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

#[derive(Debug, Clone, Deserialize)]
pub struct Thermo {
    pub delta_h_kcal: f64,
    pub delta_s_cal: f64,
}

#[derive(Debug, Clone)]
pub struct Parameters {
    pub parameter_set_id: String,
    pub citation: String,
    pub doi: String,
    pub temperature_k: f64,
    pub gas_constant_cal: f64,
    pub sodium_molar: f64,
    pub strand_concentration_molar: f64,
    pub salt_dg37_kcal_per_ln_na: f64,
    pub stacks: BTreeMap<String, Thermo>,
    pub init_gc: Thermo,
    pub init_at: Thermo,
    pub symmetry: Thermo,
    pub digest: String,
}

pub fn bundled() -> Result<Parameters, ParameterError> {
    parse(BUNDLED)
}

pub fn parse(text: &str) -> Result<Parameters, ParameterError> {
    let value: Value = serde_json::from_str(text).map_err(|err| ParameterError(err.to_string()))?;
    let digest = value
        .get("digest")
        .and_then(|item| item.as_str())
        .ok_or_else(|| ParameterError("primer parameter document has no digest".into()))?
        .to_string();
    let computed = Digest::of_document_omitting(&value, &["digest", "created_at"])
        .map_err(|err| ParameterError(err.to_string()))?;
    if computed.as_str() != digest {
        return Err(ParameterError(format!(
            "CF-PRIM-005 primer parameter digest {digest} does not match {computed}"
        )));
    }
    let raw: Raw = serde_json::from_value(value).map_err(|err| ParameterError(err.to_string()))?;
    if raw.stacks.len() != 10 {
        return Err(ParameterError(format!(
            "CF-PRIM-005 expected 10 nearest-neighbour stacks, found {}",
            raw.stacks.len()
        )));
    }
    Ok(Parameters {
        parameter_set_id: raw.parameter_set_id,
        citation: raw.provenance.citation,
        doi: raw.provenance.doi.unwrap_or_default(),
        temperature_k: raw.temperature_k,
        gas_constant_cal: raw.gas_constant_cal_per_k_mol,
        sodium_molar: raw.sodium_molar,
        strand_concentration_molar: raw.strand_concentration_molar,
        salt_dg37_kcal_per_ln_na: raw.salt_dg37_kcal_per_ln_na,
        stacks: raw.stacks,
        init_gc: raw.initiation.gc,
        init_at: raw.initiation.at,
        symmetry: raw.symmetry_correction,
        digest,
    })
}

#[derive(Debug, Deserialize)]
struct Raw {
    parameter_set_id: String,
    provenance: RawProvenance,
    temperature_k: f64,
    gas_constant_cal_per_k_mol: f64,
    sodium_molar: f64,
    strand_concentration_molar: f64,
    salt_dg37_kcal_per_ln_na: f64,
    stacks: BTreeMap<String, Thermo>,
    initiation: RawInit,
    symmetry_correction: Thermo,
}

#[derive(Debug, Deserialize)]
struct RawProvenance {
    citation: String,
    doi: Option<String>,
}

#[derive(Debug, Deserialize)]
struct RawInit {
    gc: Thermo,
    at: Thermo,
}
