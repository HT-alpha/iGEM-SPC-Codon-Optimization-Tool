//! Thermodynamic parameters as input data.
//!
//! Ticket P2-07a. Contract: `schemas/v1/structure_parameters.schema.json`.
//! `DEP-STRU-001`.
//!
//! ## Why there is no table in this file
//!
//! `AMD-009` verified the licence position of the folding *code* and prescribed
//! a permissively-licensed kernel. It said nothing about the nearest-neighbour
//! parameter set that kernel needs to compute anything, and the published sets
//! are embedded in tools whose licences are exactly what the amendment rejects.
//! A kernel shipping a table would either carry values someone invented — which
//! the no-fabricated-constants rule forbids — or redistribute data whose terms
//! nobody here has verified.
//!
//! So this module holds a *reader*, not a table. The set is user-supplied, it
//! carries its own provenance, and its digest goes in the landscape lockfile so
//! two runs cannot silently use different thermodynamics. A run with none is
//! refused.
//!
//! ## What a schema can and cannot do here
//!
//! No schema can tell a real parameter set from a well-formed invented one; the
//! numbers are just numbers. What it *can* do is refuse a set that declines to
//! say where it came from, which is why provenance is required and why the
//! checks below are about citation, licence and source digest rather than about
//! whether any particular energy looks plausible. Judging plausibility would be
//! the fabrication problem wearing a lab coat.

use std::collections::BTreeMap;

use cf_types::{code, CfError, Digest, Result};
use serde_json::Value;

/// Where a parameter set came from.
#[derive(Clone, Debug, PartialEq)]
pub struct Provenance {
    pub citation: String,
    pub doi: Option<String>,
    /// Digest of the file the values were transcribed from.
    pub source_digest: Digest,
    pub licence: String,
    /// The supplier's claim. Recorded, never acted on by vendoring.
    pub redistributable: bool,
    pub transcribed_by: Option<String>,
}

/// A nearest-neighbour parameter set.
#[derive(Clone, Debug, PartialEq)]
pub struct StructureParameters {
    pub parameter_set_id: String,
    pub provenance: Provenance,
    pub temperature_k: f64,
    stack: BTreeMap<String, f64>,
    hairpin_loop: BTreeMap<usize, f64>,
    terminal_au_penalty: f64,
    digest: Digest,
}

impl StructureParameters {
    /// Read a `structure_parameters` document.
    ///
    /// Every refusal below names what is missing and why it matters, because
    /// the person hitting it is being asked to go and find data — a message
    /// that only says "invalid" leaves them guessing at which of a dozen fields
    /// the problem is.
    pub fn from_json(document: &Value) -> Result<Self> {
        let map = document.as_object().ok_or_else(|| {
            CfError::new(
                code::STRU_BAD_PARAMETERS,
                "a structure parameter set must be a JSON object",
            )
        })?;

        let text = |field: &str| -> Result<String> {
            map.get(field)
                .and_then(Value::as_str)
                .map(str::to_string)
                .ok_or_else(|| {
                    CfError::new(
                        code::STRU_BAD_PARAMETERS,
                        format!("structure parameters have no {field}"),
                    )
                })
        };

        let provenance_map = map
            .get("provenance")
            .and_then(Value::as_object)
            .ok_or_else(|| {
                CfError::new(
                    code::STRU_BAD_PARAMETERS,
                    "structure parameters have no provenance block. The schema cannot tell \
                     a real parameter set from a well-formed invented one, so refusing a \
                     set that does not say where it came from is the only check available",
                )
            })?;

        let provenance_text = |field: &str, why: &str| -> Result<String> {
            provenance_map
                .get(field)
                .and_then(Value::as_str)
                .filter(|value| !value.trim().is_empty())
                .map(str::to_string)
                .ok_or_else(|| {
                    CfError::new(
                        code::STRU_BAD_PARAMETERS,
                        format!("structure parameter provenance has no {field}: {why}"),
                    )
                })
        };

        let provenance = Provenance {
            citation: provenance_text(
                "citation",
                "a reader with no network needs the full reference, not a bare identifier",
            )?,
            doi: provenance_map
                .get("doi")
                .and_then(Value::as_str)
                .map(str::to_string),
            source_digest: Digest::parse(&provenance_text(
                "source_digest",
                "without it a transcription cannot be re-checked against what it came from",
            )?)?,
            licence: provenance_text(
                "licence",
                "AMD-009 exists because a licence nobody checked turned out to forbid the \
                 distribution the project had assumed",
            )?,
            redistributable: provenance_map
                .get("redistributable")
                .and_then(Value::as_bool)
                .ok_or_else(|| {
                    CfError::new(
                        code::STRU_BAD_PARAMETERS,
                        "structure parameter provenance does not state whether the set is \
                         redistributable; CodonForge records the claim and never vendors \
                         on the strength of it, but an unstated claim is not a claim",
                    )
                })?,
            transcribed_by: provenance_map
                .get("transcribed_by")
                .and_then(Value::as_str)
                .map(str::to_string),
        };

        let temperature_k = map
            .get("temperature_k")
            .and_then(Value::as_f64)
            .filter(|value| *value > 0.0)
            .ok_or_else(|| {
                CfError::new(
                    code::STRU_BAD_PARAMETERS,
                    "structure parameters have no positive temperature_k; extrapolating \
                     from an unrecorded reference temperature is how a fold becomes \
                     confidently wrong",
                )
            })?;

        let stack: BTreeMap<String, f64> = map
            .get("stack")
            .and_then(Value::as_object)
            .ok_or_else(|| {
                CfError::new(
                    code::STRU_BAD_PARAMETERS,
                    "structure parameters have no stack table",
                )
            })?
            .iter()
            .filter_map(|(key, value)| value.as_f64().map(|v| (key.clone(), v)))
            .collect();
        if stack.is_empty() {
            return Err(CfError::new(
                code::STRU_BAD_PARAMETERS,
                "the stack table is empty; a parameter set with no stacking energies \
                 cannot fold anything",
            ));
        }

        let hairpin_loop = map
            .get("hairpin_loop")
            .and_then(Value::as_object)
            .map(|table| {
                table
                    .iter()
                    .filter_map(|(key, value)| Some((key.parse::<usize>().ok()?, value.as_f64()?)))
                    .collect()
            })
            .unwrap_or_default();

        Ok(Self {
            parameter_set_id: text("parameter_set_id")?,
            provenance,
            temperature_k,
            stack,
            hairpin_loop,
            terminal_au_penalty: map
                .get("terminal_au_penalty")
                .and_then(Value::as_f64)
                .unwrap_or(0.0),
            digest: Digest::parse(&text("digest")?)?,
        })
    }

    /// The digest recorded in the landscape lockfile.
    pub fn digest(&self) -> &Digest {
        &self.digest
    }

    /// Recompute the digest from the document, to check the one it carries.
    ///
    /// A document whose `digest` field disagrees with its content is either
    /// edited or corrupt, and either way the lockfile would record a value that
    /// does not identify what was used.
    pub fn verify_digest(document: &Value) -> Result<Digest> {
        let stated = document
            .get("digest")
            .and_then(Value::as_str)
            .ok_or_else(|| {
                CfError::new(
                    code::STRU_BAD_PARAMETERS,
                    "structure parameters carry no digest, so nothing identifies them in \
                     the lockfile",
                )
            })?;
        let computed = Digest::of_document_omitting(document, &["digest", "created_at"])?;
        if computed.as_str() != stated {
            return Err(CfError::new(
                code::STRU_BAD_PARAMETERS,
                format!(
                    "the parameter set states digest {stated} but its content digests to \
                     {computed}; the lockfile would record a value that does not identify \
                     what was used"
                ),
            ));
        }
        Ok(computed)
    }

    /// A stacking energy, or a refusal naming the pair.
    ///
    /// No default and no interpolation. An absent stack is a gap in the
    /// supplied set, and substituting zero would assert that a stack
    /// contributes nothing — a fabricated constant in the place it does the
    /// most damage, since stacking is most of a structure's stability.
    pub fn stack(&self, pair: &str) -> Result<f64> {
        self.stack.get(pair).copied().ok_or_else(|| {
            CfError::new(
                code::STRU_BAD_PARAMETERS,
                format!(
                    "the supplied parameter set has no stacking energy for {pair}; \
                     substituting zero would assert that this stack contributes nothing"
                ),
            )
        })
    }

    pub fn hairpin(&self, length: usize) -> Result<f64> {
        self.hairpin_loop.get(&length).copied().ok_or_else(|| {
            CfError::new(
                code::STRU_BAD_PARAMETERS,
                format!(
                    "the supplied parameter set has no hairpin initiation energy for a \
                         loop of {length}"
                ),
            )
        })
    }

    pub fn terminal_au_penalty(&self) -> f64 {
        self.terminal_au_penalty
    }

    pub fn stack_keys(&self) -> impl Iterator<Item = &String> {
        self.stack.keys()
    }
}

/// The parameter set a job is running under, or a refusal.
///
/// The refusal is the point. `CF-STRU-001` is what a user sees when they ask
/// for a structure term without supplying parameters, and it names the gap
/// rather than falling back to a table that does not exist.
#[derive(Clone, Debug, Default)]
pub struct ParameterSource {
    parameters: Option<StructureParameters>,
}

impl ParameterSource {
    pub fn none() -> Self {
        Self { parameters: None }
    }

    pub fn from(parameters: StructureParameters) -> Self {
        Self {
            parameters: Some(parameters),
        }
    }

    pub fn require(&self) -> Result<&StructureParameters> {
        self.parameters.as_ref().ok_or_else(|| {
            CfError::new(
                code::STRU_NO_PARAMETERS,
                "this build ships no thermodynamic parameter table and none was supplied. \
                 AMD-009 rules out vendoring the published sets, and inventing values is \
                 worse than refusing: supply a structure_parameters document, whose digest \
                 is then recorded in the landscape lockfile",
            )
        })
    }

    /// The lockfile entry for this source, or `None` when nothing is pinned.
    pub fn lockfile_digest(&self) -> Option<&Digest> {
        self.parameters.as_ref().map(StructureParameters::digest)
    }
}
