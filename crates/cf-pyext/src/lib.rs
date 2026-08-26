//! PyO3 glue producing the `codonforge._core` extension module.
//!
//! Ticket P1-01.
//!
//! THIS FILE IS THE ENTIRE RUST/PYTHON BOUNDARY. ADR-0001 requires that no call
//! cross it more than once per generation, and the way that requirement is kept
//! honest is by keeping the crossing surface small enough to read in one
//! sitting. `tools/check_module_layering.py --check` fails the build if any
//! other crate acquires a pyo3 dependency.
//!
//! Everything crosses as a JSON string rather than as a converted Python
//! object. That is deliberate on two counts: it keeps the binding free of a
//! type-conversion layer that would need its own contract tests, and it makes
//! every boundary payload trivially dumpable — so a failure in the search
//! engine is reproducible from a captured batch, which is the property the
//! batch-yield protocol depends on.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use cf_types::{canonicalize, CodonSequence, Digest, GeneticCode};

fn to_py_err(err: cf_types::CfError) -> PyErr {
    // The Python side gets the RFC 9457 problem detail, not a flattened
    // message, so an error crossing the boundary keeps its code.
    PyValueError::new_err(err.to_problem_details().to_string())
}

/// Version of the compiled core.
#[pyfunction]
fn version() -> &'static str {
    cf_types::VERSION
}

/// NCBI translation table ids that are actually defined.
///
/// Exposed because the range 1..=33 has six retired ids, and a Python caller
/// validating user input needs the real set rather than the range.
#[pyfunction]
fn defined_translation_tables() -> Vec<u8> {
    cf_types::tables::defined_ids()
}

/// Translate a coding sequence.
#[pyfunction]
fn translate(cds: &str, translation_table: u8) -> PyResult<String> {
    let seq = CodonSequence::from_cds(cds, translation_table).map_err(to_py_err)?;
    Ok(seq.protein_string())
}

/// Build a `sequence.schema.json` document from a coding sequence.
///
/// Returns a JSON string. The Python contract tests validate the result against
/// the schema, which is how the Rust and Python views of the contract are kept
/// from drifting (`TST-LAPI-001`).
#[pyfunction]
#[pyo3(signature = (cds, translation_table, five_prime_flank=None, three_prime_flank=None, sequence_id=None))]
fn sequence_document(
    cds: &str,
    translation_table: u8,
    five_prime_flank: Option<&str>,
    three_prime_flank: Option<&str>,
    sequence_id: Option<&str>,
) -> PyResult<String> {
    let mut seq = CodonSequence::from_cds(cds, translation_table).map_err(to_py_err)?;
    if let Some(flank) = five_prime_flank {
        seq = seq.with_five_prime_flank(flank).map_err(to_py_err)?;
    }
    if let Some(flank) = three_prime_flank {
        seq = seq.with_three_prime_flank(flank).map_err(to_py_err)?;
    }
    if let Some(id) = sequence_id {
        seq = seq.with_sequence_id(id);
    }
    let doc = seq.to_json().map_err(to_py_err)?;
    Ok(doc.to_string())
}

/// Re-parse a `sequence.schema.json` document, re-checking every invariant.
///
/// Returns the round-tripped document, so a Python test can assert that a
/// document survives Rust -> JSON -> Rust -> JSON unchanged.
#[pyfunction]
fn sequence_round_trip(document_json: &str) -> PyResult<String> {
    let value: serde_json::Value =
        serde_json::from_str(document_json).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let seq = CodonSequence::from_json(&value).map_err(to_py_err)?;
    Ok(seq.to_json().map_err(to_py_err)?.to_string())
}

/// RFC 8785 canonical form of a JSON document.
#[pyfunction]
fn canonical_json(document_json: &str) -> PyResult<String> {
    let value: serde_json::Value =
        serde_json::from_str(document_json).map_err(|e| PyValueError::new_err(e.to_string()))?;
    canonicalize(&value).map_err(to_py_err)
}

/// sha256 digest over the RFC 8785 canonical form, optionally omitting a field.
#[pyfunction]
#[pyo3(signature = (document_json, omit_field=None))]
fn digest_json(document_json: &str, omit_field: Option<&str>) -> PyResult<String> {
    let value: serde_json::Value =
        serde_json::from_str(document_json).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let digest = match omit_field {
        Some(field) => Digest::of_document(&value, field),
        None => Digest::of_json(&value),
    }
    .map_err(to_py_err)?;
    Ok(digest.to_string())
}

/// Synonymous alternatives for a codon under a translation table.
#[pyfunction]
fn synonyms(codon: &str, translation_table: u8) -> PyResult<Vec<String>> {
    let table = GeneticCode::get(translation_table).map_err(to_py_err)?;
    let codon = cf_types::Codon::parse(codon).map_err(to_py_err)?;
    Ok(table
        .synonyms(codon)
        .iter()
        .map(|c| c.to_string_upper())
        .collect())
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__doc__", "CodonForge compiled core (cf-pyext).")?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(defined_translation_tables, m)?)?;
    m.add_function(wrap_pyfunction!(translate, m)?)?;
    m.add_function(wrap_pyfunction!(sequence_document, m)?)?;
    m.add_function(wrap_pyfunction!(sequence_round_trip, m)?)?;
    m.add_function(wrap_pyfunction!(canonical_json, m)?)?;
    m.add_function(wrap_pyfunction!(digest_json, m)?)?;
    m.add_function(wrap_pyfunction!(synonyms, m)?)?;
    Ok(())
}
