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

use cf_constraints::{ConstraintChecker, ConstraintSet, ConstructView, KmerMasking};
use cf_landscape_api::{LandscapeQuery, LandscapeResponse};
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

/// Re-parse a `landscape_query.schema.json` document and re-emit it.
///
/// `TST-LAPI-001`. Every field the contract models survives; a field the
/// contract does not model is the interesting case, and the Python test asserts
/// on which ones are preserved verbatim (`construct_context`) and which are
/// refused (`outputs.gradients`).
#[pyfunction]
fn landscape_query_round_trip(document_json: &str) -> PyResult<String> {
    let value: serde_json::Value =
        serde_json::from_str(document_json).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let query = LandscapeQuery::from_json(&value).map_err(to_py_err)?;
    query.to_canonical_json().map_err(to_py_err)
}

/// Re-parse a `landscape_response.schema.json` document and re-emit it.
#[pyfunction]
fn landscape_response_round_trip(document_json: &str) -> PyResult<String> {
    let value: serde_json::Value =
        serde_json::from_str(document_json).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let response = LandscapeResponse::from_json(&value).map_err(to_py_err)?;
    response.to_canonical_json().map_err(to_py_err)
}

fn parse_masking(masking: &str) -> PyResult<KmerMasking> {
    match masking {
        "assembled" => Ok(KmerMasking::Assembled),
        "segmented_by_region" => Ok(KmerMasking::SegmentedByRegion),
        other => Err(PyValueError::new_err(format!(
            "unknown k-mer masking {other:?}; expected assembled or segmented_by_region"
        ))),
    }
}

/// Feasibility of many constructs against one constraint set.
///
/// Batched on purpose. Compiling the constraint set builds an Aho-Corasick
/// automaton, and doing that once per sequence would dominate the measurement;
/// batching also keeps this entry point shaped like the batch-yield protocol
/// the search engine uses, rather than inviting a per-move call.
#[pyfunction]
#[pyo3(signature = (constraint_set_json, constructs, cds_spans, masking="assembled"))]
fn constraint_feasible_batch(
    constraint_set_json: &str,
    constructs: Vec<String>,
    cds_spans: Vec<(usize, usize)>,
    masking: &str,
) -> PyResult<Vec<bool>> {
    let value: serde_json::Value = serde_json::from_str(constraint_set_json)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let set = ConstraintSet::from_json(&value).map_err(to_py_err)?;
    let checker = ConstraintChecker::compile(&set).map_err(to_py_err)?;
    let masking = parse_masking(masking)?;

    Ok(constructs
        .iter()
        .map(|construct| {
            let view =
                ConstructView::new(construct.as_bytes(), cds_spans.clone()).with_masking(masking);
            checker.is_feasible(&view)
        })
        .collect())
}

/// Full violation report for one construct, as JSON.
#[pyfunction]
#[pyo3(signature = (constraint_set_json, construct, cds_spans, masking="assembled"))]
fn constraint_violations(
    constraint_set_json: &str,
    construct: &str,
    cds_spans: Vec<(usize, usize)>,
    masking: &str,
) -> PyResult<String> {
    let value: serde_json::Value = serde_json::from_str(constraint_set_json)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let set = ConstraintSet::from_json(&value).map_err(to_py_err)?;
    let checker = ConstraintChecker::compile(&set).map_err(to_py_err)?;
    let view =
        ConstructView::new(construct.as_bytes(), cds_spans).with_masking(parse_masking(masking)?);

    let report: Vec<serde_json::Value> = checker
        .check(&view)
        .into_iter()
        .map(|v| {
            serde_json::json!({
                "constraint_id": v.constraint_id,
                "class": format!("{:?}", v.class).to_lowercase(),
                "start_nt": v.start_nt,
                "end_nt": v.end_nt,
                "detail": v.detail,
                "actionable": v.actionable,
            })
        })
        .collect();
    Ok(serde_json::json!({
        "violations": report,
        "undetectable": checker.undetectable(),
    })
    .to_string())
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
    m.add_function(wrap_pyfunction!(constraint_feasible_batch, m)?)?;
    m.add_function(wrap_pyfunction!(constraint_violations, m)?)?;
    m.add_function(wrap_pyfunction!(landscape_query_round_trip, m)?)?;
    m.add_function(wrap_pyfunction!(landscape_response_round_trip, m)?)?;
    Ok(())
}
