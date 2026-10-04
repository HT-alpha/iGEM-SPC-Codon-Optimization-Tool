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

// Clippy 1.98 flags `?` on every `PyResult<T>` as converting `PyErr` into
// itself (`useless_conversion`). That is PyO3's type alias, not a real
// conversion we added. Local 1.83 does not emit this.
#![allow(clippy::useless_conversion)]

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use cf_constraints::legacy::{self, LegacyScenario};
use cf_constraints::{ConstraintChecker, ConstraintSet, ConstructView, KmerMasking};
use cf_landscape_api::{LandscapeQuery, LandscapeResponse};
use cf_search::{
    BaselineBudget, CodonWeights, Genotype, JobCheckpoint, OperatorMix, SearchRng, SegmentRequest,
};
use cf_types::{canonicalize, CodonSequence, Digest, GeneticCode};

fn to_py_err(err: cf_types::CfError) -> PyErr {
    // The Python side gets the RFC 9457 problem detail, not a flattened
    // message, so an error crossing the boundary keeps its code.
    PyValueError::new_err(err.to_problem_details().to_string())
}

/// Version of the compiled core.
mod pipeline;

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
#[pyo3(signature = (document_json, omit_field=None, include_volatile=false))]
fn digest_json(
    document_json: &str,
    omit_field: Option<&str>,
    include_volatile: bool,
) -> PyResult<String> {
    let value: serde_json::Value =
        serde_json::from_str(document_json).map_err(|e| PyValueError::new_err(e.to_string()))?;

    let mut omit: Vec<&str> = omit_field.into_iter().collect();
    if !include_volatile {
        omit.extend(cf_types::VOLATILE_DIGEST_FIELDS);
    }

    let digest = if omit.is_empty() {
        Digest::of_json(&value)
    } else {
        Digest::of_document_omitting(&value, &omit)
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

/// The legacy incremental prefix predicate.
///
/// `P1-27b`. Exposed so `python/tests/test_p1_27b_baseline_fidelity.py` can
/// compare it against the block inside `codon_6.0.py`'s solver loop, prefix by
/// prefix. That comparison is where a transcription error would actually
/// live — the two windowed checks and their `+ 2` slack are easy to paraphrase
/// into something that agrees on most inputs.
#[pyfunction]
fn legacy_prefix_accepts(
    prefix: &str,
    motifs: Vec<String>,
    k: usize,
    max_allowed_repeats: u32,
    backbone_counts: std::collections::HashMap<String, u32>,
) -> bool {
    let scenario = LegacyScenario {
        motifs: motifs.into_iter().collect(),
        k,
        max_allowed_repeats,
        backbone_counts,
    };
    legacy::prefix_is_acceptable(prefix.as_bytes(), &scenario)
}

/// Run the ported backtracker over a coding sequence.
///
/// Returns the solved coding sequence, or `None` when every restart exhausted
/// its backtrack budget — which is a timeout, not a proof that no solution
/// exists.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
fn baseline_backtracking(
    cds: &str,
    translation_table: u8,
    motifs: Vec<String>,
    k: usize,
    max_allowed_repeats: u32,
    backbone_counts: std::collections::HashMap<String, u32>,
    seed: u64,
) -> PyResult<Option<String>> {
    let sequence = CodonSequence::from_cds(cds, translation_table).map_err(to_py_err)?;
    let genotype = Genotype::new(sequence);
    let scenario = LegacyScenario {
        motifs: motifs.into_iter().collect(),
        k,
        max_allowed_repeats,
        backbone_counts,
    };
    let mut rng = SearchRng::new(seed);
    let report = cf_search::baseline_backtracking(
        &genotype,
        &scenario,
        &CodonWeights::uniform(),
        // The prototype's own constants.
        BaselineBudget::from_prototype(1_000, 2_000),
        &mut rng,
    )
    .map_err(to_py_err)?;
    Ok(report
        .solution
        .map(|g| g.codons().iter().map(|c| c.to_string_upper()).collect()))
}

/// Advance a search by a whole number of generations.
///
/// `P1-29`. The single crossing the job lifecycle needs: Python owns the
/// budget, the checkpoint store and the cancellation flag; Rust owns the
/// search. ADR-0001 caps crossings at one per generation, so a segment of `n`
/// generations is comfortably inside the limit, and the state returned *is* the
/// checkpoint — there is no second representation to keep in sync.
///
/// `checkpoint_json` is `null` to start a fresh run from `cds`.
#[pyfunction]
#[pyo3(signature = (
    checkpoint_json,
    cds,
    translation_table,
    five_prime_flank,
    three_prime_flank,
    constraint_set_json,
    codon_weights,
    tabu_leak,
    root_seed,
    generations,
))]
#[allow(clippy::too_many_arguments)]
fn search_run_segment(
    checkpoint_json: Option<&str>,
    cds: &str,
    translation_table: u8,
    five_prime_flank: &str,
    three_prime_flank: &str,
    constraint_set_json: &str,
    codon_weights: Vec<f64>,
    tabu_leak: f64,
    root_seed: u64,
    generations: u32,
) -> PyResult<String> {
    let sequence = CodonSequence::from_cds(cds, translation_table)
        .and_then(|s| s.with_five_prime_flank(five_prime_flank))
        .and_then(|s| s.with_three_prime_flank(three_prime_flank))
        .map_err(to_py_err)?;
    let start = Genotype::new(sequence);

    let set: serde_json::Value = serde_json::from_str(constraint_set_json)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let checker = ConstraintChecker::compile(&ConstraintSet::from_json(&set).map_err(to_py_err)?)
        .map_err(to_py_err)?;

    let mut table = [0.0f64; 64];
    for (slot, value) in table.iter_mut().zip(codon_weights.iter()) {
        *slot = *value;
    }
    let weights = if codon_weights.is_empty() {
        CodonWeights::uniform()
    } else {
        CodonWeights::from_relative_adaptiveness(table).map_err(to_py_err)?
    };

    let checkpoint = match checkpoint_json {
        None | Some("null") => JobCheckpoint::start(&start),
        Some(text) => {
            let value: serde_json::Value =
                serde_json::from_str(text).map_err(|e| PyValueError::new_err(e.to_string()))?;
            JobCheckpoint::from_json(&value).map_err(to_py_err)?
        }
    };

    let mix = OperatorMix::without_context_model();
    let report = cf_search::run_segment(
        &checkpoint,
        &SegmentRequest {
            checker: &checker,
            mix: &mix,
            weights: &weights,
            leak: tabu_leak,
            region: 0,
            masking: KmerMasking::default(),
            root_seed,
        },
        generations,
    )
    .map_err(to_py_err)?;

    Ok(serde_json::json!({
        "checkpoint": report.checkpoint.to_json().map_err(to_py_err)?,
        "generations_run": report.generations_run,
        "warnings": report.warnings.iter().map(|w| {
            let mut problem = w.to_problem_details();
            problem["stage"] = serde_json::json!("operators");
            problem["subject"] = serde_json::json!(operator_subject(&w.detail));
            problem
        }).collect::<Vec<_>>(),
    })
    .to_string())
}

/// The operator a `CF-SRCH-001` warning is about.
///
/// The detail names it first, by construction in `UnavailableOperator::warning`.
/// Extracting it here rather than restructuring the error type keeps the
/// `stage`/`subject` contract in one place — the boundary that emits the
/// document — instead of spreading provenance concerns through `cf-search`.
pub(crate) fn operator_subject(detail: &str) -> String {
    detail
        .strip_prefix("mutation operator ")
        .and_then(|rest| rest.split_whitespace().next())
        .unwrap_or("")
        .to_string()
}

/// Run tier 1 seeding. `P1-30`.
#[pyfunction]
#[pyo3(signature = (request_json, beam, seeds))]
fn search_tier1_seed(request_json: &str, beam: usize, seeds: usize) -> PyResult<String> {
    with_tier_inputs(request_json, |inputs| {
        pipeline::run_tier1(inputs, beam, seeds)
    })
}

/// Run tier 2 block coordinate ascent. `P1-30`.
#[pyfunction]
#[pyo3(signature = (request_json, block, sweeps))]
fn search_tier2_polish(request_json: &str, block: usize, sweeps: u32) -> PyResult<String> {
    with_tier_inputs(request_json, |inputs| {
        pipeline::run_tier2(inputs, block, sweeps)
    })
}

/// The optimization direction each requested term declares. `P1-30`.
///
/// Exposed so the job layer can populate `resolved_objectives[].direction` and
/// `pareto_front.objective_sense` from the providers rather than from a
/// default, and so that a `minimize` term is never silently reported as
/// maximized.
#[pyfunction]
fn term_directions(request_json: &str) -> PyResult<String> {
    with_tier_inputs(request_json, pipeline::declared_directions)
}

/// The role each requested term declares. `P2-11`.
///
/// Mirrors `term_directions`. Without it the Python layer would need its own
/// copy of the term-to-role mapping, which is the duplicated-rule pattern this
/// project has now found five times.
#[pyfunction]
fn term_roles(request_json: &str) -> PyResult<String> {
    with_tier_inputs(request_json, pipeline::declared_roles)
}

/// The capability each requested term declares. `P2-11`.
#[pyfunction]
fn term_capabilities(request_json: &str) -> PyResult<String> {
    with_tier_inputs(request_json, pipeline::declared_capabilities)
}

/// Score a sequence on each requested term. `P1-30`.
#[pyfunction]
fn term_scores(request_json: &str, cds: &str) -> PyResult<String> {
    with_tier_inputs(request_json, |inputs| pipeline::score_terms(inputs, cds))
}

/// Run tier 3 for a segment. `P2-04`.
///
/// The job's execution path for the population loop, replacing the
/// mutation-only driver. Whole generations per crossing, so ADR-0001's
/// one-crossing-per-generation cap holds, and the state returned is the
/// checkpoint.
#[pyfunction]
#[pyo3(signature = (
    request_json,
    checkpoint_json,
    population,
    generations,
    stagnation_window,
    l0_scorer=None,
    l0_axis=0,
    l0_budget=0,
))]
// The signature *is* the Python API. Grouping these into a struct would move
// the argument list into a dict the caller has to build, which is worse at the
// call site and no clearer here.
#[allow(clippy::too_many_arguments)]
fn search_tier3_segment(
    request_json: &str,
    checkpoint_json: Option<&str>,
    population: usize,
    generations: u32,
    stagnation_window: usize,
    l0_scorer: Option<Bound<'_, PyAny>>,
    l0_axis: usize,
    l0_budget: usize,
) -> PyResult<String> {
    let checkpoint: Option<serde_json::Value> = match checkpoint_json {
        None | Some("null") => None,
        Some(text) => {
            Some(serde_json::from_str(text).map_err(|e| PyValueError::new_err(e.to_string()))?)
        }
    };
    let scorer = l0_scorer.map(|callable| pipeline::PyL0Scorer::new(callable, l0_axis, l0_budget));
    // Drop a stash left by an earlier call on this thread before this one runs.
    pipeline::clear_escaped_scorer_error();
    let result = with_tier_inputs(request_json, |inputs| {
        pipeline::run_tier3_segment(
            inputs,
            checkpoint.as_ref(),
            population,
            generations,
            stagnation_window,
            scorer
                .as_ref()
                .map(|s| s as &dyn cf_search::tier3::L0Scorer),
        )
    });
    // `DEP-ORCH-004`. `with_tier_inputs` already ran the scorer's `CfError`
    // through `to_py_err`, which publishes every `CfError` as a problem
    // document and does not consult the stash. When the stash is set, that
    // document is the sentinel (`CF-LAND-012` with no problem behind it).
    //
    // Any caller on this path MUST call `take_escaped_scorer_error` before
    // treating that error as real. Skipping it returns the sentinel, and
    // `Job.run` records a forged `CF-LAND-012`. This is the only production
    // caller, and the only recovery: do not add another inside `to_py_err`.
    if let Some(escaped) = pipeline::take_escaped_scorer_error() {
        return Err(escaped);
    }
    result
}

/// The RNG root tier 3 and the generation loop derive from. `DEP-DET-003`.
///
/// Exposed so the job layer passes a namespaced root rather than the job's own
/// seed, which is what tiers 1 and 2 already receive.
#[pyfunction]
fn tier3_root_seed(root_seed: u64) -> u64 {
    pipeline::tier3_seed(root_seed)
}

/// Exact hypervolume of a minimized front. `P2-01`, exposed for `P2-15`/`P2-19`.
///
/// The reference point is the caller's. This function does not invent one:
/// the search loop's own nadir-plus-one is a within-run stagnation signal and
/// is not comparable across runs.
#[pyfunction]
fn indicator_hypervolume(points: Vec<Vec<f64>>, reference: Vec<f64>) -> PyResult<f64> {
    cf_metrics::hypervolume(&points, &reference)
        .map_err(|err| PyValueError::new_err(err.to_string()))
}

/// Das–Dennis directions, in the crate's fixed order. `P2-01`.
#[pyfunction]
fn indicator_das_dennis(objectives: usize, divisions: usize) -> PyResult<Vec<Vec<f64>>> {
    cf_metrics::das_dennis(objectives, divisions)
        .map_err(|err| PyValueError::new_err(err.to_string()))
}

/// The scalarized objective the tiers maximize, for one sequence. `P1-30`.
#[pyfunction]
fn scalarized_score(request_json: &str, cds: &str) -> PyResult<f64> {
    let text = with_tier_inputs(request_json, |inputs| {
        Ok(serde_json::json!(pipeline::scalarized_score(inputs, cds)?))
    })?;
    text.parse::<f64>()
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// The same scalarization, plus which axes entered it. `TST-CLI-012`.
///
/// One loop produces the score and the two axis lists. The front's
/// `recommended.informed_by` / `not_informed_by` are this document, not a
/// second pass over capabilities.
#[pyfunction]
fn scalarized_score_report(request_json: &str, cds: &str) -> PyResult<String> {
    with_tier_inputs(request_json, |inputs| {
        pipeline::scalarized_score_report(inputs, cds)
    })
}

/// Decode the shared tier request document and hand it to a runner.
fn with_tier_inputs(
    request_json: &str,
    run: impl FnOnce(&pipeline::TierInputs<'_>) -> cf_types::Result<serde_json::Value>,
) -> PyResult<String> {
    let request: serde_json::Value =
        serde_json::from_str(request_json).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let field = |name: &str| -> PyResult<&serde_json::Value> {
        request
            .get(name)
            .ok_or_else(|| PyValueError::new_err(format!("tier request is missing {name}")))
    };
    let text = |name: &str| -> PyResult<&str> {
        field(name)?
            .as_str()
            .ok_or_else(|| PyValueError::new_err(format!("{name} must be a string")))
    };

    let (cds, five, three) = (
        text("cds")?,
        text("five_prime_flank")?,
        text("three_prime_flank")?,
    );
    let sequence = CodonSequence::from_cds(
        cds,
        field("translation_table")?.as_u64().unwrap_or(11) as u8,
    )
    .and_then(|s| s.with_five_prime_flank(five))
    .and_then(|s| s.with_three_prime_flank(three))
    .map_err(to_py_err)?;

    let constraint_set = ConstraintSet::from_json(field("constraint_set")?).map_err(to_py_err)?;
    let checker = ConstraintChecker::compile(&constraint_set).map_err(to_py_err)?;

    let terms: Vec<String> = field("terms")?
        .as_array()
        .ok_or_else(|| PyValueError::new_err("terms must be an array"))?
        .iter()
        .filter_map(|v| v.as_str().map(str::to_string))
        .collect();
    let weights: Vec<f64> = field("weights")?
        .as_array()
        .ok_or_else(|| PyValueError::new_err("weights must be an array"))?
        .iter()
        .filter_map(serde_json::Value::as_f64)
        .collect();

    let value = run(&pipeline::TierInputs {
        sequence: &sequence,
        checker: &checker,
        constraint_set: &constraint_set,
        host_profile: field("host_profile")?,
        terms: &terms,
        weights: &weights,
        host_id: text("host_id")?,
        host_digest: text("host_digest")?,
        root_seed: field("root_seed")?.as_u64().unwrap_or(0),
    })
    .map_err(to_py_err)?;
    Ok(value.to_string())
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
    m.add_function(wrap_pyfunction!(legacy_prefix_accepts, m)?)?;
    m.add_function(wrap_pyfunction!(baseline_backtracking, m)?)?;
    m.add_function(wrap_pyfunction!(search_run_segment, m)?)?;
    m.add_function(wrap_pyfunction!(search_tier1_seed, m)?)?;
    m.add_function(wrap_pyfunction!(search_tier2_polish, m)?)?;
    m.add_function(wrap_pyfunction!(term_directions, m)?)?;
    m.add_function(wrap_pyfunction!(term_roles, m)?)?;
    m.add_function(wrap_pyfunction!(term_capabilities, m)?)?;
    m.add_function(wrap_pyfunction!(term_scores, m)?)?;
    m.add_function(wrap_pyfunction!(indicator_hypervolume, m)?)?;
    m.add_function(wrap_pyfunction!(indicator_das_dennis, m)?)?;
    m.add_function(wrap_pyfunction!(scalarized_score, m)?)?;
    m.add_function(wrap_pyfunction!(scalarized_score_report, m)?)?;
    m.add_function(wrap_pyfunction!(tier3_root_seed, m)?)?;
    m.add_function(wrap_pyfunction!(search_tier3_segment, m)?)?;
    m.add_function(wrap_pyfunction!(fit_potts_context, m)?)?;
    m.add_function(wrap_pyfunction!(score_potts_context, m)?)?;
    m.add_function(wrap_pyfunction!(score_codon_usage_baseline, m)?)?;
    Ok(())
}

fn json_strings(value: &serde_json::Value, key: &str) -> PyResult<Vec<String>> {
    value
        .get(key)
        .and_then(|item| item.as_array())
        .ok_or_else(|| PyValueError::new_err(format!("{key} must be a list of sequences")))?
        .iter()
        .map(|item| {
            item.as_str()
                .map(str::to_string)
                .ok_or_else(|| PyValueError::new_err(format!("{key} contains a non-string")))
        })
        .collect()
}

/// Fit the synonymous Potts model. One crossing for the whole fit.
///
/// The dwell head is not on this surface. A request for it is not answered
/// here with zeros.
#[pyfunction]
fn fit_potts_context(request_json: &str) -> PyResult<String> {
    let request: serde_json::Value =
        serde_json::from_str(request_json).map_err(|err| PyValueError::new_err(err.to_string()))?;
    let sequences = json_strings(&request, "sequences")?;
    let window = request
        .get("window")
        .and_then(|item| item.as_u64())
        .ok_or_else(|| PyValueError::new_err("window must be an integer"))? as usize;
    let lambda_h = request
        .get("lambda_h")
        .and_then(|item| item.as_f64())
        .unwrap_or(1e-3);
    let lambda_j = request
        .get("lambda_j")
        .and_then(|item| item.as_f64())
        .unwrap_or(1e-3);
    let max_iterations = request
        .get("max_iterations")
        .and_then(|item| item.as_u64())
        .unwrap_or(25) as u32;
    let report = cf_epistasis::fit_potts(&sequences, window, lambda_h, lambda_j, max_iterations)
        .map_err(PyValueError::new_err)?;
    Ok(serde_json::json!({
        "window": report.model.window,
        "h": report.model.h,
        "couplings": report.model.couplings,
        "nats_per_codon": report.nats_per_codon,
        "n_positions": report.n_positions,
        "iterations": report.iterations,
        "gradient_norm": report.gradient_norm,
        "converged": report.converged,
        "lambda_h": report.lambda_h,
        "lambda_j": report.lambda_j,
    })
    .to_string())
}

#[pyfunction]
fn score_potts_context(request_json: &str) -> PyResult<String> {
    let request: serde_json::Value =
        serde_json::from_str(request_json).map_err(|err| PyValueError::new_err(err.to_string()))?;
    let sequences = json_strings(&request, "sequences")?;
    let window = request
        .get("window")
        .and_then(|item| item.as_u64())
        .ok_or_else(|| PyValueError::new_err("window must be an integer"))? as usize;
    let h = request
        .get("h")
        .and_then(|item| item.as_array())
        .ok_or_else(|| PyValueError::new_err("h must be an array"))?
        .iter()
        .map(|item| {
            item.as_f64()
                .ok_or_else(|| PyValueError::new_err("h contains a non-number"))
        })
        .collect::<PyResult<Vec<f64>>>()?;
    let couplings = request
        .get("couplings")
        .and_then(|item| item.as_array())
        .ok_or_else(|| PyValueError::new_err("couplings must be an array"))?
        .iter()
        .map(|item| {
            item.as_f64()
                .ok_or_else(|| PyValueError::new_err("couplings contains a non-number"))
        })
        .collect::<PyResult<Vec<f64>>>()?;
    let model = cf_epistasis::PottsModel {
        window,
        h,
        couplings,
    };
    let (nats, positions) =
        cf_epistasis::score_naturalness(&sequences, &model).map_err(PyValueError::new_err)?;
    Ok(serde_json::json!({
        "nats_per_codon": nats,
        "n_positions": positions,
    })
    .to_string())
}

#[pyfunction]
fn score_codon_usage_baseline(request_json: &str) -> PyResult<String> {
    let request: serde_json::Value =
        serde_json::from_str(request_json).map_err(|err| PyValueError::new_err(err.to_string()))?;
    let train = json_strings(&request, "train")?;
    let held_out = json_strings(&request, "held_out")?;
    let baseline =
        cf_epistasis::score_codon_usage(&train, &held_out).map_err(PyValueError::new_err)?;
    Ok(serde_json::json!({
        "nats_per_codon": baseline.nats_per_codon,
        "n_scored": baseline.n_scored,
        "n_unseen": baseline.n_unseen,
    })
    .to_string())
}
