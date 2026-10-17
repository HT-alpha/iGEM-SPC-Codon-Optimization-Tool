//! Tier assembly: the registry, tier 1 and tier 2, behind one crossing each.
//!
//! Ticket P1-30, Rust half.
//!
//! ## Why this lives in `cf-pyext` and not in `cf-search`
//!
//! Running a tier needs a `TermRegistry` holding concrete providers, and the
//! concrete providers are `cf-landscape`'s. `cf-search` depends on
//! `cf-landscape-api` and deliberately not on `cf-landscape` — it never names a
//! concrete predictor, which is what keeps the search engine testable against
//! mock terms and stops a term implementation from becoming a search
//! dependency. Something has to know both, and `cf-pyext` already does: it is
//! the assembly point by construction, being the one crate the job layer calls.
//!
//! ## One crossing per tier
//!
//! ADR-0001 caps crossings at one per generation. Tier 1 runs once and tier 2
//! runs a bounded number of sweeps, so each is a single call. The job layer
//! sequences them; neither reaches back into Python.

use cf_constraints::{explain_deviations, resolve_feasible_seed, ConstraintChecker};
use cf_landscape::terms::constraints::ViolationCountTerm;
use cf_landscape::terms::structure::LocalAccessibilityTerm;
use cf_landscape::terms::usage::{CaiTerm, CodonUsageTable, GcContentTerm};
use cf_landscape::ConstructLayout;
use cf_landscape_api::provider::EvalContext;
use cf_landscape_api::provider::TermRegistry;
use cf_landscape_api::query::Budget;
use cf_landscape_api::Capability;
use cf_landscape_api::HostContext;
use cf_search::tier3::{
    seed_population, Individual, L0Scorer, ObjectiveEvaluator, Selection, Tier3Checkpoint,
    Tier3Config, Tier3Request,
};
use cf_search::{
    tier1_seed, tier2_block_ca, CodonWeights, Genotype, OperatorMix, SearchRng, Tier1Config,
    Tier1Request, Tier2Config, Tier2Plan, Tier2Request,
};
use cf_types::error::lookup;
use cf_types::Codon;
use cf_types::{code, CfError, CodonSequence, Digest, HostId, Result, TermId};
use pyo3::types::PyAnyMethods;
use serde_json::{json, Value};
use std::cell::RefCell;

/// The folding span the structure term declares, in nucleotides.
///
/// Must match `codonforge.orchestrator.structure_diag.SPAN_NT`. Two constants
/// for one quantity is the pattern this project keeps finding, so
/// `TST-STRU-007` asserts they agree rather than trusting a comment — the Rust
/// side needs it to declare an interaction radius before any Python has run,
/// and the Python side needs it to call the engine.
pub const STRUCTURE_SPAN_NT: usize = 40;

/// Build a registry holding exactly the terms a job asked for.
///
/// Refuses an unknown term rather than skipping it. A job that asked for four
/// objectives and silently optimized three is the failure this whole warning
/// apparatus exists to make impossible, and it would be perverse to introduce
/// it here.
pub fn build_registry(
    terms: &[String],
    host_profile: &Value,
    sequence: &CodonSequence,
    constraint_set: &cf_constraints::spec::ConstraintSet,
) -> Result<TermRegistry> {
    let mut registry = TermRegistry::new();
    for term in terms {
        match term.as_str() {
            "usage.cai" => {
                let table =
                    CodonUsageTable::from_host_profile(host_profile, sequence.table().id())?;
                // `register_default` uses `max_role`. A lower role here would
                // make Control A/2 (test plan §1a) newly exercisable.
                registry.register_default(Box::new(CaiTerm::new(table)?))?;
            }
            "usage.gc_content" => {
                registry.register_default(Box::new(GcContentTerm::new()))?;
            }
            "structure.local_accessibility" => {
                // Contract only. The values arrive through the batch path
                // (P2-05a); this registration is what makes the term's
                // capability, role, direction and radius readable from the one
                // place every consumer already looks.
                registry
                    .register_default(Box::new(LocalAccessibilityTerm::new(STRUCTURE_SPAN_NT)?))?;
            }
            "constraint.violation_count" => {
                // Compiled afresh rather than shared: `ConstraintChecker` is not
                // `Clone`, and the term owns its checker for the life of the
                // registry. Compiling the same set twice is deterministic, so
                // this costs a little work and no correctness.
                registry.register_default(Box::new(ViolationCountTerm::new(
                    ConstraintChecker::compile(constraint_set)?,
                    ConstructLayout::from_sequence(sequence),
                )?))?;
            }
            other => {
                return Err(CfError::new(
                    code::LAND_TERM_UNKNOWN,
                    format!(
                        "no built-in provider supplies {other}; this build ships \
                         usage.cai, usage.gc_content, constraint.violation_count and \
                         structure.local_accessibility, and a job is refused rather \
                         than run against fewer terms than it asked for"
                    ),
                ));
            }
        }
    }
    Ok(registry)
}

/// The optimization direction each term declares, asked of the provider.
///
/// Asked, not tabulated. A table here would be a second declaration of
/// something the provider already states, and the two would drift — with the
/// symptom being a `pareto_front.objective_sense` that disagrees with the
/// search that produced it, which is unfalsifiable from the document alone.
///
/// Note what this is *not* for. Tier 1 and tier 2 both apply
/// `provider.direction()` themselves as a sense multiplier on the weighted
/// sum, so a caller's weights are magnitudes and must not be signed by the
/// caller as well. Applying the sign twice inverts the objective while every
/// document still reads correctly, which is the precise failure this direction
/// is carried to prevent.
pub fn declared_directions(inputs: &TierInputs<'_>) -> Result<Value> {
    let registry = build_registry(
        inputs.terms,
        inputs.host_profile,
        inputs.sequence,
        inputs.constraint_set,
    )?;
    let mut out = serde_json::Map::new();
    for term in inputs.terms {
        let id = TermId::new(term)?;
        let provider = registry.get(&id).ok_or_else(|| {
            CfError::new(code::LAND_TERM_UNKNOWN, format!("{term} is not registered"))
        })?;
        out.insert(
            term.clone(),
            json!(match provider.provider.direction() {
                cf_types::Direction::Maximize => "maximize",
                cf_types::Direction::Minimize => "minimize",
            }),
        );
    }
    Ok(Value::Object(out))
}

/// The role each requested term declares. `P2-11`.
///
/// Sourced from `Registration.role` — the **effective** role the registry holds,
/// not `TermProvider::max_role`, which is a ceiling. The two differ: a provider
/// declaring `Diagnostic` as its maximum is registered as diagnostic, while one
/// permitting `ParetoAxis` may still be registered lower. Reporting the ceiling
/// would tell a caller what a term *could* be rather than what it *was*, and the
/// first fixture written against this caught exactly that.
///
/// The distinction matters most for the case that motivated exposing it: a
/// `diagnostic` term is computed and reported but never ranked on, and a caller
/// could previously only infer that from its *absence* from `objective_order`.
/// An inference is not a contract — it holds until someone adds a second reason
/// for a term to be missing.
pub fn declared_roles(inputs: &TierInputs<'_>) -> Result<Value> {
    let registry = build_registry(
        inputs.terms,
        inputs.host_profile,
        inputs.sequence,
        inputs.constraint_set,
    )?;
    let mut out = serde_json::Map::new();
    for term in inputs.terms {
        let id = TermId::new(term)?;
        let provider = registry.get(&id).ok_or_else(|| {
            CfError::new(code::LAND_TERM_UNKNOWN, format!("{term} is not registered"))
        })?;
        out.insert(term.clone(), json!(provider.role.as_str()));
    }
    Ok(Value::Object(out))
}

/// The capability each requested term declares. `P2-11`.
///
/// `TierPlan.plan` says in as many words that capabilities "come from the
/// registry rather than from configuration, because a capability a caller
/// asserts is not a capability" — and then `design()` asserted `L1` for every
/// term. That went unnoticed while every built-in was in fact `L1`. The first
/// `L0` term made it visible immediately: the plan routed an out-of-core term
/// into tier 1, which asked it for an in-core value it must not supply.
pub fn declared_capabilities(inputs: &TierInputs<'_>) -> Result<Value> {
    let registry = build_registry(
        inputs.terms,
        inputs.host_profile,
        inputs.sequence,
        inputs.constraint_set,
    )?;
    let mut out = serde_json::Map::new();
    for term in inputs.terms {
        let id = TermId::new(term)?;
        let provider = registry.get(&id).ok_or_else(|| {
            CfError::new(code::LAND_TERM_UNKNOWN, format!("{term} is not registered"))
        })?;
        out.insert(
            term.clone(),
            json!(match provider.provider.capability() {
                Capability::L0 => "L0",
                Capability::L1 => "L1",
                Capability::L2 => "L2",
            }),
        );
    }
    Ok(Value::Object(out))
}

/// Score a sequence on each requested term.
///
/// The front's `objective_vector` has to be the value the providers actually
/// produce for the sequence being reported. Recomputing it here from the same
/// registry the tiers used is the only way to say that honestly; carrying a
/// number forward from a search internal would report the surrogate the search
/// optimized rather than the objective the document claims.
pub fn score_terms(inputs: &TierInputs<'_>, cds: &str) -> Result<Value> {
    // Rebuilt from the reported cds rather than mutated in place: the protein,
    // the flanks and the lock set all have to come from the same source the
    // job ran under, and a partial copy is how those drift apart.
    let sequence = CodonSequence::from_cds(cds, inputs.sequence.table().id())?
        .with_five_prime_flank(
            std::str::from_utf8(inputs.sequence.five_prime_flank()).unwrap_or(""),
        )?
        .with_three_prime_flank(
            std::str::from_utf8(inputs.sequence.three_prime_flank()).unwrap_or(""),
        )?;
    let registry = build_registry(
        inputs.terms,
        inputs.host_profile,
        &sequence,
        inputs.constraint_set,
    )?;
    let host = host_context(inputs)?;
    let ctx = EvalContext::new(&sequence, &host, Budget::default(), 0);

    let mut out = serde_json::Map::new();
    for term in inputs.terms {
        let id = TermId::new(term)?;
        let registration = registry.get(&id).ok_or_else(|| {
            CfError::new(code::LAND_TERM_UNKNOWN, format!("{term} is not registered"))
        })?;
        // An out-of-core term has no in-core value, and `null` is the honest
        // report: the schema's `unknown_terms` says the same thing at the front
        // level. Calling `score` here would hit `CF-STRU-006`, which is a
        // refusal aimed at a caller reaching for a value it must not have -
        // not at the ordinary path asking what can be known in core.
        if registration.provider.capability() == Capability::L0 {
            out.insert(term.clone(), Value::Null);
            continue;
        }
        out.insert(term.clone(), json!(registration.provider.score(&ctx)?));
    }
    Ok(Value::Object(out))
}

/// One scalarization: the number `best()` ranks on, and the axes that made it.
///
/// The two lists are not a second computation. They are a byproduct of the
/// same loop that produces the score: a term that contributed a number is
/// `informed_by`, a term that had none (an L0 / out-of-core term, today) is
/// `not_informed_by`. A caller that built those lists from capabilities
/// instead would drift the moment the skip condition here changed.
#[derive(Debug)]
pub struct Scalarization {
    pub score: f64,
    pub informed_by: Vec<String>,
    pub not_informed_by: Vec<String>,
}

/// The scalarized objective tier 1 and tier 2 maximize, for one sequence.
///
/// The same formula the tiers use, including reading `provider.direction()` for
/// the sense — so a caller comparing two sequences with this is comparing them
/// the way the search did, rather than by a second convention that could
/// disagree with it.
pub fn scalarize(inputs: &TierInputs<'_>, cds: &str) -> Result<Scalarization> {
    let scores = score_terms(inputs, cds)?;
    let directions = declared_directions(inputs)?;
    let mut total = 0.0;
    let mut informed_by = Vec::new();
    let mut not_informed_by = Vec::new();
    for (index, term) in inputs.terms.iter().enumerate() {
        // Out-of-core terms are absent from the scalar, and that is a real
        // limitation rather than an oversight: this formula is what tiers 1
        // and 2 maximize, and those tiers evaluate in core. Including a term
        // they never saw would make the archive's ordering disagree with the
        // search that produced the candidates. The full objective vector,
        // which does carry the out-of-core axis, is what tier 3 ranks on.
        let Some(value) = scores.get(term).and_then(Value::as_f64) else {
            not_informed_by.push(term.clone());
            continue;
        };
        let sense = match directions.get(term).and_then(Value::as_str) {
            Some("maximize") => 1.0,
            _ => -1.0,
        };
        total += sense * inputs.weights.get(index).copied().unwrap_or(1.0) * value;
        informed_by.push(term.clone());
    }
    Ok(Scalarization {
        score: total,
        informed_by,
        not_informed_by,
    })
}

pub fn scalarized_score(inputs: &TierInputs<'_>, cds: &str) -> Result<f64> {
    Ok(scalarize(inputs, cds)?.score)
}

pub fn scalarized_score_report(inputs: &TierInputs<'_>, cds: &str) -> Result<Value> {
    let scalarization = scalarize(inputs, cds)?;
    Ok(json!({
        "score": scalarization.score,
        "informed_by": scalarization.informed_by,
        "not_informed_by": scalarization.not_informed_by,
    }))
}

pub struct TierInputs<'a> {
    pub sequence: &'a CodonSequence,
    pub checker: &'a ConstraintChecker,
    pub constraint_set: &'a cf_constraints::spec::ConstraintSet,
    pub host_profile: &'a Value,
    pub terms: &'a [String],
    pub weights: &'a [f64],
    pub host_id: &'a str,
    pub host_digest: &'a str,
    pub root_seed: u64,
}

fn host_context(inputs: &TierInputs<'_>) -> Result<HostContext> {
    Ok(HostContext::new(
        HostId::new(inputs.host_id)?,
        Digest::parse(inputs.host_digest)?,
    ))
}

fn term_ids(terms: &[String]) -> Result<Vec<TermId>> {
    terms.iter().map(TermId::new).collect()
}

/// Run tier 1 and return the seeds it produced.
pub fn run_tier1(inputs: &TierInputs<'_>, beam: usize, seeds: usize) -> Result<Value> {
    let base = Genotype::new(inputs.sequence.clone());
    let registry = build_registry(
        inputs.terms,
        inputs.host_profile,
        inputs.sequence,
        inputs.constraint_set,
    )?;
    let host = host_context(inputs)?;
    let config = Tier1Config::default().with_beam(beam).with_seeds(seeds);
    let mut rng = SearchRng::new(stage_seed(inputs.root_seed, STREAM_TIER1));

    let report = tier1_seed(
        &Tier1Request {
            base: &base,
            registry: &registry,
            terms: &term_ids(inputs.terms)?,
            weights: inputs.weights,
            host: &host,
            checker: inputs.checker,
            config: &config,
        },
        &mut rng,
    )?;

    Ok(json!({
        "seeds": report.seeds.iter().map(|g| g.sequence().to_json()).collect::<Result<Vec<_>>>()?,
        "surrogate": report.surrogate,
        "states_pruned": report.states_pruned,
        "post_check_rejected": report.post_check_rejected,
        "delta_tables": report.delta_tables,
        "rebases": report.rebases,
    }))
}

/// Polish a seed with tier 2.
pub fn run_tier2(inputs: &TierInputs<'_>, block: usize, sweeps: u32) -> Result<Value> {
    let seed = Genotype::new(inputs.sequence.clone());
    let registry = build_registry(
        inputs.terms,
        inputs.host_profile,
        inputs.sequence,
        inputs.constraint_set,
    )?;
    let host = host_context(inputs)?;
    let mut config = Tier2Config::from_planning(block);
    config.max_sweeps = sweeps;

    let plan = Tier2Plan::check(&registry, &term_ids(inputs.terms)?, inputs.weights, &config)?;

    let report = tier2_block_ca(
        &Tier2Request {
            seed: &seed,
            registry: &registry,
            host: &host,
            checker: inputs.checker,
            config: &config,
        },
        &plan,
        stage_seed(inputs.root_seed, STREAM_TIER2),
    )?;

    Ok(json!({
        "incumbent": report.incumbent.sequence().to_json()?,
        "objective": report.objective,
        "sweeps": report.sweeps,
        "blocks_committed": report.blocks_committed,
        "blocks_rejected": report.blocks_rejected,
        "delta_tables": report.delta_tables,
        "warnings": plan.warnings().iter().chain(report.warnings.iter()).map(|w| {
            let mut problem = w.to_problem_details();
            problem["stage"] = json!("tier2");
            problem["subject"] = json!(w.detail.split_whitespace().next().unwrap_or(""));
            problem
        }).collect::<Vec<_>>(),
        "notes": report.notes,
    }))
}

/// One synonymous seed after tier 1 refuses a constraint the DP cannot carry.
///
/// `DEP-SRCH-015`. The report is always returned, feasible or not. The job
/// decides whether to continue. A feasible seed carries `CF-SRCH-017` so the
/// page can say the sequence came from repair. An infeasible seed carries
/// `CF-SRCH-006` whose detail names the rules still broken and what to change.
pub fn run_seed_resolve(inputs: &TierInputs<'_>) -> Result<Value> {
    let weights = codon_adaptiveness(inputs.host_profile);
    let report = resolve_feasible_seed(inputs.sequence, inputs.checker, &weights)?;
    let mut document = json!({
        "cds": report.sequence.cds_string(),
        "sequence": report.sequence.to_json()?,
        "feasible": report.feasible,
        "cycles": report.cycles,
        "edits": report.edits,
        "residual": report.residual.iter().map(|item| json!({
            "constraint_id": item.constraint_id,
            "kind": item.kind,
            "count": item.count,
        })).collect::<Vec<_>>(),
        "unsearched": unsearched_warnings(inputs)?,
    });
    if report.feasible {
        document["warning"] = stamped(
            CfError::new(
                code::SRCH_SEED_REPAIR,
                format!(
                    "tier 1 refused a constraint the dynamic program cannot carry; \
                     seed repair committed {} synonymous edit(s) over {} cycle(s), \
                     and later tiers start from that sequence",
                    report.edits, report.cycles
                ),
            ),
            "tier1",
            "seed_repair",
        );
    } else {
        document["problem"] = stamped(
            CfError::new(code::SRCH_TIER1_INFEASIBLE, infeasible_seed_detail(&report)),
            "tier1",
            "seed_repair",
        );
    }
    Ok(document)
}

/// Rules this version does not evaluate. A warning, not a refusal: treating
/// them as satisfied would report a design the checker never looked at.
pub fn run_unsearched(inputs: &TierInputs<'_>) -> Result<Value> {
    Ok(json!({ "warnings": unsearched_warnings(inputs)? }))
}

/// Why each editable codon differs from the most-frequent synonym.
pub fn run_explain_deviations(inputs: &TierInputs<'_>, cds: &str) -> Result<Value> {
    let sequence = CodonSequence::from_cds(cds, inputs.sequence.table().id())?
        .with_five_prime_flank(
            std::str::from_utf8(inputs.sequence.five_prime_flank()).map_err(|err| {
                CfError::new(
                    code::SEQ_INVALID_BASE,
                    format!("five-prime flank is not UTF-8: {err}"),
                )
            })?,
        )?
        .with_three_prime_flank(
            std::str::from_utf8(inputs.sequence.three_prime_flank()).map_err(|err| {
                CfError::new(
                    code::SEQ_INVALID_BASE,
                    format!("three-prime flank is not UTF-8: {err}"),
                )
            })?,
        )?;
    let weights = codon_adaptiveness(inputs.host_profile);
    let deviations = explain_deviations(&sequence, inputs.checker, &weights)?;
    Ok(json!(deviations
        .into_iter()
        .map(|item| json!({
            "codon_index": item.codon_index,
            "codon": item.codon,
            "weight": item.weight,
            "baseline": item.baseline,
            "baseline_weight": item.baseline_weight,
            "constraints_if_reverted": item.constraints_if_reverted,
        }))
        .collect::<Vec<_>>()))
}

fn codon_adaptiveness(host_profile: &Value) -> [f64; 64] {
    // A missing codon stays at 1.0. Filling the gap with 0 would make every
    // unlisted codon destroy CAI, which is a property of an incomplete table
    // rather than of the sequence.
    let mut weights = [1.0; 64];
    let Some(table) = host_profile
        .pointer("/codon_usage/table")
        .and_then(Value::as_object)
    else {
        return weights;
    };
    for (label, entry) in table {
        let Ok(codon) = Codon::parse(label) else {
            continue;
        };
        if let Some(weight) = entry.get("relative_adaptiveness").and_then(Value::as_f64) {
            weights[codon.index_usize()] = weight;
        }
    }
    weights
}

fn unsearched_warnings(inputs: &TierInputs<'_>) -> Result<Vec<Value>> {
    Ok(cf_constraints::seed::unsearched_rules(inputs.checker)
        .into_iter()
        .map(|rule| {
            stamped(
                CfError::new(
                    code::FEAS_UNSEARCHED,
                    format!("{}: {}", rule.id, rule.reason),
                ),
                "constraints",
                &rule.id,
            )
        })
        .collect())
}

fn stamped(error: CfError, stage: &str, subject: &str) -> Value {
    let mut problem = error.to_problem_details();
    problem["stage"] = json!(stage);
    problem["subject"] = json!(subject);
    problem
}

fn infeasible_seed_detail(report: &cf_constraints::SeedReport) -> String {
    let broken = if report.residual.is_empty() {
        "the checker still rejects the sequence".to_string()
    } else {
        report
            .residual
            .iter()
            .map(|item| {
                format!(
                    "{} ({}) still has {} violation(s)",
                    item.constraint_id, item.kind, item.count
                )
            })
            .collect::<Vec<_>>()
            .join("; ")
    };
    format!(
        "no synonymous sequence satisfied the hard constraints after seed repair. {broken}. {}",
        repair_advice(report)
    )
}

fn repair_advice(report: &cf_constraints::SeedReport) -> String {
    let mut parts = Vec::new();
    let kinds: std::collections::BTreeSet<&str> = report
        .residual
        .iter()
        .map(|item| item.kind.as_str())
        .collect();
    if kinds.contains("gc_window") {
        parts.push("widen the GC window band or shorten the window");
    }
    if kinds.contains("kmer_repeat") {
        parts.push("raise the repeat allowance or shorten k");
    }
    if kinds.contains("restriction_site")
        || kinds.contains("forbidden_motif")
        || kinds.contains("custom_regex")
    {
        parts.push("drop a site, motif, or pattern the synonymous codons cannot avoid");
    }
    if kinds.contains("homopolymer") {
        parts.push("raise the homopolymer run limit");
    }
    if kinds.contains("codon_ban") {
        parts.push("drop a codon ban that removes every synonym");
    }
    if parts.is_empty() {
        "nothing in the residual names a rule class this repair pass knows how to describe"
            .to_string()
    } else {
        format!("what you can change: {}", parts.join("; "))
    }
}

/// Per-stage root seeds, derived from the job's root seed.
///
/// `DEP-DET-002`. Every stage previously derived its streams from the job root
/// seed directly, and their index spaces overlap: the generation loop keys on
/// `0, 1, 2, …` and tier 2's `stream_index` keys on
/// `sweep · 2B + colour · B + block`, which also starts at zero. So tier 2's
/// first block and the first generation drew *the same numbers*. That is not a
/// determinism failure — both are still reproducible — but it makes the stages
/// statistically dependent in a way nothing declared or intended, and it is the
/// kind of coupling that is found much later by someone asking why two stages
/// keep making the same choice.
///
/// Each stage now gets a root seed drawn from a distinct, reserved index, so
/// the families of streams below them cannot overlap. The tags are far from
/// zero so no plausible generation count reaches them.
fn stage_seed(root_seed: u64, tag: u64) -> u64 {
    SearchRng::stream(root_seed, tag).next_u64()
}

const STREAM_TIER1: u64 = u64::MAX - 1;
const STREAM_TIER2: u64 = u64::MAX - 2;
/// `DEP-DET-003`. The generation loop derived from the job root seed directly,
/// keying generations `0, 1, 2, …` — the same index space a population loop's
/// per-individual streams naturally use. No collision existed while tiers 1
/// and 2 were the only namespaced stages; one appears the moment tier 3
/// sub-seeds per individual, which `P2-02` does.
pub const STREAM_TIER3: u64 = u64::MAX - 3;

/// The root seed the generation loop and tier 3 share.
pub fn tier3_seed(root_seed: u64) -> u64 {
    stage_seed(root_seed, STREAM_TIER3)
}

// --------------------------------------------------------------------------- //
// Tier 3: the population loop, in the job's execution path
// --------------------------------------------------------------------------- //

/// Run tier 3 for a segment, resuming from a checkpoint document.
///
/// `P2-04`. This is what replaces the mutation loop in the job's execution
/// path. The loop is the same driver `P1-18b` established — whole generations
/// per crossing, so ADR-0001's one-crossing-per-generation cap holds — and the
/// state returned is the checkpoint, so there is no second representation.
///
/// Objectives are evaluated through the same registry tiers 1 and 2 use, with
/// each term's declared direction applied once so everything below is
/// minimization. That is deliberately not a second evaluation path: a job whose
/// tiers optimized one thing while tier 3 ranked another would be the
/// direction bug of `P1-30` in a new place.
pub fn run_tier3_segment(
    inputs: &TierInputs<'_>,
    checkpoint: Option<&Value>,
    population: usize,
    generations: u32,
    stagnation_window: usize,
    l0: Option<&dyn L0Scorer>,
) -> Result<Value> {
    let registry = build_registry(
        inputs.terms,
        inputs.host_profile,
        inputs.sequence,
        inputs.constraint_set,
    )?;
    let host = host_context(inputs)?;
    let checker = ConstraintChecker::compile(inputs.constraint_set)?;
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();

    // DEP-SRCH-012: the evaluator supplies the in-core axes and tier 3 appends
    // the out-of-core one. Splitting here rather than letting the evaluator try
    // and fail keeps the refusal in `LocalAccessibilityTerm::score` meaningful:
    // it fires when something reaches for an in-core value it must not have,
    // not on every ordinary generation.
    let in_core: Vec<String> = inputs
        .terms
        .iter()
        .filter(|term| {
            TermId::new(term.as_str())
                .ok()
                .and_then(|id| registry.get(&id))
                .map(|reg| reg.provider.capability() != Capability::L0)
                .unwrap_or(true)
        })
        .cloned()
        .collect();

    let evaluator = RegistryEvaluator {
        registry: &registry,
        host: &host,
        terms: term_ids(&in_core)?,
        senses: term_senses(&registry, &in_core)?,
    };

    let mut config = Tier3Config::from_planning(population);
    config.max_generations = u32::MAX;
    config.stagnation_window = stagnation_window;

    let request = Tier3Request {
        checker: &checker,
        mix: &mix,
        weights: &weights,
        evaluator: &evaluator,
        config: &config,
        root_seed: tier3_seed(inputs.root_seed),
        // No L0 term reaches the registry by default: the structure term ships
        // as a diagnostic (`DEP-STRU-003`) and the LM is Phase 3. `P2-05a`
        // builds the path and proves it with a real crossing; wiring a term
        // onto it permanently is `P2-11`.
        l0,
    };

    let start = match checkpoint {
        Some(document) if !document.is_null() => checkpoint_from_json(document, &request)?,
        _ => {
            let seed = Genotype::new(inputs.sequence.clone());
            seed_population(&[seed], &request)?
        }
    };
    let report = cf_search::tier3::run_generations(&start, &request, generations)?;

    Ok(json!({
        "checkpoint": checkpoint_to_json(&report.checkpoint)?,
        "front": report.front.iter().map(|m| json!({
            "cds": m.genotype.sequence().cds_string(),
            "objectives": m.objectives,
            "feasible": m.is_feasible(),
        })).collect::<Vec<_>>(),
        "termination": report.termination.as_summary(),
        "terminated": report.termination.is_final(),
        "selection": match report.selection {
            Selection::CrowdingDistance => "crowding_distance",
            Selection::ReferenceDirections => "reference_directions",
        },
        "generations_run": report.generations_run,
        "stranded": report.stranded,
        "hypervolume_history": report.checkpoint.hypervolume_history,
        "batch": report.checkpoint.batch,
        "evaluations": report.checkpoint.evaluations,
        "warnings": report.warnings.iter().map(|w| {
            let mut problem = w.to_problem_details();
            problem["stage"] = json!("operators");
            problem["subject"] = json!(crate::operator_subject(&w.detail));
            problem
        }).collect::<Vec<_>>(),
    }))
}

/// Evaluates objectives through the shared registry, direction applied once.
struct RegistryEvaluator<'a> {
    registry: &'a TermRegistry,
    host: &'a HostContext,
    terms: Vec<TermId>,
    senses: Vec<f64>,
}

impl ObjectiveEvaluator for RegistryEvaluator<'_> {
    fn evaluate(&self, genotype: &Genotype) -> Result<Vec<Option<f64>>> {
        let sequence = genotype.sequence();
        let ctx = EvalContext::new(sequence, self.host, Budget::default(), 0);
        let mut out = Vec::with_capacity(self.terms.len());
        for (term, sense) in self.terms.iter().zip(&self.senses) {
            let registration = self.registry.get(term).ok_or_else(|| {
                CfError::new(code::LAND_TERM_UNKNOWN, format!("{term} is not registered"))
            })?;
            // `sense` turns a maximized term into a minimized one. Everything
            // below `cf-metrics` minimizes, and the sign convention lives in
            // exactly one place per §3.
            out.push(Some(sense * registration.provider.score(&ctx)?));
        }
        Ok(out)
    }

    fn objective_count(&self) -> usize {
        self.terms.len()
    }
}

fn term_senses(registry: &TermRegistry, terms: &[String]) -> Result<Vec<f64>> {
    terms
        .iter()
        .map(|term| {
            let id = TermId::new(term)?;
            let registration = registry.get(&id).ok_or_else(|| {
                CfError::new(code::LAND_TERM_UNKNOWN, format!("{term} is not registered"))
            })?;
            Ok(match registration.provider.direction() {
                cf_types::Direction::Minimize => 1.0,
                cf_types::Direction::Maximize => -1.0,
            })
        })
        .collect()
}

fn checkpoint_to_json(checkpoint: &Tier3Checkpoint) -> Result<Value> {
    let mut members = Vec::with_capacity(checkpoint.members.len());
    for member in &checkpoint.members {
        let mut encoded = json!({
            "sequence": member.genotype.sequence().to_json()?,
            "objectives": member.objectives,
            "violation": member.violation,
        });
        if let Some(inherited) = &member.inherited {
            encoded["inherited"] = inherited.sequence().to_json()?;
        }
        members.push(encoded);
    }
    Ok(json!({
        "generation": checkpoint.generation,
        "members": members,
        "tabu": checkpoint.tabu.iter().map(|(region, index, codon, remaining)| json!({
            "region": region,
            "codon_index": index,
            "codon": codon.to_string_upper(),
            "remaining": remaining,
        })).collect::<Vec<_>>(),
        "evaluations": checkpoint.evaluations,
        "repairs_attempted": checkpoint.repairs_attempted,
        "repairs_succeeded": checkpoint.repairs_succeeded,
        "hypervolume_history": checkpoint.hypervolume_history,
        "batch": checkpoint.batch,
    }))
}

fn checkpoint_from_json(document: &Value, request: &Tier3Request<'_>) -> Result<Tier3Checkpoint> {
    let bad = |what: &str| {
        CfError::new(
            code::ORCH_CHECKPOINT_INVALID,
            format!("tier 3 checkpoint is missing or malformed: {what}"),
        )
    };
    let map = document.as_object().ok_or_else(|| bad("not an object"))?;
    let mut members = Vec::new();
    for entry in map
        .get("members")
        .and_then(Value::as_array)
        .ok_or_else(|| bad("members"))?
    {
        let entry = entry.as_object().ok_or_else(|| bad("member"))?;
        let sequence =
            CodonSequence::from_json(entry.get("sequence").ok_or_else(|| bad("sequence"))?)?;
        let inherited = match entry.get("inherited") {
            Some(value) if !value.is_null() => {
                Some(Genotype::new(CodonSequence::from_json(value)?))
            }
            _ => None,
        };
        members.push(Individual {
            genotype: Genotype::new(sequence),
            inherited,
            // Restored, never recomputed: an `UNKNOWN` is a fact about a
            // budget already spent, and re-running evaluation could turn it
            // into a value and change the ranking.
            objectives: entry
                .get("objectives")
                .and_then(Value::as_array)
                .ok_or_else(|| bad("objectives"))?
                .iter()
                .map(|v| v.as_f64())
                .collect(),
            violation: entry
                .get("violation")
                .and_then(Value::as_f64)
                .ok_or_else(|| bad("violation"))?,
        });
    }
    // The tabu ledger travels with everything else. Dropping it made a resumed
    // run re-propose exactly the moves repair had just reverted - the thrash
    // the ledger exists to break, restarted by the checkpoint mechanism. This
    // is the fourth thing to be found missing from a checkpoint here, and the
    // first found by the test rather than by a control.
    let mut tabu = Vec::new();
    for entry in map
        .get("tabu")
        .and_then(Value::as_array)
        .ok_or_else(|| bad("tabu"))?
    {
        let entry = entry.as_object().ok_or_else(|| bad("tabu entry"))?;
        tabu.push((
            entry
                .get("region")
                .and_then(Value::as_u64)
                .ok_or_else(|| bad("region"))? as usize,
            entry
                .get("codon_index")
                .and_then(Value::as_u64)
                .ok_or_else(|| bad("codon_index"))? as usize,
            cf_types::Codon::parse(
                entry
                    .get("codon")
                    .and_then(Value::as_str)
                    .ok_or_else(|| bad("codon"))?,
            )?,
            entry
                .get("remaining")
                .and_then(Value::as_u64)
                .ok_or_else(|| bad("remaining"))? as u32,
        ));
    }

    let _ = request;
    Ok(Tier3Checkpoint {
        generation: map
            .get("generation")
            .and_then(Value::as_u64)
            .ok_or_else(|| bad("generation"))? as u32,
        members,
        tabu,
        evaluations: map.get("evaluations").and_then(Value::as_u64).unwrap_or(0),
        repairs_attempted: map
            .get("repairs_attempted")
            .and_then(Value::as_u64)
            .unwrap_or(0),
        repairs_succeeded: map
            .get("repairs_succeeded")
            .and_then(Value::as_u64)
            .unwrap_or(0),
        hypervolume_history: map
            .get("hypervolume_history")
            .and_then(Value::as_array)
            .map(|a| a.iter().filter_map(Value::as_f64).collect())
            .unwrap_or_default(),
        batch: map.get("batch").filter(|v| !v.is_null()).cloned(),
    })
}

/// An L0 scorer that is a Python callable. `P2-05a`.
///
/// This is where `ADR-0001`'s boundary actually is: PyO3, not an OS process.
/// `DEP-LAND-004` deferred the out-of-process worker and that stays deferred —
/// the transport changes only *where* a scorer runs, and the crossing budget is
/// about how often Rust hands control to Python, which is measurable today and
/// is what the rule caps.
///
/// The counter that proves the rule cannot live here. A count kept on the Rust
/// side would tally this type's own bookkeeping and could read one while the
/// callable was invoked sixty-four times. It has to be kept by the callable,
/// which is the only thing that knows it was called.
pub struct PyL0Scorer<'py> {
    callable: pyo3::Bound<'py, pyo3::PyAny>,
    axis: usize,
    budget: usize,
}

impl<'py> PyL0Scorer<'py> {
    pub fn new(callable: pyo3::Bound<'py, pyo3::PyAny>, axis: usize, budget: usize) -> Self {
        Self {
            callable,
            axis,
            budget,
        }
    }
}

impl L0Scorer for PyL0Scorer<'_> {
    fn axis(&self) -> usize {
        self.axis
    }

    fn budget(&self) -> usize {
        self.budget
    }

    fn score_batch(&self, genotypes: &[&Genotype]) -> Result<Vec<Option<f64>>> {
        let sequences: Vec<String> = genotypes
            .iter()
            .map(|g| g.sequence().cds_string())
            .collect();
        // A raised scorer is not a batch-shape mismatch. `call1` failing used
        // to be relabeled `CF-LAND-012`, so a missing engine (`CF-STRU-003`)
        // and a programming error arrived as the same document, and the
        // orchestrator could not tell them apart. A registered `code`
        // attribute keeps its own code. Anything else is stashed and
        // re-raised by `search_tier3_segment` unchanged. `DEP-ORCH-004`.
        let returned = self.callable.call1((sequences,)).map_err(scorer_raised)?;
        // A raised exception is classified by its `code` attribute, and that
        // decision happens only in `scorer_raised` above. This arm is the
        // other failure: `call1` returned, and the value is not
        // `Vec<Option<f64>>`. That is a shape mismatch, so it is
        // `CF-LAND-012` unconditionally. PyO3 builds the `extract` error
        // itself; a `code` attribute on the returned value, or one that
        // happens to sit on that conversion error, is not consulted. Routing
        // this arm through `scorer_raised` would reclassify a bad return as
        // whatever that attribute said, or stash it and crash.
        let values: Vec<Option<f64>> = returned
            .extract()
            .map_err(|e| CfError::new(code::LAND_BATCH_UNEXPECTED, e.to_string()))?;
        // No length check here on purpose. `BatchProtocol::absorb` already
        // refuses a mismatch, and a second check in this type would be a
        // mechanism nobody could test: removing either one changes nothing,
        // so neither is covered. The negative control for (f) went silent for
        // exactly this reason, which is the fifth fixture-design rule doing
        // its job on code written the same afternoon it was applied elsewhere.
        Ok(values)
    }
}

thread_local! {
    /// The real scorer exception, beside the sentinel `CfError`.
    ///
    /// `L0Scorer::score_batch` (in `cf-search`) returns `CfError`. This crate
    /// is the only one allowed to name `PyErr`, so the exception cannot travel
    /// inside that `CfError`. `scorer_raised` stashes it here and returns a
    /// sentinel whose code is `CF-LAND-012`.
    ///
    /// Any caller on this path MUST call [`take_escaped_scorer_error`] before
    /// treating that `CfError` as real. `to_py_err` publishes every `CfError`
    /// as a problem document and does not look at this stash, so skipping the
    /// take silently republishes the sentinel as `CF-LAND-012`.
    /// `search_tier3_segment` is the only caller, and it takes.
    /// `test_p1_29_a_non_problem_scorer_error_still_crashes` is that call's
    /// regression: drop the take and the test records the sentinel instead of
    /// raising.
    ///
    /// A type that cannot be published (the `Observation::from_exact` shape)
    /// does not fit. The trait erases the value back to `CfError`, and
    /// `CfError::to_problem_details` is how every registered code crosses this
    /// boundary. A second recovery inside `to_py_err` would run first, empty
    /// the stash, and make the take dead — a control that deleted only the
    /// take would then stay silent. One recovery, at the take.
    static ESCAPED_SCORER_ERROR: RefCell<Option<pyo3::PyErr>> = const { RefCell::new(None) };
}

pub fn clear_escaped_scorer_error() {
    ESCAPED_SCORER_ERROR.with(|slot| *slot.borrow_mut() = None);
}

/// The scorer exception `score_batch` refused to relabel as `CF-LAND-012`.
///
/// Callers on this path MUST call this before treating the returned `CfError`
/// as real. The sentinel is a `CF-LAND-012`-shaped `CfError` with no problem
/// behind it; `to_py_err` will publish that shape. Skipping this take is how
/// the sentinel becomes a forged `CF-LAND-012`.
pub fn take_escaped_scorer_error() -> Option<pyo3::PyErr> {
    ESCAPED_SCORER_ERROR.with(|slot| slot.borrow_mut().take())
}

fn stash_escaped_scorer_error(err: pyo3::PyErr) {
    ESCAPED_SCORER_ERROR.with(|slot| *slot.borrow_mut() = Some(err));
}

/// A scorer that raised. Keep a registered CodonForge code; do not invent one.
fn scorer_raised(err: pyo3::PyErr) -> CfError {
    if let Some(coded) = problem_code_from_python(&err) {
        return coded;
    }
    stash_escaped_scorer_error(err);
    // Sentinel only. Code `CF-LAND-012`, and not a problem document. The
    // caller MUST `take_escaped_scorer_error` before publishing this value.
    // Treating it as a real `CfError` republishes the sentinel as
    // `CF-LAND-012`.
    CfError::new(
        code::LAND_BATCH_UNEXPECTED,
        "the L0 scorer raised an exception that is not a CodonForge problem",
    )
}

fn problem_code_from_python(err: &pyo3::PyErr) -> Option<CfError> {
    pyo3::Python::with_gil(|py| {
        let value = err.value_bound(py);
        let code: String = value.getattr("code").ok()?.extract().ok()?;
        let entry = lookup(&code)?;
        let detail: String = match value.getattr("detail") {
            Ok(detail) => detail.extract().unwrap_or_else(|_| err.to_string()),
            Err(_) => err.to_string(),
        };
        Some(CfError::new(entry.code, detail))
    })
}
