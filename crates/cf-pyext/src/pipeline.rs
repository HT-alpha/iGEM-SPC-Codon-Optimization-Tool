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

use cf_constraints::ConstraintChecker;
use cf_landscape::terms::constraints::ViolationCountTerm;
use cf_landscape::terms::usage::{CaiTerm, CodonUsageTable, GcContentTerm};
use cf_landscape::ConstructLayout;
use cf_landscape_api::provider::EvalContext;
use cf_landscape_api::provider::TermRegistry;
use cf_landscape_api::query::Budget;
use cf_landscape_api::HostContext;
use cf_search::tier3::{
    seed_population, Individual, ObjectiveEvaluator, Selection, Tier3Checkpoint, Tier3Config,
    Tier3Request,
};
use cf_search::{
    tier1_seed, tier2_block_ca, CodonWeights, Genotype, OperatorMix, SearchRng, Tier1Config,
    Tier1Request, Tier2Config, Tier2Plan, Tier2Request,
};
use cf_types::{code, CfError, CodonSequence, Digest, HostId, Result, TermId};
use serde_json::{json, Value};

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
                registry.register_default(Box::new(CaiTerm::new(table)?))?;
            }
            "usage.gc_content" => {
                registry.register_default(Box::new(GcContentTerm::new()))?;
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
                        "no built-in provider supplies {other}; Phase 1 ships usage.cai, \
                         usage.gc_content and constraints.violation_count, and a job is \
                         refused rather than run against fewer terms than it asked for"
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
        out.insert(term.clone(), json!(registration.provider.score(&ctx)?));
    }
    Ok(Value::Object(out))
}

/// The scalarized objective tier 1 and tier 2 maximize, for one sequence.
///
/// The same formula the tiers use, including reading `provider.direction()` for
/// the sense — so a caller comparing two sequences with this is comparing them
/// the way the search did, rather than by a second convention that could
/// disagree with it.
pub fn scalarized_score(inputs: &TierInputs<'_>, cds: &str) -> Result<f64> {
    let scores = score_terms(inputs, cds)?;
    let directions = declared_directions(inputs)?;
    let mut total = 0.0;
    for (index, term) in inputs.terms.iter().enumerate() {
        let value = scores.get(term).and_then(Value::as_f64).ok_or_else(|| {
            CfError::new(code::LAND_TERM_UNKNOWN, format!("{term} did not score"))
        })?;
        let sense = match directions.get(term).and_then(Value::as_str) {
            Some("maximize") => 1.0,
            _ => -1.0,
        };
        total += sense * inputs.weights.get(index).copied().unwrap_or(1.0) * value;
    }
    Ok(total)
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

    let evaluator = RegistryEvaluator {
        registry: &registry,
        host: &host,
        terms: term_ids(inputs.terms)?,
        senses: term_senses(&registry, inputs.terms)?,
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
        "hypervolume_history": report.checkpoint.hypervolume_history,
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
        members.push(json!({
            "sequence": member.genotype.sequence().to_json()?,
            "objectives": member.objectives,
            "violation": member.violation,
        }));
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
        members.push(Individual {
            genotype: Genotype::new(sequence),
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
    })
}
