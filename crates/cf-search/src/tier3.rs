//! Tier 3: the NSGA-II population loop.
//!
//! Ticket P2-02. Specification: `03-nsga-loop-and-operators.md` §4, §7.
//!
//! ## What this is and is not
//!
//! This is the `M ≤ 3` path: constraint-domination, fast non-dominated
//! sorting, and crowding-distance selection. Reference-direction niching is
//! `P2-03`, surrogate-assisted evaluation is `P2-06`, and front certification
//! is `P2-07`. §7's pseudocode calls all three; each is *absent* here rather
//! than stubbed, because a stub that returns a plausible number is
//! indistinguishable from a working component until something depends on it.
//!
//! ## Why the output is a front and not a winner
//!
//! Tiers 1 and 2 produce a single incumbent because they optimize a
//! scalarization. Tier 3 produces a *set* — the mutually non-dominated
//! individuals — and collapsing it to a scalar winner would throw away the
//! Pareto structure the whole tier exists to produce. `CF-ORCH-007`'s
//! no-regression guarantee is extended to fronts rather than applied by
//! reduction; see `CandidateArchive` on the Python side.
//!
//! ## Determinism
//!
//! Every random draw derives from `SearchRng::stream(root, index)` on a
//! structural index — the generation, and within a generation the offspring
//! pair. Tier 3 has its own reserved root (`DEP-DET-003`), because its natural
//! per-individual indices `0..N` are exactly the index space the generation
//! loop already uses off the job root.

use cf_landscape_api::batch::{BatchProtocol, BatchResults};

use crate::rationing::{Observation, Rationer};

use cf_constraints::repair::{ActionableViolation, RepairTarget};
use cf_constraints::{ConstraintChecker, ConstructView, KmerMasking, TabuLedger};
use cf_metrics::{
    associate, crowding_distance, das_dennis, das_dennis_two_layer, fast_nondominated_sort,
    niching_select, normalize, placeable, Ranked,
};
use cf_types::{code, CfError, Codon, Result};

use crate::context::CodonWeights;
use crate::genotype::Genotype;
use crate::operators::crossover::{BlockUniform, CrossoverOperator};
use crate::operators::{MutationContext, OperatorMix};
use crate::rng::SearchRng;

/// One individual: a genotype plus everything ranking needs about it.
#[derive(Clone, Debug)]
pub struct Individual {
    pub genotype: Genotype,
    /// Minimized objective values. `None` is `UNKNOWN`, never imputed.
    pub objectives: Vec<Option<f64>>,
    /// Total *actionable* violation. `DEP-SRCH-008`.
    pub violation: f64,
    /// Genome transmitted to offspring when it is not the scored sequence.
    ///
    /// `None` is Lamarckian: repair replaces the genotype, and that genotype
    /// is what was scored. `Some` is Baldwinian: `genotype` is the repaired
    /// sequence that was scored, and this is the unrepaired genome the next
    /// generation breeds from.
    pub inherited: Option<Genotype>,
}

impl Individual {
    pub fn ranked(&self) -> Ranked {
        Ranked {
            objectives: self.objectives.clone(),
            violation: self.violation,
        }
    }

    pub fn is_feasible(&self) -> bool {
        self.violation <= 0.0
    }
}

/// Scores an individual. Supplied by the caller so `cf-search` never names a
/// concrete predictor — the same seam tiers 1 and 2 use.
pub trait ObjectiveEvaluator {
    /// Minimized objective values for a genotype. `None` for `UNKNOWN`.
    fn evaluate(&self, genotype: &Genotype) -> Result<Vec<Option<f64>>>;

    /// How many objectives `evaluate` returns.
    fn objective_count(&self) -> usize;
}

/// Which environmental-selection rule the splitting front uses.
///
/// §2's switch. Crowding distance degrades badly as `M` grows — nearly every
/// individual becomes non-dominated, so selection pressure comes almost
/// entirely from the diversity operator, and a sum of per-axis spacings stops
/// discriminating once most points sit on the front's boundary. Reference-
/// direction niching supplies structured coverage instead.
///
/// This decides *survival* only. It is deliberately not a ranking: introducing
/// a second notion of "best" here is how the archive's `front()` and `best()`
/// would come to disagree, and they are the only two answers the system gives.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Selection {
    /// NSGA-II crowding distance. `M ≤ 3`.
    CrowdingDistance,
    /// NSGA-III reference-direction niching. `M ≥ 4`.
    ReferenceDirections,
}

impl Selection {
    /// §2's rule, in one place so the switch cannot be spelled twice.
    pub fn for_objectives(objectives: usize) -> Self {
        if objectives >= 4 {
            Self::ReferenceDirections
        } else {
            Self::CrowdingDistance
        }
    }
}

/// When a run stops, and why.
///
/// `P2-04`. §7 terminates on hypervolume stagnation by default "because it is
/// the only signal that tracks the thing the user cares about"; the generation
/// and budget caps are backstops, and saying which one fired is the difference
/// between "converged" and "ran out of time".
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Termination {
    /// Hypervolume stopped improving over the stagnation window.
    HypervolumeStagnation,
    /// The generation cap was reached.
    GenerationLimit,
    /// The caller's segment ended; the run can continue.
    SegmentComplete,
}

impl Termination {
    /// The `pareto_front.summary.termination_reason` value.
    pub fn as_summary(&self) -> &'static str {
        match self {
            Self::HypervolumeStagnation => "hypervolume_stagnation",
            Self::GenerationLimit => "generation_limit",
            Self::SegmentComplete => "generation_limit",
        }
    }

    pub fn is_final(&self) -> bool {
        !matches!(self, Self::SegmentComplete)
    }
}

/// Whether a hypervolume history has stopped improving.
///
/// Relative rather than absolute, because hypervolume has no natural scale:
/// its magnitude depends on the reference point and the objective units, so an
/// absolute epsilon would mean "converged" for one job and "barely started"
/// for another.
///
/// A window shorter than `window` is never stagnant — a run has to have had
/// the chance to improve before failing to is evidence of anything.
pub fn stagnated(history: &[f64], window: usize, relative_epsilon: f64) -> bool {
    if window == 0 || history.len() <= window {
        return false;
    }
    let recent = &history[history.len() - window - 1..];
    let baseline = recent[0];
    let latest = recent[recent.len() - 1];
    let scale = baseline.abs().max(f64::MIN_POSITIVE);
    (latest - baseline) / scale <= relative_epsilon
}

#[derive(Clone, Copy, Debug)]
pub struct Tier3Config {
    /// Population size `N`.
    pub population: usize,
    pub max_generations: u32,
    /// Generations of no relative hypervolume improvement before stopping.
    ///
    /// §7's default is 30 at `N = 212`. Stated by the caller rather than
    /// defaulted here for a small population, because a window longer than a
    /// short run can never fire and one shorter than the search's own noise
    /// stops it before it has begun.
    pub stagnation_window: usize,
    pub stagnation_epsilon: f64,
    /// Crossover probability, `p_c` in §7.
    pub crossover_rate: f64,
    pub tabu_leak: f64,
    pub region: usize,
    pub masking: KmerMasking,
    /// `G_s` in §6 — generations between surrogate refits.
    pub surrogate_refit_period: u64,
}

impl Tier3Config {
    /// `N` and `p_c` from §7's worked loop.
    pub fn from_planning(population: usize) -> Self {
        Self {
            population,
            max_generations: 20,
            stagnation_window: 8,
            stagnation_epsilon: 1e-3,
            crossover_rate: 0.9,
            tabu_leak: 0.0,
            region: 0,
            masking: KmerMasking::default(),
            surrogate_refit_period: 10,
        }
    }
}

/// A resumable point in a tier-3 run.
///
/// Everything the reported front depends on, per the rule three Phase 1
/// omissions established. In particular the *objective vectors* travel with
/// the population: recomputing them on resume would be correct only if
/// evaluation were pure and cheap, and it is neither — a partial evaluation
/// (`UNKNOWN`) is a fact about a budget that has already been spent, and
/// re-running it could turn an `UNKNOWN` into a value and change the ranking.
#[derive(Clone, Debug)]
pub struct Tier3Checkpoint {
    pub generation: u32,
    pub members: Vec<Individual>,
    pub tabu: Vec<(usize, usize, Codon, u32)>,
    pub evaluations: u64,
    pub repairs_attempted: u64,
    pub repairs_succeeded: u64,
    /// One entry per generation actually run, across every segment.
    ///
    /// Travels in the checkpoint because termination depends on it: a resumed
    /// run that started its history empty could not stagnate until it had run
    /// the whole window again, so a job checkpointed near convergence would
    /// keep going where an uninterrupted one stopped.
    pub hypervolume_history: Vec<f64>,
    /// The batch protocol's snapshot, when an L0 term is running.
    ///
    /// `P2-05a`. Carries the deferral queue and the stranded set, so a resumed
    /// run still knows which individuals the budget never reached and
    /// `CF-SRCH-016` survives an interruption.
    pub batch: Option<serde_json::Value>,
}

#[derive(Clone, Debug)]
pub struct Tier3Report {
    pub checkpoint: Tier3Checkpoint,
    /// Which selection rule ran. Reported rather than inferred, so a test can
    /// assert the `M ≥ 4` switch happened instead of guessing from behaviour.
    pub selection: Selection,
    /// The non-dominated set, best front first, in a stable order.
    pub front: Vec<Individual>,
    /// Why the run stopped.
    pub termination: Termination,
    pub generations_run: u32,
    pub warnings: Vec<CfError>,
    /// Individuals the L0 budget never reached. Drives `CF-SRCH-016`.
    pub stranded: Vec<usize>,
}

pub struct Tier3Request<'a> {
    pub checker: &'a ConstraintChecker,
    pub mix: &'a OperatorMix,
    pub weights: &'a CodonWeights,
    pub evaluator: &'a dyn ObjectiveEvaluator,
    pub config: &'a Tier3Config,
    /// Already namespaced to tier 3 by the caller. `DEP-DET-003`.
    pub root_seed: u64,
    /// The out-of-core term, scored once per generation. `P2-05a`.
    ///
    /// `None` when every objective is in-core, which is the case tier 3 ran
    /// under until this ticket: `evaluate_one` computed everything per
    /// individual and `BatchProtocol` was never reached.
    pub l0: Option<&'a dyn L0Scorer>,
}

/// An objective too expensive to evaluate per move.
///
/// `ADR-0001` caps boundary crossings at one per generation, so this trait has
/// exactly one method and it takes the whole batch. A per-item method would
/// invite a caller to loop, and a loop over `N = 64` is a 64x budget violation
/// that nothing would report — the run would simply be slow, and the cause
/// would be invisible.
pub trait L0Scorer {
    /// Which objective index this scorer fills.
    fn axis(&self) -> usize;

    /// The evaluation budget for one generation (`E_max`).
    fn budget(&self) -> usize;

    /// Score every genotype handed over, in order.
    ///
    /// `None` for an item the scorer could not evaluate — `UNKNOWN`, never a
    /// placeholder. Must return one value per item: the protocol refuses a
    /// length mismatch rather than aligning by position, because a shifted
    /// alignment attributes one sequence's score to another and every value
    /// still looks plausible.
    fn score_batch(&self, genotypes: &[&Genotype]) -> Result<Vec<Option<f64>>>;
}

impl std::fmt::Debug for Tier3Request<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Tier3Request")
            .field("population", &self.config.population)
            .field("root_seed", &self.root_seed)
            .finish()
    }
}

/// Total *actionable* violation for a genotype.
///
/// `DEP-SRCH-008`, reconciling §4 with `DEP-CONS-002`/`DEP-CONS-004`. A
/// violation lying wholly inside immutable sequence cannot be removed by any
/// synonymous assignment; charging it here would rank the population by a
/// quantity none of it can influence, and since constraint-domination makes
/// any feasible individual beat any infeasible one, that quantity would drive
/// the entire selection.
pub fn actionable_violation(genotype: &Genotype, checker: &ConstraintChecker) -> f64 {
    let construct = assemble(genotype);
    let start = genotype.cds_offset_nt();
    let span = (start, start + genotype.len() * 3);
    let view = ConstructView::new(&construct, vec![span]).with_editable_spans(
        cf_constraints::repair::editable_spans(&[(start, genotype.sequence().clone())]),
    );
    checker
        .check(&view)
        .iter()
        .filter(|violation| violation.actionable)
        .count() as f64
}

fn assemble(genotype: &Genotype) -> Vec<u8> {
    let sequence = genotype.sequence();
    let mut out = sequence.five_prime_flank().to_vec();
    for codon in sequence.codons() {
        out.extend_from_slice(&codon.to_bytes());
    }
    out.extend_from_slice(sequence.three_prime_flank());
    out
}

/// Score a genotype into an `Individual`.
pub fn evaluate_one(genotype: Genotype, request: &Tier3Request<'_>) -> Result<Individual> {
    let mut objectives = request.evaluator.evaluate(&genotype)?;
    // An L0 term is an *additional* objective: the evaluator supplies the
    // in-core axes and the scorer supplies one more. The slot is created here,
    // as `UNKNOWN`, rather than by the batch step — a slot that appeared only
    // once a batch had been scored would make the vector's width depend on
    // whether the budget happened to reach an individual.
    if let Some(scorer) = request.l0 {
        let axis = scorer.axis();
        if axis < request.evaluator.objective_count() {
            return Err(CfError::new(
                code::SRCH_TIER3_CONFIG,
                format!(
                    "the L0 scorer claims axis {axis}, which the evaluator already \
                     fills; an out-of-core term overwriting an in-core one would \
                     silently discard the in-core value"
                ),
            ));
        }
        objectives.resize(axis + 1, None);
    }
    if objectives.len() != objective_width(request) {
        return Err(CfError::new(
            code::SRCH_TIER3_CONFIG,
            format!(
                "evaluator returned {} objectives but tier 3 expects {}",
                objectives.len(),
                objective_width(request)
            ),
        ));
    }
    let violation = actionable_violation(&genotype, request.checker);
    Ok(Individual {
        genotype,
        objectives,
        violation,
        inherited: None,
    })
}

fn breeding_genotype(individual: &Individual) -> &Genotype {
    individual
        .inherited
        .as_ref()
        .unwrap_or(&individual.genotype)
}

/// Start a run from a seed population.
pub fn seed_population(seeds: &[Genotype], request: &Tier3Request<'_>) -> Result<Tier3Checkpoint> {
    if seeds.is_empty() {
        return Err(CfError::new(
            code::SRCH_TIER3_CONFIG,
            "tier 3 needs at least one seed; an empty population has no front to report",
        ));
    }
    let mut members = Vec::new();
    // Cycled rather than padded with copies of the best: duplicates would be
    // suppressed by `DEDUP` anyway, and starting from every distinct seed is
    // what tier 1 produced them for.
    for index in 0..request.config.population {
        members.push(evaluate_one(seeds[index % seeds.len()].clone(), request)?);
    }
    Ok(Tier3Checkpoint {
        generation: 0,
        members,
        tabu: Vec::new(),
        evaluations: request.config.population as u64,
        repairs_attempted: 0,
        repairs_succeeded: 0,
        batch: None,
        hypervolume_history: Vec::new(),
    })
}

/// Advance `generations` complete generations. `M ≤ 3` path.
pub fn run_generations(
    checkpoint: &Tier3Checkpoint,
    request: &Tier3Request<'_>,
    generations: u32,
) -> Result<Tier3Report> {
    let objectives = objective_width(request);
    let selection = Selection::for_objectives(objectives);

    let mut ledger = TabuLedger::new(request.config.tabu_leak);
    for (region, index, codon, remaining) in &checkpoint.tabu {
        ledger.record(*region, *index, *codon, *remaining);
    }

    let mut members = checkpoint.members.clone();
    let mut history = checkpoint.hypervolume_history.clone();
    // Both travel in the checkpoint. The protocol carries the deferral queue
    // and the stranded set behind `CF-SRCH-016`; losing either on resume would
    // make a long job's warning quietly disappear at exactly the point it is
    // most likely to be interrupted.
    let mut protocol = match &checkpoint.batch {
        Some(document) => BatchProtocol::resume(document)?,
        None => BatchProtocol::new(),
    };
    let mut rationer = Rationer::new(
        request.l0.map(|s| s.budget()).unwrap_or(0),
        request.config.surrogate_refit_period,
    );
    let mut termination = Termination::SegmentComplete;
    let mut ran = 0u32;
    let mut evaluations = checkpoint.evaluations;
    let mut repairs_attempted = checkpoint.repairs_attempted;
    let mut repairs_succeeded = checkpoint.repairs_succeeded;

    let first = checkpoint.generation;
    // The cap is a bound on the *run*, not on the segment, and it is applied
    // before a generation runs rather than after. Checking it afterwards let a
    // resumed segment run one generation past it: an uninterrupted request for
    // nine against a cap of eight stopped at eight, while the same run split
    // 8 + 1 ran the eighth, resumed, and ran a ninth. Both were internally
    // consistent and they disagreed with each other, which is the exact shape
    // of bug `TST-SRCH-016` exists to catch.
    let last = (first + generations).min(request.config.max_generations);
    if last <= first {
        termination = Termination::GenerationLimit;
    }
    for generation in first..last {
        // A stream per absolute generation, so generation k sees the same
        // randomness whether it is the first of a resumed segment or the
        // fiftieth of an uninterrupted run.
        let mut rng = SearchRng::stream(request.root_seed, u64::from(generation));

        // Baldwinian genomes stay unrepaired, so the same individual is
        // repaired again at the start of the next generation. Lamarckian
        // genomes were replaced, so this loop does not see them. The repair
        // count is the evidence the two arms took different paths.
        if !repairs_are_lamarckian(request) {
            let due: Vec<usize> = members
                .iter()
                .enumerate()
                .filter(|(_, member)| {
                    member
                        .inherited
                        .as_ref()
                        .is_some_and(|genome| actionable_violation(genome, request.checker) > 0.0)
                })
                .map(|(index, _)| index)
                .collect();
            for index in due {
                let genome = members[index].inherited.clone().expect("filtered above");
                let mut repair_rng = SearchRng::stream(
                    request.root_seed,
                    10_000_000 + u64::from(generation) * 1_000 + index as u64,
                );
                let choice = repair_child(genome, request, &mut ledger, &mut repair_rng)?;
                repairs_attempted += choice.attempted;
                repairs_succeeded += choice.succeeded;
                if choice.succeeded > 0 {
                    let mut scored = evaluate_one(choice.phenotype, request)?;
                    scored.inherited = members[index].inherited.clone();
                    members[index] = scored;
                }
            }
        }

        let mut offspring: Vec<Individual> = Vec::new();
        let mut pair = 0u64;
        while offspring.len() < request.config.population {
            // A sub-stream per offspring pair, derived from the generation
            // stream's seed rather than drawn from it, so a pair sees the same
            // randomness however many pairs were bred before it — the property
            // parallel breeding would need and the one an iteration counter
            // would break.
            let mut pair_rng = SearchRng::stream(rng.next_u64(), pair);
            pair += 1;

            let a = tournament(&members, &mut pair_rng);
            let b = tournament(&members, &mut pair_rng);

            let bred = breed(
                breeding_genotype(&members[a]),
                breeding_genotype(&members[b]),
                request,
                &ledger,
                &mut pair_rng,
            )?;
            for child in bred {
                let choice = repair_child(child, request, &mut ledger, &mut pair_rng)?;
                repairs_attempted += choice.attempted;
                repairs_succeeded += choice.succeeded;
                if offspring.len() < request.config.population {
                    let mut individual = evaluate_one(choice.phenotype, request)?;
                    if choice.succeeded > 0 && !repairs_are_lamarckian(request) {
                        individual.inherited = Some(choice.breed);
                    }
                    offspring.push(individual);
                    evaluations += 1;
                }
            }
        }

        let mut combined = members;
        combined.extend(offspring);

        // The one crossing this generation gets. `P2-05a`.
        //
        // Placed before selection, not after: selection ranks on the objective
        // vector, so scoring afterwards would rank the population on a
        // dimension that was still `UNKNOWN` when the ranking happened, and the
        // scores would describe survivors chosen without them.
        if let Some(scorer) = request.l0 {
            score_l0_batch(
                &mut combined,
                scorer,
                &mut protocol,
                &mut rationer,
                &mut rng,
            )?;
        }

        members = select(combined, request.config.population, selection);
        ledger.tick();
        ran += 1;

        history.push(front_hypervolume(&members));
        if stagnated(
            &history,
            request.config.stagnation_window,
            request.config.stagnation_epsilon,
        ) {
            termination = Termination::HypervolumeStagnation;
            break;
        }
        if generation + 1 >= request.config.max_generations {
            termination = Termination::GenerationLimit;
            break;
        }
    }

    let front = nondominated(&members);
    Ok(Tier3Report {
        checkpoint: Tier3Checkpoint {
            generation: first + ran,
            members,
            tabu: ledger.entries().collect(),
            evaluations,
            repairs_attempted,
            repairs_succeeded,
            hypervolume_history: history,
            batch: request.l0.map(|_| protocol.snapshot()),
        },
        front,
        termination,
        selection,
        stranded: request.l0.map(|_| protocol.stranded()).unwrap_or_default(),
        generations_run: ran,
        warnings: request.mix.warnings(),
    })
}

/// Binary tournament under constraint-domination.
///
/// A tie goes to the lower index rather than to the RNG. The draw already
/// happened when the two contestants were chosen; breaking the tie randomly
/// again would consume randomness in a way that depends on how many ties
/// occurred, which makes the stream position a function of the population's
/// contents rather than of the generation.
fn tournament(members: &[Individual], rng: &mut SearchRng) -> usize {
    let a = rng.below(members.len());
    let b = rng.below(members.len());
    match cf_metrics::constraint_dominance(&members[a].ranked(), &members[b].ranked()) {
        cf_metrics::Dominance::Left => a,
        cf_metrics::Dominance::Right => b,
        cf_metrics::Dominance::Neither => a.min(b),
    }
}

fn breed(
    a: &Genotype,
    b: &Genotype,
    request: &Tier3Request<'_>,
    ledger: &TabuLedger,
    rng: &mut SearchRng,
) -> Result<Vec<Genotype>> {
    let crossover = BlockUniform::default();
    let mut children = if rng.next_f64() < request.config.crossover_rate {
        let (first, second) = crossover.cross(a, b, rng)?;
        vec![first, second]
    } else {
        vec![a.clone(), b.clone()]
    };

    for child in children.iter_mut() {
        let construct = assemble(child);
        let start = child.cds_offset_nt();
        let view = ConstructView::new(&construct, vec![(start, start + child.len() * 3)])
            .with_editable_spans(cf_constraints::repair::editable_spans(&[(
                start,
                child.sequence().clone(),
            )]));
        let found = request.checker.check(&view);
        let actionable: Vec<ActionableViolation<'_>> =
            found.iter().filter_map(ActionableViolation::new).collect();

        let ctx = MutationContext::new(child, request.weights)
            .with_violations(&actionable)
            .with_tabu(ledger)
            .with_region(request.config.region);
        let proposal = request.mix.propose(&ctx, rng)?.compacted(child);
        *child = child.with_move(&proposal)?;
    }
    Ok(children)
}

struct RepairedChild {
    /// Sequence that was scored. The repaired form when repair succeeded.
    phenotype: Genotype,
    /// Sequence the next generation mutates. The unrepaired genome when
    /// Baldwinian repair succeeded; otherwise the phenotype.
    breed: Genotype,
    attempted: u64,
    succeeded: u64,
}

/// The policy flag, read in one place so a control can silence it.
fn repairs_are_lamarckian(request: &Tier3Request<'_>) -> bool {
    request.checker.set().repair_policy.lamarckian
}

fn repair_child(
    genotype: Genotype,
    request: &Tier3Request<'_>,
    ledger: &mut TabuLedger,
    rng: &mut SearchRng,
) -> Result<RepairedChild> {
    if actionable_violation(&genotype, request.checker) <= 0.0 {
        return Ok(RepairedChild {
            phenotype: genotype.clone(),
            breed: genotype,
            attempted: 0,
            succeeded: 0,
        });
    }

    let mut target = RepairTarget::new(
        assemble(&genotype),
        vec![(genotype.cds_offset_nt(), genotype.sequence().clone())],
    )?;
    let mut repair_rng = cf_constraints::RepairRng::new(rng.next_u64());
    let outcome = cf_constraints::repair::repair(
        &mut target,
        request.checker,
        request.config.masking,
        &NoObjective,
        ledger,
        &mut repair_rng,
    )?;

    if outcome.residual.is_empty() {
        let phenotype = Genotype::new(target.regions()[0].1.clone());
        // Lamarckian writes the repair back. Baldwinian scores `phenotype`
        // and breeds from the genome repair did not replace, so the same
        // violation is repaired again next generation.
        let breed = if repairs_are_lamarckian(request) {
            phenotype.clone()
        } else {
            genotype
        };
        Ok(RepairedChild {
            phenotype,
            breed,
            attempted: 1,
            succeeded: 1,
        })
    } else {
        // Repair could not clear every blocking violation. The individual is
        // kept, not discarded: §4 is explicit that infeasible individuals
        // survive as stepping stones across a constraint barrier, and
        // constraint-domination is what applies the pressure instead.
        Ok(RepairedChild {
            phenotype: genotype.clone(),
            breed: genotype,
            attempted: 1,
            succeeded: 0,
        })
    }
}

/// Repair with no objective preference, since tier 3 ranks rather than scores.
#[derive(Debug)]
struct NoObjective;

impl cf_constraints::repair::ObjectiveDelta for NoObjective {
    fn delta(&self, _region: usize, _codon_index: usize, _codon: Codon) -> Option<f64> {
        // `None`, not `Some(0.0)`: "no preference" and "every fix is equally
        // good" are different statements, and repair's tie-breaking reads them
        // differently. Tier 3 ranks rather than scores, so it genuinely has no
        // preference to express here.
        None
    }
}

/// Environmental selection: whole fronts in order, then the splitting front.
///
/// The two strategies differ only in how the splitting front is cut. Front
/// order itself is non-domination and is the same for both, which is what
/// keeps `M ≥ 4` from being a different algorithm rather than a different
/// diversity operator.
fn select(combined: Vec<Individual>, keep: usize, selection: Selection) -> Vec<Individual> {
    let ranked: Vec<Ranked> = combined.iter().map(Individual::ranked).collect();
    let fronts = fast_nondominated_sort(&ranked);

    let mut chosen: Vec<usize> = Vec::with_capacity(keep);
    for front in &fronts {
        if chosen.len() + front.len() <= keep {
            chosen.extend(front.iter().copied());
            continue;
        }
        let slots = keep - chosen.len();
        let cut = match selection {
            Selection::CrowdingDistance => by_crowding(&ranked, front, slots),
            Selection::ReferenceDirections => {
                by_reference_directions(&ranked, &chosen, front, slots)
            }
        };
        chosen.extend(cut);
        break;
    }

    chosen.sort_unstable();
    chosen.into_iter().map(|i| combined[i].clone()).collect()
}

/// NSGA-II: take the least crowded, ties by index.
fn by_crowding(ranked: &[Ranked], front: &[usize], slots: usize) -> Vec<usize> {
    let distance = crowding_distance(ranked, front);
    let mut order: Vec<usize> = (0..front.len()).collect();
    order.sort_by(|&a, &b| {
        distance[b]
            .partial_cmp(&distance[a])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| front[a].cmp(&front[b]))
    });
    order.iter().take(slots).map(|&p| front[p]).collect()
}

/// NSGA-III: normalize, associate with reference directions, fill the
/// emptiest niches first.
///
/// Members with an `UNKNOWN` term cannot be placed in objective space at all,
/// so they are not associated with a direction — associating them would be
/// inventing the coordinate. They are taken last, in index order, which is §4's
/// "back of its feasibility class" applied to the splitting front. They are not
/// *dropped*: an individual whose provider missed a budget has not been shown
/// to be bad, only unmeasured.
fn by_reference_directions(
    ranked: &[Ranked],
    already: &[usize],
    front: &[usize],
    slots: usize,
) -> Vec<usize> {
    let (placeable_admitted, _) = placeable(ranked, already);
    let (placeable_front, unplaceable) = placeable(ranked, front);

    let objectives = ranked.first().map(|r| r.objectives.len()).unwrap_or(0);
    let directions = reference_directions(objectives, front.len().max(1));

    // Normalized over the admitted members and the splitting front together:
    // the hyperplane is a property of the generation, not of the slice being
    // cut, and normalizing the slice alone would rescale every generation
    // differently.
    let considered: Vec<usize> = placeable_admitted
        .iter()
        .chain(placeable_front.iter())
        .copied()
        .collect();
    let raw: Vec<Vec<f64>> = considered
        .iter()
        .map(|&i| {
            ranked[i]
                .objectives
                .iter()
                .map(|v| v.unwrap_or(0.0))
                .collect()
        })
        .collect();
    let normalization = normalize(&raw, None);
    let points: Vec<Vec<f64>> = raw.iter().map(|p| normalization.apply(p)).collect();

    let associations_by_position = associate(&points, &directions);
    // Re-indexed onto the population, so `niching_select` speaks the same
    // indices as everything else here.
    let mut associations = vec![
        cf_metrics::Association {
            reference: 0,
            distance: f64::INFINITY,
        };
        ranked.len()
    ];
    for (position, &index) in considered.iter().enumerate() {
        associations[index] = associations_by_position[position];
    }

    let mut cut = niching_select(
        &associations,
        &placeable_admitted,
        &placeable_front,
        directions.len(),
        slots,
    );
    for index in unplaceable {
        if cut.len() >= slots {
            break;
        }
        cut.push(index);
    }
    cut.sort_unstable();
    cut
}

/// Reference directions sized to the population, per §3's table.
///
/// `p` is the smallest number of divisions giving at least as many directions
/// as there are slots to fill, and `M ≥ 6` uses the two-layer construction
/// because a single layer then puts almost every direction on the simplex
/// boundary, leaving the interior — where balanced trade-offs live —
/// uncovered.
fn reference_directions(objectives: usize, target: usize) -> Vec<Vec<f64>> {
    if objectives == 0 {
        return vec![vec![1.0]];
    }
    for divisions in 1..=30usize {
        let single = das_dennis(objectives, divisions).unwrap_or_default();
        if single.len() >= target {
            if objectives >= 6 {
                let inner = divisions.saturating_sub(2).max(1);
                if let Ok(two_layer) = das_dennis_two_layer(objectives, divisions, inner) {
                    return two_layer;
                }
            }
            return single;
        }
    }
    das_dennis(objectives, 30).unwrap_or_else(|_| vec![vec![1.0; objectives]])
}

/// Hypervolume of the feasible non-dominated set, against a nadir-derived
/// reference point.
///
/// The reference point is the per-axis worst of the population plus one, so it
/// dominates every member and the volume is positive and comparable *within* a
/// run. It is deliberately not comparable *across* runs: a comparison like
/// `EXP-05`'s needs a fixed reference point supplied by the benchmark harness
/// (`P2-19`), and using this one there would compare two differently-scaled
/// numbers. Here it is only a stagnation signal.
fn front_hypervolume(members: &[Individual]) -> f64 {
    let feasible: Vec<&Individual> = members
        .iter()
        .filter(|m| m.is_feasible() && m.objectives.iter().all(Option::is_some))
        .collect();
    if feasible.is_empty() {
        return 0.0;
    }
    let objectives = feasible[0].objectives.len();
    let points: Vec<Vec<f64>> = feasible
        .iter()
        .map(|m| m.objectives.iter().map(|v| v.unwrap()).collect())
        .collect();
    let reference: Vec<f64> = (0..objectives)
        .map(|axis| {
            points
                .iter()
                .map(|p| p[axis])
                .fold(f64::NEG_INFINITY, f64::max)
                + 1.0
        })
        .collect();
    cf_metrics::hypervolume(&points, &reference).unwrap_or(0.0)
}

/// The mutually non-dominated members, in ascending member order.
pub fn nondominated(members: &[Individual]) -> Vec<Individual> {
    let ranked: Vec<Ranked> = members.iter().map(Individual::ranked).collect();
    match fast_nondominated_sort(&ranked).first() {
        Some(front) => front.iter().map(|&i| members[i].clone()).collect(),
        None => Vec::new(),
    }
}

/// Score the L0 axis for one generation: rank, cut, cross once, absorb.
///
/// `P2-05a`. This is the function that joins three pieces built separately and
/// never connected — `P2-06`'s ranker, `P2-05`'s protocol, and the generation
/// loop — and the order matters at every step:
///
/// * the **ranker** decides priority and nothing else (`DEP-SRCH-011`);
/// * the **protocol** applies `E_max`, puts last generation's deferred items at
///   the head, and refuses a second crossing in the same generation;
/// * the **scorer** is called exactly once, with the chosen genotypes.
///
/// Individuals not in the batch keep `None` on that axis. `UNKNOWN` is a
/// first-class value in constraint-domination, so they are ranked honestly
/// rather than imputed — which is what makes rationing sound rather than a
/// silent approximation.
fn score_l0_batch(
    population: &mut [Individual],
    scorer: &dyn L0Scorer,
    protocol: &mut BatchProtocol,
    rationer: &mut Rationer,
    rng: &mut SearchRng,
) -> Result<()> {
    let axis = scorer.axis();
    if population.is_empty() {
        return Ok(());
    }

    // Rank on what is known so far: the in-core objectives, with the L0 axis
    // dropped rather than filled in. Ranking on a placeholder would let the
    // placeholder choose who gets evaluated.
    let known: Vec<Vec<f64>> = population
        .iter()
        .map(|individual| {
            individual
                .objectives
                .iter()
                .enumerate()
                .filter(|(index, _)| *index != axis)
                .filter_map(|(_, value)| *value)
                .collect()
        })
        .collect();
    let violations: Vec<f64> = population.iter().map(|i| i.violation).collect();

    let ranked = rationer.rank(&known, &violations, rng);
    let (_, chosen) = protocol.step(&ranked.ranked, scorer.budget())?;

    if !chosen.is_empty() {
        let genotypes: Vec<&Genotype> = chosen
            .iter()
            .map(|index| &population[*index].genotype)
            .collect();

        // The crossing.
        let values = scorer.score_batch(&genotypes)?;

        protocol.absorb(&BatchResults {
            generation: protocol.generation(),
            values: values.clone(),
        })?;

        for (index, value) in chosen.iter().zip(&values) {
            population[*index].objectives[axis] = *value;
            if let Some(exact) = value {
                rationer.absorb(Observation::from_exact(known[*index].clone(), *exact));
            }
        }
    }

    protocol.advance()?;
    Ok(())
}

/// How many objectives an individual carries, in-core plus out-of-core.
///
/// `Selection::for_objectives` reads this rather than the evaluator's own
/// count: an L0 term is a real dimension of the front, and a run that switched
/// from crowding distance to reference directions only when the in-core count
/// crossed the threshold would pick the wrong rule for a job whose fourth axis
/// happens to be the expensive one.
fn objective_width(request: &Tier3Request<'_>) -> usize {
    match request.l0 {
        Some(scorer) => request.evaluator.objective_count().max(scorer.axis() + 1),
        None => request.evaluator.objective_count(),
    }
}
