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

use cf_constraints::repair::{ActionableViolation, RepairTarget};
use cf_constraints::{ConstraintChecker, ConstructView, KmerMasking, TabuLedger};
use cf_metrics::{crowding_distance, fast_nondominated_sort, Ranked};
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

#[derive(Clone, Copy, Debug)]
pub struct Tier3Config {
    /// Population size `N`.
    pub population: usize,
    pub max_generations: u32,
    /// Crossover probability, `p_c` in §7.
    pub crossover_rate: f64,
    pub tabu_leak: f64,
    pub region: usize,
    pub masking: KmerMasking,
}

impl Tier3Config {
    /// `N` and `p_c` from §7's worked loop.
    pub fn from_planning(population: usize) -> Self {
        Self {
            population,
            max_generations: 20,
            crossover_rate: 0.9,
            tabu_leak: 0.0,
            region: 0,
            masking: KmerMasking::default(),
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
}

#[derive(Clone, Debug)]
pub struct Tier3Report {
    pub checkpoint: Tier3Checkpoint,
    /// The non-dominated set, best front first, in a stable order.
    pub front: Vec<Individual>,
    pub generations_run: u32,
    pub warnings: Vec<CfError>,
}

pub struct Tier3Request<'a> {
    pub checker: &'a ConstraintChecker,
    pub mix: &'a OperatorMix,
    pub weights: &'a CodonWeights,
    pub evaluator: &'a dyn ObjectiveEvaluator,
    pub config: &'a Tier3Config,
    /// Already namespaced to tier 3 by the caller. `DEP-DET-003`.
    pub root_seed: u64,
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
    let objectives = request.evaluator.evaluate(&genotype)?;
    if objectives.len() != request.evaluator.objective_count() {
        return Err(CfError::new(
            code::SRCH_TIER3_CONFIG,
            format!(
                "evaluator returned {} objectives but declares {}",
                objectives.len(),
                request.evaluator.objective_count()
            ),
        ));
    }
    let violation = actionable_violation(&genotype, request.checker);
    Ok(Individual {
        genotype,
        objectives,
        violation,
    })
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
    })
}

/// Advance `generations` complete generations. `M ≤ 3` path.
pub fn run_generations(
    checkpoint: &Tier3Checkpoint,
    request: &Tier3Request<'_>,
    generations: u32,
) -> Result<Tier3Report> {
    if request.evaluator.objective_count() > 3 {
        return Err(CfError::new(
            code::SRCH_TIER3_CONFIG,
            format!(
                "tier 3's NSGA-II path handles at most 3 objectives, not {}; \
                 crowding distance degrades badly as M grows and reference-direction \
                 niching (P2-03) is what handles it",
                request.evaluator.objective_count()
            ),
        ));
    }

    let mut ledger = TabuLedger::new(request.config.tabu_leak);
    for (region, index, codon, remaining) in &checkpoint.tabu {
        ledger.record(*region, *index, *codon, *remaining);
    }

    let mut members = checkpoint.members.clone();
    let mut evaluations = checkpoint.evaluations;
    let mut repairs_attempted = checkpoint.repairs_attempted;
    let mut repairs_succeeded = checkpoint.repairs_succeeded;

    let first = checkpoint.generation;
    for generation in first..first + generations {
        // A stream per absolute generation, so generation k sees the same
        // randomness whether it is the first of a resumed segment or the
        // fiftieth of an uninterrupted run.
        let mut rng = SearchRng::stream(request.root_seed, u64::from(generation));

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
                &members[a].genotype,
                &members[b].genotype,
                request,
                &ledger,
                &mut pair_rng,
            )?;
            for child in bred {
                let (repaired, attempted, succeeded) =
                    repair_child(child, request, &mut ledger, &mut pair_rng)?;
                repairs_attempted += attempted;
                repairs_succeeded += succeeded;
                if offspring.len() < request.config.population {
                    offspring.push(evaluate_one(repaired, request)?);
                    evaluations += 1;
                }
            }
        }

        let mut combined = members;
        combined.extend(offspring);
        members = select(combined, request.config.population);
        ledger.tick();
    }

    let front = nondominated(&members);
    Ok(Tier3Report {
        checkpoint: Tier3Checkpoint {
            generation: first + generations,
            members,
            tabu: ledger.entries().collect(),
            evaluations,
            repairs_attempted,
            repairs_succeeded,
        },
        front,
        generations_run: generations,
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

fn repair_child(
    genotype: Genotype,
    request: &Tier3Request<'_>,
    ledger: &mut TabuLedger,
    rng: &mut SearchRng,
) -> Result<(Genotype, u64, u64)> {
    if actionable_violation(&genotype, request.checker) <= 0.0 {
        return Ok((genotype, 0, 0));
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
        // Lamarckian: the genotype is replaced by its repaired form, which is
        // `EXIT-003`'s standing default until `EXP-05` decides otherwise. The
        // Baldwinian alternative evaluates the repaired form while inheriting
        // the original, and choosing between them is what `P1-34` exists for.
        Ok((Genotype::new(target.regions()[0].1.clone()), 1, 1))
    } else {
        // Repair could not clear every blocking violation. The individual is
        // kept, not discarded: §4 is explicit that infeasible individuals
        // survive as stepping stones across a constraint barrier, and
        // constraint-domination is what applies the pressure instead.
        Ok((genotype, 1, 0))
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

/// NSGA-II environmental selection: fronts in order, crowding within the
/// splitting front.
fn select(combined: Vec<Individual>, keep: usize) -> Vec<Individual> {
    let ranked: Vec<Ranked> = combined.iter().map(Individual::ranked).collect();
    let fronts = fast_nondominated_sort(&ranked);

    let mut chosen: Vec<usize> = Vec::with_capacity(keep);
    for front in &fronts {
        if chosen.len() + front.len() <= keep {
            chosen.extend(front.iter().copied());
            continue;
        }
        // The splitting front: take the least crowded, ties by index so the
        // result is a function of the population rather than of the sort's
        // internals.
        let distance = crowding_distance(&ranked, front);
        let mut order: Vec<usize> = (0..front.len()).collect();
        order.sort_by(|&a, &b| {
            distance[b]
                .partial_cmp(&distance[a])
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| front[a].cmp(&front[b]))
        });
        for &position in order.iter().take(keep - chosen.len()) {
            chosen.push(front[position]);
        }
        break;
    }

    chosen.sort_unstable();
    chosen.into_iter().map(|i| combined[i].clone()).collect()
}

/// The mutually non-dominated members, in ascending member order.
pub fn nondominated(members: &[Individual]) -> Vec<Individual> {
    let ranked: Vec<Ranked> = members.iter().map(Individual::ranked).collect();
    match fast_nondominated_sort(&ranked).first() {
        Some(front) => front.iter().map(|&i| members[i].clone()).collect(),
        None => Vec::new(),
    }
}
