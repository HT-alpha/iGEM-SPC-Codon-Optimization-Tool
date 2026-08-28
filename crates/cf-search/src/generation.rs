//! The generation loop, and where `TABU_TICK` goes.
//!
//! Ticket P1-18b. Specification:
//! [`06-constraint-repair.md`](../../../planning/03-algorithms/06-constraint-repair.md)
//! (the anti-thrashing ledger) and
//! [`03-nsga-loop-and-operators.md`](../../../planning/03-algorithms/03-nsga-loop-and-operators.md)
//! §5 (the operators it drives).
//!
//! ## The cycle this exists to break
//!
//! Mutation proposes a move into a violation; repair reverts it; mutation
//! proposes the same move next generation. Nothing is wrong with either
//! component — the loop is the defect, and it is invisible to a unit test of
//! either half. `P1-18a` shipped the ledger that records a reverted move;
//! this closes the loop by giving the ledger a clock.
//!
//! ## Why the tick is at the *end* of a generation
//!
//! An entry recorded with tenure `T` must survive exactly `T` generations. It
//! is recorded during generation `g`, so it has to be live for the whole of
//! `g` and for `T − 1` generations after. Ticking at the end of `g` decrements
//! it once for `g` itself, leaving `T − 1`. Ticking at the *start* of a
//! generation would decrement the entry before it had ever been consulted,
//! expiring it a generation early; ticking twice, or after the mutation step
//! but before repair, shifts the window the other way. Off-by-one in either
//! direction is a real bug with no symptom other than thrashing that resumes
//! slightly sooner or later than it should, which is why `TST-CONS-005` scripts
//! several generations rather than checking one tick.
//!
//! ## What this deliberately does not own
//!
//! Four things a real job needs are `P1-29`'s, and this type takes them as
//! inputs rather than inventing them:
//!
//! * the **termination budget** — here a generation count, in a job a compute
//!   budget with a tier plan;
//! * the **RNG plan** — this derives per-generation streams from one root seed
//!   via [`SearchRng::stream`], which `P1-29` must seed from the lockfile's
//!   `rng.root_seed` rather than from anywhere else;
//! * **checkpoint/resume** — [`GenerationReport`] carries enough to describe a
//!   run but nothing here serializes it;
//! * **cooperative cancel** — the loop runs to its budget.
//!
//! And one obligation inherited from `P1-23`: [`OperatorMix::warnings`] returns
//! `CF-SRCH-001` for every operator the build cannot provide, and the
//! renormalization that follows is a real change in how the search explores.
//! This loop collects them into [`GenerationReport::warnings`] because it is the
//! first execution context with a lifetime — but it cannot discharge the
//! obligation, because the destination is `job_status.warnings` and that
//! document is `P1-29`'s.

use cf_constraints::repair::{repair, NoObjective, ObjectiveDelta};
use cf_constraints::{
    ActionableViolation, ConstraintChecker, ConstructView, GiveUp, KmerMasking, RepairRng,
    RepairTarget, TabuLedger,
};
use cf_types::{CfError, Result};

use crate::context::CodonWeights;
use crate::genotype::Genotype;
use crate::operators::{MutationContext, OperatorMix};
use crate::rng::SearchRng;

/// How the loop runs.
#[derive(Clone, Debug)]
pub struct GenerationConfig {
    pub generations: u32,
    /// Probability that a tabu entry is ignored anyway.
    ///
    /// Keeps the ledger from permanently fencing off a region that later
    /// becomes reachable. Supplied rather than defaulted: it is a search-tuning
    /// decision with no neutral value, and picking one here would bury it.
    pub tabu_leak: f64,
    pub region: usize,
    /// K-mer masking, forwarded to repair. Legacy parity per EXIT-001/EXIT-002;
    /// not this ticket's default to change.
    pub masking: KmerMasking,
    /// Absolute index of the first generation this call runs.
    ///
    /// Non-zero when resuming. Streams are derived from the *absolute* index,
    /// so generation `k` sees the same randomness whether it is the first
    /// generation of a resumed segment or the fiftieth of an uninterrupted run.
    /// A counter that restarted at zero on resume would make a resumed job
    /// diverge from an uninterrupted one while both looked reproducible —
    /// each is reproducible against itself.
    pub first_generation: u32,
}

impl GenerationConfig {
    /// Tenure comes from the constraint set's repair policy, not from here.
    ///
    /// `repair` reads `RepairPolicy::tabu_tenure_generations` when it records a
    /// reverted edit, so a second copy on this struct would be a second source
    /// of truth for the same number — and the loop's tick would be calibrated
    /// against the wrong one the first time they diverged.
    pub fn tabu_tenure(checker: &ConstraintChecker) -> u32 {
        checker.set().repair_policy.tabu_tenure_generations
    }
}

/// What one generation did.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GenerationRecord {
    pub generation: u32,
    /// Codon indices the mutation operator proposed.
    pub proposed: Vec<usize>,
    /// Whether the proposal left the individual infeasible.
    pub needed_repair: bool,
    /// Codon indices repair reverted or changed.
    pub repaired: Vec<usize>,
    /// Entries live in the ledger at the end of the generation, after the tick.
    pub tabu_live: usize,
    /// Whether the individual was feasible when the generation ended.
    pub feasible: bool,
}

/// What a run did.
#[derive(Clone, Debug)]
pub struct GenerationReport {
    pub incumbent: Genotype,
    pub history: Vec<GenerationRecord>,
    /// Warnings a job must surface. See the module note: collected here,
    /// delivered by `P1-29`.
    pub warnings: Vec<CfError>,
    pub repairs_attempted: u64,
    pub repairs_succeeded: u64,
    pub gave_up: Vec<GiveUp>,
    /// Live tabu entries when the run stopped, for a checkpoint to carry.
    ///
    /// The ledger holds state across generations by design, so a resume that
    /// dropped it would re-propose exactly the moves repair had just reverted.
    pub tabu_snapshot: Vec<(usize, usize, cf_types::Codon, u32)>,
}

impl GenerationReport {
    /// Generations in which a given codon index was proposed.
    ///
    /// The observable `TST-CONS-005` is written against: a suppressed move must
    /// be absent for exactly `tabu_tenure` generations, and "exactly" is only
    /// checkable against the full history.
    pub fn generations_proposing(&self, codon_index: usize) -> Vec<u32> {
        self.history
            .iter()
            .filter(|record| record.proposed.contains(&codon_index))
            .map(|record| record.generation)
            .collect()
    }
}

/// Everything one run of the cycle needs.
///
/// A struct rather than eight parameters, matching `Tier1Request` and
/// `Tier2Request`: the `&`-shaped arguments in the middle were positionally
/// interchangeable, and a caller that swapped two of them would find out at
/// runtime if at all.
#[derive(Clone, Copy)]
pub struct GenerationRequest<'a> {
    pub start: &'a Genotype,
    pub checker: &'a ConstraintChecker,
    pub mix: &'a OperatorMix,
    pub weights: &'a CodonWeights,
    pub config: &'a GenerationConfig,
    pub root_seed: u64,
}

impl std::fmt::Debug for GenerationRequest<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GenerationRequest")
            .field("generations", &self.config.generations)
            .field("first_generation", &self.config.first_generation)
            .field("root_seed", &self.root_seed)
            .finish()
    }
}

/// Run the mutation/repair cycle.
///
/// One generation is: propose → apply → check → repair if needed → **tick**.
pub fn run_generations(
    request: &GenerationRequest<'_>,
    resume_from: Option<TabuLedger>,
) -> Result<GenerationReport> {
    run_generations_with_objective(request, &NoObjective, resume_from)
}

/// As [`run_generations`], with an objective for repair to prefer cheap fixes.
///
/// The objective arrives as a trait object rather than as a landscape handle:
/// `cf-constraints` inverts that dependency deliberately (see `ObjectiveDelta`),
/// and the loop has no business re-coupling them.
pub fn run_generations_with_objective(
    request: &GenerationRequest<'_>,
    objective: &dyn ObjectiveDelta,
    resume_from: Option<TabuLedger>,
) -> Result<GenerationReport> {
    let GenerationRequest {
        start,
        checker,
        mix,
        weights,
        config,
        root_seed,
    } = *request;
    let mut incumbent = start.clone();
    let mut ledger = resume_from.unwrap_or_else(|| TabuLedger::new(config.tabu_leak));
    let mut history = Vec::with_capacity(config.generations as usize);
    let mut repairs_attempted = 0u64;
    let mut repairs_succeeded = 0u64;
    let mut gave_up = Vec::new();

    let first = config.first_generation;
    for generation in first..first + config.generations {
        // A stream per generation, derived rather than drawn, so generation `g`
        // sees the same stream whatever ran before it — including whether
        // anything ran before it at all. `SearchRng::stream` goes through
        // `SearchRng::new`, which is injective: the property the `seed | 1`
        // defect violated.
        let mut rng = SearchRng::stream(root_seed, u64::from(generation));

        let construct = assemble(&incumbent);
        let view = ConstructView::new(&construct, vec![cds_span(&incumbent)]);
        let found = checker.check(&view);
        let actionable: Vec<ActionableViolation<'_>> =
            found.iter().filter_map(ActionableViolation::new).collect();

        let ctx = MutationContext::new(&incumbent, weights)
            .with_violations(&actionable)
            .with_tabu(&ledger)
            .with_region(config.region);
        let proposal = mix.propose(&ctx, &mut rng)?.compacted(&incumbent);
        let proposed: Vec<usize> = proposal.touched();

        let mutated = incumbent.with_move(&proposal)?;
        let construct = assemble(&mutated);
        let view = ConstructView::new(&construct, vec![cds_span(&mutated)]);
        let after = checker.check(&view);
        let blocking: Vec<ActionableViolation<'_>> = after
            .iter()
            .filter(|v| v.class == cf_constraints::spec::ConstraintClass::Blocking)
            .filter_map(ActionableViolation::new)
            .collect();

        let needed_repair = !blocking.is_empty();
        let mut repaired = Vec::new();
        let mut next = mutated.clone();

        if needed_repair {
            repairs_attempted += 1;
            let mut target = RepairTarget::new(
                assemble(&mutated),
                vec![(mutated.cds_offset_nt(), mutated.sequence().clone())],
            )?;
            let mut repair_rng = RepairRng::new(rng.next_u64());
            let outcome = repair(
                &mut target,
                checker,
                config.masking,
                objective,
                &mut ledger,
                &mut repair_rng,
            )?;
            repaired = outcome.edits.iter().map(|e| e.codon_index).collect();
            gave_up.extend(outcome.gave_up.iter().cloned());
            if outcome.residual.is_empty() {
                repairs_succeeded += 1;
                next = Genotype::new(target.regions()[0].1.clone());
            } else {
                // Repair could not clear every blocking violation, so the
                // proposal is abandoned and the incumbent stands. Keeping the
                // infeasible individual would let the same violation drive the
                // next generation's `violation_directed` draw — the thrash the
                // ledger exists to break, restarted by the loop itself.
                next = incumbent.clone();
            }
        }

        incumbent = next;

        // END OF GENERATION. Everything above ran with the ledger as it stood
        // when the generation began; the tick belongs here and nowhere else.
        ledger.tick();

        let construct = assemble(&incumbent);
        let view = ConstructView::new(&construct, vec![cds_span(&incumbent)]);
        history.push(GenerationRecord {
            generation,
            proposed,
            needed_repair,
            repaired,
            tabu_live: ledger.len(),
            feasible: checker.is_feasible(&view),
        });
    }

    Ok(GenerationReport {
        incumbent,
        history,
        warnings: mix.warnings(),
        repairs_attempted,
        repairs_succeeded,
        gave_up,
        tabu_snapshot: ledger.entries().collect(),
    })
}

fn assemble(g: &Genotype) -> Vec<u8> {
    let sequence = g.sequence();
    let mut out = sequence.five_prime_flank().to_vec();
    for codon in sequence.codons() {
        out.extend_from_slice(&codon.to_bytes());
    }
    out.extend_from_slice(sequence.three_prime_flank());
    out
}

fn cds_span(g: &Genotype) -> (usize, usize) {
    let start = g.cds_offset_nt();
    (start, start + g.len() * 3)
}
