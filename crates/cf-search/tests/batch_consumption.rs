//! `P2-05a` — tier 3 consumes the batch protocol.
//!
//! `TST-SRCH-023`. The ticket joins three pieces that were built separately and
//! never connected: `P2-06`'s ranker, `P2-05`'s protocol, and the generation
//! loop. Each was tested against itself; nothing had ever run the three
//! together, which is precisely the gap this exists to close.
//!
//! ## The sixth fixture-design rule, applied to a crossing
//!
//! Agreement is necessary and not sufficient here in a specific way: a scorer
//! that ignores the genotypes handed to it and returns cached values in the
//! right *shape* would satisfy any check comparing "batch result" against
//! "per-item result", because both would be the same wrong numbers. So the
//! scorer below **records exactly which genotypes it was given on every call**,
//! and the tests assert against that record rather than only against the
//! objective values that come back.
//!
//! ## What the counter has to be, and why
//!
//! `ADR-0001` caps crossings at one per generation. A counter living inside the
//! Rust loop would count the loop's own bookkeeping and could read zero while
//! the boundary was crossed sixty-four times. The counter therefore lives in
//! the *scorer* — the only thing that actually gets called — and counts its own
//! invocations. That is why `CallLog` is shared state rather than a return
//! value.

use std::cell::RefCell;

use cf_constraints::spec::{Constraint, ConstraintClass, ConstraintSet, Rule, Scope, Strand};
use cf_constraints::{ConstraintChecker, KmerMasking};
use cf_search::tier3::{
    run_generations, seed_population, L0Scorer, ObjectiveEvaluator, Tier3Config, Tier3Request,
};
use cf_search::{CodonWeights, Genotype, OperatorMix};
use cf_types::{CodonSequence, Result};

const CDS: &str = "ATGCTGCGTAGCGGCGTGGCGCCGACCAAACTGCGTAGCTAA";

fn checker_for() -> ConstraintChecker {
    ConstraintChecker::compile(&ConstraintSet {
        schema_version: Some("1.0.0".into()),
        constraint_set_id: "cst_p2_05a".into(),
        name: None,
        derived_from: vec![],
        constraints: vec![Constraint {
            id: "no_cgcagc".into(),
            class: ConstraintClass::Blocking,
            scope: Scope::Construct,
            cistron_id: None,
            weight: 1.0,
            enabled: true,
            rule: Rule::ForbiddenMotif {
                motif: "CGCAGC".into(),
                strand: Strand::Both,
                max_occurrences: 0,
            },
        }],
        repair_policy: Default::default(),
        digest: None,
    })
    .unwrap()
}

fn toy_weights() -> CodonWeights {
    CodonWeights::uniform()
}

/// What the scorer was actually asked, on every call.
#[derive(Default)]
struct CallLog {
    /// One entry per call: the sequences handed over, in order.
    calls: Vec<Vec<String>>,
}

impl CallLog {
    fn crossings(&self) -> usize {
        self.calls.len()
    }

    fn all_scored(&self) -> Vec<String> {
        self.calls.iter().flatten().cloned().collect()
    }
}

/// A scorer that records its arguments and derives its answer from them.
///
/// The value returned is a pure function of the sequence, so a scorer that
/// returned stale values for the wrong genotypes produces numbers that do not
/// match the sequences it claims to have scored — which the tests can then
/// detect. A scorer returning a constant would pass a shape check and fail
/// here.
struct RecordingScorer {
    axis: usize,
    budget: usize,
    log: RefCell<CallLog>,
}

impl RecordingScorer {
    fn new(axis: usize, budget: usize) -> Self {
        Self {
            axis,
            budget,
            log: RefCell::new(CallLog::default()),
        }
    }

    /// GC fraction. Cheap, deterministic, and a function of the whole sequence,
    /// so it cannot be produced without having actually seen the genotype.
    fn score_of(sequence: &str) -> f64 {
        let gc = sequence.chars().filter(|c| *c == 'G' || *c == 'C').count();
        gc as f64 / sequence.len().max(1) as f64
    }
}

impl L0Scorer for RecordingScorer {
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
        let values = sequences.iter().map(|s| Some(Self::score_of(s))).collect();
        self.log.borrow_mut().calls.push(sequences);
        Ok(values)
    }
}

/// Two in-core objectives plus an axis the L0 scorer fills.
struct TwoPlusOne;

impl ObjectiveEvaluator for TwoPlusOne {
    fn evaluate(&self, genotype: &Genotype) -> Result<Vec<Option<f64>>> {
        let sequence = genotype.sequence().cds_string();
        let a = sequence.chars().filter(|c| *c == 'A').count() as f64;
        let t = sequence.chars().filter(|c| *c == 'T').count() as f64;
        // Two in-core axes. The L0 axis is *not* declared here: tier 3 appends
        // it, because an evaluator that already filled the slot would be
        // overwritten by the scorer and would lose its value silently.
        Ok(vec![Some(a), Some(t)])
    }

    fn objective_count(&self) -> usize {
        2
    }
}

const L0_AXIS: usize = 2;

fn request<'a>(
    checker: &'a ConstraintChecker,
    mix: &'a OperatorMix,
    weights: &'a CodonWeights,
    evaluator: &'a TwoPlusOne,
    config: &'a Tier3Config,
    scorer: Option<&'a RecordingScorer>,
) -> Tier3Request<'a> {
    Tier3Request {
        checker,
        mix,
        weights,
        evaluator,
        config,
        root_seed: 0x5eed,
        l0: scorer.map(|s| s as &dyn L0Scorer),
    }
}

fn config(population: usize, generations: u32) -> Tier3Config {
    let mut config = Tier3Config::from_planning(population);
    config.max_generations = generations;
    // Long enough that stagnation cannot end the run early and make a
    // crossing-count assertion pass by running fewer generations than asked.
    config.stagnation_window = 1_000;
    config.masking = KmerMasking::Assembled;
    config
}

fn toy_seeds() -> Vec<Genotype> {
    vec![Genotype::new(CodonSequence::from_cds(CDS, 1).expect("cds"))]
}

fn run(
    population: usize,
    generations: u32,
    budget: usize,
) -> (RecordingScorer, cf_search::tier3::Tier3Report) {
    let checker = checker_for();
    let mix = OperatorMix::without_context_model();
    let weights = toy_weights();
    let evaluator = TwoPlusOne;
    let cfg = config(population, generations);
    let scorer = RecordingScorer::new(L0_AXIS, budget);

    let seed_request = request(&checker, &mix, &weights, &evaluator, &cfg, Some(&scorer));
    let start = seed_population(&toy_seeds(), &seed_request).expect("seed");
    let report = run_generations(&start, &seed_request, generations).expect("run");
    (scorer, report)
}

// --------------------------------------------------------------------------- //
// crossing count
// --------------------------------------------------------------------------- //

/// Exactly one crossing per generation. Not "about one", not "at most one".
#[test]
fn tst_srch_023_one_crossing_per_generation_exactly() {
    for generations in [1u32, 3, 7] {
        let (scorer, report) = run(8, generations, 4);
        assert_eq!(
            report.generations_run, generations,
            "the run must actually run the generations asked for, or a crossing \
             count of zero would pass by never generating"
        );
        assert_eq!(
            scorer.log.borrow().crossings(),
            generations as usize,
            "ADR-0001 caps the boundary at one crossing per generation; {} \
             generations produced {} crossings",
            generations,
            scorer.log.borrow().crossings()
        );
    }
}

/// A per-item loop would show up here as N crossings, not one.
#[test]
fn tst_srch_023_the_crossing_count_does_not_scale_with_population() {
    let small = run(8, 4, 4).0.log.borrow().crossings();
    let large = run(24, 4, 4).0.log.borrow().crossings();
    assert_eq!(
        small, large,
        "crossings must be a function of generations alone; scaling with \
         population is the per-item loop this ticket exists to prevent"
    );
}

// --------------------------------------------------------------------------- //
// argument fidelity — the sixth rule
// --------------------------------------------------------------------------- //

/// The scorer is handed real, distinct genotypes — not a placeholder batch.
#[test]
fn tst_srch_023_the_batch_contains_the_genotypes_that_were_ranked() {
    let (scorer, _) = run(12, 3, 5);
    let log = scorer.log.borrow();

    assert!(log.crossings() > 0);
    for call in &log.calls {
        assert!(!call.is_empty(), "an empty batch is a wasted crossing");
        assert!(
            call.len() <= 5,
            "a batch of {} exceeds the E_max of 5, so the budget is not being applied",
            call.len()
        );
        for sequence in call {
            assert_eq!(
                sequence.len() % 3,
                0,
                "a batch member that is not codon-aligned is not a real genotype"
            );
        }
    }

    // Real variety, not the same sequence handed over repeatedly. A batch path
    // that always passed the incumbent would satisfy every count above.
    let distinct: std::collections::BTreeSet<String> = log.all_scored().into_iter().collect();
    assert!(
        distinct.len() > 1,
        "every batched genotype was identical; the batch is not carrying the \
         population it claims to"
    );
}

/// Changing which genotypes are batched changes which scores come back.
///
/// This is the argument-fidelity check. A scorer whose answer did not depend on
/// its input — or a loop that scored the wrong individuals — would produce the
/// same objective values under a different population, and every count-based
/// assertion above would still pass.
#[test]
fn tst_srch_023_different_genotypes_produce_different_scores() {
    let (_scorer, report) = run(12, 3, 5);

    // The recorded scores must be reproducible from the recorded sequences.
    // If the loop attached one sequence's score to another individual, this is
    // where it shows: the value on the axis would not match the genotype.
    let mut checked = 0;
    for individual in &report.checkpoint.members {
        if let Some(value) = individual.objectives[L0_AXIS] {
            let sequence = individual.genotype.sequence().cds_string();
            assert!(
                (value - RecordingScorer::score_of(&sequence)).abs() < 1e-12,
                "an individual carries an L0 score that its own sequence does not \
                 produce; the batch results were attached to the wrong genotypes"
            );
            checked += 1;
        }
    }
    assert!(
        checked > 0,
        "no individual carries an L0 score, so this test verified nothing"
    );
}

// --------------------------------------------------------------------------- //
// plausibility — the sixth rule's second half
// --------------------------------------------------------------------------- //

/// Unbatched individuals stay UNKNOWN. `None` is a value, not a gap to fill.
#[test]
fn tst_srch_023_individuals_outside_the_budget_stay_unknown() {
    let (_, report) = run(16, 2, 3);

    let scored = report
        .checkpoint
        .members
        .iter()
        .filter(|i| i.objectives[L0_AXIS].is_some())
        .count();
    let unknown = report.checkpoint.members.len() - scored;

    assert!(
        unknown > 0,
        "a population of 16 at E_max 3 over 2 generations cannot have every \
         member scored; if it does, something is imputing values"
    );
    assert!(
        scored > 0,
        "and some must be scored, or the batch path never ran"
    );
}

/// A scorer that declines an item leaves it UNKNOWN rather than defaulting.
#[test]
fn tst_srch_023_a_declined_item_is_unknown_not_zero() {
    struct Declining {
        axis: usize,
    }
    impl L0Scorer for Declining {
        fn axis(&self) -> usize {
            self.axis
        }
        fn budget(&self) -> usize {
            4
        }
        fn score_batch(&self, genotypes: &[&Genotype]) -> Result<Vec<Option<f64>>> {
            Ok(vec![None; genotypes.len()])
        }
    }

    let checker = checker_for();
    let mix = OperatorMix::without_context_model();
    let weights = toy_weights();
    let evaluator = TwoPlusOne;
    let cfg = config(8, 2);
    let declining = Declining { axis: L0_AXIS };
    let req = Tier3Request {
        checker: &checker,
        mix: &mix,
        weights: &weights,
        evaluator: &evaluator,
        config: &cfg,
        root_seed: 0x5eed,
        l0: Some(&declining as &dyn L0Scorer),
    };
    let start = seed_population(&toy_seeds(), &req).expect("seed");
    let report = run_generations(&start, &req, 2).expect("run");

    assert!(
        report
            .checkpoint
            .members
            .iter()
            .all(|i| i.objectives[L0_AXIS].is_none()),
        "a scorer that returned UNKNOWN for everything must leave every axis \
         UNKNOWN; a zero here would be a fabricated objective value"
    );
}

// --------------------------------------------------------------------------- //
// rationing and CF-SRCH-016 through the integrated path
// --------------------------------------------------------------------------- //

/// The stranded set is populated by a real run, not by a unit fixture.
///
/// `P2-05`/`P2-06` proved the protocol and the ranker against themselves. This
/// is the first time an individual is stranded by the budget of a search that
/// actually ran, which is what `CF-SRCH-016` claims to report.
#[test]
fn tst_srch_023_the_budget_strands_individuals_and_the_run_reports_them() {
    let (_, report) = run(20, 3, 2);

    assert!(
        !report.stranded.is_empty(),
        "20 individuals at E_max 2 over 3 generations must strand some; if the \
         report is empty the warning can never fire on a real job"
    );
    assert!(
        report.stranded.iter().all(|index| *index < 20 + 20),
        "a stranded index must refer to a real population slot"
    );
}

/// A budget large enough for everyone strands nobody.
///
/// The negative half: without it, a report that always listed something would
/// pass the test above and would warn on every job.
#[test]
fn tst_srch_023_a_sufficient_budget_strands_nobody() {
    let (_, report) = run(6, 3, 64);
    assert!(
        report.stranded.is_empty(),
        "a budget larger than the population stranded {:?}",
        report.stranded
    );
}

/// The protocol's state survives the checkpoint, so a resume keeps the warning.
#[test]
fn tst_srch_023_the_batch_state_travels_in_the_checkpoint() {
    let (_, report) = run(20, 2, 2);
    let batch = report
        .checkpoint
        .batch
        .as_ref()
        .expect("an L0 run must checkpoint its protocol state");

    let stranded = batch
        .get("stranded")
        .and_then(|v| v.as_array())
        .expect("the snapshot must carry the stranded set");
    assert!(
        !stranded.is_empty(),
        "the stranded set was dropped from the checkpoint; CF-SRCH-016 would \
         vanish on resume, which is when a long job most needs it"
    );

    // Through text, per the serialization fixture rule.
    let text = serde_json::to_string(batch).expect("serialize");
    let reparsed: serde_json::Value = serde_json::from_str(&text).expect("parse");
    assert_eq!(
        &reparsed, batch,
        "the snapshot must survive its own transport"
    );
}

/// A run with no L0 term neither crosses nor checkpoints a protocol.
#[test]
fn tst_srch_023_a_run_without_an_l0_term_is_unchanged() {
    let checker = checker_for();
    let mix = OperatorMix::without_context_model();
    let weights = toy_weights();
    let evaluator = TwoPlusOne;
    let cfg = config(8, 3);
    let req = request(&checker, &mix, &weights, &evaluator, &cfg, None);

    let start = seed_population(&toy_seeds(), &req).expect("seed");
    let report = run_generations(&start, &req, 3).expect("run");

    assert!(
        report.checkpoint.batch.is_none(),
        "a run with no out-of-core term should carry no protocol state"
    );
    assert!(report.stranded.is_empty());
}

/// The axis guard has a test, because a guard without one is decoration.
///
/// Added when the negative control for it went silent: the fixture above
/// declares two objectives and an axis of 2, so `axis < objective_count` was
/// never true and the guard was never reached. Removing it changed nothing,
/// which is indistinguishable from it not existing.
#[test]
fn tst_srch_023_a_scorer_may_not_claim_an_axis_the_evaluator_fills() {
    let checker = checker_for();
    let mix = OperatorMix::without_context_model();
    let weights = toy_weights();
    let evaluator = TwoPlusOne; // declares 2 objectives, indices 0 and 1
    let cfg = config(6, 1);

    // Axis 1 is in-core. Scoring it would overwrite the evaluator's value with
    // no trace, and the objective would silently become something else.
    let overlapping = RecordingScorer::new(1, 4);
    let req = request(&checker, &mix, &weights, &evaluator, &cfg, Some(&overlapping));

    let error = seed_population(&toy_seeds(), &req)
        .expect_err("an overlapping axis must be refused");
    assert_eq!(error.code, "CF-SRCH-014");
    assert!(
        error.detail.contains("axis 1"),
        "the refusal must name the axis: {}",
        error.detail
    );

    // And the non-overlapping case still works, so this is a guard rather than
    // a blanket refusal.
    let appended = RecordingScorer::new(2, 4);
    let ok = request(&checker, &mix, &weights, &evaluator, &cfg, Some(&appended));
    assert!(seed_population(&toy_seeds(), &ok).is_ok());
}

/// The L0 axis influences selection, which is only true if it is scored first.
///
/// Also added after a silent control. Scoring after `select` leaves every
/// individual's L0 value `UNKNOWN` at the moment survival is decided, so the
/// term is computed, reported, and completely without effect — a decorative
/// objective that costs a boundary crossing per generation.
///
/// Two scorers with opposite preferences over the same seed and the same RNG:
/// if the axis reaches selection, the surviving populations differ. If it does
/// not, selection sees identical inputs both times and the survivors match.
#[test]
fn tst_srch_023_the_l0_axis_actually_reaches_selection() {
    struct Signed {
        sign: f64,
    }
    impl L0Scorer for Signed {
        fn axis(&self) -> usize {
            L0_AXIS
        }
        fn budget(&self) -> usize {
            64 // large enough that every individual is scored
        }
        fn score_batch(&self, genotypes: &[&Genotype]) -> Result<Vec<Option<f64>>> {
            Ok(genotypes
                .iter()
                .map(|g| Some(self.sign * RecordingScorer::score_of(&g.sequence().cds_string())))
                .collect())
        }
    }

    let survivors = |sign: f64| -> Vec<String> {
        let checker = checker_for();
        let mix = OperatorMix::without_context_model();
        let weights = toy_weights();
        let evaluator = TwoPlusOne;
        let cfg = config(12, 4);
        let scorer = Signed { sign };
        let req = Tier3Request {
            checker: &checker,
            mix: &mix,
            weights: &weights,
            evaluator: &evaluator,
            config: &cfg,
            root_seed: 0x5eed,
            l0: Some(&scorer as &dyn L0Scorer),
        };
        let start = seed_population(&toy_seeds(), &req).expect("seed");
        let report = run_generations(&start, &req, 4).expect("run");
        let mut out: Vec<String> = report
            .checkpoint
            .members
            .iter()
            .map(|i| i.genotype.sequence().cds_string())
            .collect();
        out.sort();
        out
    };

    assert_ne!(
        survivors(1.0),
        survivors(-1.0),
        "reversing the L0 objective's direction did not change who survived, so \
         the axis is UNKNOWN at the moment selection ranks and the term is \
         decorative - computed, reported, and without effect"
    );
}
