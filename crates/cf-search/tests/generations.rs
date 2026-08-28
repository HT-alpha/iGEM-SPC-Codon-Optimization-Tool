//! `TST-CONS-005` — anti-thrashing across generations, and the sub-seed
//! determinism check `GATE-DET-001` will later depend on.
//!
//! Ticket P1-18b.
//!
//! The point of the ledger is not that `is_tabu` answers correctly — `P1-18a`
//! established that in isolation. It is that a move repair reverted stays out
//! of the operator's reach for exactly `tabu_tenure_generations`, which is a
//! property of the *loop*, invisible to a unit test of either half, and
//! sensitive to where the tick sits by exactly one generation in each
//! direction.

mod support;

use cf_constraints::spec::{
    Constraint, ConstraintClass, ConstraintSet, RepairPolicy, Rule, Scope, Strand,
};
use cf_constraints::ConstraintChecker;
use cf_search::generation::{run_generations, GenerationConfig, GenerationRequest};
use cf_search::{CodonWeights, Genotype, OperatorMix, SearchRng};
use cf_types::CodonSequence;

use support::weights;

/// `M K W *`.
///
/// Exactly one mutable position — lysine at codon 1, with two synonyms — and
/// two single-codon families around it. That is what makes the cycle *scripted*
/// rather than merely likely: the operator has one position and one alternative
/// to propose, so every generation's proposal is determined, and the only thing
/// that can vary is whether the ledger suppresses it.
const CDS: &str = "ATGAAATGGTAA";
const FIVE_PRIME: &str = "GGAGGACATATG";
const THREE_PRIME: &str = "CATCACCATCAC";

fn base() -> Genotype {
    Genotype::new(
        CodonSequence::from_cds(CDS, 1)
            .unwrap()
            .with_five_prime_flank(FIVE_PRIME)
            .unwrap()
            .with_three_prime_flank(THREE_PRIME)
            .unwrap(),
    )
}

/// Forbids `AAGTGG`, which exists only when lysine is `AAG`.
///
/// So the single available mutation is exactly the one that creates a
/// violation, and repair's only available fix is to undo it. A thrash by
/// construction.
fn thrashing_constraints(tenure: u32) -> ConstraintSet {
    ConstraintSet {
        schema_version: Some("1.0.0".into()),
        constraint_set_id: "cst_thrash".into(),
        name: None,
        derived_from: vec![],
        constraints: vec![Constraint {
            id: "no_aagtgg".into(),
            class: ConstraintClass::Blocking,
            scope: Scope::Construct,
            cistron_id: None,
            weight: 1.0,
            enabled: true,
            rule: Rule::ForbiddenMotif {
                motif: "AAGTGG".into(),
                strand: Strand::Sense,
                max_occurrences: 0,
            },
        }],
        repair_policy: RepairPolicy {
            tabu_tenure_generations: tenure,
            ..RepairPolicy::default()
        },
        digest: None,
    }
}

fn config(leak: f64) -> GenerationConfig {
    GenerationConfig {
        generations: 20,
        tabu_leak: leak,
        region: 0,
        masking: cf_constraints::KmerMasking::default(),
        first_generation: 0,
    }
}

fn request<'a>(
    base: &'a Genotype,
    checker: &'a ConstraintChecker,
    mix: &'a OperatorMix,
    weights: &'a CodonWeights,
    config: &'a GenerationConfig,
    root_seed: u64,
) -> GenerationRequest<'a> {
    GenerationRequest {
        start: base,
        checker,
        mix,
        weights,
        config,
        root_seed,
    }
}

// --------------------------------------------------------------------------- //
// TST-CONS-005
// --------------------------------------------------------------------------- //

/// A scripted thrash is broken, and it recurs on exactly the tenure period.
///
/// Both halves matter. "Broken" alone would be satisfied by a ledger that
/// suppressed the move forever, which is a different bug — it fences off a
/// region of the space permanently. The period is what shows the tick is
/// placed correctly: one generation early or late and the recurrence lands on
/// a different set of generations, which is precisely what the negative control
/// injects.
#[test]
fn tst_cons_005_a_scripted_thrash_is_broken_within_the_tenure() {
    let tenure = 4u32;
    let checker = ConstraintChecker::compile(&thrashing_constraints(tenure)).unwrap();
    assert_eq!(GenerationConfig::tabu_tenure(&checker), tenure);

    let base = base();
    assert_eq!(
        base.mutable_positions(),
        vec![1],
        "the fixture must have exactly one mutable position, or the script is not scripted"
    );

    let mix = OperatorMix::without_context_model();
    let weights = weights();
    // Leak zero: the leak is a separate mechanism with its own unit test, and
    // leaving it on would make the recurrence pattern probabilistic.
    let report = run_generations(
        &request(&base, &checker, &mix, &weights, &config(0.0), 20260826),
        None,
    )
    .unwrap();

    let proposing = report.generations_proposing(1);
    assert!(
        proposing.contains(&0),
        "generation 0 must propose the thrashing move; history: {:?}",
        report.history
    );
    assert!(
        report.history[0].needed_repair,
        "the proposal was supposed to create a violation"
    );

    // The recurrence is periodic with period `tenure`, starting at 0.
    let expected: Vec<u32> = (0..config(0.0).generations)
        .step_by(tenure as usize)
        .collect();
    assert_eq!(
        proposing, expected,
        "the suppressed move did not recur on the tenure period. history: {:?}",
        report.history
    );

    // Which means it was suppressed for `tenure - 1` consecutive generations
    // after each revert — broken, and not broken forever.
    let gaps: Vec<u32> = proposing.windows(2).map(|w| w[1] - w[0]).collect();
    assert!(
        gaps.iter().all(|g| *g == tenure),
        "gaps between recurrences were {gaps:?}, expected all {tenure}"
    );
    assert!(report.repairs_attempted >= 4);
    assert_eq!(report.repairs_succeeded, report.repairs_attempted);
}

/// A different tenure moves the period, which is what makes the period a
/// measurement rather than a coincidence.
#[test]
fn tst_cons_005_the_period_follows_the_configured_tenure() {
    let base = base();
    let mix = OperatorMix::without_context_model();
    let weights = weights();

    for tenure in [2u32, 3, 5, 6] {
        let checker = ConstraintChecker::compile(&thrashing_constraints(tenure)).unwrap();
        let report = run_generations(
            &request(&base, &checker, &mix, &weights, &config(0.0), 4242),
            None,
        )
        .unwrap();
        let proposing = report.generations_proposing(1);
        let expected: Vec<u32> = (0..20).step_by(tenure as usize).collect();
        assert_eq!(
            proposing, expected,
            "tenure {tenure} produced {proposing:?}"
        );
    }
}

/// The individual ends feasible in every generation.
///
/// Repair either clears the violation or the proposal is abandoned; a loop that
/// carried an infeasible incumbent forward would let the same violation drive
/// the next generation's `violation_directed` draw — restarting the thrash from
/// the loop rather than from the operator.
#[test]
fn tst_cons_005_the_incumbent_is_never_left_infeasible() {
    let checker = ConstraintChecker::compile(&thrashing_constraints(4)).unwrap();
    let report = run_generations(
        &request(
            &base(),
            &checker,
            &OperatorMix::without_context_model(),
            &weights(),
            &config(0.0),
            7,
        ),
        None,
    )
    .unwrap();
    assert!(report.history.iter().all(|record| record.feasible));
}

/// The leak lets a suppressed move through sometimes, and it is off by default
/// in this fixture so the period test above stays deterministic.
#[test]
fn tst_cons_005_the_leak_admits_suppressed_moves() {
    let checker = ConstraintChecker::compile(&thrashing_constraints(8)).unwrap();
    let base = base();
    let mix = OperatorMix::without_context_model();
    let weights = weights();

    let sealed = run_generations(
        &request(&base, &checker, &mix, &weights, &config(0.0), 99),
        None,
    )
    .unwrap();
    // Not 1.0: `TabuLedger::new` refuses that, because a leak of 1 disables the
    // ledger rather than loosening it, and a config that silently disabled the
    // anti-thrashing mechanism is worse than one that fails.
    let leaky = run_generations(
        &request(&base, &checker, &mix, &weights, &config(0.99), 99),
        None,
    )
    .unwrap();

    let sealed_count = sealed.generations_proposing(1).len();
    let leaky_count = leaky.generations_proposing(1).len();
    assert_eq!(
        sealed_count, 3,
        "tenure 8 over 20 generations should recur at 0, 8 and 16"
    );
    assert!(
        leaky_count >= 15,
        "a 0.99 leak suppressed {} of 20 generations; it should suppress about one",
        20 - leaky_count
    );
    assert!(
        leaky_count > sealed_count * 3,
        "the leak made no difference ({leaky_count} against {sealed_count}), so it is \
         not being consulted"
    );
}

/// The loop collects the operator-availability warnings `P1-23` produced.
///
/// It cannot deliver them — the destination is `job_status.warnings`, which is
/// `P1-29`'s document — but it is the first execution context with a lifetime,
/// so the plumbing exists here rather than being invented later.
#[test]
fn p1_18b_the_loop_carries_the_operator_warnings_forward() {
    let checker = ConstraintChecker::compile(&thrashing_constraints(4)).unwrap();
    let report = run_generations(
        &request(
            &base(),
            &checker,
            &OperatorMix::without_context_model(),
            &weights(),
            &config(0.0),
            1,
        ),
        None,
    )
    .unwrap();

    assert_eq!(report.warnings.len(), 2);
    for warning in &report.warnings {
        assert_eq!(warning.code, "CF-SRCH-001");
        assert!(warning.detail.contains("P3-05"));
    }
}

// --------------------------------------------------------------------------- //
// determinism — brought forward from P1-33
// --------------------------------------------------------------------------- //

/// Per-generation sub-seeds do not collide, over a grid.
///
/// `SearchRng::new`'s `seed | 1` defect collapsed every even seed onto its odd
/// successor and was found only because a test happened to compare two seeds.
/// The loop derives a stream per generation, so it introduces a second place
/// the same class of bug could appear. `GATE-DET-001` and `P1-33` are where
/// this would otherwise first be tested, which is far too late — a run whose
/// generations shared a stream would explore a fraction of the space while
/// reporting a full one.
#[test]
fn p1_18b_generation_sub_seeds_do_not_collide() {
    let mut seen = std::collections::HashMap::new();
    for root in 0..64u64 {
        for generation in 0..64u64 {
            let draw = SearchRng::stream(root, generation).next_u64();
            if let Some(previous) = seen.insert(draw, (root, generation)) {
                panic!(
                    "streams ({root}, {generation}) and {previous:?} produce the same first \
                     draw; a run's generations would share randomness"
                );
            }
        }
    }
    assert_eq!(seen.len(), 64 * 64);
}

/// A whole run is reproducible from its root seed, and different roots diverge.
#[test]
fn p1_18b_a_run_is_reproducible_from_its_root_seed() {
    let checker = ConstraintChecker::compile(&thrashing_constraints(4)).unwrap();
    let base = base();
    let mix = OperatorMix::without_context_model();
    let weights = weights();

    let run = |seed: u64| {
        let report = run_generations(
            &request(&base, &checker, &mix, &weights, &config(0.5), seed),
            None,
        )
        .unwrap();
        (report.incumbent.hash(), report.history)
    };

    assert_eq!(run(1234), run(1234));

    // A leak of 0.5 makes the run seed-sensitive; without that the fixture is
    // deterministic regardless of seed and the divergence assertion would be
    // vacuous.
    let mut distinct = std::collections::BTreeSet::new();
    for seed in 0..32u64 {
        distinct.insert(format!("{:?}", run(seed).1));
    }
    assert!(
        distinct.len() > 4,
        "only {} distinct histories across 32 root seeds; the sub-seeding is not \
         reaching the operators",
        distinct.len()
    );
}
