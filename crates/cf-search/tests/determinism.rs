//! TST-SRCH-004 and GATE-DET-001: the determinism harness.
//!
//! Ticket P1-33.
//!
//! The harness would be worthless without something parallel to point it at,
//! so `P1-33` made tier 2's same-colour block sweep actually run on `rayon`
//! first. `p1_33_the_harness_catches_a_scheduling_dependent_run` is the proof
//! that pointing it at that sweep means something: it builds a deliberately
//! order-dependent computation and requires the harness to reject it.

use cf_constraints::spec::{Constraint, ConstraintClass, ConstraintSet, Rule, Scope, Strand};
use cf_constraints::{ConstraintChecker, KmerMasking};
use cf_landscape_api::provider::TermRegistry;
use cf_landscape_api::HostContext;
use cf_search::determinism::{
    compare_splits, compare_thread_counts, ordered_sum, standard_splits, with_threads,
    GATE_THREAD_COUNTS,
};
use cf_search::job::{run_segment, JobCheckpoint, SegmentRequest};
use cf_search::tier2::{tier2_block_ca, Tier2Config, Tier2Plan, Tier2Request};
use cf_search::{CodonWeights, Genotype, OperatorMix};
use cf_types::{CodonSequence, Digest, HostId, TermId};

mod support;
use support::{Additive, Calls};

const CDS: &str = "ATGCTGCGTAGCGGCGTGGCGCCGACCAAACTGCGTAGCGGCGTGGCGCCGACCAAACTGCGTAGCGGCGTGGCGTAA";
const ROOT_SEED: u64 = 20260828;

fn checker() -> ConstraintChecker {
    ConstraintChecker::compile(&ConstraintSet {
        schema_version: Some("1.0.0".into()),
        constraint_set_id: "cst_det".into(),
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

fn base() -> Genotype {
    Genotype::new(
        CodonSequence::from_cds(CDS, 1)
            .unwrap()
            .with_five_prime_flank("TTGACAATTAATCATCGGCTCGTATAATGTGTGGA")
            .unwrap()
            .with_three_prime_flank("AAAGCCCGGATCCTCTAGAGGG")
            .unwrap(),
    )
}

fn host() -> HostContext {
    HostContext::new(
        HostId::new("escherichia_coli:511145/BL21").unwrap(),
        Digest::of_bytes(b"det-fixture"),
    )
}

/// Tier 2 over enough blocks that a thread pool has something to schedule.
fn tier2_run() -> String {
    let base = base();
    let checker = checker();
    let host = host();
    let mut registry = TermRegistry::new();
    registry
        .register_default(Box::new(Additive::new("usage.cai", Calls::default(), 7)))
        .unwrap();

    let config = Tier2Config::from_planning(3);
    let terms = vec![TermId::new("usage.cai").unwrap()];
    let plan = Tier2Plan::check(&registry, &terms, &[1.0], &config).unwrap();

    let report = tier2_block_ca(
        &Tier2Request {
            seed: &base,
            registry: &registry,
            host: &host,
            checker: &checker,
            config: &config,
        },
        &plan,
        ROOT_SEED,
    )
    .unwrap();

    // The whole result, not a summary. A harness that compared a scalar could
    // agree on the scalar while the sequences differed.
    format!(
        "{}|{:.17}|{}|{}|{:?}",
        report.incumbent.sequence().cds_string(),
        report.objective,
        report.blocks_committed,
        report.blocks_rejected,
        report.surrogate_trace,
    )
}

/// `TST-SRCH-004` / `GATE-DET-001`: 1, 4 and 16 threads agree.
#[test]
fn tst_srch_004_tier_2_is_bit_identical_across_thread_counts() {
    // The fixture must give the pool something to do, or 16 threads and 1
    // thread are the same execution and the comparison is vacuous.
    let blocks = (CDS.len() / 3).div_ceil(3);
    assert!(
        blocks >= 8,
        "only {blocks} blocks; a thread-count comparison needs enough same-colour \
         blocks that scheduling can actually differ"
    );

    compare_thread_counts(&GATE_THREAD_COUNTS, tier2_run)
        .expect("tier 2's answer must not depend on the size of the thread pool");
}

/// The harness rejects output that varies with the pool size.
///
/// Without this, `tst_srch_004` proves only that the harness runs. The
/// mechanism is tested directly rather than through a deliberate data race: a
/// racy control would fail only *probably*, and a control that usually fires is
/// not evidence, it is a flake waiting to be re-run until it passes.
#[test]
fn p1_33_the_harness_rejects_output_that_varies_with_the_pool() {
    let pool_dependent = || rayon::current_num_threads();

    let error = compare_thread_counts(&GATE_THREAD_COUNTS, pool_dependent).expect_err(
        "the harness accepted output that differs between thread counts; it would \
         accept the same defect in tier 2",
    );
    assert_eq!(error.code, "CF-SRCH-013");
    assert!(error.detail.contains("GATE-DET-001"), "{}", error.detail);
}

/// Summation order changes the answer, so a completion-ordered reduction would.
///
/// The arithmetic half of the same hazard. Tier 2 avoids it by collecting
/// proposals into an index-ordered vector and committing in block order, so no
/// sum ever follows completion order — but that is only worth doing if the
/// order is capable of changing the result, which is what this pins.
///
/// The magnitudes matter and the obvious choice does not work: a first version
/// of this test alternated `1e16` and `1.0`, where `1e16 + 1.0 == 1e16`
/// exactly. Every ordering absorbed the small values identically and the test
/// passed for a reason unrelated to ordering. `1e16, 1, -1e16, 1` cancels
/// instead of absorbing, so the sum genuinely depends on where the
/// cancellation falls.
#[test]
fn p1_33_a_fixed_order_reduction_is_stable_where_the_order_matters() {
    let mut values = Vec::new();
    for _ in 0..512 {
        values.extend_from_slice(&[1e16, 1.0, -1e16, 1.0]);
    }
    let reversed: Vec<f64> = values.iter().copied().rev().collect();

    assert_ne!(
        ordered_sum(&values),
        ordered_sum(&reversed),
        "these magnitudes are not order-sensitive, so this proves nothing"
    );

    // The same values in the same order give the same answer at any pool size:
    // a fixed order is what makes the reduction thread-count independent.
    let sum_at =
        |threads: usize| with_threads(threads, || format!("{:.17}", ordered_sum(&values))).unwrap();
    assert_eq!(sum_at(1), sum_at(4));
    assert_eq!(sum_at(4), sum_at(16));
}

/// `P1-29`'s split fixture, reused through the harness rather than rebuilt.
#[test]
fn p1_33_segment_splits_agree_with_an_uninterrupted_run() {
    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let request = SegmentRequest {
        checker: &checker,
        mix: &mix,
        weights: &weights,
        leak: 0.0,
        region: 0,
        masking: KmerMasking::default(),
        root_seed: ROOT_SEED,
    };

    compare_splits(
        12,
        &standard_splits(),
        || JobCheckpoint::start(&base()),
        |state, n| run_segment(&state, &request, n).map(|report| report.checkpoint),
        |state| state.to_json().unwrap(),
    )
    .expect("a resumed run must equal an uninterrupted one of the same seed plan");
}

/// The harness rejects a split that does diverge.
#[test]
fn p1_33_the_harness_catches_a_split_that_diverges() {
    // `advance` that restarts its stream counter on each call — the defect
    // `P1-29`'s absolute generation index exists to prevent.
    let outcome = compare_splits(
        12,
        &standard_splits(),
        || 0u64,
        |state, n| Ok(state + u64::from(n) * u64::from(n)),
        |state| *state,
    );
    let error = outcome.expect_err("a non-additive advance must not look transparent");
    assert_eq!(error.code, "CF-SRCH-013");
    assert!(error.detail.contains("resuming is not transparent"));
}

/// A split that is not the same total is refused, not silently compared.
#[test]
fn p1_33_a_split_that_is_not_the_same_run_is_refused() {
    let outcome = compare_splits(
        12,
        &[vec![6, 5]],
        || 0u64,
        |state, n| Ok(state + u64::from(n)),
        |state| *state,
    );
    let error = outcome.expect_err("11 steps is not 12 steps");
    assert!(
        error.detail.contains("not the same run"),
        "{}",
        error.detail
    );
}

/// Thread count does not change how often the landscape is consulted.
///
/// A pool that duplicated work would still agree on the answer, so equality
/// alone cannot see it — and the cost budget in §7 is stated per half-sweep,
/// not per thread.
#[test]
fn p1_33_thread_count_does_not_change_the_landscape_call_count() {
    let counts: Vec<(u64, u64, u64)> = GATE_THREAD_COUNTS
        .iter()
        .map(|&threads| {
            with_threads(threads, || {
                let base = base();
                let checker = checker();
                let host = host();
                let calls = Calls::default();
                let mut registry = TermRegistry::new();
                registry
                    .register_default(Box::new(Additive::new("usage.cai", calls.clone(), 7)))
                    .unwrap();

                let config = Tier2Config::from_planning(3);
                let terms = vec![TermId::new("usage.cai").unwrap()];
                let plan = Tier2Plan::check(&registry, &terms, &[1.0], &config).unwrap();
                tier2_block_ca(
                    &Tier2Request {
                        seed: &base,
                        registry: &registry,
                        host: &host,
                        checker: &checker,
                        config: &config,
                    },
                    &plan,
                    ROOT_SEED,
                )
                .unwrap();
                (calls.scores(), calls.tables(), calls.rebases())
            })
            .unwrap()
        })
        .collect();

    assert_eq!(
        counts[0], counts[1],
        "1 and 4 threads consulted the landscape a different number of times"
    );
    assert_eq!(counts[1], counts[2]);
    assert!(counts[0].1 > 0, "the fixture never requested a delta table");
}
