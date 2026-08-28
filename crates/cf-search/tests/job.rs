//! The segment runner: resume determinism and checkpoint round-tripping.
//!
//! Ticket P1-29, Rust half.
//!
//! These exist because a negative control did not fire. The resume-determinism
//! claim was tested only from Python, so injecting a stream-counter reset into
//! `cf-search` broke nothing that `cargo test -p cf-search` could see. A defect
//! in this crate should be catchable by this crate's tests.

use cf_constraints::spec::{Constraint, ConstraintClass, ConstraintSet, Rule, Scope, Strand};
use cf_constraints::{ConstraintChecker, KmerMasking, TabuLedger};
use cf_search::job::{run_segment, JobCheckpoint, SegmentRequest};
use cf_search::{CodonWeights, Genotype, OperatorMix};
use cf_types::CodonSequence;

const CDS: &str = "ATGCTGCGTAGCGGCGTGGCGCCGACCAAACTGCGTAGCTAA";
const ROOT_SEED: u64 = 20260828;

fn checker() -> ConstraintChecker {
    // Both strands: on the sense strand alone this motif never arises here, so
    // repair never fires and every claim about the tabu ledger is vacuous.
    ConstraintChecker::compile(&ConstraintSet {
        schema_version: Some("1.0.0".into()),
        constraint_set_id: "cst_job".into(),
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

fn request<'a>(
    checker: &'a ConstraintChecker,
    mix: &'a OperatorMix,
    weights: &'a CodonWeights,
) -> SegmentRequest<'a> {
    SegmentRequest {
        checker,
        mix,
        weights,
        leak: 0.0,
        region: 0,
        masking: KmerMasking::default(),
        root_seed: ROOT_SEED,
    }
}

fn run(generations: &[u32]) -> JobCheckpoint {
    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let request = request(&checker, &mix, &weights);

    let mut checkpoint = JobCheckpoint::start(&base());
    for count in generations {
        checkpoint = run_segment(&checkpoint, &request, *count)
            .unwrap()
            .checkpoint;
    }
    checkpoint
}

fn as_json(checkpoint: &JobCheckpoint) -> serde_json::Value {
    checkpoint.to_json().unwrap()
}

/// Twelve generations at once, or in pieces, is the same twelve generations.
///
/// Every split of the same total, because a run that resumed correctly at one
/// boundary and not another would pass a single-split test.
#[test]
fn p1_29_a_segmented_run_equals_an_uninterrupted_one() {
    let uninterrupted = run(&[12]);

    // The fixture has to exercise repair, or the ledger is empty and its
    // travelling in the checkpoint is untested.
    assert!(
        uninterrupted.repairs_attempted > 0,
        "repair never fired, so this proves nothing about resume"
    );

    for split in [
        vec![6, 6],
        vec![4, 4, 4],
        vec![3, 3, 3, 3],
        vec![1, 11],
        vec![11, 1],
        vec![1; 12],
        vec![5, 2, 5],
    ] {
        let pieces = run(&split);
        assert_eq!(
            as_json(&pieces),
            as_json(&uninterrupted),
            "a run split {split:?} diverged from an uninterrupted run of the same seed plan"
        );
    }
}

/// The ledger must cross the checkpoint, and the fixture must prove it matters.
#[test]
fn p1_29_the_tabu_ledger_travels_in_the_checkpoint() {
    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let request = request(&checker, &mix, &weights);

    // Find a boundary where the ledger is actually live.
    let mut checkpoint = JobCheckpoint::start(&base());
    let mut live_at = None;
    for generation in 1..=12u32 {
        checkpoint = run_segment(&checkpoint, &request, 1).unwrap().checkpoint;
        if !checkpoint.tabu.is_empty() && live_at.is_none() {
            live_at = Some(generation);
        }
    }
    let live_at = live_at.expect("no generation left a live tabu entry; fixture is inadequate");

    // Re-run to that boundary and confirm the entry survives serialization.
    let stopped = run(&[live_at]);
    assert!(!stopped.tabu.is_empty());
    let round_tripped = JobCheckpoint::from_json(&as_json(&stopped)).unwrap();
    assert_eq!(round_tripped.tabu, stopped.tabu);
    assert_eq!(as_json(&round_tripped), as_json(&stopped));
}

/// A zero-generation segment is a no-op that still reports warnings.
#[test]
fn p1_29_an_empty_segment_changes_nothing_but_still_warns() {
    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let start = JobCheckpoint::start(&base());

    let report = run_segment(&start, &request(&checker, &mix, &weights), 0).unwrap();
    assert_eq!(report.generations_run, 0);
    assert_eq!(as_json(&report.checkpoint), as_json(&start));
    assert_eq!(
        report.warnings.len(),
        2,
        "the operator warnings are properties of the mix, not of doing work"
    );
}

/// A malformed checkpoint is refused rather than partly believed.
#[test]
fn p1_29_a_malformed_checkpoint_is_refused() {
    let good = as_json(&run(&[2]));
    for missing in ["generation", "incumbent", "tabu", "history"] {
        let mut broken = good.clone();
        broken.as_object_mut().unwrap().remove(missing);
        let error = JobCheckpoint::from_json(&broken)
            .expect_err("a checkpoint missing a required field must not load");
        assert_eq!(error.code, "CF-ORCH-001");
        assert!(error.detail.contains(missing), "{}", error.detail);
    }
}

/// The ledger snapshot is ordered, so a checkpoint is a function of the state.
///
/// `HashMap` iteration order is randomly seeded per process in Rust, so an
/// unsorted snapshot would make two runs from identical state serialize
/// differently — and would surface as an intermittent resume-comparison
/// failure, which is the worst way to find out. Asserted directly rather than
/// left to chance: entries inserted in opposite orders must snapshot the same.
#[test]
fn p1_29_the_tabu_snapshot_is_ordered() {
    let entries = [
        (0usize, 7usize, cf_types::Codon::parse("CTG").unwrap(), 3u32),
        (0, 2, cf_types::Codon::parse("GCG").unwrap(), 5),
        (0, 2, cf_types::Codon::parse("AAA").unwrap(), 1),
        (1, 0, cf_types::Codon::parse("TTT").unwrap(), 2),
    ];

    let mut forward = TabuLedger::new(0.0);
    for (region, index, codon, remaining) in entries {
        forward.record(region, index, codon, remaining);
    }
    let mut backward = TabuLedger::new(0.0);
    for (region, index, codon, remaining) in entries.iter().rev() {
        backward.record(*region, *index, *codon, *remaining);
    }

    let a: Vec<_> = forward.entries().collect();
    let b: Vec<_> = backward.entries().collect();
    assert_eq!(a, b, "insertion order changed the snapshot");
    assert_eq!(a.len(), entries.len());

    let mut sorted = a.clone();
    sorted.sort_by_key(|(region, index, codon, _)| (*region, *index, codon.index()));
    assert_eq!(a, sorted, "the snapshot is not in a defined order");
}
