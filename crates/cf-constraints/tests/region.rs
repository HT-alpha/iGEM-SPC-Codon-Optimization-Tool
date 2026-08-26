//! `TST-CONS-009` — the bounded-neighbourhood check agrees with the full check.
//!
//! Ticket P1-21 (support for `constraint.violation_count`'s exact delta).
//!
//! `check_region` exists so a landscape term can compute an exact delta in
//! `O(W)` instead of `O(L)`. "Exact" is the entire value of the method: a
//! windowed check that is merely usually right would produce deltas that are
//! usually right, and a search steered by usually-right deltas converges
//! confidently on the wrong answer. So the test is a differential against the
//! whole-construct check over randomized inputs, not a handful of examples.

use cf_constraints::spec::{
    Constraint, ConstraintClass, ConstraintSet, HomopolymerBase, RepairPolicy, Rule, Scope, Strand,
};
use cf_constraints::{ConstraintChecker, ConstructView, Violation};

fn set_with(rules: Vec<(&str, Rule)>) -> ConstraintSet {
    ConstraintSet {
        schema_version: Some("1.0.0".into()),
        constraint_set_id: "cst_region_test".into(),
        name: None,
        derived_from: vec![],
        constraints: rules
            .into_iter()
            .map(|(id, rule)| Constraint {
                id: id.into(),
                class: ConstraintClass::Blocking,
                scope: Scope::Construct,
                cistron_id: None,
                weight: 1.0,
                enabled: true,
                rule,
            })
            .collect(),
        repair_policy: RepairPolicy::default(),
        digest: None,
    }
}

fn local_rules() -> ConstraintSet {
    set_with(vec![
        (
            "no_ecori",
            Rule::ForbiddenMotif {
                motif: "GAATTC".into(),
                strand: Strand::Both,
                max_occurrences: 0,
            },
        ),
        (
            "no_bamhi",
            Rule::RestrictionSite {
                enzyme: "BamHI".into(),
                site: "GGATCC".into(),
                strand: Strand::Both,
                max_occurrences: 0,
            },
        ),
        (
            "runs",
            Rule::Homopolymer {
                base: HomopolymerBase::Any,
                max_run: 4,
            },
        ),
        (
            "gc",
            Rule::GcWindow {
                window_nt: 30,
                min: 0.25,
                max: 0.75,
            },
        ),
    ])
}

/// A deterministic pseudo-random construct generator.
///
/// Biased toward runs and toward the two restriction sites, because a uniform
/// random sequence almost never contains either and the test would then be
/// comparing two empty lists a thousand times.
fn construct(seed: u64, len: usize) -> Vec<u8> {
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut out = Vec::with_capacity(len);
    while out.len() < len {
        match next() % 16 {
            0 => out.extend_from_slice(b"GAATTC"),
            1 => out.extend_from_slice(b"GGATCC"),
            2 => {
                let base = b"ACGT"[(next() % 4) as usize];
                let run = 3 + (next() % 6) as usize;
                out.extend(std::iter::repeat(base).take(run));
            }
            3 => out.extend_from_slice(b"GCGCGCGCGCGC"),
            4 => out.extend_from_slice(b"ATATATATATAT"),
            _ => out.push(b"ACGT"[(next() % 4) as usize]),
        }
    }
    out.truncate(len);
    out
}

fn intersecting(all: &[Violation], start: usize, end: usize) -> Vec<Violation> {
    let mut hit: Vec<Violation> = all
        .iter()
        .filter(|v| v.start_nt < end && start < v.end_nt)
        .cloned()
        .collect();
    hit.sort_by(|a, b| {
        (a.start_nt, a.end_nt, &a.constraint_id).cmp(&(b.start_nt, b.end_nt, &b.constraint_id))
    });
    hit
}

#[test]
fn tst_cons_009_region_check_equals_the_full_check_restricted() {
    let checker = ConstraintChecker::compile(&local_rules()).unwrap();
    assert_eq!(
        checker.influence_radius_nt(),
        Some(30),
        "the radius is the widest local reach: the 30 nt GC window"
    );

    let mut compared = 0usize;
    let mut nonempty = 0usize;

    for seed in 1..=60u64 {
        let dna = construct(seed, 240);
        let view = ConstructView::new(&dna, vec![(30, 210)]);
        let all = checker.check(&view);

        for start in (0..dna.len()).step_by(3) {
            let end = (start + 3).min(dna.len());
            let regional = checker
                .check_region(&view, start, end)
                .expect("this set is local, so a neighbourhood exists");
            let expected = intersecting(&all, start, end);
            assert_eq!(
                regional, expected,
                "seed {seed}, region [{start},{end}) disagreed with the full check"
            );
            compared += 1;
            if !expected.is_empty() {
                nonempty += 1;
            }
        }
    }

    // Without this the test could pass by comparing empty lists forever.
    assert!(
        nonempty * 10 > compared,
        "fixture is too clean to exercise the path: only {nonempty} of {compared} \
         regions contained a violation"
    );
}

/// A k-mer repeat rule makes the reach unbounded, and that is reported rather
/// than approximated.
#[test]
fn tst_cons_009_kmer_repeat_has_no_bounded_neighbourhood() {
    let with_kmer = set_with(vec![(
        "repeats",
        Rule::KmerRepeat {
            k: 8,
            max_occurrences: 2,
            count_overlapping: true,
            canonical_strand: true,
        },
    )]);
    let checker = ConstraintChecker::compile(&with_kmer).unwrap();
    assert_eq!(checker.influence_radius_nt(), None);

    let dna = construct(7, 120);
    let view = ConstructView::new(&dna, vec![(0, 120)]);
    assert!(
        checker.check_region(&view, 30, 33).is_none(),
        "there is no neighbourhood that bounds a global occurrence count"
    );
}

/// A homopolymer run longer than the fixed radius is still seen whole.
///
/// The neighbourhood is grown to the run's own ends. A version that offset by
/// the radius and stopped would report the run truncated — or, if the run were
/// long enough to escape the window entirely on both sides, would report a
/// violation whose recorded span was wrong.
#[test]
fn tst_cons_009_a_run_longer_than_the_radius_is_not_truncated() {
    let runs_only = set_with(vec![(
        "runs",
        Rule::Homopolymer {
            base: HomopolymerBase::Any,
            max_run: 4,
        },
    )]);
    let checker = ConstraintChecker::compile(&runs_only).unwrap();
    assert_eq!(
        checker.influence_radius_nt(),
        Some(0),
        "a homopolymer rule contributes no fixed radius; the growth step covers it"
    );

    // 80 explicit As plus the leading A of the trailing `ACGT` — an 81 nt run
    // spanning 4..85, far wider than the radius this rule set contributes.
    let mut dna = b"ACGT".to_vec();
    dna.extend(std::iter::repeat(b'A').take(80));
    dna.extend_from_slice(b"ACGT");
    let view = ConstructView::new(&dna, vec![(0, dna.len())]);

    let middle = 4 + 40;
    let regional = checker.check_region(&view, middle, middle + 3).unwrap();
    let expected = intersecting(&checker.check(&view), middle, middle + 3);
    assert_eq!(regional, expected);
    assert_eq!(regional.len(), 1);
    assert_eq!((regional[0].start_nt, regional[0].end_nt), (4, 85));
}
