//! `TST-CONS-010` — a violation confined to a locked span is not actionable.
//!
//! Found by `P1-25`'s audit, and the fourth time the `DEP-CONS-002` boundary
//! has produced a defect from a new direction.
//!
//! `Violation::actionable` is documented as "whether any editable position is
//! involved, and therefore whether the search engine could do anything about
//! it". It was computed against the whole coding region, which is coarser: a
//! forbidden motif lying entirely inside a *locked* span satisfied it. The
//! consequences compounded across three layers, each of which looked correct on
//! its own —
//!
//! * the checker called every candidate infeasible, for something no candidate
//!   could change;
//! * the feasibility analyser's check 1, which exists to report exactly that
//!   case once at submission, uses the same predicate and stayed silent;
//! * tier 1's beam therefore emptied and reported `CF-SRCH-006`, a claim about
//!   the design space rather than about the input.
//!
//! So the case fell through all three. This test pins the corrected behaviour
//! at the layer that defines it.

use cf_constraints::feasibility::{analyse, CheckId, DesignSpace};
use cf_constraints::spec::{
    Constraint, ConstraintClass, ConstraintSet, RepairPolicy, Rule, Scope, Strand,
};
use cf_constraints::{ConstraintChecker, ConstructView, KmerMasking};
use cf_types::{CodonSequence, LockReason};

/// `M K W K L *`, with lysine-tryptophan-lysine locked.
///
/// The locked codons spell `AAATGGAAA`, which is the forbidden motif — so the
/// violation is real, is inside the coding region, and is untouchable.
const CDS: &str = "ATGAAATGGAAACTGTAA";
const FLANK: &str = "GGAGGACATACC";

fn locked_sequence() -> CodonSequence {
    let mut seq = CodonSequence::from_cds(CDS, 1)
        .unwrap()
        .with_five_prime_flank(FLANK)
        .unwrap();
    seq.add_lock(1, 4, LockReason::UserLocked).unwrap();
    seq
}

fn constraint_set() -> ConstraintSet {
    ConstraintSet {
        schema_version: Some("1.0.0".into()),
        constraint_set_id: "cst_locked".into(),
        name: None,
        derived_from: vec![],
        constraints: vec![Constraint {
            id: "no_aaatggaaa".into(),
            class: ConstraintClass::Blocking,
            scope: Scope::Construct,
            cistron_id: None,
            weight: 1.0,
            enabled: true,
            rule: Rule::ForbiddenMotif {
                motif: "AAATGGAAA".into(),
                strand: Strand::Sense,
                max_occurrences: 0,
            },
        }],
        repair_policy: RepairPolicy::default(),
        digest: None,
    }
}

fn construct(sequence: &CodonSequence) -> Vec<u8> {
    let mut out = sequence.five_prime_flank().to_vec();
    for codon in sequence.codons() {
        out.extend_from_slice(&codon.to_bytes());
    }
    out.extend_from_slice(sequence.three_prime_flank());
    out
}

#[test]
fn tst_cons_010_a_motif_inside_a_locked_span_is_reported_but_not_actionable() {
    let sequence = locked_sequence();
    let checker = ConstraintChecker::compile(&constraint_set()).unwrap();
    let bytes = construct(&sequence);
    let cds = (FLANK.len(), FLANK.len() + CDS.len());

    // Without lock information the violation looks actionable — the behaviour
    // every caller had before, and correct for a construct with no locks.
    let coarse = ConstructView::new(&bytes, vec![cds]);
    let found = checker.check(&coarse);
    assert_eq!(
        found.len(),
        1,
        "the fixture must contain exactly one violation"
    );
    assert!(found[0].actionable);
    assert!(!checker.is_feasible(&coarse));

    // With it, the same violation is reported and is not actionable.
    let editable = ConstructView::new(&bytes, vec![cds])
        .with_editable_spans(vec![(cds.0, cds.0 + 3), (cds.0 + 12, cds.1)]);
    let found = checker.check(&editable);
    assert_eq!(found.len(), 1, "the violation is still reported");
    assert!(
        !found[0].actionable,
        "a motif spelled entirely by locked codons cannot be acted on"
    );
    assert!(
        checker.is_feasible(&editable),
        "no candidate can fix it, so it must not make every candidate infeasible"
    );
    assert_eq!(checker.unfixable(&editable).len(), 1);
}

/// And the analyser now reports it once, which is where it belongs.
#[test]
fn tst_cons_010_the_analyser_reports_it_as_check_1() {
    let sequence = locked_sequence();
    let checker = ConstraintChecker::compile(&constraint_set()).unwrap();
    let bytes = construct(&sequence);
    let space = DesignSpace::new(&bytes, vec![(FLANK.len(), &sequence)]);

    let report = analyse(&space, &checker, KmerMasking::Assembled).unwrap();
    let locked: Vec<_> = report
        .findings
        .iter()
        .filter(|f| f.check == CheckId::LockedSpan)
        .collect();
    assert_eq!(
        locked.len(),
        1,
        "check 1 must fire; findings were {:?}",
        report.findings
    );
    assert_eq!(locked[0].code, "CF-FEAS-001");
    assert!(locked[0].detail.contains("no_aaatggaaa"));
    assert!(locked[0]
        .suggested_relaxation
        .as_ref()
        .is_some_and(|s| s.contains("immutable sequence")));
}

/// The fallback is exact for a construct with no locks, so nothing that
/// predates the field changes behaviour.
#[test]
fn tst_cons_010_an_unlocked_construct_is_unaffected() {
    let sequence = CodonSequence::from_cds(CDS, 1)
        .unwrap()
        .with_five_prime_flank(FLANK)
        .unwrap();
    let checker = ConstraintChecker::compile(&constraint_set()).unwrap();
    let bytes = construct(&sequence);
    let cds = (FLANK.len(), FLANK.len() + CDS.len());

    let without = ConstructView::new(&bytes, vec![cds]);
    // `from_cds` locks the start and stop codons only, so everything between
    // them is editable and the motif overlaps it.
    let editable = cf_constraints::repair::editable_spans(&[(FLANK.len(), sequence.clone())]);
    let with = ConstructView::new(&bytes, vec![cds]).with_editable_spans(editable);

    assert_eq!(checker.check(&without), checker.check(&with));
    assert!(checker.check(&with)[0].actionable);
}
