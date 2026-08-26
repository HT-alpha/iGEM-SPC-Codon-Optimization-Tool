//! `TST-CONS-006` — each of the five feasibility checks fires on a
//! purpose-built infeasible input and stays silent on a feasible one.
//!
//! Ticket P1-19. One test per check, plus `feasible_set_produces_no_finding`,
//! which is the test that would catch an analyser that reported everything as
//! infeasible — the failure mode a feasibility gate is most likely to have and
//! least likely to notice.

use cf_constraints::feasibility::{analyse, CheckId, DesignSpace};
use cf_constraints::spec::{
    Constraint, ConstraintClass, ConstraintSet, HomopolymerBase, RepairPolicy, Rule, Scope, Strand,
};
use cf_constraints::{ConstraintChecker, KmerMasking};
use cf_types::CodonSequence;

fn set_with(rules: Vec<(&str, Rule)>) -> ConstraintSet {
    ConstraintSet {
        schema_version: Some("1.0.0".into()),
        constraint_set_id: "cst_feas_test".into(),
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

fn analyse_one(
    construct: &str,
    regions: &[(usize, &CodonSequence)],
    set: &ConstraintSet,
) -> cf_constraints::FeasibilityReport {
    let checker = ConstraintChecker::compile(set).unwrap();
    let space = DesignSpace::new(construct.as_bytes(), regions.to_vec());
    analyse(&space, &checker, KmerMasking::Assembled).unwrap()
}

/// A short, fully-editable coding region: M-K-R-I-L-stop.
fn demo_cds() -> CodonSequence {
    CodonSequence::from_cds("ATGAAACGTATTCTGTAA", 11).unwrap()
}

// --------------------------------------------------------------------------- //

/// Check 1: a blocking violation confined to immutable sequence.
///
/// This is where `DEP-CONS-002` routes non-actionable violations. The
/// reasoning is not re-derived here: `Violation::actionable` is the mechanism
/// and `DEP-CONS-002` is the source of truth.
#[test]
fn tst_cons_006_check_1_locked_span() {
    let cds = demo_cds();
    // GGCC twice, both outside the coding region at 6..24, so no synonymous
    // choice can move either one. GGCC is used rather than ACGT because the
    // demo CDS itself contains ACGT at offset 5 — which would make the repeat
    // actionable and route it away from this check.
    //   index  0 1 2 3 4 5 | 6 .. 23 (CDS) | 24 25 26 27 28 29
    //          G G C C T T |               |  G  G  C  C  A  A
    let construct = format!("GGCCTT{}GGCCAA", cds.cds_string());
    assert!(!cds.cds_string().contains("GGCC"));
    let set = set_with(vec![(
        "c_kmer",
        Rule::KmerRepeat {
            k: 4,
            max_occurrences: 1,
            count_overlapping: true,
            canonical_strand: false,
        },
    )]);

    let report = analyse_one(&construct, &[(6, &cds)], &set);
    let findings: Vec<_> = report.findings_for(CheckId::LockedSpan).collect();
    assert_eq!(findings.len(), 1, "{:?}", report.findings);
    assert_eq!(findings[0].code, "CF-FEAS-001");
    assert_eq!(findings[0].constraint_id.as_deref(), Some("c_kmer"));
    assert!(findings[0].detail.contains("GGCC"));
    assert!(findings[0].suggested_relaxation.is_some());
    assert!(
        !findings[0].exact,
        "the locked-span check is a necessary condition, not a proof"
    );
}

/// Check 2: pigeonhole over k-mers.
#[test]
fn tst_cons_006_check_2_counting_bound() {
    let cds = demo_cds();
    let construct = cds.cds_string();
    // 18 nt gives 15 overlapping 4-mers; 4^4 = 256 keys at 1 each is ample, so
    // this must NOT fire. k=4 with max_occurrences 0 makes capacity 0.
    let ok = set_with(vec![(
        "c_kmer",
        Rule::KmerRepeat {
            k: 4,
            max_occurrences: 1,
            count_overlapping: true,
            canonical_strand: false,
        },
    )]);
    assert_eq!(
        analyse_one(&construct, &[(0, &cds)], &ok)
            .findings_for(CheckId::CountingBound)
            .count(),
        0
    );

    // Capacity zero: no sequence at all can satisfy it.
    let impossible = set_with(vec![(
        "c_kmer",
        Rule::KmerRepeat {
            k: 4,
            max_occurrences: 0,
            count_overlapping: true,
            canonical_strand: false,
        },
    )]);
    let report = analyse_one(&construct, &[(0, &cds)], &impossible);
    let findings: Vec<_> = report.findings_for(CheckId::CountingBound).collect();
    assert_eq!(findings.len(), 1, "{:?}", report.findings);
    assert_eq!(findings[0].code, "CF-FEAS-002");
    assert!(findings[0].exact, "the counting bound is a proof");
    let relaxation = findings[0].suggested_relaxation.as_ref().unwrap();
    assert!(
        relaxation.contains("max_occurrences"),
        "must name what to relax and by how much: {relaxation}"
    );
}

/// Check 3: a position whose every synonym violates something, in a way no
/// other position could rescue.
#[test]
fn tst_cons_006_check_3_per_position_exhaustion() {
    // Leucine has six codons under table 11: TTA TTG CTT CTC CTA CTG. Ban all
    // six, and the leucine position has no admissible choice. The ban is a
    // codon_ban rule, whose violation span is exactly the codon — so it is
    // confined to that one editable position and the check is sound.
    let cds = demo_cds();
    let construct = cds.cds_string();
    let set = set_with(vec![(
        "c_ban_leu",
        Rule::CodonBan {
            codons: ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]
                .iter()
                .map(|s| s.to_string())
                .collect(),
            reason: Some("test".into()),
        },
    )]);

    let report = analyse_one(&construct, &[(0, &cds)], &set);
    let findings: Vec<_> = report.findings_for(CheckId::PerPosition).collect();
    assert_eq!(findings.len(), 1, "{:?}", report.findings);
    assert_eq!(findings[0].code, "CF-FEAS-003");
    // Codon 4 is the leucine.
    assert_eq!(findings[0].span_nt, Some((12, 15)));
    assert!(!findings[0].exact);

    // Banning only some of the six leaves a choice, so nothing fires.
    let partial = set_with(vec![(
        "c_ban_leu",
        Rule::CodonBan {
            codons: ["TTA", "TTG", "CTT"]
                .iter()
                .map(|s| s.to_string())
                .collect(),
            reason: None,
        },
    )]);
    assert_eq!(
        analyse_one(&construct, &[(0, &cds)], &partial)
            .findings_for(CheckId::PerPosition)
            .count(),
        0
    );
}

/// Check 4: exact product-automaton reachability (`DEP-CONS-003`).
#[test]
fn tst_cons_006_check_4_automaton_reachability_is_exact() {
    // Arginine has six codons: CGT CGC CGA CGG AGA AGG. Forbid every one as a
    // literal motif, and no assignment can pass that position. Because the
    // check is exact reachability rather than a beam, the finding is a proof.
    let cds = demo_cds();
    let construct = cds.cds_string();
    let rules: Vec<(&str, Rule)> = vec![
        ("c_m0", motif("CGT")),
        ("c_m1", motif("CGC")),
        ("c_m2", motif("CGA")),
        ("c_m3", motif("CGG")),
        ("c_m4", motif("AGA")),
        ("c_m5", motif("AGG")),
    ];
    let report = analyse_one(&construct, &[(0, &cds)], &set_with(rules));
    let findings: Vec<_> = report
        .findings_for(CheckId::AutomatonReachability)
        .collect();
    assert_eq!(findings.len(), 1, "{:?}", report.findings);
    assert_eq!(findings[0].code, "CF-FEAS-004");
    assert!(
        findings[0].exact,
        "exact reachability proves unsatisfiability; a beam could not"
    );
    // Codon 2 is the arginine.
    assert_eq!(findings[0].span_nt, Some((6, 9)));

    // Leave one arginine codon available and the region becomes reachable.
    let rules: Vec<(&str, Rule)> = vec![
        ("c_m0", motif("CGT")),
        ("c_m1", motif("CGC")),
        ("c_m2", motif("CGA")),
        ("c_m3", motif("CGG")),
        ("c_m4", motif("AGA")),
    ];
    assert_eq!(
        analyse_one(&construct, &[(0, &cds)], &set_with(rules))
            .findings_for(CheckId::AutomatonReachability)
            .count(),
        0,
        "one surviving synonym is enough"
    );
}

/// The exact check must not report a false positive where a beam would.
///
/// A beam of width 1 driven by a bad heuristic can lose the only surviving path
/// and report unsatisfiability. Exact reachability keeps every reachable state,
/// so a satisfiable set is never reported infeasible.
#[test]
fn tst_cons_006_check_4_no_false_positive_on_a_narrow_survivor() {
    // A long region where exactly one assignment at each arginine survives.
    let cds =
        CodonSequence::from_cds("ATGCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTCGTTAA", 11).unwrap();
    let construct = cds.cds_string();
    // Forbid five of the six arginine codons: CGT is the sole survivor, and it
    // must be chosen at all fifteen positions.
    let rules: Vec<(&str, Rule)> = vec![
        ("c_m1", motif("CGC")),
        ("c_m2", motif("CGA")),
        ("c_m3", motif("CGG")),
        ("c_m4", motif("AGA")),
        ("c_m5", motif("AGG")),
    ];
    let report = analyse_one(&construct, &[(0, &cds)], &set_with(rules));
    assert_eq!(
        report.findings_for(CheckId::AutomatonReachability).count(),
        0,
        "a single surviving path must be found: {:?}",
        report.findings
    );
}

/// Homopolymer limits are regular, so they participate in check 4.
#[test]
fn tst_cons_006_check_4_includes_homopolymers() {
    // M-W-G. Tryptophan has exactly one codon, TGG, and every glycine codon
    // begins GG, so a W immediately followed by a G always produces a run of
    // four Gs. No synonymous assignment avoids it.
    //
    //   A T G | T G G | G G x | T A A
    //               \--GGGG--/
    let cds = CodonSequence::from_cds("ATGTGGGGTTAA", 11).unwrap();
    assert_eq!(cds.protein_string(), "MWG*");
    let construct = cds.cds_string();

    let impossible = set_with(vec![(
        "c_hp",
        Rule::Homopolymer {
            base: HomopolymerBase::Any,
            max_run: 3,
        },
    )]);
    let report = analyse_one(&construct, &[(0, &cds)], &impossible);
    assert_eq!(
        report.findings_for(CheckId::AutomatonReachability).count(),
        1,
        "a forced four-G run must be proved unreachable: {:?}",
        report.findings
    );

    // Allow runs of four and the same region becomes reachable.
    let ok = set_with(vec![(
        "c_hp",
        Rule::Homopolymer {
            base: HomopolymerBase::Any,
            max_run: 4,
        },
    )]);
    assert_eq!(
        analyse_one(&construct, &[(0, &cds)], &ok)
            .findings_for(CheckId::AutomatonReachability)
            .count(),
        0
    );
}

/// Check 5: a GC band no synonymous assignment can reach.
#[test]
fn tst_cons_006_check_5_gc_attainability() {
    let cds = demo_cds();
    let construct = cds.cds_string();

    // A band the region can reach.
    let ok = set_with(vec![(
        "c_gc",
        Rule::GcWindow {
            window_nt: 18,
            min: 0.0,
            max: 1.0,
        },
    )]);
    assert_eq!(
        analyse_one(&construct, &[(0, &cds)], &ok)
            .findings_for(CheckId::GcAttainability)
            .count(),
        0
    );

    // A band it cannot: this 6-residue peptide cannot reach 95% GC.
    let impossible = set_with(vec![(
        "c_gc",
        Rule::GcWindow {
            window_nt: 18,
            min: 0.95,
            max: 1.0,
        },
    )]);
    let report = analyse_one(&construct, &[(0, &cds)], &impossible);
    let findings: Vec<_> = report.findings_for(CheckId::GcAttainability).collect();
    assert_eq!(findings.len(), 1, "{:?}", report.findings);
    assert_eq!(findings[0].code, "CF-FEAS-005");
    assert!(findings[0].exact, "attainability is exact for one window");
    let relaxation = findings[0].suggested_relaxation.as_ref().unwrap();
    assert!(
        relaxation.contains("widen the band"),
        "must report the attainable band: {relaxation}"
    );
}

/// The test that catches an analyser which reports everything as infeasible.
#[test]
fn tst_cons_006_feasible_set_produces_no_finding() {
    let cds = demo_cds();
    let construct = format!("GGATCCTTT{}TTTGGATCC", cds.cds_string());
    let set = set_with(vec![
        (
            "c_bsai",
            Rule::RestrictionSite {
                enzyme: "BsaI".into(),
                site: "GGTCTC".into(),
                strand: Strand::Both,
                max_occurrences: 0,
            },
        ),
        (
            "c_hp",
            Rule::Homopolymer {
                base: HomopolymerBase::Any,
                max_run: 6,
            },
        ),
        (
            "c_kmer",
            Rule::KmerRepeat {
                k: 8,
                max_occurrences: 2,
                count_overlapping: true,
                canonical_strand: false,
            },
        ),
        (
            "c_gc",
            Rule::GcWindow {
                window_nt: 20,
                min: 0.1,
                max: 0.9,
            },
        ),
    ]);

    let report = analyse_one(&construct, &[(9, &cds)], &set);
    assert!(
        report.feasible,
        "a satisfiable set must produce no finding: {:?}",
        report.findings
    );
    assert_eq!(report.checks_run.len(), 5, "all five checks must run");
    assert!(
        report.necessary_conditions_only,
        "the report must state that passing is not a guarantee"
    );
}

/// All five checks are wired, distinguishable, and correctly classified.
#[test]
fn tst_cons_006_all_five_checks_are_accounted_for() {
    assert_eq!(CheckId::ALL.len(), 5);
    let names: Vec<&str> = CheckId::ALL.iter().map(|c| c.as_str()).collect();
    assert_eq!(
        names,
        [
            "locked_span",
            "counting_bound",
            "per_position",
            "automaton_reachability",
            "gc_attainability"
        ]
    );
    let codes: Vec<&str> = CheckId::ALL.iter().map(|c| c.error_code()).collect();
    assert_eq!(
        codes,
        [
            "CF-FEAS-001",
            "CF-FEAS-002",
            "CF-FEAS-003",
            "CF-FEAS-004",
            "CF-FEAS-005"
        ]
    );
    // Three exact, two necessary-condition. Getting this wrong would let a
    // suspicion be acted on as a proof.
    let exact: Vec<&str> = CheckId::ALL
        .iter()
        .filter(|c| c.is_exact())
        .map(|c| c.as_str())
        .collect();
    assert_eq!(
        exact,
        [
            "counting_bound",
            "automaton_reachability",
            "gc_attainability"
        ]
    );
}

fn motif(literal: &str) -> Rule {
    Rule::ForbiddenMotif {
        motif: literal.into(),
        strand: Strand::Sense,
        max_occurrences: 0,
    }
}
