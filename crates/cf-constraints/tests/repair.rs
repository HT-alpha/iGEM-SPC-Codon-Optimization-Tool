//! `TST-CONS-002`, `TST-CONS-003`, `TST-CONS-004` — the repair operator.
//!
//! Ticket P1-17.

use std::collections::BTreeSet;

use cf_constraints::repair::{
    declared_bound, partition, repair, ActionableViolation, GiveUpReason, NoObjective,
    ObjectiveDelta, RepairRng, RepairTarget,
};
use cf_constraints::spec::{
    Constraint, ConstraintClass, ConstraintSet, HomopolymerBase, RepairPolicy, Rule, Scope, Strand,
};
use cf_constraints::{ConstraintChecker, KmerMasking, TabuLedger};
use cf_types::{Codon, CodonSequence};

fn set_with(rules: Vec<(&str, Rule)>, policy: RepairPolicy) -> ConstraintSet {
    ConstraintSet {
        schema_version: Some("1.0.0".into()),
        constraint_set_id: "cst_repair_test".into(),
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
        repair_policy: policy,
        digest: None,
    }
}

fn site(enzyme: &str, seq: &str) -> Rule {
    Rule::RestrictionSite {
        enzyme: enzyme.into(),
        site: seq.into(),
        strand: Strand::Both,
        max_occurrences: 0,
    }
}

/// `M K R I L *` with a HindIII site planted across codons 2–3.
///
/// AAGCTT spans `AAG|CTT`, so either codon can break it: K has AAA/AAG and L
/// has six codons. That gives a violation with several minimal fixes, which is
/// what the diversity half of `TST-CONS-002` needs.
fn planted_target() -> (RepairTarget, ConstraintSet) {
    //          0     1     2     3     4     5
    let cds = CodonSequence::from_cds("ATGAAGCTTATTCTGTAA", 11).unwrap();
    assert_eq!(cds.protein_string(), "MKLIL*");
    let construct = format!("TTTTTT{}TTTTTT", cds.cds_string());
    let target = RepairTarget::new(construct.into_bytes(), vec![(6, cds)]).unwrap();
    let set = set_with(
        vec![("c_hindiii", site("HindIII", "AAGCTT"))],
        RepairPolicy::default(),
    );
    (target, set)
}

#[test]
fn tst_cons_002_repairs_a_planted_site_minimally() {
    let (mut target, set) = planted_target();
    let checker = ConstraintChecker::compile(&set).unwrap();
    let before = target.proteins();

    assert!(
        !checker.is_feasible(&target.view(KmerMasking::Assembled)),
        "the fixture must start infeasible"
    );

    let mut tabu = TabuLedger::default();
    let mut rng = RepairRng::new(1);
    let outcome = repair(
        &mut target,
        &checker,
        KmerMasking::Assembled,
        &NoObjective,
        &mut tabu,
        &mut rng,
    )
    .unwrap();

    assert!(outcome.gave_up.is_empty(), "{:?}", outcome.gave_up);
    assert!(outcome.residual.is_empty(), "{:?}", outcome.residual);
    assert_eq!(
        outcome.edits.len(),
        1,
        "a one-codon fix exists, so repair must not use more: {:?}",
        outcome.edits
    );
    assert!(checker.is_feasible(&target.view(KmerMasking::Assembled)));
    assert_eq!(target.proteins(), before, "the protein must be untouched");
}

/// The diversity half: randomized tie-breaking must not collapse onto one fix.
#[test]
fn tst_cons_002_randomized_tie_breaking_yields_distinct_fixes() {
    let mut distinct: BTreeSet<(usize, String)> = BTreeSet::new();

    for seed in 0..1_000u64 {
        let (mut target, set) = planted_target();
        let checker = ConstraintChecker::compile(&set).unwrap();
        let mut tabu = TabuLedger::default();
        let mut rng = RepairRng::new(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(1));
        let outcome = repair(
            &mut target,
            &checker,
            KmerMasking::Assembled,
            &NoObjective,
            &mut tabu,
            &mut rng,
        )
        .unwrap();
        assert!(outcome.residual.is_empty());
        for edit in &outcome.edits {
            distinct.insert((edit.codon_index, edit.to.to_string_upper()));
        }
    }

    assert!(
        distinct.len() > 1,
        "deterministic tie-breaking would collapse the population onto one \
         repair; got only {distinct:?}"
    );
}

/// Objective awareness: with a delta table, the cheapest fix wins.
#[test]
fn tst_cons_002_objective_aware_tie_breaking() {
    struct PreferCodon(Codon);
    impl ObjectiveDelta for PreferCodon {
        fn delta(&self, _region: usize, _codon_index: usize, to: Codon) -> Option<f64> {
            Some(if to == self.0 { -10.0 } else { 1.0 })
        }
    }

    // Steer repair toward a specific leucine codon and check it is chosen every
    // time, despite the randomized third component of the ordering key.
    let preferred = Codon::parse("CTC").unwrap();
    for seed in 0..50u64 {
        let (mut target, set) = planted_target();
        let checker = ConstraintChecker::compile(&set).unwrap();
        let mut tabu = TabuLedger::default();
        let mut rng = RepairRng::new(seed + 1);
        let outcome = repair(
            &mut target,
            &checker,
            KmerMasking::Assembled,
            &PreferCodon(preferred),
            &mut tabu,
            &mut rng,
        )
        .unwrap();
        assert!(outcome.residual.is_empty());
        assert!(
            outcome.edits.iter().any(|e| e.to == preferred),
            "objective cost must outrank the random tie-break: {:?}",
            outcome.edits
        );
        assert!(outcome.objective_cost < 0.0);
    }
}

/// `TST-CONS-003`: termination and bounds, including conflicting violations.
#[test]
fn tst_cons_003_terminates_within_its_declared_bound() {
    // Two constraints that cannot both be satisfied at the leucine: every
    // leucine codon is banned, and a homopolymer rule blocks the rest.
    let cds = CodonSequence::from_cds("ATGAAGCTTATTCTGTAA", 11).unwrap();
    let construct = format!("TTTTTT{}TTTTTT", cds.cds_string());
    let set = set_with(
        vec![
            ("c_hindiii", site("HindIII", "AAGCTT")),
            (
                "c_ban_all_leu",
                Rule::CodonBan {
                    codons: ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]
                        .iter()
                        .map(|s| s.to_string())
                        .collect(),
                    reason: Some("adversarial".into()),
                },
            ),
            (
                "c_ban_all_lys",
                Rule::CodonBan {
                    codons: ["AAA", "AAG"].iter().map(|s| s.to_string()).collect(),
                    reason: Some("adversarial".into()),
                },
            ),
        ],
        RepairPolicy::default(),
    );
    let checker = ConstraintChecker::compile(&set).unwrap();
    let mut target = RepairTarget::new(construct.into_bytes(), vec![(6, cds)]).unwrap();
    let before = target.proteins();

    let mut tabu = TabuLedger::default();
    let mut rng = RepairRng::new(7);
    let outcome = repair(
        &mut target,
        &checker,
        KmerMasking::Assembled,
        &NoObjective,
        &mut tabu,
        &mut rng,
    )
    .unwrap();

    // It terminated, it did not pretend to succeed, and it stayed in budget.
    assert!(
        !outcome.gave_up.is_empty(),
        "mutually conflicting rules must be given up on"
    );
    assert!(
        outcome.work.within_bound(),
        "evaluated {} candidates against a declared bound of {}",
        outcome.work.candidates_evaluated,
        outcome.work.declared_bound
    );
    assert_eq!(target.proteins(), before);
}

/// The give-up reason distinguishes a proof from a timeout.
#[test]
fn tst_cons_003_give_up_reason_distinguishes_proof_from_budget() {
    // Depth 1 only, so the search is exhaustive and failure is a proof.
    let policy = RepairPolicy {
        max_depth: 1,
        max_escalations: 0,
        ..RepairPolicy::default()
    };
    let cds = CodonSequence::from_cds("ATGAAGCTTATTCTGTAA", 11).unwrap();
    let construct = format!("TTTTTT{}TTTTTT", cds.cds_string());
    let set = set_with(
        vec![
            ("c_hindiii", site("HindIII", "AAGCTT")),
            (
                "c_ban_all_leu",
                Rule::CodonBan {
                    codons: ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]
                        .iter()
                        .map(|s| s.to_string())
                        .collect(),
                    reason: None,
                },
            ),
            (
                "c_ban_all_lys",
                Rule::CodonBan {
                    codons: ["AAA", "AAG"].iter().map(|s| s.to_string()).collect(),
                    reason: None,
                },
            ),
        ],
        policy,
    );
    let checker = ConstraintChecker::compile(&set).unwrap();
    let mut target = RepairTarget::new(construct.into_bytes(), vec![(6, cds)]).unwrap();
    let outcome = repair(
        &mut target,
        &checker,
        KmerMasking::Assembled,
        &NoObjective,
        &mut TabuLedger::default(),
        &mut RepairRng::new(3),
    )
    .unwrap();

    assert!(!outcome.gave_up.is_empty());
    assert!(
        outcome
            .gave_up
            .iter()
            .all(|g| g.reason == GiveUpReason::ProvablyNoFixInWindow),
        "an exhaustive depth-1 search that finds nothing is a proof, not a \
         timeout: {:?}",
        outcome.gave_up
    );
}

/// `TST-CONS-004`: protein and lock invariance.
///
/// The 10^6 cases are substitution attempts through `RepairTarget::substitute`,
/// which is the operator's only mutation path — so covering it covers every
/// edit repair can make. Full repair runs are exercised separately below, at a
/// count the full re-check cost allows.
#[test]
fn tst_cons_004_protein_and_lock_invariance() {
    let cds = CodonSequence::from_cds("ATGAAACGTATTCTGGATGAAGGCGCAACCGTTAGCCTGTAA", 11).unwrap();
    let construct = format!("GGGGGG{}CCCCCC", cds.cds_string());
    let expected_protein = cds.protein_string();
    let codon_count = cds.codon_count();
    let mut target = RepairTarget::new(construct.into_bytes(), vec![(6, cds)]).unwrap();

    let mut state = 0x2545_F491_4F6C_DD1Du64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };

    let mut accepted = 0u64;
    let mut rejected_locked = 0u64;
    let mut rejected_non_synonymous = 0u64;

    for _ in 0..1_000_000u64 {
        let index = (next() as usize) % codon_count;
        let codon = Codon::from_index((next() % 64) as u8).unwrap();
        match target.substitute(0, index, codon) {
            Ok(_) => accepted += 1,
            Err(e) if e.code == cf_types::code::SEQ_POSITION_LOCKED => rejected_locked += 1,
            Err(e) if e.code == cf_types::code::SEQ_NOT_SYNONYMOUS => rejected_non_synonymous += 1,
            Err(e) => panic!("unexpected error {e}"),
        }
        assert_eq!(
            target.regions()[0].1.protein_string(),
            expected_protein,
            "protein changed"
        );
    }

    // The counts are asserted so a `substitute` that silently accepted
    // everything could not pass by never erroring.
    assert!(accepted > 0, "some substitutions must succeed");
    assert!(
        rejected_locked > 0,
        "locked positions must be attempted and refused"
    );
    assert!(
        rejected_non_synonymous > 0,
        "non-synonymous codons must be attempted and refused"
    );

    // The start and stop codons are unchanged after a million attempts.
    let final_cds = target.regions()[0].1.cds_string();
    assert!(final_cds.starts_with("ATG"));
    assert!(final_cds.ends_with("TAA"));
}

#[test]
fn tst_cons_004_full_repair_runs_preserve_the_protein() {
    for seed in 0..500u64 {
        let (mut target, set) = planted_target();
        let checker = ConstraintChecker::compile(&set).unwrap();
        let before = target.proteins();
        repair(
            &mut target,
            &checker,
            KmerMasking::Assembled,
            &NoObjective,
            &mut TabuLedger::default(),
            &mut RepairRng::new(seed + 1),
        )
        .unwrap();
        assert_eq!(target.proteins(), before, "seed {seed}");
    }
}

/// `DEP-CONS-002`: the actionable boundary is a type, not a docstring.
#[test]
fn dep_cons_002_non_actionable_violations_are_not_constructible_for_repair() {
    // GGCC twice, both outside the coding region: a violation no candidate can
    // fix, which must route to the analyser rather than to repair.
    let cds = CodonSequence::from_cds("ATGAAACGTATTCTGTAA", 11).unwrap();
    let construct = format!("GGCCTT{}GGCCAA", cds.cds_string());
    let set = set_with(
        vec![(
            "c_kmer",
            Rule::KmerRepeat {
                k: 4,
                max_occurrences: 1,
                count_overlapping: true,
                canonical_strand: false,
            },
        )],
        RepairPolicy::default(),
    );
    let checker = ConstraintChecker::compile(&set).unwrap();
    let target = RepairTarget::new(construct.into_bytes(), vec![(6, cds)]).unwrap();

    let violations = checker.check(&target.view(KmerMasking::Assembled));
    assert_eq!(violations.len(), 1);
    assert!(!violations[0].actionable);

    // The type refuses to wrap it.
    assert!(ActionableViolation::new(&violations[0]).is_none());
    let (actionable, for_analyser) = partition(&violations);
    assert!(actionable.is_empty());
    assert_eq!(for_analyser.len(), 1);

    // And repair leaves it alone rather than thrashing on it.
    let mut target = target;
    let outcome = repair(
        &mut target,
        &checker,
        KmerMasking::Assembled,
        &NoObjective,
        &mut TabuLedger::default(),
        &mut RepairRng::new(11),
    )
    .unwrap();
    assert!(outcome.edits.is_empty(), "{:?}", outcome.edits);
    assert!(
        outcome.gave_up.is_empty(),
        "it was never repair's to attempt"
    );
    assert_eq!(outcome.work.candidates_evaluated, 0);
}

/// Repair honours the tabu ledger and feeds it.
#[test]
fn unit_repair_records_reverted_codons_in_the_tabu_ledger() {
    let (mut target, set) = planted_target();
    let checker = ConstraintChecker::compile(&set).unwrap();
    let mut tabu = TabuLedger::default();
    let outcome = repair(
        &mut target,
        &checker,
        KmerMasking::Assembled,
        &NoObjective,
        &mut tabu,
        &mut RepairRng::new(5),
    )
    .unwrap();

    assert_eq!(outcome.edits.len(), 1);
    let edit = outcome.edits[0];
    assert!(
        tabu.is_tabu(edit.region, edit.codon_index, edit.from),
        "the codon repair moved away from must be discouraged"
    );
    assert_eq!(
        tabu.remaining(edit.region, edit.codon_index, edit.from),
        set.repair_policy.tabu_tenure_generations
    );
}

#[test]
fn unit_declared_bound_is_computable_before_the_run() {
    let policy = RepairPolicy::default();
    let bound = declared_bound(&policy, 2, 6);
    assert!(bound > 0);
    // Monotone in every dimension a reviewer would expect.
    assert!(declared_bound(&policy, 4, 6) > bound);
    assert!(declared_bound(&policy, 2, 8) > bound);
    let deeper = RepairPolicy {
        max_depth: 5,
        ..policy.clone()
    };
    assert!(declared_bound(&deeper, 2, 6) > bound);
}

#[test]
fn unit_repair_target_rejects_an_inconsistent_region() {
    let cds = CodonSequence::from_cds("ATGAAACGTATTCTGTAA", 11).unwrap();
    // Offset points at the wrong bytes.
    assert!(
        RepairTarget::new(b"TTTTTTTTTTTTTTTTTTTTTTTT".to_vec(), vec![(0, cds.clone())]).is_err()
    );
    // Region runs past the end.
    assert!(RepairTarget::new(b"ATGAAA".to_vec(), vec![(0, cds)]).is_err());
}

#[test]
fn unit_homopolymer_repair() {
    // A run of seven As inside the coding region, fixable by a synonymous edit:
    // K is AAA/AAG, so AAG breaks the run.
    let cds = CodonSequence::from_cds("ATGAAAAAAATTCTGTAA", 11).unwrap();
    let construct = cds.cds_string();
    let set = set_with(
        vec![(
            "c_hp",
            Rule::Homopolymer {
                base: HomopolymerBase::Any,
                max_run: 5,
            },
        )],
        RepairPolicy::default(),
    );
    let checker = ConstraintChecker::compile(&set).unwrap();
    let mut target = RepairTarget::new(construct.into_bytes(), vec![(0, cds)]).unwrap();
    assert!(!checker.is_feasible(&target.view(KmerMasking::Assembled)));

    let outcome = repair(
        &mut target,
        &checker,
        KmerMasking::Assembled,
        &NoObjective,
        &mut TabuLedger::default(),
        &mut RepairRng::new(13),
    )
    .unwrap();
    assert!(outcome.residual.is_empty(), "{:?}", outcome.residual);
    assert!(checker.is_feasible(&target.view(KmerMasking::Assembled)));
}
