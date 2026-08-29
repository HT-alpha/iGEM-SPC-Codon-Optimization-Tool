//! TST-SRCH-007, TST-SRCH-016 and TST-SRCH-017: the NSGA-II population loop.
//!
//! Ticket P2-02.
//!
//! Two of these exist because the sequencing audit found `P2-02`'s only
//! declared acceptance criterion was `TST-SRCH-007`, a comparator truth table.
//! A loop that ranks correctly and loses its population on resume, or that
//! reports one member of a tied front, would have passed.

use cf_constraints::spec::{Constraint, ConstraintClass, ConstraintSet, Rule, Scope, Strand};
use cf_constraints::{ConstraintChecker, KmerMasking};
use cf_metrics::{constraint_dominance, Dominance, Ranked};
use cf_search::tier3::{
    nondominated, run_generations, seed_population, Individual, ObjectiveEvaluator, Tier3Config,
    Tier3Request,
};
use cf_search::{CodonWeights, Genotype, OperatorMix};
use cf_types::{CodonSequence, Result};

const CDS: &str = "ATGCTGCGTAGCGGCGTGGCGCCGACCAAACTGCGTAGCTAA";
const ROOT_SEED: u64 = 20260829;

fn checker() -> ConstraintChecker {
    ConstraintChecker::compile(&ConstraintSet {
        schema_version: Some("1.0.0".into()),
        constraint_set_id: "cst_tier3".into(),
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

/// Two conflicting objectives, both minimized.
///
/// GC fraction against a codon-index sum. They genuinely conflict, so the
/// population develops a real trade-off front rather than converging to one
/// point — a fixture with a single optimum would make every claim about front
/// structure vacuous, which is the failure mode `SYN-01`'s degenerate
/// configuration and the tier-2 determinism fixture both had to be rebuilt for.
#[derive(Debug)]
struct Conflicting;

impl ObjectiveEvaluator for Conflicting {
    fn evaluate(&self, genotype: &Genotype) -> Result<Vec<Option<f64>>> {
        let codons = genotype.codons();
        let gc: usize = codons
            .iter()
            .flat_map(|c| c.to_bytes())
            .filter(|b| *b == b'G' || *b == b'C')
            .count();
        let index: usize = codons.iter().map(|c| c.index_usize()).sum();
        Ok(vec![
            Some(gc as f64 / (codons.len() * 3) as f64),
            Some(index as f64 / codons.len() as f64),
        ])
    }
    fn objective_count(&self) -> usize {
        2
    }
}

/// Every individual scores identically, so the whole population is one
/// mutually non-dominated front.
#[derive(Debug)]
struct AllTied;

impl ObjectiveEvaluator for AllTied {
    fn evaluate(&self, _genotype: &Genotype) -> Result<Vec<Option<f64>>> {
        Ok(vec![Some(1.0), Some(1.0)])
    }
    fn objective_count(&self) -> usize {
        2
    }
}

fn config(population: usize) -> Tier3Config {
    let mut config = Tier3Config::from_planning(population);
    config.max_generations = 8;
    config
}

fn request<'a>(
    checker: &'a ConstraintChecker,
    mix: &'a OperatorMix,
    weights: &'a CodonWeights,
    evaluator: &'a dyn ObjectiveEvaluator,
    config: &'a Tier3Config,
) -> Tier3Request<'a> {
    Tier3Request {
        checker,
        mix,
        weights,
        evaluator,
        config,
        root_seed: ROOT_SEED,
    }
}

// --------------------------------------------------------------------------- //
// TST-SRCH-007 — constraint-domination truth table
// --------------------------------------------------------------------------- //

#[test]
fn tst_srch_007_constraint_domination_matches_debs_definition() {
    let feasible = |a: f64, b: f64| Ranked::feasible(vec![Some(a), Some(b)]);
    let infeasible = |a: f64, b: f64, cv: f64| Ranked::infeasible(vec![Some(a), Some(b)], cv);
    let unknown = |a: f64| Ranked::feasible(vec![Some(a), None]);

    // 1. Feasible beats infeasible, whatever the objectives say.
    assert_eq!(
        constraint_dominance(&feasible(9.0, 9.0), &infeasible(0.0, 0.0, 1.0)),
        Dominance::Left,
        "an infeasible individual with better objectives must not win"
    );
    assert_eq!(
        constraint_dominance(&infeasible(0.0, 0.0, 1.0), &feasible(9.0, 9.0)),
        Dominance::Right
    );

    // 2. Two infeasible compare on total violation, and never tie ambiguously.
    assert_eq!(
        constraint_dominance(&infeasible(9.0, 9.0, 1.0), &infeasible(0.0, 0.0, 2.0)),
        Dominance::Left
    );
    assert_eq!(
        constraint_dominance(&infeasible(0.0, 0.0, 2.0), &infeasible(9.0, 9.0, 2.0)),
        Dominance::Neither,
        "equal violation is incomparable, not a Pareto comparison"
    );

    // 3. Any UNKNOWN makes two feasible individuals incomparable. Never
    //    imputed: a mean or worst case would bias the front invisibly.
    assert_eq!(
        constraint_dominance(&unknown(0.0), &feasible(9.0, 9.0)),
        Dominance::Neither,
        "an unknown term was imputed rather than left incomparable"
    );
    assert_eq!(
        constraint_dominance(&feasible(0.0, 0.0), &unknown(9.0)),
        Dominance::Neither
    );
    assert_eq!(
        constraint_dominance(&unknown(0.0), &unknown(9.0)),
        Dominance::Neither
    );

    // 4. Otherwise Pareto.
    assert_eq!(
        constraint_dominance(&feasible(1.0, 1.0), &feasible(2.0, 2.0)),
        Dominance::Left
    );
    assert_eq!(
        constraint_dominance(&feasible(1.0, 5.0), &feasible(5.0, 1.0)),
        Dominance::Neither
    );
    assert_eq!(
        constraint_dominance(&feasible(1.0, 1.0), &feasible(1.0, 1.0)),
        Dominance::Neither
    );
    assert_eq!(
        constraint_dominance(&feasible(1.0, 1.0), &feasible(1.0, 2.0)),
        Dominance::Left,
        "weak dominance with one strict improvement must dominate"
    );
}

#[test]
fn tst_srch_007_an_unknown_never_hides_infeasibility() {
    // Feasibility is decided before unknowns are consulted, so an infeasible
    // individual cannot escape the comparison by having a missing term.
    let hiding = Ranked::infeasible(vec![Some(0.0), None], 3.0);
    let honest = Ranked::feasible(vec![Some(9.0), Some(9.0)]);
    assert_eq!(constraint_dominance(&honest, &hiding), Dominance::Left);
}

// --------------------------------------------------------------------------- //
// TST-SRCH-017 — Pareto ties are reported whole
// --------------------------------------------------------------------------- //

#[test]
fn tst_srch_017_a_tied_front_is_reported_whole_not_collapsed() {
    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let evaluator = AllTied;
    let config = config(8);
    let request = request(&checker, &mix, &weights, &evaluator, &config);

    let start = seed_population(&[base()], &request).unwrap();
    let report = run_generations(&start, &request, 3).unwrap();

    assert_eq!(
        report.front.len(),
        report.checkpoint.members.len(),
        "every member ties, so the whole population is one non-dominated front; \
         reporting fewer collapses the Pareto structure tier 3 exists to produce"
    );
    assert!(report.front.len() > 1);
}

#[test]
fn tst_srch_017_the_reported_front_is_independent_of_member_order() {
    let evaluator = Conflicting;
    let population: Vec<Individual> = [
        (vec![Some(1.0), Some(5.0)], 0.0),
        (vec![Some(5.0), Some(1.0)], 0.0),
        (vec![Some(3.0), Some(3.0)], 0.0),
        (vec![Some(9.0), Some(9.0)], 0.0),
        (vec![Some(0.0), Some(0.0)], 2.0),
    ]
    .into_iter()
    .map(|(objectives, violation)| Individual {
        genotype: base(),
        objectives,
        violation,
    })
    .collect();
    let _ = evaluator;

    let forward: Vec<Vec<Option<f64>>> = nondominated(&population)
        .iter()
        .map(|i| i.objectives.clone())
        .collect();

    let mut reversed = population.clone();
    reversed.reverse();
    let mut backward: Vec<Vec<Option<f64>>> = nondominated(&reversed)
        .iter()
        .map(|i| i.objectives.clone())
        .collect();
    backward.reverse();

    assert_eq!(
        forward, backward,
        "the reported front depends on member order"
    );
    assert_eq!(
        forward.len(),
        3,
        "expected the three mutually non-dominated members"
    );
    // The infeasible member is excluded even though its objectives dominate.
    assert!(!forward.contains(&vec![Some(0.0), Some(0.0)]));
}

#[test]
fn tst_srch_017_a_conflicting_objective_pair_produces_a_real_front() {
    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let evaluator = Conflicting;
    let config = config(16);
    let request = request(&checker, &mix, &weights, &evaluator, &config);

    let start = seed_population(&[base()], &request).unwrap();
    let report = run_generations(&start, &request, 6).unwrap();

    assert!(
        report.front.len() > 1,
        "the fixture converged to a single point, so nothing here tests front structure"
    );
    // Mutually non-dominated by construction, checked rather than assumed.
    for (i, a) in report.front.iter().enumerate() {
        for (j, b) in report.front.iter().enumerate() {
            if i == j {
                continue;
            }
            assert_eq!(
                constraint_dominance(&a.ranked(), &b.ranked()),
                Dominance::Neither,
                "a reported front member dominates another"
            );
        }
    }
}

// --------------------------------------------------------------------------- //
// TST-SRCH-016 — checkpoint completeness
// --------------------------------------------------------------------------- //

fn run_split(splits: &[u32]) -> (Vec<Vec<Option<f64>>>, u64, usize) {
    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let evaluator = Conflicting;
    let config = config(12);
    let request = request(&checker, &mix, &weights, &evaluator, &config);

    let mut checkpoint = seed_population(&[base()], &request).unwrap();
    let mut front = Vec::new();
    for &count in splits {
        let report = run_generations(&checkpoint, &request, count).unwrap();
        checkpoint = report.checkpoint;
        front = report.front;
    }
    (
        front.iter().map(|i| i.objectives.clone()).collect(),
        checkpoint.evaluations,
        checkpoint.tabu.len(),
    )
}

#[test]
fn tst_srch_016_a_resumed_population_matches_an_uninterrupted_one() {
    let whole = run_split(&[9]);
    for split in [
        vec![1, 8],
        vec![3, 3, 3],
        vec![8, 1],
        vec![1, 1, 1, 1, 1, 1, 1, 1, 1],
        vec![4, 2, 3],
    ] {
        assert_eq!(
            run_split(&split),
            whole,
            "a run split {split:?} diverged from an uninterrupted one"
        );
    }
}

#[test]
fn tst_srch_016_the_checkpoint_carries_the_objective_vectors() {
    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let evaluator = Conflicting;
    let config = config(10);
    let request = request(&checker, &mix, &weights, &evaluator, &config);

    let start = seed_population(&[base()], &request).unwrap();
    let report = run_generations(&start, &request, 4).unwrap();

    assert_eq!(report.checkpoint.members.len(), 10);
    for member in &report.checkpoint.members {
        assert_eq!(member.objectives.len(), 2);
        assert!(member.objectives.iter().all(Option::is_some));
    }
    assert!(report.checkpoint.evaluations >= 10);
    assert_eq!(report.checkpoint.generation, 4);
}

#[test]
fn tst_srch_016_the_generation_index_is_absolute_not_per_segment() {
    // Two segments must not both key their randomness on generation 0. This is
    // the defect P1-29 fixed for the mutation loop and P2-02 could reintroduce.
    let (front_a, _, _) = run_split(&[2, 2]);
    let (front_b, _, _) = run_split(&[4]);
    assert_eq!(front_a, front_b);
}

// --------------------------------------------------------------------------- //
// refusals
// --------------------------------------------------------------------------- //

#[test]
fn p2_02_more_than_three_objectives_is_refused_not_approximated() {
    #[derive(Debug)]
    struct FourWay;
    impl ObjectiveEvaluator for FourWay {
        fn evaluate(&self, _g: &Genotype) -> Result<Vec<Option<f64>>> {
            Ok(vec![Some(1.0); 4])
        }
        fn objective_count(&self) -> usize {
            4
        }
    }

    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let evaluator = FourWay;
    let config = config(4);
    let request = request(&checker, &mix, &weights, &evaluator, &config);

    let start = seed_population(&[base()], &request).unwrap();
    let error = run_generations(&start, &request, 1)
        .expect_err("M > 3 must be refused; crowding distance degrades badly as M grows");
    assert_eq!(error.code, "CF-SRCH-014");
    assert!(error.detail.contains("P2-03"), "{}", error.detail);
}

#[test]
fn p2_02_an_empty_seed_set_is_refused() {
    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let evaluator = Conflicting;
    let config = config(4);
    let request = request(&checker, &mix, &weights, &evaluator, &config);

    let error = seed_population(&[], &request).expect_err("an empty population has no front");
    assert_eq!(error.code, "CF-SRCH-014");
}

#[test]
fn p2_02_infeasible_individuals_are_kept_not_discarded() {
    // §4: infeasible individuals survive as stepping stones across a
    // constraint barrier. Discarding them is the standard way to make a
    // tightly-constrained codon problem return nothing.
    let population: Vec<Individual> = [(vec![Some(5.0), Some(5.0)], 1.0)]
        .into_iter()
        .map(|(objectives, violation)| Individual {
            genotype: base(),
            objectives,
            violation,
        })
        .collect();
    assert_eq!(nondominated(&population).len(), 1);
    assert!(!population[0].is_feasible());
}

#[test]
fn p2_02_the_masking_default_is_unchanged() {
    // EXIT-001/EXIT-002: legacy parity stays until an explicit Phase 1 exit
    // decision, and a new tier is exactly where a default quietly drifts.
    assert_eq!(
        Tier3Config::from_planning(8).masking,
        KmerMasking::default()
    );
}

/// A checker whose motif sits in the immutable 5' flank.
///
/// So the violation is real, detected, and *non-actionable*: no synonymous
/// assignment can remove it.
fn flank_checker() -> ConstraintChecker {
    ConstraintChecker::compile(&ConstraintSet {
        schema_version: Some("1.0.0".into()),
        constraint_set_id: "cst_flank".into(),
        name: None,
        derived_from: vec![],
        constraints: vec![Constraint {
            id: "no_ttgaca".into(),
            class: ConstraintClass::Blocking,
            scope: Scope::Construct,
            cistron_id: None,
            weight: 1.0,
            enabled: true,
            rule: Rule::ForbiddenMotif {
                motif: "TTGACAATT".into(),
                strand: Strand::Sense,
                max_occurrences: 0,
            },
        }],
        repair_policy: Default::default(),
        digest: None,
    })
    .unwrap()
}

#[test]
fn p2_02_a_violation_in_immutable_sequence_does_not_make_an_individual_infeasible() {
    // DEP-SRCH-008. The motif is in the promoter, which no codon choice can
    // change. Charging it to the individual would make the entire population
    // infeasible, and since constraint-domination makes any feasible
    // individual beat any infeasible one, the whole ranking would then be
    // driven by a quantity nothing in the population can influence.
    let checker = flank_checker();
    let genotype = base();

    // The fixture must actually contain the violation, or this asserts nothing.
    let construct: Vec<u8> = {
        let s = genotype.sequence();
        let mut out = s.five_prime_flank().to_vec();
        for codon in s.codons() {
            out.extend_from_slice(&codon.to_bytes());
        }
        out.extend_from_slice(s.three_prime_flank());
        out
    };
    let start = genotype.cds_offset_nt();
    let coarse =
        cf_constraints::ConstructView::new(&construct, vec![(start, start + genotype.len() * 3)]);
    assert_eq!(
        checker.check(&coarse).len(),
        1,
        "the fixture must plant a detectable violation in the flank"
    );

    assert_eq!(
        cf_search::tier3::actionable_violation(&genotype, &checker),
        0.0,
        "a violation no synonymous assignment can remove was counted against the individual"
    );
}

/// Forbids the only codon for tryptophan, so no assignment is feasible.
fn unsatisfiable_checker() -> ConstraintChecker {
    ConstraintChecker::compile(&ConstraintSet {
        schema_version: Some("1.0.0".into()),
        constraint_set_id: "cst_unsat".into(),
        name: None,
        derived_from: vec![],
        constraints: vec![Constraint {
            id: "no_tgg".into(),
            class: ConstraintClass::Blocking,
            scope: Scope::Cds,
            cistron_id: None,
            weight: 1.0,
            enabled: true,
            rule: Rule::ForbiddenMotif {
                motif: "TGG".into(),
                strand: Strand::Sense,
                max_occurrences: 0,
            },
        }],
        repair_policy: Default::default(),
        digest: None,
    })
    .unwrap()
}

fn tryptophan_base() -> Genotype {
    Genotype::new(
        CodonSequence::from_cds("ATGTGGCTGCGTAGCGGCGTGTAA", 1)
            .unwrap()
            .with_five_prime_flank("GGAGGACATACC")
            .unwrap(),
    )
}

#[test]
fn p2_02_an_unsatisfiable_population_is_ranked_rather_than_emptied() {
    // §4: infeasible individuals are kept, not discarded. They survive as
    // stepping stones across a constraint barrier, and discarding them is the
    // standard way to make a tightly-constrained codon design return nothing.
    let checker = unsatisfiable_checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let evaluator = Conflicting;
    let config = config(8);
    let request = request(&checker, &mix, &weights, &evaluator, &config);

    let start = seed_population(&[tryptophan_base()], &request).unwrap();
    // Tryptophan has one codon, so every individual violates and none can be
    // repaired - the fixture must confirm that rather than assume it.
    assert!(
        start.members.iter().all(|m| !m.is_feasible()),
        "the fixture is satisfiable, so it cannot test what happens when nothing is"
    );

    let report = run_generations(&start, &request, 3).unwrap();

    assert_eq!(
        report.checkpoint.members.len(),
        8,
        "the population was emptied by discarding infeasible individuals"
    );
    assert!(
        !report.front.is_empty(),
        "an infeasible run reported no front at all"
    );
    assert!(report.front.iter().all(|m| !m.is_feasible()));
}
