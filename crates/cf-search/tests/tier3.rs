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
    nondominated, run_generations, seed_population, stagnated, Individual, ObjectiveEvaluator,
    Selection, Termination, Tier3Config, Tier3Request,
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
        l0: None,
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
        inherited: None,
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

// --------------------------------------------------------------------------- //
// TST-SRCH-018 — the M >= 4 switch
// --------------------------------------------------------------------------- //

/// `M` conflicting objectives, so the front has real structure at any width.
#[derive(Debug)]
struct ManyWay(usize);

impl ObjectiveEvaluator for ManyWay {
    fn evaluate(&self, genotype: &Genotype) -> Result<Vec<Option<f64>>> {
        let codons = genotype.codons();
        // Each objective reads a different residue class, so they conflict
        // rather than being transforms of one another - a fixture whose
        // objectives all agreed would make every member non-dominated for the
        // wrong reason and say nothing about niching.
        Ok((0..self.0)
            .map(|axis| {
                let total: usize = codons
                    .iter()
                    .enumerate()
                    .filter(|(index, _)| index % self.0 == axis)
                    .map(|(_, c)| c.index_usize())
                    .sum();
                Some(total as f64 / codons.len() as f64)
            })
            .collect())
    }
    fn objective_count(&self) -> usize {
        self.0
    }
}

fn run_with(objectives: usize, generations: u32) -> cf_search::tier3::Tier3Report {
    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let evaluator = ManyWay(objectives);
    let config = config(12);
    let request = request(&checker, &mix, &weights, &evaluator, &config);
    let start = seed_population(&[base()], &request).unwrap();
    run_generations(&start, &request, generations).unwrap()
}

#[test]
fn tst_srch_018_the_selection_rule_switches_at_four_objectives() {
    // The rule lives in one place, so it cannot be spelled differently in two.
    assert_eq!(Selection::for_objectives(1), Selection::CrowdingDistance);
    assert_eq!(Selection::for_objectives(3), Selection::CrowdingDistance);
    assert_eq!(Selection::for_objectives(4), Selection::ReferenceDirections);
    assert_eq!(Selection::for_objectives(7), Selection::ReferenceDirections);

    // And the run reports which one it used, rather than leaving a caller to
    // infer it from behaviour.
    assert_eq!(run_with(2, 3).selection, Selection::CrowdingDistance);
    assert_eq!(run_with(5, 3).selection, Selection::ReferenceDirections);
}

#[test]
fn tst_srch_018_a_five_objective_run_keeps_its_population_and_front() {
    let report = run_with(5, 4);

    assert_eq!(report.checkpoint.members.len(), 12, "niching lost members");
    assert!(!report.front.is_empty());
    // Mutually non-dominated by construction, checked rather than assumed.
    for (i, a) in report.front.iter().enumerate() {
        for (j, b) in report.front.iter().enumerate() {
            if i != j {
                assert_eq!(
                    constraint_dominance(&a.ranked(), &b.ranked()),
                    Dominance::Neither
                );
            }
        }
    }
}

#[test]
fn tst_srch_018_niching_selection_is_reproducible() {
    // The whole point of DEP-SRCH-009: with no randomness in the tie-breaks,
    // two runs of the same seed plan agree, and so do two splits of one run.
    let whole = run_with(5, 4);
    let a: Vec<Vec<Option<f64>>> = whole
        .checkpoint
        .members
        .iter()
        .map(|m| m.objectives.clone())
        .collect();
    let again = run_with(5, 4);
    let b: Vec<Vec<Option<f64>>> = again
        .checkpoint
        .members
        .iter()
        .map(|m| m.objectives.clone())
        .collect();
    assert_eq!(a, b);
}

#[test]
fn tst_srch_018_a_single_objective_run_still_agrees_with_a_scalar_best() {
    // front()/best() agreement, at the tier-3 level. With one objective the
    // non-dominated set is the minimum, so the front must be exactly the
    // members attaining it - a second ranking notion introduced by NSGA-III
    // would show up here as a front that disagrees.
    #[derive(Debug)]
    struct OneWay;
    impl ObjectiveEvaluator for OneWay {
        fn evaluate(&self, genotype: &Genotype) -> Result<Vec<Option<f64>>> {
            let total: usize = genotype.codons().iter().map(|c| c.index_usize()).sum();
            Ok(vec![Some(total as f64)])
        }
        fn objective_count(&self) -> usize {
            1
        }
    }

    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let evaluator = OneWay;
    let config = config(10);
    let request = request(&checker, &mix, &weights, &evaluator, &config);
    let start = seed_population(&[base()], &request).unwrap();
    let report = run_generations(&start, &request, 4).unwrap();

    let best = report
        .checkpoint
        .members
        .iter()
        .filter(|m| m.is_feasible())
        .map(|m| m.objectives[0].unwrap())
        .fold(f64::INFINITY, f64::min);
    for member in &report.front {
        assert_eq!(
            member.objectives[0].unwrap(),
            best,
            "the front contains a member that is not the scalar minimum"
        );
    }
    assert!(!report.front.is_empty());
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
            inherited: None,
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

// --------------------------------------------------------------------------- //
// TST-SRCH-019 — termination
// --------------------------------------------------------------------------- //

#[test]
fn tst_srch_019_stagnation_is_relative_and_needs_a_full_window() {
    // A window longer than the history is never stagnant: a run has to have
    // had the chance to improve before failing to is evidence of anything.
    assert!(!stagnated(&[1.0, 1.0, 1.0], 5, 1e-3));
    assert!(!stagnated(&[], 3, 1e-3));
    assert!(
        !stagnated(&[1.0; 4], 0, 1e-3),
        "a zero window must not fire"
    );

    // Flat over the window: stagnant.
    assert!(stagnated(&[1.0, 1.0, 1.0, 1.0, 1.0], 3, 1e-3));

    // Still climbing: not stagnant.
    assert!(!stagnated(&[1.0, 2.0, 3.0, 4.0, 5.0], 3, 1e-3));

    // Relative, not absolute. The same absolute gain is stagnation at a large
    // scale and progress at a small one, which is why hypervolume - whose
    // magnitude depends on the reference point and the objective units -
    // cannot use an absolute epsilon.
    assert!(stagnated(&[1.0e6, 1.0e6, 1.0e6, 1.0e6 + 1.0], 3, 1e-3));
    assert!(!stagnated(&[1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3 + 1.0], 3, 1e-3));
}

/// A constant objective, so the population cannot improve at all.
///
/// The stagnation condition is *constructed*, not waited for. A fixture that
/// ran a real search and hoped it stalled would depend on where a particular
/// random walk went - which is how a motif fixture in this repository broke
/// twice when an unrelated RNG change shifted the trajectory.
#[derive(Debug)]
struct Constant;

impl ObjectiveEvaluator for Constant {
    fn evaluate(&self, _genotype: &Genotype) -> Result<Vec<Option<f64>>> {
        Ok(vec![Some(1.0), Some(1.0)])
    }
    fn objective_count(&self) -> usize {
        2
    }
}

fn run_until_termination(
    evaluator: &dyn ObjectiveEvaluator,
    window: usize,
    max_generations: u32,
) -> cf_search::tier3::Tier3Report {
    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let mut config = config(8);
    config.stagnation_window = window;
    config.max_generations = max_generations;
    let request = request(&checker, &mix, &weights, evaluator, &config);
    let start = seed_population(&[base()], &request).unwrap();
    run_generations(&start, &request, max_generations).unwrap()
}

#[test]
fn tst_srch_019_a_run_that_cannot_improve_terminates_on_stagnation() {
    let report = run_until_termination(&Constant, 3, 50);

    assert_eq!(
        report.termination,
        Termination::HypervolumeStagnation,
        "a run with a constant objective ran to its generation cap instead of stopping"
    );
    assert!(report.termination.is_final());
    assert!(
        report.generations_run < 50,
        "termination fired only at the cap, which is not termination"
    );
    assert_eq!(report.termination.as_summary(), "hypervolume_stagnation");
}

#[test]
fn tst_srch_019_a_run_that_is_still_improving_reaches_its_generation_cap() {
    // The other half: a fixture where stagnation must *not* fire, so the test
    // above cannot pass by a condition that always holds. The window is longer
    // than the run, so stagnation is unreachable by construction.
    let report = run_until_termination(&Conflicting, 100, 5);

    assert_eq!(report.termination, Termination::GenerationLimit);
    assert_eq!(report.generations_run, 5);
}

#[test]
fn tst_srch_019_the_history_has_one_entry_per_generation_run() {
    let report = run_until_termination(&Conflicting, 100, 6);
    assert_eq!(report.checkpoint.hypervolume_history.len(), 6);
    assert_eq!(report.generations_run, 6);
    assert!(report
        .checkpoint
        .hypervolume_history
        .iter()
        .all(|v| v.is_finite()));
}

#[test]
fn tst_srch_019_the_history_travels_in_the_checkpoint() {
    // Termination depends on it, so a resumed run that started with an empty
    // history could not stagnate until it had run a whole window again - a job
    // checkpointed near convergence would keep going where an uninterrupted
    // one stopped.
    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let mut config = config(8);
    config.stagnation_window = 3;
    config.max_generations = 50;
    let evaluator = Constant;
    let request = request(&checker, &mix, &weights, &evaluator, &config);

    let start = seed_population(&[base()], &request).unwrap();
    let first = run_generations(&start, &request, 2).unwrap();
    assert_eq!(first.termination, Termination::SegmentComplete);
    assert_eq!(first.checkpoint.hypervolume_history.len(), 2);

    let resumed = run_generations(&first.checkpoint, &request, 10).unwrap();
    assert_eq!(
        resumed.termination,
        Termination::HypervolumeStagnation,
        "the resumed run did not inherit the history and could not stagnate"
    );
    assert!(resumed.checkpoint.hypervolume_history.len() > 2);
}

#[test]
fn tst_srch_019_a_segment_that_has_not_terminated_says_so() {
    // SegmentComplete is not a termination: the caller may continue. Reporting
    // it as a generation limit would tell a job it had converged when it had
    // only reached the end of a checkpoint interval.
    let report = run_until_termination(&Conflicting, 100, 20);
    assert!(report.termination.is_final());

    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let mut config = config(8);
    config.stagnation_window = 100;
    config.max_generations = 20;
    let evaluator = Conflicting;
    let request = request(&checker, &mix, &weights, &evaluator, &config);
    let start = seed_population(&[base()], &request).unwrap();
    let partial = run_generations(&start, &request, 3).unwrap();

    assert_eq!(partial.termination, Termination::SegmentComplete);
    assert!(!partial.termination.is_final());
}

/// Four objectives, the last always `UNKNOWN`.
///
/// So every member of every splitting front is unplaceable in objective space
/// and niching has to fall back on them rather than merely rank them last.
/// A fixture where *some* members are placeable cannot test this: the
/// placeable ones fill every slot and the fallback branch never runs, which is
/// what a first version of this test found - four generations in, not one
/// member carried an unknown and the assertion that the fixture bound caught
/// it.
#[derive(Debug)]
struct AlwaysUnknown;

impl ObjectiveEvaluator for AlwaysUnknown {
    fn evaluate(&self, genotype: &Genotype) -> Result<Vec<Option<f64>>> {
        let codons = genotype.codons();
        let total: usize = codons.iter().map(|c| c.index_usize()).sum();
        Ok(vec![
            Some(total as f64 / codons.len() as f64),
            Some((total % 7) as f64),
            Some((total % 5) as f64),
            None,
        ])
    }
    fn objective_count(&self) -> usize {
        4
    }
}

#[test]
fn tst_srch_018_nsga3_keeps_members_it_cannot_place_in_objective_space() {
    // §4: an UNKNOWN term is never imputed, and an individual carrying one is
    // not shown to be bad - only unmeasured. Dropping such members would leave
    // the population unable to refill itself the moment a provider missed a
    // budget for everyone, which is exactly when a search must not collapse.
    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let evaluator = AlwaysUnknown;
    let mut config = config(12);
    config.stagnation_window = 1000;
    let request = request(&checker, &mix, &weights, &evaluator, &config);

    let start = seed_population(&[base()], &request).unwrap();
    assert!(
        start
            .members
            .iter()
            .all(|m| m.objectives.iter().any(Option::is_none)),
        "the fixture must actually produce unplaceable members"
    );

    let report = run_generations(&start, &request, 3).unwrap();

    assert_eq!(report.selection, Selection::ReferenceDirections);
    assert_eq!(
        report.checkpoint.members.len(),
        12,
        "selection dropped members it could not place in objective space, so a \
         population whose provider missed its budget would collapse"
    );
    assert!(!report.front.is_empty());
}

#[test]
fn tst_srch_018_nsga3_prefers_a_placeable_member_to_an_unplaceable_one() {
    // The other half of the rule: unplaceable members are ranked *last*, not
    // treated as equals. With a mixed splitting front the placeable members
    // fill the slots first.
    let ranked = vec![
        Ranked::feasible(vec![Some(1.0), Some(1.0), Some(1.0), None]),
        Ranked::feasible(vec![Some(2.0), Some(2.0), Some(2.0), Some(2.0)]),
    ];
    let (known, unknown) = cf_metrics::placeable(&ranked, &[0, 1]);
    assert_eq!(known, vec![1]);
    assert_eq!(unknown, vec![0]);
}

#[test]
fn tst_srch_019_a_resumed_segment_shorter_than_the_window_still_stagnates() {
    // The history has to be inherited, not restarted. With a window of 3 and a
    // resumed segment of 2 generations, stagnation is unreachable from an empty
    // history - so a run that terminates here can only have done so on the
    // history it carried. The earlier test used a long resumed segment, which
    // could stagnate on its own and so could not tell the two apart.
    let checker = checker();
    let mix = OperatorMix::without_context_model();
    let weights = CodonWeights::uniform();
    let mut config = config(8);
    config.stagnation_window = 3;
    config.max_generations = 100;
    let evaluator = Constant;
    let request = request(&checker, &mix, &weights, &evaluator, &config);

    let start = seed_population(&[base()], &request).unwrap();
    let first = run_generations(&start, &request, 3).unwrap();
    assert_eq!(first.termination, Termination::SegmentComplete);
    assert_eq!(first.checkpoint.hypervolume_history.len(), 3);

    let resumed = run_generations(&first.checkpoint, &request, 2).unwrap();
    assert_eq!(
        resumed.termination,
        Termination::HypervolumeStagnation,
        "a two-generation segment terminated on a three-generation window only if \
         it inherited the history"
    );
    assert_eq!(resumed.generations_run, 1);
}
