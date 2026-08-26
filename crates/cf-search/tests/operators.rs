//! `TST-SRCH-006` and the `P1-23` interop properties.
//!
//! The acceptance criterion is protein and lock invariance after 10⁶ randomized
//! operator applications. The rest of this file covers the properties `P1-24`,
//! `P1-18b` and `P1-27b` will rely on once they start, which the audit found
//! were not exercised by the declared criterion alone.

use std::sync::Arc;

use cf_constraints::spec::{
    Constraint, ConstraintClass, ConstraintSet, HomopolymerBase, RepairPolicy, Rule, Scope, Strand,
};
use cf_constraints::{ActionableViolation, ConstraintChecker, ConstructView, TabuLedger};
use cf_landscape_api::{Edit, Target};
use cf_search::context::ContextModel;
use cf_search::{
    BlockUniform, CodonWeights, CrossoverOperator, Genotype, Move, MoveKind, MutationContext,
    MutationOperator, OnePoint, OperatorMix, Population, PottsGuided, SearchRng, TwoPoint,
    UniformSynonymous, ViolationDirected,
};
use cf_types::{Codon, CodonSequence, Digest, GeneticCode, LockReason, Result};

// --------------------------------------------------------------------------- //
// fixtures
// --------------------------------------------------------------------------- //

/// `M L R S G V A P T K L R S *` — three six-fold families and several
/// four-fold ones, so an operator has real choice nearly everywhere.
const CDS: &str = "ATGCTGCGTAGCGGCGTGGCGCCGACCAAACTGCGTAGCTAA";

const FIVE_PRIME: &str = "TTGACAATTAATCATCGGCTCGTATAATGTGTGGA";
const THREE_PRIME: &str = "AAAGCCCGGATCCTCTAGAGGG";

fn sequence() -> CodonSequence {
    CodonSequence::from_cds(CDS, 1)
        .unwrap()
        .with_five_prime_flank(FIVE_PRIME)
        .unwrap()
        .with_three_prime_flank(THREE_PRIME)
        .unwrap()
}

/// The same sequence with codons 2..5 locked.
fn locked_sequence() -> CodonSequence {
    let mut seq = sequence();
    seq.add_lock(2, 5, LockReason::UserLocked).unwrap();
    seq
}

fn genotype() -> Genotype {
    Genotype::new(sequence())
}

/// A non-uniform weight table, so a weighted draw is distinguishable from a
/// uniform one.
fn weights() -> CodonWeights {
    let mut table = [0.0f64; 64];
    for (index, slot) in table.iter_mut().enumerate() {
        *slot = ((index as f64 * 17.0) % 23.0 + 1.0) / 23.0;
    }
    let weights = CodonWeights::from_relative_adaptiveness(table).unwrap();
    assert!(!weights.is_uniform());
    weights
}

/// A context model that is genuinely context-dependent.
///
/// Not a Potts model — it is a stand-in for one so the seam can be exercised
/// before `P3-05`. It conditions on the neighbouring codons, which is enough to
/// distinguish a real conditional from a context-free draw and therefore enough
/// to catch a `PottsGuided` that ignored its model.
#[derive(Debug)]
struct NeighbourModel;

impl ContextModel for NeighbourModel {
    fn name(&self) -> &str {
        "test_neighbour"
    }

    fn window(&self) -> usize {
        3
    }

    fn conditional(
        &self,
        sequence: &CodonSequence,
        index: usize,
        temperature: f64,
        out: &mut Vec<(Codon, f64)>,
    ) -> Result<()> {
        out.clear();
        let left = index
            .checked_sub(1)
            .map(|i| u32::from(sequence.codons()[i].index()))
            .unwrap_or(0);
        let right = sequence
            .codons()
            .get(index + 1)
            .map(|c| u32::from(c.index()))
            .unwrap_or(0);
        let context = f64::from(left ^ right);

        let mut total = 0.0;
        for codon in sequence.synonyms_at(index) {
            let affinity = ((f64::from(codon.index()) * 0.37 + context * 0.11).sin() + 1.5)
                .powf(1.0 / temperature);
            out.push((*codon, affinity));
            total += affinity;
        }
        for entry in out.iter_mut() {
            entry.1 /= total;
        }
        Ok(())
    }
}

fn constraint_set() -> ConstraintSet {
    ConstraintSet {
        schema_version: Some("1.0.0".into()),
        constraint_set_id: "cst_p1_23".into(),
        name: None,
        derived_from: vec![],
        constraints: vec![
            Constraint {
                id: "no_bamhi".into(),
                class: ConstraintClass::Blocking,
                scope: Scope::Construct,
                cistron_id: None,
                weight: 1.0,
                enabled: true,
                rule: Rule::RestrictionSite {
                    enzyme: "BamHI".into(),
                    site: "GGATCC".into(),
                    strand: Strand::Both,
                    max_occurrences: 0,
                },
            },
            Constraint {
                id: "no_ggcgc".into(),
                class: ConstraintClass::Blocking,
                scope: Scope::Construct,
                cistron_id: None,
                weight: 1.0,
                enabled: true,
                rule: Rule::ForbiddenMotif {
                    motif: "GGCGC".into(),
                    strand: Strand::Both,
                    max_occurrences: 0,
                },
            },
            Constraint {
                id: "runs".into(),
                class: ConstraintClass::Blocking,
                scope: Scope::Construct,
                cistron_id: None,
                weight: 1.0,
                enabled: true,
                rule: Rule::Homopolymer {
                    base: HomopolymerBase::Any,
                    max_run: 3,
                },
            },
        ],
        repair_policy: RepairPolicy::default(),
        digest: None,
    }
}

fn assembled(g: &Genotype) -> Vec<u8> {
    let mut out = FIVE_PRIME.as_bytes().to_vec();
    for codon in g.codons() {
        out.extend_from_slice(&codon.to_bytes());
    }
    out.extend_from_slice(THREE_PRIME.as_bytes());
    out
}

// --------------------------------------------------------------------------- //
// TST-SRCH-006 — the acceptance criterion
// --------------------------------------------------------------------------- //

/// Protein and lock invariance after 10⁶ randomized operator applications.
///
/// Runs against a genotype with a locked span, because an operator that ignored
/// locks would still preserve the protein and would pass a protein-only check.
/// Every mutation operator participates, including the two that need a context
/// model, so the seam is covered rather than skipped.
#[test]
fn tst_srch_006_protein_and_lock_invariance_under_a_million_applications() {
    let start = Genotype::new(locked_sequence());
    let protein = start.sequence().protein_string();
    let locked: Vec<Codon> = (2..5).map(|i| start.codons()[i]).collect();
    let weights = weights();
    let mix = OperatorMix::with_context_model(Arc::new(NeighbourModel));
    assert_eq!(mix.len(), 4, "all four operators participate");

    let checker = ConstraintChecker::compile(&constraint_set()).unwrap();
    let mut rng = SearchRng::new(20260826);
    let mut current = start.clone();
    let mut applications = 0u64;
    let mut changes = 0u64;

    // 10^6 applications, per the criterion. Each is an operator proposal plus
    // its application, so the loop count is the application count.
    while applications < 1_000_000 {
        let construct = assembled(&current);
        let view = ConstructView::new(
            &construct,
            vec![(FIVE_PRIME.len(), FIVE_PRIME.len() + CDS.len())],
        );
        let found = checker.check(&view);
        let actionable: Vec<ActionableViolation<'_>> =
            found.iter().filter_map(ActionableViolation::new).collect();

        let ctx = MutationContext::new(&current, &weights).with_violations(&actionable);
        let mv = mix.propose(&ctx, &mut rng).unwrap();
        let before = current.codons().to_vec();
        current.apply(&mv).unwrap();
        applications += 1;
        if current.codons() != before.as_slice() {
            changes += 1;
        }

        // Checked every application, not at the end: a violation that was
        // introduced and then reverted would be invisible to a final check.
        assert_eq!(
            current.sequence().protein_string(),
            protein,
            "protein changed after {applications} applications"
        );
        for (offset, codon) in locked.iter().enumerate() {
            assert_eq!(
                current.codons()[2 + offset],
                *codon,
                "locked codon {} changed after {applications} applications",
                2 + offset
            );
        }
    }

    assert!(
        changes * 4 > applications,
        "only {changes} of {applications} applications changed anything; the invariance \
         assertions are close to vacuous"
    );
    assert!(current.hash_is_consistent());
    assert!(!current.same_codons_as(&start), "the walk went nowhere");
}

/// Crossover preserves both invariants too.
///
/// Separate from the mutation walk because crossover is the operator most
/// likely to write a codon directly — it is copying slices between parents —
/// and it is the one §5.2's pseudocode guards with a "restore locked positions
/// afterwards" pass.
#[test]
fn tst_srch_006_crossover_preserves_protein_and_locks() {
    let weights = weights();
    let mix = OperatorMix::without_context_model();
    let mut rng = SearchRng::new(4242);

    let mut a = Genotype::new(locked_sequence());
    let mut b = Genotype::new(locked_sequence());
    let protein = a.sequence().protein_string();
    let locked: Vec<Codon> = (2..5).map(|i| a.codons()[i]).collect();

    // Diverge the parents first, or crossover has nothing to exchange.
    for _ in 0..200 {
        let ctx = MutationContext::new(&a, &weights);
        let mv = mix.propose(&ctx, &mut rng).unwrap();
        a.apply(&mv).unwrap();
        let ctx = MutationContext::new(&b, &weights);
        let mv = mix.propose(&ctx, &mut rng).unwrap();
        b.apply(&mv).unwrap();
    }
    assert!(!a.same_codons_as(&b));

    let operators: Vec<Box<dyn CrossoverOperator>> = vec![
        Box::new(BlockUniform::default()),
        Box::new(OnePoint),
        Box::new(TwoPoint),
    ];
    let mut exchanges = 0usize;
    for operator in &operators {
        for _ in 0..20_000 {
            let (x, y) = operator.cross(&a, &b, &mut rng).unwrap();
            for child in [&x, &y] {
                assert_eq!(child.sequence().protein_string(), protein);
                for (offset, codon) in locked.iter().enumerate() {
                    assert_eq!(child.codons()[2 + offset], *codon);
                }
                assert!(child.hash_is_consistent());
            }
            if !x.same_codons_as(&a) {
                exchanges += 1;
            }
        }
    }
    assert!(
        exchanges > 10_000,
        "only {exchanges} crossovers exchanged anything"
    );
}

// --------------------------------------------------------------------------- //
// interop with the exact-delta path
// --------------------------------------------------------------------------- //

/// An operator move is the shape the landscape's delta path consumes.
///
/// Not "can be converted into" — is. `Move` carries `cf_landscape_api::Edit`
/// values, so `to_delta_target` moves them rather than translating them. If the
/// genotype had recorded only its final codons, a caller would have had to
/// recover the edit list with an `O(L)` diff, which would put the query back on
/// the expensive path by construction — the failure this asserts against.
#[test]
fn tst_srch_010_every_operator_move_is_a_delta_query() {
    let weights = weights();
    let mix = OperatorMix::with_context_model(Arc::new(NeighbourModel));
    let mut rng = SearchRng::new(77);
    let base = genotype();
    let digest = base.sequence().digest().unwrap();

    let mut singles = 0usize;
    let mut kinds = std::collections::BTreeSet::new();

    for _ in 0..20_000 {
        let ctx = MutationContext::new(&base, &weights);
        let mv = mix.propose(&ctx, &mut rng).unwrap().compacted(&base);
        if mv.is_empty() {
            continue;
        }
        kinds.insert(mv.kind());
        if mv.is_single_substitution() {
            singles += 1;
        }

        let target = mv.to_delta_target(digest.clone(), None).unwrap();
        match target {
            Target::Delta {
                edits, base_digest, ..
            } => {
                assert_eq!(base_digest, digest);
                assert_eq!(edits.as_slice(), mv.edits());
                // Every edit is applicable to the base it was proposed against.
                assert!(base.with_move(&mv).is_ok());
            }
            other => panic!("expected a delta target, got {}", other.kind()),
        }
    }

    assert!(
        singles > 5_000,
        "only {singles} single-substitution moves; tier 2's accept/reject shape is barely covered"
    );
    assert!(
        kinds.contains(&MoveKind::UniformSynonymous) && kinds.contains(&MoveKind::PottsGuided),
        "the mix produced only {kinds:?}"
    );
}

/// A crossover child can be scored as a delta against its parent.
///
/// Block-uniform crossover exchanges whole blocks, so the edit count is
/// usually large — but the edit list exists, which is what lets `P1-24` decide
/// per child whether a delta or a rescore is cheaper instead of having the
/// choice made for it.
#[test]
fn tst_srch_010_a_crossover_child_carries_its_edit_list() {
    let weights = weights();
    let mix = OperatorMix::without_context_model();
    let mut rng = SearchRng::new(5150);

    let a = genotype();
    let mut b = genotype();
    for _ in 0..50 {
        let ctx = MutationContext::new(&b, &weights);
        let mv = mix.propose(&ctx, &mut rng).unwrap();
        b.apply(&mv).unwrap();
    }

    let (child, _) = BlockUniform::default().cross(&a, &b, &mut rng).unwrap();
    let mv = a.diff(&child).unwrap();
    assert!(a.with_move(&mv).unwrap().same_codons_as(&child));
    if !mv.is_empty() {
        assert!(mv.to_delta_target(Digest::of_bytes(b"base"), None).is_ok());
    }
}

// --------------------------------------------------------------------------- //
// the actionable/non-actionable boundary
// --------------------------------------------------------------------------- //

/// `violation_directed` only ever proposes edits inside the coding region.
///
/// The construct's 3' flank contains a BamHI site, which is a real blocking
/// violation and entirely non-actionable: it lies in the backbone, so no codon
/// index maps to it. `DEP-CONS-002` says such a violation belongs to the
/// feasibility analyser, not to the search, and the operator's signature
/// enforces it — `ActionableViolation` cannot be built from it.
#[test]
fn p1_23_violation_directed_cannot_touch_the_backbone() {
    let checker = ConstraintChecker::compile(&constraint_set()).unwrap();
    let g = genotype();
    let construct = assembled(&g);
    let view = ConstructView::new(
        &construct,
        vec![(FIVE_PRIME.len(), FIVE_PRIME.len() + CDS.len())],
    );
    let found = checker.check(&view);

    let backbone: Vec<_> = found.iter().filter(|v| !v.actionable).collect();
    assert!(
        !backbone.is_empty(),
        "fixture has no non-actionable violation, so it cannot show the boundary"
    );
    for violation in &backbone {
        assert!(
            ActionableViolation::new(violation).is_none(),
            "a backbone-confined violation must not become an ActionableViolation"
        );
        assert!(
            g.codons_in_span(violation.start_nt, violation.end_nt)
                .is_empty(),
            "a backbone violation maps to no codon"
        );
    }

    // And with only non-actionable violations in hand, the operator declines.
    let actionable: Vec<ActionableViolation<'_>> =
        found.iter().filter_map(ActionableViolation::new).collect();
    let weights = weights();
    let ctx = MutationContext::new(&g, &weights).with_violations(&[]);
    assert!(!ViolationDirected.applicable(&ctx));

    // When there is an actionable one, every proposed edit is inside the CDS
    // and unlocked.
    if !actionable.is_empty() {
        let ctx = MutationContext::new(&g, &weights).with_violations(&actionable);
        let mut rng = SearchRng::new(3);
        let mutable = g.mutable_positions();
        for _ in 0..5_000 {
            let mv = ViolationDirected.propose(&ctx, &mut rng).unwrap();
            for edit in mv.edits() {
                assert!(mutable.contains(&edit.codon_index));
            }
        }
    }
}

/// A locked codon inside a violation span is not proposed.
#[test]
fn p1_23_a_locked_position_inside_a_violation_span_is_skipped() {
    let g = Genotype::new(locked_sequence());
    let weights = weights();
    let mutable = g.mutable_positions();
    assert!(!mutable.contains(&3), "codon 3 is inside the locked span");

    // A violation covering the locked span and one codon either side.
    let (lo, _) = g.codon_span_nt(1);
    let (_, hi) = g.codon_span_nt(5);
    let violation = cf_constraints::Violation {
        constraint_id: "synthetic".into(),
        class: ConstraintClass::Blocking,
        start_nt: lo,
        end_nt: hi,
        detail: "covers the locked span".into(),
        actionable: true,
    };
    let actionable = [ActionableViolation::new(&violation).unwrap()];
    let ctx = MutationContext::new(&g, &weights).with_violations(&actionable);

    let mut rng = SearchRng::new(11);
    let mut seen = std::collections::BTreeSet::new();
    for _ in 0..20_000 {
        let mv = ViolationDirected.propose(&ctx, &mut rng).unwrap();
        for edit in mv.edits() {
            seen.insert(edit.codon_index);
            assert!(
                !(2..5).contains(&edit.codon_index),
                "proposed an edit at locked codon {}",
                edit.codon_index
            );
        }
    }
    assert!(
        seen.len() >= 2,
        "only {} distinct positions proposed; the assertion is weak",
        seen.len()
    );
}

// --------------------------------------------------------------------------- //
// operator availability
// --------------------------------------------------------------------------- //

/// Without a context model, two operators are absent and say so.
#[test]
fn p1_23_missing_operators_are_reported_not_approximated() {
    let mix = OperatorMix::without_context_model();
    assert_eq!(
        mix.names(),
        vec!["uniform_synonymous", "violation_directed"]
    );

    let unavailable: Vec<&str> = mix.unavailable().iter().map(|u| u.name).collect();
    assert_eq!(unavailable, vec!["potts_guided", "block_gibbs"]);

    let warnings = mix.warnings();
    assert_eq!(warnings.len(), 2);
    for warning in &warnings {
        assert_eq!(warning.code, "CF-SRCH-001");
        assert!(warning.detail.contains("P3-05"), "{}", warning.detail);
        assert!(warning.detail.contains("renormalized"));
    }

    // With a model, nothing is missing.
    let full = OperatorMix::with_context_model(Arc::new(NeighbourModel));
    assert!(full.unavailable().is_empty());
    assert!(full.warnings().is_empty());
}

/// A feasible individual in a violation-only mix gets an error, not an empty
/// move.
#[test]
fn p1_23_an_inapplicable_mix_fails_loudly() {
    let mut mix = OperatorMix::empty();
    mix.register(Box::new(ViolationDirected), 1.0);
    let g = genotype();
    let weights = weights();
    let ctx = MutationContext::new(&g, &weights);
    let mut rng = SearchRng::new(1);

    let err = mix.propose(&ctx, &mut rng).unwrap_err();
    assert_eq!(err.code, "CF-SRCH-002");

    assert_eq!(
        OperatorMix::empty()
            .propose(&ctx, &mut rng)
            .unwrap_err()
            .code,
        "CF-SRCH-002"
    );
}

/// `potts_guided` actually consults its model.
///
/// A `PottsGuided` that ignored the conditional and drew uniformly would pass
/// every invariance test in this file, so the distribution is checked directly:
/// the model's conditional is strongly non-uniform, and the observed frequencies
/// have to follow it.
#[test]
fn p1_23_potts_guided_samples_from_the_model() {
    let g = genotype();
    let weights = weights();
    let ctx = MutationContext::new(&g, &weights);
    let operator = PottsGuided::new(Arc::new(NeighbourModel)).with_temperature(1.0);

    // Leucine at codon 1: six synonyms, so the conditional has room to differ.
    let index = 1;
    assert_eq!(
        GeneticCode::get(1)
            .unwrap()
            .synonyms(g.codons()[index])
            .len(),
        6
    );
    let mut expected = Vec::new();
    NeighbourModel
        .conditional(g.sequence(), index, 1.0, &mut expected)
        .unwrap();
    let current = g.codons()[index];
    expected.retain(|(c, _)| *c != current);
    let total: f64 = expected.iter().map(|(_, p)| p).sum();

    let mut counts = std::collections::BTreeMap::new();
    let mut rng = SearchRng::new(999);
    let mut draws = 0usize;
    while draws < 60_000 {
        let mv = operator.propose(&ctx, &mut rng).unwrap();
        for edit in mv.edits() {
            if edit.codon_index == index {
                *counts.entry(edit.codon).or_insert(0usize) += 1;
                draws += 1;
            }
        }
    }

    for (codon, probability) in &expected {
        let observed = *counts.get(codon).unwrap_or(&0) as f64 / draws as f64;
        let target = probability / total;
        assert!(
            (observed - target).abs() < 0.02,
            "codon {} observed at {observed:.4}, model says {target:.4}",
            codon.to_string_upper()
        );
    }
}

// --------------------------------------------------------------------------- //
// tabu consultation (the part of P1-18 that does not need P1-18b)
// --------------------------------------------------------------------------- //

/// A tabu edit is suppressed; the leak lets it through at the stated rate.
///
/// The generation clock that calls `TABU_TICK` is `P1-18b`, so the ledger is
/// advanced by hand here. What is tested is the consultation itself, which is
/// the part of the contract `P1-23` owns.
#[test]
fn p1_23_the_mutation_operator_consults_the_tabu_ledger() {
    let g = genotype();
    let weights = weights();
    let index = 1;
    let current = g.codons()[index];
    let synonyms: Vec<Codon> = g
        .sequence()
        .synonyms_at(index)
        .iter()
        .copied()
        .filter(|c| *c != current)
        .collect();
    assert_eq!(synonyms.len(), 5);

    // No leak: a suppressed codon must never appear.
    let mut ledger = TabuLedger::new(0.0);
    for codon in synonyms.iter().take(4) {
        ledger.record(0, index, *codon, 10);
    }
    let ctx = MutationContext::new(&g, &weights).with_tabu(&ledger);
    let mut rng = SearchRng::new(31337);

    let mut seen = std::collections::BTreeSet::new();
    for _ in 0..20_000 {
        let mv = UniformSynonymous.propose(&ctx, &mut rng).unwrap();
        for edit in mv.edits() {
            if edit.codon_index == index {
                seen.insert(edit.codon);
            }
        }
    }
    assert_eq!(
        seen.len(),
        1,
        "expected only the untabooed synonym, saw {seen:?}"
    );
    assert!(seen.contains(&synonyms[4]));

    // The ledger ages out; after enough ticks the suppression lifts.
    for _ in 0..10 {
        ledger.tick();
    }
    assert!(!ledger.is_tabu(0, index, synonyms[0]));

    let ctx = MutationContext::new(&g, &weights).with_tabu(&ledger);
    let mut seen_after = std::collections::BTreeSet::new();
    for _ in 0..20_000 {
        let mv = UniformSynonymous.propose(&ctx, &mut rng).unwrap();
        for edit in mv.edits() {
            if edit.codon_index == index {
                seen_after.insert(edit.codon);
            }
        }
    }
    assert_eq!(seen_after.len(), 5, "every synonym is available again");
}

/// A fully tabu position yields nothing rather than a suppressed edit.
#[test]
fn p1_23_a_fully_suppressed_position_proposes_nothing_there() {
    let g = genotype();
    let weights = weights();
    let index = 1;
    let mut ledger = TabuLedger::new(0.0);
    for codon in g.sequence().synonyms_at(index) {
        ledger.record(0, index, *codon, 5);
    }
    let ctx = MutationContext::new(&g, &weights).with_tabu(&ledger);
    let mut rng = SearchRng::new(8);

    for _ in 0..20_000 {
        let mv = ViolationDirected.propose(&ctx, &mut rng).unwrap();
        assert!(mv.is_empty() || mv.edits().iter().all(|e| e.codon_index != index));
    }
}

// --------------------------------------------------------------------------- //
// duplicate suppression and seed injection
// --------------------------------------------------------------------------- //

#[test]
fn p1_23_duplicates_are_suppressed_and_collisions_are_confirmed() {
    let weights = weights();
    let mix = OperatorMix::without_context_model();
    let mut rng = SearchRng::new(606);
    let parent = genotype();

    let mut population = Population::new();
    assert!(population.insert(parent.clone()));
    assert!(
        !population.insert(parent.clone()),
        "an exact copy is refused"
    );

    let mut forced = 0usize;
    for _ in 0..2_000 {
        let ctx = MutationContext::new(&parent, &weights);
        let offspring = population
            .breed_distinct(&parent, &mix, &ctx, &mut rng)
            .unwrap();
        if offspring.forced {
            forced += 1;
        }
        population.insert(offspring.genotype);
    }
    assert!(population.len() > 1_000, "population barely grew");

    // Every member is distinct by full comparison, not just by hash.
    for (i, a) in population.members().iter().enumerate() {
        for b in population.members().iter().skip(i + 1) {
            assert!(!a.same_codons_as(b));
        }
    }
    assert!(forced < 2_000, "every single offspring needed forcing");
}

#[test]
fn p1_23_seed_injection_respects_the_cap_and_reports_what_it_did() {
    let weights = weights();
    let mix = OperatorMix::without_context_model();
    let mut rng = SearchRng::new(1234);

    let base = genotype();
    let mut seeds = Vec::new();
    for n in 0..10 {
        let mut seed = base.clone();
        for _ in 0..(n + 1) {
            let ctx = MutationContext::new(&seed, &weights);
            let mv = mix.propose(&ctx, &mut rng).unwrap();
            seed.apply(&mv).unwrap();
        }
        seeds.push(seed);
    }

    let mut population = Population::new();
    let report = population
        .inject_seeds(&seeds, 20, 0.2, &mix, &weights, &mut rng)
        .unwrap();

    assert_eq!(report.offered, 10);
    assert_eq!(report.accepted, 4, "20 * 0.2 = 4");
    assert_eq!(report.capped, 6);
    assert_eq!(report.mutated_copies, 8, "two mutated copies per seed");
    assert_eq!(population.len(), report.accepted + report.mutated_copies);

    let bad = population.inject_seeds(&seeds, 20, 1.5, &mix, &weights, &mut rng);
    assert_eq!(bad.unwrap_err().code, "CF-SRCH-004");
}

// --------------------------------------------------------------------------- //
// determinism
// --------------------------------------------------------------------------- //

/// The same seed produces the same walk.
///
/// `GATE-DET-001` eventually demands bit-identical output across thread counts.
/// That is `P1-33`, but it is only reachable if the operators are deterministic
/// given their stream, so the property is pinned here rather than discovered to
/// be missing later.
#[test]
fn p1_23_operators_are_deterministic_given_a_seed() {
    let weights = weights();
    let walk = |seed: u64| {
        let mix = OperatorMix::with_context_model(Arc::new(NeighbourModel));
        let mut rng = SearchRng::new(seed);
        let mut g = genotype();
        let mut trace = Vec::new();
        for _ in 0..2_000 {
            let ctx = MutationContext::new(&g, &weights);
            let mv = mix.propose(&ctx, &mut rng).unwrap();
            trace.push((mv.kind(), mv.touched()));
            g.apply(&mv).unwrap();
        }
        (g.hash(), trace)
    };

    let (hash_a, trace_a) = walk(2026);
    let (hash_b, trace_b) = walk(2026);
    assert_eq!(hash_a, hash_b);
    assert_eq!(trace_a, trace_b);

    let (hash_c, _) = walk(2027);
    assert_ne!(hash_a, hash_c, "different seeds must diverge");
}

/// A move applied and reverted returns the genotype exactly, hash included.
#[test]
fn p1_23_a_move_is_reversible() {
    let mut g = genotype();
    let before = g.clone();
    let mv = Move::new(
        MoveKind::Explicit,
        vec![
            Edit::new(1, Codon::parse("CTC").unwrap()),
            Edit::new(3, Codon::parse("TCT").unwrap()),
        ],
    );
    g.apply(&mv).unwrap();
    let back = g.diff(&before).unwrap();
    g.apply(&back).unwrap();
    assert!(g.same_codons_as(&before));
    assert_eq!(g.hash(), before.hash());
}
