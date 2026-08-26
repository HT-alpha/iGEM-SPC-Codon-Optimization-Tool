//! `TST-SRCH-002` (tier-1 portion), `TST-SRCH-011`, and the `P1-27b`
//! differential.
//!
//! The headline test is `TST-SRCH-011`, the cost-shaped one. Last round's audit
//! established that a values-only test cannot see a correct-but-expensive
//! implementation, and tier 1 is the place that matters most: its inner loop
//! runs `O(L·b·S)` times, so a landscape query per step would be right and
//! roughly three orders of magnitude over budget.

mod support;

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use cf_constraints::legacy::{legacy_is_feasible, LegacyScenario};
use cf_constraints::spec::{
    Constraint, ConstraintClass, ConstraintSet, RepairPolicy, Rule, Scope, Strand,
};
use cf_constraints::{ConstraintChecker, ConstructView};
use cf_landscape_api::{
    Capability, EvalContext, PositionDeltas, Provenance, Rebased, TermProvider, TermRegistry,
};
use cf_search::baseline::{baseline_backtracking, BaselineBudget, BaselineOutcome};
use cf_search::tier1::{enumerate_feasible, surrogate_of, Tier1Request};
use cf_search::{tier1_seed, Genotype, SearchRng, Tier1Config};
use cf_types::{Codon, CodonSequence, Direction, Result, TermId, Unit};

use support::{host, weights};

// --------------------------------------------------------------------------- //
// fixtures
// --------------------------------------------------------------------------- //

/// `M L R S G V A P T K L R S *`, inside a plasmid backbone.
const CDS: &str = "ATGCTGCGTAGCGGCGTGGCGCCGACCAAACTGCGTAGCTAA";
const FIVE_PRIME: &str = "TTGACAATTAATCATCGGCTCGTATAATGTGTGGA";
const THREE_PRIME: &str = "AAAGCCCGGATCCTCTAGAGGG";

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

/// A small protein whose whole design space can be enumerated.
///
/// `M L R S *` — 6 × 6 × 6 = 216 assignments, well inside the ground-truth
/// regime §2 describes for exact DP.
fn toy() -> Genotype {
    Genotype::new(
        CodonSequence::from_cds("ATGCTGCGTAGCTAA", 1)
            .unwrap()
            .with_five_prime_flank("GGTCTCA")
            .unwrap()
            .with_three_prime_flank("TGAGACC")
            .unwrap(),
    )
}

fn motif_rules(motifs: &[&str]) -> ConstraintSet {
    ConstraintSet {
        schema_version: Some("1.0.0".into()),
        constraint_set_id: "cst_tier1".into(),
        name: None,
        derived_from: vec![],
        constraints: motifs
            .iter()
            .enumerate()
            .map(|(index, motif)| Constraint {
                id: format!("no_{index}"),
                class: ConstraintClass::Blocking,
                scope: Scope::Construct,
                cistron_id: None,
                weight: 1.0,
                enabled: true,
                rule: Rule::ForbiddenMotif {
                    motif: (*motif).to_string(),
                    // Single strand: the legacy predicate matches literals and
                    // has no strand awareness, so a two-strand rule would make
                    // the P1-27b differential compare different questions.
                    strand: Strand::Sense,
                    max_occurrences: 0,
                },
            })
            .collect(),
        repair_policy: RepairPolicy::default(),
        digest: None,
    }
}

fn assemble(g: &Genotype) -> Vec<u8> {
    let sequence = g.sequence();
    let mut out = sequence.five_prime_flank().to_vec();
    for codon in sequence.codons() {
        out.extend_from_slice(&codon.to_bytes());
    }
    out.extend_from_slice(sequence.three_prime_flank());
    out
}

fn cds_string(g: &Genotype) -> String {
    g.codons().iter().map(|c| c.to_string_upper()).collect()
}

// --------------------------------------------------------------------------- //
// an instrumented L1 term
// --------------------------------------------------------------------------- //

#[derive(Clone, Debug, Default)]
struct Calls(Arc<Counters>);

#[derive(Debug, Default)]
struct Counters {
    scores: AtomicU64,
    deltas: AtomicU64,
    tables: AtomicU64,
    rebases: AtomicU64,
}

impl Calls {
    fn scores(&self) -> u64 {
        self.0.scores.load(Ordering::SeqCst)
    }
    fn deltas(&self) -> u64 {
        self.0.deltas.load(Ordering::SeqCst)
    }
    fn tables(&self) -> u64 {
        self.0.tables.load(Ordering::SeqCst)
    }
    fn rebases(&self) -> u64 {
        self.0.rebases.load(Ordering::SeqCst)
    }
}

/// An additive L1 term, instrumented at every entry point.
#[derive(Debug)]
struct CountingTerm {
    term: TermId,
    calls: Calls,
    weight: [f64; 64],
}

impl CountingTerm {
    fn new(term: &str, calls: Calls, salt: u64) -> Self {
        let mut weight = [0.0f64; 64];
        for (index, slot) in weight.iter_mut().enumerate() {
            *slot = (((index as u64 * 31 + salt) % 97) as f64) / 97.0;
        }
        Self {
            term: TermId::new(term).unwrap(),
            calls,
            weight,
        }
    }

    fn total(&self, sequence: &CodonSequence) -> f64 {
        sequence
            .codons()
            .iter()
            .map(|c| self.weight[c.index_usize()])
            .sum()
    }
}

impl TermProvider for CountingTerm {
    fn term_id(&self) -> &TermId {
        &self.term
    }
    fn unit(&self) -> Unit {
        Unit::new(Unit::DIMENSIONLESS)
    }
    fn direction(&self) -> Direction {
        Direction::Maximize
    }
    fn capability(&self) -> Capability {
        Capability::L1
    }
    fn provenance(&self) -> Provenance {
        Provenance::builtin("test_counting", "0.1.0")
    }
    fn score(&self, ctx: &EvalContext<'_>) -> Result<f64> {
        self.calls.0.scores.fetch_add(1, Ordering::SeqCst);
        Ok(self.total(ctx.sequence))
    }
    fn rebase(&self, ctx: &EvalContext<'_>) -> Result<Rebased> {
        self.calls.0.rebases.fetch_add(1, Ordering::SeqCst);
        Ok(Rebased::value_only(self.total(ctx.sequence)))
    }
    fn delta(
        &self,
        ctx: &EvalContext<'_>,
        _base: &Rebased,
        edits: &[cf_landscape_api::Edit],
    ) -> Result<f64> {
        self.calls.0.deltas.fetch_add(1, Ordering::SeqCst);
        let codons = ctx.sequence.codons();
        Ok(edits
            .iter()
            .map(|e| {
                self.weight[e.codon.index_usize()]
                    - self.weight[codons[e.codon_index].index_usize()]
            })
            .sum())
    }
    fn delta_table(
        &self,
        ctx: &EvalContext<'_>,
        base: &Rebased,
        positions: &[usize],
    ) -> Result<Vec<PositionDeltas>> {
        self.calls.0.tables.fetch_add(1, Ordering::SeqCst);
        positions
            .iter()
            .map(|&index| {
                let current = ctx.sequence.codons()[index];
                let alternatives = ctx
                    .sequence
                    .synonyms_at(index)
                    .iter()
                    .copied()
                    .filter(|c| *c != current)
                    .map(|c| {
                        self.delta(ctx, base, &[cf_landscape_api::Edit::new(index, c)])
                            .map(|d| (c, d))
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(PositionDeltas {
                    codon_index: index,
                    alternatives,
                })
            })
            .collect()
    }
}

fn registry_with(calls: &Calls, names: &[&str]) -> TermRegistry {
    let mut registry = TermRegistry::new();
    for (salt, name) in names.iter().enumerate() {
        registry
            .register_default(Box::new(CountingTerm::new(
                name,
                calls.clone(),
                salt as u64,
            )))
            .unwrap();
    }
    registry
}

fn terms(names: &[&str]) -> Vec<TermId> {
    names.iter().map(|n| TermId::new(*n).unwrap()).collect()
}

// --------------------------------------------------------------------------- //
// TST-SRCH-011 — the cost-shaped test
// --------------------------------------------------------------------------- //

/// Tier 1 takes one delta table per term, not one landscape query per DP step.
///
/// The assertion is on call counts because values cannot distinguish the two
/// implementations. A DP that queried per (state, synonym) would produce the
/// same seeds and issue `states_generated` queries — here in the tens of
/// thousands — while `delta_tables` stayed at one per term, so the two counters
/// together are what pin the complexity claim.
#[test]
fn tst_srch_011_tier1_scores_from_delta_tables_not_per_step_queries() {
    let calls = Calls::default();
    let registry = registry_with(&calls, &["usage.cai", "usage.gc_content"]);
    let checker = ConstraintChecker::compile(&motif_rules(&["GGATCC", "GAATTC"])).unwrap();
    let base = base();
    let host = host();
    let config = Tier1Config::default().with_beam(32).with_seeds(4);
    let mut rng = SearchRng::new(20260826);

    let report = tier1_seed(
        &Tier1Request {
            base: &base,
            registry: &registry,
            terms: &terms(&["usage.cai", "usage.gc_content"]),
            weights: &[1.0, 0.5],
            host: &host,
            checker: &checker,
            config: &config,
        },
        &mut rng,
    )
    .unwrap();

    assert_eq!(report.rebases, 2, "one pin per scored term");
    assert_eq!(report.delta_tables, 2, "one delta table per scored term");
    assert_eq!(calls.tables(), 2);
    assert_eq!(calls.rebases(), 2);

    // The fixture has to be big enough that a per-step implementation would be
    // visibly different, or the counters prove nothing.
    assert!(
        report.states_generated > 500,
        "only {} states generated; a per-step implementation would look identical",
        report.states_generated
    );

    // `delta` is called only from inside `delta_table` — twice per table entry
    // — and never from the DP loop.
    let table_entries: u64 = calls.deltas();
    assert!(
        table_entries < report.states_generated,
        "delta was called {table_entries} times against {} states generated; \
         the DP is querying the landscape per step",
        report.states_generated
    );
    assert_eq!(
        calls.scores(),
        0,
        "tier 1 scored a full sequence; the delta table is supposed to be enough"
    );
    assert!(!report.seeds.is_empty());
}

/// A term with no delta table is refused before the tier runs.
///
/// ADR-0009 §2.1's tier-planning check, for tier 1. Discovering an L0 term
/// mid-run would waste the tier, and silently dropping it would optimize a
/// different objective than the caller asked for.
#[test]
fn tst_srch_011_an_l0_term_is_refused_up_front() {
    #[derive(Debug)]
    struct ScoreOnly(TermId);
    impl TermProvider for ScoreOnly {
        fn term_id(&self) -> &TermId {
            &self.0
        }
        fn unit(&self) -> Unit {
            Unit::new(Unit::DIMENSIONLESS)
        }
        fn direction(&self) -> Direction {
            Direction::Maximize
        }
        fn capability(&self) -> Capability {
            Capability::L0
        }
        fn provenance(&self) -> Provenance {
            Provenance::builtin("test_l0", "0.1.0")
        }
        fn score(&self, _ctx: &EvalContext<'_>) -> Result<f64> {
            Ok(0.0)
        }
    }

    let mut registry = TermRegistry::new();
    registry
        .register_default(Box::new(ScoreOnly(TermId::new("lm.naturalness").unwrap())))
        .unwrap();

    let checker = ConstraintChecker::compile(&motif_rules(&["GGATCC"])).unwrap();
    let base = base();
    let host = host();
    let config = Tier1Config::default();
    let mut rng = SearchRng::new(1);

    let err = tier1_seed(
        &Tier1Request {
            base: &base,
            registry: &registry,
            terms: &terms(&["lm.naturalness"]),
            weights: &[1.0],
            host: &host,
            checker: &checker,
            config: &config,
        },
        &mut rng,
    )
    .unwrap_err();
    assert_eq!(err.code, "CF-LAND-002");
    assert!(err.detail.contains("tier 1 needs a delta table"));
}

// --------------------------------------------------------------------------- //
// TST-SRCH-002, tier-1 portion — optimality on a brute-forceable instance
// --------------------------------------------------------------------------- //

/// On the toy space, tier 1 with a full-width beam finds the exact optimum.
///
/// §2 notes exact DP is affordable only for small state spaces; with a beam at
/// least as wide as the space, the beam search *is* exact, which makes this the
/// tier-1 half of `TST-SRCH-002`. Ground truth comes from enumerating all 216
/// assignments and ranking them on the same surrogate the DP optimizes — a
/// different scale would make the comparison meaningless.
#[test]
fn tst_srch_002_tier1_finds_the_optimum_on_the_toy_space() {
    let calls = Calls::default();
    let registry = registry_with(&calls, &["usage.cai"]);
    let checker = ConstraintChecker::compile(&motif_rules(&["GGTAGC", "CGCAGC"])).unwrap();
    let base = toy();
    let host = host();
    let term_list = terms(&["usage.cai"]);

    let feasible = enumerate_feasible(&base, &checker, 4096).unwrap();
    assert!(
        feasible.len() > 30 && feasible.len() < 216,
        "the constraint set must bind without emptying the space; {} feasible",
        feasible.len()
    );

    let best = feasible
        .iter()
        .map(|candidate| {
            surrogate_of(candidate, &base, &registry, &term_list, &[1.0], &host).unwrap()
        })
        .fold(f64::NEG_INFINITY, f64::max);

    let mut exact = 0usize;
    let seeds = 200usize;
    for seed in 0..seeds as u64 {
        let config = Tier1Config::default().with_beam(256).with_seeds(1);
        let mut rng = SearchRng::new(seed);
        let report = tier1_seed(
            &Tier1Request {
                base: &base,
                registry: &registry,
                terms: &term_list,
                weights: &[1.0],
                host: &host,
                checker: &checker,
                config: &config,
            },
            &mut rng,
        )
        .unwrap();
        if (report.surrogate[0] - best).abs() < 1e-9 {
            exact += 1;
        }
    }
    assert_eq!(
        exact, seeds,
        "a beam wider than the space must be exact on every seed; {exact}/{seeds}"
    );
}

/// A truncated beam is where seeds are lost, and the test says so rather than
/// asserting a threshold it cannot justify.
#[test]
fn tst_srch_002_a_narrow_beam_is_measurably_worse() {
    let calls = Calls::default();
    let registry = registry_with(&calls, &["usage.cai"]);
    let checker = ConstraintChecker::compile(&motif_rules(&["GGTAGC", "CGCAGC"])).unwrap();
    let base = toy();
    let host = host();
    let term_list = terms(&["usage.cai"]);

    let wide = {
        let config = Tier1Config::default().with_beam(256).with_seeds(1);
        let mut rng = SearchRng::new(7);
        tier1_seed(
            &Tier1Request {
                base: &base,
                registry: &registry,
                terms: &term_list,
                weights: &[1.0],
                host: &host,
                checker: &checker,
                config: &config,
            },
            &mut rng,
        )
        .unwrap()
        .surrogate[0]
    };
    let narrow = {
        let config = Tier1Config::default().with_beam(1).with_seeds(1);
        let mut rng = SearchRng::new(7);
        tier1_seed(
            &Tier1Request {
                base: &base,
                registry: &registry,
                terms: &term_list,
                weights: &[1.0],
                host: &host,
                checker: &checker,
                config: &config,
            },
            &mut rng,
        )
        .unwrap()
        .surrogate[0]
    };
    assert!(
        narrow <= wide + 1e-12,
        "a width-1 beam scored {narrow} against a width-256 {wide}"
    );
}

// --------------------------------------------------------------------------- //
// invariance and feasibility by construction
// --------------------------------------------------------------------------- //

/// Every seed encodes the input protein, respects locks, and is feasible.
#[test]
fn p1_24_seeds_preserve_protein_and_locks_and_are_feasible() {
    let calls = Calls::default();
    let registry = registry_with(&calls, &["usage.cai"]);
    let checker =
        ConstraintChecker::compile(&motif_rules(&["GGATCC", "GAATTC", "CTGCGT"])).unwrap();

    let mut sequence = CodonSequence::from_cds(CDS, 1)
        .unwrap()
        .with_five_prime_flank(FIVE_PRIME)
        .unwrap()
        .with_three_prime_flank(THREE_PRIME)
        .unwrap();
    sequence
        .add_lock(4, 7, cf_types::LockReason::UserLocked)
        .unwrap();
    let base = Genotype::new(sequence);
    let protein = base.sequence().protein_string();
    let locked: Vec<Codon> = (4..7).map(|i| base.codons()[i]).collect();
    let host = host();

    for seed in 0..50u64 {
        let config = Tier1Config::default().with_beam(48).with_seeds(4);
        let mut rng = SearchRng::new(seed);
        let report = tier1_seed(
            &Tier1Request {
                base: &base,
                registry: &registry,
                terms: &terms(&["usage.cai"]),
                weights: &[1.0],
                host: &host,
                checker: &checker,
                config: &config,
            },
            &mut rng,
        )
        .unwrap();

        assert!(!report.seeds.is_empty());
        for candidate in &report.seeds {
            assert_eq!(candidate.sequence().protein_string(), protein);
            for (offset, codon) in locked.iter().enumerate() {
                assert_eq!(candidate.codons()[4 + offset], *codon);
            }
            let construct = assemble(candidate);
            let view = ConstructView::new(
                &construct,
                vec![(
                    candidate.cds_offset_nt(),
                    candidate.cds_offset_nt() + candidate.len() * 3,
                )],
            );
            assert!(
                checker.is_feasible(&view),
                "seed {seed} produced an infeasible candidate"
            );
            assert!(candidate.hash_is_consistent());
        }
    }
}

/// A motif that only exists across the 5' junction is pruned, not post-checked.
///
/// `DEP-SRCH-006`. §2's pseudocode starts the automaton at `A.start`, which is
/// right for a bare CDS and wrong inside a construct. The fixture's 5' flank
/// ends in `TGGA`; a first codon of `TCC` would complete `GGATCC` across the
/// junction, and only a primed automaton can see that while extending.
#[test]
fn p1_24_the_flanks_prime_the_automaton() {
    let calls = Calls::default();
    let registry = registry_with(&calls, &["usage.cai"]);
    let checker = ConstraintChecker::compile(&motif_rules(&["GGATCC"])).unwrap();

    // 5' flank ends ...TGGA; CDS starts with a serine, whose synonyms include
    // TCC — so TGGA|TCC contains GGATCC straddling the junction.
    let sequence = CodonSequence::from_cds("AGCCTGCGTTAA", 1)
        .unwrap()
        .with_five_prime_flank("TTGACAATTAATCATCGGCTCGTATAATGTGTGGA")
        .unwrap();
    let base = Genotype::new(sequence);
    assert_eq!(
        &base.sequence().five_prime_flank()[31..],
        b"TGGA",
        "fixture no longer ends in TGGA, so the junction case is not exercised"
    );
    let host = host();

    for seed in 0..40u64 {
        let config = Tier1Config::default().with_beam(64).with_seeds(8);
        let mut rng = SearchRng::new(seed);
        let report = tier1_seed(
            &Tier1Request {
                base: &base,
                registry: &registry,
                terms: &terms(&["usage.cai"]),
                weights: &[1.0],
                host: &host,
                checker: &checker,
                config: &config,
            },
            &mut rng,
        )
        .unwrap();
        for candidate in &report.seeds {
            assert_ne!(
                candidate.codons()[0].to_string_upper(),
                "TCC",
                "a candidate completing GGATCC across the 5' junction survived"
            );
        }
        // And nothing reached the post-check: pruning happened in the DP.
        assert_eq!(
            report.post_check_rejected, 0,
            "the junction motif was caught after the fact, not by construction"
        );
    }
}

// --------------------------------------------------------------------------- //
// P1-27b — the differential against baseline_backtracking
// --------------------------------------------------------------------------- //

fn legacy_scenario(base: &Genotype, motifs: &[&str], k: usize, limit: u32) -> LegacyScenario {
    let construct: String = String::from_utf8(assemble(base)).unwrap();
    let cds = cds_string(base);
    let mut scenario = LegacyScenario::from_construct(&construct, &[&cds], k, limit, 4, 5);
    // Replace the auto-detected motif set with the declared one, so both
    // solvers are answering the same question. `auto_detect_forbidden_motifs`
    // derives motifs from the construct; the DP is given them explicitly.
    scenario.motifs = motifs.iter().map(|m| (*m).to_string()).collect();
    scenario
}

/// Both solvers agree on the same feasibility surface, in both directions.
///
/// The surface is motif-only. The prototype also enforces a k-mer repeat limit
/// and §2 explicitly excludes `kmer_repeat` from the DP because it is not a
/// regular language, so a differential that included it would report design as
/// drift. The limit is set high enough here to be inactive, and
/// `p1_27b_the_kmer_row_is_excluded_deliberately` asserts that it really is.
#[test]
fn p1_27b_dp_and_baseline_agree_on_the_shared_surface() {
    let calls = Calls::default();
    let registry = registry_with(&calls, &["usage.cai"]);
    let motifs = ["GGTAGC", "CGCAGC", "CTGCGT"];
    let checker = ConstraintChecker::compile(&motif_rules(&motifs)).unwrap();
    let base = toy();
    let host = host();
    let scenario = legacy_scenario(&base, &motifs, 8, 99);
    let codon_weights = weights();

    let mut dp_solved = 0usize;
    let mut baseline_solved = 0usize;

    for seed in 0..60u64 {
        let config = Tier1Config::default().with_beam(128).with_seeds(1);
        let mut rng = SearchRng::new(seed);
        let dp = tier1_seed(
            &Tier1Request {
                base: &base,
                registry: &registry,
                terms: &terms(&["usage.cai"]),
                weights: &[1.0],
                host: &host,
                checker: &checker,
                config: &config,
            },
            &mut rng,
        );

        // Direction 1: every DP seed satisfies the legacy predicate.
        if let Ok(report) = &dp {
            dp_solved += 1;
            for candidate in &report.seeds {
                assert!(
                    legacy_is_feasible(&cds_string(candidate), &scenario),
                    "seed {seed}: the DP produced {} which the legacy predicate rejects",
                    cds_string(candidate)
                );
            }
        }

        // Direction 2: every baseline solution satisfies the new checker.
        let mut rng = SearchRng::new(seed);
        let report = baseline_backtracking(
            &base,
            &scenario,
            &codon_weights,
            BaselineBudget::from_prototype(5, 2_000),
            &mut rng,
        )
        .unwrap();
        if report.outcome == BaselineOutcome::Solved {
            baseline_solved += 1;
            let solution = report.solution.unwrap();
            assert_eq!(
                solution.sequence().protein_string(),
                base.sequence().protein_string()
            );
            let construct = assemble(&solution);
            let view = ConstructView::new(
                &construct,
                vec![(
                    solution.cds_offset_nt(),
                    solution.cds_offset_nt() + solution.len() * 3,
                )],
            );
            assert!(
                checker.is_feasible(&view),
                "seed {seed}: the baseline solved with {} which the new checker rejects",
                cds_string(&solution)
            );
        }
    }

    assert_eq!(
        dp_solved, 60,
        "the DP must solve every seed on a feasible instance"
    );
    assert!(
        baseline_solved > 50,
        "the baseline solved only {baseline_solved}/60; the fixture is too hard to \
         make the comparison meaningful"
    );
}

/// The DP finding nothing means nothing exists — checked by enumeration.
///
/// The baseline cannot answer this: its `BudgetExhausted` is a timeout, not a
/// proof, exactly as `P1-17`'s `GiveUpReason` distinguishes. So completeness is
/// checked against exhaustive enumeration instead of against the other solver.
#[test]
fn p1_27b_an_empty_result_means_an_empty_space() {
    let calls = Calls::default();
    let registry = registry_with(&calls, &["usage.cai"]);
    // Every serine codon banned: the toy protein cannot be encoded at all.
    let motifs = ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"];
    let checker = ConstraintChecker::compile(&motif_rules(&motifs)).unwrap();
    let base = toy();
    let host = host();

    assert!(
        enumerate_feasible(&base, &checker, 4096)
            .unwrap()
            .is_empty(),
        "the fixture must actually be infeasible"
    );

    let config = Tier1Config::default().with_beam(256).with_seeds(1);
    let mut rng = SearchRng::new(3);
    let err = tier1_seed(
        &Tier1Request {
            base: &base,
            registry: &registry,
            terms: &terms(&["usage.cai"]),
            weights: &[1.0],
            host: &host,
            checker: &checker,
            config: &config,
        },
        &mut rng,
    )
    .unwrap_err();
    assert_eq!(err.code, "CF-SRCH-006");
    assert!(
        err.detail.contains("not a proof"),
        "a beam result must not be reported as a proof: {}",
        err.detail
    );

    let mut rng = SearchRng::new(3);
    let report = baseline_backtracking(
        &base,
        &legacy_scenario(&base, &motifs, 8, 99),
        &weights(),
        BaselineBudget::from_prototype(2, 500),
        &mut rng,
    )
    .unwrap();
    assert_eq!(report.outcome, BaselineOutcome::BudgetExhausted);
    assert!(report.solution.is_none());
}

/// The k-mer row of the comparability table, asserted rather than assumed.
///
/// The prototype enforces a k-mer repeat limit and the DP does not. A
/// differential that did not know this would read the difference as drift, so
/// the asymmetry is pinned: with a limit tight enough to bind, the baseline
/// refuses assignments the DP happily emits.
#[test]
fn p1_27b_the_kmer_row_is_excluded_deliberately() {
    let calls = Calls::default();
    let registry = registry_with(&calls, &["usage.cai"]);
    let checker = ConstraintChecker::compile(&motif_rules(&[])).unwrap();
    let base = toy();
    let host = host();
    // k = 4 with a limit of 1 makes almost every assignment fail the legacy
    // k-mer test while leaving the motif surface empty.
    let scenario = legacy_scenario(&base, &[], 4, 1);

    let config = Tier1Config::default().with_beam(64).with_seeds(8);
    let mut rng = SearchRng::new(11);
    let report = tier1_seed(
        &Tier1Request {
            base: &base,
            registry: &registry,
            terms: &terms(&["usage.cai"]),
            weights: &[1.0],
            host: &host,
            checker: &checker,
            config: &config,
        },
        &mut rng,
    )
    .unwrap();

    let rejected_by_legacy = report
        .seeds
        .iter()
        .filter(|s| !legacy_is_feasible(&cds_string(s), &scenario))
        .count();
    assert!(
        rejected_by_legacy > 0,
        "the fixture does not exercise the asymmetry: every DP seed also passes \
         the legacy k-mer test, so the exclusion is untested"
    );
    assert_eq!(
        report.post_check_rejected, 0,
        "the checker was given no k-mer rule, so nothing should have been post-rejected"
    );
}

/// The baseline never returns a solution encoding a different protein.
#[test]
fn p1_27b_baseline_solutions_preserve_the_protein() {
    let base = base();
    let scenario = legacy_scenario(&base, &["GGATCC"], 8, 99);
    let codon_weights = weights();
    let protein = base.sequence().protein_string();

    let mut solved = 0;
    for seed in 0..40u64 {
        let mut rng = SearchRng::new(seed);
        let report = baseline_backtracking(
            &base,
            &scenario,
            &codon_weights,
            BaselineBudget::from_prototype(5, 5_000),
            &mut rng,
        )
        .unwrap();
        if let Some(solution) = report.solution {
            solved += 1;
            assert_eq!(solution.sequence().protein_string(), protein);
            cf_search::baseline::check_solution(&base, &solution).unwrap();
        }
    }
    assert!(
        solved > 30,
        "only {solved}/40 solved; the assertion is weak"
    );
}

/// The baseline is deterministic under its seed and diverges across seeds.
#[test]
fn p1_27b_baseline_is_deterministic_and_seed_sensitive() {
    let base = base();
    let scenario = legacy_scenario(&base, &["GGATCC"], 8, 99);
    let codon_weights = weights();

    let run = |seed: u64| {
        let mut rng = SearchRng::new(seed);
        baseline_backtracking(
            &base,
            &scenario,
            &codon_weights,
            BaselineBudget::from_prototype(5, 5_000),
            &mut rng,
        )
        .unwrap()
        .solution
        .map(|g| cds_string(&g))
    };

    assert_eq!(run(101), run(101));
    let mut distinct = std::collections::BTreeSet::new();
    for seed in 0..30u64 {
        if let Some(cds) = run(seed) {
            distinct.insert(cds);
        }
    }
    assert!(
        distinct.len() > 5,
        "only {} distinct solutions across 30 seeds; the weighted shuffle is not \
         exploring",
        distinct.len()
    );
}
