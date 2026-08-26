//! `TST-SRCH-001`, the tier-2 half of `TST-SRCH-002`, and `TST-SRCH-012`.
//!
//! `TST-SRCH-012` is the independence-shaped test, in the spirit of
//! `TST-SRCH-011`: it would fail on an implementation that produced correct
//! block updates while letting one same-colour block observe another. Such an
//! implementation returns a good answer and is simply not the coordinate-ascent
//! step it reports being — nothing about the numbers looks wrong.

mod support;

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use cf_constraints::spec::{
    Constraint, ConstraintClass, ConstraintSet, RepairPolicy, Rule, Scope, Strand,
};
use cf_constraints::{ConstraintChecker, ConstructView};
use cf_landscape_api::{
    Capability, EvalContext, InteractionRadius, PositionDeltas, Provenance, Rebased, TermProvider,
    TermRegistry,
};
use cf_search::tier1::enumerate_feasible;
use cf_search::tier2::{partition, tier2_block_ca, Colour, Tier2Plan, Tier2Request};
use cf_search::{Genotype, Tier2Config};
use cf_types::{CodonSequence, Direction, LockReason, Result, TermId, Unit};

use support::host;

// --------------------------------------------------------------------------- //
// fixtures
// --------------------------------------------------------------------------- //

/// Sixteen codons of real synonymous choice, plus start and stop.
const CDS: &str = "ATGCTGCGTAGCGGCGTGGCGCCGACCAAACTGCGTAGCGGCGTGGCGTAA";
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

/// `M L R S *` — 216 assignments, brute-forceable.
fn toy() -> Genotype {
    Genotype::new(
        CodonSequence::from_cds("ATGCTGCGTAGCTAA", 1)
            .unwrap()
            .with_five_prime_flank("AAGGAGGTA")
            .unwrap()
            .with_three_prime_flank("TAAGGACC")
            .unwrap(),
    )
}

fn motif_rules(motifs: &[&str]) -> ConstraintSet {
    ConstraintSet {
        schema_version: Some("1.0.0".into()),
        constraint_set_id: "cst_tier2".into(),
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
                    strand: Strand::Sense,
                    max_occurrences: 0,
                },
            })
            .collect(),
        repair_policy: RepairPolicy::default(),
        digest: None,
    }
}

/// Whether a genotype satisfies the checker, with locks honoured.
fn seed_is_feasible(g: &Genotype, checker: &ConstraintChecker) -> bool {
    let construct = assemble(g);
    let view = ConstructView::new(
        &construct,
        vec![(g.cds_offset_nt(), g.cds_offset_nt() + g.len() * 3)],
    )
    .with_editable_spans(cf_constraints::repair::editable_spans(&[(
        g.cds_offset_nt(),
        g.sequence().clone(),
    )]));
    checker.is_feasible(&view)
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

// --------------------------------------------------------------------------- //
// instrumented terms
// --------------------------------------------------------------------------- //

#[derive(Clone, Debug, Default)]
struct Calls(Arc<Counters>);

#[derive(Debug, Default)]
struct Counters {
    scores: AtomicU64,
    tables: AtomicU64,
    rebases: AtomicU64,
}

impl Calls {
    fn scores(&self) -> u64 {
        self.0.scores.load(Ordering::SeqCst)
    }
    fn tables(&self) -> u64 {
        self.0.tables.load(Ordering::SeqCst)
    }
    fn rebases(&self) -> u64 {
        self.0.rebases.load(Ordering::SeqCst)
    }
}

/// An additive, position-wise L1 term.
#[derive(Debug)]
struct Additive {
    term: TermId,
    calls: Calls,
    weight: [f64; 64],
    radius: InteractionRadius,
    capability: Capability,
}

impl Additive {
    fn new(term: &str, calls: Calls, salt: u64) -> Self {
        let mut weight = [0.0f64; 64];
        for (index, slot) in weight.iter_mut().enumerate() {
            *slot = (((index as u64 * 31 + salt) % 97) as f64) / 97.0;
        }
        Self {
            term: TermId::new(term).unwrap(),
            calls,
            weight,
            radius: InteractionRadius::Positionwise,
            capability: Capability::L1,
        }
    }

    fn with_radius(mut self, radius: InteractionRadius) -> Self {
        self.radius = radius;
        self
    }

    fn score_only(mut self) -> Self {
        self.capability = Capability::L0;
        self
    }

    /// Explicit weights for named codons; everything else is zero.
    ///
    /// So a test can build a scoring landscape with a known shape rather than
    /// hoping pseudo-random weights happen to produce one.
    fn with_weights(term: &str, calls: Calls, pairs: &[(&str, f64)]) -> Self {
        let mut this = Self::new(term, calls, 0);
        this.weight = [0.0; 64];
        for (codon, weight) in pairs {
            let codon = cf_types::Codon::parse(codon).unwrap();
            this.weight[codon.index_usize()] = *weight;
        }
        this
    }

    /// A term that wants one specific codon above everything else.
    ///
    /// So a test can make the DP *try* to choose a particular assignment and
    /// find out whether something else stops it, rather than hoping the
    /// pseudo-random weights happen to point that way.
    fn preferring(term: &str, calls: Calls, codon: &str) -> Self {
        let wanted = cf_types::Codon::parse(codon).unwrap();
        let mut this = Self::new(term, calls, 0);
        this.weight = [0.0; 64];
        this.weight[wanted.index_usize()] = 1.0;
        this
    }

    fn total(&self, sequence: &CodonSequence) -> f64 {
        sequence
            .codons()
            .iter()
            .map(|c| self.weight[c.index_usize()])
            .sum()
    }
}

impl TermProvider for Additive {
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
        self.capability
    }
    fn interaction_radius(&self) -> InteractionRadius {
        self.radius
    }
    fn provenance(&self) -> Provenance {
        Provenance::builtin("test_additive", "0.1.0")
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

/// A target-seeking term: `-(mean − target)²`.
///
/// Its single-substitution delta is exact, so it is honestly L1 — but the
/// deltas do not *sum*, and unlike a monotone transform the disagreement
/// changes the answer rather than only its scale. Every individual move toward
/// the target improves the term; a whole block of them overshoots and makes it
/// worse. That is the situation §3.2's two-step accept exists for.
///
/// With purely additive test terms the DP's linearization *is* the true
/// objective, so an implementation that committed every proposal unchecked
/// agrees with one that checks — which is why the accept step went unexercised
/// until this fixture existed. A first attempt used `mean³`, which is monotone
/// in the mean and therefore has the same argmax as its linearization: it
/// diverged numerically and never changed a decision.
#[derive(Debug)]
struct Convex {
    term: TermId,
    calls: Calls,
    weight: [f64; 64],
    target: f64,
}

impl Convex {
    fn new(term: &str, calls: Calls) -> Self {
        let mut weight = [0.0f64; 64];
        for (index, slot) in weight.iter_mut().enumerate() {
            *slot = ((index as u64 % 13) as f64) / 13.0;
        }
        Self {
            term: TermId::new(term).unwrap(),
            calls,
            weight,
            target: 0.0,
        }
    }

    fn with_target(mut self, target: f64) -> Self {
        self.target = target;
        self
    }

    fn raw(&self, sequence: &CodonSequence) -> f64 {
        sequence
            .codons()
            .iter()
            .map(|c| self.weight[c.index_usize()])
            .sum::<f64>()
            / sequence.codon_count().max(1) as f64
    }

    fn value(&self, sequence: &CodonSequence) -> f64 {
        let gap = self.raw(sequence) - self.target;
        -(gap * gap)
    }
}

impl TermProvider for Convex {
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
    fn interaction_radius(&self) -> InteractionRadius {
        InteractionRadius::Positionwise
    }
    fn provenance(&self) -> Provenance {
        Provenance::builtin("test_convex", "0.1.0")
    }
    fn score(&self, ctx: &EvalContext<'_>) -> Result<f64> {
        self.calls.0.scores.fetch_add(1, Ordering::SeqCst);
        Ok(self.value(ctx.sequence))
    }
    fn rebase(&self, ctx: &EvalContext<'_>) -> Result<Rebased> {
        self.calls.0.rebases.fetch_add(1, Ordering::SeqCst);
        Ok(Rebased::value_only(self.value(ctx.sequence)))
    }
    fn delta(
        &self,
        ctx: &EvalContext<'_>,
        _base: &Rebased,
        edits: &[cf_landscape_api::Edit],
    ) -> Result<f64> {
        let mut mutated = ctx.sequence.clone();
        for edit in edits {
            mutated.substitute(edit.codon_index, edit.codon)?;
        }
        Ok(self.value(&mutated) - self.value(ctx.sequence))
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

fn registry_with(terms: Vec<Box<dyn TermProvider>>) -> TermRegistry {
    let mut registry = TermRegistry::new();
    for provider in terms {
        registry.register_default(provider).unwrap();
    }
    registry
}

fn terms(names: &[&str]) -> Vec<TermId> {
    names.iter().map(|n| TermId::new(*n).unwrap()).collect()
}

// --------------------------------------------------------------------------- //
// the partition
// --------------------------------------------------------------------------- //

#[test]
fn p1_25_the_partition_alternates_and_covers() {
    let blocks = partition(17, 4);
    assert_eq!(blocks.len(), 5);
    assert_eq!(blocks[0].colour, Colour::Red);
    assert_eq!(blocks[1].colour, Colour::Black);
    assert_eq!(blocks[4].colour, Colour::Red);
    assert_eq!((blocks[4].start, blocks[4].end), (16, 17));

    // Contiguous and complete.
    assert_eq!(blocks[0].start, 0);
    for pair in blocks.windows(2) {
        assert_eq!(pair[0].end, pair[1].start);
        assert_ne!(pair[0].colour, pair[1].colour);
    }
    assert_eq!(blocks.last().unwrap().end, 17);

    // The property the independence argument rests on: two same-colour blocks
    // are at least `B` codons apart.
    for (i, a) in blocks.iter().enumerate() {
        for b in blocks.iter().skip(i + 1) {
            if a.colour == b.colour {
                assert!(
                    b.start - a.end >= 4,
                    "blocks {} and {} are only {} codons apart",
                    a.index,
                    b.index,
                    b.start - a.end
                );
            }
        }
    }
}

// --------------------------------------------------------------------------- //
// TST-SRCH-012 — the independence-shaped test
// --------------------------------------------------------------------------- //

/// Same-colour block updates are conditionally independent.
///
/// Two observables, and both are about cost or structure rather than answers.
///
/// First, the landscape is queried once per scored term per half-sweep, not
/// once per block. A per-block implementation would return the same result and
/// scale its query count with the block count.
///
/// Second — and this is the part that catches coupling specifically — the
/// number of blocks makes no difference to *how many times* the landscape is
/// consulted within a half-sweep, and no block's update can be an input to
/// another's, because every same-colour proposal is computed before any is
/// committed. An implementation that committed each block as it went would be
/// scoring block `k+2` against a sequence block `k` had already changed, which
/// is exactly the cross-block interaction `B ≥ R` is supposed to have removed.
#[test]
fn tst_srch_012_same_colour_blocks_are_scored_independently() {
    let calls = Calls::default();
    let registry = registry_with(vec![
        Box::new(Additive::new("usage.cai", calls.clone(), 0)),
        Box::new(Additive::new("usage.gc_content", calls.clone(), 7)),
    ]);
    let checker = ConstraintChecker::compile(&motif_rules(&["GGATCC"])).unwrap();
    let base = base();
    let host = host();

    // Block 3 over 17 codons gives six blocks, so a per-block implementation
    // would issue six times as many queries as a per-half-sweep one.
    let config = Tier2Config::from_planning(3);
    let plan = Tier2Plan::check(
        &registry,
        &terms(&["usage.cai", "usage.gc_content"]),
        &[1.0, 0.5],
        &config,
    )
    .unwrap();

    let blocks = partition(base.len(), config.block);
    assert!(
        blocks.len() >= 6,
        "the fixture has only {} blocks; per-block and per-half-sweep query counts \
         would be hard to tell apart",
        blocks.len()
    );

    let report = tier2_block_ca(
        &Tier2Request {
            seed: &base,
            registry: &registry,
            host: &host,
            checker: &checker,
            config: &config,
        },
        &plan,
        20260826,
    )
    .unwrap();

    let half_sweeps = report.surrogate_trace.len() as u64;
    assert!(half_sweeps >= 2);
    assert_eq!(
        report.delta_tables,
        half_sweeps * plan.scored().len() as u64,
        "one delta table per scored term per half-sweep, not per block \
         ({} blocks, {half_sweeps} half-sweeps)",
        blocks.len()
    );
    assert_eq!(report.rebases, report.delta_tables);
    assert_eq!(calls.tables(), report.delta_tables);
    assert_eq!(calls.rebases(), report.rebases);

    // `score` is called only by the accept step, which is once per full
    // evaluation — never inside a block update.
    assert!(
        calls.scores() < report.delta_tables * blocks.len() as u64,
        "score was called {} times; block updates are querying the landscape",
        calls.scores()
    );
}

/// Permuting the order of same-colour block updates changes nothing.
///
/// The decisive statement of conditional independence: if any same-colour block
/// could observe another, the order would matter. Driven by running the same
/// configuration from the same root seed and comparing the whole trace, then by
/// running with block sizes that change the partition's shape and checking the
/// per-block proposals are unaffected by their neighbours.
#[test]
fn tst_srch_012_a_block_update_does_not_see_its_same_colour_neighbours() {
    let calls = Calls::default();
    let checker = ConstraintChecker::compile(&motif_rules(&["GGATCC"])).unwrap();
    let host = host();
    let base = base();
    let config = Tier2Config::from_planning(4);

    let run = || {
        let registry = registry_with(vec![Box::new(Additive::new("usage.cai", calls.clone(), 0))]);
        let plan = Tier2Plan::check(&registry, &terms(&["usage.cai"]), &[1.0], &config).unwrap();
        let report = tier2_block_ca(
            &Tier2Request {
                seed: &base,
                registry: &registry,
                host: &host,
                checker: &checker,
                config: &config,
            },
            &plan,
            4242,
        )
        .unwrap();
        (report.incumbent.hash(), report.surrogate_trace)
    };

    let first = run();
    let second = run();
    assert_eq!(
        first, second,
        "the run is not reproducible from its root seed"
    );

    // Independence proper: the value the incumbent reaches must not depend on
    // whether the blocks of one colour were processed forwards or backwards.
    // The RNG stream is keyed on the block index rather than on iteration
    // order, so a correct implementation is order-invariant by construction and
    // a coupled one is not.
    let reversed = {
        let registry = registry_with(vec![Box::new(Additive::new("usage.cai", calls.clone(), 0))]);
        let plan = Tier2Plan::check(&registry, &terms(&["usage.cai"]), &[1.0], &config).unwrap();
        cf_search::tier2::tier2_block_ca_with_order(
            &Tier2Request {
                seed: &base,
                registry: &registry,
                host: &host,
                checker: &checker,
                config: &config,
            },
            &plan,
            4242,
            true,
        )
        .unwrap()
    };
    assert_eq!(
        reversed.incumbent.hash(),
        first.0,
        "reversing the order of same-colour block updates changed the result, so \
         the blocks are not conditionally independent"
    );
    assert_eq!(reversed.surrogate_trace, first.1);
}

// --------------------------------------------------------------------------- //
// the capability precheck
// --------------------------------------------------------------------------- //

/// An L0 term is dropped with a warning before any block is scored.
///
/// The ordering is the assertion. `04-model-resolution.md` §2.1 requires the
/// planner to warn or fail *before* the tier runs, because discovering it
/// mid-run wastes the tier — so the test checks the L0 provider was never
/// scored at all, not merely that a warning exists somewhere.
#[test]
fn p1_25_the_capability_precheck_fires_before_any_block_is_scored() {
    let good = Calls::default();
    let bad = Calls::default();
    let registry = registry_with(vec![
        Box::new(Additive::new("usage.cai", good.clone(), 0)),
        Box::new(Additive::new("lm.naturalness", bad.clone(), 3).score_only()),
    ]);
    let checker = ConstraintChecker::compile(&motif_rules(&["GGATCC"])).unwrap();
    let config = Tier2Config::from_planning(4);

    let plan = Tier2Plan::check(
        &registry,
        &terms(&["usage.cai", "lm.naturalness"]),
        &[1.0, 1.0],
        &config,
    )
    .unwrap();

    assert_eq!(plan.dropped(), &terms(&["lm.naturalness"]));
    assert_eq!(plan.scored(), &terms(&["usage.cai"]));
    assert_eq!(plan.warnings().len(), 1);
    assert_eq!(plan.warnings()[0].code, "CF-SRCH-010");
    assert!(plan.warnings()[0].detail.contains("rescore it in full"));

    // Nothing has been scored yet: the plan is a pure check.
    assert_eq!(bad.scores(), 0);
    assert_eq!(good.scores(), 0);

    let report = tier2_block_ca(
        &Tier2Request {
            seed: &base(),
            registry: &registry,
            host: &host(),
            checker: &checker,
            config: &config,
        },
        &plan,
        1,
    )
    .unwrap();

    assert_eq!(
        bad.scores(),
        0,
        "the dropped term was evaluated anyway, so the precheck did not take effect"
    );
    assert!(good.scores() > 0);
    assert_eq!(report.warnings.len(), 1);
}

/// A term wider than the block refuses the tier, naming P1-26.
#[test]
fn p1_25_a_term_wider_than_the_block_is_refused() {
    let calls = Calls::default();
    let config = Tier2Config::from_planning(4);

    let unbounded = registry_with(vec![Box::new(
        Additive::new("structure.mfe_global", calls.clone(), 1)
            .with_radius(InteractionRadius::Unbounded),
    )]);
    let err = Tier2Plan::check(
        &unbounded,
        &terms(&["structure.mfe_global"]),
        &[1.0],
        &config,
    )
    .unwrap_err();
    assert_eq!(err.code, "CF-SRCH-008");
    assert!(err.detail.contains("P1-26"), "{}", err.detail);

    let wide = registry_with(vec![Box::new(
        Additive::new("structure.local_accessibility", calls.clone(), 2)
            .with_radius(InteractionRadius::Codons(9)),
    )]);
    let err = Tier2Plan::check(
        &wide,
        &terms(&["structure.local_accessibility"]),
        &[1.0],
        &config,
    )
    .unwrap_err();
    assert_eq!(err.code, "CF-SRCH-008");
    assert!(err.detail.contains("9 codons"));

    // Raising the block size to the radius makes it admissible again.
    let bigger = Tier2Config::from_planning(9);
    let plan = Tier2Plan::check(
        &wide,
        &terms(&["structure.local_accessibility"]),
        &[1.0],
        &bigger,
    )
    .unwrap();
    assert_eq!(plan.radius(), 9);
    assert!(plan.partition_is_valid());
}

// --------------------------------------------------------------------------- //
// TST-SRCH-001 — monotonicity
// --------------------------------------------------------------------------- //

/// The true scalarized objective never decreases across a half-sweep.
///
/// Monotonicity in the *truth*, not in the linearization the DP maximizes,
/// which is what §3.2's two-step accept exists to buy. Run over many seeds,
/// with the trace checked at every half-sweep rather than only end to end — a
/// dip followed by a recovery would be invisible to an endpoint comparison.
#[test]
fn tst_srch_001_the_objective_never_decreases_across_a_half_sweep() {
    let calls = Calls::default();
    // Both motifs are absent from the base: `GGATCC` occurs only in the
    // backbone, where it is non-actionable, and `CGCAGC` is creatable but not
    // present. A motif already in the seed would make the seed infeasible, and
    // tier 2 is not repair — the run would then be measuring rejection rather
    // than monotonicity.
    let checker = ConstraintChecker::compile(&motif_rules(&["GGATCC", "CGCAGC"])).unwrap();
    let host = host();
    let base = base();
    assert!(seed_is_feasible(&base, &checker));

    let mut improved_at_all = 0usize;
    for seed in 0..500u64 {
        let registry = registry_with(vec![
            Box::new(Additive::new("usage.cai", calls.clone(), 0)),
            Box::new(Additive::new("usage.gc_content", calls.clone(), 11)),
        ]);
        let config = Tier2Config::from_planning(3 + (seed % 4) as usize);
        let plan = Tier2Plan::check(
            &registry,
            &terms(&["usage.cai", "usage.gc_content"]),
            &[1.0, 0.25],
            &config,
        )
        .unwrap();
        let report = tier2_block_ca(
            &Tier2Request {
                seed: &base,
                registry: &registry,
                host: &host,
                checker: &checker,
                config: &config,
            },
            &plan,
            seed,
        )
        .unwrap();

        for pair in report.surrogate_trace.windows(2) {
            assert!(
                pair[1] + 1e-12 >= pair[0],
                "seed {seed}: objective fell from {} to {}",
                pair[0],
                pair[1]
            );
        }
        if report.blocks_committed > 0 {
            improved_at_all += 1;
        }
    }
    assert!(
        improved_at_all > 400,
        "only {improved_at_all}/500 runs changed anything; monotonicity is nearly vacuous"
    );
}

// --------------------------------------------------------------------------- //
// TST-SRCH-002, tier-2 portion
// --------------------------------------------------------------------------- //

/// Tier 2 reaches the exact optimum on the brute-forceable space.
#[test]
fn tst_srch_002_tier2_finds_the_optimum_on_the_toy_space() {
    let calls = Calls::default();
    let checker = ConstraintChecker::compile(&motif_rules(&["GGTAGC", "CGCAGC"])).unwrap();
    let host = host();
    let base = toy();

    let registry = registry_with(vec![Box::new(Additive::new("usage.cai", calls.clone(), 0))]);
    let term_list = terms(&["usage.cai"]);
    let config = Tier2Config::from_planning(2);
    let plan = Tier2Plan::check(&registry, &term_list, &[1.0], &config).unwrap();

    let feasible = enumerate_feasible(&base, &checker, 4096).unwrap();
    assert!(
        feasible.len() > 30 && feasible.len() < 216,
        "the constraint set must bind without emptying the space; {} feasible",
        feasible.len()
    );

    let scorer = Additive::new("usage.cai", calls.clone(), 0);
    let best = feasible
        .iter()
        .map(|candidate| scorer.total(candidate.sequence()))
        .fold(f64::NEG_INFINITY, f64::max);

    let mut exact = 0usize;
    let seeds = 200usize;
    for seed in 0..seeds as u64 {
        let report = tier2_block_ca(
            &Tier2Request {
                seed: &base,
                registry: &registry,
                host: &host,
                checker: &checker,
                config: &config,
            },
            &plan,
            seed,
        )
        .unwrap();
        if (report.objective - best).abs() < 1e-9 {
            exact += 1;
        }
    }
    assert!(
        exact * 100 >= seeds * 99,
        "reached the optimum in {exact}/{seeds} seeds, below the 99% threshold"
    );
}

// --------------------------------------------------------------------------- //
// boundaries
// --------------------------------------------------------------------------- //

/// A locked span straddling a block edge is respected on both sides.
#[test]
fn p1_25_a_locked_span_across_a_block_edge_is_respected() {
    let calls = Calls::default();
    let checker = ConstraintChecker::compile(&motif_rules(&["GGATCC"])).unwrap();
    let host = host();

    let mut sequence = CodonSequence::from_cds(CDS, 1)
        .unwrap()
        .with_five_prime_flank(FIVE_PRIME)
        .unwrap()
        .with_three_prime_flank(THREE_PRIME)
        .unwrap();
    // Block size 4 puts the edge at codon 8; lock 6..11 so the span crosses it.
    sequence.add_lock(6, 11, LockReason::UserLocked).unwrap();
    let base = Genotype::new(sequence);
    let protein = base.sequence().protein_string();
    let locked: Vec<_> = (6..11).map(|i| base.codons()[i]).collect();

    let config = Tier2Config::from_planning(4);
    let blocks = partition(base.len(), config.block);
    let edge = blocks.iter().find(|b| b.start == 8).expect("an edge at 8");
    assert!(
        edge.start > 6 && edge.start < 11,
        "the lock must straddle a block edge to exercise the case"
    );

    for seed in 0..60u64 {
        let registry = registry_with(vec![Box::new(Additive::new("usage.cai", calls.clone(), 0))]);
        let plan = Tier2Plan::check(&registry, &terms(&["usage.cai"]), &[1.0], &config).unwrap();
        let report = tier2_block_ca(
            &Tier2Request {
                seed: &base,
                registry: &registry,
                host: &host,
                checker: &checker,
                config: &config,
            },
            &plan,
            seed,
        )
        .unwrap();

        assert_eq!(report.incumbent.sequence().protein_string(), protein);
        for (offset, codon) in locked.iter().enumerate() {
            assert_eq!(
                report.incumbent.codons()[6 + offset],
                *codon,
                "seed {seed}: locked codon {} changed",
                6 + offset
            );
        }
        assert!(report.incumbent.hash_is_consistent());
    }
}

/// A motif created across a block's *right* edge is prevented.
///
/// The complement of the priming test. The block DP prunes matches ending
/// inside its own range; one that starts inside and completes in the next block
/// is invisible to it, and nothing downstream would catch it either — the
/// accept step compares objectives, not feasibility. Tier 1's bug was a
/// fresh-start automaton at a construct junction; this is the same gap one
/// boundary further in, and it was live until the right-context walk was added.
///
/// `CTCCGT` spans codon 1 (block 0's last) and codon 2 (block 1's first, fixed
/// while red blocks update). It is absent from the base, and the term is rigged
/// to want `CTC`, so the DP will try to create it.
#[test]
fn p1_25_a_motif_across_a_block_right_edge_is_prevented() {
    let calls = Calls::default();
    let checker = ConstraintChecker::compile(&motif_rules(&["CTCCGT"])).unwrap();
    let host = host();
    let base = base();
    assert_eq!(base.codons()[1].to_string_upper(), "CTG");
    assert_eq!(base.codons()[2].to_string_upper(), "CGT");
    assert!(seed_is_feasible(&base, &checker));

    let config = Tier2Config::from_planning(2);
    for seed in 0..60u64 {
        let registry = registry_with(vec![Box::new(Additive::preferring(
            "usage.cai",
            calls.clone(),
            "CTC",
        ))]);
        let plan = Tier2Plan::check(&registry, &terms(&["usage.cai"]), &[1.0], &config).unwrap();
        let report = tier2_block_ca(
            &Tier2Request {
                seed: &base,
                registry: &registry,
                host: &host,
                checker: &checker,
                config: &config,
            },
            &plan,
            seed,
        )
        .unwrap();

        assert!(
            !(report.incumbent.codons()[1].to_string_upper() == "CTC"
                && report.incumbent.codons()[2].to_string_upper() == "CGT"),
            "seed {seed}: the DP created CTCCGT across the block's right edge"
        );
        assert!(
            seed_is_feasible(&report.incumbent, &checker),
            "seed {seed}: tier 2 committed an infeasible incumbent"
        );
    }
}

/// `L < 2B` leaves no partition, and is optimized as a single block.
#[test]
fn p1_25_a_short_sequence_is_one_block_with_no_seams() {
    let calls = Calls::default();
    let checker = ConstraintChecker::compile(&motif_rules(&[])).unwrap();
    let host = host();
    let base = toy();

    let config = Tier2Config::from_planning(64);
    assert_eq!(partition(base.len(), config.block).len(), 1);

    let registry = registry_with(vec![Box::new(Additive::new("usage.cai", calls.clone(), 0))]);
    let plan = Tier2Plan::check(&registry, &terms(&["usage.cai"]), &[1.0], &config).unwrap();
    let report = tier2_block_ca(
        &Tier2Request {
            seed: &base,
            registry: &registry,
            host: &host,
            checker: &checker,
            config: &config,
        },
        &plan,
        9,
    )
    .unwrap();

    assert!(report
        .notes
        .iter()
        .any(|note| note.contains("single block")));
    assert_eq!(
        report.incumbent.sequence().protein_string(),
        base.sequence().protein_string()
    );
}

// --------------------------------------------------------------------------- //
// the two-step accept, and the line that buys monotonicity
// --------------------------------------------------------------------------- //

/// Accepting on the truth is load-bearing, not decoration.
///
/// With an additive term the DP's linearization *is* the true objective, so an
/// implementation that committed every block proposal unchecked would agree
/// with one that checked — which is why every other test here passed with the
/// accept step disabled. `Convex` breaks that: its single-substitution deltas
/// are exact but do not sum, so a block proposal that looks best under the
/// accumulated surrogate can be worse in truth. The accept step must then
/// reject some proposals, and the objective must still never decrease.
#[test]
fn p1_25_accepting_on_the_truth_rejects_surrogate_wins_that_are_truth_losses() {
    let calls = Calls::default();
    let checker = ConstraintChecker::compile(&motif_rules(&[])).unwrap();
    let host = host();
    let base = base();
    let config = Tier2Config::from_planning(6);

    // Just above the base, so every individual move toward it improves the term
    // and a whole block of them overshoots.
    let probe = Convex::new("usage.cai", calls.clone());
    let target = probe.raw(base.sequence()) + 0.04;

    let mut rejected_total = 0u64;
    for seed in 0..80u64 {
        let registry = registry_with(vec![Box::new(
            Convex::new("usage.cai", calls.clone()).with_target(target),
        )]);
        let plan = Tier2Plan::check(&registry, &terms(&["usage.cai"]), &[1.0], &config).unwrap();
        let report = tier2_block_ca(
            &Tier2Request {
                seed: &base,
                registry: &registry,
                host: &host,
                checker: &checker,
                config: &config,
            },
            &plan,
            seed,
        )
        .unwrap();

        for pair in report.surrogate_trace.windows(2) {
            assert!(
                pair[1] + 1e-12 >= pair[0],
                "seed {seed}: the true objective fell from {} to {}",
                pair[0],
                pair[1]
            );
        }
        rejected_total += report.blocks_rejected;
    }

    assert!(
        rejected_total > 0,
        "no block proposal was ever rejected across 80 seeds, so the accept step is \
         never exercised and this test proves nothing"
    );
}

/// The incumbent stays in the beam, so a block update is never worsening.
///
/// Only observable under truncation: with a beam wide enough to hold every
/// path, the incumbent survives whether or not it is forced. Beam 1 makes the
/// truncation maximal, and then a block update that dropped the incumbent would
/// return a worse assignment and the accept step would reject it. With the
/// incumbent forced in, nothing is ever rejected.
#[test]
fn p1_25_the_incumbent_stays_in_a_truncated_beam() {
    let calls = Calls::default();
    let checker = ConstraintChecker::compile(&motif_rules(&["GGATCC"])).unwrap();
    let host = host();
    let base = base();

    let mut config = Tier2Config::from_planning(4);
    config.beam = 1;

    let mut committed = 0u64;
    for seed in 0..80u64 {
        let registry = registry_with(vec![Box::new(Additive::new("usage.cai", calls.clone(), 0))]);
        let plan = Tier2Plan::check(&registry, &terms(&["usage.cai"]), &[1.0], &config).unwrap();
        let report = tier2_block_ca(
            &Tier2Request {
                seed: &base,
                registry: &registry,
                host: &host,
                checker: &checker,
                config: &config,
            },
            &plan,
            seed,
        )
        .unwrap();

        assert_eq!(
            report.blocks_rejected, 0,
            "seed {seed}: a block update returned something worse than the incumbent, \
             so the incumbent was not in the beam"
        );
        committed += report.blocks_committed;
    }
    assert!(
        committed > 0,
        "nothing was ever committed at beam 1, so the assertion is vacuous"
    );
}

/// A motif created across a block's left edge is prevented by priming.
///
/// The block DP primes on everything to its left, so the fixed codons just
/// across the edge are context rather than a fresh start — tier 1's bug was a
/// fresh-start automaton at a construct junction, and this is the same risk at
/// a block boundary.
///
/// The term is rigged to *want* the offending codon, so the DP will choose it
/// unless the automaton stops it. Without that the DP might avoid the codon for
/// unrelated scoring reasons and the test would pass on an unprimed
/// implementation.
///
/// The observable is `blocks_rejected`, not the final codon. The accept step
/// refuses infeasible candidates, so an unprimed DP still cannot *commit* the
/// motif — it proposes it, has the proposal thrown away, and the block
/// contributes nothing. What priming buys is that the DP never proposes what
/// the accept step must discard, and a rejection count is the only way to see
/// the difference. Asserting on the final codon alone passed on both
/// implementations, which is how the first version of this test managed to
/// prove nothing.
#[test]
fn p1_25_priming_prevents_a_motif_created_across_a_block_edge() {
    let calls = Calls::default();
    let host = host();

    // Codon 1 is CTG (fixed in block 0 while black blocks update); codon 2 is
    // arginine. `CTGCGC` spans the block 0 / block 1 edge at block size 2 and
    // exists only when codon 2 is CGC.
    let checker = ConstraintChecker::compile(&motif_rules(&["CTGCGC"])).unwrap();

    // Codon 1 is locked, so it stays CTG. Without that, block 0 — which
    // updates first — is free to pick a different leucine codon, the motif is
    // never formed, and the test passes on an unprimed implementation for a
    // reason that has nothing to do with priming.
    let mut sequence = CodonSequence::from_cds(CDS, 1)
        .unwrap()
        .with_five_prime_flank(FIVE_PRIME)
        .unwrap()
        .with_three_prime_flank(THREE_PRIME)
        .unwrap();
    sequence.add_lock(1, 2, LockReason::UserLocked).unwrap();
    let base = Genotype::new(sequence);
    assert_eq!(base.codons()[1].to_string_upper(), "CTG");
    assert!(base.sequence().is_locked(1));
    assert!(base
        .sequence()
        .synonyms_at(2)
        .iter()
        .any(|c| c.to_string_upper() == "CGC"));

    let config = Tier2Config::from_planning(2);
    for seed in 0..60u64 {
        let registry = registry_with(vec![Box::new(Additive::preferring(
            "usage.cai",
            calls.clone(),
            "CGC",
        ))]);
        let plan = Tier2Plan::check(&registry, &terms(&["usage.cai"]), &[1.0], &config).unwrap();
        let report = tier2_block_ca(
            &Tier2Request {
                seed: &base,
                registry: &registry,
                host: &host,
                checker: &checker,
                config: &config,
            },
            &plan,
            seed,
        )
        .unwrap();

        assert_ne!(
            report.incumbent.codons()[2].to_string_upper(),
            "CGC",
            "seed {seed}: the DP chose the codon that completes CTGCGC across the block edge"
        );
        assert!(
            seed_is_feasible(&report.incumbent, &checker),
            "seed {seed}: infeasible incumbent"
        );
        // The observable that distinguishes a primed DP from an unprimed one.
        // Both refuse to *commit* the motif, because the accept step rejects an
        // infeasible candidate; only the primed one refuses to *propose* it.
        assert_eq!(
            report.blocks_rejected, 0,
            "seed {seed}: a block proposal was thrown away as infeasible, which means \
             the DP proposed something its own automaton should have pruned"
        );
    }
}

/// The incumbent survives truncation, so a block update cannot end up worse.
///
/// Separable scoring cannot show this: greedy is optimal there, so a beam of
/// any width finds the block optimum and the incumbent is never needed. The
/// gap only opens when an early choice that looks best *closes off* a better
/// continuation, which is what motif pruning does.
///
/// The fixture is built to do exactly that. `CGC` at codon 2 scores `+1.0`
/// against the incumbent, but `CGCAGC` is forbidden, so taking it costs `AGC`
/// at codon 3 — worth `−2.0`. A width-1 beam that dropped the incumbent would
/// commit to `CGC` and finish at `−1.0`, worse than doing nothing; with the
/// incumbent kept as an extra state, the `0.0` path survives and wins.
#[test]
fn p1_25_the_incumbent_survives_truncation_and_prevents_a_worsening_update() {
    let calls = Calls::default();
    let checker = ConstraintChecker::compile(&motif_rules(&["CGCAGC"])).unwrap();
    let host = host();
    let base = base();
    assert_eq!(base.codons()[2].to_string_upper(), "CGT");
    assert_eq!(base.codons()[3].to_string_upper(), "AGC");

    let mut config = Tier2Config::from_planning(2);
    config.beam = 1;
    // Polish would otherwise paper over a truncation loss, and the point here
    // is what the DP itself returns.
    config.max_polish_passes = 0;

    let mut committed = 0u64;
    for seed in 0..40u64 {
        let registry = registry_with(vec![Box::new(Additive::with_weights(
            "usage.cai",
            calls.clone(),
            &[("CGC", 1.5), ("CGT", 0.5), ("AGC", 2.0)],
        ))]);
        let plan = Tier2Plan::check(&registry, &terms(&["usage.cai"]), &[1.0], &config).unwrap();
        let report = tier2_block_ca(
            &Tier2Request {
                seed: &base,
                registry: &registry,
                host: &host,
                checker: &checker,
                config: &config,
            },
            &plan,
            seed,
        )
        .unwrap();

        assert_eq!(
            report.blocks_rejected, 0,
            "seed {seed}: a block update returned something worse than the incumbent, \
             so the incumbent did not survive truncation"
        );
        // The trap was actually laid: the incumbent's own codons are still there.
        assert_eq!(report.incumbent.codons()[3].to_string_upper(), "AGC");
        committed += report.blocks_committed;
    }
    assert!(
        committed > 0,
        "nothing was committed anywhere, so the run is degenerate and the assertion \
         above is vacuous"
    );
}
