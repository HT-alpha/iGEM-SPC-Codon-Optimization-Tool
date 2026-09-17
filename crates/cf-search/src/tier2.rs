//! Tier 2 — red-black block-local coordinate ascent.
//!
//! Ticket P1-25. Specification:
//! [`02-block-local-coordinate-ascent.md`](../../../planning/03-algorithms/02-block-local-coordinate-ascent.md)
//! §3.
//!
//! ## The partition, and what actually makes it sound
//!
//! Codon positions are cut into contiguous blocks of `B` and two-coloured.
//! Holding every black block fixed, the red blocks are conditionally
//! independent and can be optimized simultaneously with nothing to reconcile —
//! `DEP-SRCH-001`, the departure that removes the architecture proposal's
//! overlap-reconciliation step.
//!
//! §3.1 states the condition as `B ≥ W`. That is the condition only if `W` is
//! the largest interaction radius among the scored terms, and §4 exists
//! precisely because a term can have a larger one. So [`Tier2Plan::check`]
//! refuses `B < R` against radii the *providers* declare, not against a
//! configured `W`. A term that declares [`InteractionRadius::Unbounded`], or a
//! radius above `B`, is `CF-SRCH-008`, naming `P1-26`: the overlapping-block
//! fallback is the answer for that case, it is not built, and a silent
//! assumption that blocks are independent would be worse than a refusal.
//! `DEP-SRCH-014`: that is the only `B ≥ R` check. `CF-SRCH-009` is
//! configuration malformation (a zero block, a term/weight count mismatch),
//! not a second copy of the same inequality. `tier2_block_ca` used to repeat
//! `B ≥ R` under `009`; the orchestrator mapped only `008`, so the copy leaked
//! as a raw `ValueError`.
//!
//! **The partition is over codon positions, not editable positions.** §3.1's
//! prose says "partition the editable positions" while its diagram shows codon
//! ranges `[0,B) [B,2B) …`, and the two differ on a sequence with locks. Codon
//! positions are used, because the independence argument is about *codon*
//! distance: two positions in different red blocks are then at least `B` codons
//! apart directly, with no need to reason about how many locked positions sit
//! between them. A block that happens to contain only locked codons is a no-op.
//!
//! ## Surrogate and truth
//!
//! The achievement scalarizing function is not additive over positions, so it
//! cannot drive a chain DP. §3.2's two-step accept applies: linearize at the
//! incumbent to get an additive surrogate the DP can use, then **accept on the
//! truth** — evaluate the real `ASF` of the block proposal and take it only if
//! it does not worsen. Because the incumbent is always in the beam, "do
//! nothing" is always available, so acceptance is always possible and the
//! sweep is monotone in the true scalarized objective rather than in the
//! linearization.
//!
//! Tier 1 proposes and does not certify. Tier 2 makes an accept/reject decision
//! per block, so the truth is evaluated at the point the decision is made, not
//! afterwards.

use std::collections::HashMap;

use cf_constraints::{ConstraintChecker, ConstructView};
use cf_landscape_api::{
    Budget, BudgetTier, Capability, EvalContext, HostContext, InteractionRadius, Rebased,
    TermRegistry,
};
use cf_types::{code, CfError, Codon, Direction, Result, TermId};
use rayon::prelude::*;

use crate::genotype::{Genotype, Move, MoveKind};
use crate::rng::SearchRng;

/// Block colour.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub enum Colour {
    Red,
    Black,
}

impl Colour {
    pub fn as_str(self) -> &'static str {
        match self {
            Colour::Red => "red",
            Colour::Black => "black",
        }
    }

    fn index(self) -> u64 {
        match self {
            Colour::Red => 0,
            Colour::Black => 1,
        }
    }
}

/// One block of the partition.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Block {
    pub index: usize,
    pub colour: Colour,
    /// Half-open codon range.
    pub start: usize,
    pub end: usize,
}

impl Block {
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.start >= self.end
    }
}

/// How tier 2 is configured.
#[derive(Clone, Debug)]
pub struct Tier2Config {
    /// Block size `B`, in codons.
    pub block: usize,
    /// Beam width for the block DP.
    pub beam: usize,
    pub max_sweeps: u32,
    /// 1-opt polish passes after each block DP, per §3.3.
    pub max_polish_passes: u32,
}

impl Tier2Config {
    /// `B` and `b` from §3.4's worked cost model.
    ///
    /// Not defaults in the `Default` sense — a caller states them — but the
    /// values the cost analysis in §7 was computed against, so a caller with no
    /// reason to differ has a stated one to start from.
    pub fn from_planning(block: usize) -> Self {
        Self {
            block,
            beam: 32,
            max_sweeps: 20,
            max_polish_passes: 3,
        }
    }
}

/// The checked plan for a tier 2 run.
///
/// Built before any block is scored. Every refusal and every warning below is
/// produced here, which is what "the planner drops the term with a warning, or
/// fails" in [`04-model-resolution.md`](../../../planning/04-api/04-model-resolution.md)
/// §2.1 requires — discovering an L0 term mid-run wastes the whole tier.
/// `B ≥ R` is decided here as `CF-SRCH-008` (`DEP-SRCH-014`); a plan that
/// exists already satisfies it.
#[derive(Debug)]
pub struct Tier2Plan {
    /// Terms tier 2 will score, in the caller's order.
    scored: Vec<TermId>,
    weights: Vec<f64>,
    /// Terms dropped because they cannot be scored incrementally.
    dropped: Vec<TermId>,
    warnings: Vec<CfError>,
    radius: usize,
    block: usize,
}

impl Tier2Plan {
    /// Check a term set against a block size, before the tier starts.
    pub fn check(
        registry: &TermRegistry,
        terms: &[TermId],
        weights: &[f64],
        config: &Tier2Config,
    ) -> Result<Self> {
        if terms.len() != weights.len() {
            return Err(CfError::new(
                code::SRCH_TIER2_CONFIG,
                format!(
                    "tier 2 was given {} terms and {} weights",
                    terms.len(),
                    weights.len()
                ),
            ));
        }
        if config.block == 0 {
            return Err(CfError::new(
                code::SRCH_TIER2_CONFIG,
                "block size must be at least one codon",
            ));
        }

        let mut scored = Vec::new();
        let mut kept_weights = Vec::new();
        let mut dropped = Vec::new();
        let mut warnings = Vec::new();
        let mut radius = 0usize;

        for (term, weight) in terms.iter().zip(weights) {
            let registration = registry.resolve(std::slice::from_ref(term))?[0];
            let provider = registration.provider.as_ref();

            if !provider.capability().offers_exact_delta() {
                // §2.1: drop it with a warning rather than paying a full
                // rescore per move. Dropping changes the objective, so it is
                // said out loud; the alternative reading of §2.1 is to fail,
                // which a caller can do by treating the warning as fatal.
                dropped.push(term.clone());
                warnings.push(CfError::new(
                    code::SRCH_TIER2_TERM_DROPPED,
                    format!(
                        "{term}'s provider is {}, so every tier 2 move would rescore it in \
                         full. The term is dropped from tier 2 and the tier optimizes without \
                         it; supply an L1 provider to keep it",
                        provider.capability()
                    ),
                ));
                continue;
            }

            match provider.interaction_radius() {
                InteractionRadius::Unbounded => {
                    return Err(CfError::new(
                        code::SRCH_TIER2_UNBOUNDED_TERM,
                        format!(
                            "{term} declares an unbounded interaction radius, so same-colour \
                             blocks are not conditionally independent under it and the \
                             red-black partition does not apply. The overlapping-block \
                             reconciliation of §4 is the answer and arrives with P1-26; drop \
                             the term from tier 2 until then"
                        ),
                    ));
                }
                other => {
                    let codons = other.codons().expect("bounded");
                    if codons > config.block {
                        return Err(CfError::new(
                            code::SRCH_TIER2_UNBOUNDED_TERM,
                            format!(
                                "{term} couples positions up to {codons} codons apart, wider \
                                 than the block size {}. Same-colour blocks would not be \
                                 independent; raise the block size to at least {codons} or \
                                 wait for P1-26",
                                config.block
                            ),
                        ));
                    }
                    radius = radius.max(codons);
                }
            }

            scored.push(term.clone());
            kept_weights.push(*weight);
        }

        Ok(Self {
            scored,
            weights: kept_weights,
            dropped,
            warnings,
            radius,
            block: config.block,
        })
    }

    pub fn scored(&self) -> &[TermId] {
        &self.scored
    }

    pub fn dropped(&self) -> &[TermId] {
        &self.dropped
    }

    /// Warnings a job must surface, as `P1-23`'s operator warnings are.
    pub fn warnings(&self) -> &[CfError] {
        &self.warnings
    }

    /// The largest declared interaction radius among the scored terms.
    pub fn radius(&self) -> usize {
        self.radius
    }

    /// The block size `check` admitted. The run partitions against this, not
    /// against a request that might have drifted.
    pub fn block(&self) -> usize {
        self.block
    }

    /// Whether the partition is valid: `B ≥ R`.
    ///
    /// An invariant of any plan that `check` returned, not a second gate.
    /// `tier2_block_ca` debug-asserts it; it does not refuse again.
    pub fn partition_is_valid(&self) -> bool {
        self.block >= self.radius
    }
}

/// Partition codon positions into two-coloured contiguous blocks.
///
/// Returned in index order. A caller optimizes one colour at a time.
pub fn partition(codon_count: usize, block: usize) -> Vec<Block> {
    let block = block.max(1);
    let mut blocks = Vec::new();
    let mut start = 0usize;
    while start < codon_count {
        let end = (start + block).min(codon_count);
        let index = blocks.len();
        blocks.push(Block {
            index,
            colour: if index % 2 == 0 {
                Colour::Red
            } else {
                Colour::Black
            },
            start,
            end,
        });
        start = end;
    }
    blocks
}

/// What a tier 2 run did.
#[derive(Clone, Debug)]
pub struct Tier2Report {
    pub incumbent: Genotype,
    /// True `ASF` of the incumbent — the objective the accept step used.
    pub objective: f64,
    pub sweeps: u32,
    pub blocks_committed: u64,
    pub blocks_rejected: u64,
    /// Delta tables requested. One per scored term per half-sweep, never per
    /// block: block updates read from a table, they do not query the landscape.
    pub delta_tables: u64,
    pub rebases: u64,
    /// Per-half-sweep surrogate values, for the monotonicity property.
    pub surrogate_trace: Vec<f64>,
    pub warnings: Vec<CfError>,
    pub notes: Vec<String>,
}

/// Everything one tier 2 run needs.
#[derive(Clone, Debug)]
pub struct Tier2Request<'a> {
    pub seed: &'a Genotype,
    pub registry: &'a TermRegistry,
    pub host: &'a HostContext,
    pub checker: &'a ConstraintChecker,
    pub config: &'a Tier2Config,
}

/// Per-position, per-synonym surrogate increments, indexed by codon.
type Increments = Vec<HashMap<u8, f64>>;

/// Run tier 2.
pub fn tier2_block_ca(
    request: &Tier2Request<'_>,
    plan: &Tier2Plan,
    root_seed: u64,
) -> Result<Tier2Report> {
    tier2_block_ca_with_order(request, plan, root_seed, false)
}

/// As [`tier2_block_ca`], with the same-colour blocks visited in reverse.
///
/// Exists for `TST-SRCH-012` and for nothing else. Conditional independence
/// says the visiting order of same-colour blocks cannot affect the result, and
/// the only way to assert that is to vary it. A caller has no reason to reverse
/// the order; a test has every reason. Each block's RNG stream is keyed on its
/// index rather than on iteration order, so a correct implementation is
/// order-invariant by construction and a coupled one is not.
pub fn tier2_block_ca_with_order(
    request: &Tier2Request<'_>,
    plan: &Tier2Plan,
    root_seed: u64,
    reverse_blocks: bool,
) -> Result<Tier2Report> {
    let Tier2Request {
        seed,
        registry,
        host,
        checker,
        config,
    } = *request;

    // `DEP-SRCH-014`: `B ≥ R` was already decided by `Tier2Plan::check` as
    // `CF-SRCH-008`. Repeating it here under `CF-SRCH-009` was the same
    // inequality twice, and only `008` was mapped onto `job_status.problem`.
    debug_assert!(
        plan.partition_is_valid(),
        "Tier2Plan::check is the sole B≥R gate; a plan that fails it must not reach the run"
    );
    debug_assert_eq!(
        plan.block(),
        config.block,
        "the run must use the block size check() already admitted"
    );

    let mut notes = Vec::new();
    let block = plan.block();
    let blocks = partition(seed.len(), block);
    if blocks.len() < 2 {
        // §4: `L < 2W` leaves no useful partition, so the whole sequence is one
        // block. No reconciliation is needed either — there are no seams.
        notes.push(format!(
            "sequence is {} codons against a block size of {}; optimized as a single block",
            seed.len(),
            block
        ));
    }

    let mut incumbent = seed.clone();
    let mut objective = evaluate(&incumbent, registry, plan, host)?;
    let mut best = (incumbent.clone(), objective);

    let mut blocks_committed = 0u64;
    let mut blocks_rejected = 0u64;
    let mut delta_tables = 0u64;
    let mut rebases = 0u64;
    let mut surrogate_trace = Vec::new();
    let mut sweeps = 0u32;

    for sweep in 0..config.max_sweeps {
        sweeps = sweep + 1;
        let mut changed = false;

        for colour in [Colour::Red, Colour::Black] {
            // One delta table per scored term, for the whole half-sweep. Every
            // block reads from it, so scoring during a half-sweep is read-only
            // and same-colour blocks cannot observe each other.
            let (increments, pinned, tables) = build_increments(&incumbent, registry, plan, host)?;
            rebases += pinned;
            delta_tables += tables;

            let mut same_colour: Vec<&Block> =
                blocks.iter().filter(|b| b.colour == colour).collect();
            if reverse_blocks {
                same_colour.reverse();
            }

            // Same-colour blocks are conditionally independent by construction
            // (§3.1, `B ≥ R`), so they are proposed in parallel. `rayon` is
            // what ADR-0001 named for exactly this sweep.
            //
            // Determinism does not depend on the scheduling. Each block's
            // stream is keyed on `(sweep, colour, block index)` rather than on
            // iteration order, so a block sees the same randomness however
            // many others ran first or on which thread; the results are
            // collected into an index-ordered vector rather than a completion
            // -ordered one; and the commit below is sorted by block index. Two
            // thread counts therefore produce the same bytes, which is what
            // `GATE-DET-001` requires and what `TST-SRCH-004` measures.
            let proposals: Vec<(Block, Vec<Codon>)> = same_colour
                .par_iter()
                .map(|block| {
                    let mut rng = SearchRng::stream(
                        root_seed,
                        stream_index(sweep, colour, block.index, blocks.len()),
                    );
                    update_block(&incumbent, block, &increments, checker, config, &mut rng)
                        .map(|proposal| (**block, proposal))
                })
                .collect::<Result<Vec<_>>>()?;
            let mut proposals = proposals;

            // Committed after every same-colour block has been proposed, so no
            // block's update can be an input to another's. Sorted by index so
            // the commit order is the partition's order whatever order the
            // proposals were computed in.
            proposals.sort_by_key(|(block, _)| block.index);
            for (block, assignment) in proposals {
                let candidate = with_block(&incumbent, &block, &assignment)?;
                if candidate.same_codons_as(&incumbent) {
                    continue;
                }
                // Never commit an infeasible incumbent. The block DP prunes
                // within its own reach and the right-context walk covers its
                // trailing edge, so this should never fire — which is exactly
                // why it is here rather than assumed: the accept step compares
                // objectives, and an infeasible candidate can score well.
                if !block_is_feasible(&incumbent, &block, &assignment, checker) {
                    blocks_rejected += 1;
                    continue;
                }
                // Accept on the truth (§3.2 step 2), not on the linearization
                // the DP maximized.
                let candidate_objective = evaluate(&candidate, registry, plan, host)?;
                if candidate_objective + 1e-12 >= objective {
                    incumbent = candidate;
                    objective = candidate_objective;
                    blocks_committed += 1;
                    changed = true;
                } else {
                    blocks_rejected += 1;
                }
            }

            surrogate_trace.push(objective);
        }

        if objective > best.1 {
            best = (incumbent.clone(), objective);
        }
        if !changed {
            break;
        }
    }

    Ok(Tier2Report {
        incumbent: best.0,
        objective: best.1,
        sweeps,
        blocks_committed,
        blocks_rejected,
        delta_tables,
        rebases,
        surrogate_trace,
        warnings: plan.warnings().to_vec(),
        notes,
    })
}

/// The true scalarized objective.
///
/// A weighted sum over the scored terms in minimization sense, negated so the
/// caller maximizes. Not the augmented Tchebycheff `ASF` of §3.2: that needs
/// the ideal and nadir points, which come from `cf-metrics` and the tier-3
/// population, neither of which exists yet. The two-step accept structure is
/// what matters here and is identical either way — the DP maximizes an additive
/// linearization, and the decision is taken on a freshly evaluated true value.
/// `P1-25`'s ticket records the substitution rather than implying `ASF` is
/// implemented.
fn evaluate(
    candidate: &Genotype,
    registry: &TermRegistry,
    plan: &Tier2Plan,
    host: &HostContext,
) -> Result<f64> {
    let ctx = EvalContext::new(
        candidate.sequence(),
        host,
        Budget::new(BudgetTier::Standard),
        0,
    );
    let mut total = 0.0;
    for (term, weight) in plan.scored.iter().zip(&plan.weights) {
        let registration = registry.resolve(std::slice::from_ref(term))?[0];
        let provider = registration.provider.as_ref();
        let sense = match provider.direction() {
            Direction::Maximize => 1.0,
            Direction::Minimize => -1.0,
        };
        total += weight * sense * provider.score(&ctx)?;
    }
    Ok(total)
}

/// One delta table per scored term, against the current incumbent.
fn build_increments(
    incumbent: &Genotype,
    registry: &TermRegistry,
    plan: &Tier2Plan,
    host: &HostContext,
) -> Result<(Increments, u64, u64)> {
    let mut increments: Increments = vec![HashMap::new(); incumbent.len()];
    let positions = incumbent.mutable_positions();
    let ctx = EvalContext::new(incumbent.sequence(), host, Budget::new(BudgetTier::Fast), 0);
    let mut rebases = 0u64;
    let mut tables = 0u64;

    for (term, weight) in plan.scored.iter().zip(&plan.weights) {
        let registration = registry.resolve(std::slice::from_ref(term))?[0];
        let provider = registration.provider.as_ref();
        let pinned: Rebased = provider.rebase(&ctx)?;
        rebases += 1;
        let table = provider.delta_table(&ctx, &pinned, &positions)?;
        tables += 1;

        let sense = match provider.direction() {
            Direction::Maximize => 1.0,
            Direction::Minimize => -1.0,
        };
        for entry in table {
            let slot = &mut increments[entry.codon_index];
            for (codon, delta) in entry.alternatives {
                *slot.entry(codon.index()).or_insert(0.0) += weight * sense * delta;
            }
            slot.entry(incumbent.codons()[entry.codon_index].index())
                .or_insert(0.0);
        }
    }
    Ok((increments, rebases, tables))
}

/// A distinct RNG stream per `(sweep, colour, block)`.
fn stream_index(sweep: u32, colour: Colour, block: usize, blocks: usize) -> u64 {
    let per_sweep = 2 * blocks.max(1) as u64;
    u64::from(sweep) * per_sweep + colour.index() * blocks.max(1) as u64 + block as u64
}

/// `UPDATE_BLOCK` (§3.3), plus the 1-opt polish.
///
/// The incumbent assignment is kept in the beam at every step, which is what
/// makes the block update non-worsening in the surrogate and therefore what
/// buys convergence.
fn update_block(
    incumbent: &Genotype,
    block: &Block,
    increments: &Increments,
    checker: &ConstraintChecker,
    config: &Tier2Config,
    rng: &mut SearchRng,
) -> Result<Vec<Codon>> {
    let automata = checker.motif_automata();
    let sequence = incumbent.sequence();

    // Prime on everything to the left of the block: the 5' flank and every
    // codon before `block.start`. Those are fixed for this half-sweep, so a
    // motif inside them is not something this block can fix.
    let mut state = [automata[0].start_state(), automata[1].start_state()];
    let mut distance = crate::tier1::Distance::none();
    state = crate::tier1::walk(&automata, state, sequence.five_prime_flank());
    for index in 0..block.start {
        state = crate::tier1::walk(&automata, state, &sequence.codons()[index].to_bytes());
        distance = distance.after_codon(false);
    }

    #[derive(Clone)]
    struct State {
        dfa: [usize; 2],
        assign: Vec<Codon>,
        score: f64,
        on_incumbent: bool,
    }

    let mut beam = vec![State {
        dfa: state,
        assign: Vec::with_capacity(block.len()),
        score: 0.0,
        on_incumbent: true,
    }];

    for (index, increment) in increments
        .iter()
        .enumerate()
        .take(block.end)
        .skip(block.start)
    {
        let current = sequence.codons()[index];
        let locked = sequence.is_locked(index);
        let options: Vec<Codon> = if locked {
            vec![current]
        } else {
            sequence.synonyms_at(index).to_vec()
        };

        let mut next: Vec<State> = Vec::with_capacity(beam.len() * options.len());
        for st in &beam {
            for codon in &options {
                let Some(dfa) = crate::tier1::step_editable(
                    &automata,
                    checker,
                    st.dfa,
                    &codon.to_bytes(),
                    !locked,
                    distance,
                ) else {
                    continue;
                };
                let mut assign = st.assign.clone();
                assign.push(*codon);
                next.push(State {
                    dfa,
                    assign,
                    score: st.score + increment.get(&codon.index()).copied().unwrap_or(0.0),
                    on_incumbent: st.on_incumbent && *codon == current,
                });
            }
        }
        distance = distance.after_codon(!locked);

        if next.is_empty() {
            // The block cannot be reassigned at all without entering a
            // forbidden motif. That is not a failure: the incumbent is a valid
            // answer, and returning it is what "do nothing is always available"
            // means.
            return Ok(sequence.codons()[block.start..block.end].to_vec());
        }

        // Keep the best `b`, and always keep the incumbent path — the one line
        // that buys monotonicity (§3.3).
        let jitter: Vec<u64> = (0..next.len()).map(|_| rng.next_u64()).collect();
        let mut order: Vec<usize> = (0..next.len()).collect();
        order.sort_by(|a, b| {
            next[*b]
                .score
                .total_cmp(&next[*a].score)
                .then_with(|| jitter[*a].cmp(&jitter[*b]))
        });
        // The incumbent path is kept *in addition* to the best `b`, not in
        // competition with them. §3.3 writes the union before the truncation,
        // which reads as the incumbent being eligible for truncation — and then
        // "block updates are made non-worsening by seeding the beam with the
        // incumbent" is not a guarantee, only a tendency. Keeping it as a
        // `b + 1`th state makes it one.
        //
        // Displacing the best state instead, which is the other way to read it,
        // is worse still: at `b = 1` the beam would then always *be* the
        // incumbent and the block could never improve.
        order.truncate(config.beam.max(1));
        if !order.iter().any(|i| next[*i].on_incumbent) {
            if let Some(position) = next.iter().position(|s| s.on_incumbent) {
                order.push(position);
            }
        }
        order.sort_unstable();
        let mut kept = Vec::with_capacity(order.len());
        for index in order.into_iter().rev() {
            kept.push(next.swap_remove(index));
        }
        beam = kept;
    }

    // Walk the fixed context to the right of the block before choosing.
    //
    // The DP prunes matches that end inside the block; it cannot see one that
    // *starts* inside and completes in the next block. Those codons are fixed
    // for this half-sweep, so walking them is cheap and exact — and without it
    // a block update can create a motif across its right edge, which nothing
    // downstream would catch: the accept step compares objectives, not
    // feasibility. Structurally the same gap as tier 1's 3' junction, one
    // boundary further in.
    let right_reach = automata[0]
        .max_pattern_len()
        .max(automata[1].max_pattern_len());
    let right_end = (block.end + right_reach.div_ceil(3)).min(sequence.codon_count());
    let survivors: Vec<State> = beam
        .iter()
        .filter(|st| {
            let mut dfa = st.dfa;
            let mut d = distance;
            for index in block.end..right_end {
                match crate::tier1::step_editable(
                    &automata,
                    checker,
                    dfa,
                    &sequence.codons()[index].to_bytes(),
                    false,
                    d,
                ) {
                    Some(next) => dfa = next,
                    None => return false,
                }
                d = d.after_codon(false);
            }
            if right_end == sequence.codon_count() {
                return crate::tier1::step_editable(
                    &automata,
                    checker,
                    dfa,
                    sequence.three_prime_flank(),
                    false,
                    d,
                )
                .is_some();
            }
            true
        })
        .cloned()
        .collect();

    // If every path dies against the right context — which happens when the
    // incumbent itself is already violating, so its own path is pruned too —
    // doing nothing is still available and is the right answer.
    let beam = if survivors.is_empty() {
        beam
    } else {
        survivors
    };

    let best = beam
        .into_iter()
        .max_by(|a, b| a.score.total_cmp(&b.score))
        .expect("the beam is non-empty");
    Ok(polish(
        incumbent,
        block,
        best.assign,
        increments,
        checker,
        config,
    ))
}

/// 1-opt polish (§3.3): try every synonym at every position, accept
/// improvements, repeat until nothing changes.
///
/// Removes residual beam-truncation loss. Each candidate is checked against the
/// constraint checker rather than the automaton, because a single substitution
/// is cheap to verify exactly and the polish is off the inner loop.
fn polish(
    incumbent: &Genotype,
    block: &Block,
    mut assignment: Vec<Codon>,
    increments: &Increments,
    checker: &ConstraintChecker,
    config: &Tier2Config,
) -> Vec<Codon> {
    let sequence = incumbent.sequence();
    for _ in 0..config.max_polish_passes {
        let mut improved = false;
        for offset in 0..assignment.len() {
            let index = block.start + offset;
            if sequence.is_locked(index) {
                continue;
            }
            let current = assignment[offset];
            let current_score = increments[index]
                .get(&current.index())
                .copied()
                .unwrap_or(0.0);
            for candidate in sequence.synonyms_at(index) {
                if *candidate == current {
                    continue;
                }
                let score = increments[index]
                    .get(&candidate.index())
                    .copied()
                    .unwrap_or(0.0);
                if score <= current_score {
                    continue;
                }
                assignment[offset] = *candidate;
                if block_is_feasible(incumbent, block, &assignment, checker) {
                    improved = true;
                    break;
                }
                assignment[offset] = current;
            }
        }
        if !improved {
            break;
        }
    }
    assignment
}

fn block_is_feasible(
    incumbent: &Genotype,
    block: &Block,
    assignment: &[Codon],
    checker: &ConstraintChecker,
) -> bool {
    let Ok(candidate) = with_block(incumbent, block, assignment) else {
        return false;
    };
    let construct = assemble(&candidate);
    let view = ConstructView::new(&construct, vec![cds_span(&candidate)])
        .with_editable_spans(editable_spans(&candidate));
    checker.is_feasible(&view)
}

/// Apply a block assignment, through `Genotype::apply`.
///
/// So protein invariance and lock respect are inherited from
/// `CodonSequence::substitute` rather than restated — a locked position whose
/// assignment differed would fail here, which is the defensive check §3.3's
/// pseudocode performs by hand.
fn with_block(incumbent: &Genotype, block: &Block, assignment: &[Codon]) -> Result<Genotype> {
    let edits = assignment
        .iter()
        .enumerate()
        .filter(|(offset, codon)| incumbent.codons()[block.start + offset] != **codon)
        .map(|(offset, codon)| cf_landscape_api::Edit::new(block.start + offset, *codon))
        .collect();
    incumbent.with_move(&Move::new(MoveKind::Explicit, edits))
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

fn cds_span(g: &Genotype) -> (usize, usize) {
    let start = g.cds_offset_nt();
    (start, start + g.len() * 3)
}

/// Editable nucleotide spans, so a violation inside a locked span is not
/// counted against a candidate (`DEP-CONS-004`).
fn editable_spans(g: &Genotype) -> Vec<(usize, usize)> {
    cf_constraints::repair::editable_spans(&[(g.cds_offset_nt(), g.sequence().clone())])
}

/// Capability of each requested term, for a caller planning a tier.
pub fn planned_radii(
    registry: &TermRegistry,
    terms: &[TermId],
) -> Vec<(TermId, Capability, InteractionRadius)> {
    terms
        .iter()
        .filter_map(|term| {
            registry.get(term).map(|r| {
                (
                    term.clone(),
                    r.provider.capability(),
                    r.provider.interaction_radius(),
                )
            })
        })
        .collect()
}
