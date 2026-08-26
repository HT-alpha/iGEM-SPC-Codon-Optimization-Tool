//! Tier 1 — constrained beam DP seeding.
//!
//! Ticket P1-24. Specification:
//! [`02-block-local-coordinate-ascent.md`](../../../planning/03-algorithms/02-block-local-coordinate-ascent.md)
//! §2.
//!
//! Gets from "nothing" to "feasible and decent" fast: a beam search over the
//! codon chain, carried on the product with the motif automaton, so a state
//! that would enter a violating state is never generated and the output is
//! feasible with respect to the regular constraints **by construction**.
//!
//! ## Where the scores come from, and why that is the whole design
//!
//! The pseudocode's inner loop evaluates `Δ ← Σ_m λ_m · delta_m(x | ctx)` for
//! every (beam state, synonym) pair — `O(L·b·S)` of them, about 120 000 on a
//! 666-codon protein at `b = 64`. If each were a landscape query, tier 1 would
//! issue 120 000 queries, every one of them correct, and would be roughly three
//! orders of magnitude slower than the budget allows. Nothing about the answer
//! would look wrong.
//!
//! So the DP takes **one delta table per term, once, against the base**, and
//! reads increments out of it. `TST-SRCH-011` asserts that on landscape call
//! counts rather than on values, because values cannot distinguish the two.
//!
//! ## What the accumulated score is, precisely
//!
//! It is an **additive surrogate**: the sum over positions of the per-position
//! delta each term reports for the chosen synonym. That is exact for a term
//! whose contribution decomposes by position — `usage.gc_content` is one — and
//! it is a surrogate for one that does not. `usage.cai` is the instructive
//! case: it is log-additive, so maximizing the surrogate maximizes CAI exactly
//! when CAI is the only term, but the surrogate is not CAI's own arithmetic
//! difference. §3.2's rule applies and is implemented in
//! [`Tier1Report::seeds`]: seeds are re-scored on the truth before anything
//! ranks them. Tier 1 proposes; it does not certify.
//!
//! ## Two scope limits, stated rather than approximated
//!
//! Constraints that are not regular languages — `kmer_repeat`, which counts
//! over the whole construct, and `gc_window` — are not enforced in the DP.
//! Encoding a global occurrence count in the DP state would multiply the state
//! space by the count bound for every monitored k-mer. §2 defers them to the
//! post-check and the repair operator, and so does this.
//!
//! A motif rule with `max_occurrences > 0` is likewise not pruned: a state that
//! has seen fewer occurrences than the allowance is legal, and distinguishing
//! them needs a counter in the state. Those patterns are left to the
//! post-check, which is why the automaton has to report *which* pattern
//! matched rather than only that one did.

use std::collections::HashMap;

use cf_constraints::{ConstraintChecker, ConstructView, MotifAutomaton};
use cf_landscape_api::{Budget, BudgetTier, EvalContext, HostContext, Rebased, TermRegistry};
use cf_types::{code, CfError, Codon, CodonSequence, Direction, Result, TermId};

use crate::genotype::{Genotype, Move, MoveKind};
use crate::rng::SearchRng;

/// How tier 1 is configured.
#[derive(Clone, Debug)]
pub struct Tier1Config {
    /// Beam width `b`.
    pub beam: usize,
    /// Maximum states sharing one trailing context, per
    /// `TOP_B_BY_SCORE_WITH_DIVERSITY`.
    pub context_cap: usize,
    /// Trailing context length used for the diversity cap, in codons.
    ///
    /// The interaction window `W`. One means context-free, which is what the
    /// Phase 1 terms are.
    pub context_window: usize,
    /// How many distinct seeds to return.
    pub seeds: usize,
}

impl Default for Tier1Config {
    fn default() -> Self {
        Self {
            // §2's worked cost model uses b = 64.
            beam: 64,
            // Not a measured constant: an explicit, conservative cap that keeps
            // no more than a quarter of the beam on one lineage. Stated here so
            // it is a decision rather than a literal buried in a loop.
            context_cap: 16,
            context_window: 1,
            seeds: 4,
        }
    }
}

impl Tier1Config {
    pub fn with_beam(mut self, beam: usize) -> Self {
        self.beam = beam.max(1);
        self
    }

    pub fn with_seeds(mut self, seeds: usize) -> Self {
        self.seeds = seeds.max(1);
        self
    }

    pub fn with_context_window(mut self, window: usize) -> Self {
        self.context_window = window.max(1);
        self
    }
}

/// What tier 1 produced, and what it cost.
#[derive(Clone, Debug)]
pub struct Tier1Report {
    /// Distinct feasible seeds, best surrogate first.
    pub seeds: Vec<Genotype>,
    /// Surrogate score of each seed, in the same order.
    pub surrogate: Vec<f64>,
    /// Landscape delta tables requested. One per scored term, always.
    pub delta_tables: usize,
    /// Providers pinned. One per scored term, always.
    pub rebases: usize,
    /// Beam states generated, before pruning and truncation.
    pub states_generated: u64,
    /// States rejected because they entered a violating automaton state.
    pub states_pruned: u64,
    /// Seeds dropped by the post-check because a non-regular constraint —
    /// a k-mer repeat or a GC window — was violated.
    ///
    /// Not a failure. §2 hands those to the repair operator on purpose; the
    /// count is reported so a caller can see when tier 1's structural guarantee
    /// is covering less of the constraint set than it looks like it is.
    pub post_check_rejected: usize,
    pub notes: Vec<String>,
}

/// A live beam state.
#[derive(Clone, Debug)]
struct BeamState {
    /// Automaton states, one per automaton, in `motif_automata` order.
    dfa: [usize; 2],
    codons: Vec<Codon>,
    score: f64,
}

/// Per-position, per-synonym surrogate increments, indexed by codon.
///
/// Built once from the delta tables. `None` means the synonym was not offered
/// for that position, which is different from an increment of zero.
type Increments = Vec<HashMap<u8, f64>>;

/// Everything one tier 1 run needs.
///
/// A struct rather than eight parameters, because the four `&[..]`-shaped
/// arguments in the middle were positionally interchangeable and a caller that
/// swapped `terms` and `weights` would get a compile error only by luck.
#[derive(Clone, Debug)]
pub struct Tier1Request<'a> {
    /// Supplies the protein, the flanks and the lock set. The DP re-chooses
    /// every unlocked codon.
    pub base: &'a Genotype,
    pub registry: &'a TermRegistry,
    pub terms: &'a [TermId],
    /// Scalarization coefficients `λ`, one per term, in the same order.
    pub weights: &'a [f64],
    pub host: &'a HostContext,
    pub checker: &'a ConstraintChecker,
    pub config: &'a Tier1Config,
}

/// Run tier 1.
pub fn tier1_seed(request: &Tier1Request<'_>, rng: &mut SearchRng) -> Result<Tier1Report> {
    let Tier1Request {
        base,
        registry,
        terms,
        weights,
        host,
        checker,
        config,
    } = *request;

    if terms.len() != weights.len() {
        return Err(CfError::new(
            code::SRCH_TIER1_CONFIG,
            format!(
                "tier 1 was given {} terms and {} weights",
                terms.len(),
                weights.len()
            ),
        ));
    }

    let mut notes = Vec::new();
    let (increments, rebases, delta_tables) =
        build_increments(base, registry, terms, weights, host, &mut notes)?;

    let automata = checker.motif_automata();
    let mut beam = vec![prime(base, &automata)];
    let mut states_generated = 0u64;
    let mut states_pruned = 0u64;
    // Distance back to the nearest byte the search can change. Undefined
    // through the 5' flank and through leading locked codons, so a motif lying
    // wholly in immutable sequence is never pruned.
    let mut distance = Distance::none();

    for (index, increment) in increments.iter().enumerate() {
        let current = base.codons()[index];
        let locked = base.sequence().is_locked(index);
        // A locked position has exactly one option, and it is the one already
        // there. Nothing else may be proposed — §2's `if i in locked` line, and
        // the reason the DP never has to consult the lock set again.
        let options: Vec<Codon> = if locked {
            vec![current]
        } else {
            base.sequence().synonyms_at(index).to_vec()
        };

        let mut next: Vec<BeamState> = Vec::with_capacity(beam.len() * options.len());
        for state in &beam {
            for codon in &options {
                states_generated += 1;
                let Some(dfa) = step_editable(
                    &automata,
                    checker,
                    state.dfa,
                    &codon.to_bytes(),
                    !locked,
                    distance,
                ) else {
                    states_pruned += 1;
                    continue;
                };
                let mut codons = state.codons.clone();
                codons.push(*codon);
                next.push(BeamState {
                    dfa,
                    codons,
                    score: state.score + increment.get(&codon.index()).copied().unwrap_or(0.0),
                });
            }
        }

        distance = distance.after_codon(!locked);

        if next.is_empty() {
            return Err(CfError::new(
                code::SRCH_TIER1_INFEASIBLE,
                format!(
                    "tier 1 beam emptied at codon {index}: every synonym enters a forbidden \
                     motif from every surviving prefix. This is a beam result, not a proof — \
                     run the feasibility analyser for one"
                ),
            ));
        }
        beam = truncate(next, config, rng);
    }

    // Walk the 3' flank so a motif straddling that junction is pruned here
    // rather than found by the post-check. DEP-SRCH-006.
    let tail = distance;
    beam.retain(|state| {
        step_editable(
            &automata,
            checker,
            state.dfa,
            base.sequence().three_prime_flank(),
            false,
            tail,
        )
        .is_some()
    });
    if beam.is_empty() {
        return Err(CfError::new(
            code::SRCH_TIER1_INFEASIBLE,
            "every tier 1 candidate creates a forbidden motif across the 3' junction",
        ));
    }

    beam.sort_by(|a, b| b.score.total_cmp(&a.score));

    let mut seeds = Vec::new();
    let mut surrogate = Vec::new();
    let mut post_check_rejected = 0usize;
    for state in &beam {
        if seeds.len() >= config.seeds {
            break;
        }
        let candidate = base.with_move(&Move::new(
            MoveKind::Explicit,
            base.codons()
                .iter()
                .zip(&state.codons)
                .enumerate()
                .filter(|(_, (a, b))| a != b)
                .map(|(i, (_, b))| cf_landscape_api::Edit::new(i, *b))
                .collect(),
        ))?;

        // The non-regular constraints the DP does not carry.
        let construct = assemble(&candidate);
        let view = ConstructView::new(&construct, vec![cds_span(&candidate)]);
        if !checker.is_feasible(&view) {
            post_check_rejected += 1;
            continue;
        }
        if seeds
            .iter()
            .any(|s: &Genotype| s.same_codons_as(&candidate))
        {
            continue;
        }
        seeds.push(candidate);
        surrogate.push(state.score);
    }

    if seeds.is_empty() {
        return Err(CfError::new(
            code::SRCH_TIER1_INFEASIBLE,
            format!(
                "all {post_check_rejected} tier 1 candidates failed the post-check on a \
                 constraint the DP cannot carry (k-mer repeat or GC window); these are \
                 repair's to fix, so widen the beam or run repair on a seed"
            ),
        ));
    }
    if post_check_rejected > 0 {
        notes.push(format!(
            "{post_check_rejected} candidate(s) rejected by the post-check: tier 1 enforces \
             regular constraints by construction, not k-mer repeats or GC windows"
        ));
    }

    Ok(Tier1Report {
        seeds,
        surrogate,
        delta_tables,
        rebases,
        states_generated,
        states_pruned,
        post_check_rejected,
        notes,
    })
}

/// One delta table per term, once, against the base.
///
/// The count is the point. A DP that queried per state would be correct and
/// three orders of magnitude slower, and only these counters can tell.
fn build_increments(
    base: &Genotype,
    registry: &TermRegistry,
    terms: &[TermId],
    weights: &[f64],
    host: &HostContext,
    notes: &mut Vec<String>,
) -> Result<(Increments, usize, usize)> {
    let mut increments: Increments = vec![HashMap::new(); base.len()];
    let positions: Vec<usize> = base.mutable_positions();
    let ctx = EvalContext::new(base.sequence(), host, Budget::new(BudgetTier::Fast), 0);

    let mut rebases = 0usize;
    let mut delta_tables = 0usize;

    for (term, weight) in terms.iter().zip(weights) {
        let registration = registry.resolve(std::slice::from_ref(term))?[0];
        let provider = registration.provider.as_ref();
        if !provider.capability().offers_exact_delta() {
            // The tier planning check ADR-0009 §2.1 describes, arriving for
            // tier 1. Discovering this mid-run would waste the whole tier, and
            // silently dropping the term would optimize a different objective
            // than the caller asked for.
            return Err(CfError::new(
                code::LAND_CAPABILITY,
                format!(
                    "tier 1 needs a delta table for {term}, but its provider is {} — \
                     drop the term from tier 1 or supply an L1 provider",
                    provider.capability()
                ),
            ));
        }

        let pinned: Rebased = provider.rebase(&ctx)?;
        rebases += 1;
        let table = provider.delta_table(&ctx, &pinned, &positions)?;
        delta_tables += 1;

        // All internal optimization is maximization of the surrogate, so a
        // Minimize term is negated exactly once, here.
        let sense = match provider.direction() {
            Direction::Maximize => 1.0,
            Direction::Minimize => -1.0,
        };
        for entry in table {
            let slot = &mut increments[entry.codon_index];
            for (codon, delta) in entry.alternatives {
                *slot.entry(codon.index()).or_insert(0.0) += weight * sense * delta;
            }
            // The codon already present has, by definition, a delta of zero.
            slot.entry(base.codons()[entry.codon_index].index())
                .or_insert(0.0);
        }
    }

    if terms.is_empty() {
        notes.push(
            "tier 1 ran with no scored terms: every feasible assignment ties, so the seeds \
             are whatever the beam happened to keep"
                .into(),
        );
    }
    Ok((increments, rebases, delta_tables))
}

/// The automaton state after consuming the 5' flank.
///
/// `DEP-SRCH-006`. §2's pseudocode starts at `A.start` with an empty prefix,
/// which is right for a bare coding sequence and wrong for one inside a
/// construct: a motif straddling the 5' junction is then invisible to the DP.
/// Priming makes the junction structural instead of a post-check finding, which
/// matters because junction semantics are exactly what `EXIT-001`/`EXIT-002`
/// are about.
fn prime(base: &Genotype, automata: &[&MotifAutomaton; 2]) -> BeamState {
    let start = [automata[0].start_state(), automata[1].start_state()];
    BeamState {
        dfa: walk(automata, start, base.sequence().five_prime_flank()),
        codons: Vec::with_capacity(base.len()),
        score: 0.0,
    }
}

/// Advance both automata without pruning.
///
/// Used for the 5' flank. A motif lying entirely inside immutable sequence is a
/// property of the *input*, not of any candidate — `DEP-CONS-002`, the same
/// distinction `Violation::actionable` draws and the reason
/// `ConstraintChecker::is_feasible` ignores backbone-confined violations.
/// Pruning on one here would make tier 1 report every candidate as infeasible
/// for something no candidate can change, which is the failure mode the
/// feasibility analyser exists to report once instead.
fn walk(automata: &[&MotifAutomaton; 2], from: [usize; 2], bytes: &[u8]) -> [usize; 2] {
    let mut state = from;
    for byte in bytes {
        for (slot, automaton) in automata.iter().enumerate() {
            if automaton.is_empty() {
                continue;
            }
            state[slot] = automaton
                .transition(state[slot], *byte)
                .unwrap_or_else(|| automaton.start_state());
        }
    }
    state
}

/// Advance over bytes, pruning only on an *actionable* match.
///
/// `reach_nt` counts nucleotides since the last position the search could
/// change. A pattern of length `m` completing here began `m - 1` bytes earlier,
/// so it overlaps something editable exactly when `m >= reach_nt` after this
/// byte is counted.
///
/// The obvious rule — prune on any match ending inside the coding region — is
/// what tier 1 shipped with, and it is wrong for the same reason the 3' flank
/// walk was wrong. A locked span inside the CDS is immutable, so a motif lying
/// entirely within one is non-actionable (`DEP-CONS-002`): under the old rule
/// every beam state hit it, the beam emptied, and tier 1 reported
/// `CF-SRCH-006` for a violation no candidate could fix. Fourth time this
/// boundary has produced a bug, each time from a different direction.
fn step_editable(
    automata: &[&MotifAutomaton; 2],
    checker: &ConstraintChecker,
    from: [usize; 2],
    bytes: &[u8],
    editable: bool,
    distance: Distance,
) -> Option<[usize; 2]> {
    let mut state = from;
    let mut distance = distance;
    for byte in bytes {
        distance = distance.advance(editable);
        for (slot, automaton) in automata.iter().enumerate() {
            if automaton.is_empty() {
                continue;
            }
            let next = automaton
                .transition(state[slot], *byte)
                .unwrap_or_else(|| automaton.start_state());
            for pattern in automaton.outputs(next) {
                if !checker.pattern_is_zero_tolerance(automaton, *pattern) {
                    continue;
                }
                if distance.covers(automaton.patterns()[*pattern].len) {
                    return None;
                }
            }
            state[slot] = next;
        }
    }
    Some(state)
}

/// Nucleotides from the current byte back to the nearest editable one.
///
/// `None` means no editable byte has been seen yet, so nothing can reach one.
#[derive(Clone, Copy, Debug)]
struct Distance(Option<usize>);

impl Distance {
    /// Nothing editable seen yet.
    fn none() -> Self {
        Self(None)
    }

    /// Consume one byte.
    fn advance(self, editable: bool) -> Self {
        if editable {
            Self(Some(0))
        } else {
            Self(self.0.map(|d| d + 1))
        }
    }

    /// Consume a whole codon, without checking anything.
    fn after_codon(self, editable: bool) -> Self {
        if editable {
            Self(Some(0))
        } else {
            Self(self.0.map(|d| d + 3))
        }
    }

    /// Whether a pattern of `length` ending at the current byte covers an
    /// editable one. It spans the last `length` bytes, so it reaches back to
    /// distance `length - 1`.
    fn covers(self, length: usize) -> bool {
        self.0.is_some_and(|d| length > d)
    }
}

/// `TOP_B_BY_SCORE_WITH_DIVERSITY`.
///
/// Keeps the best `b` states but caps how many share a trailing context, which
/// prevents the beam collapsing onto one lineage — the standard beam-search
/// failure, and the one that makes a beam behave like a greedy walk while still
/// reporting a beam width.
fn truncate(
    mut states: Vec<BeamState>,
    config: &Tier1Config,
    rng: &mut SearchRng,
) -> Vec<BeamState> {
    // Ties are broken by a draw rather than by vector order, so the beam does
    // not inherit a bias from the order synonyms happen to be listed in.
    let jitter: Vec<u64> = (0..states.len()).map(|_| rng.next_u64()).collect();
    let mut order: Vec<usize> = (0..states.len()).collect();
    order.sort_by(|a, b| {
        states[*b]
            .score
            .total_cmp(&states[*a].score)
            .then_with(|| jitter[*a].cmp(&jitter[*b]))
    });

    let mut kept = Vec::with_capacity(config.beam);
    let mut per_context: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut overflow: Vec<usize> = Vec::new();

    for index in order {
        if kept.len() >= config.beam {
            break;
        }
        let context: Vec<u8> = states[index]
            .codons
            .iter()
            .rev()
            .take(config.context_window)
            .map(|c| c.index())
            .collect();
        let seen = per_context.entry(context).or_insert(0);
        if *seen >= config.context_cap {
            overflow.push(index);
            continue;
        }
        *seen += 1;
        kept.push(index);
    }
    // If the cap left the beam short — few distinct contexts exist near the
    // start of the sequence — refill from the best of what it excluded. A cap
    // that shrank the beam would trade one beam-search failure for another.
    for index in overflow {
        if kept.len() >= config.beam {
            break;
        }
        kept.push(index);
    }

    kept.sort_unstable();
    let mut out = Vec::with_capacity(kept.len());
    for index in kept.into_iter().rev() {
        out.push(states.swap_remove(index));
    }
    out
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

/// Exhaustive enumeration of the design space, for ground truth.
///
/// §2 notes that exact DP is affordable only when the state space is small, and
/// that the synthetic suite uses it as ground truth. This is the simplest exact
/// oracle: enumerate every synonymous assignment. It is `O(∏ |Syn(aᵢ)|)` and is
/// only for test instances — a caller that reaches for it on a real protein
/// will find out immediately.
pub fn enumerate_feasible(
    base: &Genotype,
    checker: &ConstraintChecker,
    cap: usize,
) -> Result<Vec<Genotype>> {
    let positions = base.mutable_positions();
    let space: usize = positions
        .iter()
        .map(|i| base.sequence().synonyms_at(*i).len())
        .try_fold(1usize, |acc, n| acc.checked_mul(n))
        .unwrap_or(usize::MAX);
    if space > cap {
        return Err(CfError::new(
            code::SRCH_TIER1_CONFIG,
            format!("design space is {space} assignments, above the enumeration cap of {cap}"),
        ));
    }

    let mut found = Vec::new();
    let mut assignment = vec![0usize; positions.len()];
    loop {
        let edits: Vec<cf_landscape_api::Edit> = positions
            .iter()
            .zip(&assignment)
            .map(|(index, choice)| {
                cf_landscape_api::Edit::new(*index, base.sequence().synonyms_at(*index)[*choice])
            })
            .collect();
        let candidate = base.with_move(&Move::new(MoveKind::Explicit, edits))?;
        let construct = assemble(&candidate);
        let view = ConstructView::new(&construct, vec![cds_span(&candidate)]);
        if checker.is_feasible(&view) {
            found.push(candidate);
        }

        let mut carry = positions.len();
        while carry > 0 {
            carry -= 1;
            assignment[carry] += 1;
            if assignment[carry] < base.sequence().synonyms_at(positions[carry]).len() {
                break;
            }
            assignment[carry] = 0;
            if carry == 0 {
                return Ok(found);
            }
        }
        if positions.is_empty() {
            return Ok(found);
        }
    }
}

/// The surrogate a tier 1 run would assign to an arbitrary assignment.
///
/// Exposed so a test can rank an exhaustive enumeration on the *same* scale the
/// DP optimizes, which is the only way to ask whether the DP found the optimum
/// rather than merely a good answer.
pub fn surrogate_of(
    candidate: &Genotype,
    base: &Genotype,
    registry: &TermRegistry,
    terms: &[TermId],
    weights: &[f64],
    host: &HostContext,
) -> Result<f64> {
    let mut notes = Vec::new();
    let (increments, _, _) = build_increments(base, registry, terms, weights, host, &mut notes)?;
    Ok(candidate
        .codons()
        .iter()
        .enumerate()
        .map(|(index, codon)| {
            increments[index]
                .get(&codon.index())
                .copied()
                .unwrap_or(0.0)
        })
        .sum())
}

/// A codon sequence's protein, for readability in test failures.
pub fn protein_of(sequence: &CodonSequence) -> String {
    sequence.protein_string()
}
