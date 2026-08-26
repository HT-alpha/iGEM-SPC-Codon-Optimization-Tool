//! `baseline_backtracking` — the `codon_6.0.py` solver, ported.
//!
//! Ticket P1-27b. The other half of the `P1-27` split: `P1-27a` ported the
//! prototype's *feasibility predicate* into `cf-constraints` and pinned it with
//! a differential against the loaded script; this ports the *solver*.
//!
//! ## What this is an oracle for, and what it is not
//!
//! The round's standard is that an oracle must be the actual prototype, not a
//! transcription that could share a bug with the thing it checks. That is right
//! about oracles and this is not one on its own. It is a *reference solver*: an
//! independently-structured search — depth-first with restarts, against tier
//! 1's breadth-first beam — whose witnesses each side can check with the
//! other's predicate.
//!
//! It earns that role only by being demonstrably faithful, so
//! `python/tests/test_p1_27b_baseline_fidelity.py` runs it against the real
//! `codon_6.0.py` on the same inputs. Until that passes, a disagreement between
//! this and tier 1 would be ambiguous between "tier 1 is wrong" and "the port
//! is wrong", which is the failure mode `P1-27a` avoided by loading the script.
//!
//! ## The two solvers are not comparable on score
//!
//! `optimize_cds_sequence` has no objective. It returns the first assignment
//! that survives its checks, ordered by a weighted shuffle over codon
//! frequency. Tier 1 maximizes a surrogate. So the comparable surface is
//! **feasibility**, and only over the constraints both enforce:
//!
//! | constraint | prototype | tier 1 DP |
//! | --- | --- | --- |
//! | forbidden motif / restriction site | yes, over the CDS | yes, by automaton product |
//! | k-mer repeat | yes, CDS + backbone counts | no — not a regular language (§2) |
//! | GC window | no | no |
//! | homopolymer | only via `auto_detect_forbidden_motifs` | as a motif, if compiled in |
//!
//! A differential that ignored the middle row would report disagreements that
//! are design, not drift.
//!
//! ## Faithful in the ways that matter, explicit where it cannot be
//!
//! Preserved: the stack machine and its ordering, so the search visits
//! candidates in the same order; the restart budget and the per-attempt
//! backtrack budget, so it gives up in the same places; the weighted shuffle
//! over codon frequency; and both windowed incremental checks, including the
//! `+2` slack the prototype uses.
//!
//! Not preserved: the RNG. The prototype uses Python's Mersenne Twister through
//! `random.uniform`, and reproducing that stream bit-for-bit would mean
//! vendoring MT19937 to imitate a script that never promised reproducibility.
//! So the port is deterministic under its own seed and the fidelity test
//! compares *properties* — solves-or-not, and validity of the answer — over
//! many seeds, rather than sequence equality.

use cf_constraints::legacy;
use cf_types::{code, CfError, Codon, CodonSequence, GeneticCode, Result};

use crate::genotype::{Genotype, Move, MoveKind};
use crate::rng::SearchRng;

/// The prototype's budgets.
///
/// Named rather than inlined because they are the prototype's numbers, not
/// ours: changing either changes what "the baseline solves this" means.
#[derive(Clone, Copy, Debug)]
pub struct BaselineBudget {
    /// `MAX_RESTARTS` in `codon_6.0.py`.
    pub max_restarts: u32,
    /// `MAX_BACKTRACKS_PER_ATTEMPT` in `codon_6.0.py`.
    pub max_backtracks_per_attempt: u32,
}

impl BaselineBudget {
    /// Read from the prototype at port time. Overridable so a test can force
    /// the give-up path without a long wait.
    pub fn from_prototype(max_restarts: u32, max_backtracks_per_attempt: u32) -> Self {
        Self {
            max_restarts,
            max_backtracks_per_attempt,
        }
    }
}

/// Why the baseline stopped.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BaselineOutcome {
    Solved,
    /// Every restart exhausted its backtrack budget.
    ///
    /// The prototype's `Exceeded max restarts` message. It is not a proof that
    /// no solution exists, and nothing may read it as one — that distinction is
    /// the same one `P1-17`'s `GiveUpReason` draws.
    BudgetExhausted,
}

/// What the baseline did.
#[derive(Clone, Debug)]
pub struct BaselineReport {
    pub outcome: BaselineOutcome,
    pub solution: Option<Genotype>,
    pub attempts: u32,
    pub backtracks: u64,
    /// Prefix extensions tried, including rejected ones.
    pub extensions: u64,
}

/// One stack frame: the prototype's `(aa_idx, current_seq, remaining_codons)`.
struct Frame {
    position: usize,
    prefix: Vec<Codon>,
    remaining: Vec<Codon>,
}

/// Run the ported backtracker.
///
/// `base` supplies the protein and the flanks; the solver re-chooses every
/// unlocked codon, exactly as the prototype re-chooses every codon of a CDS.
/// A locked position offers only the codon already there, which the prototype
/// has no concept of — stated as a divergence rather than smuggled in, since a
/// construct with locked codons is outside what the prototype ever handled.
pub fn baseline_backtracking(
    base: &Genotype,
    scenario: &legacy::LegacyScenario,
    weights: &crate::context::CodonWeights,
    budget: BaselineBudget,
    rng: &mut SearchRng,
) -> Result<BaselineReport> {
    let code = base.sequence().table();
    let protein: Vec<u8> = base.codons().iter().map(|c| code.amino_acid(*c)).collect();

    let mut attempts = 0u32;
    let mut backtracks = 0u64;
    let mut extensions = 0u64;

    for attempt in 1..=budget.max_restarts {
        attempts = attempt;
        let mut stack: Vec<Frame> = Vec::with_capacity(protein.len() + 1);
        stack.push(Frame {
            position: 0,
            prefix: Vec::new(),
            remaining: options(base, code, 0, weights, rng),
        });
        let mut attempt_backtracks = 0u32;

        while let Some(mut frame) = stack.pop() {
            if frame.position == protein.len() {
                let edits = frame
                    .prefix
                    .iter()
                    .enumerate()
                    .filter(|(index, codon)| base.codons()[*index] != **codon)
                    .map(|(index, codon)| cf_landscape_api::Edit::new(index, *codon))
                    .collect();
                let solution = base.with_move(&Move::new(MoveKind::Explicit, edits))?;
                return Ok(BaselineReport {
                    outcome: BaselineOutcome::Solved,
                    solution: Some(solution),
                    attempts,
                    backtracks,
                    extensions,
                });
            }

            let Some(next_codon) = frame.remaining.pop() else {
                // The prototype counts a frame with no options left as one
                // backtrack and drops it.
                backtracks += 1;
                attempt_backtracks += 1;
                if attempt_backtracks > budget.max_backtracks_per_attempt {
                    break;
                }
                continue;
            };

            let mut candidate = frame.prefix.clone();
            candidate.push(next_codon);
            let position = frame.position;
            stack.push(frame);
            extensions += 1;

            let bytes: Vec<u8> = candidate.iter().flat_map(|c| c.to_bytes()).collect();
            if !legacy::prefix_is_acceptable(&bytes, scenario) {
                continue;
            }

            let next_position = position + 1;
            let remaining = if next_position < protein.len() {
                options(base, code, next_position, weights, rng)
            } else {
                Vec::new()
            };
            stack.push(Frame {
                position: next_position,
                prefix: candidate,
                remaining,
            });
        }
    }

    Ok(BaselineReport {
        outcome: BaselineOutcome::BudgetExhausted,
        solution: None,
        attempts,
        backtracks,
        extensions,
    })
}

/// The prototype's `weighted_shuffle(REV_CODON_TABLE[aa])`.
///
/// Reversed at the end, as the prototype does, because the stack pops from the
/// back: the reversal is what makes the *highest*-weighted codon tried first.
fn options(
    base: &Genotype,
    code: &GeneticCode,
    position: usize,
    weights: &crate::context::CodonWeights,
    rng: &mut SearchRng,
) -> Vec<Codon> {
    if base.sequence().is_locked(position) {
        return vec![base.codons()[position]];
    }
    let mut pool: Vec<Codon> = code.synonyms(base.codons()[position]).to_vec();
    let mut order = Vec::with_capacity(pool.len());
    while !pool.is_empty() {
        let w = weights.weights_for(&pool);
        let choice = rng.weighted(&w).unwrap_or(0);
        order.push(pool.remove(choice));
    }
    order.reverse();
    order
}

/// The protein a solution must still encode.
pub fn protein_of(sequence: &CodonSequence) -> String {
    sequence.protein_string()
}

/// Reject a baseline solution that changed the protein.
///
/// The port routes every change through `Genotype::with_move`, so this cannot
/// fail — which is the point of asserting it: if it ever does, the port grew a
/// path around `CodonSequence::substitute`.
pub fn check_solution(base: &Genotype, solution: &Genotype) -> Result<()> {
    if base.sequence().protein_string() != solution.sequence().protein_string() {
        return Err(CfError::new(
            code::SEQ_PROTEIN_MISMATCH,
            "baseline_backtracking returned a solution encoding a different protein",
        ));
    }
    Ok(())
}
