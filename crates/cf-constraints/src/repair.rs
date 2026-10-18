//! The in-loop constraint repair operator.
//!
//! Ticket P1-17. Contract: `planning/03-algorithms/06-constraint-repair.md` §3.
//! Sequenced after `P1-19`; see §5.1 of the contract for why.
//!
//! # Three properties enforced by types rather than by discipline
//!
//! **Protein invariance falls out; it is not re-implemented.** The only
//! mutation path is [`RepairTarget::substitute`], which delegates to
//! `CodonSequence::substitute` — the precondition established in `P1-06`. A
//! non-synonymous replacement is rejected with `CF-SEQ-004` and a locked
//! position with `CF-SEQ-006`, so repair inherits both guarantees the same way
//! `P1-23`'s search operators will. There is deliberately no path that writes
//! to the construct bytes directly.
//!
//! **Only actionable violations are repairable.** [`ActionableViolation`] can
//! only be constructed from a [`Violation`] whose `actionable` flag is set, and
//! [`repair`] accepts nothing else. Handing repair a backbone-confined
//! violation is therefore not expressible, rather than merely discouraged
//! (`DEP-CONS-002`). [`partition`] splits a violation list into the two kinds
//! and is the intended way to route the rest to the feasibility analyser.
//!
//! **Work is bounded by a constant known before the run.** See
//! [`WorkAccounting`], which is compared against the declared bound by
//! `TST-CONS-003` rather than being documentation.
//!
//! # The termination condition, stated
//!
//! Repair terminates for two reasons at two different scales, and both are
//! needed — see §5.1 of the contract.
//!
//! *Within a call*: depth-1 search is **exhaustive** over the candidate window,
//! so "no depth-1 fix exists" is a proof about that window, not a timeout. The
//! escalation ladder is finite, and depth ≥ 2 is beam-limited. When the beam is
//! never truncated the whole search was exhaustive, and giving up is again a
//! proof; [`GiveUpReason`] records which of the two happened, so the caller can
//! tell "impossible here" from "ran out of budget".
//!
//! *Across calls*: `P1-19` runs once at submission and rejects constraint sets
//! no assignment can satisfy, so repair is never asked to solve an impossible
//! problem once per candidate per generation. That is the difference between
//! this and the prototype's `MAX_RESTARTS = 1000`, which also terminated — it
//! just terminated uselessly, a thousand times, with nothing to report.

use cf_types::error::{code, CfError, Result};
use cf_types::{Codon, CodonSequence};

use crate::detect::{ConstraintChecker, ConstructView, KmerMasking, Violation};
use crate::spec::{ConstraintClass, RepairPolicy};
use crate::tabu::TabuLedger;

// --------------------------------------------------------------------------- //
// objective awareness, without a module edge
// --------------------------------------------------------------------------- //

/// Scalarized objective cost of a candidate substitution.
///
/// Defined here as a trait rather than taking `cf-landscape-api`'s `DeltaTable`
/// directly. `cf-constraints` and `cf-landscape-api` are both rank 1 in
/// `modules.yaml`, so an edge between them is a rank inversion that
/// `tools/check_module_layering.py` rejects — correctly, because it is exactly
/// how a cycle gets introduced. Inverting the dependency keeps repair a pure
/// function of its arguments and keeps `cf-constraints` testable with a
/// hand-written implementation and no landscape present.
pub trait ObjectiveDelta {
    /// Change in the scalarized objective from placing `to` at this position.
    /// `None` when the provider has no opinion.
    fn delta(&self, region: usize, codon_index: usize, to: Codon) -> Option<f64>;
}

/// An `ObjectiveDelta` that has no opinion, for callers with no landscape.
#[derive(Debug, Default, Clone, Copy)]
pub struct NoObjective;

impl ObjectiveDelta for NoObjective {
    fn delta(&self, _region: usize, _codon_index: usize, _to: Codon) -> Option<f64> {
        None
    }
}

// --------------------------------------------------------------------------- //
// the actionable-violation boundary
// --------------------------------------------------------------------------- //

/// A violation repair is permitted to work on.
///
/// Constructible only from a violation whose `actionable` flag is set, which
/// makes `DEP-CONS-002`'s boundary a property of the type system rather than a
/// docstring.
#[derive(Clone, Copy, Debug)]
pub struct ActionableViolation<'a>(&'a Violation);

impl<'a> ActionableViolation<'a> {
    pub fn new(violation: &'a Violation) -> Option<Self> {
        violation.actionable.then_some(Self(violation))
    }

    pub fn get(&self) -> &'a Violation {
        self.0
    }
}

/// Split violations into the repairable and the analyser's.
///
/// The second half is what `P1-19`'s locked-span check reports; nothing here
/// tries to fix it.
pub fn partition(violations: &[Violation]) -> (Vec<ActionableViolation<'_>>, Vec<&Violation>) {
    let mut actionable = Vec::new();
    let mut for_analyser = Vec::new();
    for violation in violations {
        match ActionableViolation::new(violation) {
            Some(v) => actionable.push(v),
            None => for_analyser.push(violation),
        }
    }
    (actionable, for_analyser)
}

// --------------------------------------------------------------------------- //
// the mutable target
// --------------------------------------------------------------------------- //

/// A construct plus its editable coding regions, kept in sync.
///
/// The single mutation path in the repair operator.
#[derive(Clone, Debug)]
pub struct RepairTarget {
    construct: Vec<u8>,
    regions: Vec<(usize, CodonSequence)>,
}

impl RepairTarget {
    /// Build, verifying that each region's coding sequence really is at its
    /// declared offset. A silent mismatch here would let repair "fix" bytes
    /// that belong to something else.
    pub fn new(construct: Vec<u8>, regions: Vec<(usize, CodonSequence)>) -> Result<Self> {
        for (offset, sequence) in &regions {
            let end = offset + sequence.codon_count() * 3;
            if end > construct.len() {
                return Err(CfError::new(
                    code::SEQ_INDEX_OUT_OF_RANGE,
                    format!("region at {offset} extends past the construct"),
                ));
            }
            if construct[*offset..end] != *sequence.cds_string().as_bytes() {
                return Err(CfError::new(
                    code::SEQ_PROTEIN_MISMATCH,
                    format!("region at {offset} does not match the construct bytes there"),
                ));
            }
        }
        Ok(Self { construct, regions })
    }

    pub fn construct(&self) -> &[u8] {
        &self.construct
    }

    pub fn regions(&self) -> &[(usize, CodonSequence)] {
        &self.regions
    }

    pub fn cds_spans(&self) -> Vec<(usize, usize)> {
        self.regions
            .iter()
            .map(|(offset, sequence)| (*offset, offset + sequence.codon_count() * 3))
            .collect()
    }

    pub fn view(&self, masking: KmerMasking) -> ConstructView<'_> {
        ConstructView::new(&self.construct, self.cds_spans())
            .with_masking(masking)
            .with_editable_spans(editable_spans(&self.regions))
    }

    /// The only way to change anything.
    ///
    /// Delegates to `CodonSequence::substitute`, so synonymy and lock respect
    /// are enforced there, once, for every caller.
    pub fn substitute(&mut self, region: usize, codon_index: usize, to: Codon) -> Result<Codon> {
        let (offset, sequence) = self.regions.get_mut(region).ok_or_else(|| {
            CfError::new(
                code::SEQ_INDEX_OUT_OF_RANGE,
                format!("no editable region {region}"),
            )
        })?;
        let previous = sequence.substitute(codon_index, to)?;
        let at = *offset + codon_index * 3;
        self.construct[at..at + 3].copy_from_slice(&to.to_bytes());
        Ok(previous)
    }

    /// Proteins of every region, for invariance assertions.
    pub fn proteins(&self) -> Vec<String> {
        self.regions
            .iter()
            .map(|(_, sequence)| sequence.protein_string())
            .collect()
    }

    /// `(region, codon index)` pairs overlapping a nucleotide span, expanded by
    /// `radius` codons on each side, excluding locked positions.
    fn candidates_near(
        &self,
        start_nt: usize,
        end_nt: usize,
        radius: usize,
    ) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        for (region, (offset, sequence)) in self.regions.iter().enumerate() {
            let count = sequence.codon_count();
            let lo_nt = start_nt.saturating_sub(radius * 3);
            let hi_nt = end_nt.saturating_add(radius * 3);
            for index in 0..count {
                if sequence.is_locked(index) {
                    continue;
                }
                let codon_start = offset + index * 3;
                if codon_start + 3 > lo_nt && codon_start < hi_nt {
                    out.push((region, index));
                }
            }
        }
        out
    }
}

// --------------------------------------------------------------------------- //
// outcome
// --------------------------------------------------------------------------- //

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Edit {
    pub region: usize,
    pub codon_index: usize,
    pub from: Codon,
    pub to: Codon,
}

/// Why repair stopped working on a violation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GiveUpReason {
    /// The search was exhaustive over the final window and found nothing. A
    /// proof about that window, not a timeout.
    ProvablyNoFixInWindow,
    /// The beam was truncated, so absence of a fix was not established.
    BudgetExhausted,
}

#[derive(Clone, Debug, PartialEq)]
pub struct GiveUp {
    pub constraint_id: String,
    pub span_nt: (usize, usize),
    pub reason: GiveUpReason,
    pub final_window_codons: usize,
}

/// Measured work against the bound declared by the policy.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct WorkAccounting {
    pub candidates_evaluated: u64,
    pub declared_bound: u64,
}

impl WorkAccounting {
    pub fn within_bound(&self) -> bool {
        self.candidates_evaluated <= self.declared_bound
    }
}

#[derive(Clone, Debug)]
pub struct RepairOutcome {
    pub edits: Vec<Edit>,
    /// Actionable blocking violations still present when repair stopped.
    pub residual: Vec<Violation>,
    pub gave_up: Vec<GiveUp>,
    pub objective_cost: f64,
    pub work: WorkAccounting,
}

/// Worst-case candidate evaluations, computable before the run.
pub fn declared_bound(policy: &RepairPolicy, violations: usize, max_family: usize) -> u64 {
    let escalations = policy.max_escalations as u64 + 1;
    // Candidate positions grow with the window; the widest window after all
    // escalations bounds them.
    let widest = policy.max_window_codons as u64 * (1 << policy.max_escalations);
    let positions = 2 * widest + 8;
    let depth_1 = positions * max_family as u64;
    let deeper = policy.max_depth.saturating_sub(1) as u64
        * policy.beam_width as u64
        * positions
        * max_family as u64;
    violations as u64 * escalations * (depth_1 + deeper)
}

// --------------------------------------------------------------------------- //
// the operator
// --------------------------------------------------------------------------- //

/// Deterministic PRNG. Small, seeded, dependency-free: repair's tie-breaking
/// must be reproducible from the run's seed, which rules out `rand`'s thread
/// RNG and does not justify a dependency.
#[derive(Clone, Debug)]
pub struct RepairRng(u64);

impl RepairRng {
    /// Seed the stream.
    ///
    /// Not `seed | 1`: that keeps xorshift out of its absorbing zero state but
    /// maps every even seed onto its odd successor, so half of all consecutive
    /// seed pairs produce identical tie-breaking. The splitmix64 finalizer is a
    /// bijection, which keeps distinct seeds distinct. Same fix as
    /// `cf_search::SearchRng::new`, where the aliasing was found.
    pub fn new(seed: u64) -> Self {
        let mut z = seed ^ 0x9E37_79B9_7F4A_7C15;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        Self(if z == 0 { 0x2545_F491_4F6C_DD1D } else { z })
    }

    fn next_u64(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Editable nucleotide spans for a set of `(offset, sequence)` regions.
///
/// A locked codon is immutable, so it is excluded. Contiguous runs are merged
/// so the span list stays short on a sequence with few locks.
pub fn editable_spans(regions: &[(usize, CodonSequence)]) -> Vec<(usize, usize)> {
    let mut spans: Vec<(usize, usize)> = Vec::new();
    for (offset, sequence) in regions {
        for index in 0..sequence.codon_count() {
            if sequence.is_locked(index) {
                continue;
            }
            let start = offset + index * 3;
            match spans.last_mut() {
                Some(last) if last.1 == start => last.1 = start + 3,
                _ => spans.push((start, start + 3)),
            }
        }
    }
    spans
}

/// Repair the actionable violations of `target`, in place.
pub fn repair(
    target: &mut RepairTarget,
    checker: &ConstraintChecker,
    masking: KmerMasking,
    objective: &dyn ObjectiveDelta,
    tabu: &mut TabuLedger,
    rng: &mut RepairRng,
) -> Result<RepairOutcome> {
    let policy = checker.set().repair_policy.clone();
    let max_family = target
        .regions()
        .iter()
        .map(|(_, s)| s.table().max_family_len() as usize)
        .max()
        .unwrap_or(1);

    let initial = actionable_blocking(checker, target, masking);
    let mut outcome = RepairOutcome {
        edits: Vec::new(),
        residual: Vec::new(),
        gave_up: Vec::new(),
        objective_cost: 0.0,
        work: WorkAccounting {
            candidates_evaluated: 0,
            declared_bound: declared_bound(&policy, initial.len().max(1), max_family),
        },
    };

    // Order by span, then by constraint id, so repair is deterministic.
    let mut queue = initial.clone();
    queue.sort_by(|a, b| {
        (a.start_nt, a.end_nt, &a.constraint_id).cmp(&(b.start_nt, b.end_nt, &b.constraint_id))
    });

    for violation in queue {
        // An earlier repair may already have cleared it.
        if !still_present(checker, target, masking, &violation) {
            continue;
        }

        let mut radius = policy.max_window_codons as usize;
        let mut fixed = false;
        let mut beam_was_truncated = false;

        for _ in 0..=policy.max_escalations {
            let candidates = target.candidates_near(violation.start_nt, violation.end_nt, radius);
            match search_minimal_fix(
                target,
                checker,
                masking,
                &violation,
                &candidates,
                &policy,
                objective,
                tabu,
                rng,
                &mut outcome.work.candidates_evaluated,
                &mut beam_was_truncated,
            )? {
                Some(fix) => {
                    for edit in &fix.edits {
                        let previous = target.substitute(edit.region, edit.codon_index, edit.to)?;
                        debug_assert_eq!(previous, edit.from);
                        tabu.record(
                            edit.region,
                            edit.codon_index,
                            edit.from,
                            policy.tabu_tenure_generations,
                        );
                    }
                    outcome.objective_cost += fix.cost;
                    outcome.edits.extend(fix.edits);
                    fixed = true;
                    break;
                }
                None => radius *= 2,
            }
        }

        if !fixed {
            outcome.gave_up.push(GiveUp {
                constraint_id: violation.constraint_id.clone(),
                span_nt: (violation.start_nt, violation.end_nt),
                // Exhaustive unless the beam threw candidates away.
                reason: if beam_was_truncated {
                    GiveUpReason::BudgetExhausted
                } else {
                    GiveUpReason::ProvablyNoFixInWindow
                },
                final_window_codons: radius,
            });
        }
    }

    outcome.residual = actionable_blocking(checker, target, masking);
    Ok(outcome)
}

struct Fix {
    edits: Vec<Edit>,
    cost: f64,
}

#[allow(clippy::too_many_arguments)]
fn search_minimal_fix(
    target: &RepairTarget,
    checker: &ConstraintChecker,
    masking: KmerMasking,
    violation: &Violation,
    candidates: &[(usize, usize)],
    policy: &RepairPolicy,
    objective: &dyn ObjectiveDelta,
    tabu: &TabuLedger,
    rng: &mut RepairRng,
    work: &mut u64,
    beam_was_truncated: &mut bool,
) -> Result<Option<Fix>> {
    let before = actionable_blocking(checker, target, masking).len();

    // ---- depth 1: exhaustive over the window -----------------------------
    //
    // Exhaustive is the point. "No depth-1 fix exists" is then a proof about
    // this window rather than a timeout, which is what lets `GiveUpReason`
    // distinguish the two.
    let mut best: Option<(f64, f64, Vec<Edit>)> = None;
    for &(region, codon_index) in candidates {
        let (_, sequence) = &target.regions()[region];
        let current = sequence.codons()[codon_index];
        for &alternative in sequence.synonyms_at(codon_index) {
            if alternative == current {
                continue;
            }
            if tabu.is_tabu(region, codon_index, alternative) {
                continue;
            }
            *work += 1;
            let mut probe = target.clone();
            probe.substitute(region, codon_index, alternative)?;
            // One checker pass answers both questions. Two passes returned the
            // same multiset and doubled the cost of every probe.
            let (present, remaining) = probe_status(checker, &probe, masking, violation);
            if present || remaining >= before {
                // Still the original hit, or a fix that created another.
                continue;
            }
            let cost = objective
                .delta(region, codon_index, alternative)
                .unwrap_or(0.0);
            let key = (cost, rng.next_f64());
            let edits = vec![Edit {
                region,
                codon_index,
                from: current,
                to: alternative,
            }];
            if best
                .as_ref()
                .map(|(c, r, _)| (cost, key.1) < (*c, *r))
                .unwrap_or(true)
            {
                best = Some((cost, key.1, edits));
            }
        }
    }
    if let Some((cost, _, edits)) = best {
        return Ok(Some(Fix { edits, cost }));
    }

    // ---- depth >= 2: bounded beam over edit sets --------------------------
    let mut beam: Vec<(RepairTarget, Vec<Edit>, f64)> = vec![(target.clone(), Vec::new(), 0.0)];
    for _ in 2..=policy.max_depth {
        let mut next: Vec<(RepairTarget, Vec<Edit>, f64, u64, f64)> = Vec::new();
        for (state, edits, cost) in &beam {
            for &(region, codon_index) in candidates {
                if edits
                    .iter()
                    .any(|e| e.region == region && e.codon_index == codon_index)
                {
                    continue;
                }
                let (_, sequence) = &state.regions()[region];
                let current = sequence.codons()[codon_index];
                for &alternative in sequence.synonyms_at(codon_index) {
                    if alternative == current || tabu.is_tabu(region, codon_index, alternative) {
                        continue;
                    }
                    *work += 1;
                    let mut probe = state.clone();
                    probe.substitute(region, codon_index, alternative)?;
                    let (present, remaining) = probe_status(checker, &probe, masking, violation);
                    let mut extended = edits.clone();
                    extended.push(Edit {
                        region,
                        codon_index,
                        from: current,
                        to: alternative,
                    });
                    let extended_cost = cost
                        + objective
                            .delta(region, codon_index, alternative)
                            .unwrap_or(0.0);
                    if !present && remaining < before {
                        return Ok(Some(Fix {
                            edits: extended,
                            cost: extended_cost,
                        }));
                    }
                    next.push((
                        probe,
                        extended,
                        extended_cost,
                        remaining as u64,
                        rng.next_f64(),
                    ));
                }
            }
        }
        if next.is_empty() {
            break;
        }
        // Rank by (violations remaining, objective cost, random) and truncate.
        next.sort_by(|a, b| {
            (a.3, a.2, a.4)
                .partial_cmp(&(b.3, b.2, b.4))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        if next.len() > policy.beam_width as usize {
            *beam_was_truncated = true;
            next.truncate(policy.beam_width as usize);
        }
        beam = next
            .into_iter()
            .map(|(state, edits, cost, _, _)| (state, edits, cost))
            .collect();
    }
    Ok(None)
}

fn actionable_blocking(
    checker: &ConstraintChecker,
    target: &RepairTarget,
    masking: KmerMasking,
) -> Vec<Violation> {
    checker
        .check(&target.view(masking))
        .into_iter()
        .filter(|v| v.class == ConstraintClass::Blocking && v.actionable)
        .collect()
}

/// Whether `violation` is still reported, and how many blocking actionable
/// hits the probe has. Both come from one `check`.
fn probe_status(
    checker: &ConstraintChecker,
    target: &RepairTarget,
    masking: KmerMasking,
    violation: &Violation,
) -> (bool, usize) {
    let hits = actionable_blocking(checker, target, masking);
    let present = hits.iter().any(|hit| {
        hit.constraint_id == violation.constraint_id
            && hit.start_nt == violation.start_nt
            && hit.end_nt == violation.end_nt
    });
    (present, hits.len())
}

fn still_present(
    checker: &ConstraintChecker,
    target: &RepairTarget,
    masking: KmerMasking,
    violation: &Violation,
) -> bool {
    checker.check(&target.view(masking)).iter().any(|v| {
        v.constraint_id == violation.constraint_id
            && v.start_nt == violation.start_nt
            && v.end_nt == violation.end_nt
    })
}
