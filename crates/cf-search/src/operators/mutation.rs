//! The four mutation operators.
//!
//! Ticket P1-23. Specification:
//! [`03-nsga-loop-and-operators.md`](../../../../planning/03-algorithms/03-nsga-loop-and-operators.md)
//! §5.1.
//!
//! | operator | default weight | needs |
//! | --- | --- | --- |
//! | [`UniformSynonymous`] | 0.35 | host codon weights, optional |
//! | [`PottsGuided`] | 0.35 | a [`ContextModel`] — `P3-05` |
//! | [`BlockGibbs`] | 0.20 | a [`ContextModel`] — `P3-05` |
//! | [`ViolationDirected`] | 0.10 | actionable violations |
//!
//! None of them re-checks synonymy or lock membership. They sample from
//! `Genotype::mutable_positions` and `CodonSequence::synonyms_at`, and the move
//! is applied through `CodonSequence::substitute`, which is where both checks
//! already live. A second copy of either check here would be a second
//! definition, free to drift from the first.

use std::sync::Arc;

use cf_landscape_api::Edit;
use cf_types::{Codon, Result};

use crate::context::ContextModel;
use crate::genotype::{Move, MoveKind};
use crate::operators::{MutationContext, MutationOperator};
use crate::rng::SearchRng;

/// Sample one admissible replacement at `index`, honouring the tabu ledger.
///
/// Returns `None` when every alternative is suppressed. That is a real outcome
/// — a heavily tabu position has nothing to offer — and reporting it lets the
/// caller move on to another position instead of proposing a suppressed edit
/// and having repair undo it, which is the loop `P1-18` exists to break.
fn sample_replacement(
    ctx: &MutationContext<'_>,
    index: usize,
    rng: &mut SearchRng,
) -> Option<Codon> {
    let current = ctx.genotype.codons()[index];
    let mut candidates: Vec<Codon> = ctx
        .genotype
        .sequence()
        .synonyms_at(index)
        .iter()
        .copied()
        .filter(|c| *c != current)
        .collect();

    if let Some(tabu) = ctx.tabu {
        // The roll is drawn once per position, not once per candidate, so the
        // leak is a property of the position rather than of how many synonyms
        // it happens to have.
        let roll = rng.next_f64();
        candidates.retain(|c| !tabu.should_suppress(ctx.region, index, *c, roll));
    }
    if candidates.is_empty() {
        return None;
    }

    let weights = ctx.weights.weights_for(&candidates);
    rng.weighted(&weights)
        .or_else(|| Some(rng.below(candidates.len())))
        .map(|choice| candidates[choice])
}

/// Sample from a context model's conditional at `index`.
fn sample_conditional(
    model: &dyn ContextModel,
    ctx: &MutationContext<'_>,
    index: usize,
    temperature: f64,
    scratch: &mut Vec<(Codon, f64)>,
    rng: &mut SearchRng,
) -> Result<Option<Codon>> {
    model.conditional(ctx.genotype.sequence(), index, temperature, scratch)?;
    let current = ctx.genotype.codons()[index];
    scratch.retain(|(codon, _)| *codon != current);

    if let Some(tabu) = ctx.tabu {
        let roll = rng.next_f64();
        scratch.retain(|(codon, _)| !tabu.should_suppress(ctx.region, index, *codon, roll));
    }
    if scratch.is_empty() {
        return Ok(None);
    }

    let probabilities: Vec<f64> = scratch.iter().map(|(_, p)| *p).collect();
    Ok(rng.weighted(&probabilities).map(|choice| scratch[choice].0))
}

/// Each editable position mutates independently with `p_m = 1/L_edit`.
#[derive(Debug)]
pub struct UniformSynonymous;

impl MutationOperator for UniformSynonymous {
    fn name(&self) -> &'static str {
        "uniform_synonymous"
    }

    fn kind(&self) -> MoveKind {
        MoveKind::UniformSynonymous
    }

    fn propose(&self, ctx: &MutationContext<'_>, rng: &mut SearchRng) -> Result<Move> {
        let positions = ctx.genotype.mutable_positions();
        if positions.is_empty() {
            return Ok(Move::empty(self.kind()));
        }
        let rate = 1.0 / positions.len() as f64;

        let mut edits = Vec::new();
        for index in &positions {
            if rng.next_f64() >= rate {
                continue;
            }
            if let Some(codon) = sample_replacement(ctx, *index, rng) {
                edits.push(Edit::new(*index, codon));
            }
        }
        // An expected count of 1 still leaves a sizeable chance of zero, and a
        // mutation that does nothing wastes a generation slot. One forced draw
        // keeps the operator's contract — it always proposes a change when one
        // is available — without changing the per-position rate.
        if edits.is_empty() {
            let index = positions[rng.below(positions.len())];
            if let Some(codon) = sample_replacement(ctx, index, rng) {
                edits.push(Edit::new(index, codon));
            }
        }
        Ok(Move::new(self.kind(), edits))
    }
}

/// One position, replaced from the context model's conditional.
pub struct PottsGuided {
    model: Arc<dyn ContextModel>,
    temperature: f64,
}

impl std::fmt::Debug for PottsGuided {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PottsGuided")
            .field("model", &self.model.name())
            .field("temperature", &self.temperature)
            .finish()
    }
}

impl PottsGuided {
    /// `T = 1.5` per §5.1.
    ///
    /// The temperature is not a tuning knob left at a plausible value: §5.1
    /// states why it is above 1. Sampling the conditional at `T = 1` drags the
    /// population toward host-typical sequence and erodes the diversity the
    /// multi-objective search depends on, so the conditional is deliberately
    /// flattened.
    pub const DEFAULT_TEMPERATURE: f64 = 1.5;

    pub fn new(model: Arc<dyn ContextModel>) -> Self {
        Self {
            model,
            temperature: Self::DEFAULT_TEMPERATURE,
        }
    }

    pub fn with_temperature(mut self, temperature: f64) -> Self {
        self.temperature = temperature;
        self
    }
}

impl MutationOperator for PottsGuided {
    fn name(&self) -> &'static str {
        "potts_guided"
    }

    fn kind(&self) -> MoveKind {
        MoveKind::PottsGuided
    }

    fn propose(&self, ctx: &MutationContext<'_>, rng: &mut SearchRng) -> Result<Move> {
        let positions = ctx.genotype.mutable_positions();
        if positions.is_empty() {
            return Ok(Move::empty(self.kind()));
        }
        let index = positions[rng.below(positions.len())];
        let mut scratch = Vec::new();
        let chosen = sample_conditional(
            self.model.as_ref(),
            ctx,
            index,
            self.temperature,
            &mut scratch,
            rng,
        )?;
        Ok(match chosen {
            Some(codon) => Move::new(self.kind(), vec![Edit::new(index, codon)]),
            None => Move::empty(self.kind()),
        })
    }
}

/// One Gibbs sweep through a contiguous block.
pub struct BlockGibbs {
    model: Arc<dyn ContextModel>,
    temperature: f64,
    mean_length: f64,
}

impl std::fmt::Debug for BlockGibbs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BlockGibbs")
            .field("model", &self.model.name())
            .field("mean_length", &self.mean_length)
            .finish()
    }
}

impl BlockGibbs {
    /// `Geometric(1/12)` per §5.1.
    pub const DEFAULT_MEAN_LENGTH: f64 = 12.0;

    pub fn new(model: Arc<dyn ContextModel>) -> Self {
        Self {
            model,
            temperature: PottsGuided::DEFAULT_TEMPERATURE,
            mean_length: Self::DEFAULT_MEAN_LENGTH,
        }
    }
}

impl MutationOperator for BlockGibbs {
    fn name(&self) -> &'static str {
        "block_gibbs"
    }

    fn kind(&self) -> MoveKind {
        MoveKind::BlockGibbs
    }

    fn propose(&self, ctx: &MutationContext<'_>, rng: &mut SearchRng) -> Result<Move> {
        let positions = ctx.genotype.mutable_positions();
        if positions.is_empty() {
            return Ok(Move::empty(self.kind()));
        }
        let length = rng
            .geometric(1.0 / self.mean_length, positions.len())
            .min(positions.len());
        let start = rng.below(positions.len() - length + 1);

        // The sweep is over a working copy: a Gibbs step conditions on the
        // positions already resampled, which is what makes the move coherent
        // rather than a batch of independent draws.
        let mut working = ctx.genotype.clone();
        let mut scratch = Vec::new();
        let mut edits = Vec::new();

        for offset in 0..length {
            let index = positions[start + offset];
            let stepped = MutationContext {
                genotype: &working,
                ..*ctx
            };
            let chosen = sample_conditional(
                self.model.as_ref(),
                &stepped,
                index,
                self.temperature,
                &mut scratch,
                rng,
            )?;
            if let Some(codon) = chosen {
                let step = Move::new(self.kind(), vec![Edit::new(index, codon)]);
                working.apply(&step)?;
                edits.push(Edit::new(index, codon));
            }
        }
        Ok(Move::new(self.kind(), edits))
    }
}

/// Mutate positions inside a current violation span.
///
/// Takes `ActionableViolation`, not `Violation`. `DEP-CONS-002`: a blocking
/// violation confined to immutable sequence is a property of the input rather
/// than of the candidate, and an operator that tried to "escape" one would be
/// mutating positions that have nothing to do with it. The type makes the
/// distinction impossible to skip, the same way it does for `P1-17`'s repair.
#[derive(Debug)]
pub struct ViolationDirected;

impl MutationOperator for ViolationDirected {
    fn name(&self) -> &'static str {
        "violation_directed"
    }

    fn kind(&self) -> MoveKind {
        MoveKind::ViolationDirected
    }

    /// Only for individuals that actually have somewhere to escape from.
    fn applicable(&self, ctx: &MutationContext<'_>) -> bool {
        ctx.violations.iter().any(|v| {
            let violation = v.get();
            !ctx.genotype
                .codons_in_span(violation.start_nt, violation.end_nt)
                .is_empty()
        })
    }

    fn propose(&self, ctx: &MutationContext<'_>, rng: &mut SearchRng) -> Result<Move> {
        let mutable = ctx.genotype.mutable_positions();
        let mut targets: Vec<usize> = Vec::new();
        for violation in ctx.violations {
            let span = violation.get();
            for index in ctx.genotype.codons_in_span(span.start_nt, span.end_nt) {
                // Intersected with the mutable set rather than trusted: a span
                // can overlap a locked codon at its edge, and the whole point
                // of an operator is to propose moves that can be applied.
                if mutable.contains(&index) && !targets.contains(&index) {
                    targets.push(index);
                }
            }
        }
        if targets.is_empty() {
            return Ok(Move::empty(self.kind()));
        }

        let index = targets[rng.below(targets.len())];
        Ok(match sample_replacement(ctx, index, rng) {
            Some(codon) => Move::new(self.kind(), vec![Edit::new(index, codon)]),
            None => Move::empty(self.kind()),
        })
    }
}
