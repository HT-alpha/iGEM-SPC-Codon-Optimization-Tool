//! Crossover.
//!
//! Ticket P1-23. Specification:
//! [`03-nsga-loop-and-operators.md`](../../../../planning/03-algorithms/03-nsga-loop-and-operators.md)
//! §5.2.
//!
//! ## Why every crossover here is expressed as a move
//!
//! A child is produced by applying, to one parent, the edits at the positions
//! where the donor differs. That is not a roundabout way of copying a slice —
//! it is what makes the operator inherit its invariants. Every changed position
//! goes through `CodonSequence::substitute`, so a locked position cannot be
//! crossed over even if the two parents somehow disagree there, and the
//! "restore locked positions from the parents afterwards" step in the §5.2
//! pseudocode becomes an assertion the type system makes for us rather than a
//! defensive pass that could be forgotten.
//!
//! It also means a crossover has an edit list, so a child can be scored through
//! the exact-delta path against its parent when the edit count is small.

use cf_landscape_api::Edit;
use cf_types::{code, CfError, Result};

use crate::genotype::{Genotype, Move, MoveKind};
use crate::rng::SearchRng;

/// Two parents in, two children out.
pub trait CrossoverOperator: Send + Sync + std::fmt::Debug {
    fn name(&self) -> &'static str;

    /// The positions child 1 takes from parent 2, as a boolean mask.
    ///
    /// Expressed as a mask rather than as children so that every operator
    /// shares one implementation of the invariance-preserving construction
    /// below, and a new crossover pattern is a new mask rather than a new
    /// opportunity to write codons directly.
    fn mask(&self, length: usize, rng: &mut SearchRng) -> Vec<bool>;

    fn cross(
        &self,
        first: &Genotype,
        second: &Genotype,
        rng: &mut SearchRng,
    ) -> Result<(Genotype, Genotype)> {
        if first.len() != second.len() {
            return Err(CfError::new(
                code::SRCH_INCOMPATIBLE_PARENTS,
                format!(
                    "cannot cross a {}-codon parent with a {}-codon one",
                    first.len(),
                    second.len()
                ),
            ));
        }
        let mask = self.mask(first.len(), rng);

        let mut to_first = Vec::new();
        let mut to_second = Vec::new();
        for (index, take_from_other) in mask.iter().enumerate() {
            if !take_from_other {
                continue;
            }
            let a = first.codons()[index];
            let b = second.codons()[index];
            if a == b {
                continue;
            }
            to_first.push(Edit::new(index, b));
            to_second.push(Edit::new(index, a));
        }

        Ok((
            first.with_move(&Move::new(MoveKind::Crossover, to_first))?,
            second.with_move(&Move::new(MoveKind::Crossover, to_second))?,
        ))
    }
}

/// Block-uniform crossover: the default.
///
/// §5.2's reasoning, not restated: per-codon uniform crossover breaks about
/// half of all adjacent pairs, which destroys the co-adapted context tier 2
/// spent its budget building. Block lengths with mean ≈ `B` leave most `W`-wide
/// windows intact while still mixing.
#[derive(Clone, Copy, Debug)]
pub struct BlockUniform {
    mean_block: f64,
}

impl BlockUniform {
    pub const DEFAULT_MEAN_BLOCK: f64 = 30.0;

    pub fn new(mean_block: f64) -> Self {
        Self {
            mean_block: mean_block.max(1.0),
        }
    }
}

impl Default for BlockUniform {
    fn default() -> Self {
        Self::new(Self::DEFAULT_MEAN_BLOCK)
    }
}

impl CrossoverOperator for BlockUniform {
    fn name(&self) -> &'static str {
        "block_uniform"
    }

    fn mask(&self, length: usize, rng: &mut SearchRng) -> Vec<bool> {
        let mut mask = vec![false; length];
        let mut index = 0usize;
        let mut take_from_other = rng.coin();
        while index < length {
            let run = rng.geometric(1.0 / self.mean_block, length);
            let end = (index + run).min(length);
            if take_from_other {
                mask[index..end].fill(true);
            }
            take_from_other = !take_from_other;
            index = end;
        }
        mask
    }
}

/// One-point crossover. The right choice when `W = 1`, per §5.2.
#[derive(Clone, Copy, Debug, Default)]
pub struct OnePoint;

impl CrossoverOperator for OnePoint {
    fn name(&self) -> &'static str {
        "one_point"
    }

    fn mask(&self, length: usize, rng: &mut SearchRng) -> Vec<bool> {
        let mut mask = vec![false; length];
        if length < 2 {
            return mask;
        }
        let cut = 1 + rng.below(length - 1);
        mask[cut..].fill(true);
        mask
    }
}

/// Two-point crossover.
#[derive(Clone, Copy, Debug, Default)]
pub struct TwoPoint;

impl CrossoverOperator for TwoPoint {
    fn name(&self) -> &'static str {
        "two_point"
    }

    fn mask(&self, length: usize, rng: &mut SearchRng) -> Vec<bool> {
        let mut mask = vec![false; length];
        if length < 3 {
            return mask;
        }
        let a = 1 + rng.below(length - 1);
        let b = 1 + rng.below(length - 1);
        let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
        mask[lo..hi].fill(true);
        mask
    }
}
