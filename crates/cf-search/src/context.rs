//! The context-model seam, and codon weights.
//!
//! Ticket P1-23, departure `DEP-SRCH-005`.
//!
//! Two of the four mutation operators in
//! [`03-nsga-loop-and-operators.md`](../../../planning/03-algorithms/03-nsga-loop-and-operators.md)
//! §5.1 draw their replacement from the Potts conditional: `potts_guided` at
//! `T = 1.5`, and `block_gibbs` one position at a time across a block. That
//! model is `cf-epistasis` and it arrives in `P3-05`, and `cf-search` is
//! specified as talking only to `cf-landscape-api` — so even once it exists,
//! the conditional has to reach the operators through a seam.
//!
//! [`ContextModel`] is that seam. The operators are written against it now and
//! are registered only when a model is supplied, so what ships this round is
//! two working samplers with no distribution to sample from rather than two
//! operators pretending to be context-aware. The distinction matters because a
//! `potts_guided` backed by host codon usage is `uniform_synonymous` with
//! different weights, under a name that claims otherwise — and a search whose
//! diversity behaviour silently changed would be attributed to anything but
//! that.

use cf_types::{code, CfError, Codon, CodonSequence, Result, CODON_COUNT};

/// A conditional distribution over the synonyms at one position.
///
/// Implemented by the epistatic surface in `P3-05`. The window size is
/// reported so the search can assert the tier-2 block partition satisfies
/// `B ≥ W` (`P1-25`) without knowing what model it is talking to.
pub trait ContextModel: Send + Sync {
    /// Name recorded in provenance and in the operator-availability report.
    fn name(&self) -> &str;

    /// Interaction window, in codons. `1` means context-free.
    fn window(&self) -> usize;

    /// Fill `out` with `(codon, probability)` over the synonyms at `index`.
    ///
    /// The implementation is responsible for normalizing and for applying
    /// `temperature`. `out` is cleared first and is reused across calls, so a
    /// Gibbs sweep does not allocate per position.
    fn conditional(
        &self,
        sequence: &CodonSequence,
        index: usize,
        temperature: f64,
        out: &mut Vec<(Codon, f64)>,
    ) -> Result<()>;
}

/// Per-codon sampling weights, indexed by [`Codon::index_usize`].
///
/// Carries relative adaptiveness when a host profile supplies it. The table is
/// injected as data rather than reached for through `cf-landscape`: weights are
/// values, and a `cf-search → cf-landscape` edge would buy nothing except a
/// dependency the module graph does not declare.
#[derive(Clone, Debug)]
pub struct CodonWeights([f64; CODON_COUNT]);

impl CodonWeights {
    /// Every synonym equally likely.
    ///
    /// The explicit fallback when no host table is supplied, not a stand-in for
    /// one: a uniform draw over a synonymous family is a real, statable
    /// sampling policy, whereas inventing plausible-looking adaptiveness values
    /// would not be.
    pub fn uniform() -> Self {
        Self([1.0; CODON_COUNT])
    }

    /// Relative adaptiveness from a host profile's `codon_usage.table`.
    ///
    /// Rejects a negative or non-finite weight (`CF-HOST-002`). A zero weight
    /// is accepted: it means the reference set never used that codon, and an
    /// operator that never proposes it is the correct behaviour.
    pub fn from_relative_adaptiveness(weights: [f64; CODON_COUNT]) -> Result<Self> {
        for (index, weight) in weights.iter().enumerate() {
            if !weight.is_finite() || *weight < 0.0 {
                return Err(CfError::new(
                    code::HOST_PROFILE_INVALID,
                    format!(
                        "relative adaptiveness for codon {} is {weight}, which is not a \
                         usable weight",
                        Codon::from_index(index as u8)
                            .map(|c| c.to_string_upper())
                            .unwrap_or_else(|_| index.to_string())
                    ),
                ));
            }
        }
        Ok(Self(weights))
    }

    pub fn weight(&self, codon: Codon) -> f64 {
        self.0[codon.index_usize()]
    }

    /// Weights for a slice of candidate codons, in the same order.
    pub fn weights_for(&self, codons: &[Codon]) -> Vec<f64> {
        codons.iter().map(|c| self.weight(*c)).collect()
    }

    /// Whether any weight differs from any other.
    ///
    /// Used by tests to assert a fixture is not accidentally uniform, which
    /// would make a weighted-sampling assertion pass for the wrong reason.
    pub fn is_uniform(&self) -> bool {
        let first = self.0[0];
        self.0.iter().all(|w| *w == first)
    }
}

impl Default for CodonWeights {
    fn default() -> Self {
        Self::uniform()
    }
}
