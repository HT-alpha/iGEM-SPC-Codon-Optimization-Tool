//! The L0/L1/L2 capability model of ADR-0002.
//!
//! Ticket P1-20.
//!
//! The search engine **requires** L0 and **exploits** L1. It never requires L2.
//! That sentence is the whole design, and this module's job is to make the
//! second and third clauses structural rather than aspirational:
//!
//! * L1 is exploited, not required, so [`Capability::L0`] is a legal answer and
//!   composition must have a correct path for it.
//! * L2 is never required, so [`super::TermProvider`] — the only trait
//!   composition knows about — has no gradient method at all. A gradient-only
//!   provider cannot be written against this contract, which is what makes the
//!   no-gradients decision impossible to erode by increments.

use serde::{Deserialize, Serialize};

/// What a provider can compute for one term.
///
/// Ordered: `L0 < L1 < L2`. The ordering is meaningful — a provider offering a
/// higher capability offers every lower one — and is used when planning a tier.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, Serialize, Deserialize)]
pub enum Capability {
    /// Score a full sequence. Mandatory for every provider.
    #[serde(rename = "L0")]
    L0,
    /// Exact change in the term for one synonymous substitution, and the full
    /// table of such changes over a set of positions. Available for every
    /// decomposable term.
    #[serde(rename = "L1")]
    L1,
    /// Derivative with respect to relaxed codon logits.
    ///
    /// Representable because the response schema admits it and a provider may
    /// legitimately offer it to some other consumer. Not reachable through
    /// [`super::TermProvider`], and composition treats an L2 declaration
    /// exactly as it treats L1.
    #[serde(rename = "L2")]
    L2,
}

impl Capability {
    pub fn as_str(self) -> &'static str {
        match self {
            Capability::L0 => "L0",
            Capability::L1 => "L1",
            Capability::L2 => "L2",
        }
    }

    pub fn parse(text: &str) -> Option<Self> {
        match text {
            "L0" => Some(Capability::L0),
            "L1" => Some(Capability::L1),
            "L2" => Some(Capability::L2),
            _ => None,
        }
    }

    /// Whether this provider can answer a delta query in closed form.
    ///
    /// The single question composition asks. It is deliberately not "is this
    /// L1" — L2 answers yes as well, and a future L3 would too.
    #[inline]
    pub fn offers_exact_delta(self) -> bool {
        self >= Capability::L1
    }
}

impl std::fmt::Display for Capability {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// How far a term's contribution at one position reaches.
///
/// Tier 2's red-black partition rests on same-colour blocks being
/// conditionally independent, which holds when the block size is at least the
/// largest interaction radius among the scored terms — not merely at least the
/// epistatic window `W`. `02-block-local-coordinate-ascent.md` §4 says as much
/// by keeping a fallback for "a term with an effective interaction radius
/// larger than `B`", so a radius above `W` is contemplated, not hypothetical.
///
/// Declared per provider and required, rather than defaulted. A default would
/// be a claim the provider never made, and the consequence of getting it wrong
/// is invisible: blocks scored as independent while genuinely coupled produce
/// a plausible answer that is not the coordinate-ascent step it reports being.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub enum InteractionRadius {
    /// The contribution decomposes by position: choosing a codon at `i` changes
    /// nothing at any other position. CAI and GC content are both this.
    Positionwise,
    /// Coupling reaches at most this many codons either side.
    Codons(usize),
    /// No bound. Long-range RNA base pairing is the real instance — two codons
    /// 400 nt apart can pair — and it is why §4's overlapping-block fallback
    /// exists at all.
    Unbounded,
}

impl InteractionRadius {
    /// The radius in codons, or `None` when unbounded.
    pub fn codons(self) -> Option<usize> {
        match self {
            InteractionRadius::Positionwise => Some(0),
            InteractionRadius::Codons(n) => Some(n),
            InteractionRadius::Unbounded => None,
        }
    }

    /// Whether a block of `block_size` codons isolates this term.
    pub fn fits_block(self, block_size: usize) -> bool {
        self.codons().is_some_and(|radius| radius <= block_size)
    }

    pub fn as_str(self) -> &'static str {
        match self {
            InteractionRadius::Positionwise => "positionwise",
            InteractionRadius::Codons(_) => "bounded",
            InteractionRadius::Unbounded => "unbounded",
        }
    }
}

impl std::fmt::Display for InteractionRadius {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            InteractionRadius::Positionwise => f.write_str("positionwise"),
            InteractionRadius::Codons(n) => write!(f, "{n} codons"),
            InteractionRadius::Unbounded => f.write_str("unbounded"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_capability_is_ordered_and_round_trips() {
        assert!(Capability::L0 < Capability::L1);
        assert!(Capability::L1 < Capability::L2);
        for c in [Capability::L0, Capability::L1, Capability::L2] {
            assert_eq!(Capability::parse(c.as_str()), Some(c));
        }
        assert_eq!(Capability::parse("L3"), None);
    }

    #[test]
    fn unit_exact_delta_is_a_threshold_not_an_equality() {
        assert!(!Capability::L0.offers_exact_delta());
        assert!(Capability::L1.offers_exact_delta());
        assert!(Capability::L2.offers_exact_delta());
    }
}
