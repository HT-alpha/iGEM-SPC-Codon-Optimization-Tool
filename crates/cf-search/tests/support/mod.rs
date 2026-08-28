//! Fixtures shared by the `cf-search` integration tests.
//!
//! Each test binary compiles this module separately, so a fixture used by one
//! and not another is dead code from that binary's point of view.
#![allow(dead_code)]

use cf_landscape_api::HostContext;
use cf_search::CodonWeights;
use cf_types::{Digest, HostId};

pub fn host() -> HostContext {
    HostContext::new(
        HostId::new("escherichia_coli:511145/BL21").unwrap(),
        Digest::of_bytes(b"host-profile-fixture"),
    )
}

/// A deterministic, non-uniform weight table.
///
/// Non-uniform on purpose: a uniform table would make a weighted-sampling
/// assertion pass for the wrong reason.
pub fn weights() -> CodonWeights {
    let mut table = [0.0f64; 64];
    for (index, slot) in table.iter_mut().enumerate() {
        *slot = ((index as f64 * 17.0) % 23.0 + 1.0) / 23.0;
    }
    let weights = CodonWeights::from_relative_adaptiveness(table).unwrap();
    assert!(!weights.is_uniform());
    weights
}

// --------------------------------------------------------------------------- //
// Shared mock terms.
//
// `Additive` and its call counters started in `tier2.rs`. They moved here when
// `P1-33`'s harness needed the same term: a determinism harness that built its
// own mock would be measuring a different landscape from the one the tier-2
// tests exercise, and the two could then drift apart without either failing.
// --------------------------------------------------------------------------- //

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use cf_landscape_api::provider::{EvalContext, PositionDeltas, Rebased, TermProvider};
use cf_landscape_api::{Capability, InteractionRadius, Provenance};
use cf_types::{CodonSequence, Direction, Result, TermId, Unit};

#[derive(Clone, Debug, Default)]
pub struct Calls(pub Arc<Counters>);

#[derive(Debug, Default)]
pub struct Counters {
    pub scores: AtomicU64,
    pub tables: AtomicU64,
    pub rebases: AtomicU64,
}

impl Calls {
    pub fn scores(&self) -> u64 {
        self.0.scores.load(Ordering::SeqCst)
    }
    pub fn tables(&self) -> u64 {
        self.0.tables.load(Ordering::SeqCst)
    }
    pub fn rebases(&self) -> u64 {
        self.0.rebases.load(Ordering::SeqCst)
    }
}

/// An additive, position-wise L1 term.
#[derive(Debug)]
pub struct Additive {
    term: TermId,
    calls: Calls,
    weight: [f64; 64],
    radius: InteractionRadius,
    capability: Capability,
}

impl Additive {
    pub fn new(term: &str, calls: Calls, salt: u64) -> Self {
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

    pub fn with_radius(mut self, radius: InteractionRadius) -> Self {
        self.radius = radius;
        self
    }

    pub fn score_only(mut self) -> Self {
        self.capability = Capability::L0;
        self
    }

    /// Explicit weights for named codons; everything else is zero.
    ///
    /// So a test can build a scoring landscape with a known shape rather than
    /// hoping pseudo-random weights happen to produce one.
    pub fn with_weights(term: &str, calls: Calls, pairs: &[(&str, f64)]) -> Self {
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
    pub fn preferring(term: &str, calls: Calls, codon: &str) -> Self {
        let wanted = cf_types::Codon::parse(codon).unwrap();
        let mut this = Self::new(term, calls, 0);
        this.weight = [0.0; 64];
        this.weight[wanted.index_usize()] = 1.0;
        this
    }

    pub fn total(&self, sequence: &CodonSequence) -> f64 {
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
