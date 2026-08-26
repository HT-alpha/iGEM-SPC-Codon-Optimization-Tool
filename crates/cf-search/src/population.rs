//! Duplicate suppression and seed injection.
//!
//! Ticket P1-23. Specification:
//! [`03-nsga-loop-and-operators.md`](../../../planning/03-algorithms/03-nsga-loop-and-operators.md)
//! §5.3.
//!
//! Duplicate-free populations matter more here than in a continuous MOEA
//! because the space is finite and tier 2 deliberately converges many instances
//! to similar local optima. A population that quietly filled with copies would
//! keep reporting a healthy size while searching a fraction of it.

use std::collections::HashMap;

use cf_types::{code, CfError, Result};

use crate::genotype::Genotype;
use crate::operators::{MutationContext, OperatorMix};
use crate::rng::SearchRng;

/// A set of individuals with no exact duplicates.
///
/// Indexed by the genotype's incremental hash, with a full codon comparison to
/// confirm. The confirmation makes "a distinct individual is never discarded"
/// unconditional rather than probabilistic, and it costs one slice comparison
/// on a bucket that almost always holds one entry.
///
/// It is worth being precise about what it does *not* buy, because the comment
/// this replaces overstated it. Removing the confirmation changes no observable
/// behaviour for any input a test can construct: an exhaustive search over
/// every two-position difference in the 14×64 space finds no collision in
/// `Genotype`'s hash, so the branch is unreachable. What it defends against is
/// a future weaker hash or a second mutation path leaving one stale — and the
/// property that actually catches a stale hash is
/// `Genotype::hash_is_consistent`, which `TST-SRCH-006` asserts.
#[derive(Debug, Default)]
pub struct Population {
    members: Vec<Genotype>,
    by_hash: HashMap<u64, Vec<usize>>,
}

impl Population {
    /// Regeneration attempts before forcing an extra mutation, per §5.3.
    pub const MAX_ATTEMPTS: usize = 8;

    pub fn new() -> Self {
        Self::default()
    }

    pub fn len(&self) -> usize {
        self.members.len()
    }

    pub fn is_empty(&self) -> bool {
        self.members.is_empty()
    }

    pub fn members(&self) -> &[Genotype] {
        &self.members
    }

    pub fn contains(&self, candidate: &Genotype) -> bool {
        self.by_hash.get(&candidate.hash()).is_some_and(|slots| {
            slots
                .iter()
                .any(|i| self.members[*i].same_codons_as(candidate))
        })
    }

    /// Add an individual if it is not already present.
    ///
    /// Returns whether it was added, rather than silently doing nothing, so a
    /// caller counting population size cannot drift from reality.
    pub fn insert(&mut self, candidate: Genotype) -> bool {
        if self.contains(&candidate) {
            return false;
        }
        self.by_hash
            .entry(candidate.hash())
            .or_default()
            .push(self.members.len());
        self.members.push(candidate);
        true
    }

    /// Produce an offspring distinct from every current member.
    ///
    /// Up to [`Population::MAX_ATTEMPTS`] fresh proposals, then a forced extra
    /// mutation on the last one. The forced mutation can itself collide — the
    /// space may genuinely be exhausted near the parent — so the result is
    /// returned with a flag rather than looped on forever.
    pub fn breed_distinct(
        &self,
        parent: &Genotype,
        mix: &OperatorMix,
        ctx: &MutationContext<'_>,
        rng: &mut SearchRng,
    ) -> Result<Offspring> {
        let mut last: Option<Genotype> = None;

        for attempt in 0..Self::MAX_ATTEMPTS {
            let mv = mix.propose(ctx, rng)?;
            let child = parent.with_move(&mv.compacted(parent))?;
            if !self.contains(&child) && !child.same_codons_as(parent) {
                return Ok(Offspring {
                    genotype: child,
                    attempts: attempt + 1,
                    forced: false,
                });
            }
            last = Some(child);
        }

        let base = last.unwrap_or_else(|| parent.clone());
        let forced_ctx = MutationContext {
            genotype: &base,
            ..*ctx
        };
        let extra = mix.propose(&forced_ctx, rng)?;
        let child = base.with_move(&extra.compacted(&base))?;
        Ok(Offspring {
            genotype: child,
            attempts: Self::MAX_ATTEMPTS + 1,
            forced: true,
        })
    }

    /// Inject imported seeds, capped as a fraction of the population.
    ///
    /// ADR-0007 and §5.3. Each accepted seed brings two mutated copies so
    /// crossover has material immediately — a lone seed in a population of
    /// unrelated individuals contributes nothing until it happens to be
    /// selected twice.
    ///
    /// Seeds are injected at initialization only. A seed arriving mid-run would
    /// be an unrecorded input to a job whose lockfile has already been sealed.
    pub fn inject_seeds(
        &mut self,
        seeds: &[Genotype],
        target_population: usize,
        max_fraction: f64,
        mix: &OperatorMix,
        weights: &crate::context::CodonWeights,
        rng: &mut SearchRng,
    ) -> Result<SeedReport> {
        if !(0.0..=1.0).contains(&max_fraction) {
            return Err(CfError::new(
                code::SRCH_SEED_POLICY,
                format!("max_fraction_of_population is {max_fraction}, outside [0, 1]"),
            ));
        }
        let cap = (target_population as f64 * max_fraction).floor() as usize;
        let mut report = SeedReport {
            offered: seeds.len(),
            accepted: 0,
            duplicates: 0,
            capped: seeds.len().saturating_sub(cap),
            mutated_copies: 0,
        };

        for seed in seeds.iter().take(cap) {
            if !self.insert(seed.clone()) {
                report.duplicates += 1;
                continue;
            }
            report.accepted += 1;
            for _ in 0..2 {
                let ctx = MutationContext::new(seed, weights);
                let offspring = self.breed_distinct(seed, mix, &ctx, rng)?;
                if self.insert(offspring.genotype) {
                    report.mutated_copies += 1;
                }
            }
        }
        Ok(report)
    }
}

/// An offspring and how much work it took.
#[derive(Clone, Debug)]
pub struct Offspring {
    pub genotype: Genotype,
    pub attempts: usize,
    /// Whether the extra mutation of §5.3 had to be forced.
    pub forced: bool,
}

/// What seed injection did.
///
/// Reported rather than logged: a user whose seeds were mostly rejected as
/// duplicates, or truncated by the cap, has no way to discover that otherwise,
/// and would reasonably conclude the importer had ignored them.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SeedReport {
    pub offered: usize,
    pub accepted: usize,
    pub duplicates: usize,
    pub capped: usize,
    pub mutated_copies: usize,
}
