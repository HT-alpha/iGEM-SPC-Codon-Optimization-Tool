//! The tRNA demand ledger.
//!
//! Ticket P2-12. Specification: `05-burden-penalty-and-capacity.md` §3.
//!
//! ## Two quantities, and only one of them is conserved
//!
//! The ledger tracks demand at two levels, and the difference is what
//! `TST-BURD-001` is about.
//!
//! **Amino-acid demand is conserved exactly** under any synonymous
//! substitution, because a codon consumes exactly one charged tRNA carrying its
//! amino acid whichever synonym is used. It is therefore counted as an integer
//! and is bit-identical across any number of substitutions — not "equal to
//! within tolerance", equal.
//!
//! **Anticodon-pool demand is not conserved**, and must not be: synonyms are
//! read by different pools, and moving demand between pools is the entire point
//! of codon optimization for burden. It is a weighted sum and is compared
//! against a full recomputation to `1e-12` (`TST-BURD-002`).
//!
//! Conflating the two is the mistake worth guarding against: a ledger that
//! summed only pools would show conservation failing for a correct
//! substitution, and one that summed only amino acids would show a burden
//! objective that no codon choice can move.

use std::collections::BTreeMap;

use cf_types::{code, CfError, Codon, GeneticCode, Result};

use crate::wobble::{decoders, DecodingWeights};

/// Demand accumulated over a construct.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DemandLedger {
    /// Codon counts by amino acid. Integer, so conservation is exact.
    amino_acid: BTreeMap<u8, u64>,
    /// Weighted demand by anticodon. Float, and deliberately not conserved.
    anticodon: BTreeMap<[u8; 3], f64>,
    /// Codon counts, retained so a full recomputation has something to
    /// recompute from without re-reading the sequence.
    codons: BTreeMap<Codon, u64>,
}

impl DemandLedger {
    pub fn new() -> Self {
        Self::default()
    }

    /// Accumulate a whole coding sequence.
    ///
    /// `flux` is the cistron's expression flux `f_k` from §3 — demand scales
    /// with how much of the protein is made, not just with its sequence.
    pub fn accumulate(
        &mut self,
        codons: &[Codon],
        flux: f64,
        weights: &DecodingWeights,
        table: &GeneticCode,
    ) -> Result<()> {
        for codon in codons {
            self.add(*codon, flux, weights, table)?;
        }
        Ok(())
    }

    fn add(
        &mut self,
        codon: Codon,
        flux: f64,
        weights: &DecodingWeights,
        table: &GeneticCode,
    ) -> Result<()> {
        *self.amino_acid.entry(table.amino_acid(codon)).or_insert(0) += 1;
        *self.codons.entry(codon).or_insert(0) += 1;
        for decoder in decoders(codon) {
            let weight = weights.weight(codon, &decoder)?;
            *self.anticodon.entry(decoder.anticodon).or_insert(0.0) += flux * weight;
        }
        Ok(())
    }

    fn remove(
        &mut self,
        codon: Codon,
        flux: f64,
        weights: &DecodingWeights,
        table: &GeneticCode,
    ) -> Result<()> {
        let amino = table.amino_acid(codon);
        match self.amino_acid.get_mut(&amino) {
            Some(count) if *count > 0 => *count -= 1,
            _ => {
                return Err(CfError::new(
                    code::BURD_LEDGER_UNDERFLOW,
                    format!(
                        "cannot remove codon {codon} from the ledger: no demand for amino \
                         acid {} is recorded. The ledger is being asked to un-count \
                         something it never counted, which means it has drifted from the \
                         sequence it describes",
                        amino as char
                    ),
                ));
            }
        }
        if let Some(count) = self.codons.get_mut(&codon) {
            *count -= 1;
        }
        for decoder in decoders(codon) {
            let weight = weights.weight(codon, &decoder)?;
            if let Some(demand) = self.anticodon.get_mut(&decoder.anticodon) {
                *demand -= flux * weight;
            }
        }
        Ok(())
    }

    /// Apply one synonymous substitution incrementally. §3's `DELTA_BURDEN`.
    ///
    /// Touches at most six pools — three that decoded the old codon and three
    /// that decode the new one — rather than re-reading the construct, which is
    /// what makes the burden term affordable inside the search loop.
    ///
    /// The substitution is checked for synonymy rather than assumed: a ledger
    /// updated through a non-synonymous change would silently stop describing
    /// the protein it was built for, and every conservation claim below would
    /// be false without anything reporting it.
    pub fn substitute(
        &mut self,
        from: Codon,
        to: Codon,
        flux: f64,
        weights: &DecodingWeights,
        table: &GeneticCode,
    ) -> Result<()> {
        if table.amino_acid(from) != table.amino_acid(to) {
            return Err(CfError::new(
                code::BURD_NON_SYNONYMOUS,
                format!(
                    "{from} encodes {} and {to} encodes {}; a ledger updated through a \
                     non-synonymous substitution stops describing the protein it was \
                     built for",
                    table.amino_acid(from) as char,
                    table.amino_acid(to) as char
                ),
            ));
        }
        self.remove(from, flux, weights, table)?;
        self.add(to, flux, weights, table)
    }

    /// Total codon count per amino acid. Exactly conserved under synonymy.
    pub fn amino_acid_demand(&self) -> &BTreeMap<u8, u64> {
        &self.amino_acid
    }

    /// Weighted demand per anticodon pool.
    pub fn anticodon_demand(&self) -> &BTreeMap<[u8; 3], f64> {
        &self.anticodon
    }

    pub fn codon_counts(&self) -> &BTreeMap<Codon, u64> {
        &self.codons
    }

    /// Rebuild from scratch, for comparison against the incremental path.
    ///
    /// Sums in codon order rather than in the order edits happened, so the
    /// comparison in `TST-BURD-002` is against a *canonical* accumulation
    /// rather than against whichever order the caller used. Floating-point
    /// addition is not associative, so "the full recomputation" has to mean
    /// one specific order or the test is comparing two arbitrary ones.
    pub fn recompute(
        counts: &BTreeMap<Codon, u64>,
        flux: f64,
        weights: &DecodingWeights,
        table: &GeneticCode,
    ) -> Result<Self> {
        let mut ledger = Self::new();
        for (codon, count) in counts {
            for _ in 0..*count {
                ledger.add(*codon, flux, weights, table)?;
            }
        }
        Ok(ledger)
    }

    /// Largest absolute disagreement in pool demand against another ledger.
    pub fn max_pool_difference(&self, other: &Self) -> f64 {
        let mut worst = 0.0f64;
        let pools: std::collections::BTreeSet<&[u8; 3]> = self
            .anticodon
            .keys()
            .chain(other.anticodon.keys())
            .collect();
        for pool in pools {
            let a = self.anticodon.get(pool).copied().unwrap_or(0.0);
            let b = other.anticodon.get(pool).copied().unwrap_or(0.0);
            worst = worst.max((a - b).abs());
        }
        worst
    }
}
