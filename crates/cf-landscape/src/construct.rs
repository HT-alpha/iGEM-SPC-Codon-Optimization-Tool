//! Where a coding sequence sits inside the molecule that will be synthesized.
//!
//! Ticket P1-21.
//!
//! Constraints are properties of the assembled construct, not of the CDS: a
//! restriction site straddling the junction between the insert and the backbone
//! is in the plasmid whether or not it is in the gene. So a term that counts
//! violations has to assemble the construct before it can count anything, and
//! that means knowing the flanks and the offset.
//!
//! Phase 1 handles the single-cistron case, which is what
//! [`CodonSequence`]'s own flanks already describe. A multi-cistron layout
//! comes with the construct context in `P2-xx`; this type is the seam it will
//! extend, and it fails loudly rather than guessing in the meantime.

use cf_types::CodonSequence;

/// A single coding region between two immutable flanks.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConstructLayout {
    five_prime: Vec<u8>,
    three_prime: Vec<u8>,
    codon_count: usize,
}

impl ConstructLayout {
    /// Take the layout from the sequence's own flanks.
    pub fn from_sequence(sequence: &CodonSequence) -> Self {
        Self {
            five_prime: sequence.five_prime_flank().to_vec(),
            three_prime: sequence.three_prime_flank().to_vec(),
            codon_count: sequence.codon_count(),
        }
    }

    /// Nucleotide offset of the coding region within the construct.
    pub fn cds_offset_nt(&self) -> usize {
        self.five_prime.len()
    }

    /// Half-open coding spans, in construct coordinates.
    pub fn cds_spans(&self) -> Vec<(usize, usize)> {
        let start = self.cds_offset_nt();
        vec![(start, start + self.codon_count * 3)]
    }

    /// Half-open nucleotide span of one codon, in construct coordinates.
    pub fn codon_span_nt(&self, codon_index: usize) -> (usize, usize) {
        let start = self.cds_offset_nt() + codon_index * 3;
        (start, start + 3)
    }

    /// The assembled molecule.
    pub fn assemble(&self, sequence: &CodonSequence) -> Vec<u8> {
        let mut out = Vec::with_capacity(
            self.five_prime.len() + self.codon_count * 3 + self.three_prime.len(),
        );
        out.extend_from_slice(&self.five_prime);
        for codon in sequence.codons() {
            out.extend_from_slice(&codon.to_bytes());
        }
        out.extend_from_slice(&self.three_prime);
        out
    }
}
