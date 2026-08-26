//! Overlapping k-mer occurrence counting, with O(k) incremental updates.
//!
//! Ticket P1-15. Contract: `planning/02-schemas/constraint_set.schema.json`
//! (`kmer_repeat`).
//!
//! LEGACY SEMANTICS, PRESERVED DELIBERATELY. `codon_6.0.py` counts overlapping
//! occurrences via its own `count_overlapping` helper, with a comment recording
//! that `str.count` was wrong here. It was: `str.count` is non-overlapping, so
//! it reports `AAAAA` as containing one `AAA` rather than three, and tandem
//! repeats are exactly what a synthesis vendor rejects. The behaviour is
//! carried forward and made explicit as `count_overlapping` in the schema
//! rather than inherited silently — `contract_overlapping_kmer_semantics_are_
//! explicit_in_the_set` asserts the worked example declares it.
//!
//! MASKED POSITIONS. The prototype masks coding regions out of the backbone
//! with `N` and skips any k-mer containing one. That is preserved: a k-mer
//! overlapping a masked position is not counted. It has a consequence worth
//! knowing about, recorded in `docs/legacy-parity.md` — k-mers spanning the
//! junction between a coding region and the backbone are counted by neither
//! side, so a repeat straddling that boundary is invisible.

use std::collections::HashMap;

use cf_types::error::{code, CfError, Result};

/// Two-bit packed k-mer. `k <= 32` by construction (see `ConstraintSet::validate`).
pub type PackedKmer = u64;

#[inline]
fn base_code(byte: u8) -> Option<u64> {
    match byte {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

/// Reverse complement of a packed k-mer.
fn revcomp(packed: PackedKmer, k: usize) -> PackedKmer {
    let mut out = 0u64;
    let mut value = packed;
    for _ in 0..k {
        let base = value & 3;
        out = (out << 2) | (3 - base);
        value >>= 2;
    }
    out
}

/// Pack the k-mer starting at `pos`, or `None` if it runs off the end or
/// touches a masked position.
#[inline]
pub fn pack_at(seq: &[u8], pos: usize, k: usize) -> Option<PackedKmer> {
    if pos + k > seq.len() {
        return None;
    }
    let mut packed = 0u64;
    for &byte in &seq[pos..pos + k] {
        packed = (packed << 2) | base_code(byte)?;
    }
    Some(packed)
}

/// Occurrence counts for one `k`, maintained incrementally.
#[derive(Clone, Debug)]
pub struct KmerIndex {
    k: usize,
    canonical: bool,
    counts: HashMap<PackedKmer, u32>,
    /// The sequence the index describes. Held so an update knows what changed.
    seq: Vec<u8>,
}

impl KmerIndex {
    pub fn build(seq: &[u8], k: usize, canonical: bool) -> Result<Self> {
        if !(4..=32).contains(&k) {
            return Err(CfError::new(
                code::FEAS_CONSTRAINT_INVALID,
                format!("k={k} is outside the supported range 4..=32"),
            ));
        }
        let mut index = Self {
            k,
            canonical,
            counts: HashMap::new(),
            seq: seq.to_vec(),
        };
        index.rebuild();
        Ok(index)
    }

    fn rebuild(&mut self) {
        self.counts.clear();
        if self.seq.len() < self.k {
            return;
        }
        for pos in 0..=self.seq.len() - self.k {
            if let Some(key) = self.key_at(pos) {
                *self.counts.entry(key).or_insert(0) += 1;
            }
        }
    }

    #[inline]
    fn key_at(&self, pos: usize) -> Option<PackedKmer> {
        let packed = pack_at(&self.seq, pos, self.k)?;
        Some(self.canonicalize(packed))
    }

    #[inline]
    fn canonicalize(&self, packed: PackedKmer) -> PackedKmer {
        if self.canonical {
            packed.min(revcomp(packed, self.k))
        } else {
            packed
        }
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn sequence(&self) -> &[u8] {
        &self.seq
    }

    pub fn count(&self, packed: PackedKmer) -> u32 {
        *self.counts.get(&self.canonicalize(packed)).unwrap_or(&0)
    }

    pub fn count_at(&self, pos: usize) -> u32 {
        self.key_at(pos).map(|key| self.count(key)).unwrap_or(0)
    }

    /// Every k-mer occurring more than `max_occurrences` times, with the
    /// positions at which it occurs.
    pub fn over_limit(&self, max_occurrences: u32) -> Vec<(PackedKmer, u32, Vec<usize>)> {
        let mut offenders: Vec<(PackedKmer, u32)> = self
            .counts
            .iter()
            .filter(|(_, &count)| count > max_occurrences)
            .map(|(&key, &count)| (key, count))
            .collect();
        offenders.sort();
        offenders
            .into_iter()
            .map(|(key, count)| {
                let positions = (0..self.seq.len().saturating_sub(self.k - 1))
                    .filter(|&pos| self.key_at(pos) == Some(key))
                    .collect();
                (key, count, positions)
            })
            .collect()
    }

    /// Replace `replacement.len()` bytes at `offset`, updating counts in O(k).
    ///
    /// A substitution of `n` bytes at `offset` invalidates exactly the k-mers
    /// starting in `[offset - k + 1, offset + n - 1]`, so only those are
    /// removed and re-added.
    pub fn substitute(&mut self, offset: usize, replacement: &[u8]) -> Result<()> {
        let n = replacement.len();
        if offset + n > self.seq.len() {
            return Err(CfError::new(
                code::SEQ_INDEX_OUT_OF_RANGE,
                format!(
                    "substitution at {offset}..{} is outside a sequence of {} nt",
                    offset + n,
                    self.seq.len()
                ),
            ));
        }
        let lo = offset.saturating_sub(self.k - 1);
        let hi = (offset + n).min(self.seq.len().saturating_sub(self.k - 1));

        for pos in lo..hi {
            if let Some(key) = self.key_at(pos) {
                if let Some(count) = self.counts.get_mut(&key) {
                    *count -= 1;
                    if *count == 0 {
                        self.counts.remove(&key);
                    }
                }
            }
        }
        self.seq[offset..offset + n].copy_from_slice(replacement);
        for pos in lo..hi {
            if let Some(key) = self.key_at(pos) {
                *self.counts.entry(key).or_insert(0) += 1;
            }
        }
        Ok(())
    }

    /// Total distinct k-mers tracked. Used by tests to compare against a rebuild.
    pub fn distinct(&self) -> usize {
        self.counts.len()
    }

    pub fn counts(&self) -> &HashMap<PackedKmer, u32> {
        &self.counts
    }
}

/// Render a packed k-mer back to text, for error messages.
pub fn unpack(packed: PackedKmer, k: usize) -> String {
    let mut out = vec![0u8; k];
    let mut value = packed;
    for i in (0..k).rev() {
        out[i] = b"ACGT"[(value & 3) as usize];
        value >>= 2;
    }
    String::from_utf8(out).expect("ASCII")
}

/// The prototype's `count_overlapping`, transcribed for parity testing.
///
/// Kept in the shipped crate rather than in a test so that the legacy
/// semantics have exactly one definition, and so `P1-27`'s differential test
/// and the checker cannot drift apart by construction.
pub fn count_overlapping(haystack: &[u8], needle: &[u8]) -> u32 {
    if needle.is_empty() || haystack.len() < needle.len() {
        return 0;
    }
    (0..=haystack.len() - needle.len())
        .filter(|&i| &haystack[i..i + needle.len()] == needle)
        .count() as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_overlapping_counting_matches_the_prototype() {
        // The case that motivated the prototype's own fix: str.count would say 1.
        assert_eq!(count_overlapping(b"AAAAA", b"AAA"), 3);
        assert_eq!(count_overlapping(b"ATATAT", b"ATAT"), 2);
        assert_eq!(count_overlapping(b"ACGT", b"TGCA"), 0);

        let index = KmerIndex::build(b"AAAAAAA", 4, false).unwrap();
        assert_eq!(index.count(pack_at(b"AAAA", 0, 4).unwrap()), 4);
    }

    #[test]
    fn unit_index_agrees_with_the_naive_count() {
        let seq = b"ATGAAACGTATTCTGGATGAAGGCGCAACCGTTAGCCTGTAAATGAAACGTATTCTG";
        for k in [4usize, 6, 8] {
            let index = KmerIndex::build(seq, k, false).unwrap();
            for pos in 0..=seq.len() - k {
                let expected = count_overlapping(seq, &seq[pos..pos + k]);
                assert_eq!(index.count_at(pos), expected, "k={k} pos={pos}");
            }
        }
    }

    #[test]
    fn unit_masked_positions_are_skipped() {
        // The prototype masks coding regions with N and skips k-mers touching
        // one; a k-mer spanning the mask boundary is counted by neither side.
        let index = KmerIndex::build(b"AAAANNNNAAAA", 4, false).unwrap();
        assert_eq!(index.count(pack_at(b"AAAA", 0, 4).unwrap()), 2);
        assert_eq!(index.distinct(), 1);
    }

    #[test]
    fn unit_canonical_strand_folds_reverse_complements() {
        let plain = KmerIndex::build(b"AAAATTTT", 4, false).unwrap();
        assert_eq!(plain.count(pack_at(b"AAAA", 0, 4).unwrap()), 1);
        assert_eq!(plain.count(pack_at(b"TTTT", 0, 4).unwrap()), 1);

        let canonical = KmerIndex::build(b"AAAATTTT", 4, true).unwrap();
        // AAAA and TTTT are reverse complements, so they share one bucket.
        assert_eq!(canonical.count(pack_at(b"AAAA", 0, 4).unwrap()), 2);
        assert_eq!(canonical.count(pack_at(b"TTTT", 0, 4).unwrap()), 2);
    }

    /// TST-CONS-001 (k-mer portion): incremental equals a naive rebuild.
    #[test]
    fn property_tst_cons_001_incremental_matches_rebuild() {
        // A small deterministic LCG keeps the test dependency-free and
        // reproducible; the point is coverage of edit positions, not randomness
        // quality.
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };

        let alphabet = b"ACGT";
        let base: Vec<u8> = (0..300).map(|i| alphabet[(i * 7 + i / 4) % 4]).collect();

        for &k in &[4usize, 8, 12] {
            for &canonical in &[false, true] {
                let mut index = KmerIndex::build(&base, k, canonical).unwrap();
                let mut mirror = base.clone();

                for _ in 0..2_000 {
                    let offset = (next() as usize) % (mirror.len() - 3);
                    let replacement: Vec<u8> =
                        (0..3).map(|_| alphabet[(next() as usize) % 4]).collect();
                    index.substitute(offset, &replacement).unwrap();
                    mirror[offset..offset + 3].copy_from_slice(&replacement);

                    let rebuilt = KmerIndex::build(&mirror, k, canonical).unwrap();
                    assert_eq!(
                        index.counts(),
                        rebuilt.counts(),
                        "k={k} canonical={canonical} offset={offset}: incremental diverged"
                    );
                }
            }
        }
    }

    #[test]
    fn unit_over_limit_reports_offenders_and_positions() {
        let index = KmerIndex::build(b"ACGTACGTACGT", 4, false).unwrap();
        let offenders = index.over_limit(2);
        assert_eq!(offenders.len(), 1);
        let (key, count, positions) = &offenders[0];
        assert_eq!(unpack(*key, 4), "ACGT");
        assert_eq!(*count, 3);
        assert_eq!(positions, &vec![0, 4, 8]);
        assert!(index.over_limit(3).is_empty());
    }

    #[test]
    fn unit_pack_round_trips() {
        for seq in ["ACGT", "AAAA", "TTTT", "GATTACA"] {
            let k = seq.len();
            let packed = pack_at(seq.as_bytes(), 0, k).unwrap();
            assert_eq!(unpack(packed, k), seq);
        }
        assert_eq!(pack_at(b"ACNT", 0, 4), None);
        assert_eq!(pack_at(b"ACG", 0, 4), None);
    }

    #[test]
    fn unit_bad_k_is_refused() {
        assert_eq!(
            KmerIndex::build(b"ACGT", 3, false).unwrap_err().code,
            code::FEAS_CONSTRAINT_INVALID
        );
        assert!(KmerIndex::build(b"ACGT", 33, false).is_err());
    }
}
