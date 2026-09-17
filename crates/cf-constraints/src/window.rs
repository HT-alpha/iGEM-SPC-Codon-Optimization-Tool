//! Homopolymer runs and GC-window content, with incremental updates.
//!
//! Ticket P1-16. Contract: `planning/02-schemas/constraint_set.schema.json`
//! (`homopolymer`, `gc_window`).
//!
//! Homopolymers are found by a run scan; after a substitution only a window of
//! `max_run` on each side can change, because a run crossing the edit is
//! bounded by the runs adjacent to it.
//!
//! GC content uses a Fenwick tree rather than a prefix-sum array. A prefix-sum
//! array gives O(1) window queries but O(n) updates, which is the wrong trade
//! for a structure the repair loop mutates: the Fenwick tree is O(log n) for
//! both, and n here is a construct length, so log n is about 14.

/// A maximal run of one base.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Run {
    pub base: u8,
    pub start: usize,
    pub end: usize,
}

impl Run {
    pub fn len(&self) -> usize {
        self.end - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.end == self.start
    }
}

/// Every maximal run at least `min_len` long, over the whole sequence.
///
/// Masked positions (anything outside ACGT) break runs rather than extending
/// them, matching the treatment in the k-mer index.
pub fn find_runs(seq: &[u8], min_len: usize, base_filter: Option<u8>) -> Vec<Run> {
    let mut runs = Vec::new();
    if seq.is_empty() || min_len == 0 {
        return runs;
    }
    let mut start = 0usize;
    while start < seq.len() {
        let base = seq[start].to_ascii_uppercase();
        if !matches!(base, b'A' | b'C' | b'G' | b'T') {
            start += 1;
            continue;
        }
        let mut end = start + 1;
        while end < seq.len() && seq[end].to_ascii_uppercase() == base {
            end += 1;
        }
        if end - start >= min_len && base_filter.map(|b| b == base).unwrap_or(true) {
            runs.push(Run { base, start, end });
        }
        start = end;
    }
    runs
}

/// Runs overlapping the window that a substitution at `offset..offset+n` can
/// affect.
///
/// The window extends `max_run` beyond the edit on each side: a run that ends
/// just before the edit can be extended by it, and no run longer than the
/// window can be created without also being detected inside it.
pub fn find_runs_near(
    seq: &[u8],
    offset: usize,
    n: usize,
    min_len: usize,
    base_filter: Option<u8>,
) -> Vec<Run> {
    let lo = offset.saturating_sub(min_len);
    let hi = (offset + n + min_len).min(seq.len());
    // Extend to full run boundaries so a run is never reported truncated.
    let mut lo = lo;
    while lo > 0 && seq[lo - 1].eq_ignore_ascii_case(&seq[lo]) {
        lo -= 1;
    }
    let mut hi = hi;
    while hi < seq.len() && seq[hi - 1].eq_ignore_ascii_case(&seq[hi]) {
        hi += 1;
    }
    find_runs(&seq[lo..hi], min_len, base_filter)
        .into_iter()
        .map(|run| Run {
            base: run.base,
            start: run.start + lo,
            end: run.end + lo,
        })
        .collect()
}

/// Fenwick (binary indexed) tree over a 0/1 GC indicator.
#[derive(Clone, Debug)]
pub struct GcIndex {
    tree: Vec<u32>,
    seq: Vec<u8>,
}

#[inline]
fn is_gc(byte: u8) -> bool {
    matches!(byte, b'G' | b'C' | b'g' | b'c')
}

impl GcIndex {
    pub fn build(seq: &[u8]) -> Self {
        let mut index = Self {
            tree: vec![0; seq.len() + 1],
            seq: seq.to_vec(),
        };
        for (i, &byte) in seq.iter().enumerate() {
            if is_gc(byte) {
                index.add(i, 1);
            }
        }
        index
    }

    fn add(&mut self, position: usize, delta: i64) {
        let mut i = position + 1;
        while i < self.tree.len() {
            self.tree[i] = (self.tree[i] as i64 + delta) as u32;
            i += i & i.wrapping_neg();
        }
    }

    /// GC count in `[0, end)`.
    fn prefix(&self, end: usize) -> u32 {
        let mut i = end;
        let mut total = 0u32;
        while i > 0 {
            total += self.tree[i];
            i -= i & i.wrapping_neg();
        }
        total
    }

    /// GC count in `[start, end)`.
    pub fn count(&self, start: usize, end: usize) -> u32 {
        let end = end.min(self.seq.len());
        if start >= end {
            return 0;
        }
        self.prefix(end) - self.prefix(start)
    }

    /// GC fraction in `[start, end)`, over positions that are real bases.
    ///
    /// Masked positions are excluded from the denominator: including them would
    /// report a masked backbone as AT-rich.
    pub fn fraction(&self, start: usize, end: usize) -> Option<f64> {
        let end = end.min(self.seq.len());
        if start >= end {
            return None;
        }
        let denominator = self.seq[start..end]
            .iter()
            .filter(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
            .count();
        if denominator == 0 {
            return None;
        }
        Some(self.count(start, end) as f64 / denominator as f64)
    }

    pub fn substitute(&mut self, offset: usize, replacement: &[u8]) {
        for (i, &byte) in replacement.iter().enumerate() {
            let position = offset + i;
            let before = is_gc(self.seq[position]);
            let after = is_gc(byte);
            if before != after {
                self.add(position, if after { 1 } else { -1 });
            }
            self.seq[position] = byte;
        }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    pub fn sequence(&self) -> &[u8] {
        &self.seq
    }

    /// Every window of `window_nt` whose GC fraction falls outside `[min, max]`.
    pub fn windows_outside(
        &self,
        window_nt: usize,
        min: f64,
        max: f64,
    ) -> Vec<(usize, usize, f64)> {
        let mut out = Vec::new();
        if window_nt == 0 || self.seq.len() < window_nt {
            return out;
        }
        for start in 0..=self.seq.len() - window_nt {
            if let Some(fraction) = self.fraction(start, start + window_nt) {
                if fraction < min || fraction > max {
                    out.push((start, start + window_nt, fraction));
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_runs_are_maximal_and_filtered() {
        let runs = find_runs(b"AAACCCCGT", 3, None);
        assert_eq!(
            runs,
            vec![
                Run {
                    base: b'A',
                    start: 0,
                    end: 3
                },
                Run {
                    base: b'C',
                    start: 3,
                    end: 7
                },
            ]
        );
        assert_eq!(find_runs(b"AAACCCCGT", 4, None).len(), 1);
        assert_eq!(find_runs(b"AAACCCCGT", 3, Some(b'C')).len(), 1);
        // Masked positions break runs rather than extending them.
        assert_eq!(find_runs(b"AANAA", 3, None), vec![]);
    }

    #[test]
    fn unit_run_window_finds_a_run_created_by_an_edit() {
        let mut seq = b"AAACGTAAA".to_vec();
        assert!(find_runs(&seq, 6, None).is_empty());
        seq[3..6].copy_from_slice(b"AAA");
        assert_eq!(find_runs(&seq, 6, None).len(), 1);
        // The windowed scan finds the same run.
        let near = find_runs_near(&seq, 3, 3, 6, None);
        assert_eq!(near, find_runs(&seq, 6, None));
    }

    #[test]
    fn unit_gc_fraction_and_windows() {
        let index = GcIndex::build(b"GGGGAAAA");
        assert_eq!(index.count(0, 8), 4);
        assert_eq!(index.fraction(0, 8), Some(0.5));
        assert_eq!(index.fraction(0, 4), Some(1.0));
        assert_eq!(index.fraction(4, 8), Some(0.0));

        let outside = index.windows_outside(4, 0.25, 0.75);
        assert_eq!(outside.len(), 2, "the all-G and all-A windows are outside");
    }

    #[test]
    fn unit_gc_excludes_masked_positions_from_the_denominator() {
        let index = GcIndex::build(b"GGNN");
        assert_eq!(index.fraction(0, 4), Some(1.0), "masked bases are not AT");
        assert_eq!(GcIndex::build(b"NNNN").fraction(0, 4), None);
    }

    /// TST-CONS-001 (remaining classes): incremental equals a naive rebuild.
    #[test]
    fn property_tst_cons_001_gc_and_runs_incremental_match_rebuild() {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let alphabet = b"ACGT";
        let mut seq: Vec<u8> = (0..400).map(|i| alphabet[(i * 5 + i / 3) % 4]).collect();
        let mut index = GcIndex::build(&seq);

        for _ in 0..5_000 {
            let offset = (next() as usize) % (seq.len() - 3);
            let replacement: Vec<u8> = (0..3).map(|_| alphabet[(next() as usize) % 4]).collect();

            index.substitute(offset, &replacement);
            seq[offset..offset + 3].copy_from_slice(&replacement);

            let rebuilt = GcIndex::build(&seq);
            assert_eq!(index.count(0, seq.len()), rebuilt.count(0, seq.len()));
            // Spot-check windows around the edit, where divergence would appear.
            for start in offset.saturating_sub(40)..(offset + 40).min(seq.len() - 50) {
                assert_eq!(
                    index.fraction(start, start + 50),
                    rebuilt.fraction(start, start + 50),
                    "GC window at {start} diverged after an edit at {offset}"
                );
            }
            // The windowed run scan agrees with a full scan.
            let full = find_runs(&seq, 5, None);
            let near = find_runs_near(&seq, offset, 3, 5, None);
            for run in &near {
                assert!(full.contains(run), "windowed scan invented run {run:?}");
            }
            for run in &full {
                if run.end > offset.saturating_sub(5) && run.start < offset + 3 + 5 {
                    assert!(near.contains(run), "windowed scan missed run {run:?}");
                }
            }
        }
    }
}
