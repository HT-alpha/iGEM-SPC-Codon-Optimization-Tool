//! Aho–Corasick matching over IUPAC patterns, on both strands.
//!
//! Ticket P1-14. Contract: `planning/02-schemas/constraint_set.schema.json`.
//!
//! The alphabet is four symbols, so the goto function is a dense `[usize; 4]`
//! per node and matching is a pointer chase with no branching on the character
//! — which is what makes the incremental rescan cheap enough to run inside the
//! repair loop.
//!
//! IUPAC codes are expanded into the trie rather than into separate patterns:
//! `RYN` becomes six edges at three depths instead of eight literal strings,
//! and shared prefixes stay shared. The expansion is still exponential in the
//! number of ambiguous positions, so the node count is capped and an oversized
//! pattern set is refused with `CF-FEAS-007` rather than exhausting memory.
//!
//! LEGACY PARITY NOTE (P1-27). The prototype checks motifs by scanning a tail
//! window of `len(motif) + 2` characters after appending each codon. Appending
//! three characters can create an occurrence starting at most `len - 1`
//! characters earlier, so that window is exactly right and this implementation
//! uses the same width for its incremental rescan.

use cf_types::error::{code, CfError, Result};

/// Maximum trie nodes before a pattern set is refused.
const MAX_NODES: usize = 1 << 17;

const ROOT: usize = 0;

/// Expand a IUPAC code to the concrete bases it matches.
fn iupac_bases(code: u8) -> Option<&'static [u8]> {
    Some(match code.to_ascii_uppercase() {
        b'A' => &[0],
        b'C' => &[1],
        b'G' => &[2],
        b'T' | b'U' => &[3],
        b'R' => &[0, 2],
        b'Y' => &[1, 3],
        b'S' => &[1, 2],
        b'W' => &[0, 3],
        b'K' => &[2, 3],
        b'M' => &[0, 1],
        b'B' => &[1, 2, 3],
        b'D' => &[0, 2, 3],
        b'H' => &[0, 1, 3],
        b'V' => &[0, 1, 2],
        b'N' => &[0, 1, 2, 3],
        _ => return None,
    })
}

/// Reverse complement of an IUPAC pattern.
fn reverse_complement(pattern: &str) -> Result<String> {
    let mut out = Vec::with_capacity(pattern.len());
    for byte in pattern.bytes().rev() {
        let complement = match byte.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' | b'U' => b'A',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'D' => b'H',
            b'H' => b'D',
            b'V' => b'B',
            b'N' => b'N',
            other => {
                return Err(CfError::new(
                    code::FEAS_CONSTRAINT_INVALID,
                    format!("{:?} is not a IUPAC nucleotide code", other as char),
                ));
            }
        };
        out.push(complement);
    }
    Ok(String::from_utf8(out).expect("ASCII"))
}

/// One pattern registered with the automaton.
#[derive(Clone, Debug)]
pub struct MotifPattern {
    /// Index of the owning constraint, so a hit reports which rule it violated.
    pub constraint: usize,
    pub literal: String,
    pub len: usize,
    /// True when this entry is the reverse complement of the declared pattern.
    pub reverse: bool,
}

/// A hit: the pattern index and the half-open nucleotide span it occupies.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MotifHit {
    pub pattern: usize,
    pub start: usize,
    pub end: usize,
}

#[derive(Debug, Default)]
struct Node {
    goto: [usize; 4],
    fail: usize,
    /// Patterns ending at this node, plus those reachable via fail links.
    outputs: Vec<usize>,
}

/// Multi-pattern automaton over `{A, C, G, T}`.
#[derive(Debug)]
pub struct MotifAutomaton {
    nodes: Vec<Node>,
    patterns: Vec<MotifPattern>,
    max_len: usize,
}

impl MotifAutomaton {
    pub fn empty() -> Self {
        Self {
            nodes: vec![Node::default()],
            patterns: Vec::new(),
            max_len: 0,
        }
    }

    /// Build from `(constraint index, pattern, strand)` triples.
    pub fn build(specs: &[(usize, String, bool)]) -> Result<Self> {
        let mut automaton = Self::empty();

        for (constraint, pattern, both_strands) in specs {
            automaton.add(*constraint, pattern, false)?;
            if *both_strands {
                let rc = reverse_complement(pattern)?;
                // A palindromic site matches itself; registering it twice would
                // double-count every hit.
                if !rc.eq_ignore_ascii_case(pattern) {
                    automaton.add(*constraint, &rc, true)?;
                }
            }
        }
        automaton.link();
        Ok(automaton)
    }

    fn add(&mut self, constraint: usize, pattern: &str, reverse: bool) -> Result<()> {
        if pattern.is_empty() {
            return Err(CfError::new(
                code::FEAS_CONSTRAINT_INVALID,
                "an empty motif matches everywhere",
            ));
        }
        let index = self.patterns.len();
        self.patterns.push(MotifPattern {
            constraint,
            literal: pattern.to_ascii_uppercase(),
            len: pattern.len(),
            reverse,
        });
        self.max_len = self.max_len.max(pattern.len());

        // Breadth-first over the ambiguity expansion, sharing prefixes.
        let mut frontier = vec![ROOT];
        for byte in pattern.bytes() {
            let bases = iupac_bases(byte).ok_or_else(|| {
                CfError::new(
                    code::FEAS_CONSTRAINT_INVALID,
                    format!("{:?} is not a IUPAC nucleotide code", byte as char),
                )
            })?;
            let mut next = Vec::with_capacity(frontier.len() * bases.len());
            for &node in &frontier {
                for &base in bases {
                    let mut child = self.nodes[node].goto[base as usize];
                    if child == 0 {
                        if self.nodes.len() >= MAX_NODES {
                            return Err(CfError::new(
                                code::FEAS_REGEX_TOO_LARGE,
                                format!(
                                    "motif set exceeds the {MAX_NODES}-node cap; \
                                     reduce the number of ambiguous positions"
                                ),
                            ));
                        }
                        child = self.nodes.len();
                        self.nodes.push(Node::default());
                        self.nodes[node].goto[base as usize] = child;
                    }
                    if !next.contains(&child) {
                        next.push(child);
                    }
                }
            }
            frontier = next;
        }
        for node in frontier {
            self.nodes[node].outputs.push(index);
        }
        Ok(())
    }

    /// Build fail links and propagate outputs (standard Aho–Corasick).
    fn link(&mut self) {
        let mut queue = std::collections::VecDeque::new();
        for base in 0..4 {
            let child = self.nodes[ROOT].goto[base];
            if child != 0 {
                self.nodes[child].fail = ROOT;
                queue.push_back(child);
            }
        }
        while let Some(node) = queue.pop_front() {
            let inherited = self.nodes[self.nodes[node].fail].outputs.clone();
            for pattern in inherited {
                if !self.nodes[node].outputs.contains(&pattern) {
                    self.nodes[node].outputs.push(pattern);
                }
            }
            for base in 0..4 {
                let child = self.nodes[node].goto[base];
                let fail = self.nodes[node].fail;
                if child == 0 {
                    self.nodes[node].goto[base] = self.nodes[fail].goto[base];
                } else {
                    self.nodes[child].fail = self.nodes[fail].goto[base];
                    queue.push_back(child);
                }
            }
        }
    }

    pub fn patterns(&self) -> &[MotifPattern] {
        &self.patterns
    }

    pub fn is_empty(&self) -> bool {
        self.patterns.is_empty()
    }

    /// Width of the window that must be rescanned after a 3 nt substitution.
    ///
    /// A replaced codon can complete an occurrence that began up to
    /// `max_len - 1` characters earlier and can end up to `max_len - 1`
    /// characters later, so the window spans `max_len - 1` on each side of the
    /// three changed characters.
    pub fn rescan_radius(&self) -> usize {
        self.max_len.saturating_sub(1)
    }

    /// Every occurrence in `haystack`, reported as absolute spans.
    ///
    /// `offset` shifts the reported spans, so a windowed rescan can report
    /// coordinates in the full construct.
    pub fn find_all(&self, haystack: &[u8], offset: usize, out: &mut Vec<MotifHit>) {
        if self.patterns.is_empty() {
            return;
        }
        let mut state = ROOT;
        for (i, &byte) in haystack.iter().enumerate() {
            let Some(base) = base_index(byte) else {
                // An unknown character (N in a masked backbone) cannot be part
                // of a match; reset rather than pretending it is a base.
                state = ROOT;
                continue;
            };
            state = self.nodes[state].goto[base as usize];
            for &pattern in &self.nodes[state].outputs {
                let len = self.patterns[pattern].len;
                let end = i + 1;
                if end >= len {
                    out.push(MotifHit {
                        pattern,
                        start: offset + end - len,
                        end: offset + end,
                    });
                }
            }
        }
    }
}

#[inline]
fn base_index(byte: u8) -> Option<u8> {
    match byte {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hits(automaton: &MotifAutomaton, seq: &str) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        automaton.find_all(seq.as_bytes(), 0, &mut out);
        let mut spans: Vec<(usize, usize)> = out.iter().map(|h| (h.start, h.end)).collect();
        spans.sort();
        spans
    }

    #[test]
    fn unit_literal_and_overlapping_matches() {
        let a = MotifAutomaton::build(&[(0, "AAA".into(), false)]).unwrap();
        // Overlapping occurrences are all reported; a run of 5 A's has 3.
        assert_eq!(hits(&a, "AAAAA"), vec![(0, 3), (1, 4), (2, 5)]);
        assert_eq!(hits(&a, "AACAA"), vec![]);
    }

    #[test]
    fn unit_both_strands() {
        // BsaI: GGTCTC, reverse complement GAGACC.
        let a = MotifAutomaton::build(&[(0, "GGTCTC".into(), true)]).unwrap();
        assert_eq!(a.patterns().len(), 2);
        assert_eq!(hits(&a, "AAGGTCTCAA").len(), 1);
        assert_eq!(hits(&a, "AAGAGACCAA").len(), 1);
        assert_eq!(hits(&a, "AAGGTCTAAA").len(), 0);
    }

    #[test]
    fn unit_palindrome_is_not_double_counted() {
        // EcoRI GAATTC is its own reverse complement.
        let a = MotifAutomaton::build(&[(0, "GAATTC".into(), true)]).unwrap();
        assert_eq!(a.patterns().len(), 1, "a palindrome must register once");
        assert_eq!(hits(&a, "TTGAATTCTT").len(), 1);
    }

    #[test]
    fn unit_iupac_expansion() {
        // GGTCTN matches GGTCTA/C/G/T.
        let a = MotifAutomaton::build(&[(0, "GGTCTN".into(), false)]).unwrap();
        for last in ["A", "C", "G", "T"] {
            assert_eq!(
                hits(&a, &format!("AA GGTCT{last} AA").replace(' ', "")).len(),
                1
            );
        }
        // Shared prefixes mean the trie is far smaller than 4 literals.
        let ry = MotifAutomaton::build(&[(0, "RYRY".into(), false)]).unwrap();
        assert_eq!(hits(&ry, "ACAC").len(), 1);
        assert_eq!(hits(&ry, "ACAA").len(), 0);
    }

    #[test]
    fn unit_multiple_patterns_share_the_automaton() {
        let a = MotifAutomaton::build(&[
            (0, "GGTCTC".into(), true),
            (1, "AAGCTT".into(), true),
            (2, "GCGGCCGC".into(), true),
        ])
        .unwrap();
        let mut out = Vec::new();
        a.find_all(b"GGTCTCAAGCTTGCGGCCGC", 0, &mut out);
        let mut owners: Vec<usize> = out
            .iter()
            .map(|h| a.patterns()[h.pattern].constraint)
            .collect();
        owners.sort();
        assert_eq!(owners, vec![0, 1, 2]);
    }

    #[test]
    fn unit_offset_is_applied_to_reported_spans() {
        let a = MotifAutomaton::build(&[(0, "AAA".into(), false)]).unwrap();
        let mut out = Vec::new();
        a.find_all(b"GAAAG", 100, &mut out);
        assert_eq!(
            out,
            vec![MotifHit {
                pattern: 0,
                start: 101,
                end: 104
            }]
        );
    }

    #[test]
    fn unit_non_bases_break_matches_rather_than_matching() {
        // A masked backbone uses N; a motif must not match across the mask.
        let a = MotifAutomaton::build(&[(0, "AAAA".into(), false)]).unwrap();
        assert_eq!(hits(&a, "AANAA"), vec![]);
        assert_eq!(hits(&a, "AAAA"), vec![(0, 4)]);
    }

    #[test]
    fn unit_rescan_radius_covers_a_codon_substitution() {
        let a = MotifAutomaton::build(&[(0, "GCGGCCGC".into(), false)]).unwrap();
        assert_eq!(a.rescan_radius(), 7);
    }

    #[test]
    fn unit_bad_input_is_refused() {
        assert_eq!(
            MotifAutomaton::build(&[(0, "GGZTCT".into(), false)])
                .unwrap_err()
                .code,
            code::FEAS_CONSTRAINT_INVALID
        );
        assert!(MotifAutomaton::build(&[(0, String::new(), false)]).is_err());
        // 16 fully-ambiguous positions would need 4^16 leaves.
        assert_eq!(
            MotifAutomaton::build(&[(0, "N".repeat(16), false)])
                .unwrap_err()
                .code,
            code::FEAS_REGEX_TOO_LARGE
        );
    }
}
