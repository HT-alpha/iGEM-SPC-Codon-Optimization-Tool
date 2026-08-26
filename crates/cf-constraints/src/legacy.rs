//! The `codon_6.0.py` feasibility predicate, transcribed.
//!
//! Ticket P1-27. Sequencing: this depends on `P1-14`–`P1-16` (the detectors),
//! not on `P1-17` (the repair operator) as originally declared — its acceptance
//! criterion exercises only the checker.
//!
//! WHY THIS LIVES IN THE CRATE RATHER THAN IN A TEST. The legacy semantics are
//! now a contract: the k-mer detector was built to preserve them, and the
//! differential test in `python/tests/test_p1_27_legacy_differential.py`
//! compares against the real `codon_6.0.py`. Keeping the transcription next to
//! the detector means there is one definition of what "the legacy behaviour"
//! is, and that a change to either has to confront the other.
//!
//! THE PREDICATE. `optimize_cds_sequence` accepts a codon extension when, and
//! only when:
//!
//!   1. no forbidden motif occurs in the coding sequence built so far, and
//!   2. for every k-mer window, `backbone_count[kmer] + count_overlapping(cds,
//!      kmer) <= max_allowed_repeats`.
//!
//! Its incremental checks are equivalent to evaluating those two conditions on
//! the finished sequence. For (1) a motif wholly inside an earlier prefix was
//! rejected when that prefix was built. For (2) the count a prefix sees is
//! monotone in the prefix, so the tightest check happens when the last
//! occurrence of a k-mer is created, at which point the prefix count equals the
//! whole-sequence count. So the whole-sequence predicate below is exactly what
//! the prototype enforces, and it is what the differential test compares.
//!
//! TWO LEGACY BLIND SPOTS, PRESERVED HERE AND DOCUMENTED IN
//! `docs/legacy-parity.md`:
//!
//!   * motifs are matched against the coding sequence only, never the
//!     backbone, so a restriction site straddling the CDS boundary is invisible;
//!   * backbone k-mer counts come from a construct whose coding regions are
//!     masked with `N`, and masked k-mers are skipped, so a repeat straddling
//!     the boundary is counted by neither side.
//!
//! Both are reproduced faithfully. Neither is the default in
//! [`crate::ConstraintChecker`].

use std::collections::{BTreeSet, HashMap};

use crate::kmer::count_overlapping;

/// Inputs to the legacy predicate, mirroring `optimize_cds_sequence`'s
/// arguments.
#[derive(Clone, Debug)]
pub struct LegacyScenario {
    /// Forbidden motifs, as literal uppercase strings. The prototype has no
    /// IUPAC handling and no strand awareness: it matches literals, and gets
    /// both strands only because its restriction-site list spells out each
    /// reverse complement.
    pub motifs: BTreeSet<String>,
    pub k: usize,
    pub max_allowed_repeats: u32,
    /// k-mer counts over the construct with every coding region masked out.
    pub backbone_counts: HashMap<String, u32>,
}

impl LegacyScenario {
    /// Build the scenario the way `optimize_whole_construct` does.
    pub fn from_construct(
        construct: &str,
        cds_regions: &[&str],
        k: usize,
        max_allowed_repeats: u32,
        short_motif_repeat_freq: u32,
        homopolymer_len: usize,
    ) -> Self {
        let clean: String = construct
            .split_whitespace()
            .collect::<String>()
            .to_uppercase();
        let motifs = auto_detect_forbidden_motifs(&clean, short_motif_repeat_freq, homopolymer_len);

        let mut backbone = clean.clone();
        for cds in cds_regions {
            let clean_cds: String = cds.split_whitespace().collect::<String>().to_uppercase();
            backbone = backbone.replace(&clean_cds, &"N".repeat(clean_cds.len()));
        }
        Self {
            motifs,
            k,
            max_allowed_repeats,
            backbone_counts: build_kmer_frequency_counts(&backbone, k),
        }
    }
}

/// `build_kmer_frequency_counts`: overlapping counts, skipping masked k-mers.
pub fn build_kmer_frequency_counts(sequence: &str, k: usize) -> HashMap<String, u32> {
    let bytes = sequence.as_bytes();
    let mut counts: HashMap<String, u32> = HashMap::new();
    if bytes.len() < k {
        return counts;
    }
    for start in 0..=bytes.len() - k {
        let window = &bytes[start..start + k];
        if window.contains(&b'N') {
            continue;
        }
        *counts
            .entry(String::from_utf8(window.to_vec()).expect("ASCII"))
            .or_insert(0) += 1;
    }
    counts
}

/// `auto_detect_forbidden_motifs`, including its data-dependent behaviour: a
/// motif is forbidden only if it is already present in the input construct.
pub fn auto_detect_forbidden_motifs(
    construct: &str,
    repeat_threshold: u32,
    homopolymer_len: usize,
) -> BTreeSet<String> {
    let seq: String = construct
        .split_whitespace()
        .collect::<String>()
        .to_uppercase();
    let mut motifs = BTreeSet::new();

    // The prototype's literal list. Both strands of BsaI and BsmBI are spelled
    // out; EcoRI, BamHI and HindIII are palindromic; NotI is palindromic.
    const RESTRICTION_SITES: &[&str] = &[
        "GGTCTC", "GAGACC", // BsaI
        "CGTCTC", "GAGACG", // BsmBI
        "GAATTC", "GGATCC", "AAGCTT", "GCGGCCGC", // EcoRI, BamHI, HindIII, NotI
    ];
    for site in RESTRICTION_SITES {
        if seq.contains(site) {
            motifs.insert((*site).to_string());
        }
    }

    for base in ["A", "C", "G", "T"] {
        let run = base.repeat(homopolymer_len);
        if seq.contains(&run) {
            motifs.insert(run);
        }
    }

    for k in [5usize, 6] {
        let counts = {
            let bytes = seq.as_bytes();
            let mut counts: HashMap<&[u8], u32> = HashMap::new();
            if bytes.len() >= k {
                for start in 0..=bytes.len() - k {
                    *counts.entry(&bytes[start..start + k]).or_insert(0) += 1;
                }
            }
            counts
        };
        for (kmer, count) in counts {
            if count >= repeat_threshold {
                let gc = kmer.iter().filter(|b| **b == b'G' || **b == b'C').count();
                if gc as f64 / k as f64 >= 0.8 {
                    motifs.insert(String::from_utf8(kmer.to_vec()).expect("ASCII"));
                }
            }
        }
    }
    motifs
}

/// Why the legacy predicate rejected a coding sequence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LegacyRejection {
    Motif {
        motif: String,
        at: usize,
    },
    KmerOverLimit {
        kmer: String,
        count: u32,
        limit: u32,
    },
}

/// The legacy feasibility predicate for a complete coding sequence.
pub fn legacy_violations(cds: &str, scenario: &LegacyScenario) -> Vec<LegacyRejection> {
    let cds = cds.to_uppercase();
    let bytes = cds.as_bytes();
    let mut out = Vec::new();

    for motif in &scenario.motifs {
        if motif.len() > cds.len() {
            continue;
        }
        let needle = motif.as_bytes();
        for start in 0..=bytes.len() - needle.len() {
            if &bytes[start..start + needle.len()] == needle {
                out.push(LegacyRejection::Motif {
                    motif: motif.clone(),
                    at: start,
                });
            }
        }
    }

    if bytes.len() >= scenario.k {
        let mut reported = BTreeSet::new();
        for start in 0..=bytes.len() - scenario.k {
            let kmer = &cds[start..start + scenario.k];
            let backbone = *scenario.backbone_counts.get(kmer).unwrap_or(&0);
            let within = count_overlapping(bytes, kmer.as_bytes());
            let total = backbone + within;
            if total > scenario.max_allowed_repeats && reported.insert(kmer.to_string()) {
                out.push(LegacyRejection::KmerOverLimit {
                    kmer: kmer.to_string(),
                    count: total,
                    limit: scenario.max_allowed_repeats,
                });
            }
        }
    }
    out
}

/// Whether the prototype would accept this coding sequence.
pub fn legacy_is_feasible(cds: &str, scenario: &LegacyScenario) -> bool {
    legacy_violations(cds, scenario).is_empty()
}

#[cfg(test)]
mod tests {
    use super::*;

    const INSERT: &str = "AAAGAGGAGAAAGATCTATGAAACGTATTCTGTAAGGATCCAAAAAAACCAGGCATCAAATAAAACGAAAGG";
    const CDS: &str = "ATGAAACGTATTCTGTAA";

    fn scenario() -> LegacyScenario {
        LegacyScenario::from_construct(INSERT, &[CDS], 8, 2, 6, 7)
    }

    #[test]
    fn unit_auto_detection_is_data_dependent() {
        // BamHI (GGATCC) is present in the insert, so it becomes forbidden;
        // BsaI is not present, so it does not.
        let motifs = auto_detect_forbidden_motifs(INSERT, 6, 7);
        assert!(motifs.contains("GGATCC"));
        assert!(!motifs.contains("GGTCTC"));
        // A run of seven A's is present.
        assert!(motifs.contains("AAAAAAA"));
    }

    #[test]
    fn unit_backbone_masking_skips_the_coding_region() {
        let scenario = scenario();
        // No backbone k-mer may contain N, so none may overlap the CDS.
        for kmer in scenario.backbone_counts.keys() {
            assert!(!kmer.contains('N'));
        }
    }

    #[test]
    fn unit_predicate_detects_a_planted_motif() {
        let scenario = scenario();
        assert!(legacy_is_feasible(CDS, &scenario));

        // Plant BamHI inside the coding sequence.
        let planted = "ATGAAACGTGGATCCCTGTAA";
        let violations = legacy_violations(planted, &scenario);
        assert!(violations
            .iter()
            .any(|v| matches!(v, LegacyRejection::Motif { motif, .. } if motif == "GGATCC")));
    }

    #[test]
    fn unit_predicate_counts_overlapping_kmers_against_the_backbone() {
        let mut scenario = scenario();
        scenario.k = 4;
        scenario.max_allowed_repeats = 1;
        scenario.backbone_counts.clear();
        scenario.backbone_counts.insert("ATGA".to_string(), 1);
        // One occurrence in the CDS plus one in the backbone exceeds a limit of 1.
        let violations = legacy_violations("ATGAAACGT", &scenario);
        assert!(violations.iter().any(
            |v| matches!(v, LegacyRejection::KmerOverLimit { kmer, count, .. } if kmer == "ATGA" && *count == 2)
        ));
    }
}
