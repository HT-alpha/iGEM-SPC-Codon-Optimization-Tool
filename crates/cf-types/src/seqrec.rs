//! The canonical codon-indexed, frame-safe sequence representation.
//!
//! Ticket P1-06. Contract: `planning/02-schemas/sequence.schema.json`.
//!
//! FRAME SAFETY IS STRUCTURAL, NOT VALIDATED. The coding region is a
//! `Vec<Codon>` — a vector of validated codon indices — and untranslated
//! context lives in separate flank fields that no operation edits. There is
//! therefore no representable state in which an edit shifts the reading frame:
//! the only mutating operation is [`CodonSequence::substitute`], which replaces
//! one codon with another.
//!
//! PROTEIN INVARIANCE IS A PRECONDITION, NOT A TEST. `substitute` rejects a
//! non-synonymous replacement with `CF-SEQ-004`. That makes "every search
//! operator preserves the protein" a property of the type rather than something
//! each operator has to remember, which is what invariant I3 requires.

use crate::canonical::Digest;
use crate::codon::{base_index, Codon};
use crate::error::{code, CfError, Result};
use crate::tables::GeneticCode;

/// Why a span may not be edited.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum LockReason {
    StartCodon,
    StopCodon,
    UserLocked,
    OverlappingOrf,
    AffinityTag,
    CleavageSite,
    InternalRbs,
    SignalPeptide,
    ProgrammedFrameshift,
    Other,
}

impl LockReason {
    pub fn as_str(self) -> &'static str {
        match self {
            LockReason::StartCodon => "start_codon",
            LockReason::StopCodon => "stop_codon",
            LockReason::UserLocked => "user_locked",
            LockReason::OverlappingOrf => "overlapping_orf",
            LockReason::AffinityTag => "affinity_tag",
            LockReason::CleavageSite => "cleavage_site",
            LockReason::InternalRbs => "internal_rbs",
            LockReason::SignalPeptide => "signal_peptide",
            LockReason::ProgrammedFrameshift => "programmed_frameshift",
            LockReason::Other => "other",
        }
    }

    pub fn parse(text: &str) -> Option<Self> {
        Some(match text {
            "start_codon" => LockReason::StartCodon,
            "stop_codon" => LockReason::StopCodon,
            "user_locked" => LockReason::UserLocked,
            "overlapping_orf" => LockReason::OverlappingOrf,
            "affinity_tag" => LockReason::AffinityTag,
            "cleavage_site" => LockReason::CleavageSite,
            "internal_rbs" => LockReason::InternalRbs,
            "signal_peptide" => LockReason::SignalPeptide,
            "programmed_frameshift" => LockReason::ProgrammedFrameshift,
            "other" => LockReason::Other,
            _ => return None,
        })
    }
}

/// A half-open codon span `[start, end)` that no operator may alter.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct LockedSpan {
    pub start: usize,
    pub end: usize,
    pub reason: LockReason,
    pub note: Option<String>,
}

/// The canonical sequence type.
#[derive(Clone)]
pub struct CodonSequence {
    codons: Vec<Codon>,
    table: &'static GeneticCode,
    /// Cached translation. Invariant under every exposed operation, which is
    /// why it can be cached at all.
    protein: Vec<u8>,
    includes_stop: bool,
    five_prime_flank: Vec<u8>,
    three_prime_flank: Vec<u8>,
    locked: Vec<LockedSpan>,
    /// Per-codon lock lookup, derived from `locked`. `substitute` is on the
    /// per-move path, so the check must be O(1) rather than a span scan.
    locked_mask: Vec<bool>,
    sequence_id: Option<String>,
    source: Option<serde_json::Value>,
    /// Lazily computed and invalidated on mutation: recomputing a digest on
    /// every substitution would make the hot path O(L) per move.
    digest_cache: Option<Digest>,
}

impl CodonSequence {
    pub const SCHEMA_VERSION: &'static str = "1.0.0";

    /// Parse a coding sequence.
    ///
    /// Applies the default locks: the start codon, and the terminal stop codon
    /// when one is present. Locking the start codon is not cosmetic — in tables
    /// where `GTG` and `TTG` are start codons they encode valine and leucine,
    /// so a synonymous substitution there would silently destroy initiation
    /// while preserving the protein.
    pub fn from_cds(cds: &str, translation_table: u8) -> Result<Self> {
        let table = GeneticCode::get(translation_table)?;
        let bytes = cds.as_bytes();

        if bytes.is_empty() {
            return Err(CfError::new(code::SEQ_EMPTY, "coding sequence is empty"));
        }
        if bytes.len() % 3 != 0 {
            return Err(CfError::new(
                code::SEQ_NOT_IN_FRAME,
                format!("length {} is not a multiple of 3", bytes.len()),
            ));
        }

        let mut codons = Vec::with_capacity(bytes.len() / 3);
        for chunk in bytes.chunks(3) {
            codons.push(Codon::from_bytes(chunk)?);
        }

        // Stop-codon placement: a stop is allowed only as the final codon.
        let last = *codons.last().expect("non-empty");
        let includes_stop = table.is_stop(last);
        let internal_stop_limit = if includes_stop {
            codons.len() - 1
        } else {
            codons.len()
        };
        for (i, codon) in codons[..internal_stop_limit].iter().enumerate() {
            if table.is_stop(*codon) {
                return Err(CfError::new(
                    code::SEQ_INTERNAL_STOP,
                    format!(
                        "codon {i} ({codon}) is a stop in translation table {translation_table}"
                    ),
                ));
            }
        }

        let protein = table.translate(&codons);
        let mut seq = Self {
            codons,
            table,
            protein,
            includes_stop,
            five_prime_flank: Vec::new(),
            three_prime_flank: Vec::new(),
            locked: Vec::new(),
            locked_mask: Vec::new(),
            sequence_id: None,
            source: None,
            digest_cache: None,
        };
        seq.rebuild_locked_mask();
        seq.apply_default_locks();
        Ok(seq)
    }

    fn apply_default_locks(&mut self) {
        let mut locks = vec![LockedSpan {
            start: 0,
            end: 1,
            reason: LockReason::StartCodon,
            note: None,
        }];
        if self.includes_stop {
            locks.push(LockedSpan {
                start: self.codons.len() - 1,
                end: self.codons.len(),
                reason: LockReason::StopCodon,
                note: None,
            });
        }
        self.locked = locks;
        self.rebuild_locked_mask();
    }

    fn rebuild_locked_mask(&mut self) {
        self.locked_mask = vec![false; self.codons.len()];
        for span in &self.locked {
            for i in span.start..span.end.min(self.codons.len()) {
                self.locked_mask[i] = true;
            }
        }
    }

    // ---- accessors -------------------------------------------------------

    pub fn codons(&self) -> &[Codon] {
        &self.codons
    }

    pub fn codon_count(&self) -> usize {
        self.codons.len()
    }

    pub fn table(&self) -> &'static GeneticCode {
        self.table
    }

    pub fn includes_stop(&self) -> bool {
        self.includes_stop
    }

    pub fn protein(&self) -> &[u8] {
        &self.protein
    }

    pub fn protein_string(&self) -> String {
        String::from_utf8(self.protein.clone()).expect("ASCII residues")
    }

    pub fn cds_string(&self) -> String {
        let mut out = String::with_capacity(self.codons.len() * 3);
        for codon in &self.codons {
            out.push_str(std::str::from_utf8(&codon.to_bytes()).expect("ASCII"));
        }
        out
    }

    pub fn five_prime_flank(&self) -> &[u8] {
        &self.five_prime_flank
    }

    pub fn three_prime_flank(&self) -> &[u8] {
        &self.three_prime_flank
    }

    pub fn locked_spans(&self) -> &[LockedSpan] {
        &self.locked
    }

    #[inline]
    pub fn is_locked(&self, index: usize) -> bool {
        self.locked_mask.get(index).copied().unwrap_or(true)
    }

    /// Codon indices an operator is allowed to change.
    pub fn editable_positions(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.codons.len()).filter(move |i| !self.locked_mask[*i])
    }

    /// Synonymous alternatives at a position, including the current codon.
    #[inline]
    pub fn synonyms_at(&self, index: usize) -> &[Codon] {
        self.table.synonyms(self.codons[index])
    }

    pub fn sequence_id(&self) -> Option<&str> {
        self.sequence_id.as_deref()
    }

    // ---- builders --------------------------------------------------------

    pub fn with_five_prime_flank(mut self, flank: &str) -> Result<Self> {
        self.five_prime_flank = normalize_flank(flank)?;
        self.digest_cache = None;
        Ok(self)
    }

    pub fn with_three_prime_flank(mut self, flank: &str) -> Result<Self> {
        self.three_prime_flank = normalize_flank(flank)?;
        self.digest_cache = None;
        Ok(self)
    }

    pub fn with_sequence_id(mut self, id: impl Into<String>) -> Self {
        self.sequence_id = Some(id.into());
        self.digest_cache = None;
        self
    }

    pub fn with_source(mut self, source: serde_json::Value) -> Self {
        self.source = Some(source);
        self.digest_cache = None;
        self
    }

    /// Add a lock, validating invariant I4.
    pub fn add_lock(&mut self, start: usize, end: usize, reason: LockReason) -> Result<()> {
        if start >= end || end > self.codons.len() {
            return Err(CfError::new(
                code::SEQ_LOCK_OUT_OF_RANGE,
                format!(
                    "locked span [{start}, {end}) is not within [0, {})",
                    self.codons.len()
                ),
            ));
        }
        self.locked.push(LockedSpan {
            start,
            end,
            reason,
            note: None,
        });
        self.rebuild_locked_mask();
        self.digest_cache = None;
        Ok(())
    }

    /// Remove every lock, including the defaults.
    ///
    /// Exposed for the synthetic test suite, which needs a fully editable
    /// sequence to enumerate a design space exhaustively.
    pub fn clear_locks(&mut self) {
        self.locked.clear();
        self.rebuild_locked_mask();
        self.digest_cache = None;
    }

    // ---- the one mutating operation --------------------------------------

    /// Replace the codon at `index`, returning the codon that was there.
    ///
    /// Rejects, in this order: an out-of-range index (`CF-SEQ-005`), a locked
    /// position (`CF-SEQ-006`), and a non-synonymous replacement
    /// (`CF-SEQ-004`). The third is what makes protein invariance structural.
    #[inline]
    pub fn substitute(&mut self, index: usize, new: Codon) -> Result<Codon> {
        if index >= self.codons.len() {
            return Err(CfError::new(
                code::SEQ_INDEX_OUT_OF_RANGE,
                format!(
                    "codon index {index} is outside a sequence of {} codons",
                    self.codons.len()
                ),
            ));
        }
        if self.locked_mask[index] {
            return Err(CfError::new(
                code::SEQ_POSITION_LOCKED,
                format!("codon {index} is locked and may not be edited"),
            ));
        }
        let old = self.codons[index];
        if old == new {
            return Ok(old);
        }
        if !self.table.are_synonymous(old, new) {
            return Err(CfError::new(
                code::SEQ_NOT_SYNONYMOUS,
                format!(
                    "{old} encodes {} but {new} encodes {} in translation table {}",
                    self.table.amino_acid(old) as char,
                    self.table.amino_acid(new) as char,
                    self.table.id()
                ),
            ));
        }
        self.codons[index] = new;
        self.digest_cache = None;
        Ok(old)
    }

    /// Non-mutating form, for callers that need the original.
    pub fn substituted(&self, index: usize, new: Codon) -> Result<Self> {
        let mut next = self.clone();
        next.substitute(index, new)?;
        Ok(next)
    }

    // ---- serialization ---------------------------------------------------

    /// Serialize to a document conforming to `sequence.schema.json`.
    pub fn to_json(&self) -> Result<serde_json::Value> {
        let mut map = serde_json::Map::new();
        map.insert("schema_version".into(), Self::SCHEMA_VERSION.into());
        if let Some(id) = &self.sequence_id {
            map.insert("sequence_id".into(), id.clone().into());
        }
        map.insert("alphabet".into(), "DNA".into());
        map.insert("translation_table".into(), self.table.id().into());
        map.insert("cds".into(), self.cds_string().into());
        map.insert("codon_count".into(), self.codons.len().into());
        map.insert("protein".into(), self.protein_string().into());
        map.insert("includes_stop".into(), self.includes_stop.into());
        map.insert(
            "five_prime_flank".into(),
            String::from_utf8(self.five_prime_flank.clone())
                .expect("ASCII")
                .into(),
        );
        map.insert(
            "three_prime_flank".into(),
            String::from_utf8(self.three_prime_flank.clone())
                .expect("ASCII")
                .into(),
        );
        let locks: Vec<serde_json::Value> = self
            .locked
            .iter()
            .map(|span| {
                let mut entry = serde_json::Map::new();
                entry.insert(
                    "span".into(),
                    serde_json::json!({ "start": span.start, "end": span.end }),
                );
                entry.insert("reason".into(), span.reason.as_str().into());
                if let Some(note) = &span.note {
                    entry.insert("note".into(), note.clone().into());
                }
                serde_json::Value::Object(entry)
            })
            .collect();
        map.insert("locked_codons".into(), locks.into());
        if let Some(source) = &self.source {
            map.insert("source".into(), source.clone());
        }

        let partial = serde_json::Value::Object(map.clone());
        let digest = Digest::of_json(&partial)?;
        map.insert("digest".into(), digest.to_string().into());
        Ok(serde_json::Value::Object(map))
    }

    /// The document digest (invariant I5).
    pub fn digest(&self) -> Result<Digest> {
        if let Some(cached) = &self.digest_cache {
            return Ok(cached.clone());
        }
        let doc = self.to_json()?;
        Digest::parse(doc["digest"].as_str().expect("to_json writes a digest"))
    }

    /// Parse a document, re-checking every invariant rather than trusting it.
    pub fn from_json(value: &serde_json::Value) -> Result<Self> {
        let cds = value["cds"]
            .as_str()
            .ok_or_else(|| CfError::new(code::SEQ_EMPTY, "document has no `cds` string"))?;
        let table_id = value["translation_table"].as_u64().ok_or_else(|| {
            CfError::new(
                code::SEQ_UNKNOWN_TABLE,
                "document has no `translation_table`",
            )
        })? as u8;

        let mut seq = Self::from_cds(cds, table_id)?;

        // I1: declared codon count matches.
        if let Some(declared) = value["codon_count"].as_u64() {
            if declared as usize != seq.codons.len() {
                return Err(CfError::new(
                    code::SEQ_NOT_IN_FRAME,
                    format!(
                        "declared codon_count {declared} but the cds has {} codons",
                        seq.codons.len()
                    ),
                ));
            }
        }

        // I2: declared protein matches the translation.
        if let Some(declared) = value["protein"].as_str() {
            if declared != seq.protein_string() {
                return Err(CfError::new(
                    code::SEQ_PROTEIN_MISMATCH,
                    format!(
                        "declared protein {declared:?} but the cds translates to {:?}",
                        seq.protein_string()
                    ),
                ));
            }
        }

        if let Some(flag) = value["includes_stop"].as_bool() {
            if flag != seq.includes_stop {
                return Err(CfError::new(
                    code::SEQ_STOP_PLACEMENT,
                    format!(
                        "declared includes_stop={flag} but the terminal codon {} a stop",
                        if seq.includes_stop { "is" } else { "is not" }
                    ),
                ));
            }
        }

        if let Some(flank) = value["five_prime_flank"].as_str() {
            seq.five_prime_flank = normalize_flank(flank)?;
        }
        if let Some(flank) = value["three_prime_flank"].as_str() {
            seq.three_prime_flank = normalize_flank(flank)?;
        }
        if let Some(id) = value["sequence_id"].as_str() {
            seq.sequence_id = Some(id.to_string());
        }
        if let Some(source) = value.get("source") {
            if !source.is_null() {
                seq.source = Some(source.clone());
            }
        }

        // I4: locked spans are in range.
        if let Some(locks) = value["locked_codons"].as_array() {
            let mut parsed = Vec::with_capacity(locks.len());
            for lock in locks {
                let start = lock["span"]["start"].as_u64().ok_or_else(|| {
                    CfError::new(code::SEQ_LOCK_OUT_OF_RANGE, "lock has no span.start")
                })? as usize;
                let end = lock["span"]["end"].as_u64().ok_or_else(|| {
                    CfError::new(code::SEQ_LOCK_OUT_OF_RANGE, "lock has no span.end")
                })? as usize;
                if start >= end || end > seq.codons.len() {
                    return Err(CfError::new(
                        code::SEQ_LOCK_OUT_OF_RANGE,
                        format!(
                            "locked span [{start}, {end}) is not within [0, {})",
                            seq.codons.len()
                        ),
                    ));
                }
                let reason = lock["reason"]
                    .as_str()
                    .and_then(LockReason::parse)
                    .unwrap_or(LockReason::Other);
                parsed.push(LockedSpan {
                    start,
                    end,
                    reason,
                    note: lock["note"].as_str().map(str::to_string),
                });
            }
            seq.locked = parsed;
            seq.rebuild_locked_mask();
        }

        // I5: declared digest matches the recomputed one.
        if let Some(declared) = value["digest"].as_str() {
            let recomputed = seq.digest()?;
            if declared != recomputed.as_str() {
                return Err(CfError::new(
                    code::SEQ_DIGEST_MISMATCH,
                    format!("declared {declared} but recomputed {recomputed}"),
                ));
            }
        }

        Ok(seq)
    }
}

fn normalize_flank(flank: &str) -> Result<Vec<u8>> {
    let mut out = Vec::with_capacity(flank.len());
    for byte in flank.bytes() {
        let index = base_index(byte).ok_or_else(|| {
            CfError::new(
                code::SEQ_INVALID_BASE,
                format!(
                    "flank contains {:?}, which is not A, C, G, T or U",
                    byte as char
                ),
            )
        })?;
        out.push(crate::codon::BASES[index as usize]);
    }
    Ok(out)
}

impl std::fmt::Debug for CodonSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CodonSequence({} codons, table {}, protein {})",
            self.codons.len(),
            self.table.id(),
            self.protein_string()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const DEMO: &str = "ATGAAACGTATTCTGTAA";

    fn demo() -> CodonSequence {
        CodonSequence::from_cds(DEMO, 11).unwrap()
    }

    /// TST-CORE-001 (I1, I2): length and translation agree with the declaration.
    #[test]
    fn unit_tst_core_001_i1_i2_length_and_translation() {
        let seq = demo();
        assert_eq!(seq.codon_count(), 6);
        assert_eq!(seq.cds_string().len(), 3 * seq.codon_count());
        assert_eq!(seq.protein_string(), "MKRIL*");
        assert!(seq.includes_stop());
    }

    /// TST-CORE-001 (I3): protein invariance is a precondition of `substitute`.
    #[test]
    fn unit_tst_core_001_i3_protein_invariance_is_enforced() {
        let mut seq = demo();
        let before = seq.protein_string();

        // A synonymous substitution is accepted and leaves the protein alone.
        let cgc = Codon::parse("CGC").unwrap();
        let previous = seq.substitute(2, cgc).unwrap();
        assert_eq!(previous, Codon::parse("CGT").unwrap());
        assert_eq!(seq.protein_string(), before);
        assert_eq!(seq.cds_string(), "ATGAAACGCATTCTGTAA");

        // A non-synonymous one is refused, and the sequence is untouched.
        let aaa = Codon::parse("AAA").unwrap();
        let err = seq.substitute(2, aaa).unwrap_err();
        assert_eq!(err.code, code::SEQ_NOT_SYNONYMOUS);
        assert_eq!(seq.cds_string(), "ATGAAACGCATTCTGTAA");
        assert_eq!(seq.protein_string(), before);
    }

    #[test]
    fn unit_tst_core_001_i3_holds_over_exhaustive_substitution() {
        // Every legal substitution at every editable position, checked against
        // the translation rather than against the cached protein.
        let seq = demo();
        let expected = seq.protein_string();
        let mut applied = 0usize;
        for index in seq.editable_positions() {
            for &alternative in seq.synonyms_at(index) {
                let mutated = seq.substituted(index, alternative).unwrap();
                let retranslated: String =
                    String::from_utf8(mutated.table().translate(mutated.codons())).unwrap();
                assert_eq!(
                    retranslated, expected,
                    "protein changed at {index} -> {alternative}"
                );
                applied += 1;
            }
        }
        assert!(
            applied >= 12,
            "expected several alternatives, tried {applied}"
        );
    }

    /// TST-CORE-001 (I4): locks are range-checked and enforced.
    #[test]
    fn unit_tst_core_001_i4_locks() {
        let mut seq = demo();
        // Defaults: start codon and terminal stop.
        assert!(seq.is_locked(0));
        assert!(seq.is_locked(5));
        assert!(!seq.is_locked(2));

        let ttg = Codon::parse("TTG").unwrap(); // synonymous with CTG (Leu)
        assert_eq!(
            seq.substitute(0, Codon::parse("ATG").unwrap())
                .unwrap_err()
                .code,
            code::SEQ_POSITION_LOCKED
        );

        assert_eq!(
            seq.add_lock(4, 4, LockReason::UserLocked).unwrap_err().code,
            code::SEQ_LOCK_OUT_OF_RANGE
        );
        assert_eq!(
            seq.add_lock(4, 99, LockReason::UserLocked)
                .unwrap_err()
                .code,
            code::SEQ_LOCK_OUT_OF_RANGE
        );

        assert!(seq.substitute(4, ttg).is_ok());
        seq.add_lock(4, 5, LockReason::UserLocked).unwrap();
        assert_eq!(
            seq.substitute(4, Codon::parse("CTG").unwrap())
                .unwrap_err()
                .code,
            code::SEQ_POSITION_LOCKED
        );

        seq.clear_locks();
        assert!(!seq.is_locked(0));
    }

    /// TST-CORE-001 (I5): the digest is over the document minus the digest.
    #[test]
    fn unit_tst_core_001_i5_digest_round_trip() {
        let seq = demo()
            .with_five_prime_flank("AAAGAGGAGAAAGATCT")
            .unwrap()
            .with_sequence_id("seq_demo");
        let doc = seq.to_json().unwrap();
        let declared = doc["digest"].as_str().unwrap();
        let recomputed = Digest::of_document(&doc, "digest").unwrap();
        assert_eq!(declared, recomputed.as_str());

        // Round-trip through the document form.
        let parsed = CodonSequence::from_json(&doc).unwrap();
        assert_eq!(parsed.cds_string(), seq.cds_string());
        assert_eq!(parsed.protein_string(), seq.protein_string());
        assert_eq!(parsed.five_prime_flank(), seq.five_prime_flank());
        assert_eq!(parsed.locked_spans().len(), seq.locked_spans().len());
        assert_eq!(parsed.to_json().unwrap(), doc);
    }

    #[test]
    fn unit_from_json_rejects_a_tampered_document() {
        let seq = demo();
        let mut doc = seq.to_json().unwrap();

        // Protein tampered: caught by I2 before the digest is even consulted.
        doc["protein"] = serde_json::json!("MKRIQ*");
        assert_eq!(
            CodonSequence::from_json(&doc).unwrap_err().code,
            code::SEQ_PROTEIN_MISMATCH
        );

        // Digest tampered while everything else is consistent: caught by I5.
        let mut doc = seq.to_json().unwrap();
        doc["digest"] = serde_json::json!(
            "sha256:0000000000000000000000000000000000000000000000000000000000000000"
        );
        assert_eq!(
            CodonSequence::from_json(&doc).unwrap_err().code,
            code::SEQ_DIGEST_MISMATCH
        );

        // Codon count tampered: caught by I1.
        let mut doc = seq.to_json().unwrap();
        doc["codon_count"] = serde_json::json!(7);
        assert_eq!(
            CodonSequence::from_json(&doc).unwrap_err().code,
            code::SEQ_NOT_IN_FRAME
        );
    }

    #[test]
    fn unit_frame_and_alphabet_are_rejected_at_construction() {
        assert_eq!(
            CodonSequence::from_cds("ATGAAAC", 11).unwrap_err().code,
            code::SEQ_NOT_IN_FRAME
        );
        assert_eq!(
            CodonSequence::from_cds("ATGNNNTAA", 11).unwrap_err().code,
            code::SEQ_INVALID_BASE
        );
        assert_eq!(
            CodonSequence::from_cds("", 11).unwrap_err().code,
            code::SEQ_EMPTY
        );
        assert_eq!(
            CodonSequence::from_cds("ATGTAAAAATAA", 11)
                .unwrap_err()
                .code,
            code::SEQ_INTERNAL_STOP
        );
        assert_eq!(
            CodonSequence::from_cds(DEMO, 7).unwrap_err().code,
            code::SEQ_UNKNOWN_TABLE
        );
    }

    #[test]
    fn unit_rna_input_is_transcoded() {
        let seq = CodonSequence::from_cds("AUGAAACGUAUUCUGUAA", 11).unwrap();
        assert_eq!(seq.cds_string(), DEMO);
        assert_eq!(seq.protein_string(), "MKRIL*");
    }

    #[test]
    fn unit_digest_changes_with_content_and_is_cached_consistently() {
        let seq = demo();
        let d1 = seq.digest().unwrap();
        let d2 = seq.digest().unwrap();
        assert_eq!(d1, d2);

        let mutated = seq.substituted(2, Codon::parse("CGC").unwrap()).unwrap();
        assert_ne!(seq.digest().unwrap(), mutated.digest().unwrap());
    }

    #[test]
    fn unit_table_without_stop_is_supported() {
        // A fusion partner: no terminal stop, so includes_stop is false and any
        // stop anywhere is an internal stop.
        let seq = CodonSequence::from_cds("ATGAAACGT", 11).unwrap();
        assert!(!seq.includes_stop());
        assert_eq!(seq.protein_string(), "MKR");
        assert_eq!(
            seq.locked_spans().len(),
            1,
            "only the start codon is locked"
        );
    }
}
