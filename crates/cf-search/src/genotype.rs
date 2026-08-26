//! The genotype, and the move type that connects an operator to the landscape.
//!
//! Ticket P1-23.
//!
//! ## Invariance is inherited, not re-implemented
//!
//! The design decision that makes the whole search cheap is in
//! [`03-nsga-loop-and-operators.md`](../../../planning/03-algorithms/03-nsga-loop-and-operators.md)
//! §1: because the genotype is a product space over synonymous choices, *every*
//! codon-boundary operator produces a valid individual. There is no repair, no
//! penalty and no rejection sampling for the protein constraint.
//!
//! That is only true if nothing can write a codon outside its synonymous
//! family, so [`Genotype`] exposes no mutable view of its codons. The single
//! mutation path is [`Genotype::apply`], which delegates to
//! `CodonSequence::substitute` — already the enforcement point for synonymy
//! (`CF-SEQ-004`) and for locked spans (`CF-SEQ-006`) since `P1-06`. An
//! operator therefore cannot break either invariant without first bypassing a
//! type that offers no way to bypass it.
//!
//! ## A move is the landscape's delta query, not a translation of one
//!
//! [`Move`] carries `cf_landscape_api::Edit` values — the same type
//! `Target::Delta` holds. An operator's output is already the shape the exact-
//! delta path consumes, so `P1-24`'s tier-1 DP gets the cheap path by handing
//! the move over rather than by reconstructing an edit list from a diff. A
//! genotype that only recorded its final codons would have forced exactly that
//! reconstruction, and a diff of two sequences is `O(L)` — which would have
//! put the query back on the expensive path by construction.

use cf_landscape_api::{Edit, Target};
use cf_types::{code, CfError, Codon, CodonSequence, Digest, Result};

/// Which operator produced a move. Diagnostic only; nothing branches on it.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub enum MoveKind {
    UniformSynonymous,
    PottsGuided,
    BlockGibbs,
    ViolationDirected,
    Crossover,
    /// A move assembled by a caller rather than by a registered operator.
    Explicit,
}

impl MoveKind {
    pub fn as_str(self) -> &'static str {
        match self {
            MoveKind::UniformSynonymous => "uniform_synonymous",
            MoveKind::PottsGuided => "potts_guided",
            MoveKind::BlockGibbs => "block_gibbs",
            MoveKind::ViolationDirected => "violation_directed",
            MoveKind::Crossover => "crossover",
            MoveKind::Explicit => "explicit",
        }
    }
}

/// A set of substitutions, relative to the genotype it was proposed against.
#[derive(Clone, PartialEq, Debug)]
pub struct Move {
    edits: Vec<Edit>,
    kind: MoveKind,
}

impl Move {
    pub fn new(kind: MoveKind, edits: Vec<Edit>) -> Self {
        Self { edits, kind }
    }

    pub fn empty(kind: MoveKind) -> Self {
        Self {
            edits: Vec::new(),
            kind,
        }
    }

    /// The edits, in the exact type `Target::Delta` carries.
    pub fn edits(&self) -> &[Edit] {
        &self.edits
    }

    pub fn kind(&self) -> MoveKind {
        self.kind
    }

    pub fn len(&self) -> usize {
        self.edits.len()
    }

    pub fn is_empty(&self) -> bool {
        self.edits.is_empty()
    }

    /// Whether this is the single-substitution move tier 2 accept/reject and
    /// the repair operator are built around.
    pub fn is_single_substitution(&self) -> bool {
        self.edits.len() == 1
    }

    /// Codon indices this move touches, in ascending order without duplicates.
    pub fn touched(&self) -> Vec<usize> {
        let mut indices: Vec<usize> = self.edits.iter().map(|e| e.codon_index).collect();
        indices.sort_unstable();
        indices.dedup();
        indices
    }

    /// The landscape query target for this move against a known base.
    ///
    /// The reason this method exists rather than a comment saying "the caller
    /// can build one": it makes the interop a compile-checked property. If the
    /// move type ever stopped carrying edits, this would not build.
    pub fn to_delta_target(
        &self,
        base_digest: Digest,
        cistron_id: Option<String>,
    ) -> Result<Target> {
        if self.edits.is_empty() {
            return Err(CfError::new(
                code::LAND_CAPABILITY,
                "an empty move has no delta target; score the base instead",
            ));
        }
        Ok(Target::Delta {
            base_digest,
            edits: self.edits.clone(),
            cistron_id,
        })
    }

    /// Drop edits that would not change anything.
    ///
    /// A no-op edit is not wrong, but it makes a move look larger than it is
    /// and inflates the delta query it produces.
    pub fn compacted(mut self, against: &Genotype) -> Self {
        let codons = against.codons();
        self.edits
            .retain(|e| codons.get(e.codon_index).is_some_and(|c| *c != e.codon));
        self
    }
}

/// One individual: a coding sequence plus an incrementally maintained hash.
///
/// The hash is a filter for duplicate suppression, never a decision on its own:
/// [`Genotype::same_codons_as`] confirms a match, so a collision costs one
/// regeneration rather than silently discarding a distinct offspring.
#[derive(Clone, Debug)]
pub struct Genotype {
    sequence: CodonSequence,
    hash: u64,
}

impl Genotype {
    pub fn new(sequence: CodonSequence) -> Self {
        let hash = Self::full_hash(sequence.codons());
        Self { sequence, hash }
    }

    pub fn sequence(&self) -> &CodonSequence {
        &self.sequence
    }

    pub fn codons(&self) -> &[Codon] {
        self.sequence.codons()
    }

    pub fn len(&self) -> usize {
        self.sequence.codon_count()
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.codon_count() == 0
    }

    pub fn hash(&self) -> u64 {
        self.hash
    }

    pub fn same_codons_as(&self, other: &Genotype) -> bool {
        self.codons() == other.codons()
    }

    /// Positions an operator may sample from: editable and with a real choice.
    ///
    /// A position whose amino acid has one codon is excluded even though it is
    /// unlocked, because sampling it would waste a draw on a position with no
    /// alternative — methionine and tryptophan under the standard code.
    pub fn mutable_positions(&self) -> Vec<usize> {
        self.sequence
            .editable_positions()
            .filter(|i| self.sequence.synonyms_at(*i).len() > 1)
            .collect()
    }

    /// Nucleotide offset of the coding region inside its construct.
    ///
    /// Taken from the sequence's own 5' flank rather than from a separate
    /// layout type, so there is one place the offset can come from and no way
    /// for two of them to disagree.
    pub fn cds_offset_nt(&self) -> usize {
        self.sequence.five_prime_flank().len()
    }

    /// Half-open nucleotide span of one codon, in construct coordinates.
    pub fn codon_span_nt(&self, codon_index: usize) -> (usize, usize) {
        let start = self.cds_offset_nt() + codon_index * 3;
        (start, start + 3)
    }

    /// Codon indices overlapping a nucleotide span, clipped to this sequence.
    pub fn codons_in_span(&self, start_nt: usize, end_nt: usize) -> Vec<usize> {
        let offset = self.cds_offset_nt();
        let end_of_cds = offset + self.len() * 3;
        if end_nt <= offset || start_nt >= end_of_cds {
            return Vec::new();
        }
        let first = start_nt.max(offset).saturating_sub(offset) / 3;
        let last = (end_nt.min(end_of_cds) - offset).div_ceil(3);
        (first..last.min(self.len())).collect()
    }

    /// Apply a move in place.
    ///
    /// The only way to change a genotype. Every edit goes through
    /// `CodonSequence::substitute`, so a non-synonymous replacement fails with
    /// `CF-SEQ-004` and a locked position with `CF-SEQ-006`, and neither check
    /// is restated here — restating it would create a second definition that
    /// could drift.
    ///
    /// Not atomic: a move that fails partway leaves earlier edits applied. That
    /// is deliberate rather than overlooked. Making it atomic means cloning the
    /// sequence on every application, which is `O(L)` on the hottest path in
    /// the engine, to guard against a failure that operators cannot produce —
    /// they sample from `synonyms_at` over `mutable_positions`. Callers that
    /// accept moves from outside the operator set use [`Genotype::with_move`],
    /// which clones once and is atomic by construction.
    pub fn apply(&mut self, mv: &Move) -> Result<()> {
        for edit in mv.edits() {
            let previous = self.sequence.substitute(edit.codon_index, edit.codon)?;
            self.hash ^= Self::mix(edit.codon_index, previous);
            self.hash ^= Self::mix(edit.codon_index, edit.codon);
        }
        Ok(())
    }

    /// A copy with the move applied, or an error and no change to `self`.
    pub fn with_move(&self, mv: &Move) -> Result<Self> {
        let mut next = self.clone();
        next.apply(mv)?;
        Ok(next)
    }

    /// The move that turns this genotype into `other`.
    ///
    /// `O(L)`. Provided for seeds and for imported candidates, which arrive as
    /// whole sequences with no edit history; it is not how operators produce
    /// moves, because paying `O(L)` per move is the cost the delta path exists
    /// to avoid.
    pub fn diff(&self, other: &Genotype) -> Result<Move> {
        if self.len() != other.len() {
            return Err(CfError::new(
                code::SEQ_INDEX_OUT_OF_RANGE,
                format!(
                    "cannot diff a {}-codon genotype against a {}-codon one",
                    self.len(),
                    other.len()
                ),
            ));
        }
        let edits = self
            .codons()
            .iter()
            .zip(other.codons())
            .enumerate()
            .filter(|(_, (a, b))| a != b)
            .map(|(index, (_, b))| Edit::new(index, *b))
            .collect();
        Ok(Move::new(MoveKind::Explicit, edits))
    }

    /// Position-sensitive, incrementally updatable hash.
    ///
    /// XOR of a per-`(position, codon)` mix, so a substitution updates it in
    /// `O(1)` by cancelling the old contribution. A rolling polynomial hash
    /// would have needed an `O(L)` recomputation per move, which duplicate
    /// suppression runs on every offspring.
    fn mix(index: usize, codon: Codon) -> u64 {
        let mut value = (index as u64) << 6 | u64::from(codon.index());
        value = value.wrapping_add(0x9E37_79B9_7F4A_7C15);
        value ^= value >> 30;
        value = value.wrapping_mul(0xBF58_476D_1CE4_E5B9);
        value ^= value >> 27;
        value = value.wrapping_mul(0x94D0_49BB_1331_11EB);
        value ^ (value >> 31)
    }

    fn full_hash(codons: &[Codon]) -> u64 {
        codons
            .iter()
            .enumerate()
            .fold(0u64, |acc, (index, codon)| acc ^ Self::mix(index, *codon))
    }

    /// Whether the cached hash still matches a full recomputation.
    ///
    /// Exposed for tests. The incremental update is the kind of optimization
    /// that stays correct until someone adds a second mutation path, and a
    /// stale hash shows up as a duplicate that was never suppressed — a
    /// symptom nobody traces back to here.
    pub fn hash_is_consistent(&self) -> bool {
        self.hash == Self::full_hash(self.codons())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn genotype() -> Genotype {
        Genotype::new(CodonSequence::from_cds("ATGCTGCGTAGCGGCGTGTAA", 1).unwrap())
    }

    #[test]
    fn unit_apply_updates_the_hash_incrementally() {
        let mut g = genotype();
        let before = g.hash();
        g.apply(&Move::new(
            MoveKind::Explicit,
            vec![Edit::new(1, Codon::parse("CTC").unwrap())],
        ))
        .unwrap();
        assert_ne!(g.hash(), before);
        assert!(g.hash_is_consistent());

        // And back again: the XOR construction must cancel exactly.
        g.apply(&Move::new(
            MoveKind::Explicit,
            vec![Edit::new(1, Codon::parse("CTG").unwrap())],
        ))
        .unwrap();
        assert_eq!(g.hash(), before);
    }

    #[test]
    fn unit_a_non_synonymous_move_is_refused_by_the_core() {
        let mut g = genotype();
        let err = g
            .apply(&Move::new(
                MoveKind::Explicit,
                vec![Edit::new(1, Codon::parse("ATG").unwrap())],
            ))
            .unwrap_err();
        assert_eq!(err.code, "CF-SEQ-004");
    }

    #[test]
    fn unit_mutable_positions_exclude_single_codon_families() {
        // ATG is methionine and TAA is a stop; neither has an alternative.
        let g = genotype();
        let mutable = g.mutable_positions();
        assert!(!mutable.contains(&0));
        assert!(!mutable.contains(&(g.len() - 1)));
        assert!(mutable.contains(&1));
    }

    #[test]
    fn unit_span_mapping_clips_to_the_coding_region() {
        let g = Genotype::new(
            CodonSequence::from_cds("ATGCTGCGTAGCGGCGTGTAA", 1)
                .unwrap()
                .with_five_prime_flank("AAAAAA")
                .unwrap(),
        );
        assert_eq!(g.cds_offset_nt(), 6);
        assert_eq!(g.codon_span_nt(0), (6, 9));
        // Entirely in the flank.
        assert!(g.codons_in_span(0, 6).is_empty());
        // Straddling the junction: only the coding part is returned.
        assert_eq!(g.codons_in_span(3, 12), vec![0, 1]);
    }

    #[test]
    fn unit_a_move_is_a_delta_target() {
        let mv = Move::new(
            MoveKind::UniformSynonymous,
            vec![Edit::new(1, Codon::parse("CTC").unwrap())],
        );
        assert!(mv.is_single_substitution());
        let target = mv.to_delta_target(Digest::of_bytes(b"base"), None).unwrap();
        match target {
            Target::Delta { edits, .. } => assert_eq!(edits, mv.edits()),
            other => panic!("expected a delta target, got {}", other.kind()),
        }
        assert!(Move::empty(MoveKind::Explicit)
            .to_delta_target(Digest::of_bytes(b"base"), None)
            .is_err());
    }

    #[test]
    fn unit_compaction_drops_no_op_edits() {
        let g = genotype();
        let current = g.codons()[1];
        let mv = Move::new(
            MoveKind::Explicit,
            vec![
                Edit::new(1, current),
                Edit::new(2, Codon::parse("CGC").unwrap()),
            ],
        )
        .compacted(&g);
        assert_eq!(mv.len(), 1);
        assert_eq!(mv.edits()[0].codon_index, 2);
    }

    #[test]
    fn unit_diff_round_trips_through_apply() {
        let a = genotype();
        let b = a
            .with_move(&Move::new(
                MoveKind::Explicit,
                vec![
                    Edit::new(1, Codon::parse("CTC").unwrap()),
                    Edit::new(3, Codon::parse("TCT").unwrap()),
                ],
            ))
            .unwrap();
        let mv = a.diff(&b).unwrap();
        assert_eq!(mv.len(), 2);
        assert!(a.with_move(&mv).unwrap().same_codons_as(&b));
    }
}
