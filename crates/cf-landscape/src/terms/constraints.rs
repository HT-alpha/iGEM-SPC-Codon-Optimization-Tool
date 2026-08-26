//! `constraint.violation_count`.
//!
//! Ticket P1-21.
//!
//! ## The capability is a property of the constraint set, not of this code
//!
//! Every other Phase 1 term is L1 unconditionally, because CAI and GC are
//! functions of the codons and nothing else. Violation count is not: whether an
//! exact delta is affordable depends on what the user asked for.
//!
//! A forbidden motif, a homopolymer limit and a GC window all have bounded
//! reach — editing one codon can only create or destroy a violation within a
//! computable distance. A k-mer repeat limit does not: it counts occurrences
//! over the whole construct, so changing three nucleotides here can push a
//! k-mer somewhere else across its threshold, and there is no window that
//! bounds the effect.
//!
//! So [`ViolationCountTerm::capability`] answers L1 for a local set and L0 for
//! one containing a k-mer rule, and the composition layer degrades that term —
//! only that term — to a full rescore. The alternative would be to declare L1
//! and compute the "delta" by rescoring twice, which is exact and costs `O(L)`;
//! that reports a cheap path while walking the expensive one, and the whole
//! argument for exact deltas is a cost argument.
//!
//! ## What is counted
//!
//! Actionable blocking violations. `DEP-CONS-002`: a blocking violation
//! confined to immutable sequence is a property of the input, not of the
//! candidate, so counting it would add the same constant to every candidate
//! while making all of them look infeasible. Advisory violations are excluded
//! because they do not block, and mixing the two would make one number mean two
//! things.

use cf_constraints::{ConstraintChecker, ConstructView, KmerMasking};
use cf_landscape_api::{
    Capability, Edit, EvalContext, PositionDeltas, Provenance, Rebased, TermProvider, TermRole,
};
use cf_types::{code, CfError, CodonSequence, Direction, Result, TermId, Unit};

use crate::construct::ConstructLayout;

/// Counts actionable blocking violations over the construct the CDS sits in.
#[derive(Debug)]
pub struct ViolationCountTerm {
    term: TermId,
    checker: ConstraintChecker,
    layout: ConstructLayout,
    masking: KmerMasking,
}

impl ViolationCountTerm {
    pub fn new(checker: ConstraintChecker, layout: ConstructLayout) -> Result<Self> {
        Ok(Self {
            term: TermId::new("constraint.violation_count")?,
            checker,
            layout,
            // Legacy parity, per EXIT-001/EXIT-002. Not this ticket's decision
            // to revisit.
            masking: KmerMasking::default(),
        })
    }

    pub fn with_masking(mut self, masking: KmerMasking) -> Self {
        self.masking = masking;
        self
    }

    /// Constraints the checker cannot evaluate at all.
    ///
    /// Surfaced so composition can report them: a constraint that is silently
    /// not checked is worse than one that fails loudly, and this term's value
    /// would otherwise look like a complete answer.
    pub fn undetectable(&self) -> &[String] {
        self.checker.undetectable()
    }

    fn count(&self, construct: &[u8]) -> usize {
        let view =
            ConstructView::new(construct, self.layout.cds_spans()).with_masking(self.masking);
        self.checker
            .check(&view)
            .iter()
            .filter(|v| v.class == cf_constraints::spec::ConstraintClass::Blocking && v.actionable)
            .count()
    }

    /// Actionable blocking violations intersecting a nucleotide span.
    fn count_region(&self, construct: &[u8], start: usize, end: usize) -> Option<usize> {
        let view =
            ConstructView::new(construct, self.layout.cds_spans()).with_masking(self.masking);
        Some(
            self.checker
                .check_region(&view, start, end)?
                .iter()
                .filter(|v| {
                    v.class == cf_constraints::spec::ConstraintClass::Blocking && v.actionable
                })
                .count(),
        )
    }

    /// Nucleotide span touched by a set of codon edits, in construct
    /// coordinates.
    fn touched_span(&self, edits: &[Edit], codon_count: usize) -> Result<(usize, usize)> {
        let mut lo = usize::MAX;
        let mut hi = 0usize;
        for edit in edits {
            if edit.codon_index >= codon_count {
                return Err(CfError::new(
                    code::SEQ_INDEX_OUT_OF_RANGE,
                    format!("codon {} is outside the sequence", edit.codon_index),
                ));
            }
            let (s, e) = self.layout.codon_span_nt(edit.codon_index);
            lo = lo.min(s);
            hi = hi.max(e);
        }
        Ok((lo, hi))
    }

    fn delta_from(
        &self,
        sequence: &CodonSequence,
        base_construct: &[u8],
        edits: &[Edit],
    ) -> Result<f64> {
        if edits.is_empty() {
            return Ok(0.0);
        }
        let (start, end) = self.touched_span(edits, sequence.codon_count())?;

        let Some(before) = self.count_region(base_construct, start, end) else {
            return Err(CfError::new(
                code::LAND_CAPABILITY,
                "constraint.violation_count is L0 for this constraint set: a k-mer repeat \
                 rule counts over the whole construct, so no window bounds a substitution's \
                 effect",
            ));
        };

        // The construct is patched rather than reassembled from the sequence:
        // the edits are synonymous, so they change exactly the nucleotides of
        // the codons they name and nothing else moves.
        let mut mutated = base_construct.to_vec();
        for edit in edits {
            let (lo, _) = self.layout.codon_span_nt(edit.codon_index);
            mutated[lo..lo + 3].copy_from_slice(&edit.codon.to_bytes());
        }
        let after = self
            .count_region(&mutated, start, end)
            .expect("locality was already established");

        Ok(after as f64 - before as f64)
    }

    /// The pinned construct, or a freshly assembled one.
    fn base_construct(&self, sequence: &CodonSequence, base: &Rebased) -> Vec<u8> {
        base.state::<Vec<u8>>()
            .cloned()
            .unwrap_or_else(|| self.layout.assemble(sequence))
    }
}

impl TermProvider for ViolationCountTerm {
    fn term_id(&self) -> &TermId {
        &self.term
    }
    fn unit(&self) -> Unit {
        Unit::new(Unit::DIMENSIONLESS)
    }
    fn direction(&self) -> Direction {
        Direction::Minimize
    }
    fn capability(&self) -> Capability {
        if self.checker.influence_radius_nt().is_some() {
            Capability::L1
        } else {
            Capability::L0
        }
    }
    fn max_role(&self) -> TermRole {
        TermRole::ParetoAxis
    }
    fn provenance(&self) -> Provenance {
        Provenance::builtin("constraint_native", crate::ENGINE_VERSION)
    }
    fn score(&self, ctx: &EvalContext<'_>) -> Result<f64> {
        Ok(self.count(&self.layout.assemble(ctx.sequence)) as f64)
    }
    /// Pins the assembled construct.
    ///
    /// This removes the automaton walk and the k-mer index build from the
    /// delta path, leaving one `O(L)` byte copy per delta. That copy is
    /// several orders of magnitude cheaper than what it replaces but it is not
    /// `O(W)`, and saying so is more useful than implying otherwise: see the
    /// note on `DEP-LAND-003`.
    fn rebase(&self, ctx: &EvalContext<'_>) -> Result<Rebased> {
        let construct = self.layout.assemble(ctx.sequence);
        Ok(Rebased::with_state(
            self.count(&construct) as f64,
            construct,
        ))
    }
    fn delta(&self, ctx: &EvalContext<'_>, base: &Rebased, edits: &[Edit]) -> Result<f64> {
        let construct = self.base_construct(ctx.sequence, base);
        self.delta_from(ctx.sequence, &construct, edits)
    }
    fn delta_table(
        &self,
        ctx: &EvalContext<'_>,
        base: &Rebased,
        positions: &[usize],
    ) -> Result<Vec<PositionDeltas>> {
        let construct = self.base_construct(ctx.sequence, base);
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
                        self.delta_from(ctx.sequence, &construct, &[Edit::new(index, c)])
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
