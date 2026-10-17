//! `epistasis.context_naturalness` as a ranked term.
//!
//! The value is the mean conditional log probability under a fitted
//! synonymous Potts model. Higher is better. The dwell head is not this term
//! and is not scored here.
//!
//! The delta recomputes only the positions inside the fitted window of an
//! edit (`mean_log_probability_delta`). The codon buffer is updated in place
//! for that window and restored; it is not cloned per synonym. A full rescore
//! is the check, not the implementation: the two have to agree, and the window
//! sum is what tier 1 can afford.

use cf_landscape_api::{
    Capability, Edit, EvalContext, InteractionRadius, PositionDeltas, Provenance, Rebased,
    ResolvedBy, TermProvider, TermRole,
};
use cf_types::{code, ArtifactId, CfError, Codon, Digest, Direction, Result, TermId, Unit};

use crate::potts::{mean_log_probability, PottsField, PottsModel};

/// Ranked context naturalness. Maximize. L1.
#[derive(Clone, Debug)]
pub struct ContextNaturalnessTerm {
    term: TermId,
    model: PottsModel,
    artifact_id: ArtifactId,
    digest: Digest,
}

impl ContextNaturalnessTerm {
    /// `model` must already be the verified parameter set.
    ///
    /// `artifact_id` and `digest` are recorded so a score can be tied back to
    /// the document that produced the parameters. They are not checked again
    /// here; the caller verified the artifact before constructing the term.
    pub fn new(model: PottsModel, artifact_id: &str, digest: &str) -> Result<Self> {
        let width = 64 + model.window * 4096;
        if model.h.len() != 64 || model.couplings.len() != model.window * 4096 {
            return Err(CfError::new(
                code::LAND_TERM_UNKNOWN,
                format!(
                    "Potts parameters have h {} and {} couplings; a window of {} needs {width}",
                    model.h.len(),
                    model.couplings.len(),
                    model.window
                ),
            ));
        }
        Ok(Self {
            term: TermId::new("epistasis.context_naturalness")?,
            model,
            artifact_id: ArtifactId::new(artifact_id)?,
            digest: Digest::parse(digest)?,
        })
    }

    fn codons_of(&self, sequence: &cf_types::CodonSequence) -> Vec<u8> {
        sequence
            .codons()
            .iter()
            .map(|codon| codon.index())
            .collect()
    }

    fn score_codons(&self, codons: &[u8]) -> Result<f64> {
        mean_log_probability(codons, &self.model).map_err(|detail| {
            CfError::new(
                code::LAND_TERM_UNKNOWN,
                format!("epistasis.context_naturalness could not score the sequence: {detail}"),
            )
        })
    }
}

impl TermProvider for ContextNaturalnessTerm {
    fn term_id(&self) -> &TermId {
        &self.term
    }

    fn unit(&self) -> Unit {
        Unit::new(Unit::DIMENSIONLESS)
    }

    fn direction(&self) -> Direction {
        Direction::Maximize
    }

    fn capability(&self) -> Capability {
        Capability::L1
    }

    fn max_role(&self) -> TermRole {
        TermRole::ParetoAxis
    }

    fn interaction_radius(&self) -> InteractionRadius {
        InteractionRadius::Codons(self.model.window)
    }

    fn provenance(&self) -> Provenance {
        let mut provenance = Provenance::builtin("potts_context", env!("CARGO_PKG_VERSION"));
        provenance.resolved_by = ResolvedBy::Pin;
        provenance.artifact_id = Some(self.artifact_id.clone());
        provenance.digest = Some(self.digest.clone());
        provenance
    }

    fn score(&self, ctx: &EvalContext<'_>) -> Result<f64> {
        self.score_codons(&self.codons_of(ctx.sequence))
    }

    fn rebase(&self, ctx: &EvalContext<'_>) -> Result<Rebased> {
        Ok(Rebased::value_only(self.score(ctx)?))
    }

    fn delta(&self, ctx: &EvalContext<'_>, base: &Rebased, edits: &[Edit]) -> Result<f64> {
        let _ = base;
        let mut codons = self.codons_of(ctx.sequence);
        let applied = applied_edits(ctx, &codons, edits)?;
        let field = PottsField::bind(&self.model, &codons).map_err(naturalness_delta_error)?;
        field
            .delta(&mut codons, &applied)
            .map_err(naturalness_delta_error)
    }

    fn delta_table(
        &self,
        ctx: &EvalContext<'_>,
        _base: &Rebased,
        positions: &[usize],
    ) -> Result<Vec<PositionDeltas>> {
        let mut codons = self.codons_of(ctx.sequence);
        let field = PottsField::bind(&self.model, &codons).map_err(naturalness_delta_error)?;
        positions
            .iter()
            .map(|&index| {
                if index >= codons.len() {
                    return Err(CfError::new(
                        code::SEQ_INDEX_OUT_OF_RANGE,
                        format!("codon {index} is outside the sequence"),
                    ));
                }
                let current = ctx.sequence.codons()[index];
                let alternatives = ctx
                    .sequence
                    .synonyms_at(index)
                    .iter()
                    .copied()
                    .filter(|codon| *codon != current)
                    .map(|codon| {
                        let change = field
                            .delta(&mut codons, &[(index, codon.index())])
                            .map_err(naturalness_delta_error)?;
                        Ok((codon, change))
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

fn applied_edits(ctx: &EvalContext<'_>, codons: &[u8], edits: &[Edit]) -> Result<Vec<(usize, u8)>> {
    let table = ctx.sequence.table();
    let mut applied = Vec::with_capacity(edits.len());
    for edit in edits {
        if edit.codon_index >= codons.len() {
            return Err(CfError::new(
                code::SEQ_INDEX_OUT_OF_RANGE,
                format!("codon {} is outside the sequence", edit.codon_index),
            ));
        }
        let old = Codon::from_index(codons[edit.codon_index])?;
        if table.amino_acid(old) != table.amino_acid(edit.codon) {
            return Err(CfError::new(
                code::SEQ_NOT_SYNONYMOUS,
                format!("{old} and {} do not encode the same amino acid", edit.codon),
            ));
        }
        applied.push((edit.codon_index, edit.codon.index()));
    }
    Ok(applied)
}

fn naturalness_delta_error(detail: String) -> CfError {
    CfError::new(
        code::LAND_TERM_UNKNOWN,
        format!("epistasis.context_naturalness delta failed: {detail}"),
    )
}
