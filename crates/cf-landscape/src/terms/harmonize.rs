//! Codon harmonization and synonymous-frequency matching.
//!
//! `usage.harmonize_rca` is the relative-codon-adaptiveness match described
//! by Claassens et al., PLoS ONE 2017, 12(9):e0184355
//! (doi:10.1371/journal.pone.0184355), which is the same per-codon replacement
//! Angov, Hillier, Kincaid and Lyon, PLoS ONE 2008, 3(5):e2189
//! (doi:10.1371/journal.pone.0002189) frame as matching a source usage
//! frequency with the closest host frequency. The score is the mean absolute
//! difference between the source codon's relative adaptiveness and the
//! chosen host codon's relative adaptiveness. The unconstrained minimum is
//! the closest host codon at each position. A constraint can force a
//! compromise. This is the objective class DNA Chisel names `harmonize_rca`
//! (Zulkower and Rosser, Bioinformatics 2020). It is not a claim that the
//! two implementations return the same bits.
//!
//! `usage.match_frequency` is the squared error between the sequence's
//! synonymous codon fractions and the host table's synonymous fractions,
//! summed over amino acids that have a choice. Relative adaptiveness `w`
//! recovers that fraction: `f_i = w_i / Σ w` inside the family, which is
//! the definition `w_i = f_i / max f`. Maximizing CAI piles onto the top
//! codon. This term penalizes that pile-up. It is the frequency-matching
//! objective class, not a transcription of DNA Chisel's `match_codon_usage`
//! source.
//!
//! `DEP-LAND-006`.

use cf_landscape_api::{
    Capability, Edit, EvalContext, InteractionRadius, PositionDeltas, Provenance, Rebased,
    TermProvider,
};
use cf_types::{
    code, CfError, Codon, CodonSequence, Direction, GeneticCode, Result, TermId, Unit, CODON_COUNT,
};

use super::usage::CodonUsageTable;

/// Mean absolute RCA difference. Positionwise: the source codon is fixed.
#[derive(Debug)]
pub struct HarmonizeRcaTerm {
    term: TermId,
    target: CodonUsageTable,
    /// Source relative adaptiveness of the native codon at each position.
    source_rca: Vec<f64>,
    counted: Vec<bool>,
}

#[derive(Clone, Copy, Debug)]
struct HarmState {
    sum_abs: f64,
    counted: usize,
}

impl HarmonizeRcaTerm {
    pub fn new(
        target: CodonUsageTable,
        source: CodonUsageTable,
        native: &CodonSequence,
    ) -> Result<Self> {
        let mut source_rca = Vec::with_capacity(native.codon_count());
        let mut counted = Vec::with_capacity(native.codon_count());
        for index in 0..native.codon_count() {
            let codon = native.codons()[index];
            let family = native.synonyms_at(index);
            let is_counted = !native.table().is_stop(codon) && family.len() > 1;
            if is_counted {
                let observed = family.iter().any(|member| target.weight(*member) > 0.0);
                if !observed {
                    return Err(CfError::new(
                        code::HOST_TERM_UNSUPPORTED,
                        format!(
                            "usage.harmonize_rca cannot rank position {index}: the expression \
                             host has no observed codon for that amino acid, so every synonym \
                             scores the same"
                        ),
                    ));
                }
            }
            source_rca.push(source.weight(codon));
            counted.push(is_counted);
        }
        Ok(Self {
            term: TermId::new("usage.harmonize_rca")?,
            target,
            source_rca,
            counted,
        })
    }

    fn state_of(&self, sequence: &CodonSequence) -> Result<HarmState> {
        self.require_same_length(sequence)?;
        let mut state = HarmState {
            sum_abs: 0.0,
            counted: 0,
        };
        for index in 0..sequence.codon_count() {
            if !self.counted[index] {
                continue;
            }
            state.counted += 1;
            state.sum_abs += self.deviation(sequence.codons()[index], index);
        }
        Ok(state)
    }

    fn deviation(&self, codon: Codon, index: usize) -> f64 {
        (self.target.weight(codon) - self.source_rca[index]).abs()
    }

    fn value(state: HarmState) -> f64 {
        if state.counted == 0 {
            0.0
        } else {
            state.sum_abs / state.counted as f64
        }
    }

    fn require_same_length(&self, sequence: &CodonSequence) -> Result<()> {
        if sequence.codon_count() == self.source_rca.len() {
            Ok(())
        } else {
            Err(CfError::new(
                code::SEQ_PROTEIN_MISMATCH,
                format!(
                    "usage.harmonize_rca was built on a native CDS of {} codons and scored \
                     a sequence of {} codons",
                    self.source_rca.len(),
                    sequence.codon_count()
                ),
            ))
        }
    }

    fn delta_from(
        &self,
        sequence: &CodonSequence,
        before: HarmState,
        edits: &[Edit],
    ) -> Result<f64> {
        self.require_same_length(sequence)?;
        let mut codons: Vec<Codon> = sequence.codons().to_vec();
        let mut sum = before.sum_abs;
        for edit in edits {
            if edit.codon_index >= codons.len() {
                return Err(CfError::new(
                    code::SEQ_INDEX_OUT_OF_RANGE,
                    format!("codon {} is outside the sequence", edit.codon_index),
                ));
            }
            if self.counted[edit.codon_index] {
                sum -= self.deviation(codons[edit.codon_index], edit.codon_index);
                sum += self.deviation(edit.codon, edit.codon_index);
            }
            codons[edit.codon_index] = edit.codon;
        }
        Ok(Self::value(HarmState {
            sum_abs: sum,
            counted: before.counted,
        }) - Self::value(before))
    }
}

impl TermProvider for HarmonizeRcaTerm {
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
        Capability::L1
    }
    fn interaction_radius(&self) -> InteractionRadius {
        InteractionRadius::Positionwise
    }
    fn provenance(&self) -> Provenance {
        Provenance::builtin("harmonize_rca", crate::ENGINE_VERSION)
    }
    fn score(&self, ctx: &EvalContext<'_>) -> Result<f64> {
        Ok(Self::value(self.state_of(ctx.sequence)?))
    }
    fn rebase(&self, ctx: &EvalContext<'_>) -> Result<Rebased> {
        let state = self.state_of(ctx.sequence)?;
        Ok(Rebased::with_state(Self::value(state), state))
    }
    fn delta(&self, ctx: &EvalContext<'_>, base: &Rebased, edits: &[Edit]) -> Result<f64> {
        let state = base
            .state::<HarmState>()
            .copied()
            .map(Ok)
            .unwrap_or_else(|| self.state_of(ctx.sequence))?;
        self.delta_from(ctx.sequence, state, edits)
    }
    fn delta_table(
        &self,
        ctx: &EvalContext<'_>,
        base: &Rebased,
        positions: &[usize],
    ) -> Result<Vec<PositionDeltas>> {
        let state = base
            .state::<HarmState>()
            .copied()
            .map(Ok)
            .unwrap_or_else(|| self.state_of(ctx.sequence))?;
        positions
            .iter()
            .map(|&index| {
                let current = ctx.sequence.codons()[index];
                let alternatives = ctx
                    .sequence
                    .synonyms_at(index)
                    .iter()
                    .copied()
                    .filter(|codon| *codon != current)
                    .map(|codon| {
                        self.delta_from(ctx.sequence, state, &[Edit::new(index, codon)])
                            .map(|delta| (codon, delta))
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

/// Squared error of synonymous fractions. Global: every codon of an amino
/// acid moves the fraction of every other codon of that amino acid.
#[derive(Debug)]
pub struct MatchFrequencyTerm {
    term: TermId,
    code: &'static GeneticCode,
    target_fraction: [f64; CODON_COUNT],
}

#[derive(Clone, Copy, Debug)]
struct FreqState {
    counts: [u32; CODON_COUNT],
}

impl MatchFrequencyTerm {
    pub fn new(table: CodonUsageTable, translation_table: u8) -> Result<Self> {
        let code: &'static GeneticCode = GeneticCode::get(translation_table)?;
        let mut target_fraction = [0.0f64; CODON_COUNT];
        let mut seen = [false; 256];
        for codon in Codon::all() {
            if code.is_stop(codon) {
                continue;
            }
            let aa = code.amino_acid(codon) as usize;
            if seen[aa] {
                continue;
            }
            seen[aa] = true;
            let family = code.synonyms(codon);
            if family.len() <= 1 {
                continue;
            }
            let sum: f64 = family.iter().map(|member| table.weight(*member)).sum();
            if sum == 0.0 {
                return Err(CfError::new(
                    code::HOST_TERM_UNSUPPORTED,
                    format!(
                        "usage.match_frequency has no observed codon for amino acid {}",
                        code.amino_acid(codon) as char
                    ),
                ));
            }
            for member in family {
                target_fraction[member.index_usize()] = table.weight(*member) / sum;
            }
        }
        Ok(Self {
            term: TermId::new("usage.match_frequency")?,
            code,
            target_fraction,
        })
    }

    fn tally(&self, sequence: &CodonSequence) -> [u32; CODON_COUNT] {
        let mut counts = [0u32; CODON_COUNT];
        for codon in sequence.codons() {
            counts[codon.index_usize()] += 1;
        }
        counts
    }

    fn score_counts(&self, counts: &[u32; CODON_COUNT]) -> f64 {
        let mut total = 0.0;
        let mut seen = [false; 256];
        for codon in Codon::all() {
            if self.code.is_stop(codon) {
                continue;
            }
            let aa = self.code.amino_acid(codon) as usize;
            if seen[aa] {
                continue;
            }
            seen[aa] = true;
            let family = self.code.synonyms(codon);
            if family.len() <= 1 {
                continue;
            }
            let n: u32 = family
                .iter()
                .map(|member| counts[member.index_usize()])
                .sum();
            if n == 0 {
                continue;
            }
            let width = f64::from(n);
            for member in family {
                let observed = f64::from(counts[member.index_usize()]) / width;
                let diff = observed - self.target_fraction[member.index_usize()];
                total += diff * diff;
            }
        }
        total
    }

    fn counts_after(&self, sequence: &CodonSequence, edits: &[Edit]) -> Result<FreqState> {
        let mut codons: Vec<Codon> = sequence.codons().to_vec();
        for edit in edits {
            if edit.codon_index >= codons.len() {
                return Err(CfError::new(
                    code::SEQ_INDEX_OUT_OF_RANGE,
                    format!("codon {} is outside the sequence", edit.codon_index),
                ));
            }
            codons[edit.codon_index] = edit.codon;
        }
        let mut counts = [0u32; CODON_COUNT];
        for codon in codons {
            counts[codon.index_usize()] += 1;
        }
        Ok(FreqState { counts })
    }
}

impl TermProvider for MatchFrequencyTerm {
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
        Capability::L1
    }
    /// One substitution changes the fraction of every synonym of that amino
    /// acid, wherever it sits. A block partition would treat those positions
    /// as independent.
    fn interaction_radius(&self) -> InteractionRadius {
        InteractionRadius::Unbounded
    }
    fn provenance(&self) -> Provenance {
        Provenance::builtin("match_frequency", crate::ENGINE_VERSION)
    }
    fn score(&self, ctx: &EvalContext<'_>) -> Result<f64> {
        Ok(self.score_counts(&self.tally(ctx.sequence)))
    }
    fn rebase(&self, ctx: &EvalContext<'_>) -> Result<Rebased> {
        let state = FreqState {
            counts: self.tally(ctx.sequence),
        };
        Ok(Rebased::with_state(self.score_counts(&state.counts), state))
    }
    fn delta(&self, ctx: &EvalContext<'_>, base: &Rebased, edits: &[Edit]) -> Result<f64> {
        let before = base.state::<FreqState>().copied().unwrap_or(FreqState {
            counts: self.tally(ctx.sequence),
        });
        let after = self.counts_after(ctx.sequence, edits)?;
        Ok(self.score_counts(&after.counts) - self.score_counts(&before.counts))
    }
    fn delta_table(
        &self,
        ctx: &EvalContext<'_>,
        base: &Rebased,
        positions: &[usize],
    ) -> Result<Vec<PositionDeltas>> {
        let before = base.state::<FreqState>().copied().unwrap_or(FreqState {
            counts: self.tally(ctx.sequence),
        });
        let base_score = self.score_counts(&before.counts);
        positions
            .iter()
            .map(|&index| {
                let current = ctx.sequence.codons()[index];
                let alternatives = ctx
                    .sequence
                    .synonyms_at(index)
                    .iter()
                    .copied()
                    .filter(|codon| *codon != current)
                    .map(|codon| {
                        let mut counts = before.counts;
                        let current_count = &mut counts[current.index_usize()];
                        *current_count = current_count.saturating_sub(1);
                        counts[codon.index_usize()] += 1;
                        (codon, self.score_counts(&counts) - base_score)
                    })
                    .collect();
                Ok(PositionDeltas {
                    codon_index: index,
                    alternatives,
                })
            })
            .collect()
    }
}
