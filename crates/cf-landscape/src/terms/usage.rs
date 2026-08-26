//! `usage.cai` and `usage.gc_content`.
//!
//! Ticket P1-21. Both are L1: additive or log-additive over codons, so the
//! exact change for a substitution is `O(1)` and needs no rescore.

use cf_landscape_api::{
    Capability, Edit, EvalContext, InteractionRadius, PositionDeltas, Provenance, Rebased,
    TermProvider, TermRole,
};
use cf_types::{code, CfError, Codon, CodonSequence, Direction, Result, TermId, Unit, CODON_COUNT};
use serde_json::Value;

/// Which reference set a codon usage table was built from.
///
/// Carried because CAI is *defined* against a highly-expressed reference, and
/// the host profile schema calls using the genome-wide table by mistake "a
/// classic silent error". Silent is the operative word: the computation
/// succeeds and returns a number in the right range.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum UsageSubset {
    GenomeWide,
    HighlyExpressed,
    UserSupplied,
}

impl UsageSubset {
    pub fn parse(text: &str) -> Option<Self> {
        match text {
            "genome_wide" => Some(UsageSubset::GenomeWide),
            "highly_expressed" => Some(UsageSubset::HighlyExpressed),
            "user_supplied" => Some(UsageSubset::UserSupplied),
            _ => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            UsageSubset::GenomeWide => "genome_wide",
            UsageSubset::HighlyExpressed => "highly_expressed",
            UsageSubset::UserSupplied => "user_supplied",
        }
    }
}

/// Relative adaptiveness `w` per codon, indexed by [`Codon::index`].
#[derive(Clone, Debug)]
pub struct CodonUsageTable {
    weights: [f64; CODON_COUNT],
    subset: UsageSubset,
    translation_table: u8,
}

impl CodonUsageTable {
    /// Read the `codon_usage` block of a `host_profile.schema.json` document.
    ///
    /// Two things are refused rather than worked around.
    ///
    /// A **genome-wide** table is refused for CAI (`CF-HOST-003`). CAI's
    /// reference is the highly-expressed set by definition; against a
    /// genome-wide table the quantity is still computable and no longer means
    /// what its name says.
    ///
    /// An **entirely unobserved synonymous family** is refused
    /// (`CF-HOST-003`). Individual zero weights are fine and in fact useful —
    /// a rare codon scores badly and the search moves off it, which is CAI
    /// working. But normalization sets at least one member of every observed
    /// family to `w = 1`, so an all-zero family means the family was never
    /// seen at all. Every sequence encoding that amino acid then scores
    /// exactly 0 no matter which codon it picks, so the term stops
    /// discriminating precisely where the user most needs it to.
    pub fn from_host_profile(profile: &Value, translation_table: u8) -> Result<Self> {
        let usage = profile.get("codon_usage").ok_or_else(|| {
            CfError::new(
                code::HOST_PROFILE_INVALID,
                "host profile has no codon_usage block",
            )
        })?;
        let subset = usage
            .get("subset")
            .and_then(Value::as_str)
            .and_then(UsageSubset::parse)
            .ok_or_else(|| {
                CfError::new(
                    code::HOST_PROFILE_INVALID,
                    "codon_usage.subset is missing or unrecognised",
                )
            })?;
        let table = usage
            .get("table")
            .and_then(Value::as_object)
            .ok_or_else(|| {
                CfError::new(
                    code::HOST_PROFILE_INVALID,
                    "codon_usage.table is missing or not an object",
                )
            })?;

        let mut weights = [0.0f64; CODON_COUNT];
        for (codon, entry) in table {
            let codon = Codon::parse(codon)?;
            let w = entry
                .get("relative_adaptiveness")
                .and_then(Value::as_f64)
                .ok_or_else(|| {
                    CfError::new(
                        code::HOST_PROFILE_INVALID,
                        format!("codon_usage.table.{codon} has no relative_adaptiveness"),
                    )
                })?;
            if !(0.0..=1.0).contains(&w) {
                return Err(CfError::new(
                    code::HOST_PROFILE_INVALID,
                    format!("relative_adaptiveness for {codon} is {w}, outside [0, 1]"),
                ));
            }
            weights[codon.index_usize()] = w;
        }

        Ok(Self {
            weights,
            subset,
            translation_table,
        })
    }

    pub fn subset(&self) -> UsageSubset {
        self.subset
    }

    pub fn weight(&self, codon: Codon) -> f64 {
        self.weights[codon.index_usize()]
    }

    /// Families in which every codon has `w = 0`.
    fn unobserved_families(&self) -> Result<Vec<u8>> {
        let code = cf_types::GeneticCode::get(self.translation_table)?;
        let mut dead = Vec::new();
        for codon in Codon::all() {
            let aa = code.amino_acid(codon);
            if dead.contains(&aa) {
                continue;
            }
            let family = code.synonyms(codon);
            if family.iter().all(|c| self.weights[c.index_usize()] == 0.0) {
                dead.push(aa);
            }
        }
        Ok(dead)
    }
}

/// Codon Adaptation Index.
///
/// `CAI = (Π wᵢ)^(1/n)` over counted positions. Counted positions exclude stop
/// codons and single-codon families (Met, Trp under the standard code): their
/// `w` is 1 by construction, so including them would leave the product
/// unchanged while inflating `n`, which pulls every CAI toward 1 by an amount
/// that depends on how much methionine the protein happens to contain. That is
/// the standard definition, not a CodonForge convention, and it is stated here
/// because the alternative reading is easy to implement by accident.
///
/// ## Why the state is `(sum_log, zero_count)` rather than a running product
///
/// A product over hundreds of terms underflows; a sum of logs does not. But
/// `log 0 = -∞` swallows the sum, and a zero weight is legitimate — it is how
/// a rare codon is represented. Counting zeros separately keeps both cases
/// exact: `CAI = 0` when any counted position has `w = 0`, and
/// `exp(sum_log / n)` otherwise. The delta then costs `O(1)` in either case,
/// including the case where a move is what creates or removes the zero.
#[derive(Debug)]
pub struct CaiTerm {
    term: TermId,
    table: CodonUsageTable,
}

/// The `(Σ log w, #zeros, n)` state CAI is computed from.
#[derive(Clone, Copy, Debug, PartialEq)]
struct CaiState {
    sum_log: f64,
    zeros: usize,
    counted: usize,
}

impl CaiState {
    fn value(self) -> f64 {
        if self.counted == 0 || self.zeros > 0 {
            0.0
        } else {
            (self.sum_log / self.counted as f64).exp()
        }
    }

    fn add(&mut self, w: f64) {
        self.counted += 1;
        if w == 0.0 {
            self.zeros += 1;
        } else {
            self.sum_log += w.ln();
        }
    }

    fn remove(&mut self, w: f64) {
        self.counted -= 1;
        if w == 0.0 {
            self.zeros -= 1;
        } else {
            self.sum_log -= w.ln();
        }
    }
}

impl CaiTerm {
    pub fn new(table: CodonUsageTable) -> Result<Self> {
        if table.subset == UsageSubset::GenomeWide {
            return Err(CfError::new(
                code::HOST_TERM_UNSUPPORTED,
                "usage.cai needs a highly-expressed reference set; this host profile's \
                 codon_usage.subset is genome_wide. CAI computed against genome-wide usage \
                 is a different quantity that happens to share the name",
            ));
        }
        let dead = table.unobserved_families()?;
        if !dead.is_empty() {
            let names: String = dead.iter().map(|a| *a as char).collect();
            return Err(CfError::new(
                code::HOST_TERM_UNSUPPORTED,
                format!(
                    "usage.cai is uninformative for this host: amino acid(s) {names} have no \
                     observed codon at all, so every sequence containing one scores 0 whatever \
                     it picks. Widen the reference set"
                ),
            ));
        }
        Ok(Self {
            term: TermId::new("usage.cai")?,
            table: table.clone(),
        })
    }

    /// Whether this position enters the geometric mean.
    fn counted(&self, sequence: &CodonSequence, index: usize) -> bool {
        let codon = sequence.codons()[index];
        let code = sequence.table();
        !code.is_stop(codon) && code.synonyms(codon).len() > 1
    }

    fn state(&self, sequence: &CodonSequence) -> CaiState {
        let mut state = CaiState {
            sum_log: 0.0,
            zeros: 0,
            counted: 0,
        };
        for index in 0..sequence.codon_count() {
            if self.counted(sequence, index) {
                state.add(self.table.weight(sequence.codons()[index]));
            }
        }
        state
    }

    /// Exact change from applying `edits` to `sequence`, given the base state.
    ///
    /// `O(|edits|)`, because `before` is pinned rather than recomputed. A
    /// synonymous substitution never changes which positions are counted —
    /// synonyms share an amino acid, so family size and stop-ness are
    /// invariant — which is what makes the pinned state editable in place.
    fn delta_from(
        &self,
        sequence: &CodonSequence,
        before: CaiState,
        edits: &[Edit],
    ) -> Result<f64> {
        let mut after = before;
        for edit in edits {
            if edit.codon_index >= sequence.codon_count() {
                return Err(CfError::new(
                    code::SEQ_INDEX_OUT_OF_RANGE,
                    format!("codon {} is outside the sequence", edit.codon_index),
                ));
            }
            if !self.counted(sequence, edit.codon_index) {
                continue;
            }
            after.remove(self.table.weight(sequence.codons()[edit.codon_index]));
            after.add(self.table.weight(edit.codon));
        }
        Ok(after.value() - before.value())
    }

    /// The pinned state, or a fresh one if the caller handed over someone
    /// else's `Rebased`.
    fn base_state(&self, sequence: &CodonSequence, base: &Rebased) -> CaiState {
        base.state::<CaiState>()
            .copied()
            .unwrap_or_else(|| self.state(sequence))
    }
}

impl TermProvider for CaiTerm {
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
    /// CAI is a geometric mean over independent per-codon weights, so choosing
    /// a codon at one position changes nothing at any other.
    fn interaction_radius(&self) -> InteractionRadius {
        InteractionRadius::Positionwise
    }
    fn provenance(&self) -> Provenance {
        Provenance::builtin("cai_native", crate::ENGINE_VERSION)
    }
    fn score(&self, ctx: &EvalContext<'_>) -> Result<f64> {
        Ok(self.state(ctx.sequence).value())
    }
    fn rebase(&self, ctx: &EvalContext<'_>) -> Result<Rebased> {
        let state = self.state(ctx.sequence);
        Ok(Rebased::with_state(state.value(), state))
    }
    fn delta(&self, ctx: &EvalContext<'_>, base: &Rebased, edits: &[Edit]) -> Result<f64> {
        let state = self.base_state(ctx.sequence, base);
        self.delta_from(ctx.sequence, state, edits)
    }
    fn delta_table(
        &self,
        ctx: &EvalContext<'_>,
        base: &Rebased,
        positions: &[usize],
    ) -> Result<Vec<PositionDeltas>> {
        // Pinned once for the whole table. Recomputing it per alternative was
        // the `O(L·P·S)` defect the `P1-23` audit found, where the capability
        // model promises `O(P·S)`.
        let state = self.base_state(ctx.sequence, base);
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
                        self.delta_from(ctx.sequence, state, &[Edit::new(index, c)])
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

/// GC fraction of the coding sequence.
///
/// ## Why this is capped at `diagnostic`
///
/// GC content has no direction. Too high and the transcript folds and the gene
/// is hard to synthesize; too low and it does likewise for different reasons.
/// The quantity a designer actually wants bounded is `|gc − target|`, and the
/// target is a property of the vendor and the host, not of this term. Declaring
/// `minimize` or `maximize` here would be picking a direction the biology does
/// not have, and the standing rule against fabricated defaults applies to
/// directions as much as to constants. So the raw fraction is reported, and
/// [`TermProvider::max_role`] keeps it out of any ranking; the `Minimize`
/// returned by [`TermProvider::direction`] is unreachable for ordering
/// purposes and exists only because the response schema requires the field.
///
/// A directed GC objective belongs in a separate term with a stated target.
/// Bounding GC as a *constraint* already works today: `Rule::GcWindow`.
#[derive(Debug)]
pub struct GcContentTerm {
    term: TermId,
}

impl Default for GcContentTerm {
    fn default() -> Self {
        Self::new()
    }
}

impl GcContentTerm {
    pub fn new() -> Self {
        Self {
            term: TermId::new("usage.gc_content").expect("literal term id"),
        }
    }

    fn gc_of(codon: Codon) -> u32 {
        codon
            .to_string_upper()
            .bytes()
            .filter(|b| *b == b'G' || *b == b'C')
            .count() as u32
    }

    fn counts(sequence: &CodonSequence) -> (u32, u32) {
        let gc: u32 = sequence.codons().iter().copied().map(Self::gc_of).sum();
        (gc, (sequence.codon_count() * 3) as u32)
    }
}

impl TermProvider for GcContentTerm {
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
    fn max_role(&self) -> TermRole {
        TermRole::Diagnostic
    }
    /// A fraction of independent per-base counts.
    fn interaction_radius(&self) -> InteractionRadius {
        InteractionRadius::Positionwise
    }
    fn provenance(&self) -> Provenance {
        Provenance::builtin("gc_native", crate::ENGINE_VERSION)
    }
    fn score(&self, ctx: &EvalContext<'_>) -> Result<f64> {
        let (gc, total) = Self::counts(ctx.sequence);
        Ok(if total == 0 {
            0.0
        } else {
            f64::from(gc) / f64::from(total)
        })
    }
    fn delta(&self, ctx: &EvalContext<'_>, _base: &Rebased, edits: &[Edit]) -> Result<f64> {
        // `codon_count * 3`, not a pass over the sequence. Calling `counts()`
        // here summed GC over every codon purely to obtain a denominator that
        // is known in `O(1)`.
        let total = (ctx.sequence.codon_count() * 3) as u32;
        if total == 0 {
            return Ok(0.0);
        }
        let codons = ctx.sequence.codons();
        // The GC change is an exact integer, so it is accumulated as one and
        // divided once. Computing `new/total - old/total` instead would round
        // twice and land one ULP away from the value `delta_table` reports for
        // the same move — two spellings of one quantity that disagree, which is
        // how a "harmless" floating-point difference becomes a delta test that
        // fails on a Tuesday.
        let mut change = 0i64;
        for edit in edits {
            if edit.codon_index >= codons.len() {
                return Err(CfError::new(
                    code::SEQ_INDEX_OUT_OF_RANGE,
                    format!("codon {} is outside the sequence", edit.codon_index),
                ));
            }
            change += i64::from(Self::gc_of(edit.codon))
                - i64::from(Self::gc_of(codons[edit.codon_index]));
        }
        Ok(change as f64 / f64::from(total))
    }
    fn delta_table(
        &self,
        ctx: &EvalContext<'_>,
        _base: &Rebased,
        positions: &[usize],
    ) -> Result<Vec<PositionDeltas>> {
        let total = (ctx.sequence.codon_count() * 3) as f64;
        positions
            .iter()
            .map(|&index| {
                let current = ctx.sequence.codons()[index];
                Ok(PositionDeltas {
                    codon_index: index,
                    alternatives: ctx
                        .sequence
                        .synonyms_at(index)
                        .iter()
                        .copied()
                        .filter(|c| *c != current)
                        .map(|c| {
                            (
                                c,
                                (f64::from(Self::gc_of(c)) - f64::from(Self::gc_of(current)))
                                    / total,
                            )
                        })
                        .collect(),
                })
            })
            .collect()
    }
}
