//! Test providers shared by the contract tests.
//!
//! Every provider here is deliberately trivial to score, so that a failure is
//! a failure of the contract rather than of arithmetic. Several are dishonest
//! on purpose — they are the fixtures the negative controls need, and building
//! them as first-class objects is cheaper than editing a real provider to
//! misbehave and remembering to edit it back.

use std::sync::atomic::{AtomicU64, Ordering};

use cf_landscape_api::{
    Capability, EvalContext, InteractionRadius, PositionDeltas, Provenance, Rebased, TermProvider,
    TermRole, Uncertainty, UncertaintyMethod,
};
use cf_types::{code, CfError, Direction, Result, TermId, Unit};

/// Per-codon weight used by the additive test terms.
///
/// Additive so the exact delta has a closed form a reader can check by eye:
/// swapping codon `a` for codon `b` changes the score by `w(b) - w(a)`.
pub fn weight(codon: cf_types::Codon) -> f64 {
    let bytes = codon.to_string_upper().into_bytes();
    bytes.iter().map(|b| f64::from(*b - b'A')).sum::<f64>() / 100.0
}

fn additive_score(sequence: &cf_types::CodonSequence) -> f64 {
    sequence.codons().iter().copied().map(weight).sum()
}

/// An honest L1 provider: additive, with exact deltas.
#[derive(Debug)]
pub struct AdditiveTerm {
    term: TermId,
    role_cap: TermRole,
}

impl AdditiveTerm {
    pub fn new(term: &str) -> Self {
        Self {
            term: TermId::new(term).unwrap(),
            role_cap: TermRole::ParetoAxis,
        }
    }

    /// Cap the provider at diagnostic, the way a kinetics term is capped
    /// until `EXP-06`.
    pub fn diagnostic_only(mut self) -> Self {
        self.role_cap = TermRole::Diagnostic;
        self
    }
}

impl TermProvider for AdditiveTerm {
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
    fn interaction_radius(&self) -> InteractionRadius {
        InteractionRadius::Positionwise
    }
    fn max_role(&self) -> TermRole {
        self.role_cap
    }
    fn provenance(&self) -> Provenance {
        Provenance::builtin("test_additive", "0.1.0")
    }
    fn score(&self, ctx: &EvalContext<'_>) -> Result<f64> {
        Ok(additive_score(ctx.sequence))
    }
    fn delta(
        &self,
        ctx: &EvalContext<'_>,
        _base: &Rebased,
        edits: &[cf_landscape_api::Edit],
    ) -> Result<f64> {
        let codons = ctx.sequence.codons();
        Ok(edits
            .iter()
            .map(|e| weight(e.codon) - weight(codons[e.codon_index]))
            .sum())
    }
    fn delta_table(
        &self,
        ctx: &EvalContext<'_>,
        _base: &Rebased,
        positions: &[usize],
    ) -> Result<Vec<PositionDeltas>> {
        Ok(positions
            .iter()
            .map(|&i| {
                let current = ctx.sequence.codons()[i];
                PositionDeltas {
                    codon_index: i,
                    alternatives: ctx
                        .sequence
                        .synonyms_at(i)
                        .iter()
                        .copied()
                        .filter(|c| *c != current)
                        .map(|c| (c, weight(c) - weight(current)))
                        .collect(),
                }
            })
            .collect())
    }
}

/// An honest L0 provider: scores, and refuses delta queries.
///
/// The point of this fixture is that L0 is *legal*. Composition has to have a
/// correct path for it, and that path has to cost only this term.
#[derive(Debug)]
pub struct ScoreOnlyTerm {
    term: TermId,
    pub calls: AtomicU64,
}

impl ScoreOnlyTerm {
    pub fn new(term: &str) -> Self {
        Self {
            term: TermId::new(term).unwrap(),
            calls: AtomicU64::new(0),
        }
    }
}

impl TermProvider for ScoreOnlyTerm {
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
        Capability::L0
    }
    fn interaction_radius(&self) -> InteractionRadius {
        InteractionRadius::Positionwise
    }
    fn provenance(&self) -> Provenance {
        Provenance::builtin("test_score_only", "0.1.0")
    }
    fn score(&self, ctx: &EvalContext<'_>) -> Result<f64> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        Ok(additive_score(ctx.sequence) * 2.0)
    }
}

/// Declares L1 and returns a delta that is not the change.
///
/// The fixture for negative control (a): a term that silently claims delta
/// support it does not have. It returns `0.0`, which is the value a plausible
/// stub returns and is indistinguishable from "no change" at a glance.
#[derive(Debug)]
pub struct LyingL1Term {
    term: TermId,
}

impl LyingL1Term {
    pub fn new(term: &str) -> Self {
        Self {
            term: TermId::new(term).unwrap(),
        }
    }
}

impl TermProvider for LyingL1Term {
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
    fn interaction_radius(&self) -> InteractionRadius {
        InteractionRadius::Positionwise
    }
    fn provenance(&self) -> Provenance {
        Provenance::builtin("test_lying", "0.1.0")
    }
    fn score(&self, ctx: &EvalContext<'_>) -> Result<f64> {
        Ok(additive_score(ctx.sequence))
    }
    fn delta(
        &self,
        _ctx: &EvalContext<'_>,
        _base: &Rebased,
        _edits: &[cf_landscape_api::Edit],
    ) -> Result<f64> {
        Ok(0.0)
    }
    fn delta_table(
        &self,
        ctx: &EvalContext<'_>,
        _base: &Rebased,
        positions: &[usize],
    ) -> Result<Vec<PositionDeltas>> {
        Ok(positions
            .iter()
            .map(|&i| PositionDeltas {
                codon_index: i,
                alternatives: ctx
                    .sequence
                    .synonyms_at(i)
                    .iter()
                    .copied()
                    .map(|c| (c, 0.0))
                    .collect(),
            })
            .collect())
    }
}

/// Declares L1 but has no delta table at all.
///
/// Distinct from [`LyingL1Term`]: this one forgot to implement the method
/// rather than implementing it wrongly, which is the more common way a
/// capability claim goes stale.
#[derive(Debug)]
pub struct MissingTableTerm {
    term: TermId,
}

impl MissingTableTerm {
    pub fn new(term: &str) -> Self {
        Self {
            term: TermId::new(term).unwrap(),
        }
    }
}

impl TermProvider for MissingTableTerm {
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
    fn interaction_radius(&self) -> InteractionRadius {
        InteractionRadius::Positionwise
    }
    fn provenance(&self) -> Provenance {
        Provenance::builtin("test_missing_table", "0.1.0")
    }
    fn score(&self, ctx: &EvalContext<'_>) -> Result<f64> {
        Ok(additive_score(ctx.sequence))
    }
    fn delta(
        &self,
        ctx: &EvalContext<'_>,
        _base: &Rebased,
        edits: &[cf_landscape_api::Edit],
    ) -> Result<f64> {
        let codons = ctx.sequence.codons();
        Ok(edits
            .iter()
            .map(|e| weight(e.codon) - weight(codons[e.codon_index]))
            .sum())
    }
}

/// A sampled provider whose value depends on the seed.
///
/// Uses a splitmix64 step so the dependence is real but reproducible; the test
/// cares that the same seed gives the same number, not what the number is.
#[derive(Debug)]
pub struct SampledTerm {
    term: TermId,
}

impl SampledTerm {
    pub fn new(term: &str) -> Self {
        Self {
            term: TermId::new(term).unwrap(),
        }
    }
}

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

impl TermProvider for SampledTerm {
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
        Capability::L0
    }
    fn interaction_radius(&self) -> InteractionRadius {
        InteractionRadius::Positionwise
    }
    fn uncertainty_method(&self) -> UncertaintyMethod {
        UncertaintyMethod::McDropout
    }
    fn provenance(&self) -> Provenance {
        Provenance::builtin("test_sampled", "0.1.0")
    }
    fn score(&self, ctx: &EvalContext<'_>) -> Result<f64> {
        let noise = (splitmix64(ctx.seed) >> 11) as f64 / (1u64 << 53) as f64;
        Ok(additive_score(ctx.sequence) + noise)
    }
    fn uncertainty(&self, _ctx: &EvalContext<'_>, _value: f64) -> Uncertainty {
        Uncertainty {
            method: UncertaintyMethod::McDropout,
            calibrated: true,
            sigma: Some(0.05),
            interval: None,
            coverage: Some(0.9),
            n_members: Some(8),
            // Never measured. The term is therefore ineligible for
            // probabilistic dominance, and says so instead of guessing.
            shared_fraction: None,
        }
    }
}

/// A provider whose value depends on how many times it has been called.
///
/// Nothing legitimate behaves this way. It exists to give the term-independence
/// assertion something to catch: without a fixture that violates the property,
/// the assertion would pass on every provider and prove nothing.
#[derive(Debug)]
pub struct OrderSensitiveTerm {
    term: TermId,
    calls: AtomicU64,
}

impl OrderSensitiveTerm {
    pub fn new(term: &str) -> Self {
        Self {
            term: TermId::new(term).unwrap(),
            calls: AtomicU64::new(0),
        }
    }
}

impl TermProvider for OrderSensitiveTerm {
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
        Capability::L0
    }
    fn interaction_radius(&self) -> InteractionRadius {
        InteractionRadius::Positionwise
    }
    fn provenance(&self) -> Provenance {
        Provenance::builtin("test_order_sensitive", "0.1.0")
    }
    fn score(&self, _ctx: &EvalContext<'_>) -> Result<f64> {
        Ok(self.calls.fetch_add(1, Ordering::SeqCst) as f64)
    }
}

/// A provider that always fails, for exercising the error path.
#[derive(Debug)]
pub struct FailingTerm {
    term: TermId,
}

impl FailingTerm {
    pub fn new(term: &str) -> Self {
        Self {
            term: TermId::new(term).unwrap(),
        }
    }
}

impl TermProvider for FailingTerm {
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
        Capability::L0
    }
    fn interaction_radius(&self) -> InteractionRadius {
        InteractionRadius::Positionwise
    }
    fn provenance(&self) -> Provenance {
        Provenance::builtin("test_failing", "0.1.0")
    }
    fn score(&self, _ctx: &EvalContext<'_>) -> Result<f64> {
        Err(CfError::new(
            code::BUDGET_EXCEEDED,
            "test provider always exceeds its budget",
        ))
    }
}
