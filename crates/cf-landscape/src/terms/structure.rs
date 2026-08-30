//! The windowed-accessibility term's contract.
//!
//! Ticket P2-11. `AMD-013`, `AMD-001`, `DEP-STRU-003`, `DEP-STRU-004`.
//!
//! ## A provider that declares but does not compute
//!
//! Every other built-in term computes its own value in Rust. This one cannot:
//! `AMD-013` put the folding engine out of process, so the numbers arrive from
//! a Python scorer through `P2-05a`'s batch path, once per generation.
//!
//! What is left for a provider to do is still worth doing, and is the reason
//! this type exists rather than the term being handled entirely on the Python
//! side. The registry is where a term's *contract* lives — its capability, its
//! role, its direction, its interaction radius — and every consumer reads that
//! contract from one place. A term whose contract lived only in the
//! orchestrator would need a second copy of each of those facts, which is the
//! duplicated-rule pattern this project has now found five times.
//!
//! So: [`score`](TermProvider::score) refuses. Not because scoring is
//! unimplemented, but because scoring in-core is the thing `AMD-001` and
//! `ADR-0001` between them forbid — a per-move fold would cross the boundary
//! once per candidate, roughly 19 200 times for a 64/300 run against a budget
//! of 300. `CF-STRU-006` says so by name rather than returning a plausible
//! number computed some cheaper way.

use cf_landscape_api::{
    Capability, Edit, EvalContext, InteractionRadius, Provenance, Rebased, TermProvider, TermRole,
};
use cf_types::{code, CfError, Direction, Result, TermId, Unit};

/// `structure.local_accessibility`.
#[derive(Debug)]
pub struct LocalAccessibilityTerm {
    term: TermId,
    /// The folding span, in nucleotides. Sets the interaction radius, so it is
    /// held rather than defaulted at each call.
    span_nt: usize,
}

impl LocalAccessibilityTerm {
    /// `span_nt` is the maximum base-pair span the engine is configured with.
    ///
    /// Taken as a parameter because it decides the interaction radius, and a
    /// radius guessed here would disagree with the engine that produced the
    /// numbers — tier 2's `B >= R` partition condition is checked against it,
    /// so a wrong value is not a cosmetic mismatch but a silently invalid
    /// block decomposition.
    pub fn new(span_nt: usize) -> Result<Self> {
        if span_nt == 0 {
            return Err(CfError::new(
                code::STRU_BAD_WINDOW,
                "a folding span of zero nucleotides admits no base pair, so the term \
                 would report a constant and rank nothing",
            ));
        }
        Ok(Self {
            term: TermId::new("structure.local_accessibility")?,
            span_nt,
        })
    }

    pub fn span_nt(&self) -> usize {
        self.span_nt
    }
}

impl TermProvider for LocalAccessibilityTerm {
    fn term_id(&self) -> &TermId {
        &self.term
    }

    fn unit(&self) -> Unit {
        Unit::new(Unit::DIMENSIONLESS)
    }

    /// Minimize. The value is a mean negative log accessibility, so a strongly
    /// paired region costs more than an open one.
    fn direction(&self) -> Direction {
        Direction::Minimize
    }

    /// `L0`. Scalar, whole-sequence, out of core.
    ///
    /// `DEP-STRU-003` corrected this from the `L1` the plan declared: an exact
    /// per-move delta was reachable only for the in-process kernel `AMD-013`
    /// removed.
    fn capability(&self) -> Capability {
        Capability::L0
    }

    /// A real Pareto axis, unlike the diagnostic form it shipped as.
    fn max_role(&self) -> TermRole {
        TermRole::ParetoAxis
    }

    /// A base pair can only form within the folding span, so an edit cannot
    /// change accessibility beyond it. This is the same locality claim the
    /// incremental cache rests on, and `TST-STRU-006` tests it against the real
    /// engine rather than assuming it.
    fn interaction_radius(&self) -> InteractionRadius {
        InteractionRadius::Codons(self.span_nt.div_ceil(3))
    }

    fn provenance(&self) -> Provenance {
        Provenance::builtin("local_accessibility", crate::ENGINE_VERSION)
    }

    /// Refuses, by design.
    ///
    /// The temptation this guards against is specific: a cheaper in-core proxy
    /// — GC content in a window, say — would return a number for every
    /// candidate, every downstream check would pass, and the objective would be
    /// measuring something other than what it claims. Refusing makes the
    /// out-of-core path mandatory rather than preferred.
    fn score(&self, _ctx: &EvalContext<'_>) -> Result<f64> {
        Err(CfError::new(
            code::STRU_OUT_OF_CORE,
            "structure.local_accessibility cannot be scored in core: the folding \
             engine is out of process (AMD-013), so values arrive through the \
             batch path once per generation (P2-05a). A cheaper in-core proxy \
             would return a number for every candidate and quietly be measuring \
             something else",
        ))
    }

    fn delta(&self, _ctx: &EvalContext<'_>, _base: &Rebased, _edits: &[Edit]) -> Result<f64> {
        Err(CfError::new(
            code::STRU_OUT_OF_CORE,
            "structure.local_accessibility declares L0, so it has no exact delta; \
             a caller reaching here has ignored the declared capability",
        ))
    }
}
