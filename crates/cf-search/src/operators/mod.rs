//! The variation operators and the mix that selects among them.
//!
//! Ticket P1-23.

pub mod crossover;
pub mod mutation;

use cf_constraints::{ActionableViolation, TabuLedger};
use cf_types::{code, CfError, Result};

use crate::context::{CodonWeights, ContextModel};
use crate::genotype::{Genotype, Move, MoveKind};
use crate::rng::SearchRng;

pub use crossover::{BlockUniform, CrossoverOperator, OnePoint, TwoPoint};
pub use mutation::{BlockGibbs, PottsGuided, UniformSynonymous, ViolationDirected};

/// Everything an operator may look at.
#[derive(Clone, Copy)]
pub struct MutationContext<'a> {
    pub genotype: &'a Genotype,
    pub weights: &'a CodonWeights,
    /// Actionable violations of this individual. Empty for a feasible one.
    pub violations: &'a [ActionableViolation<'a>],
    /// The tabu ledger, when a generation loop is driving.
    ///
    /// `None` in `P1-23`, because `TABU_TICK` belongs to the loop and the loop
    /// is `P1-18b`. The consultation is implemented and tested here against a
    /// hand-advanced ledger; only the clock is deferred.
    pub tabu: Option<&'a TabuLedger>,
    /// Which coding region this genotype is, for tabu bookkeeping. Zero for a
    /// single-cistron construct.
    pub region: usize,
}

impl<'a> MutationContext<'a> {
    pub fn new(genotype: &'a Genotype, weights: &'a CodonWeights) -> Self {
        const NO_VIOLATIONS: &[ActionableViolation<'static>] = &[];
        Self {
            genotype,
            weights,
            violations: NO_VIOLATIONS,
            tabu: None,
            region: 0,
        }
    }

    pub fn with_violations(mut self, violations: &'a [ActionableViolation<'a>]) -> Self {
        self.violations = violations;
        self
    }

    pub fn with_tabu(mut self, tabu: &'a TabuLedger) -> Self {
        self.tabu = Some(tabu);
        self
    }

    pub fn with_region(mut self, region: usize) -> Self {
        self.region = region;
        self
    }
}

impl std::fmt::Debug for MutationContext<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MutationContext")
            .field("codons", &self.genotype.len())
            .field("violations", &self.violations.len())
            .field("tabu", &self.tabu.is_some())
            .field("region", &self.region)
            .finish()
    }
}

/// One mutation operator.
pub trait MutationOperator: Send + Sync + std::fmt::Debug {
    fn name(&self) -> &'static str;
    fn kind(&self) -> MoveKind;

    /// Whether this operator has anything to offer for this individual.
    ///
    /// `violation_directed` is the reason this exists: §5.1 scopes it to
    /// infeasible individuals, and selecting it for a feasible one would return
    /// an empty move and waste the draw.
    fn applicable(&self, ctx: &MutationContext<'_>) -> bool {
        let _ = ctx;
        true
    }

    fn propose(&self, ctx: &MutationContext<'_>, rng: &mut SearchRng) -> Result<Move>;
}

/// An operator the plan specifies but this build cannot provide.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnavailableOperator {
    pub name: &'static str,
    pub reason: &'static str,
    pub ticket: &'static str,
}

impl UnavailableOperator {
    /// The warning a job should surface.
    ///
    /// Returned as a `CfError` with a registered code rather than as a string,
    /// so the orchestrator can put it straight into `job_status.warnings`
    /// without inventing a code at the boundary.
    pub fn warning(&self) -> CfError {
        CfError::new(
            code::SRCH_OPERATOR_UNAVAILABLE,
            format!(
                "mutation operator {} is not available in this build ({}); it arrives with {}. \
                 The remaining operators' weights are renormalized, which changes how the \
                 search explores",
                self.name, self.reason, self.ticket
            ),
        )
    }
}

/// The weighted set of mutation operators.
///
/// Weights are renormalized over the operators that are both *registered* and
/// *applicable* to the individual at hand. Both filters matter: a build without
/// a context model has two of four operators, and a feasible individual has no
/// use for `violation_directed`.
#[derive(Debug)]
pub struct OperatorMix {
    operators: Vec<(Box<dyn MutationOperator>, f64)>,
    unavailable: Vec<UnavailableOperator>,
}

impl OperatorMix {
    /// Default weights from §5.1.
    pub const W_UNIFORM: f64 = 0.35;
    pub const W_POTTS: f64 = 0.35;
    pub const W_BLOCK_GIBBS: f64 = 0.20;
    pub const W_VIOLATION: f64 = 0.10;

    pub fn empty() -> Self {
        Self {
            operators: Vec::new(),
            unavailable: Vec::new(),
        }
    }

    /// The operators available without an epistatic surface.
    ///
    /// `potts_guided` and `block_gibbs` are recorded as unavailable rather than
    /// approximated. Backing them with host codon usage would turn
    /// `potts_guided` into `uniform_synonymous` with different weights under a
    /// name that says otherwise, and the resulting change in exploration
    /// behaviour would be attributed to anything but the operator that caused
    /// it. See `DEP-SRCH-005`.
    pub fn without_context_model() -> Self {
        let mut mix = Self::empty();
        mix.register(Box::new(UniformSynonymous), Self::W_UNIFORM);
        mix.register(Box::new(ViolationDirected), Self::W_VIOLATION);
        mix.unavailable.push(UnavailableOperator {
            name: "potts_guided",
            reason: "no epistatic context model is registered",
            ticket: "P3-05",
        });
        mix.unavailable.push(UnavailableOperator {
            name: "block_gibbs",
            reason: "no epistatic context model is registered",
            ticket: "P3-05",
        });
        mix
    }

    /// The full four-operator mix of §5.1.
    pub fn with_context_model(model: std::sync::Arc<dyn ContextModel>) -> Self {
        let mut mix = Self::empty();
        mix.register(Box::new(UniformSynonymous), Self::W_UNIFORM);
        mix.register(Box::new(PottsGuided::new(model.clone())), Self::W_POTTS);
        mix.register(Box::new(BlockGibbs::new(model)), Self::W_BLOCK_GIBBS);
        mix.register(Box::new(ViolationDirected), Self::W_VIOLATION);
        mix
    }

    pub fn register(&mut self, operator: Box<dyn MutationOperator>, weight: f64) {
        self.operators.push((operator, weight.max(0.0)));
    }

    pub fn len(&self) -> usize {
        self.operators.len()
    }

    pub fn is_empty(&self) -> bool {
        self.operators.is_empty()
    }

    pub fn names(&self) -> Vec<&'static str> {
        self.operators.iter().map(|(op, _)| op.name()).collect()
    }

    /// Operators the plan specifies that this build does not have.
    pub fn unavailable(&self) -> &[UnavailableOperator] {
        &self.unavailable
    }

    /// Warnings a job must surface before running.
    ///
    /// Nothing in `P1-23` can surface them — there is no job. `P1-29` owns the
    /// wiring into `job_status.warnings`; this method is what it will call.
    pub fn warnings(&self) -> Vec<CfError> {
        self.unavailable.iter().map(|u| u.warning()).collect()
    }

    /// Propose a move, choosing an operator by renormalized weight.
    pub fn propose(&self, ctx: &MutationContext<'_>, rng: &mut SearchRng) -> Result<Move> {
        if self.operators.is_empty() {
            return Err(CfError::new(
                code::SRCH_NO_OPERATOR,
                "no mutation operator is registered",
            ));
        }
        let weights: Vec<f64> = self
            .operators
            .iter()
            .map(|(op, w)| if op.applicable(ctx) { *w } else { 0.0 })
            .collect();

        let Some(choice) = rng.weighted(&weights) else {
            // Every registered operator declined. Real for a feasible
            // individual in a mix that is only `violation_directed`, and a
            // silent empty move would look like a mutation that happened to
            // change nothing.
            return Err(CfError::new(
                code::SRCH_NO_OPERATOR,
                "no registered mutation operator applies to this individual",
            ));
        };
        self.operators[choice].0.propose(ctx, rng)
    }
}
