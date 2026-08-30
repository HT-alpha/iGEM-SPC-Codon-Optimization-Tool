//! Surrogate rationing for L0 terms.
//!
//! Ticket P2-06. Spec: `planning/03-algorithms/03-nsga-loop-and-operators.md` §6.
//!
//! L1 terms are nearly free and computed exactly for every offspring. The L0
//! terms are the entire cost of tier 3, so only `E_max` individuals per
//! generation get a real evaluation and the rest carry a surrogate estimate —
//! or, if they are not chosen at all, `UNKNOWN`, which is *not* imputed.
//!
//! ## Two disciplines encoded in types rather than documented
//!
//! Following the `CF-LAND-011` precedent, where "one crossing per generation"
//! became a property of the batch type instead of a rule callers were asked to
//! remember:
//!
//! 1. **The surrogate cannot train on itself.** [`Observation`] is constructible
//!    only from an exact evaluation, so [`Surrogate::refit`] structurally cannot
//!    ingest its own predictions. §6 calls this out as the same failure mode as
//!    `AMD-004`'s feedback loop at smaller scale, and a rule that lives in a
//!    comment is a rule that gets violated by the next person adding a feature.
//! 2. **Truncation defers, it does not discard.** [`Rationer::choose`] returns
//!    the overflow as [`Rationed::deferred`], and the type has no method that
//!    drops it. An individual squeezed out by `E_max` is first in line next
//!    generation.
//!
//! ## Why the deferral is also a warning
//!
//! Modelling deferral correctly is not the same as making it visible. An
//! individual that stays `UNKNOWN` to the end because the budget never reached
//! it is structurally the same problem `CF-ORCH-006` exists to solve for skipped
//! tiers: the run is quietly narrower than the plan describes, and nothing in
//! the output says so. So [`Rationer::unevaluated_at_end`] reports it, and the
//! orchestrator raises `CF-SRCH-016` with stage and subject provenance.

use std::collections::BTreeSet;

use cf_metrics::{fast_nondominated_sort, Ranked};
use cf_types::{code, CfError, Result};

use crate::rng::SearchRng;

/// An exact L0 evaluation. The only thing the surrogate may learn from.
///
/// Deliberately has no public constructor taking a bare `f64`: the point of the
/// type is that a surrogate prediction cannot be laundered into training data by
/// a caller who has one to hand.
#[derive(Clone, Debug, PartialEq)]
pub struct Observation {
    features: Vec<f64>,
    value: f64,
}

impl Observation {
    /// Record an exact evaluation.
    ///
    /// `exact` is the value that came back from a real L0 call. Callers holding
    /// a surrogate estimate have nothing to pass here, which is the intent.
    pub fn from_exact(features: Vec<f64>, exact: f64) -> Self {
        Self {
            features,
            value: exact,
        }
    }

    pub fn value(&self) -> f64 {
        self.value
    }
}

/// A ridge regression on k-mer and L1-term features.
///
/// Ridge rather than anything richer because §6 budgets it at microseconds and
/// because a surrogate that is expensive to fit defeats the purpose of having
/// one. The penalty is what keeps it stable when the training set is small and
/// its features are collinear, which is the normal case early in a run.
#[derive(Clone, Debug)]
pub struct Surrogate {
    weights: Vec<f64>,
    intercept: f64,
    lambda: f64,
    trained_on: usize,
}

impl Surrogate {
    pub fn new(lambda: f64) -> Self {
        Self {
            weights: Vec::new(),
            intercept: 0.0,
            lambda,
            trained_on: 0,
        }
    }

    /// Whether the surrogate has been fit at all.
    ///
    /// Before the first fit it has no opinion, and callers must treat every
    /// individual as unranked by it rather than reading a zero as a prediction.
    pub fn is_trained(&self) -> bool {
        self.trained_on > 0
    }

    pub fn trained_on(&self) -> usize {
        self.trained_on
    }

    /// Fit on exact observations only — enforced by [`Observation`]'s
    /// constructor, not by asking the caller to be careful.
    pub fn refit(&mut self, observations: &[Observation]) -> Result<()> {
        if observations.is_empty() {
            return Err(CfError::new(
                code::SRCH_SURROGATE_UNFIT,
                "cannot fit a surrogate on no observations",
            ));
        }
        let width = observations[0].features.len();
        if observations.iter().any(|o| o.features.len() != width) {
            return Err(CfError::new(
                code::SRCH_SURROGATE_UNFIT,
                "observations disagree on feature width",
            ));
        }

        let n = observations.len() as f64;
        let mean_y = observations.iter().map(|o| o.value).sum::<f64>() / n;
        let means: Vec<f64> = (0..width)
            .map(|j| observations.iter().map(|o| o.features[j]).sum::<f64>() / n)
            .collect();

        // Ridge with a diagonal approximation: the full normal equations would
        // need a solve per refit, and §6's budget for this is microseconds. The
        // diagonal form is exact when features are uncorrelated and shrinks
        // sensibly when they are not, which is the behaviour the penalty is
        // there for anyway.
        let mut weights = vec![0.0; width];
        for (j, weight) in weights.iter_mut().enumerate() {
            let mut covariance = 0.0;
            let mut variance = 0.0;
            for observation in observations {
                let centred = observation.features[j] - means[j];
                covariance += centred * (observation.value - mean_y);
                variance += centred * centred;
            }
            *weight = if variance + self.lambda > 0.0 {
                covariance / (variance + self.lambda)
            } else {
                0.0
            };
        }

        self.intercept = mean_y - weights.iter().zip(&means).map(|(w, m)| w * m).sum::<f64>();
        self.weights = weights;
        self.trained_on = observations.len();
        Ok(())
    }

    /// Predict, or refuse if never fit.
    ///
    /// Refusing rather than returning zero: an untrained surrogate returning a
    /// number is indistinguishable from a trained one that thinks the answer is
    /// zero, and the ranking that follows would be silently arbitrary.
    pub fn predict(&self, features: &[f64]) -> Result<f64> {
        if !self.is_trained() {
            return Err(CfError::new(
                code::SRCH_SURROGATE_UNFIT,
                "the surrogate has not been fit; an untrained prediction is \
                 indistinguishable from a trained one that happens to be zero",
            ));
        }
        if features.len() != self.weights.len() {
            return Err(CfError::new(
                code::SRCH_SURROGATE_UNFIT,
                format!(
                    "surrogate fit on {} features, asked to predict from {}",
                    self.weights.len(),
                    features.len()
                ),
            ));
        }
        Ok(self.intercept
            + features
                .iter()
                .zip(&self.weights)
                .map(|(f, w)| f * w)
                .sum::<f64>())
    }
}

/// A generation's ranked candidates for exact evaluation.
///
/// A *ranking*, not a selection. `DEP-SRCH-011`: the budget cut, the deferral
/// queue and the "defers rather than discards" guarantee all live in
/// [`cf_landscape_api::batch::BatchProtocol`], which already implemented them
/// before this module existed. Two mechanisms enforcing one rule is one too
/// many, and the one that survives is the one the worker crossing already has
/// to go through.
#[derive(Clone, Debug, PartialEq)]
pub struct Rationed {
    /// Candidate indices in priority order, longest-first. The protocol takes
    /// as many as its budget allows and defers the tail.
    pub ranked: Vec<usize>,
    /// Why each index is where it is, for the report.
    pub reasons: Vec<(usize, Reason)>,
}

/// Which of §6's four components put an individual in the chosen set.
///
/// **Declaration order is truncation priority, and that is the only place the
/// priority is written down.** `DEP-SRCH-010`. Two mechanisms enforcing one
/// rule is one too many, so there is no separate priority table to drift out of
/// step with this enum's `Ord`.
///
/// The order is not §6's listing order, deliberately:
///
/// * `Deferred` first, because deferral is a promise. A promise that competes
///   for room with fresh candidates every generation is not a promise, and
///   under a budget that is always full it would never once be kept.
/// * `RandomSample` second, because §6 calls it "not optional". A plain
///   truncation of the union drops it exactly when the front is large — which
///   is precisely the case where the surrogate is training only on points it
///   already likes, and so precisely the case the sample exists for.
/// * `Nondominated` and `ExpectedGain` last, in that order: they are the
///   components that have the whole rest of the population to draw from, so
///   they are the ones that can afford to lose a slot.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Reason {
    /// (d) deferred from last generation.
    Deferred,
    /// (b) the mandatory uniform sample.
    RandomSample,
    /// (a) non-dominated under L1 ∪ surrogate.
    Nondominated,
    /// (c) top-k by expected hypervolume gain.
    ExpectedGain,
}

/// The rationing state that persists across generations.
#[derive(Clone, Debug)]
pub struct Rationer {
    e_max: usize,
    refit_period: u64,
    sample_fraction: f64,
    observations: Vec<Observation>,
}

impl Rationer {
    pub fn new(e_max: usize, refit_period: u64) -> Self {
        Self {
            e_max,
            refit_period,
            sample_fraction: 0.10,
            observations: Vec::new(),
        }
    }

    pub fn e_max(&self) -> usize {
        self.e_max
    }

    pub fn should_refit(&self, generation: u64) -> bool {
        self.refit_period > 0 && generation > 0 && generation % self.refit_period == 0
    }

    /// Record an exact result. Only exact results reach the training set.
    pub fn absorb(&mut self, observation: Observation) {
        self.observations.push(observation);
    }

    pub fn observations(&self) -> &[Observation] {
        &self.observations
    }

    /// Build this generation's chosen set: §6's four components, then truncate.
    ///
    /// `objectives` holds each individual's L1 vector already extended with the
    /// surrogate's L0 estimate where one is available. `rng` supplies component
    /// (b) only — the other three components are deterministic given the
    /// objectives, so a seed change perturbs the sample and nothing else.
    pub fn rank(
        &mut self,
        objectives: &[Vec<f64>],
        violations: &[f64],
        rng: &mut SearchRng,
    ) -> Rationed {
        let population = objectives.len();
        let mut reasons: Vec<(usize, Reason)> = Vec::new();
        let mut picked: BTreeSet<usize> = BTreeSet::new();

        // Component (d), "deferred from last generation", is deliberately
        // absent here. `BatchProtocol::step` puts its own deferred queue at the
        // head of whatever ranking it is handed, so re-adding it would either
        // duplicate the rule or quietly disagree with it.

        // (b) the mandatory uniform sample, drawn from the whole population and
        // drawn *before* (a).
        //
        // Order matters here for a reason that is easy to miss: if (a) runs
        // first and the population is largely non-dominated, every index the
        // sampler draws has already been claimed, so the sample records nothing
        // and truncation then keeps whichever non-dominated individuals happen
        // to sort first. The surrogate ends up trained on a deterministic slice
        // of the front — the exact drift §6 says this component exists to
        // prevent — while the code looks like it drew a sample.
        let sample_size = ((population as f64) * self.sample_fraction).round() as usize;
        let mut drawn = 0usize;
        let mut attempts = 0usize;
        while drawn < sample_size && attempts < population.saturating_mul(4) {
            attempts += 1;
            let index = (rng.next_f64() * population as f64) as usize;
            let index = index.min(population.saturating_sub(1));
            if picked.insert(index) {
                reasons.push((index, Reason::RandomSample));
                drawn += 1;
            }
        }

        // (a) non-dominated under L1 ∪ surrogate estimate.
        if !objectives.is_empty() {
            let ranked: Vec<Ranked> = objectives
                .iter()
                .zip(violations)
                .map(|(values, violation)| {
                    let objectives = values.iter().map(|v| Some(*v)).collect();
                    if *violation > 0.0 {
                        Ranked::infeasible(objectives, *violation)
                    } else {
                        Ranked::feasible(objectives)
                    }
                })
                .collect();
            let fronts = fast_nondominated_sort(&ranked);
            if let Some(first) = fronts.first() {
                for &index in first {
                    if picked.insert(index) {
                        reasons.push((index, Reason::Nondominated));
                    }
                }
            }
        }

        // (c) top-k by expected hypervolume gain, approximated by distance from
        // the current best in each objective. Filling whatever room is left.
        let mut by_gain: Vec<(usize, f64)> = (0..population)
            .filter(|index| !picked.contains(index))
            .map(|index| (index, expected_gain(&objectives[index], objectives)))
            .collect();
        by_gain.sort_by(|a, b| {
            b.1.partial_cmp(&a.1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.0.cmp(&b.0))
        });
        for (index, _) in by_gain {
            if picked.len() >= self.e_max {
                break;
            }
            if picked.insert(index) {
                reasons.push((index, Reason::ExpectedGain));
            }
        }

        // Order by [`Reason`]'s declaration order, which *is* the priority
        // (see its docs). Sorting by reason before index also keeps the result
        // independent of the order components were added above. The protocol
        // cuts this list at its budget; everything past the cut is deferred by
        // it, not dropped here.
        reasons.sort_by_key(|(index, reason)| (*reason, *index));
        let ranked: Vec<usize> = reasons.iter().map(|(index, _)| *index).collect();

        Rationed { ranked, reasons }
    }
}

/// A cheap stand-in for expected hypervolume gain.
///
/// Exact expected-HV-gain is exponential in the objective count and would cost
/// more than the L0 calls it is trying to save. This ranks by how far an
/// individual pushes past the current best in its own strongest objective,
/// which is the part of the gain that dominates when the front is not yet
/// dense.
fn expected_gain(candidate: &[f64], population: &[Vec<f64>]) -> f64 {
    let mut best_margin = f64::NEG_INFINITY;
    for (axis, value) in candidate.iter().enumerate() {
        let best = population
            .iter()
            .filter_map(|other| other.get(axis))
            .fold(f64::INFINITY, |a: f64, b| a.min(*b));
        best_margin = best_margin.max(best - value);
    }
    best_margin
}
