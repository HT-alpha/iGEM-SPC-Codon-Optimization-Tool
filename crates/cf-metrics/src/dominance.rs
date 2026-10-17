//! Non-dominated sorting and crowding distance.
//!
//! Ticket P2-01. Specification: `03-nsga-loop-and-operators.md` §4.
//!
//! ## Everything here is minimization
//!
//! Objective *sense* belongs to the provider (`TermProvider::direction`) and is
//! applied once, on the way in. A metric that also knew about senses would be a
//! second place the direction is interpreted, and this codebase has paid for a
//! direction interpreted twice: the job layer signed weights the tiers had
//! already signed, and a `minimize` term came out maximized with every document
//! still reading correctly.
//!
//! ## `UNKNOWN` is a value, not a missing number
//!
//! A term whose provider missed its latency budget is `None`, and §4 is
//! explicit that it must never be imputed: an individual with any unknown term
//! is *incomparable* to every other feasible individual rather than being
//! guessed at. Imputing a mean or a worst case biases the front in a way
//! nothing in the output reveals.

use std::cmp::Ordering;

/// The result of comparing two individuals.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Dominance {
    /// The left individual is better.
    Left,
    /// The right individual is better.
    Right,
    /// Neither dominates: mutually non-dominated, equal, or incomparable.
    Neither,
}

/// What ranking needs to know about an individual.
///
/// Deliberately not a genotype: `cf-metrics` ranks numbers and knows nothing
/// about sequences, which is what lets the same code rank a synthetic fixture
/// and a real population.
#[derive(Clone, Debug, PartialEq)]
pub struct Ranked {
    /// Objective values, all minimized. `None` is `UNKNOWN`.
    pub objectives: Vec<Option<f64>>,
    /// Total actionable constraint violation. Zero means feasible.
    ///
    /// `DEP-SRCH-008`: *actionable*. A violation no synonymous assignment can
    /// remove is not the individual's fault, and charging it here would rank
    /// the population by a quantity none of it can influence.
    pub violation: f64,
}

impl Ranked {
    pub fn feasible(objectives: Vec<Option<f64>>) -> Self {
        Self {
            objectives,
            violation: 0.0,
        }
    }

    pub fn infeasible(objectives: Vec<Option<f64>>, violation: f64) -> Self {
        Self {
            objectives,
            violation,
        }
    }

    pub fn is_feasible(&self) -> bool {
        self.violation <= 0.0
    }

    pub fn has_unknown(&self) -> bool {
        self.objectives.iter().any(Option::is_none)
    }
}

/// Deb's constraint-domination, with `UNKNOWN` as a first-class value.
///
/// The order of the tests is the definition and not a convenience:
///
/// 1. A feasible individual beats an infeasible one, whatever the objectives
///    say. This is what drives a population towards feasibility while keeping
///    infeasible individuals as stepping stones across a constraint barrier —
///    discarding them is the standard way to make a tightly-constrained codon
///    problem return nothing.
/// 2. Two infeasible individuals compare on total violation alone, which is
///    scalar and therefore never incomparable.
/// 3. Two feasible individuals with any unknown term are incomparable.
/// 4. Otherwise, Pareto comparison.
pub fn constraint_dominance(left: &Ranked, right: &Ranked) -> Dominance {
    match (left.is_feasible(), right.is_feasible()) {
        (true, false) => return Dominance::Left,
        (false, true) => return Dominance::Right,
        (false, false) => {
            return match left.violation.partial_cmp(&right.violation) {
                Some(Ordering::Less) => Dominance::Left,
                Some(Ordering::Greater) => Dominance::Right,
                _ => Dominance::Neither,
            };
        }
        (true, true) => {}
    }

    // The one place `UNKNOWN` is handled. `pareto_dominance` below takes known
    // values, so this screen cannot be bypassed and cannot be duplicated: an
    // earlier version also returned `Neither` on an unknown *inside* the Pareto
    // comparison, which made the two mechanisms redundant — and a negative
    // control that removed this screen changed nothing, because the other one
    // silently covered for it. Two implementations of one rule is one too many
    // even when they agree.
    if left.has_unknown() || right.has_unknown() {
        return Dominance::Neither;
    }
    // The values are known. Comparing them in place is the same `<` walk as
    // copying both vectors out first; the copy was one allocation per
    // comparison, and a 1 600-pair primer front does about 1.3 million of
    // those.
    pareto_dominance(&left.objectives, &right.objectives)
}

/// Pareto comparison over fully-known, minimized objectives.
fn pareto_dominance(left: &[Option<f64>], right: &[Option<f64>]) -> Dominance {
    let mut left_better = false;
    let mut right_better = false;

    for (a, b) in left.iter().zip(right.iter()) {
        let a = a.expect("screened above");
        let b = b.expect("screened above");
        if a < b {
            left_better = true;
        } else if b < a {
            right_better = true;
        }
        if left_better && right_better {
            return Dominance::Neither;
        }
    }

    match (left_better, right_better) {
        (true, false) => Dominance::Left,
        (false, true) => Dominance::Right,
        _ => Dominance::Neither,
    }
}

/// Fast non-dominated sort. Returns fronts, best first, as index lists.
///
/// `O(M · N²)`, which §4 chooses deliberately: at `M=5, N=212` that is about
/// `2.2 × 10⁵` comparisons per generation, three to four orders of magnitude
/// below evaluation cost, so the asymptotically better algorithms would buy
/// nothing and cost clarity.
///
/// Each front is returned in ascending index order, so the result is a function
/// of the input set rather than of the traversal — the same discipline the
/// tier-2 commit follows, and for the same reason.
pub fn fast_nondominated_sort(population: &[Ranked]) -> Vec<Vec<usize>> {
    let n = population.len();
    let mut dominated_by: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut domination_count = vec![0usize; n];
    let mut fronts: Vec<Vec<usize>> = Vec::new();
    let mut current: Vec<usize> = Vec::new();

    for i in 0..n {
        for j in (i + 1)..n {
            match constraint_dominance(&population[i], &population[j]) {
                Dominance::Left => {
                    dominated_by[i].push(j);
                    domination_count[j] += 1;
                }
                Dominance::Right => {
                    dominated_by[j].push(i);
                    domination_count[i] += 1;
                }
                Dominance::Neither => {}
            }
        }
    }

    for (index, count) in domination_count.iter().enumerate() {
        if *count == 0 {
            current.push(index);
        }
    }

    while !current.is_empty() {
        let mut next: Vec<usize> = Vec::new();
        for &index in &current {
            for &other in &dominated_by[index] {
                domination_count[other] -= 1;
                if domination_count[other] == 0 {
                    next.push(other);
                }
            }
        }
        next.sort_unstable();
        fronts.push(std::mem::take(&mut current));
        current = next;
    }

    fronts
}

/// Crowding distance within one front, indexed as `front` is.
///
/// Two edge cases carry the whole behaviour, and `P2-01` shipped without a test
/// for either because none was specified (`TST-METR-004` now covers them):
///
/// * **Boundary individuals take infinity**, so the extremes of a front are
///   never crowded out. Without it a front collapses towards its middle over
///   generations and the extreme trade-offs — often the interesting designs —
///   are lost.
/// * **A front whose members share an objective value has zero range in it.**
///   Normalizing by that range divides by zero. The term contributes nothing
///   instead, which is the correct reading: an objective on which everyone
///   agrees says nothing about who is crowded.
pub fn crowding_distance(population: &[Ranked], front: &[usize]) -> Vec<f64> {
    let count = front.len();
    let mut distance = vec![0.0f64; count];
    if count == 0 {
        return distance;
    }
    let objectives = population[front[0]].objectives.len();

    for objective in 0..objectives {
        // An unknown cannot be ordered, so it cannot contribute a distance.
        // Treated as absent rather than as an extreme: an individual is not
        // crowded or uncrowded on a term nobody measured for it.
        let mut order: Vec<usize> = (0..count)
            .filter(|&i| population[front[i]].objectives[objective].is_some())
            .collect();
        if order.is_empty() {
            continue;
        }
        if order.len() == 1 {
            // The only member with a known value on this objective is both of
            // its boundaries, so it is never crowded out. Skipping instead
            // would leave a one-member front at distance zero and let
            // selection discard the only candidate it has.
            distance[order[0]] = f64::INFINITY;
            continue;
        }
        let value = |i: usize| population[front[i]].objectives[objective].unwrap();
        // Ties broken on the front index so the sort is total: two individuals
        // with the same value must not swap depending on the sort's internals.
        order.sort_by(|&a, &b| {
            value(a)
                .partial_cmp(&value(b))
                .unwrap_or(Ordering::Equal)
                .then_with(|| front[a].cmp(&front[b]))
        });

        let low = value(order[0]);
        let high = value(order[order.len() - 1]);
        distance[order[0]] = f64::INFINITY;
        distance[order[order.len() - 1]] = f64::INFINITY;

        let range = high - low;
        if range <= 0.0 {
            continue;
        }
        for position in 1..order.len() - 1 {
            let index = order[position];
            if distance[index].is_infinite() {
                continue;
            }
            distance[index] += (value(order[position + 1]) - value(order[position - 1])) / range;
        }
    }

    distance
}
