//! NSGA-III: normalization, association, and niching selection.
//!
//! Ticket P2-03. Specification: `03-nsga-loop-and-operators.md` §3.
//!
//! ## Why this replaces crowding distance rather than supplementing it
//!
//! §2's reason, not restated: as `M` grows, nearly every individual becomes
//! non-dominated, so selection pressure comes almost entirely from the
//! diversity operator — and crowding distance is a poor one in high dimension,
//! because a sum of per-axis spacings stops discriminating once most points are
//! on the boundary of the front. Reference directions supply structured
//! coverage instead: each direction is a niche, and selection fills the
//! emptiest niches first.
//!
//! ## Everything here is deterministic
//!
//! `DEP-SRCH-009`. Deb's algorithm makes two random choices — which tied
//! least-occupied reference point to fill, and which member of an occupied
//! niche to take. Both are deterministic here: the lowest reference index, and
//! the minimum perpendicular distance. `GATE-DET-001` requires bit-identical
//! output across thread counts, and a rule with no randomness in it cannot
//! depend on scheduling, where a correctly-seeded one still has to be proved
//! not to.

use crate::dominance::Ranked;

/// The normalized objective space for one generation.
#[derive(Clone, Debug, PartialEq)]
pub struct Normalization {
    /// Running per-axis minimum. Never regresses across generations.
    pub ideal: Vec<f64>,
    /// Per-axis intercepts of the hyperplane through the extreme points.
    pub intercepts: Vec<f64>,
    /// Whether a degeneracy guard replaced a computed intercept.
    pub degenerate: bool,
}

impl Normalization {
    /// Translate and scale a point into `[0, 1]`-ish normalized space.
    pub fn apply(&self, objectives: &[f64]) -> Vec<f64> {
        objectives
            .iter()
            .zip(&self.ideal)
            .zip(&self.intercepts)
            .map(|((value, ideal), intercept)| {
                let span = intercept - ideal;
                // A zero span means every member shares this axis. Contributing
                // zero says "this axis distinguishes nobody", which is true;
                // dividing would produce NaN and poison every distance.
                //
                // Defensive rather than load-bearing for `normalize`'s own
                // output: its degenerate fallback already sets an intercept of
                // `ideal + 1` when the per-axis maximum is zero, so a span of
                // zero cannot arrive from there. This guard is for a
                // `Normalization` a caller built directly, and the distinction
                // is worth stating because a negative control removing it does
                // not fire through `normalize`.
                if span.abs() <= f64::EPSILON {
                    0.0
                } else {
                    (value - ideal) / span
                }
            })
            .collect()
    }
}

/// The standard NSGA-III normalization, with §3's two degeneracy guards.
///
/// `previous_ideal` carries the running minimum forward: the ideal point is
/// defined over every generation seen, not over the current one, so a
/// generation that happens to lose its best member on an axis does not move
/// the origin back and rescale the whole population.
pub fn normalize(points: &[Vec<f64>], previous_ideal: Option<&[f64]>) -> Normalization {
    let objectives = points.first().map(Vec::len).unwrap_or(0);
    if objectives == 0 || points.is_empty() {
        return Normalization {
            ideal: Vec::new(),
            intercepts: Vec::new(),
            degenerate: true,
        };
    }

    let mut ideal = vec![f64::INFINITY; objectives];
    for point in points {
        for (axis, value) in point.iter().enumerate() {
            ideal[axis] = ideal[axis].min(*value);
        }
    }
    if let Some(previous) = previous_ideal {
        for (axis, value) in previous.iter().enumerate().take(objectives) {
            ideal[axis] = ideal[axis].min(*value);
        }
    }

    // Extreme points: for each axis, the member minimizing an achievement
    // scalarizing function weighted almost entirely towards that axis. The
    // 1e-6 off-axis weight is the standard formulation and is what keeps the
    // ASF from being indifferent among all points sharing the axis value.
    let translated: Vec<Vec<f64>> = points
        .iter()
        .map(|p| p.iter().zip(&ideal).map(|(v, i)| v - i).collect())
        .collect();

    let mut extremes: Vec<Vec<f64>> = Vec::with_capacity(objectives);
    for axis in 0..objectives {
        let mut best = 0usize;
        let mut best_value = f64::INFINITY;
        for (index, point) in translated.iter().enumerate() {
            let asf = point
                .iter()
                .enumerate()
                .map(|(a, v)| if a == axis { *v } else { v / 1e-6 })
                .fold(f64::NEG_INFINITY, f64::max);
            // Strictly less, so the lowest index wins a tie. `DEP-SRCH-009`.
            if asf < best_value {
                best_value = asf;
                best = index;
            }
        }
        extremes.push(translated[best].clone());
    }

    let per_axis_max: Vec<f64> = (0..objectives)
        .map(|axis| {
            translated
                .iter()
                .map(|p| p[axis])
                .fold(f64::NEG_INFINITY, f64::max)
        })
        .collect();

    // The one validity test. An intercept must be finite and positive: a
    // singular extreme-point matrix produces an infinity or a NaN here, a
    // degenerate front produces a non-positive value, and §3's two guards have
    // the same answer for both — fall back to the per-axis maximum.
    match solve_intercepts(&extremes) {
        Some(intercepts)
            if intercepts
                .iter()
                .all(|v| v.is_finite() && *v > f64::EPSILON) =>
        {
            Normalization {
                intercepts: intercepts.iter().zip(&ideal).map(|(a, i)| a + i).collect(),
                ideal,
                degenerate: false,
            }
        }
        // §3's guards, as one branch because they have one answer: a singular
        // extreme-point matrix and a non-positive intercept both mean the
        // hyperplane is not usable, and both fall back to the per-axis maximum
        // of the current front. A degenerate front — two perfectly correlated
        // objectives, say — otherwise divides by ~0 and every distance becomes
        // meaningless in a way nothing downstream can detect.
        _ => Normalization {
            intercepts: per_axis_max
                .iter()
                .zip(&ideal)
                .map(|(m, i)| if *m > 0.0 { m + i } else { i + 1.0 })
                .collect(),
            ideal,
            degenerate: true,
        },
    }
}

/// Solve `E · a = 1` for the intercepts, or `None` if `E` is singular.
///
/// Gaussian elimination with partial pivoting. Small `M` — at most a handful of
/// objectives — so the direct method is right and an iterative one would be
/// pointless machinery.
fn solve_intercepts(extremes: &[Vec<f64>]) -> Option<Vec<f64>> {
    let n = extremes.len();
    if n == 0 || extremes.iter().any(|row| row.len() != n) {
        return None;
    }
    let mut matrix: Vec<Vec<f64>> = extremes.to_vec();
    let mut rhs = vec![1.0f64; n];

    for column in 0..n {
        let pivot = (column..n).max_by(|a, b| {
            matrix[*a][column]
                .abs()
                .partial_cmp(&matrix[*b][column].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        // No zero-pivot guard here, deliberately. A singular matrix divides by
        // zero, which in IEEE arithmetic gives an infinity or a NaN rather than
        // a trap, and both fail the single validity test below. Checking here
        // as well would be two mechanisms for one rule: they agreed, so
        // removing either changed nothing and neither was testable — which is
        // exactly how the `UNKNOWN` screen went unnoticed in `cf-metrics`
        // until a negative control refused to fire.
        matrix.swap(column, pivot);
        rhs.swap(column, pivot);

        for row in (column + 1)..n {
            let factor = matrix[row][column] / matrix[column][column];
            for k in column..n {
                matrix[row][k] -= factor * matrix[column][k];
            }
            rhs[row] -= factor * rhs[column];
        }
    }

    let mut solution = vec![0.0f64; n];
    for row in (0..n).rev() {
        let mut value = rhs[row];
        for (k, solved) in solution.iter().enumerate().take(n).skip(row + 1) {
            value -= matrix[row][k] * solved;
        }
        solution[row] = value / matrix[row][row];
    }
    // The intercept on axis `i` is `1 / a_i`.
    solution
        .iter()
        .map(|a| {
            if a.abs() <= 1e-12 {
                None
            } else {
                Some(1.0 / a)
            }
        })
        .collect()
}

/// Which reference direction a member is associated with, and how far off it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Association {
    pub reference: usize,
    pub distance: f64,
}

/// Perpendicular distance from a normalized point to a reference direction.
pub fn perpendicular_distance(point: &[f64], direction: &[f64]) -> f64 {
    let norm_squared: f64 = direction.iter().map(|d| d * d).sum();
    if norm_squared <= f64::EPSILON {
        return point.iter().map(|p| p * p).sum::<f64>().sqrt();
    }
    let dot: f64 = point.iter().zip(direction).map(|(p, d)| p * d).sum();
    let scale = dot / norm_squared;
    point
        .iter()
        .zip(direction)
        .map(|(p, d)| (p - scale * d).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Associate each normalized point with its nearest reference direction.
///
/// Ties go to the lowest reference index, so association is a function of the
/// point set. `DEP-SRCH-009`.
pub fn associate(points: &[Vec<f64>], directions: &[Vec<f64>]) -> Vec<Association> {
    points
        .iter()
        .map(|point| {
            let mut best = Association {
                reference: 0,
                distance: f64::INFINITY,
            };
            for (index, direction) in directions.iter().enumerate() {
                let distance = perpendicular_distance(point, direction);
                if distance < best.distance {
                    best = Association {
                        reference: index,
                        distance,
                    };
                }
            }
            best
        })
        .collect()
}

/// Deb's niching selection (Algorithm 4), made deterministic.
///
/// `already` holds the members admitted from earlier fronts, `candidates` the
/// splitting front, both as indices into the same association list. Returns the
/// `slots` chosen candidates, in ascending index order so the caller cannot
/// depend on selection order.
pub fn niching_select(
    associations: &[Association],
    already: &[usize],
    candidates: &[usize],
    reference_count: usize,
    slots: usize,
) -> Vec<usize> {
    let mut counts = vec![0usize; reference_count.max(1)];
    for &index in already {
        counts[associations[index].reference] += 1;
    }

    let mut remaining: Vec<usize> = candidates.to_vec();
    remaining.sort_unstable();
    let mut chosen: Vec<usize> = Vec::with_capacity(slots);
    // A reference point with no remaining candidate is excluded rather than
    // repeatedly re-chosen; without this the loop spins on an empty niche.
    let mut exhausted = vec![false; counts.len()];

    while chosen.len() < slots && !remaining.is_empty() {
        let niche = counts
            .iter()
            .enumerate()
            .filter(|(index, _)| !exhausted[*index])
            // `min_by_key` on `(count, index)` gives the lowest index among
            // ties, which is where Deb samples. `DEP-SRCH-009`.
            .min_by_key(|(index, count)| (**count, *index))
            .map(|(index, _)| index);

        let Some(niche) = niche else { break };

        let pick = remaining
            .iter()
            .copied()
            .filter(|&member| associations[member].reference == niche)
            // Minimum perpendicular distance, then lowest index. Deb takes a
            // random member when the niche is already occupied; taking the
            // closest is deterministic and, for an empty niche, is Deb's own
            // rule anyway.
            .min_by(|a, b| {
                associations[*a]
                    .distance
                    .partial_cmp(&associations[*b].distance)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.cmp(b))
            });

        match pick {
            Some(member) => {
                chosen.push(member);
                remaining.retain(|&m| m != member);
                counts[niche] += 1;
            }
            None => exhausted[niche] = true,
        }
    }

    chosen.sort_unstable();
    chosen
}

/// Normalized objective vectors for a set of fully-known feasible members.
///
/// Members with an `UNKNOWN` term are excluded and reported separately: they
/// cannot be placed in objective space at all, so associating them with a
/// direction would be inventing the coordinate. §4 puts them at the back of
/// their feasibility class, which is the caller's job — this only says which
/// they are.
pub fn placeable(population: &[Ranked], subset: &[usize]) -> (Vec<usize>, Vec<usize>) {
    let mut known = Vec::new();
    let mut unknown = Vec::new();
    for &index in subset {
        if population[index].has_unknown() {
            unknown.push(index);
        } else {
            known.push(index);
        }
    }
    (known, unknown)
}
