//! Quality indicators: hypervolume, IGD+, ε-indicator, Das–Dennis directions.
//!
//! Ticket P2-01. Specification: `03-nsga-loop-and-operators.md` §3.
//!
//! All minimization, all reference points given rather than inferred. An
//! indicator that picked its own reference point would make two runs
//! incomparable whenever their fronts differed in extent, which is exactly when
//! a comparison matters.

use std::cmp::Ordering;

use crate::{MetricsError, Result};

/// Exact hypervolume of a minimized front with respect to `reference`.
///
/// WFG's recursive formulation: the hypervolume of a set is the sum over its
/// members of that member's *exclusive* contribution, where each contribution
/// is computed against the set limited to the region the member dominates. It
/// is exponential in the worst case and exact, which is the right trade here —
/// §8 puts non-dominated sorting three to four orders of magnitude below
/// evaluation cost, and the indicator is computed once per generation over a
/// front of a few hundred points, not per individual.
///
/// Points not dominating the reference contribute nothing rather than a
/// negative volume.
pub fn hypervolume(front: &[Vec<f64>], reference: &[f64]) -> Result<f64> {
    for point in front {
        if point.len() != reference.len() {
            return Err(MetricsError::dimension(format!(
                "point has {} objectives against a {}-objective reference point",
                point.len(),
                reference.len()
            )));
        }
    }

    let dominating: Vec<Vec<f64>> = front
        .iter()
        .filter(|p| p.iter().zip(reference).all(|(v, r)| v < r))
        .cloned()
        .collect();

    Ok(wfg(&filter_dominated(&dominating), reference))
}

fn wfg(front: &[Vec<f64>], reference: &[f64]) -> f64 {
    // Sorted so the *recursion* is a function of the point set rather than of
    // arrival order. This is defensive, not load-bearing, and the distinction
    // is worth stating because it is easy to assume otherwise: WFG's
    // exclusive-contribution decomposition is order-independent by
    // construction, and removing this sort leaves the result bit-identical
    // even on deliberately adversarial float magnitudes. What it buys is that
    // the intermediate decomposition is stable too, so a later short-circuit
    // or tolerance cannot make the total depend on arrival order without
    // someone noticing. `TST-METR-003` asserts the property directly rather
    // than trusting this line for it.
    let mut ordered = front.to_vec();
    ordered.sort_by(|a, b| {
        for (x, y) in a.iter().zip(b.iter()) {
            match x.partial_cmp(y) {
                Some(Ordering::Equal) | None => continue,
                Some(other) => return other,
            }
        }
        Ordering::Equal
    });

    let mut total = 0.0;
    for (index, point) in ordered.iter().enumerate() {
        total += exclusive(point, &ordered[index + 1..], reference);
    }
    total
}

/// The volume `point` dominates that the rest of the set does not.
fn exclusive(point: &[f64], rest: &[Vec<f64>], reference: &[f64]) -> f64 {
    let inclusive: f64 = point
        .iter()
        .zip(reference)
        .map(|(v, r)| (r - v).max(0.0))
        .product();
    if rest.is_empty() {
        return inclusive;
    }
    let limited = filter_dominated(&limit_set(point, rest));
    inclusive - wfg(&limited, reference)
}

/// The rest of the set, clipped into the region `point` dominates.
fn limit_set(point: &[f64], rest: &[Vec<f64>]) -> Vec<Vec<f64>> {
    rest.iter()
        .map(|other| {
            other
                .iter()
                .zip(point)
                .map(|(o, p)| o.max(*p))
                .collect::<Vec<f64>>()
        })
        .collect()
}

/// Drop points dominated by another in the set, and exact duplicates.
fn filter_dominated(points: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let mut kept: Vec<Vec<f64>> = Vec::new();
    for (index, point) in points.iter().enumerate() {
        let dominated = points.iter().enumerate().any(|(other_index, other)| {
            if other_index == index {
                return false;
            }
            let weakly = other.iter().zip(point).all(|(o, p)| o <= p);
            let strictly = other.iter().zip(point).any(|(o, p)| o < p);
            // A duplicate is dropped once: the earlier copy survives.
            (weakly && strictly) || (weakly && other == point && other_index < index)
        });
        if !dominated {
            kept.push(point.clone());
        }
    }
    kept
}

/// IGD+ from a reference front to an approximation, minimized.
///
/// The `+` matters: the modified distance measures only the components in which
/// the approximation is *worse* than the reference point. Plain IGD rewards an
/// approximation for overshooting past the true front, so a front that is
/// wrong in a useful direction scores better than one that is right.
///
/// Zero exactly when every reference point is weakly dominated by some member
/// of the approximation, which is what `TST-METR-003` asserts.
pub fn igd_plus(approximation: &[Vec<f64>], reference: &[Vec<f64>]) -> Result<f64> {
    if reference.is_empty() {
        return Err(MetricsError::dimension(
            "IGD+ needs a non-empty reference front",
        ));
    }
    if approximation.is_empty() {
        return Ok(f64::INFINITY);
    }

    let mut total = 0.0;
    for point in reference {
        let mut nearest = f64::INFINITY;
        for candidate in approximation {
            if candidate.len() != point.len() {
                return Err(MetricsError::dimension(
                    "approximation and reference have different objective counts",
                ));
            }
            let distance: f64 = candidate
                .iter()
                .zip(point)
                .map(|(a, r)| (a - r).max(0.0).powi(2))
                .sum::<f64>()
                .sqrt();
            nearest = nearest.min(distance);
        }
        total += nearest;
    }
    Ok(total / reference.len() as f64)
}

/// Additive ε-indicator: the smallest `ε` such that the approximation shifted
/// by `ε` weakly dominates every reference point.
pub fn epsilon_additive(approximation: &[Vec<f64>], reference: &[Vec<f64>]) -> Result<f64> {
    if approximation.is_empty() || reference.is_empty() {
        return Err(MetricsError::dimension(
            "the epsilon indicator needs two non-empty fronts",
        ));
    }
    let mut worst = f64::NEG_INFINITY;
    for point in reference {
        let mut best = f64::INFINITY;
        for candidate in approximation {
            if candidate.len() != point.len() {
                return Err(MetricsError::dimension(
                    "approximation and reference have different objective counts",
                ));
            }
            let needed = candidate
                .iter()
                .zip(point)
                .map(|(a, r)| a - r)
                .fold(f64::NEG_INFINITY, f64::max);
            best = best.min(needed);
        }
        worst = worst.max(best);
    }
    Ok(worst)
}

/// Das–Dennis reference directions on the `M`-simplex with `p` divisions.
///
/// Returns `C(M+p-1, p)` points, each summing to 1. Generated recursively in a
/// fixed order, so two runs produce the same directions in the same positions —
/// niching associates individuals with directions by index, so a reordering
/// would silently change which individuals compete.
pub fn das_dennis(objectives: usize, divisions: usize) -> Result<Vec<Vec<f64>>> {
    if objectives == 0 {
        return Err(MetricsError::dimension(
            "reference directions need at least one objective",
        ));
    }
    let mut out = Vec::new();
    let mut current = vec![0usize; objectives];
    generate(&mut out, &mut current, 0, divisions, divisions);
    Ok(out)
}

fn generate(
    out: &mut Vec<Vec<f64>>,
    current: &mut Vec<usize>,
    index: usize,
    remaining: usize,
    divisions: usize,
) {
    if index == current.len() - 1 {
        current[index] = remaining;
        out.push(
            current
                .iter()
                .map(|&c| c as f64 / divisions as f64)
                .collect(),
        );
        return;
    }
    for take in 0..=remaining {
        current[index] = take;
        generate(out, current, index + 1, remaining - take, divisions);
    }
}

/// The two-layer construction for `M ≥ 4`.
///
/// §3's reason: a single layer at a usable `p` puts almost every direction on
/// the boundary of the simplex in high dimension, leaving the interior — where
/// balanced trade-offs live — uncovered. The inner layer is a shrunken copy
/// pulled towards the centroid.
pub fn das_dennis_two_layer(
    objectives: usize,
    outer: usize,
    inner: usize,
) -> Result<Vec<Vec<f64>>> {
    let mut directions = das_dennis(objectives, outer)?;
    if inner == 0 {
        return Ok(directions);
    }
    let centre = 1.0 / objectives as f64;
    for direction in das_dennis(objectives, inner)? {
        directions.push(
            direction
                .iter()
                .map(|value| centre + 0.5 * (value - centre))
                .collect(),
        );
    }
    Ok(directions)
}
