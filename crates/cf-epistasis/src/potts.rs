//! Synonymous-restricted Potts pseudo-likelihood.
//!
//! The objective is concave in the parameters (`01-epistatic-surface-fit.md`
//! §2.4), so the fit starts at zero and uses L-BFGS. After every accepted
//! step the parameters are projected onto the zero-sum gauge inside each
//! synonymous family. `μ = 0`: the group-lasso path is a different optimizer
//! and is not this fit.

use cf_types::{Codon, GeneticCode};

const HISTORY: usize = 8;
const TRANSLATION_TABLE: u8 = 11;

/// One fitted context model. `h` is 64 fields. `couplings` is `W` matrices
/// of 64×64 in row-major order, offset `d` at `d * 4096`.
#[derive(Clone, Debug)]
pub struct PottsModel {
    pub window: usize,
    pub h: Vec<f64>,
    pub couplings: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct FitReport {
    pub model: PottsModel,
    pub nats_per_codon: f64,
    pub n_positions: u64,
    pub iterations: u32,
    pub gradient_norm: f64,
    pub converged: bool,
    pub lambda_h: f64,
    pub lambda_j: f64,
}

#[derive(Clone, Debug)]
pub struct UsageBaseline {
    pub nats_per_codon: f64,
    pub n_scored: u64,
    pub n_unseen: u64,
}

struct Families {
    /// Synonymous codons of each codon, length 1–8.
    members: Vec<Vec<u8>>,
}

impl Families {
    fn bacterial() -> Self {
        let table = GeneticCode::get(TRANSLATION_TABLE).expect("table 11 is vendored");
        let mut members = Vec::with_capacity(64);
        for index in 0..64 {
            let codon = Codon::from_index(index as u8).expect("0..64");
            members.push(table.synonyms(codon).iter().map(|c| c.index()).collect());
        }
        Self { members }
    }
}

fn parse_genes(sequences: &[String]) -> Result<Vec<Vec<u8>>, String> {
    let mut genes = Vec::with_capacity(sequences.len());
    for (which, sequence) in sequences.iter().enumerate() {
        let upper = sequence.trim().to_ascii_uppercase();
        if upper.len() % 3 != 0 {
            return Err(format!("sequence {which} is not a whole number of codons"));
        }
        let bytes = upper.as_bytes();
        let mut codons = Vec::with_capacity(bytes.len() / 3);
        for chunk in bytes.chunks(3) {
            let codon = Codon::from_bytes(chunk).map_err(|err| err.detail.to_string())?;
            codons.push(codon.index());
        }
        genes.push(codons);
    }
    Ok(genes)
}

fn is_family_leader(index: usize, family: &[u8]) -> bool {
    family.len() > 1 && *family.iter().min().expect("family") == index as u8
}

fn project(theta: &mut [f64], window: usize, families: &Families) {
    for index in 0..64 {
        let family = &families.members[index];
        if !is_family_leader(index, family) {
            // The mean is subtracted once per family, from its lowest codon
            // index, so a codon and its synonyms are not each shifted.
            continue;
        }
        let mean = family.iter().map(|&c| theta[c as usize]).sum::<f64>() / family.len() as f64;
        for &codon in family {
            theta[codon as usize] -= mean;
        }
    }
    for offset in 0..window {
        let base = 64 + offset * 4096;
        for _ in 0..6 {
            for column in 0..64 {
                for index in 0..64 {
                    let family = &families.members[index];
                    if !is_family_leader(index, family) {
                        continue;
                    }
                    let mean = family
                        .iter()
                        .map(|&row| theta[base + row as usize * 64 + column])
                        .sum::<f64>()
                        / family.len() as f64;
                    for &row in family {
                        theta[base + row as usize * 64 + column] -= mean;
                    }
                }
            }
            for row in 0..64 {
                for index in 0..64 {
                    let family = &families.members[index];
                    if !is_family_leader(index, family) {
                        continue;
                    }
                    let mean = family
                        .iter()
                        .map(|&column| theta[base + row * 64 + column as usize])
                        .sum::<f64>()
                        / family.len() as f64;
                    for &column in family {
                        theta[base + row * 64 + column as usize] -= mean;
                    }
                }
            }
        }
    }
}

/// Negative penalized pseudo-likelihood and its gradient.
fn evaluate(
    genes: &[Vec<u8>],
    theta: &[f64],
    window: usize,
    lambda_h: f64,
    lambda_j: f64,
    families: &Families,
) -> (f64, Vec<f64>, u64) {
    let n_param = theta.len();
    let mut gradient = vec![0.0; n_param];
    let mut negative = 0.0;
    let mut positions = 0u64;
    let h = &theta[..64];
    for gene in genes {
        let length = gene.len();
        for i in 0..length {
            let observed = gene[i] as usize;
            let family = &families.members[observed];
            let size = family.len();
            if size <= 1 {
                continue;
            }
            positions += 1;
            let mut field = [0.0; 8];
            let mut observed_slot = 0usize;
            for (slot, &codon) in family.iter().enumerate() {
                let codon = codon as usize;
                if codon == observed {
                    observed_slot = slot;
                }
                let mut value = h[codon];
                for distance in 1..=window {
                    let matrix = 64 + (distance - 1) * 4096;
                    if i + distance < length {
                        let neighbor = gene[i + distance] as usize;
                        value += theta[matrix + codon * 64 + neighbor];
                    }
                    if i >= distance {
                        let neighbor = gene[i - distance] as usize;
                        value += theta[matrix + neighbor * 64 + codon];
                    }
                }
                field[slot] = value;
            }
            let mut peak = f64::NEG_INFINITY;
            for value in field.iter().take(size) {
                peak = peak.max(*value);
            }
            let mut partition = 0.0;
            let mut probability = [0.0; 8];
            for (value, field_value) in probability.iter_mut().zip(field.iter()).take(size) {
                *value = (field_value - peak).exp();
                partition += *value;
            }
            for value in probability.iter_mut().take(size) {
                *value /= partition;
            }
            negative -= probability[observed_slot].ln();
            for slot in 0..size {
                let codon = family[slot] as usize;
                let residual = probability[slot] - if slot == observed_slot { 1.0 } else { 0.0 };
                gradient[codon] += residual;
                for distance in 1..=window {
                    let matrix = 64 + (distance - 1) * 4096;
                    if i + distance < length {
                        let neighbor = gene[i + distance] as usize;
                        gradient[matrix + codon * 64 + neighbor] += residual;
                    }
                    if i >= distance {
                        let neighbor = gene[i - distance] as usize;
                        gradient[matrix + neighbor * 64 + codon] += residual;
                    }
                }
            }
        }
    }
    for index in 0..64 {
        negative += lambda_h * theta[index] * theta[index];
        gradient[index] += 2.0 * lambda_h * theta[index];
    }
    for index in 64..n_param {
        negative += lambda_j * theta[index] * theta[index];
        gradient[index] += 2.0 * lambda_j * theta[index];
    }
    (negative, gradient, positions)
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn lbfgs_direction(gradient: &[f64], history: &[(Vec<f64>, Vec<f64>)]) -> Vec<f64> {
    let mut q = gradient.to_vec();
    let mut alpha = vec![0.0; history.len()];
    for (index, (step, change)) in history.iter().enumerate().rev() {
        let rho = 1.0 / dot(change, step);
        alpha[index] = rho * dot(step, &q);
        for (slot, value) in q.iter_mut().enumerate() {
            *value -= alpha[index] * change[slot];
        }
    }
    if let Some((step, change)) = history.last() {
        let scale = dot(step, change) / dot(change, change);
        for value in &mut q {
            *value *= scale;
        }
    }
    for (index, (step, change)) in history.iter().enumerate() {
        let rho = 1.0 / dot(change, step);
        let beta = rho * dot(change, &q);
        for (slot, value) in q.iter_mut().enumerate() {
            *value += step[slot] * (alpha[index] - beta);
        }
    }
    for value in &mut q {
        *value = -*value;
    }
    q
}

/// Fit `W` on `sequences`. `lambda_h` and `lambda_j` are numerical
/// stabilizers on the squared parameters, recorded with the artifact.
/// They are not measured expression rates.
pub fn fit_potts(
    sequences: &[String],
    window: usize,
    lambda_h: f64,
    lambda_j: f64,
    max_iterations: u32,
) -> Result<FitReport, String> {
    fit_potts_from(sequences, window, lambda_h, lambda_j, max_iterations, None)
}

/// `initial` continues a previous fit. It is the concatenated `h` and
/// couplings, already in the zero-sum gauge. L-BFGS history starts empty.
pub fn fit_potts_from(
    sequences: &[String],
    window: usize,
    lambda_h: f64,
    lambda_j: f64,
    max_iterations: u32,
    initial: Option<&[f64]>,
) -> Result<FitReport, String> {
    if window == 0 || window > 20 {
        return Err(format!("window {window} is outside 1..=20"));
    }
    if lambda_h < 0.0 || lambda_j < 0.0 {
        return Err("regularization weights are negative".to_string());
    }
    let genes = parse_genes(sequences)?;
    if genes.is_empty() {
        return Err("the fit was given no sequences".to_string());
    }
    let families = Families::bacterial();
    let n_param = 64 + window * 4096;
    let mut theta = vec![0.0; n_param];
    if let Some(start) = initial {
        if start.len() != n_param {
            return Err(format!(
                "initial parameters have length {}, window {window} needs {n_param}",
                start.len()
            ));
        }
        theta.copy_from_slice(start);
    }
    let mut history: Vec<(Vec<f64>, Vec<f64>)> = Vec::new();
    let mut previous_theta: Option<Vec<f64>> = None;
    let mut previous_gradient: Option<Vec<f64>> = None;
    let mut iterations = 0u32;
    let mut gradient_norm: f64;
    let mut converged = false;
    let mut positions = 0u64;
    let mut value: f64;
    for _ in 0..max_iterations {
        project(&mut theta, window, &families);
        let (next_value, mut gradient, counted) =
            evaluate(&genes, &theta, window, lambda_h, lambda_j, &families);
        // The likelihood is unchanged by a gauge shift. Projecting the
        // gradient drops that null-space component so a step is not spent
        // on a direction the following projection would cancel.
        project(&mut gradient, window, &families);
        value = next_value;
        positions = counted;
        iterations += 1;
        gradient_norm = dot(&gradient, &gradient).sqrt();
        let scale = (counted as f64).sqrt().max(1.0);
        if gradient_norm / scale < 1e-4 {
            converged = true;
            break;
        }
        if let (Some(before), Some(before_gradient)) = (&previous_theta, &previous_gradient) {
            let step: Vec<f64> = theta
                .iter()
                .zip(before)
                .map(|(now, old)| now - old)
                .collect();
            let change: Vec<f64> = gradient
                .iter()
                .zip(before_gradient)
                .map(|(now, old)| now - old)
                .collect();
            if dot(&change, &step) > 1e-8 {
                history.push((step, change));
                if history.len() > HISTORY {
                    history.remove(0);
                }
            }
        }
        let direction = lbfgs_direction(&gradient, &history);
        let slope = dot(&gradient, &direction);
        if slope >= 0.0 {
            break;
        }
        let mut alpha = 1.0;
        let mut accepted = false;
        for _ in 0..16 {
            let mut trial = theta.clone();
            for (slot, step) in trial.iter_mut().enumerate() {
                *step += alpha * direction[slot];
            }
            project(&mut trial, window, &families);
            let (trial_value, _, _) =
                evaluate(&genes, &trial, window, lambda_h, lambda_j, &families);
            if trial_value <= value + 1e-4 * alpha * slope {
                previous_theta = Some(theta);
                previous_gradient = Some(gradient);
                theta = trial;
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if !accepted {
            break;
        }
    }
    project(&mut theta, window, &families);
    let (final_value, mut final_gradient, counted) =
        evaluate(&genes, &theta, window, lambda_h, lambda_j, &families);
    project(&mut final_gradient, window, &families);
    gradient_norm = dot(&final_gradient, &final_gradient).sqrt();
    let nats = if counted == 0 {
        0.0
    } else {
        // The penalized objective is not the reported likelihood. Recompute
        // the mean log probability with the penalties off.
        let (plain, _, _) = evaluate(&genes, &theta, window, 0.0, 0.0, &families);
        -plain / counted as f64
    };
    let _ = final_value;
    Ok(FitReport {
        model: PottsModel {
            window,
            h: theta[..64].to_vec(),
            couplings: theta[64..].to_vec(),
        },
        nats_per_codon: nats,
        n_positions: positions,
        iterations,
        gradient_norm,
        converged,
        lambda_h,
        lambda_j,
    })
}

pub fn score_naturalness(sequences: &[String], model: &PottsModel) -> Result<(f64, u64), String> {
    let genes = parse_genes(sequences)?;
    let families = Families::bacterial();
    let mut theta = model.h.clone();
    theta.extend_from_slice(&model.couplings);
    if theta.len() != 64 + model.window * 4096 {
        return Err("couplings do not match the window".to_string());
    }
    let (negative, _, counted) = evaluate(&genes, &theta, model.window, 0.0, 0.0, &families);
    if counted == 0 {
        return Err("no position had a synonymous choice".to_string());
    }
    Ok((-negative / counted as f64, counted))
}

/// Mean log probability of one codon assignment under `model`.
///
/// The same quantity `score_naturalness` reports for a DNA string: higher
/// (less negative) is closer to the training conditional. Positions whose
/// synonymous family has one codon do not enter the mean.
pub fn mean_log_probability(codons: &[u8], model: &PottsModel) -> Result<f64, String> {
    let (theta, families) = prepared(model)?;
    let (total, counted) =
        sum_log_probability(codons, &theta, model.window, &families, 0..codons.len());
    if counted == 0 {
        return Err("no position had a synonymous choice".to_string());
    }
    Ok(total / counted as f64)
}

/// Exact change in `mean_log_probability` between two equal-length assignments.
///
/// Only positions within `window` of a changed codon are recomputed. A
/// synonymous edit does not change how many positions are counted, so the
/// change in the mean is the change in that local sum divided by the count.
pub fn mean_log_probability_delta(
    before: &[u8],
    after: &[u8],
    model: &PottsModel,
) -> Result<f64, String> {
    if before.len() != after.len() {
        return Err("the two assignments have different lengths".to_string());
    }
    let (theta, families) = prepared(model)?;
    let (_, total_count) =
        sum_log_probability(before, &theta, model.window, &families, 0..before.len());
    delta_of_slices(before, after, &theta, model.window, &families, total_count)
}

/// Same change as [`mean_log_probability_delta`], writing `edits` into
/// `codons` for the after-sum and restoring them afterwards.
///
/// `total_count` is the number of positions with a synonymous choice on the
/// unedited assignment. A synonymous edit does not change it, so the caller
/// can compute it once per sequence instead of once per synonym.
fn delta_in_place(
    codons: &mut [u8],
    edits: &[(usize, u8)],
    theta: &[f64],
    window: usize,
    families: &Families,
    total_count: u64,
) -> Result<f64, String> {
    let mut saved: Vec<(usize, u8)> = Vec::with_capacity(edits.len());
    for &(index, codon) in edits {
        if index >= codons.len() {
            restore(codons, &saved);
            return Err(format!("codon {index} is outside the sequence"));
        }
        if saved.iter().all(|(seen, _)| *seen != index) {
            saved.push((index, codons[index]));
        }
        codons[index] = codon;
    }
    let mut changed = Vec::new();
    for &(index, previous) in &saved {
        if codons[index] != previous {
            changed.push(index);
        }
    }
    changed.sort_unstable();
    if changed.is_empty() {
        restore(codons, &saved);
        return Ok(0.0);
    }
    restore(codons, &saved);
    let affected = affected_positions(codons.len(), &changed, window);
    let (before_sum, before_count) =
        sum_log_probability(codons, theta, window, families, affected.iter().copied());
    for &(index, codon) in edits {
        codons[index] = codon;
    }
    let (after_sum, after_count) =
        sum_log_probability(codons, theta, window, families, affected.iter().copied());
    restore(codons, &saved);
    if total_count == 0 || before_count != after_count {
        return Err("no position had a synonymous choice".to_string());
    }
    Ok((after_sum - before_sum) / total_count as f64)
}

/// Theta, families, and the synonymous-position count for one assignment.
///
/// Built once per delta table. Each synonym writes its codon into the shared
/// buffer, reads the window sum, and restores the codon.
pub(crate) struct PottsField {
    theta: Vec<f64>,
    families: Families,
    window: usize,
    total_count: u64,
}

impl PottsField {
    pub(crate) fn bind(model: &PottsModel, codons: &[u8]) -> Result<Self, String> {
        let (theta, families) = prepared(model)?;
        let (_, total_count) =
            sum_log_probability(codons, &theta, model.window, &families, 0..codons.len());
        Ok(Self {
            theta,
            families,
            window: model.window,
            total_count,
        })
    }

    pub(crate) fn delta(&self, codons: &mut [u8], edits: &[(usize, u8)]) -> Result<f64, String> {
        delta_in_place(
            codons,
            edits,
            &self.theta,
            self.window,
            &self.families,
            self.total_count,
        )
    }
}

fn delta_of_slices(
    before: &[u8],
    after: &[u8],
    theta: &[f64],
    window: usize,
    families: &Families,
    total_count: u64,
) -> Result<f64, String> {
    let mut changed = Vec::new();
    for (index, (left, right)) in before.iter().zip(after).enumerate() {
        if left != right {
            changed.push(index);
        }
    }
    if changed.is_empty() {
        return Ok(0.0);
    }
    let affected = affected_positions(before.len(), &changed, window);
    let (before_sum, before_count) =
        sum_log_probability(before, theta, window, families, affected.iter().copied());
    let (after_sum, after_count) =
        sum_log_probability(after, theta, window, families, affected.iter().copied());
    if total_count == 0 || before_count != after_count {
        return Err("no position had a synonymous choice".to_string());
    }
    Ok((after_sum - before_sum) / total_count as f64)
}

fn affected_positions(length: usize, changed: &[usize], window: usize) -> Vec<usize> {
    let mut affected = Vec::new();
    for index in 0..length {
        if changed.iter().any(|&site| index.abs_diff(site) <= window) {
            affected.push(index);
        }
    }
    affected
}

fn restore(codons: &mut [u8], saved: &[(usize, u8)]) {
    for &(index, previous) in saved {
        codons[index] = previous;
    }
}

fn prepared(model: &PottsModel) -> Result<(Vec<f64>, Families), String> {
    let mut theta = model.h.clone();
    theta.extend_from_slice(&model.couplings);
    if theta.len() != 64 + model.window * 4096 {
        return Err("couplings do not match the window".to_string());
    }
    Ok((theta, Families::bacterial()))
}

fn sum_log_probability(
    gene: &[u8],
    theta: &[f64],
    window: usize,
    families: &Families,
    positions: impl IntoIterator<Item = usize>,
) -> (f64, u64) {
    let mut total = 0.0;
    let mut counted = 0u64;
    for index in positions {
        if let Some(value) = log_probability_at(gene, index, theta, window, families) {
            total += value;
            counted += 1;
        }
    }
    (total, counted)
}

fn log_probability_at(
    gene: &[u8],
    index: usize,
    theta: &[f64],
    window: usize,
    families: &Families,
) -> Option<f64> {
    let observed = *gene.get(index)? as usize;
    let family = families.members.get(observed)?;
    let size = family.len();
    if size <= 1 {
        return None;
    }
    let h = &theta[..64];
    let mut field = [0.0; 8];
    let mut observed_slot = 0usize;
    for (slot, &codon) in family.iter().enumerate() {
        let codon = codon as usize;
        if codon == observed {
            observed_slot = slot;
        }
        let mut value = h[codon];
        for distance in 1..=window {
            let matrix = 64 + (distance - 1) * 4096;
            if index + distance < gene.len() {
                let neighbor = gene[index + distance] as usize;
                value += theta[matrix + codon * 64 + neighbor];
            }
            if index >= distance {
                let neighbor = gene[index - distance] as usize;
                value += theta[matrix + neighbor * 64 + codon];
            }
        }
        field[slot] = value;
    }
    let mut peak = f64::NEG_INFINITY;
    for value in field.iter().take(size) {
        peak = peak.max(*value);
    }
    let mut partition = 0.0;
    for value in field.iter().take(size) {
        partition += (value - peak).exp();
    }
    Some(field[observed_slot] - peak - partition.ln())
}

pub fn score_codon_usage(train: &[String], held_out: &[String]) -> Result<UsageBaseline, String> {
    let train_genes = parse_genes(train)?;
    let held_genes = parse_genes(held_out)?;
    let families = Families::bacterial();
    let mut counts = [0u64; 64];
    for gene in &train_genes {
        for &codon in gene {
            counts[codon as usize] += 1;
        }
    }
    let mut total = 0.0;
    let mut scored = 0u64;
    let mut unseen = 0u64;
    for gene in &held_genes {
        for &codon in gene {
            let family = &families.members[codon as usize];
            if family.len() <= 1 {
                continue;
            }
            let family_count: u64 = family.iter().map(|&member| counts[member as usize]).sum();
            let own = counts[codon as usize];
            if family_count == 0 || own == 0 {
                unseen += 1;
                continue;
            }
            total += (own as f64 / family_count as f64).ln();
            scored += 1;
        }
    }
    if scored == 0 {
        return Err("the codon-usage baseline scored no held-out position".to_string());
    }
    Ok(UsageBaseline {
        nats_per_codon: total / scored as f64,
        n_scored: scored,
        n_unseen: unseen,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_context_bias_beats_the_codon_usage_baseline() {
        // Leucine after AAA is TTA. Leucine after anything else is CTG.
        // The independent baseline cannot see the neighbour.
        let mut train = Vec::new();
        let mut held = Vec::new();
        for copy in 0..40 {
            let mut sequence = String::new();
            for site in 0..24 {
                let leucine = if site % 2 == 0 { "TTA" } else { "CTG" };
                let lysine = if site % 2 == 0 { "AAA" } else { "AAG" };
                sequence.push_str(lysine);
                sequence.push_str(leucine);
            }
            if copy < 30 {
                train.push(sequence);
            } else {
                held.push(sequence);
            }
        }
        let report = fit_potts(&train, 1, 1e-4, 1e-4, 12).expect("fit");
        let (model_nats, _) = score_naturalness(&held, &report.model).expect("score");
        let baseline = score_codon_usage(&train, &held).expect("baseline");
        assert!(
            model_nats > baseline.nats_per_codon + 0.05,
            "model {model_nats} did not beat baseline {}",
            baseline.nats_per_codon
        );
    }

    #[test]
    fn the_gauge_sums_are_zero_inside_each_family() {
        let sequences = vec!["ATGAAACGTATGAAACGTATGAAACGT".to_string(); 8];
        let report = fit_potts(&sequences, 1, 1e-3, 1e-3, 4).expect("fit");
        let families = Families::bacterial();
        for index in 0..64 {
            let family = &families.members[index];
            if !is_family_leader(index, family) {
                continue;
            }
            let sum: f64 = family.iter().map(|&c| report.model.h[c as usize]).sum();
            assert!(sum.abs() < 1e-8, "family {index} sums to {sum}");
        }
    }

    #[test]
    fn the_local_delta_matches_a_rescore() {
        let gct = Codon::from_bytes(b"GCT").unwrap().index() as usize;
        let mut model = PottsModel {
            window: 1,
            h: vec![0.0; 64],
            couplings: vec![0.0; 4096],
        };
        model.h[gct] = 2.0;
        // A neighbor coupling so the edit is not a single-site field.
        let atg = Codon::from_bytes(b"ATG").unwrap().index() as usize;
        model.couplings[atg * 64 + gct] = 0.5;
        let before = parse_genes(&["ATGGCTGCTTAA".to_string()]).unwrap();
        let after = parse_genes(&["ATGGCCGCTTAA".to_string()]).unwrap();
        let delta = mean_log_probability_delta(&before[0], &after[0], &model).unwrap();
        let expected = mean_log_probability(&after[0], &model).unwrap()
            - mean_log_probability(&before[0], &model).unwrap();
        assert!(
            (delta - expected).abs() < 1e-12,
            "delta {delta} rescore {expected}"
        );
        assert!(delta.abs() > 1e-6, "the edit did not move the score");
    }

    #[test]
    fn inplace_delta_is_bit_identical_to_the_slice_delta() {
        let families = Families::bacterial();
        let mut state = 0xA5A5_F00D_1234_5678u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for case in 0..48 {
            let window = 1 + (next() as usize % 3);
            let mut model = PottsModel {
                window,
                h: vec![0.0; 64],
                couplings: vec![0.0; window * 4096],
            };
            for _ in 0..24 {
                let index = next() as usize % 64;
                model.h[index] = (next() % 200) as f64 / 25.0 - 4.0;
            }
            for _ in 0..80 {
                let index = next() as usize % model.couplings.len();
                model.couplings[index] = (next() % 200) as f64 / 40.0 - 2.5;
            }
            // A coupling at exactly `window` so a short window is a different sum.
            let anchor = Codon::from_bytes(b"GCT").unwrap().index() as usize;
            let neighbor = Codon::from_bytes(b"CGT").unwrap().index() as usize;
            let matrix = (window - 1) * 4096;
            model.couplings[matrix + anchor * 64 + neighbor] = 1.25;
            let length = 12;
            let mut codons = Vec::with_capacity(length);
            for position in 0..length {
                codons.push(if position % 5 == 0 {
                    anchor as u8
                } else if position % 5 == window.min(4) {
                    neighbor as u8
                } else {
                    (next() % 64) as u8
                });
            }
            let mut edits = Vec::new();
            for _ in 0..2 {
                let index = next() as usize % length;
                let family = &families.members[codons[index] as usize];
                if family.len() <= 1 {
                    continue;
                }
                let codon = family[next() as usize % family.len()];
                edits.push((index, codon));
            }
            if edits.is_empty() {
                continue;
            }
            let mut after = codons.clone();
            for &(index, codon) in &edits {
                after[index] = codon;
            }
            let (theta, prepared_families) = prepared(&model).unwrap();
            let (_, total_count) =
                sum_log_probability(&codons, &theta, window, &prepared_families, 0..codons.len());
            let expected = mean_log_probability_delta(&codons, &after, &model).unwrap();
            let mut working = codons.clone();
            let got = delta_in_place(
                &mut working,
                &edits,
                &theta,
                window,
                &prepared_families,
                total_count,
            )
            .unwrap();
            assert_eq!(working, codons, "case {case} left the codon buffer edited");
            assert_eq!(
                got.to_bits(),
                expected.to_bits(),
                "case {case} window {window} inplace {got} slice {expected}"
            );
        }
    }
}
