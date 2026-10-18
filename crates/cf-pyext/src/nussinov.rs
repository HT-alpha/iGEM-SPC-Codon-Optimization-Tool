//! Opt-in Nussinov maximum matching.
//!
//! Nussinov and Jacobson, PNAS 1980, 77:6309–6313. A pair scores 1 when it is
//! A–U, G–C or G–U. The hairpin must leave at least three unpaired bases.
//! There is no stacking term and no Turner parameter. `temperature_c` is
//! accepted by the Python protocol and does not enter this score.
//!
//! The returned value is an unpaired fraction in [0, 1], one per window of
//! `window_nt`, the same length ViennaRNA's binding returns. It is not an
//! unpaired probability and not a free energy. `DEP-STRU-006`.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyfunction]
pub fn nussinov_accessibility(
    sequence: &str,
    window_nt: usize,
    span_nt: usize,
) -> PyResult<Vec<f64>> {
    if window_nt < 1 {
        return Err(PyValueError::new_err(
            "CF-STRU-004 window_nt must be at least 1",
        ));
    }
    let seq = normalize(sequence)?;
    if window_nt > seq.len() {
        return Err(PyValueError::new_err(format!(
            "CF-STRU-004 window_nt {window_nt} exceeds the sequence length {}",
            seq.len()
        )));
    }
    if span_nt < window_nt {
        return Err(PyValueError::new_err(format!(
            "CF-STRU-004 span_nt {span_nt} is smaller than window_nt {window_nt}"
        )));
    }
    let n = seq.len();
    let count = n - window_nt + 1;
    if n < span_nt {
        let paired = matching(&seq);
        return Ok((0..count)
            .map(|start| unpaired_fraction(&paired, start, window_nt))
            .collect());
    }
    let mut cache: Vec<Option<Vec<bool>>> = vec![None; n - span_nt + 1];
    Ok((0..count)
        .map(|start| {
            let span_start = span_start_for(start, start + window_nt, span_nt, n);
            let paired = cache[span_start]
                .get_or_insert_with(|| matching(&seq[span_start..span_start + span_nt]));
            unpaired_fraction(paired, start - span_start, window_nt)
        })
        .collect())
}

fn span_start_for(window_start: usize, window_end: usize, span_nt: usize, n: usize) -> usize {
    let mid = window_start + (window_end - window_start) / 2;
    let ideal = mid.saturating_sub(span_nt / 2);
    ideal.min(n - span_nt)
}

fn unpaired_fraction(paired: &[bool], start: usize, window_nt: usize) -> f64 {
    let unpaired = paired[start..start + window_nt]
        .iter()
        .filter(|paired| !**paired)
        .count();
    unpaired as f64 / window_nt as f64
}

fn normalize(sequence: &str) -> PyResult<Vec<u8>> {
    let mut out = Vec::with_capacity(sequence.len());
    for (index, base) in sequence.chars().enumerate() {
        let byte = match base.to_ascii_uppercase() {
            'A' => b'A',
            'C' => b'C',
            'G' => b'G',
            'T' | 'U' => b'U',
            other => {
                return Err(PyValueError::new_err(format!(
                    "CF-STRU-004 base {other} at position {index} is not A, C, G, T or U"
                )));
            }
        };
        out.push(byte);
    }
    if out.is_empty() {
        return Err(PyValueError::new_err("CF-STRU-004 sequence is empty"));
    }
    Ok(out)
}

fn pairs(left: u8, right: u8) -> bool {
    matches!(
        (left, right),
        (b'A', b'U') | (b'U', b'A') | (b'G', b'C') | (b'C', b'G') | (b'G', b'U') | (b'U', b'G')
    )
}

/// Maximum-cardinality matching. On a tie the leftmost partner is kept.
fn matching(seq: &[u8]) -> Vec<bool> {
    let n = seq.len();
    if n == 0 {
        return Vec::new();
    }
    let mut best = vec![vec![0u16; n]; n];
    let mut partner: Vec<Vec<Option<u16>>> = vec![vec![None; n]; n];
    for len in 1..n {
        for i in 0..n - len {
            let j = i + len;
            let mut score = best[i][j - 1];
            let mut choice = None;
            if j >= i + 4 {
                for k in i..j {
                    if j - k < 4 || !pairs(seq[k], seq[j]) {
                        continue;
                    }
                    let left = if k > i { best[i][k - 1] } else { 0 };
                    let between = if k + 1 < j { best[k + 1][j - 1] } else { 0 };
                    let candidate = left + between + 1;
                    if candidate > score {
                        score = candidate;
                        choice = Some(k as u16);
                    }
                }
            }
            best[i][j] = score;
            partner[i][j] = choice;
        }
    }
    let mut paired = vec![false; n];
    trace(&partner, &mut paired, 0, n - 1);
    paired
}

fn trace(partner: &[Vec<Option<u16>>], paired: &mut [bool], i: usize, j: usize) {
    if i >= j {
        return;
    }
    match partner[i][j] {
        None => trace(partner, paired, i, j - 1),
        Some(k_raw) => {
            let k = k_raw as usize;
            paired[k] = true;
            paired[j] = true;
            if k > i {
                trace(partner, paired, i, k - 1);
            }
            if k + 1 < j {
                trace(partner, paired, k + 1, j - 1);
            }
        }
    }
}
