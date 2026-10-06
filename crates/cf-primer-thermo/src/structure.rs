//! Two-state hairpin and dimer ΔG.
//!
//! A hairpin is the most stable Watson–Crick stem of at least one stack, with
//! a loop of at least three unpaired bases. Table 2 has no hairpin-loop
//! initiation, so none is added. A dimer is the most stable contiguous
//! Watson–Crick alignment of the two oligos, including overhangs. The 3′
//! figure is that minimum restricted to alignments that pair a 3′ terminal base.

use crate::duplex::{complement, perfect_duplex};
use crate::parameters::{ParameterError, Parameters};

#[derive(Debug, Clone, PartialEq)]
pub struct Fold {
    pub delta_g_kcal: f64,
    pub diagram: String,
}

pub fn hairpin(parameters: &Parameters, seq: &str) -> Result<Fold, ParameterError> {
    let bytes = seq.as_bytes();
    let mut best: Option<(f64, String)> = None;
    if bytes.len() >= 2 + 3 + 2 {
        for stem in 2..=(bytes.len() - 3) / 2 {
            let loop_start = stem;
            let loop_end = bytes.len() - stem;
            if loop_end < loop_start + 3 {
                continue;
            }
            let paired = (0..stem)
                .all(|index| complement(bytes[index]) == Some(bytes[bytes.len() - 1 - index]));
            if !paired {
                continue;
            }
            let stem_seq = &seq[..stem];
            let duplex = perfect_duplex(parameters, stem_seq)?;
            let diagram = hairpin_diagram(stem_seq, loop_end - loop_start);
            match &best {
                Some((delta_g, _)) if duplex.delta_g37_kcal >= *delta_g => {}
                _ => best = Some((duplex.delta_g37_kcal, diagram)),
            }
        }
    }
    Ok(match best {
        Some((delta_g_kcal, diagram)) => Fold {
            delta_g_kcal,
            diagram,
        },
        None => Fold {
            delta_g_kcal: 0.0,
            diagram: "no Watson-Crick stem of two or more pairs".into(),
        },
    })
}

fn hairpin_diagram(stem: &str, loop_len: usize) -> String {
    let bottom: String = stem
        .bytes()
        .rev()
        .map(|base| complement(base).unwrap() as char)
        .collect();
    let marks = "|".repeat(stem.len());
    format!("5'-{stem}\n   {marks}\n3'-{bottom}  loop {loop_len} nt")
}

#[derive(Debug, Clone, PartialEq)]
pub struct Dimer {
    pub delta_g_kcal: f64,
    pub three_prime_delta_g_kcal: f64,
    pub alignment: String,
}

pub fn dimer(parameters: &Parameters, left: &str, right: &str) -> Result<Dimer, ParameterError> {
    let a = left.as_bytes();
    let b: Vec<u8> = right.as_bytes().iter().rev().copied().collect();
    let mut best = 0.0;
    let mut best_align = "no Watson-Crick stack".to_string();
    let mut best_3p = 0.0;
    // One register is one offset between the two strands. Every Watson-Crick
    // stack in that register belongs to the same duplex, so initiation is
    // charged once. Scoring each contiguous run as its own duplex charges
    // initiation again and drops the other helix.
    let a_len = a.len() as isize;
    let b_len = b.len() as isize;
    for offset in -(b_len - 1)..a_len {
        let mut paired = vec![false; a.len()];
        for i in 0..a.len() {
            let j = i as isize - offset;
            if (0..b_len).contains(&j) && complement(a[i]) == Some(b[j as usize]) {
                paired[i] = true;
            }
        }
        let mut runs: Vec<(usize, usize)> = Vec::new();
        let mut index = 0;
        while index < paired.len() {
            if !paired[index] {
                index += 1;
                continue;
            }
            let start = index;
            while index < paired.len() && paired[index] {
                index += 1;
            }
            if index - start >= 2 {
                runs.push((start, index));
            }
        }
        if runs.is_empty() {
            continue;
        }
        let mut delta_h = 0.0;
        let mut delta_s = 0.0;
        for (start, end) in &runs {
            let duplex = perfect_duplex(parameters, &left[*start..*end])?;
            // perfect_duplex includes initiation and, for a palindromic
            // slice, the symmetry term. Subtract those and re-apply once
            // for the whole register below.
            let (init_h, init_s) = initiation_totals(
                parameters,
                left.as_bytes()[*start],
                left.as_bytes()[*end - 1],
            );
            delta_h += duplex.delta_h_kcal - init_h;
            delta_s += duplex.delta_s_cal - init_s;
            if duplex.self_complementary {
                delta_h -= parameters.symmetry.delta_h_kcal;
                delta_s -= parameters.symmetry.delta_s_cal;
            }
        }
        let (init_h, init_s) = initiation_totals(
            parameters,
            left.as_bytes()[runs[0].0],
            left.as_bytes()[runs[runs.len() - 1].1 - 1],
        );
        delta_h += init_h;
        delta_s += init_s;
        if runs.len() == 1 {
            let only = &left[runs[0].0..runs[0].1];
            if crate::duplex::is_self_complementary(only) {
                delta_h += parameters.symmetry.delta_h_kcal;
                delta_s += parameters.symmetry.delta_s_cal;
            }
        }
        let delta_g = delta_h - parameters.temperature_k * (delta_s / 1000.0);
        let alignment = register_alignment(left, &b, offset, &paired);
        if delta_g < best {
            best = delta_g;
            best_align = alignment.clone();
        }
        let pairs_left_3p = paired[a.len() - 1];
        let pairs_right_3p = offset <= 0 && ((-offset) as usize) < b.len() && {
            // b[0] is the 3' end of `right`. It is paired when some i has i - offset == 0.
            paired
                .iter()
                .enumerate()
                .any(|(i, is_pair)| *is_pair && i as isize - offset == 0)
        };
        if (pairs_left_3p || pairs_right_3p) && delta_g < best_3p {
            best_3p = delta_g;
        }
        let _ = alignment;
    }
    Ok(Dimer {
        delta_g_kcal: best,
        three_prime_delta_g_kcal: best_3p,
        alignment: best_align,
    })
}

fn initiation_totals(parameters: &Parameters, five: u8, three: u8) -> (f64, f64) {
    let five_prime = if five == b'G' || five == b'C' {
        &parameters.init_gc
    } else {
        &parameters.init_at
    };
    let three_prime = if three == b'G' || three == b'C' {
        &parameters.init_gc
    } else {
        &parameters.init_at
    };
    (
        five_prime.delta_h_kcal + three_prime.delta_h_kcal,
        five_prime.delta_s_cal + three_prime.delta_s_cal,
    )
}

fn register_alignment(left: &str, reversed_right: &[u8], offset: isize, paired: &[bool]) -> String {
    let mut marks = String::new();
    let mut bottom = String::new();
    for (index, base_paired) in paired.iter().enumerate() {
        let j = index as isize - offset;
        if (0..reversed_right.len() as isize).contains(&j) {
            bottom.push(reversed_right[j as usize] as char);
            marks.push(if *base_paired { '|' } else { '.' });
        }
    }
    format!("{left}\n{marks}\n{bottom}")
}
