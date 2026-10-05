//! Two-state hairpin and dimer ΔG.
//!
//! A hairpin is the most stable Watson–Crick stem of at least one stack, with
//! a loop of at least three unpaired bases. Table 2 has no hairpin-loop
//! initiation, so none is added. A dimer is the most stable contiguous
//! Watson–Crick alignment of the two oligos, including overhangs. The 3′
//! figure is that minimum restricted to alignments that pair a 3′ terminal base.

use crate::duplex::{complement, perfect_duplex, Duplex};
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
    let mut best_3p_align = "no 3' Watson-Crick stack".to_string();
    for start_a in 0..a.len() {
        for start_b in 0..b.len() {
            let mut length = 0;
            while start_a + length < a.len()
                && start_b + length < b.len()
                && complement(a[start_a + length]) == Some(b[start_b + length])
            {
                length += 1;
            }
            if length < 2 {
                continue;
            }
            let slice = &left[start_a..start_a + length];
            let duplex: Duplex = perfect_duplex(parameters, slice)?;
            let alignment = alignment_text(left, right, start_a, start_b, length);
            if duplex.delta_g37_kcal < best {
                best = duplex.delta_g37_kcal;
                best_align = alignment.clone();
            }
            let pairs_left_3p = start_a + length == a.len();
            let pairs_right_3p = start_b == 0;
            if (pairs_left_3p || pairs_right_3p) && duplex.delta_g37_kcal < best_3p {
                best_3p = duplex.delta_g37_kcal;
                best_3p_align = alignment;
            }
        }
    }
    let _ = best_3p_align;
    Ok(Dimer {
        delta_g_kcal: best,
        three_prime_delta_g_kcal: best_3p,
        alignment: best_align,
    })
}

fn alignment_text(
    left: &str,
    right: &str,
    start_a: usize,
    start_b_rev: usize,
    length: usize,
) -> String {
    let right_rev: String = right.chars().rev().collect();
    let pad = start_a.max(start_b_rev);
    let mut top = " ".repeat(pad - start_a);
    top.push_str(left);
    let mut bottom = " ".repeat(pad - start_b_rev);
    bottom.push_str(&right_rev);
    let mut marks = " ".repeat(pad);
    marks.push_str(&"|".repeat(length));
    format!("{top}\n{marks}\n{bottom}")
}
