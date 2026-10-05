//! Watson–Crick duplex ΔH, ΔS, ΔG°37 and Tm from SantaLucia 1998 Eqs. 2 and 3.
//!
//! Mismatch nearest-neighbour parameters are not in Table 2, so a mismatch
//! ends a duplex. Dangling-end parameters are not in Table 2 and are not
//! added.

use crate::parameters::{ParameterError, Parameters, Thermo};

#[derive(Debug, Clone, PartialEq)]
pub struct Duplex {
    pub delta_h_kcal: f64,
    pub delta_s_cal: f64,
    pub delta_g37_kcal: f64,
    pub tm_celsius: f64,
    pub self_complementary: bool,
}

pub fn complement(base: u8) -> Option<u8> {
    match base {
        b'A' => Some(b'T'),
        b'T' => Some(b'A'),
        b'C' => Some(b'G'),
        b'G' => Some(b'C'),
        _ => None,
    }
}

pub fn reverse_complement(seq: &str) -> Result<String, ParameterError> {
    let mut out = Vec::with_capacity(seq.len());
    for base in seq.bytes().rev() {
        out.push(complement(base).ok_or_else(|| {
            ParameterError(format!(
                "CF-PRIM-006 base {} is not A, C, G or T",
                base as char
            ))
        })?);
    }
    Ok(String::from_utf8(out).expect("ascii"))
}

pub fn is_self_complementary(seq: &str) -> bool {
    reverse_complement(seq).ok().as_deref() == Some(seq)
}

fn stack_key(left: u8, right: u8) -> Result<String, ParameterError> {
    let bottom_left = complement(left).ok_or_else(|| ParameterError("bad base".into()))?;
    let bottom_right = complement(right).ok_or_else(|| ParameterError("bad base".into()))?;
    let forward = format!(
        "{}{}/{}{}",
        left as char, right as char, bottom_left as char, bottom_right as char
    );
    Ok(forward)
}

fn lookup(
    parameters: &Parameters,
    left: u8,
    right: u8,
) -> Result<&Thermo, ParameterError> {
    let forward = stack_key(left, right)?;
    if let Some(thermo) = parameters.stacks.get(&forward) {
        return Ok(thermo);
    }
    let flipped_left = complement(right).ok_or_else(|| ParameterError("bad base".into()))?;
    let flipped_right = complement(left).ok_or_else(|| ParameterError("bad base".into()))?;
    let reverse = stack_key(flipped_left, flipped_right)?;
    parameters
        .stacks
        .get(&reverse)
        .ok_or_else(|| ParameterError(format!("CF-PRIM-005 no stack for {forward} or {reverse}")))
}

fn initiation(parameters: &Parameters, base: u8) -> Thermo {
    if base == b'G' || base == b'C' {
        parameters.init_gc.clone()
    } else {
        parameters.init_at.clone()
    }
}

/// Perfect duplex of `top` with its Watson–Crick complement.
pub fn perfect_duplex(parameters: &Parameters, top: &str) -> Result<Duplex, ParameterError> {
    let bytes = top.as_bytes();
    if bytes.len() < 2 {
        return Err(ParameterError(
            "CF-PRIM-006 a duplex needs at least two paired bases".into(),
        ));
    }
    let mut delta_h = 0.0;
    let mut delta_s = 0.0;
    for pair in bytes.windows(2) {
        let thermo = lookup(parameters, pair[0], pair[1])?;
        delta_h += thermo.delta_h_kcal;
        delta_s += thermo.delta_s_cal;
    }
    for end in [bytes[0], bytes[bytes.len() - 1]] {
        let init = initiation(parameters, end);
        delta_h += init.delta_h_kcal;
        delta_s += init.delta_s_cal;
    }
    let self_complementary = is_self_complementary(top);
    if self_complementary {
        delta_h += parameters.symmetry.delta_h_kcal;
        delta_s += parameters.symmetry.delta_s_cal;
    }
    // Eq. 7 / Eq. 8. N = length − 1 (SantaLucia 1998: an 8-bp duplex has N = 7).
    let sodium = parameters.sodium_molar;
    if (sodium - 1.0).abs() > 1e-12 {
        let n = (bytes.len() - 1) as f64;
        let delta_g_salt = parameters.salt_dg37_kcal_per_ln_na * n * sodium.ln();
        let delta_s_salt_cal = -1000.0 * delta_g_salt / parameters.temperature_k;
        delta_s += delta_s_salt_cal;
    }
    let delta_g37 = delta_h - parameters.temperature_k * (delta_s / 1000.0);
    let ct = parameters.strand_concentration_molar;
    let concentration = if self_complementary { ct } else { ct / 4.0 };
    let tm_k = (delta_h * 1000.0) / (delta_s + parameters.gas_constant_cal * concentration.ln());
    Ok(Duplex {
        delta_h_kcal: delta_h,
        delta_s_cal: delta_s,
        delta_g37_kcal: delta_g37,
        tm_celsius: tm_k - 273.15,
        self_complementary,
    })
}

pub fn gc_percent(seq: &str) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .bytes()
        .filter(|base| *base == b'G' || *base == b'C')
        .count();
    100.0 * gc as f64 / seq.len() as f64
}

pub fn longest_homopolymer(seq: &str) -> usize {
    let mut best = 0;
    let mut run = 0;
    let mut previous = 0;
    for base in seq.bytes() {
        if base == previous {
            run += 1;
        } else {
            run = 1;
            previous = base;
        }
        best = best.max(run);
    }
    best
}

pub fn gc_clamp(seq: &str) -> bool {
    matches!(seq.as_bytes().last(), Some(b'G' | b'C'))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DgZone {
    Accept,
    Flag,
    Refuse,
}

/// IDT's published cutoffs, as used by this suite: below −9 kcal/mol refuses,
/// below −4 kcal/mol flags. A value of exactly −4 or −9 stays on the milder side.
pub fn dg_zone(delta_g: f64) -> DgZone {
    if delta_g < -9.0 {
        DgZone::Refuse
    } else if delta_g < -4.0 {
        DgZone::Flag
    } else {
        DgZone::Accept
    }
}
