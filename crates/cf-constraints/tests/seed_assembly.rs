//! Seed resolution against the validation assembly set.
//!
//! These are the four sequences the v2.5.0 note showed failing CF-SRCH-006.
//! The test is the acceptance check for DEP-SRCH-015: a synonymous sequence
//! that the checker accepts, or a residual that names the rule still broken.

use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use cf_constraints::{resolve_feasible_seed, ConstraintChecker, ConstraintSet};
use cf_types::codon::Codon;
use cf_types::seqrec::CodonSequence;
use sha2::{Digest, Sha256};

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn fasta(name: &str) -> String {
    let text = fs::read_to_string(
        repo_root()
            .join("benchmarks/validation/sequences")
            .join(name),
    )
    .unwrap_or_else(|err| panic!("read {name}: {err}"));
    text.lines()
        .filter(|line| !line.starts_with('>'))
        .collect::<String>()
        .chars()
        .filter(|ch| !ch.is_whitespace())
        .collect::<String>()
        .to_ascii_uppercase()
}

fn synthetic_weights() -> [f64; 64] {
    // The orchestrator fixture's non-uniform table. Tie-breaks differ from the
    // ribosomal weights, which is how a one-codon search got stuck on the
    // elastin repeat.
    let bases = *b"TCAG";
    let mut weights = [1.0; 64];
    for (i, a) in bases.iter().enumerate() {
        for (j, b) in bases.iter().enumerate() {
            for (k, c) in bases.iter().enumerate() {
                let codon = Codon::from_bytes(&[*a, *b, *c]).unwrap();
                let index = i * 16 + j * 4 + k;
                let raw = 0.05 + 0.95 * ((index * 17) % 23) as f64 / 22.0;
                weights[codon.index_usize()] = (raw * 10_000.0).round() / 10_000.0;
            }
        }
    }
    weights
}

fn weights() -> [f64; 64] {
    let text =
        fs::read_to_string(repo_root().join("examples/mg1655_ribosomal_usage.json")).unwrap();
    let value: serde_json::Value = serde_json::from_str(&text).unwrap();
    let mut weights = [1.0; 64];
    let table = value["codon_usage"]["table"].as_object().unwrap();
    for (codon, entry) in table {
        let parsed = Codon::parse(codon).unwrap();
        weights[parsed.index_usize()] = entry["relative_adaptiveness"].as_f64().unwrap();
    }
    weights
}

fn checker() -> ConstraintChecker {
    let text =
        fs::read_to_string(repo_root().join("benchmarks/validation/constraints_assembly.json"))
            .unwrap();
    let value: serde_json::Value = serde_json::from_str(&text).unwrap();
    let set = ConstraintSet::from_json(&value).unwrap();
    ConstraintChecker::compile(&set).unwrap()
}

fn assert_repaired_with(name: &str, weights: &[f64; 64]) -> String {
    let cds = fasta(name);
    let sequence = CodonSequence::from_cds(&cds, 11).unwrap();
    let compiled = checker();
    let started = Instant::now();
    let report = resolve_feasible_seed(&sequence, &compiled, weights).unwrap();
    eprintln!(
        "{name}: feasible={} edits={} cycles={} elapsed={:.3}s residual={:?}",
        report.feasible,
        report.edits,
        report.cycles,
        started.elapsed().as_secs_f64(),
        report.residual
    );
    assert!(
        report.feasible,
        "{name} stayed infeasible: {:?}",
        report.residual
    );
    assert_ne!(report.sequence.cds_string(), cds);
    assert_eq!(report.sequence.protein(), sequence.protein());
    report.sequence.cds_string()
}

fn assert_repaired(name: &str) {
    let repaired = assert_repaired_with(name, &weights());
    // Ribosomal weights. Feasibility still holds if an accept/reject flips and
    // a different synonymous seed is returned, so the sequence itself is pinned.
    let expected = match name {
        "elp_vpgvg.fa" => "067af67421dbeb2dc933198ea718fd2d9b7e1b11a2e6c9e725dbf0ac2f5c4da9",
        "yagF.fa" => "720ee3bdcc76e18b75afadf7a801238941aa1b79052f0574ab12a5033e45cab9",
        "waaS.fa" => "11cdd2fb04a48e29ff284107103ddab0c83e741341b7bf9926fff06d8b90627d",
        "yeeJ.fa" => "686adf6a7ed03ac848270c94afbb9f7bbfd7873dbf926bc03390001d80eb9b4c",
        other => panic!("no pinned seed for {other}"),
    };
    let mut hasher = Sha256::new();
    hasher.update(repaired.as_bytes());
    let digest = format!("{:x}", hasher.finalize());
    assert_eq!(digest, expected, "{name} seed changed");
}

#[test]
fn seed_repairs_the_elastin_repeat() {
    assert_repaired("elp_vpgvg.fa");
}

#[test]
fn seed_repairs_elastin_under_the_orchestrator_weight_table() {
    assert_repaired_with("elp_vpgvg.fa", &synthetic_weights());
}

#[test]
fn seed_repairs_yagf() {
    assert_repaired("yagF.fa");
}

#[test]
fn seed_repairs_waas() {
    assert_repaired("waaS.fa");
}

#[test]
fn seed_repairs_yeej() {
    assert_repaired("yeeJ.fa");
}
