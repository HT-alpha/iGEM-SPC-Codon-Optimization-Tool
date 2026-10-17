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

fn geometric_cai(cds: &str, weights: &[f64; 64]) -> f64 {
    let sequence = CodonSequence::from_cds(cds, 11).unwrap();
    let mut sum = 0.0;
    let mut counted = 0usize;
    for (index, codon) in sequence.codons().iter().copied().enumerate() {
        if sequence.table().is_stop(codon) || sequence.synonyms_at(index).len() <= 1 {
            continue;
        }
        let weight = weights[codon.index_usize()];
        if weight == 0.0 {
            return 0.0;
        }
        sum += weight.ln();
        counted += 1;
    }
    if counted == 0 {
        1.0
    } else {
        (sum / counted as f64).exp()
    }
}

fn assert_repaired(name: &str) {
    let table = weights();
    let repaired = assert_repaired_with(name, &table);
    // Ribosomal weights. The digest pins the sequence after the conflict climb
    // and the four-start keep (`DEP-SRCH-021`, `DEP-SRCH-022`). The CAI floor
    // is that sequence's CAI, which is at least the v2.9.3 seed-7 assembly
    // number. waaS's floor is above 0.9220206722508919: stopping at the first
    // feasible start (the weight-ordered one) returns 0.9212217911296222.
    let (expected, floor) = match name {
        "elp_vpgvg.fa" => (
            "848d0378fa464d9a65a26c663646614a82e082abf4ed14e5d3344497eec5c910",
            0.4642192192077335,
        ),
        "yagF.fa" => (
            "773523e6566cfdcf0996ac51123ffd332c78e5a4703ed40062949e20644cb514",
            0.9453398274862173,
        ),
        "waaS.fa" => (
            "d227b303a3c477cea7f5ed35c3a0fbcb783cb9169fa38a3737e0fe7699419f16",
            0.9259163095438141,
        ),
        "yeeJ.fa" => (
            "72bd74389f99f3aba3dfcb5422c05c91db606adaacb65f67a5a80f881dc6e540",
            0.9727484991594185,
        ),
        "rpsL.fa" => (
            "0624d04f50821da0804a0cf91fa59b9a6bde2ba8ef6156bb79e9680ba09eedfc",
            0.9866866273575934,
        ),
        other => panic!("no pinned seed for {other}"),
    };
    let mut hasher = Sha256::new();
    hasher.update(repaired.as_bytes());
    let digest = format!("{:x}", hasher.finalize());
    let cai = geometric_cai(&repaired, &table);
    eprintln!("{name}: cai={cai:.16} sha256={digest}");
    assert_eq!(digest, expected, "{name} seed changed");
    assert!(cai + 1e-12 >= floor, "{name} CAI {cai} fell below {floor}");
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

#[test]
fn seed_repairs_rpsl() {
    assert_repaired("rpsL.fa");
}
