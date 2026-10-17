//! `burden.growth_cost` moves when the codon mix moves, and its delta matches
//! a rescore.
//!
//! A constant `g` would satisfy the unit and the direction and still be useless.
//! The piled assignment and the host mix are the two designs the term exists
//! to separate.

use cf_burden::{growth_cost, power_mean, GrowthCostTerm, AGGREGATION_P, SATURATION_N};
use cf_landscape_api::{
    verify_capability, Budget, BudgetTier, EvalContext, HostContext, TermProvider,
};
use cf_types::{CodonSequence, Digest, HostId};
use serde_json::json;

fn profile(weight: f64) -> serde_json::Value {
    let mut table = serde_json::Map::new();
    for codon in cf_types::Codon::all() {
        table.insert(
            codon.to_string_upper(),
            json!({"relative_adaptiveness": weight}),
        );
    }
    json!({"codon_usage": {"subset": "highly_expressed", "table": table}})
}

fn host() -> HostContext {
    HostContext::new(
        HostId::new("escherichia_coli:511145/BL21").unwrap(),
        Digest::of_bytes(b"burden-penalty-fixture"),
    )
}

fn score(term: &GrowthCostTerm, cds: &str) -> f64 {
    let sequence = CodonSequence::from_cds(cds, 11).unwrap();
    let host = host();
    let ctx = EvalContext::new(&sequence, &host, Budget::new(BudgetTier::Standard), 0);
    term.score(&ctx).unwrap()
}

#[test]
fn power_mean_of_equals_is_that_value_and_saturation_at_one_is_one_half() {
    assert!((power_mean(&[1.0, 1.0, 1.0, 1.0], AGGREGATION_P) - 1.0).abs() < 1e-12);
    assert!((growth_cost(1.0, SATURATION_N, 1.0) - 0.5).abs() < 1e-12);
}

#[test]
fn a_piled_assignment_costs_more_than_the_host_mix() {
    let term = GrowthCostTerm::from_host_profile(&profile(1.0)).unwrap();
    // Eight alanines, one codon versus two of each synonym. Met and the stop
    // are single-codon or stop and do not enter the penalty.
    let piled = score(&term, "ATGGCTGCTGCTGCTGCTGCTGCTGCTGCTTAA");
    let mixed = score(&term, "ATGGCTGCCGCAGCGGCTGCCGCAGCGTAA");
    assert!(
        piled > mixed,
        "piling alanine onto GCT scored {piled}, the even mix scored {mixed}"
    );
    assert!(
        (mixed - 0.5).abs() < 1e-9,
        "an even mix against a flat table is u = 1 and g = 0.5, got {mixed}"
    );
}

#[test]
fn the_delta_matches_a_rescore() {
    let term = GrowthCostTerm::from_host_profile(&profile(1.0)).unwrap();
    let sequence = CodonSequence::from_cds("ATGGCTGCTGCTGCTGCTGCTGCTGCTTAA", 11).unwrap();
    let host = host();
    let ctx = EvalContext::new(&sequence, &host, Budget::new(BudgetTier::Standard), 0);
    verify_capability(&term, &ctx).unwrap();
}

#[test]
fn a_profile_with_no_supply_prior_is_refused() {
    let error = GrowthCostTerm::from_host_profile(&json!({})).unwrap_err();
    assert_eq!(error.code, "CF-BURD-004");
}
