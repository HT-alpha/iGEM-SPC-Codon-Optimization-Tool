//! `TST-BURD-005`/`TST-BURD-006` — capacity estimators.
//!
//! Ticket P2-13. Spec: `planning/03-algorithms/05-burden-penalty-and-capacity.md` §2.
//!
//! The sequencing audit found `TST-BURD-005` as written tests *honesty* — unit
//! `dimensionless`, `calibrated: false`, `low_confidence_host` attached — and
//! that an estimator returning a plausible constant would pass every one of
//! those checks. That is a real gap: the semantics would be impeccable and the
//! number meaningless.
//!
//! So `TST-BURD-005` keeps the semantic checks, and `TST-BURD-006` adds the one
//! the audit asked for: the estimate must **move with** the host's real
//! provenance data, following P2-12's precedent of testing against
//! `trna.gene_copy_number` and `tai.s_values` rather than a synthetic stand-in
//! that would pass regardless.

use std::collections::BTreeMap;

use cf_burden::{
    decoders, fit_capacity_tier1, fit_capacity_tier2, fit_capacity_tier3, DataCompleteness,
    DecodingWeights, ThetaUnit,
};
use cf_types::Codon;

/// Native demand derived from real provenance inputs, exactly as the
/// orchestrator will build it.
///
/// Takes `gene_copy_number` per anticodon and `s_values` for wobble decoding,
/// so varying either changes the result — which is the property under test.
fn native_demand_from_provenance(
    codon_counts: &[(&str, f64)],
    s_values: &BTreeMap<String, f64>,
) -> BTreeMap<[u8; 3], f64> {
    let weights = DecodingWeights::from_s_values(s_values.clone()).expect("s_values");
    let mut demand: BTreeMap<[u8; 3], f64> = BTreeMap::new();
    for (codon_text, count) in codon_counts {
        let codon = Codon::parse(codon_text).expect("codon");
        for decoder in decoders(codon) {
            let weight = weights.weight(codon, &decoder).expect("weight");
            if weight > 0.0 {
                *demand.entry(decoder.anticodon).or_insert(0.0) += count * weight;
            }
        }
    }
    demand
}

fn s_values() -> BTreeMap<String, f64> {
    // Wobble penalties, every pairing named explicitly. P2-12 refuses to
    // default one, so an incomplete table here is a refusal rather than a
    // silently-wrong weight — which is how this fixture was caught being
    // incomplete on its first run.
    [
        ("A:T", 0.0),
        ("T:A", 0.0),
        ("G:C", 0.0),
        ("C:G", 0.0),
        ("G:T", 0.41),
        ("T:G", 0.68),
        ("I:T", 0.0),
        ("I:C", 0.28),
        ("I:A", 0.9999),
    ]
    .into_iter()
    .map(|(key, value)| (key.to_string(), value))
    .collect()
}

fn counts() -> Vec<(&'static str, f64)> {
    vec![
        ("GCT", 120.0),
        ("GCC", 300.0),
        ("GCA", 80.0),
        ("GCG", 210.0),
        ("CTG", 400.0),
        ("CTT", 60.0),
    ]
}

/// `TST-BURD-005` — tier-3 semantics. The honesty requirement of ADR-0004.
#[test]
fn tst_burd_005_tier_3_reports_a_mismatch_score_not_a_growth_prediction() {
    let demand = native_demand_from_provenance(&counts(), &s_values());
    let estimate = fit_capacity_tier3(&demand).expect("tier 3");

    assert_eq!(estimate.tier(), 3);
    assert_eq!(estimate.unit(), ThetaUnit::FractionOfNativeDemand);
    assert_eq!(
        estimate.reported_unit(),
        "dimensionless",
        "DEP-BURD-003: reporting the same unit at all three tiers would be the single \
         most misleading thing DNA Design Suite could do"
    );
    assert!(
        !estimate.calibrated(),
        "tier 3 rests on no measurement in this host"
    );
    assert!(
        estimate.rationale_tags().contains(&"low_confidence_host"),
        "the decision layer's tag must be attached at the source, not by a caller who \
         might forget"
    );
    assert_eq!(estimate.theta(), 1.0, "theta = 1 by construction at tier 3");
    assert_eq!(estimate.g_max(), 1.0, "so growth_cost lands in [0,1)");

    let total: f64 = estimate.supply().values().sum();
    assert!(
        (total - 1.0).abs() < 1e-12,
        "supply shares must normalize: got {total}"
    );
}

/// Tier 1 is the only tier that may call itself calibrated.
#[test]
fn tst_burd_005_only_a_measured_tier_claims_calibration() {
    let demand = native_demand_from_provenance(&counts(), &s_values());

    let tier3 = fit_capacity_tier3(&demand).unwrap();
    let tier2 = fit_capacity_tier2(20_000.0, 16.0, &demand).unwrap();
    let panel: Vec<(f64, f64)> = (0..24)
        .map(|i| {
            let expression = 1.0 * 1.3_f64.powi(i);
            (expression, 0.6 * expression / (expression + 50.0))
        })
        .collect();
    let tier1 = fit_capacity_tier1(&panel, &demand).unwrap();

    assert!(!tier3.calibrated());
    assert!(
        !tier2.calibrated(),
        "tier 2 derives capacity from physiology, not from a measured burden curve; \
         calling it calibrated would claim a measurement nobody took"
    );
    assert!(tier1.calibrated());
    assert_eq!(tier1.unit(), ThetaUnit::AaPerCellPerSecond);
    assert_eq!(tier2.unit(), ThetaUnit::RibosomeOccupancyFraction);
}

/// `TST-BURD-006` — the estimate tracks the host's real data.
///
/// This is the check the audit asked for. Every assertion below would fail
/// against an estimator that returned a fixed profile, which is exactly what
/// `TST-BURD-005` alone cannot detect.
#[test]
fn tst_burd_006_the_estimate_moves_with_gene_copy_number() {
    let base = native_demand_from_provenance(&counts(), &s_values());
    let base_estimate = fit_capacity_tier3(&base).unwrap();

    // A host that uses GCC far more heavily. Its supply profile must shift
    // toward whatever decodes GCC.
    let mut shifted_counts = counts();
    shifted_counts[1].1 = 3_000.0;
    let shifted = native_demand_from_provenance(&shifted_counts, &s_values());
    let shifted_estimate = fit_capacity_tier3(&shifted).unwrap();

    assert_ne!(
        base_estimate.supply(),
        shifted_estimate.supply(),
        "changing the host's codon demand must change its supply profile; an estimator \
         that ignores gene_copy_number would pass every semantic check and be useless"
    );

    let moved = base_estimate
        .supply()
        .iter()
        .filter(|(anticodon, share)| {
            shifted_estimate
                .supply()
                .get(*anticodon)
                .is_some_and(|other| (other - *share).abs() > 1e-6)
        })
        .count();
    assert!(
        moved > 0,
        "at least one anticodon's share must move measurably"
    );
}

/// And with the tAI s-values, which is the other half of P2-12's provenance.
#[test]
fn tst_burd_006_the_estimate_moves_with_the_s_values() {
    let base = fit_capacity_tier3(&native_demand_from_provenance(&counts(), &s_values())).unwrap();

    // Forbid G:T wobble entirely. Codons previously decoded by a wobbling tRNA
    // must redistribute, so the supply profile must change.
    let mut strict = s_values();
    strict.insert("G:T".into(), 0.0);
    let strict_estimate =
        fit_capacity_tier3(&native_demand_from_provenance(&counts(), &strict)).unwrap();

    assert_ne!(
        base.supply(),
        strict_estimate.supply(),
        "the wobble weights feed the demand that feeds the supply; an estimator that \
         ignored tai.s_values would return the same profile for a host where G:T \
         wobble does not happen"
    );
}

/// No data means a refusal, not a default host.
#[test]
fn tst_burd_006_a_host_with_no_demand_data_is_refused() {
    let empty: BTreeMap<[u8; 3], f64> = BTreeMap::new();
    let error = fit_capacity_tier3(&empty).expect_err("no demand data must be refused");
    assert_eq!(error.code, "CF-BURD-004");
}

/// Tier 1's panel requirements are enforced, not advisory.
///
/// Fixture rule: both branches exercised. The passing case above uses 24
/// constructs over three decades; these are the two ways to fail.
#[test]
fn tst_burd_006_an_inadequate_panel_does_not_buy_tier_1() {
    let demand = native_demand_from_provenance(&counts(), &s_values());

    let too_few: Vec<(f64, f64)> = (0..12).map(|i| (1.0 * 2.0_f64.powi(i), 0.3)).collect();
    let error = fit_capacity_tier1(&too_few, &demand).expect_err("12 constructs is not 20");
    assert_eq!(error.code, "CF-BURD-004");

    // Enough constructs, but all crowded into half a decade: the saturating
    // region is never observed, so theta is unidentifiable.
    let narrow: Vec<(f64, f64)> = (0..24).map(|i| (100.0 + i as f64, 0.3)).collect();
    let error = fit_capacity_tier1(&narrow, &demand)
        .expect_err("a panel spanning no dynamic range cannot identify theta");
    assert_eq!(error.code, "CF-BURD-004");
}

/// Tier selection reads the host's data, and never claims more than it has.
#[test]
fn tst_burd_006_tier_selection_never_exceeds_the_data() {
    let nothing = DataCompleteness::default();
    assert_eq!(
        nothing.tier(),
        3,
        "the no-data case is where most hosts land"
    );

    let ribosomes = DataCompleteness {
        has_ribosome_count: true,
        ..Default::default()
    };
    assert_eq!(ribosomes.tier(), 2);

    let panel = DataCompleteness {
        has_burden_panel: true,
        has_ribosome_count: true,
        ..Default::default()
    };
    assert_eq!(panel.tier(), 1);

    // Having ribosome counts does not buy tier 1 without a burden panel.
    assert_eq!(
        DataCompleteness {
            has_ribosome_count: true,
            has_trna_gene_copies: true,
            has_codon_usage: true,
            has_burden_panel: false,
        }
        .tier(),
        2,
        "the tier is the highest the data supports, never higher"
    );
}
