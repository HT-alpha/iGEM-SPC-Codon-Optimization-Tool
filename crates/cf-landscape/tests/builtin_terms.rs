//! `P1-21` — the three built-in terms, and the delta test applied to each.
//!
//! The acceptance criterion is a `TST-EPI-002`-style delta test: for every
//! synonymous alternative at every editable position, the claimed exact delta
//! equals the change a full rescore observes. That is the property the whole
//! capability model rests on, so it is checked exhaustively over the fixture
//! rather than sampled.

mod support;

use cf_constraints::ConstraintChecker;
use cf_landscape::terms::usage::UsageSubset;
use cf_landscape::{CaiTerm, CodonUsageTable, GcContentTerm, ViolationCountTerm};
use cf_landscape_api::{verify_capability, Capability, Edit, TermProvider, TermRole};
use cf_types::{Codon, CodonSequence, GeneticCode};

use support::{
    ctx, flanked_sequence, global_constraints, host, layout, local_constraints, sequence,
    usage_table, usage_table_json,
};

/// Every claimed delta equals the change a rescore observes.
///
/// Returns how many alternatives were checked so a caller can assert the
/// fixture actually exercised the path — a term whose sequence had no
/// synonymous choice would pass this loop without comparing anything.
fn assert_deltas_are_exact(provider: &dyn TermProvider, base: &CodonSequence) -> usize {
    let host = host();
    let base_ctx = ctx(base, &host);
    let pinned = provider.rebase(&base_ctx).unwrap();
    let baseline = provider.score(&base_ctx).unwrap();
    assert_eq!(
        pinned.value().to_bits(),
        baseline.to_bits(),
        "{}: the pinned base disagrees with a fresh score, so every delta against it is offset",
        provider.term_id()
    );
    let mut checked = 0;

    for index in base.editable_positions() {
        let current = base.codons()[index];
        for alternative in base.synonyms_at(index).iter().copied() {
            if alternative == current {
                continue;
            }
            let mutated = base.substituted(index, alternative).unwrap();
            let truth = provider.score(&ctx(&mutated, &host)).unwrap() - baseline;
            let claimed = provider
                .delta(&base_ctx, &pinned, &[Edit::new(index, alternative)])
                .unwrap();
            assert!(
                (claimed - truth).abs() <= 1e-12 + 1e-9 * truth.abs().max(claimed.abs()),
                "{}: delta at codon {index} -> {} was {claimed}, rescoring says {truth}",
                provider.term_id(),
                alternative.to_string_upper()
            );
            checked += 1;
        }
    }
    checked
}

/// The delta table agrees with the single-edit delta at every entry.
fn assert_table_matches_delta(provider: &dyn TermProvider, base: &CodonSequence) {
    let host = host();
    let base_ctx = ctx(base, &host);
    let pinned = provider.rebase(&base_ctx).unwrap();
    let positions: Vec<usize> = base.editable_positions().collect();
    let table = provider
        .delta_table(&base_ctx, &pinned, &positions)
        .unwrap();

    assert_eq!(table.len(), positions.len());
    for entry in table {
        for (codon, delta) in entry.alternatives {
            let single = provider
                .delta(&base_ctx, &pinned, &[Edit::new(entry.codon_index, codon)])
                .unwrap();
            assert_eq!(
                delta.to_bits(),
                single.to_bits(),
                "{}: table and single-edit delta disagree at codon {}",
                provider.term_id(),
                entry.codon_index
            );
        }
    }
}

// --------------------------------------------------------------------------- //
// usage.cai
// --------------------------------------------------------------------------- //

#[test]
fn p1_21_cai_deltas_are_exact() {
    let term = CaiTerm::new(usage_table("highly_expressed")).unwrap();
    let base = sequence();
    let checked = assert_deltas_are_exact(&term, &base);
    assert!(
        checked >= 20,
        "fixture exercised only {checked} alternatives; the delta test would be near-vacuous"
    );
    assert_table_matches_delta(&term, &base);
    verify_capability(&term, &ctx(&base, &host())).unwrap();
}

/// CAI is the geometric mean over *counted* positions: no stops, no
/// single-codon families.
///
/// Checked against an independently computed value rather than against the
/// implementation's own state, so a change to which positions are counted
/// fails here instead of silently redefining the term.
#[test]
fn p1_21_cai_excludes_stops_and_single_codon_families() {
    let table = usage_table("highly_expressed");
    let term = CaiTerm::new(table.clone()).unwrap();
    let seq = sequence();
    let code = GeneticCode::get(1).unwrap();

    let mut logs = Vec::new();
    for codon in seq.codons().iter().copied() {
        if code.is_stop(codon) || code.synonyms(codon).len() == 1 {
            continue;
        }
        logs.push(table.weight(codon).ln());
    }
    let expected = (logs.iter().sum::<f64>() / logs.len() as f64).exp();

    let value = term.score(&ctx(&seq, &host())).unwrap();
    assert!((value - expected).abs() < 1e-12, "{value} != {expected}");

    // The fixture must actually contain an excluded position, or the assertion
    // above would hold for the wrong reason.
    assert!(
        seq.codons()
            .iter()
            .any(|c| code.is_stop(*c) || code.synonyms(*c).len() == 1),
        "fixture has no stop and no single-codon family, so it cannot show the exclusion"
    );
    assert!(logs.len() < seq.codon_count());
}

/// A genome-wide table is refused rather than used.
///
/// The host profile schema calls this "a classic silent error", and silent is
/// the operative word: the arithmetic succeeds and the answer looks reasonable.
#[test]
fn p1_21_cai_refuses_a_genome_wide_reference() {
    let err = CaiTerm::new(usage_table("genome_wide")).unwrap_err();
    assert_eq!(err.code, "CF-HOST-003");
    assert!(err.detail.contains("genome_wide"));

    assert!(CaiTerm::new(usage_table("highly_expressed")).is_ok());
    assert!(CaiTerm::new(usage_table("user_supplied")).is_ok());
}

/// An amino acid with no observed codon at all makes CAI uninformative, and
/// that is reported rather than returned as zero.
#[test]
fn p1_21_cai_refuses_an_entirely_unobserved_family() {
    let mut document = usage_table_json("highly_expressed");
    let code = GeneticCode::get(1).unwrap();
    let tryptophan_free_family = code.synonyms(Codon::parse("TGT").unwrap()).to_vec();
    for codon in &tryptophan_free_family {
        document["codon_usage"]["table"][codon.to_string_upper()]["relative_adaptiveness"] =
            serde_json::json!(0.0);
    }

    let table = CodonUsageTable::from_host_profile(&document, 1).unwrap();
    let err = CaiTerm::new(table).unwrap_err();
    assert_eq!(err.code, "CF-HOST-003");
    assert!(err.detail.contains('C'), "the message names the amino acid");
}

/// A single zero weight is legitimate: the codon scores badly and the search
/// moves off it. That is CAI doing its job, not a failure.
#[test]
fn p1_21_cai_tolerates_an_individual_zero_weight() {
    let mut document = usage_table_json("highly_expressed");
    document["codon_usage"]["table"]["CTA"]["relative_adaptiveness"] = serde_json::json!(0.0);
    let table = CodonUsageTable::from_host_profile(&document, 1).unwrap();
    let term = CaiTerm::new(table).unwrap();

    let base = sequence();
    let with_zero = base.substituted(1, Codon::parse("CTA").unwrap()).unwrap();
    assert_eq!(term.score(&ctx(&with_zero, &host())).unwrap(), 0.0);
    assert!(term.score(&ctx(&base, &host())).unwrap() > 0.0);

    // And the delta remains exact across the boundary in both directions.
    assert_deltas_are_exact(&term, &with_zero);
}

/// The pinned state is load-bearing, not decorative.
///
/// Correctness alone cannot show this: an implementation that ignores the pin
/// and recomputes returns the same numbers and differs only in costing `O(L)`
/// per alternative instead of `O(1)` — the `DEP-LAND-003` defect. So the pin is
/// poisoned instead. A `Rebased` built from a *different* sequence is invalid
/// input, and composition is the only thing that constructs one, but it makes
/// the question observable: an implementation that consults the pin produces
/// deltas consistent with the wrong base, and one that recomputes produces
/// correct deltas and thereby convicts itself.
#[test]
fn p1_21_cai_uses_its_pinned_state_rather_than_recomputing() {
    let term = CaiTerm::new(usage_table("highly_expressed")).unwrap();
    let host = host();

    let here = sequence();
    let elsewhere = here
        .substituted(4, Codon::parse("GGT").unwrap())
        .unwrap()
        .substituted(5, Codon::parse("GTA").unwrap())
        .unwrap()
        .substituted(6, Codon::parse("GCT").unwrap())
        .unwrap();

    let pin_here = term.rebase(&ctx(&here, &host)).unwrap();
    let pin_elsewhere = term.rebase(&ctx(&elsewhere, &host)).unwrap();
    assert_ne!(
        pin_here.value().to_bits(),
        pin_elsewhere.value().to_bits(),
        "the two pins must differ or the assertion below proves nothing"
    );

    let honest = term
        .delta_table(&ctx(&here, &host), &pin_here, &[1, 2])
        .unwrap();
    let poisoned = term
        .delta_table(&ctx(&here, &host), &pin_elsewhere, &[1, 2])
        .unwrap();
    assert_ne!(
        honest, poisoned,
        "delta_table ignored its pinned state and recomputed the base"
    );

    let honest_single = term
        .delta(
            &ctx(&here, &host),
            &pin_here,
            &[Edit::new(1, Codon::parse("CTC").unwrap())],
        )
        .unwrap();
    let poisoned_single = term
        .delta(
            &ctx(&here, &host),
            &pin_elsewhere,
            &[Edit::new(1, Codon::parse("CTC").unwrap())],
        )
        .unwrap();
    assert_ne!(
        honest_single.to_bits(),
        poisoned_single.to_bits(),
        "delta ignored its pinned state and recomputed the base"
    );

    // And the honest path is still the correct one.
    let baseline = term.score(&ctx(&here, &host)).unwrap();
    let mutated = here.substituted(1, Codon::parse("CTC").unwrap()).unwrap();
    let truth = term.score(&ctx(&mutated, &host)).unwrap() - baseline;
    assert!((honest_single - truth).abs() < 1e-12);
}

#[test]
fn p1_21_cai_subset_is_carried() {
    assert_eq!(
        usage_table("highly_expressed").subset(),
        UsageSubset::HighlyExpressed
    );
}

// --------------------------------------------------------------------------- //
// usage.gc_content
// --------------------------------------------------------------------------- //

#[test]
fn p1_21_gc_deltas_are_exact() {
    let term = GcContentTerm::new();
    let base = sequence();
    assert!(assert_deltas_are_exact(&term, &base) >= 20);
    assert_table_matches_delta(&term, &base);
    verify_capability(&term, &ctx(&base, &host())).unwrap();
}

#[test]
fn p1_21_gc_value_is_the_cds_fraction() {
    let seq = sequence();
    let expected = support::CDS
        .bytes()
        .filter(|b| *b == b'G' || *b == b'C')
        .count() as f64
        / support::CDS.len() as f64;
    assert!((GcContentTerm::new().score(&ctx(&seq, &host())).unwrap() - expected).abs() < 1e-12);
}

/// GC content cannot be a Pareto axis.
///
/// The quantity has no direction — too high and the transcript folds, too low
/// and it does so differently — so ranking on it would mean asserting one the
/// biology does not have.
#[test]
fn p1_21_gc_is_capped_at_diagnostic() {
    assert_eq!(GcContentTerm::new().max_role(), TermRole::Diagnostic);
}

// --------------------------------------------------------------------------- //
// constraint.violation_count
// --------------------------------------------------------------------------- //

#[test]
fn p1_21_violation_count_is_l1_for_a_local_set_and_deltas_are_exact() {
    let checker = ConstraintChecker::compile(&local_constraints()).unwrap();
    let term = ViolationCountTerm::new(checker, layout()).unwrap();
    assert_eq!(term.capability(), Capability::L1);

    let base = flanked_sequence();
    let checked = assert_deltas_are_exact(&term, &base);
    assert!(checked >= 20);
    assert_table_matches_delta(&term, &base);
    verify_capability(&term, &ctx(&base, &host())).unwrap();
}

/// The fixture must contain a real violation, or the delta test compares zero
/// with zero at every position.
#[test]
fn p1_21_violation_count_fixture_actually_violates_something() {
    let checker = ConstraintChecker::compile(&local_constraints()).unwrap();
    let term = ViolationCountTerm::new(checker, layout()).unwrap();
    let base = flanked_sequence();
    assert!(
        term.score(&ctx(&base, &host())).unwrap() > 0.0,
        "the fixture construct satisfies every constraint, so nothing is being counted"
    );
}

/// A k-mer repeat rule counts over the whole construct, so no window bounds a
/// substitution's effect and the term is honestly L0.
///
/// Declaring L1 and computing the "delta" by rescoring twice would be exact and
/// would cost `O(L)` — a cheap path advertised while the expensive one is
/// walked, which defeats the argument the capability model is made of.
#[test]
fn p1_21_violation_count_degrades_to_l0_for_a_global_set() {
    let checker = ConstraintChecker::compile(&global_constraints()).unwrap();
    let term = ViolationCountTerm::new(checker, layout()).unwrap();
    assert_eq!(term.capability(), Capability::L0);

    let base = flanked_sequence();
    let host = host();
    let base_ctx = ctx(&base, &host);
    let pinned = term.rebase(&base_ctx).unwrap();
    let err = term
        .delta(
            &base_ctx,
            &pinned,
            &[Edit::new(1, Codon::parse("CTC").unwrap())],
        )
        .unwrap_err();
    assert_eq!(err.code, "CF-LAND-002");
    assert!(err.detail.contains("k-mer repeat"));

    // The declared capability still matches the surface: L0 promises a score
    // and nothing else, and that is exactly what it delivers.
    verify_capability(&term, &ctx(&base, &host)).unwrap();
}
