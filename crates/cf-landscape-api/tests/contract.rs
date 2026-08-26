//! `TST-LAPI-001`, `TST-LAPI-002`, `TST-LAPI-003`.
//!
//! Ticket P1-20.

mod support;

use cf_landscape_api::{
    verify_capability, Budget, BudgetTier, CachePolicy, Capability, ConstructRef, DeltaEntry,
    Deltas, Diagnostics, Edit, EvalContext, HostContext, LandscapeQuery, LandscapeResponse,
    ModelPin, ObjectiveAxis, ObjectiveVector, Outputs, Provenance, Status, Target, TermProvider,
    TermRegistry, TermRole, TermValue, Uncertainty,
};
use cf_types::{ArtifactId, Codon, CodonSequence, Digest, Direction, HostId, TermId, Unit};

use support::{
    AdditiveTerm, FailingTerm, LyingL1Term, MissingTableTerm, OrderSensitiveTerm, SampledTerm,
    ScoreOnlyTerm,
};

// --------------------------------------------------------------------------- //
// fixtures
// --------------------------------------------------------------------------- //

/// A CDS with several redundant amino acids, so there are real synonymous moves.
///
/// `M L R S G V *` — leucine, arginine and serine each have six codons, so the
/// delta table is non-trivial and a provider cannot pass by returning an empty
/// alternative list.
const CDS: &str = "ATGCTGCGTAGCGGCGTGTAA";

fn sequence() -> CodonSequence {
    CodonSequence::from_cds(CDS, 1).unwrap()
}

fn host() -> HostContext {
    HostContext::new(
        HostId::new("escherichia_coli:511145/BL21").unwrap(),
        Digest::of_bytes(b"host-profile-fixture"),
    )
}

fn terms(names: &[&str]) -> Vec<TermId> {
    names.iter().map(|n| TermId::new(*n).unwrap()).collect()
}

fn term(name: &str) -> TermId {
    TermId::new(name).unwrap()
}

fn ctx<'a>(sequence: &'a CodonSequence, host: &'a HostContext, seed: u64) -> EvalContext<'a> {
    EvalContext::new(sequence, host, Budget::new(BudgetTier::Standard), seed)
}

// --------------------------------------------------------------------------- //
// TST-LAPI-001 — document round-trip
// --------------------------------------------------------------------------- //

/// Every golden query document survives Rust -> JSON -> Rust -> JSON unchanged.
///
/// The Python leg of `TST-LAPI-001` lives in
/// `python/tests/test_contract_landscape_schema.py`, which validates the same
/// documents against the schema and re-canonicalizes them.
#[test]
fn tst_lapi_001_query_documents_round_trip() {
    let base = Digest::of_bytes(b"base-sequence-document");

    let full = LandscapeQuery::new(
        "q_01ARZ3NDEKTSV4RRFFQ69G5FAV",
        Target::Full {
            sequence: Box::new(sequence()),
            cistron_id: Some("cds_reporter".into()),
        },
        host().with_condition("cond_m9_glucose_37c"),
        terms(&[
            "usage.cai",
            "usage.gc_content",
            "constraint.violation_count",
        ]),
        Budget::new(BudgetTier::Standard).with_latency_ms(250),
    )
    .with_outputs(Outputs::default().with_uncertainty())
    .with_cache_policy(CachePolicy::namespaced(
        "sha256:".to_owned() + &"a".repeat(64),
    ))
    .with_seed(20260826);

    let delta = LandscapeQuery::new(
        "q_01ARZ3NDEKTSV4RRFFQ69G5FAW",
        Target::Delta {
            base_digest: base.clone(),
            edits: vec![
                Edit::new(1, Codon::parse("CTC").unwrap()),
                Edit::new(3, Codon::parse("TCT").unwrap()),
            ],
            cistron_id: None,
        },
        host(),
        terms(&["usage.cai"]),
        Budget::new(BudgetTier::Fast),
    );

    let mut table = LandscapeQuery::new(
        "q_01ARZ3NDEKTSV4RRFFQ69G5FAX",
        Target::DeltaTable {
            base_digest: base.clone(),
            positions: vec![1, 2, 3, 4],
            cistron_id: None,
        },
        host(),
        terms(&["epistasis.dwell_deviation", "usage.gc_content"]),
        Budget::new(BudgetTier::Exact),
    )
    .with_outputs(Outputs::default().with_deltas());
    // Pinned to a term that actually has an artifact. CAI is `builtin` — it
    // has no artifact to pin, and a fixture that pinned one would be modelling
    // something the resolution ladder cannot produce.
    table.model_pins = vec![(
        term("epistasis.dwell_deviation"),
        ModelPin {
            artifact_id: ArtifactId::new("cf:epistasis_dwell:escherichia_coli:511145/BL21:1.2.0")
                .unwrap(),
            digest: Digest::of_bytes(b"dwell-artifact"),
        },
    )];
    table.construct_context = Some(ConstructRef::Digest(Digest::of_bytes(b"construct")));

    for query in [full, delta, table] {
        let once = query.to_canonical_json().unwrap();
        let reparsed = LandscapeQuery::from_json(&query.to_json().unwrap()).unwrap();
        let twice = reparsed.to_canonical_json().unwrap();
        assert_eq!(
            once, twice,
            "query {} did not survive a round trip",
            query.query_id
        );
    }
}

#[test]
fn tst_lapi_001_response_documents_round_trip() {
    let axes = vec![
        ObjectiveAxis::new(term("usage.cai"), TermRole::ParetoAxis).unwrap(),
        ObjectiveAxis::new(term("usage.gc_content"), TermRole::ParetoAxis).unwrap(),
    ];
    let response = LandscapeResponse {
        query_id: "q_01ARZ3NDEKTSV4RRFFQ69G5FAV".into(),
        status: Status::Partial,
        terms: vec![
            (
                term("usage.cai"),
                TermValue::new(
                    0.7421,
                    Unit::new(Unit::DIMENSIONLESS),
                    Direction::Maximize,
                    Capability::L1,
                    Provenance::builtin("cai_native", "0.1.0"),
                )
                .with_uncertainty(Uncertainty::deterministic()),
            ),
            (
                term("usage.gc_content"),
                TermValue::new(
                    0.5238,
                    Unit::new(Unit::DIMENSIONLESS),
                    Direction::Minimize,
                    Capability::L1,
                    Provenance::builtin("gc_native", "0.1.0"),
                ),
            ),
        ],
        objective_vector: Some(
            ObjectiveVector::new(
                axes,
                vec![Direction::Maximize, Direction::Minimize],
                vec![Some(0.7421), None],
            )
            .unwrap(),
        ),
        deltas: Some(Deltas {
            base_digest: Digest::of_bytes(b"base-sequence-document"),
            exact: true,
            entries: vec![DeltaEntry {
                codon_index: 1,
                alternatives: vec![cf_landscape_api::AlternativeDelta {
                    codon: Some(Codon::parse("CTC").unwrap()),
                    delta: vec![(term("usage.cai"), -0.0125)],
                }],
            }],
        }),
        diagnostics: Diagnostics {
            cache_hits: 12,
            cache_misses: 3,
            evaluations: 3,
            full_rescores: 1,
            rebases: 1,
            terms_deferred: terms(&["usage.gc_content"]),
            engine_notes: vec!["deferred by latency budget".into()],
        },
        problem: None,
    };

    let once = response.to_canonical_json().unwrap();
    let twice = LandscapeResponse::from_json(&response.to_json().unwrap())
        .unwrap()
        .to_canonical_json()
        .unwrap();
    assert_eq!(once, twice);
}

/// A response carrying a gradient block is refused, not silently dropped.
///
/// Dropping would let an L2 path exist and be invisible, which is the erosion
/// ADR-0002 guards against.
#[test]
fn tst_lapi_001_gradient_block_is_refused_on_parse() {
    let mut document = LandscapeResponse {
        query_id: "q_01ARZ3NDEKTSV4RRFFQ69G5FAV".into(),
        status: Status::Ok,
        terms: Vec::new(),
        objective_vector: None,
        deltas: None,
        diagnostics: Diagnostics::default(),
        problem: None,
    }
    .to_json()
    .unwrap();
    document["gradients"] = serde_json::json!([{
        "term_id": "lm.naturalness",
        "parameterization": "relaxed_codon_logits",
        "shape": [7, 6],
        "values": [[0.0]],
    }]);

    let err = LandscapeResponse::from_json(&document).unwrap_err();
    assert_eq!(err.code, "CF-LAND-002");
    assert!(err.detail.contains("ADR-0002"));
}

// --------------------------------------------------------------------------- //
// TST-LAPI-002 — a capability claim is checked by calling
// --------------------------------------------------------------------------- //

#[test]
fn tst_lapi_002_honest_providers_verify() {
    let seq = sequence();
    let host = host();
    let ctx = ctx(&seq, &host, 0);

    verify_capability(&AdditiveTerm::new("usage.cai"), &ctx).expect("honest L1 verifies");
    verify_capability(&ScoreOnlyTerm::new("usage.gc_content"), &ctx).expect("honest L0 verifies");
}

/// An L1 claim without a *correct* delta fails.
///
/// Both failure shapes are covered, because they are different mistakes: a
/// provider that implements the method wrongly, and one that never implemented
/// it. A check that only caught the second would pass the first, and the first
/// is the one that corrupts a search silently.
#[test]
fn tst_lapi_002_false_l1_claims_fail() {
    let seq = sequence();
    let host = host();
    let ctx = ctx(&seq, &host, 0);

    let lying =
        verify_capability(&LyingL1Term::new("epistasis.dwell_deviation"), &ctx).unwrap_err();
    assert_eq!(lying.code, "CF-LAND-003");
    assert!(
        lying.detail.contains("rescoring gives"),
        "the failure must name the disagreement, got: {}",
        lying.detail
    );

    let missing = verify_capability(&MissingTableTerm::new("lm.naturalness"), &ctx).unwrap_err();
    assert_eq!(missing.code, "CF-LAND-002");
}

/// Verification refuses to pass vacuously on a sequence with no moves.
///
/// `ATGTGGTAA` is Met-Trp-Stop: methionine and tryptophan have one codon each,
/// so there is no synonymous substitution to check a delta against. Reporting
/// that is the difference between "verified" and "found nothing to verify".
#[test]
fn tst_lapi_002_unverifiable_probe_is_reported_not_passed() {
    let seq = CodonSequence::from_cds("ATGTGGTAA", 1).unwrap();
    let host = host();
    let ctx = ctx(&seq, &host, 0);

    let err = verify_capability(&AdditiveTerm::new("usage.cai"), &ctx).unwrap_err();
    assert_eq!(err.code, "CF-LAND-002");
    assert!(err.detail.contains("no editable position"));
}

// --------------------------------------------------------------------------- //
// TST-LAPI-003 — the three invariants composition depends on
// --------------------------------------------------------------------------- //

/// (a) Resolution order is the caller's order, not the registry's.
///
/// Registered in one order and requested in three others. A registry that
/// returned its own iteration order would pass with one of them and fail the
/// rest — which is exactly how this defect survives casual testing.
#[test]
fn tst_lapi_003_resolution_order_is_the_callers_order() {
    let mut registry = TermRegistry::new();
    for name in [
        "usage.gc_content",
        "constraint.violation_count",
        "usage.cai",
    ] {
        registry
            .register_default(Box::new(AdditiveTerm::new(name)))
            .unwrap();
    }

    for requested in [
        vec![
            "usage.cai",
            "usage.gc_content",
            "constraint.violation_count",
        ],
        vec![
            "constraint.violation_count",
            "usage.cai",
            "usage.gc_content",
        ],
        vec!["usage.gc_content", "usage.cai"],
    ] {
        let wanted = terms(&requested);
        let resolved: Vec<TermId> = registry
            .resolve(&wanted)
            .unwrap()
            .iter()
            .map(|r| r.provider.term_id().clone())
            .collect();
        assert_eq!(
            resolved, wanted,
            "resolve() must preserve the caller's order"
        );
    }

    let missing = registry.resolve(&terms(&["usage.cai", "structure.local_accessibility"]));
    assert_eq!(missing.unwrap_err().code, "CF-LAND-001");
}

/// (b) The same seed gives the same answer; a different seed does not.
///
/// The second half matters as much as the first: a provider that ignored the
/// seed entirely would satisfy determinism trivially while making the `seed`
/// field a decoration.
#[test]
fn tst_lapi_003_evaluation_is_deterministic_under_a_seed() {
    let seq = sequence();
    let host = host();
    let sampled = SampledTerm::new("lm.naturalness");

    let a = sampled.score(&ctx(&seq, &host, 4242)).unwrap();
    let b = sampled.score(&ctx(&seq, &host, 4242)).unwrap();
    assert_eq!(
        a.to_bits(),
        b.to_bits(),
        "same seed must give the same bits"
    );

    let c = sampled.score(&ctx(&seq, &host, 4243)).unwrap();
    assert_ne!(
        a.to_bits(),
        c.to_bits(),
        "a provider that ignores the seed makes the field a decoration"
    );

    // A deterministic provider is seed-independent, which is why the query
    // only demands a seed when the provider is sampled.
    let additive = AdditiveTerm::new("usage.cai");
    assert_eq!(
        additive.score(&ctx(&seq, &host, 1)).unwrap().to_bits(),
        additive.score(&ctx(&seq, &host, 2)).unwrap().to_bits()
    );
}

/// (c) Evaluating a term set in any permutation gives identical values.
#[test]
fn tst_lapi_003_terms_are_independent_of_evaluation_order() {
    let seq = sequence();
    let host = host();
    let ctx = ctx(&seq, &host, 11);

    let providers: Vec<Box<dyn TermProvider>> = vec![
        Box::new(AdditiveTerm::new("usage.cai")),
        Box::new(AdditiveTerm::new("usage.gc_content")),
        Box::new(ScoreOnlyTerm::new("constraint.violation_count")),
    ];

    let mut reference: Vec<(String, u64)> = Vec::new();
    for p in &providers {
        reference.push((
            p.term_id().as_str().to_owned(),
            p.score(&ctx).unwrap().to_bits(),
        ));
    }
    reference.sort();

    for permutation in [[2usize, 0, 1], [1, 2, 0], [2, 1, 0]] {
        let mut observed: Vec<(String, u64)> = Vec::new();
        for i in permutation {
            let p = &providers[i];
            observed.push((
                p.term_id().as_str().to_owned(),
                p.score(&ctx).unwrap().to_bits(),
            ));
        }
        observed.sort();
        assert_eq!(
            observed, reference,
            "evaluation order changed a term's value"
        );
    }
}

/// The fixture that proves the independence assertion can fail.
///
/// Without this, `tst_lapi_003_terms_are_independent_of_evaluation_order`
/// would pass on any set of providers and would be evidence of nothing.
#[test]
fn tst_lapi_003_independence_assertion_catches_a_stateful_provider() {
    let seq = sequence();
    let host = host();
    let ctx = ctx(&seq, &host, 0);

    let stateful = OrderSensitiveTerm::new("meta.call_count");
    let first = stateful.score(&ctx).unwrap();
    let second = stateful.score(&ctx).unwrap();
    assert_ne!(
        first, second,
        "the negative-control fixture must actually be order-sensitive"
    );
}

// --------------------------------------------------------------------------- //
// role gate
// --------------------------------------------------------------------------- //

/// A provider capped at diagnostic cannot be registered as an axis or a bound.
///
/// This is `AMD-002` enforced at the point of registration rather than at the
/// point of ranking. By the time a diagnostic term reached a ranking it would
/// already have changed the answer.
#[test]
fn tst_lapi_003_a_diagnostic_capped_provider_cannot_be_ranked_on() {
    let mut registry = TermRegistry::new();

    let as_axis = registry.register(
        Box::new(AdditiveTerm::new("kinetics.pause_concordance").diagnostic_only()),
        TermRole::ParetoAxis,
    );
    assert_eq!(as_axis.unwrap_err().code, "CF-LAND-004");

    let as_bound = registry.register(
        Box::new(AdditiveTerm::new("kinetics.pause_concordance").diagnostic_only()),
        TermRole::EpsilonConstraint,
    );
    assert_eq!(
        as_bound.unwrap_err().code,
        "CF-LAND-004",
        "an epsilon constraint changes feasibility, so it changes the ranking too"
    );

    registry
        .register(
            Box::new(AdditiveTerm::new("kinetics.pause_concordance").diagnostic_only()),
            TermRole::Diagnostic,
        )
        .expect("diagnostic is the role it is permitted to hold");
    assert_eq!(
        registry.role_of(&term("kinetics.pause_concordance")),
        Some(TermRole::Diagnostic)
    );
}

#[test]
fn tst_lapi_003_duplicate_registration_is_refused() {
    let mut registry = TermRegistry::new();
    registry
        .register_default(Box::new(AdditiveTerm::new("usage.cai")))
        .unwrap();
    let again = registry.register_default(Box::new(ScoreOnlyTerm::new("usage.cai")));
    assert_eq!(again.unwrap_err().code, "CF-LAND-005");
}

#[test]
fn tst_lapi_003_a_failing_provider_surfaces_its_code() {
    let seq = sequence();
    let host = host();
    let err = FailingTerm::new("structure.local_accessibility")
        .score(&ctx(&seq, &host, 0))
        .unwrap_err();
    assert_eq!(err.code, "CF-BUDGET-003");
}
