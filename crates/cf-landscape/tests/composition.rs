//! `TST-LAND-001` … `TST-LAND-005` — composition, cache, and degradation.
//!
//! Ticket P1-22.

mod support;

use std::collections::{BTreeSet, HashMap};

use cf_constraints::ConstraintChecker;
use cf_landscape::compose::HostCapability;
use cf_landscape::{CaiTerm, Composer, EvaluationCache, GcContentTerm, ViolationCountTerm};
use cf_landscape_api::{
    Budget, BudgetTier, CachePolicy, Capability, Edit, LandscapeQuery, Outputs, Status, Target,
    TermProvider, TermRegistry, TermRole,
};
use cf_types::{Codon, Digest, TermId};

use support::{
    ctx, flanked_sequence, global_constraints, host, layout, local_constraints, usage_table,
};

const LOCKFILE_A: &str = "sha256:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
const LOCKFILE_B: &str = "sha256:bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb";

fn terms(names: &[&str]) -> Vec<TermId> {
    names.iter().map(|n| TermId::new(*n).unwrap()).collect()
}

fn term(name: &str) -> TermId {
    TermId::new(name).unwrap()
}

/// A registry whose violation counter is L1 (local constraint set).
fn local_registry() -> TermRegistry {
    let mut registry = TermRegistry::new();
    registry
        .register_default(Box::new(
            CaiTerm::new(usage_table("highly_expressed")).unwrap(),
        ))
        .unwrap();
    registry
        .register_default(Box::new(GcContentTerm::new()))
        .unwrap();
    registry
        .register_default(Box::new(
            ViolationCountTerm::new(
                ConstraintChecker::compile(&local_constraints()).unwrap(),
                layout(),
            )
            .unwrap(),
        ))
        .unwrap();
    registry
}

/// The same registry, except the violation counter is L0 because the
/// constraint set contains a k-mer repeat rule.
fn mixed_registry() -> TermRegistry {
    let mut registry = TermRegistry::new();
    registry
        .register_default(Box::new(
            CaiTerm::new(usage_table("highly_expressed")).unwrap(),
        ))
        .unwrap();
    registry
        .register_default(Box::new(GcContentTerm::new()))
        .unwrap();
    registry
        .register_default(Box::new(
            ViolationCountTerm::new(
                ConstraintChecker::compile(&global_constraints()).unwrap(),
                layout(),
            )
            .unwrap(),
        ))
        .unwrap();
    registry
}

fn full_query(id: &str, requested: &[&str]) -> LandscapeQuery {
    LandscapeQuery::new(
        id,
        Target::Full {
            sequence: Box::new(flanked_sequence()),
            cistron_id: None,
        },
        host(),
        terms(requested),
        Budget::new(BudgetTier::Standard),
    )
}

fn delta_query(id: &str, base: &Digest, requested: &[&str], edits: Vec<Edit>) -> LandscapeQuery {
    LandscapeQuery::new(
        id,
        Target::Delta {
            base_digest: base.clone(),
            edits,
            cistron_id: None,
        },
        host(),
        terms(requested),
        Budget::new(BudgetTier::Standard),
    )
}

const ALL: &[&str] = &[
    "usage.cai",
    "usage.gc_content",
    "constraint.violation_count",
];

// --------------------------------------------------------------------------- //
// TST-LAND-001 — cache correctness and namespace isolation
// --------------------------------------------------------------------------- //

#[test]
fn tst_land_001_a_cache_hit_equals_a_cold_evaluation() {
    let mut warm = Composer::new(local_registry(), LOCKFILE_A);
    let cold = &mut Composer::new(local_registry(), LOCKFILE_A);

    let first = warm
        .evaluate(&full_query("q_01ARZ3NDEKTSV4RRFFQ69G5FAV", ALL))
        .unwrap();
    let second = warm
        .evaluate(&full_query("q_01ARZ3NDEKTSV4RRFFQ69G5FAV", ALL))
        .unwrap();
    let fresh = cold
        .evaluate(&full_query("q_01ARZ3NDEKTSV4RRFFQ69G5FAV", ALL))
        .unwrap();

    assert_eq!(first.diagnostics.cache_hits, 0);
    assert_eq!(first.diagnostics.cache_misses, 3);
    assert_eq!(
        second.diagnostics.cache_hits, 3,
        "the second pass is all hits"
    );
    assert_eq!(second.diagnostics.evaluations, 0);

    for name in ALL {
        let t = term(name);
        let a = first.get(&t).unwrap().value.to_bits();
        let b = second.get(&t).unwrap().value.to_bits();
        let c = fresh.get(&t).unwrap().value.to_bits();
        assert_eq!(a, b, "{name}: cache hit differs from the cold value");
        assert_eq!(a, c, "{name}: a fresh composer differs");
    }
    assert_eq!(
        first.to_canonical_json().unwrap(),
        fresh.to_canonical_json().unwrap()
    );
}

/// Two lockfiles cannot see each other's entries.
///
/// The check is not that the values differ — with these built-ins they do not,
/// because the terms have no artifacts. It is that the second composer records
/// a *miss*: a shared key space would have handed it the first one's entry
/// regardless of what the numbers happened to be.
#[test]
fn tst_land_001_lockfile_digest_partitions_the_cache() {
    let mut composer = Composer::new(local_registry(), LOCKFILE_A);
    let query = full_query("q_01ARZ3NDEKTSV4RRFFQ69G5FAV", ALL);
    composer.evaluate(&query).unwrap();

    let warmed = composer.cache().len();
    assert_eq!(warmed, 3);

    let mut other =
        Composer::new(local_registry(), LOCKFILE_B).with_cache(EvaluationCache::new(LOCKFILE_B));
    let response = other.evaluate(&query).unwrap();
    assert_eq!(
        response.diagnostics.cache_misses, 3,
        "a differently-namespaced composer must miss on every term"
    );

    // And directly: the same question under two namespaces is two keys.
    let base = Digest::of_bytes(b"base");
    let a = EvaluationCache::new(LOCKFILE_A).key(&base, &[], &term("usage.cai"));
    let b = EvaluationCache::new(LOCKFILE_B).key(&base, &[], &term("usage.cai"));
    assert_ne!(a, b);
}

/// `cache_policy.read = false` forces a recomputation that still agrees.
#[test]
fn tst_land_001_cache_bypass_agrees_with_the_cached_answer() {
    let mut composer = Composer::new(local_registry(), LOCKFILE_A);
    let cached = composer
        .evaluate(&full_query("q_01ARZ3NDEKTSV4RRFFQ69G5FAV", ALL))
        .unwrap();

    let bypass =
        full_query("q_01ARZ3NDEKTSV4RRFFQ69G5FAV", ALL).with_cache_policy(CachePolicy::bypass());
    let recomputed = composer.evaluate(&bypass).unwrap();

    assert_eq!(recomputed.diagnostics.evaluations, 3);
    for name in ALL {
        assert_eq!(
            cached.get(&term(name)).unwrap().value.to_bits(),
            recomputed.get(&term(name)).unwrap().value.to_bits()
        );
    }
}

/// Edit order does not change a delta query's cache identity.
#[test]
fn tst_land_001_edit_order_is_not_part_of_the_key() {
    let cache = EvaluationCache::new(LOCKFILE_A);
    let base = Digest::of_bytes(b"base");
    let t = term("usage.cai");
    let forward = [
        Edit::new(1, Codon::parse("CTC").unwrap()),
        Edit::new(4, Codon::parse("GGT").unwrap()),
    ];
    let reversed = [forward[1], forward[0]];
    assert_eq!(
        cache.key(&base, &forward, &t),
        cache.key(&base, &reversed, &t)
    );
}

// --------------------------------------------------------------------------- //
// TST-LAND-002 — partial propagates as UNKNOWN
// --------------------------------------------------------------------------- //

/// A term the host cannot support is refused at submission, not dropped.
///
/// ADR-0004 makes this concrete: a host with only a genome assembly has no
/// tRNA data, so `burden.growth_cost` is unsupportable. Ingest records that in
/// the profile; this is the boundary where the record becomes an answer.
#[test]
fn tst_land_002_an_unsupported_term_is_refused_with_cf_host_003() {
    let mut composer = Composer::new(local_registry(), LOCKFILE_A).with_host_capability(
        HostCapability::supporting(["usage.cai".to_owned(), "usage.gc_content".to_owned()]),
    );

    let err = composer
        .evaluate(&full_query("q_01ARZ3NDEKTSV4RRFFQ69G5FAV", ALL))
        .unwrap_err();
    assert_eq!(err.code, "CF-HOST-003");
    assert!(err.detail.contains("constraint.violation_count"));

    composer
        .evaluate(&full_query(
            "q_01ARZ3NDEKTSV4RRFFQ69G5FAV",
            &["usage.cai", "usage.gc_content"],
        ))
        .expect("supported terms still compose");
}

/// A deferred term is absent from `terms`, listed in `terms_deferred`, and
/// `None` in the objective vector — never imputed.
#[test]
fn tst_land_002_a_deferred_term_becomes_unknown_and_is_never_imputed() {
    let mut registry = local_registry();
    registry
        .register(
            Box::new(support_provider::OverBudget::new("lm.naturalness")),
            TermRole::ParetoAxis,
        )
        .unwrap();
    let mut composer = Composer::new(registry, LOCKFILE_A);

    let requested = ["usage.cai", "lm.naturalness", "constraint.violation_count"];
    let response = composer
        .evaluate(&full_query("q_01ARZ3NDEKTSV4RRFFQ69G5FAV", &requested))
        .unwrap();

    assert_eq!(response.status, Status::Partial);
    assert!(response.get(&term("lm.naturalness")).is_none());
    assert_eq!(response.deferred(), &[term("lm.naturalness")]);

    let vector = response.objective_vector.as_ref().unwrap();
    assert_eq!(
        vector.order(),
        &[
            term("usage.cai"),
            term("lm.naturalness"),
            term("constraint.violation_count")
        ],
        "a deferred axis keeps its slot; dropping it would change the front's dimensionality"
    );
    assert_eq!(vector.values()[1], None);
    assert!(vector.has_unknown());

    // The surviving terms are unaffected: deferral is per term.
    assert!(response.get(&term("usage.cai")).is_some());
    assert!(response.get(&term("constraint.violation_count")).is_some());
}

// --------------------------------------------------------------------------- //
// TST-LAND-003 — no predictor crate in the tree
// --------------------------------------------------------------------------- //

/// `cf-landscape` depends on no predictor crate, transitively.
///
/// Computed from the manifests rather than by shelling out to `cargo tree`, so
/// it runs in the same conditions as every other test and cannot be skipped by
/// a missing toolchain.
#[test]
fn tst_land_003_no_predictor_crate_is_reachable() {
    const PREDICTORS: &[&str] = &[
        "cf-epistasis",
        "cf-biophys-structure",
        "cf-biophys-kinetics",
        "cf-biophys-lm",
        "cf-burden",
    ];

    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf();

    let mut direct: HashMap<String, BTreeSet<String>> = HashMap::new();
    for entry in std::fs::read_dir(&root).unwrap() {
        let dir = entry.unwrap().path();
        let manifest = dir.join("Cargo.toml");
        if !manifest.is_file() {
            continue;
        }
        let text = std::fs::read_to_string(&manifest).unwrap();
        let name = dir.file_name().unwrap().to_string_lossy().to_string();
        let deps = text
            .lines()
            .filter_map(|line| line.trim().strip_prefix("cf-"))
            .filter_map(|rest| rest.split_whitespace().next())
            .filter(|token| token.ends_with('=') || !token.contains('='))
            .map(|token| format!("cf-{}", token.trim_end_matches('=').trim()))
            .collect();
        direct.insert(name, deps);
    }

    let mut reachable = BTreeSet::new();
    let mut frontier = vec!["cf-landscape".to_string()];
    while let Some(next) = frontier.pop() {
        for dep in direct.get(&next).into_iter().flatten() {
            if reachable.insert(dep.clone()) {
                frontier.push(dep.clone());
            }
        }
    }

    assert!(
        reachable.contains("cf-landscape-api"),
        "the closure computation found nothing, so the assertion below is vacuous: {reachable:?}"
    );
    for predictor in PREDICTORS {
        assert!(
            !reachable.contains(*predictor),
            "cf-landscape reaches {predictor}; the swappable-layer property is gone. \
             Reachable: {reachable:?}"
        );
    }
    assert!(
        reachable.contains("cf-constraints"),
        "cf-constraints is expected and permitted here; see modules.yaml"
    );
}

// --------------------------------------------------------------------------- //
// TST-LAND-004 — degradation is per term, not per query
// --------------------------------------------------------------------------- //

/// One L0 term does not drag the others onto the rescore path.
///
/// This is the test the sequencing audit added. `TST-LAND-001`, `002` and
/// `003` all pass under a composition that responds to a single L0 term by
/// rescoring everything: the numbers come out right and only the cost changes.
/// `full_rescores` is what makes the difference observable.
#[test]
fn tst_land_004_a_single_l0_term_does_not_force_a_global_rescore() {
    let registry = mixed_registry();
    assert_eq!(
        registry
            .get(&term("constraint.violation_count"))
            .unwrap()
            .provider
            .capability(),
        Capability::L0,
        "the fixture must actually be mixed"
    );
    assert_eq!(
        registry
            .get(&term("usage.cai"))
            .unwrap()
            .provider
            .capability(),
        Capability::L1
    );

    let mut composer = Composer::new(registry, LOCKFILE_A);
    let base = composer.remember_base(&flanked_sequence()).unwrap();

    let query = delta_query(
        "q_01ARZ3NDEKTSV4RRFFQ69G5FAW",
        &base,
        ALL,
        vec![Edit::new(1, Codon::parse("CTC").unwrap())],
    );
    let response = composer.evaluate(&query).unwrap();

    assert_eq!(
        response.diagnostics.full_rescores, 1,
        "exactly the one L0 term should have been rescored, not all three"
    );
    assert_eq!(response.diagnostics.evaluations, 3);

    // And the values are right whichever path each term took.
    let mutated = flanked_sequence()
        .substituted(1, Codon::parse("CTC").unwrap())
        .unwrap();
    let host = host();
    let expected_cai = CaiTerm::new(usage_table("highly_expressed"))
        .unwrap()
        .score(&ctx(&mutated, &host))
        .unwrap();
    let expected_violations = ViolationCountTerm::new(
        ConstraintChecker::compile(&global_constraints()).unwrap(),
        layout(),
    )
    .unwrap()
    .score(&ctx(&mutated, &host))
    .unwrap();

    let got_cai = response.get(&term("usage.cai")).unwrap().value;
    assert!(
        (got_cai - expected_cai).abs() < 1e-12,
        "delta path gave {got_cai}, rescoring gives {expected_cai}"
    );
    assert_eq!(
        response
            .get(&term("constraint.violation_count"))
            .unwrap()
            .value,
        expected_violations
    );
}

/// With every term L1, nothing is rescored at all.
///
/// The other half of the claim: `full_rescores == 1` above is only meaningful
/// if the counter can reach zero.
#[test]
fn tst_land_004_an_all_l1_registry_rescores_nothing() {
    let mut composer = Composer::new(local_registry(), LOCKFILE_A);
    let base = composer.remember_base(&flanked_sequence()).unwrap();
    let response = composer
        .evaluate(&delta_query(
            "q_01ARZ3NDEKTSV4RRFFQ69G5FAW",
            &base,
            ALL,
            vec![Edit::new(1, Codon::parse("CTC").unwrap())],
        ))
        .unwrap();
    assert_eq!(response.diagnostics.full_rescores, 0);
}

/// A delta query naming a base the session has forgotten fails with the code
/// the contract specifies, rather than being silently promoted to a full query.
#[test]
fn tst_land_004_an_evicted_base_is_reported() {
    let mut composer = Composer::new(local_registry(), LOCKFILE_A);
    let base = composer.remember_base(&flanked_sequence()).unwrap();
    composer.forget_bases();

    let err = composer
        .evaluate(&delta_query(
            "q_01ARZ3NDEKTSV4RRFFQ69G5FAW",
            &base,
            &["usage.cai"],
            vec![Edit::new(1, Codon::parse("CTC").unwrap())],
        ))
        .unwrap_err();
    assert_eq!(err.code, "CF-LAND-010");
}

/// A delta table carries exact deltas for the L1 terms and omits the L0 term,
/// with `exact: false` recording that the table is not complete.
#[test]
fn tst_land_004_a_delta_table_omits_an_l0_term_rather_than_zeroing_it() {
    let mut composer = Composer::new(mixed_registry(), LOCKFILE_A);
    let base = composer.remember_base(&flanked_sequence()).unwrap();

    let query = LandscapeQuery::new(
        "q_01ARZ3NDEKTSV4RRFFQ69G5FAX",
        Target::DeltaTable {
            base_digest: base,
            positions: vec![1, 2, 3],
            cistron_id: None,
        },
        host(),
        terms(ALL),
        Budget::new(BudgetTier::Exact),
    )
    .with_outputs(Outputs::default().with_deltas());

    let response = composer.evaluate(&query).unwrap();
    let deltas = response.deltas.as_ref().unwrap();
    assert!(
        !deltas.exact,
        "one requested term could not contribute, so the table is not exact"
    );
    assert_eq!(deltas.entries.len(), 3);

    for entry in &deltas.entries {
        assert!(!entry.alternatives.is_empty());
        for alternative in &entry.alternatives {
            assert!(
                alternative
                    .get(&term("constraint.violation_count"))
                    .is_none(),
                "an L0 term must be absent, not present as 0.0 — 0.0 reads as 'no change'"
            );
            assert!(alternative.get(&term("usage.cai")).is_some());
            assert!(alternative.get(&term("usage.gc_content")).is_some());
        }
    }
}

// --------------------------------------------------------------------------- //
// TST-LAND-005 — a diagnostic term cannot change a ranking
// --------------------------------------------------------------------------- //

/// AMD-002, enforced by construction.
///
/// `usage.gc_content` caps itself at diagnostic, so it is evaluated, reported,
/// and absent from the objective vector. Requesting it alongside the axes
/// leaves the vector identical to the vector without it — not merely equal in
/// value, but the same order and the same length.
#[test]
fn tst_land_005_a_diagnostic_term_is_reported_and_never_ranked() {
    let mut composer = Composer::new(local_registry(), LOCKFILE_A);

    let without = composer
        .evaluate(&full_query(
            "q_01ARZ3NDEKTSV4RRFFQ69G5FAV",
            &["usage.cai", "constraint.violation_count"],
        ))
        .unwrap();
    let with = composer
        .evaluate(&full_query("q_01ARZ3NDEKTSV4RRFFQ69G5FAV", ALL))
        .unwrap();

    assert!(
        with.get(&term("usage.gc_content")).is_some(),
        "the diagnostic term is still computed and reported"
    );

    let a = without.objective_vector.as_ref().unwrap();
    let b = with.objective_vector.as_ref().unwrap();
    assert_eq!(a.order(), b.order());
    assert_eq!(a.values(), b.values());
    assert_eq!(a.sense(), b.sense());
    assert!(
        !b.order().contains(&term("usage.gc_content")),
        "a diagnostic term reached the objective vector"
    );
}

/// The folding-kinetics case specifically, since that is the term `AMD-002` is
/// about.
///
/// The fixture is deliberately loud: the diagnostic term reports a value far
/// outside every other term's range and declares `Minimize`, so if it reached
/// the objective vector it would not merely perturb the ranking, it would
/// dominate it. The composed vector is byte-identical with and without it.
///
/// `EXP-06` is what would promote this term to an axis. Until it runs, the
/// promotion is not a configuration change: the provider caps itself, so
/// registering it as an axis fails.
#[test]
fn tst_land_005_folding_kinetics_is_inert_until_exp_06() {
    let mut registry = local_registry();
    registry
        .register_default(Box::new(support_provider::LoudDiagnostic::new(
            "kinetics.pause_concordance",
        )))
        .unwrap();
    assert_eq!(
        registry.role_of(&term("kinetics.pause_concordance")),
        Some(TermRole::Diagnostic),
        "the provider caps itself; nothing in a job configuration can raise it"
    );

    let mut composer = Composer::new(registry, LOCKFILE_A);
    let axes_only = ["usage.cai", "constraint.violation_count"];
    let with_kinetics = [
        "usage.cai",
        "kinetics.pause_concordance",
        "constraint.violation_count",
    ];

    let before = composer
        .evaluate(&full_query("q_01ARZ3NDEKTSV4RRFFQ69G5FAV", &axes_only))
        .unwrap();
    let after = composer
        .evaluate(&full_query("q_01ARZ3NDEKTSV4RRFFQ69G5FAV", &with_kinetics))
        .unwrap();

    let reported = after.get(&term("kinetics.pause_concordance")).unwrap();
    assert_eq!(
        reported.value, -1.0e6,
        "the diagnostic value is still reported"
    );

    let a = before.objective_vector.as_ref().unwrap();
    let b = after.objective_vector.as_ref().unwrap();
    assert_eq!(a.order(), b.order());
    assert_eq!(a.values(), b.values());
    assert_eq!(a.sense(), b.sense());

    // And a second registry that tries to rank on it cannot be built at all.
    let mut refusing = TermRegistry::new();
    assert_eq!(
        refusing
            .register(
                Box::new(support_provider::LoudDiagnostic::new(
                    "kinetics.pause_concordance"
                )),
                TermRole::ParetoAxis,
            )
            .unwrap_err()
            .code,
        "CF-LAND-004"
    );
}

/// The type gate, independent of any composer.
#[test]
fn tst_land_005_a_diagnostic_provider_cannot_be_registered_as_an_axis() {
    let mut registry = TermRegistry::new();
    let refused = registry.register(Box::new(GcContentTerm::new()), TermRole::ParetoAxis);
    assert_eq!(refused.unwrap_err().code, "CF-LAND-004");
}

/// A term nobody registered fails by name; it is not skipped.
#[test]
fn tst_land_005_an_unregistered_term_fails_loudly() {
    let mut composer = Composer::new(local_registry(), LOCKFILE_A);
    let err = composer
        .evaluate(&full_query(
            "q_01ARZ3NDEKTSV4RRFFQ69G5FAV",
            &["usage.cai", "epistasis.dwell_deviation"],
        ))
        .unwrap_err();
    assert_eq!(err.code, "CF-LAND-001");
    assert!(err.detail.contains("epistasis.dwell_deviation"));
}

mod support_provider {
    use cf_landscape_api::{Capability, EvalContext, Provenance, TermProvider, TermRole};
    use cf_types::{code, CfError, Direction, Result, TermId, Unit};

    /// A diagnostic-capped provider with a value that would dominate any
    /// ranking it entered.
    ///
    /// The magnitude is the test: a diagnostic term whose value sat inside the
    /// other terms' range could leak into the objective vector and change
    /// nothing detectable, and the assertion would pass on a broken gate.
    #[derive(Debug)]
    pub struct LoudDiagnostic {
        term: TermId,
    }

    impl LoudDiagnostic {
        pub fn new(term: &str) -> Self {
            Self {
                term: TermId::new(term).unwrap(),
            }
        }
    }

    impl TermProvider for LoudDiagnostic {
        fn term_id(&self) -> &TermId {
            &self.term
        }
        fn unit(&self) -> Unit {
            Unit::new(Unit::DIMENSIONLESS)
        }
        fn direction(&self) -> Direction {
            Direction::Minimize
        }
        fn capability(&self) -> Capability {
            Capability::L0
        }
        fn max_role(&self) -> TermRole {
            TermRole::Diagnostic
        }
        fn provenance(&self) -> Provenance {
            Provenance::builtin("test_loud_diagnostic", "0.1.0")
        }
        fn score(&self, _ctx: &EvalContext<'_>) -> Result<f64> {
            Ok(-1.0e6)
        }
    }

    /// A provider that always misses its budget.
    ///
    /// Deterministic on purpose: driving the deferral path with a real wall
    /// clock would make the test flaky on a loaded machine, and a flaky test
    /// for an invariant this important gets muted.
    #[derive(Debug)]
    pub struct OverBudget {
        term: TermId,
    }

    impl OverBudget {
        pub fn new(term: &str) -> Self {
            Self {
                term: TermId::new(term).unwrap(),
            }
        }
    }

    impl TermProvider for OverBudget {
        fn term_id(&self) -> &TermId {
            &self.term
        }
        fn unit(&self) -> Unit {
            Unit::new(Unit::NATS_PER_CODON)
        }
        fn direction(&self) -> Direction {
            Direction::Maximize
        }
        fn capability(&self) -> Capability {
            Capability::L0
        }
        fn provenance(&self) -> Provenance {
            Provenance::builtin("test_over_budget", "0.1.0")
        }
        fn score(&self, _ctx: &EvalContext<'_>) -> Result<f64> {
            Err(CfError::new(
                code::BUDGET_EXCEEDED,
                "provider exceeded its latency budget",
            ))
        }
    }
}
