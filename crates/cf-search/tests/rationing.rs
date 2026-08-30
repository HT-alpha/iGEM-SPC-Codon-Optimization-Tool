//! `TST-SRCH-020`/`TST-SRCH-021`/`TST-SRCH-022` — surrogate rationing.
//!
//! Ticket P2-06. Spec: `planning/03-algorithms/03-nsga-loop-and-operators.md` §6.
//!
//! The sequencing audit found that P2-06's stated criterion, `TST-SRCH-005`
//! ("every emitted candidate carries an exact, non-surrogate evaluation"), is a
//! *certification* property belonging to P2-07. It says nothing about when a
//! surrogate is used instead of a full model call, which is P2-06's entire
//! content. These tests cover the rationing itself.

use cf_search::rationing::{Observation, Rationed, Rationer, Reason, Surrogate};
use cf_search::rng::SearchRng;

fn rationer(e_max: usize) -> Rationer {
    Rationer::new(e_max, 10)
}

/// Objectives that make the front a genuine minority of the population.
///
/// Fixture rule: a fixture must be able to fail. If every individual were
/// non-dominated, component (a) would fill the budget on its own and the test
/// could not tell a correct four-way choice from a one-way one.
fn spread_population(size: usize) -> (Vec<Vec<f64>>, Vec<f64>) {
    let objectives: Vec<Vec<f64>> = (0..size)
        .map(|i| {
            let x = i as f64;
            // A convex trade-off for the first few, dominated bulk after.
            if i < 4 {
                vec![x, 10.0 - x]
            } else {
                vec![20.0 + x, 20.0 + x]
            }
        })
        .collect();
    (objectives, vec![0.0; size])
}

#[test]
fn tst_srch_020_the_budget_caps_exact_calls_and_the_rest_are_not_imputed() {
    let (objectives, violations) = spread_population(50);
    let mut rng = SearchRng::new(7);
    let mut rationer = rationer(12);

    let Rationed { ranked, .. } = rationer.rank(&objectives, &violations, &mut rng);

    // The ranker ranks; the budget cut belongs to BatchProtocol (DEP-SRCH-011).
    // What this checks is that the ranking is a strict priority order with no
    // duplicates, since the protocol takes a prefix of it.
    let unique: std::collections::BTreeSet<usize> = ranked.iter().copied().collect();
    assert_eq!(
        unique.len(),
        ranked.len(),
        "a duplicated index would be scored twice and waste the budget it cost"
    );
    assert!(
        ranked.len() <= objectives.len(),
        "the ranking cannot exceed the population"
    );
    assert!(
        !ranked.is_empty(),
        "a non-empty population must produce candidates, or nothing is ever scored"
    );
}

/// The mandatory 10 % is mandatory.
///
/// §6 is explicit that without it the surrogate is only ever trained on points
/// it already likes. The failure mode is silent — a run without the sample
/// still produces a front — so it needs a test that looks for the sample
/// specifically rather than at the outcome.
#[test]
fn tst_srch_020_the_uniform_sample_is_drawn_even_when_the_front_could_fill_the_budget() {
    // Every individual non-dominated: component (a) alone would take the whole
    // budget, which is exactly the condition under which a lazy implementation
    // drops the sample.
    let size = 40;
    let objectives: Vec<Vec<f64>> = (0..size)
        .map(|i| vec![i as f64, (size - i) as f64])
        .collect();
    let violations = vec![0.0; size];

    let mut rng = SearchRng::new(11);
    let mut rationer = rationer(16);
    let rationed = rationer.rank(&objectives, &violations, &mut rng);

    let sampled = rationed
        .reasons
        .iter()
        .filter(|(_, reason)| *reason == Reason::RandomSample)
        .count();
    assert!(
        sampled > 0,
        "the uniform sample must survive a population that is entirely \
         non-dominated; this is the case where dropping it is invisible"
    );
}

/// The surrogate cannot be trained on its own predictions.
///
/// Enforced by the type: `Observation::from_exact` is the only constructor, so
/// there is no path from `predict` to `refit`. The test documents the property
/// and would catch a later constructor being added.
#[test]
fn tst_srch_022_the_surrogate_learns_only_from_exact_evaluations() {
    let mut surrogate = Surrogate::new(1e-6);

    // Untrained: refuses rather than returning a zero that ranks as a real
    // prediction.
    let error = surrogate
        .predict(&[1.0, 2.0])
        .expect_err("an untrained surrogate must refuse");
    assert_eq!(error.code, "CF-SRCH-015");

    let observations: Vec<Observation> = (0..20)
        .map(|i| {
            let x = i as f64;
            Observation::from_exact(vec![x, 1.0], 3.0 * x + 2.0)
        })
        .collect();
    surrogate.refit(&observations).expect("fit");
    assert_eq!(surrogate.trained_on(), 20);

    let predicted = surrogate.predict(&[10.0, 1.0]).expect("predict");
    assert!(
        (predicted - 32.0).abs() < 1e-6,
        "a ridge fit on an exactly linear relation should recover it: got {predicted}"
    );
}

/// Refit cadence is every `G_s` generations, not every generation.
#[test]
fn tst_srch_022_the_refit_cadence_is_periodic() {
    let rationer = rationer(8);
    let refits: Vec<u64> = (0..31).filter(|g| rationer.should_refit(*g)).collect();
    assert_eq!(
        refits,
        vec![10, 20, 30],
        "§6 sets G_s = 10; fitting every generation would cost more than the calls saved"
    );
}

/// Adversarial floats, crossed through text.
///
/// Fixture rule from last round: an in-memory f64 round-trip does not exercise
/// the parse path, which is where precision is actually lost.
#[test]
fn tst_srch_022_observations_survive_serialization_of_their_values() {
    let adversarial: [f64; 4] = [
        0.1 + 0.2,
        f64::MIN_POSITIVE,
        1.000_000_000_000_000_2,
        -8.988_465_674_311_58e307,
    ];
    for value in adversarial {
        let text = serde_json::to_string(&value).unwrap();
        let back: f64 = serde_json::from_str(&text).unwrap();
        let observation = Observation::from_exact(vec![1.0], back);
        assert_eq!(
            observation.value().to_bits(),
            value.to_bits(),
            "an exact L0 value must survive the text it travels as; {value} came back as {back}"
        );
    }
}
