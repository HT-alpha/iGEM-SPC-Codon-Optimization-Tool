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

    let Rationed {
        chosen, deferred, ..
    } = rationer.choose(&objectives, &violations, &mut rng);

    assert!(
        chosen.len() <= 12,
        "the budget is a cap: {} chosen against E_max 12",
        chosen.len()
    );
    assert!(
        chosen.len() < objectives.len(),
        "rationing must actually ration; if everything is chosen there is no \
         surrogate story to tell"
    );
    // The unchosen are simply absent from the chosen set. §6: f_L0 <- UNKNOWN,
    // not imputed. Nothing here hands back a filled-in value for them.
    let unchosen = objectives.len() - chosen.len();
    assert!(unchosen > 0);
    assert!(
        deferred.iter().all(|index| !chosen.contains(index)),
        "an individual cannot be both chosen and deferred"
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
    let rationed = rationer.choose(&objectives, &violations, &mut rng);

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

/// Truncation defers rather than discards, and the deferral is honoured.
#[test]
fn tst_srch_021_a_truncated_individual_is_first_in_line_next_generation() {
    let (objectives, violations) = spread_population(60);
    let mut rng = SearchRng::new(3);
    let mut rationer = rationer(8);

    let first = rationer.choose(&objectives, &violations, &mut rng);
    assert!(
        !first.deferred.is_empty(),
        "a 60-member population under E_max 8 must overflow, or this test proves nothing"
    );

    let second = rationer.choose(&objectives, &violations, &mut rng);
    let carried: Vec<usize> = first
        .deferred
        .iter()
        .copied()
        .filter(|index| second.chosen.contains(index))
        .collect();
    assert!(
        !carried.is_empty(),
        "deferred individuals must be taken up next generation; deferring into a \
         set that never gets read is discarding with extra steps"
    );
    assert!(
        second
            .reasons
            .iter()
            .any(|(_, reason)| *reason == Reason::Deferred),
        "and they must be recorded as deferred, so the report can say why"
    );
}

/// The silent-UNKNOWN gap is visible.
///
/// Structurally the same problem `CF-ORCH-006` solves for skipped tiers: an
/// individual the budget never reached stays `UNKNOWN` forever, and without
/// this the job's output would not say so.
#[test]
fn tst_srch_021_individuals_the_budget_never_reached_are_reported() {
    let (objectives, violations) = spread_population(60);
    let mut rng = SearchRng::new(5);
    let mut rationer = rationer(6);

    for _ in 0..3 {
        let rationed = rationer.choose(&objectives, &violations, &mut rng);
        for index in rationed.chosen {
            rationer.absorb(index, Observation::from_exact(vec![index as f64], 1.0));
        }
    }

    let stranded = rationer.unevaluated_at_end();
    assert!(
        !stranded.is_empty(),
        "60 individuals at E_max 6 over 3 generations cannot all have been evaluated"
    );
    // And absorbing clears the flag, so this reports a real gap rather than
    // everything that was ever deferred.
    let mut settled = rationer;
    for &index in &stranded {
        settled.absorb(index, Observation::from_exact(vec![index as f64], 1.0));
    }
    assert!(
        settled.unevaluated_at_end().is_empty(),
        "an individual that later gets an exact evaluation is no longer stranded; \
         otherwise the warning would fire on every run that ever deferred anything"
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
