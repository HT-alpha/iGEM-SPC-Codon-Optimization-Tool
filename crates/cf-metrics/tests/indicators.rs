//! TST-METR-001 through TST-METR-004: the indicator math.
//!
//! Ticket P2-01.
//!
//! `TST-METR-001` and `TST-METR-003` are checked against brute force rather
//! than against a second closed form, because a second closed form written by
//! the same hand agrees with the first for the same reasons it is wrong.

use cf_metrics::{
    constraint_dominance, crowding_distance, das_dennis, das_dennis_two_layer, epsilon_additive,
    fast_nondominated_sort, hypervolume, igd_plus, Dominance, Ranked,
};

/// splitmix64, so the fixtures replay exactly.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    /// A coarse grid on purpose: ties are the interesting case, and a uniform
    /// double never produces one.
    fn coarse(&mut self, levels: u64) -> f64 {
        (self.next() % levels) as f64
    }
}

fn feasible(values: &[f64]) -> Ranked {
    Ranked::feasible(values.iter().map(|v| Some(*v)).collect())
}

// --------------------------------------------------------------------------- //
// TST-METR-001 — sorting and hypervolume against brute force
// --------------------------------------------------------------------------- //

/// Brute-force front assignment: rank 0 is undominated, rank k is undominated
/// once ranks below k are removed.
fn brute_force_fronts(population: &[Ranked]) -> Vec<Vec<usize>> {
    let mut remaining: Vec<usize> = (0..population.len()).collect();
    let mut fronts = Vec::new();
    while !remaining.is_empty() {
        let front: Vec<usize> = remaining
            .iter()
            .copied()
            .filter(|&i| {
                !remaining.iter().any(|&j| {
                    j != i
                        && constraint_dominance(&population[j], &population[i]) == Dominance::Left
                })
            })
            .collect();
        assert!(
            !front.is_empty(),
            "a non-empty set with no undominated member"
        );
        remaining.retain(|i| !front.contains(i));
        fronts.push(front);
    }
    fronts
}

#[test]
fn tst_metr_001_nondominated_sorting_agrees_with_brute_force() {
    let mut rng = Rng(20260829);
    let mut with_ties = 0;

    for _ in 0..10_000 {
        let objectives = 2 + rng.below(4); // M in 2..=5
        let size = 1 + rng.below(30); // N in 1..=30
        let population: Vec<Ranked> = (0..size)
            .map(|_| {
                let values: Vec<f64> = (0..objectives).map(|_| rng.coarse(4)).collect();
                if rng.below(4) == 0 {
                    Ranked::infeasible(
                        values.iter().map(|v| Some(*v)).collect(),
                        rng.coarse(3) + 1.0,
                    )
                } else {
                    feasible(&values)
                }
            })
            .collect();

        let fast = fast_nondominated_sort(&population);
        let brute = brute_force_fronts(&population);
        assert_eq!(fast.len(), brute.len());
        for (a, b) in fast.iter().zip(brute.iter()) {
            let mut a = a.clone();
            let mut b = b.clone();
            a.sort_unstable();
            b.sort_unstable();
            assert_eq!(a, b, "front membership differs from brute force");
        }
        if fast.first().map(|f| f.len()).unwrap_or(0) > 1 {
            with_ties += 1;
        }
    }

    // The coarse grid must actually produce mutually non-dominated sets, or
    // this is 10 000 instances of a total order. `SYN-01`'s degenerate
    // configuration exists for the same reason.
    assert!(
        with_ties > 5000,
        "only {with_ties} of 10 000 fixtures had a multi-member first front; \
         the generator is producing near-total orders and the sort is barely tested"
    );
}

#[test]
fn tst_metr_001_hypervolume_agrees_with_a_grid_count() {
    // A Monte-Carlo-free brute force: count dominated cells on an integer grid
    // and compare with the exact volume. Only sound because the fixtures are
    // integer-valued, which is why the generator is coarse.
    let mut rng = Rng(4242);
    for _ in 0..300 {
        let objectives = 2 + rng.below(2); // M in 2..=3
        let size = 1 + rng.below(6);
        let reference = vec![5.0; objectives];
        let front: Vec<Vec<f64>> = (0..size)
            .map(|_| (0..objectives).map(|_| rng.coarse(5)).collect())
            .collect();

        let exact = hypervolume(&front, &reference).unwrap();

        let mut cells = 0u64;
        let side = 5usize;
        let total = side.pow(objectives as u32);
        for cell in 0..total {
            let mut coordinate = Vec::with_capacity(objectives);
            let mut rest = cell;
            for _ in 0..objectives {
                coordinate.push((rest % side) as f64);
                rest /= side;
            }
            if front
                .iter()
                .any(|p| p.iter().zip(&coordinate).all(|(v, c)| *v <= *c))
            {
                cells += 1;
            }
        }
        assert!(
            (exact - cells as f64).abs() < 1e-9,
            "hypervolume {exact} against {cells} dominated unit cells for {front:?}"
        );
    }
}

// --------------------------------------------------------------------------- //
// TST-METR-003 — indicator properties
// --------------------------------------------------------------------------- //

#[test]
fn tst_metr_003_hypervolume_is_invariant_to_candidate_ordering() {
    let mut rng = Rng(77);
    for _ in 0..500 {
        let objectives = 2 + rng.below(2);
        let size = 2 + rng.below(6);
        let reference = vec![6.0; objectives];
        let mut front: Vec<Vec<f64>> = (0..size)
            .map(|_| (0..objectives).map(|_| rng.coarse(5)).collect())
            .collect();

        let first = hypervolume(&front, &reference).unwrap();
        for _ in 0..4 {
            // A deterministic shuffle: rotate and swap, so a failure replays.
            let shift = 1 + rng.below(front.len());
            front.rotate_left(shift);
            let index = rng.below(front.len());
            front.swap(0, index);
            let again = hypervolume(&front, &reference).unwrap();
            assert_eq!(
                first.to_bits(),
                again.to_bits(),
                "hypervolume changed with candidate order: {first} vs {again}"
            );
        }
    }
}

#[test]
fn tst_metr_003_hypervolume_is_monotone_under_front_improvement() {
    let reference = vec![10.0, 10.0];
    let front = vec![vec![4.0, 6.0], vec![6.0, 4.0]];
    let base = hypervolume(&front, &reference).unwrap();

    // Improving one point strictly increases the volume.
    let improved = vec![vec![3.0, 6.0], vec![6.0, 4.0]];
    assert!(hypervolume(&improved, &reference).unwrap() > base);

    // Adding a dominated point changes nothing.
    let padded = vec![vec![4.0, 6.0], vec![6.0, 4.0], vec![7.0, 8.0]];
    assert_eq!(
        hypervolume(&padded, &reference).unwrap().to_bits(),
        base.to_bits()
    );

    // A point outside the reference contributes nothing rather than a negative.
    let outside = vec![vec![4.0, 6.0], vec![6.0, 4.0], vec![20.0, 20.0]];
    assert_eq!(
        hypervolume(&outside, &reference).unwrap().to_bits(),
        base.to_bits()
    );
}

#[test]
fn tst_metr_003_igd_plus_is_zero_exactly_when_the_front_is_reached() {
    let reference = vec![vec![1.0, 3.0], vec![2.0, 2.0], vec![3.0, 1.0]];

    assert_eq!(igd_plus(&reference, &reference).unwrap(), 0.0);

    // Strictly better than the reference is still zero: IGD+ measures only the
    // components in which the approximation is worse. Plain IGD would charge
    // for the distance and reward a front for stopping exactly on target.
    let better: Vec<Vec<f64>> = reference
        .iter()
        .map(|p| p.iter().map(|v| v - 0.5).collect())
        .collect();
    assert_eq!(igd_plus(&better, &reference).unwrap(), 0.0);

    let worse = vec![vec![4.0, 4.0]];
    assert!(igd_plus(&worse, &reference).unwrap() > 0.0);
}

#[test]
fn tst_metr_003_the_epsilon_indicator_is_zero_on_itself_and_positive_when_worse() {
    let reference = vec![vec![1.0, 3.0], vec![3.0, 1.0]];
    assert_eq!(epsilon_additive(&reference, &reference).unwrap(), 0.0);

    let shifted: Vec<Vec<f64>> = reference
        .iter()
        .map(|p| p.iter().map(|v| v + 2.0).collect())
        .collect();
    assert!((epsilon_additive(&shifted, &reference).unwrap() - 2.0).abs() < 1e-12);
}

// --------------------------------------------------------------------------- //
// TST-METR-002 — Das-Dennis directions
// --------------------------------------------------------------------------- //

fn binomial(n: usize, k: usize) -> usize {
    (0..k).fold(1usize, |acc, i| acc * (n - i) / (i + 1))
}

#[test]
fn tst_metr_002_das_dennis_counts_sum_and_are_unique() {
    for objectives in 2..=5usize {
        for divisions in 1..=6usize {
            let directions = das_dennis(objectives, divisions).unwrap();
            assert_eq!(
                directions.len(),
                binomial(objectives + divisions - 1, divisions),
                "M={objectives} p={divisions}"
            );
            for direction in &directions {
                assert_eq!(direction.len(), objectives);
                let sum: f64 = direction.iter().sum();
                assert!((sum - 1.0).abs() < 1e-12, "direction sums to {sum}");
                assert!(direction.iter().all(|v| *v >= 0.0));
            }
            let mut seen: Vec<Vec<u64>> = directions
                .iter()
                .map(|d| d.iter().map(|v| (v * 1e9).round() as u64).collect())
                .collect();
            seen.sort();
            let count = seen.len();
            seen.dedup();
            assert_eq!(seen.len(), count, "duplicate directions at M={objectives}");
        }
    }
}

#[test]
fn tst_metr_002_the_two_layer_construction_covers_the_interior() {
    let objectives = 5;
    let single = das_dennis(objectives, 3).unwrap();
    let two_layer = das_dennis_two_layer(objectives, 3, 2).unwrap();

    let interior = |set: &[Vec<f64>]| set.iter().filter(|d| d.iter().all(|v| *v > 1e-9)).count();

    // A single layer at a usable p puts almost everything on the boundary in
    // high dimension, which is the reason §3 gives for the two-layer scheme.
    assert_eq!(
        interior(&single),
        0,
        "the single layer already covers the interior"
    );
    assert!(
        interior(&two_layer) > 0,
        "the two-layer construction covered no interior point"
    );
    assert!(two_layer.len() > single.len());
    for direction in &two_layer {
        assert!((direction.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }
}

// --------------------------------------------------------------------------- //
// TST-METR-004 — crowding distance, which had no acceptance test
// --------------------------------------------------------------------------- //

#[test]
fn tst_metr_004_boundary_individuals_take_infinity() {
    let population = vec![
        feasible(&[1.0, 5.0]),
        feasible(&[2.0, 3.0]),
        feasible(&[5.0, 1.0]),
    ];
    let distance = crowding_distance(&population, &[0, 1, 2]);

    assert!(distance[0].is_infinite(), "the low extreme was crowded out");
    assert!(
        distance[2].is_infinite(),
        "the high extreme was crowded out"
    );
    assert!(distance[1].is_finite() && distance[1] > 0.0);
}

#[test]
fn tst_metr_004_a_degenerate_objective_contributes_zero_rather_than_dividing_by_zero() {
    // Every member shares the second objective, so its range is zero.
    let population = vec![
        feasible(&[1.0, 7.0]),
        feasible(&[2.0, 7.0]),
        feasible(&[3.0, 7.0]),
    ];
    let distance = crowding_distance(&population, &[0, 1, 2]);

    assert!(
        distance.iter().all(|d| !d.is_nan()),
        "a zero range produced NaN"
    );
    assert!(distance[0].is_infinite() && distance[2].is_infinite());
    // The middle takes its spacing from the first objective only.
    assert!(distance[1].is_finite());
    assert!((distance[1] - 1.0).abs() < 1e-12, "got {}", distance[1]);
}

#[test]
fn tst_metr_004_a_fully_degenerate_front_is_all_infinite_and_not_nan() {
    let population = vec![feasible(&[4.0, 4.0]), feasible(&[4.0, 4.0])];
    let distance = crowding_distance(&population, &[0, 1]);
    assert!(distance.iter().all(|d| !d.is_nan()));
    assert!(distance.iter().all(|d| d.is_infinite()));
}

#[test]
fn tst_metr_004_crowding_distance_is_invariant_to_input_ordering() {
    let population = vec![
        feasible(&[1.0, 5.0]),
        feasible(&[2.0, 3.0]),
        feasible(&[3.0, 2.0]),
        feasible(&[5.0, 1.0]),
    ];
    let forward = crowding_distance(&population, &[0, 1, 2, 3]);
    let reversed = crowding_distance(&population, &[3, 2, 1, 0]);

    for (index, value) in forward.iter().enumerate() {
        let mirrored = reversed[3 - index];
        assert_eq!(
            value.to_bits(),
            mirrored.to_bits(),
            "crowding distance changed with front order at {index}"
        );
    }
}

#[test]
fn tst_metr_004_an_unknown_term_neither_crowds_nor_uncrowds() {
    let population = vec![
        Ranked::feasible(vec![Some(1.0), None]),
        Ranked::feasible(vec![Some(2.0), Some(3.0)]),
        Ranked::feasible(vec![Some(3.0), Some(1.0)]),
    ];
    let distance = crowding_distance(&population, &[0, 1, 2]);
    assert!(distance.iter().all(|d| !d.is_nan()));
    // Two known values in objective 1 means both are boundaries there.
    assert!(distance[1].is_infinite());
}

#[test]
fn tst_metr_004_a_single_member_front_is_infinite() {
    let population = vec![feasible(&[1.0, 1.0])];
    assert!(crowding_distance(&population, &[0])[0].is_infinite());
}
