//! TST-METR-005 and TST-METR-006: NSGA-III normalization and niching.
//!
//! Ticket P2-03.
//!
//! These exist because the sequencing audit found `P2-03`'s two declared
//! acceptance criteria were both already green before a line of NSGA-III was
//! written — `TST-SRCH-007` from `P2-02`'s comparator and `TST-METR-002` from
//! `P2-01`'s Das–Dennis. Normalization, the degeneracy guards, association and
//! niching had no criterion at all, so the ticket could have been closed
//! without any of them.

use cf_metrics::{
    associate, das_dennis, niching_select, normalize, perpendicular_distance, placeable,
    Association, Normalization, Ranked,
};

// --------------------------------------------------------------------------- //
// TST-METR-005 — normalization
// --------------------------------------------------------------------------- //

#[test]
fn tst_metr_005_the_ideal_point_is_the_running_minimum() {
    let first = normalize(&[vec![4.0, 9.0], vec![9.0, 4.0]], None);
    assert_eq!(first.ideal, vec![4.0, 4.0]);

    // A later generation that lost its best member on an axis must not move
    // the origin back: the ideal point is defined over every generation seen,
    // and rescaling the whole population because one member died would make
    // distances incomparable between generations.
    let later = normalize(&[vec![7.0, 9.0], vec![9.0, 7.0]], Some(&first.ideal));
    assert_eq!(later.ideal, vec![4.0, 4.0]);

    // And it does fall when a genuinely better member appears.
    let better = normalize(&[vec![1.0, 9.0], vec![9.0, 2.0]], Some(&later.ideal));
    assert_eq!(better.ideal, vec![1.0, 2.0]);
}

#[test]
fn tst_metr_005_extreme_points_give_the_hyperplane_intercepts() {
    // A clean two-objective front on the line x + y = 10, translated so the
    // ideal is the origin: the intercepts are the axis crossings, 10 and 10.
    let points = vec![vec![0.0, 10.0], vec![5.0, 5.0], vec![10.0, 0.0]];
    let normalization = normalize(&points, None);

    assert!(
        !normalization.degenerate,
        "a well-formed front was called degenerate"
    );
    assert_eq!(normalization.ideal, vec![0.0, 0.0]);
    for intercept in &normalization.intercepts {
        assert!((intercept - 10.0).abs() < 1e-9, "got {intercept}");
    }

    // The extremes map to the unit axes and the middle to the middle.
    assert_eq!(normalization.apply(&[0.0, 10.0]), vec![0.0, 1.0]);
    assert_eq!(normalization.apply(&[10.0, 0.0]), vec![1.0, 0.0]);
    let mid = normalization.apply(&[5.0, 5.0]);
    assert!((mid[0] - 0.5).abs() < 1e-9 && (mid[1] - 0.5).abs() < 1e-9);
}

#[test]
fn tst_metr_005_a_singular_extreme_point_matrix_falls_back() {
    // Two perfectly correlated objectives: every point is on the line y = x,
    // so the extreme-point matrix is singular and there is no hyperplane. §3's
    // first guard. Without it the solve divides by ~0 and every subsequent
    // distance is meaningless in a way nothing downstream can detect.
    let points = vec![vec![1.0, 1.0], vec![2.0, 2.0], vec![3.0, 3.0]];
    let normalization = normalize(&points, None);

    assert!(
        normalization.degenerate,
        "a singular matrix was not detected"
    );
    assert!(normalization.intercepts.iter().all(|v| v.is_finite()));
    let applied = normalization.apply(&[2.0, 2.0]);
    assert!(applied.iter().all(|v| v.is_finite() && !v.is_nan()));
}

#[test]
fn tst_metr_005_a_degenerate_axis_normalizes_to_zero_rather_than_nan() {
    // Every member shares the second objective, so its span is zero.
    let points = vec![vec![1.0, 7.0], vec![3.0, 7.0]];
    let normalization = normalize(&points, None);
    let applied = normalization.apply(&[2.0, 7.0]);

    assert!(
        applied.iter().all(|v| !v.is_nan()),
        "a zero span produced NaN"
    );
    assert_eq!(applied[1], 0.0, "a degenerate axis must distinguish nobody");
}

#[test]
fn tst_metr_005_normalization_is_invariant_to_point_order() {
    let points = vec![vec![0.0, 10.0], vec![5.0, 5.0], vec![10.0, 0.0]];
    let mut reversed = points.clone();
    reversed.reverse();
    assert_eq!(normalize(&points, None), normalize(&reversed, None));
}

#[test]
fn tst_metr_005_an_empty_population_is_reported_degenerate_not_panicked() {
    let normalization = normalize(&[], None);
    assert!(normalization.degenerate);
    assert!(normalization.intercepts.is_empty());
}

// --------------------------------------------------------------------------- //
// TST-METR-006 — association and niching
// --------------------------------------------------------------------------- //

#[test]
fn tst_metr_006_a_point_on_a_direction_has_zero_perpendicular_distance() {
    let direction = vec![0.5, 0.5];
    assert!(perpendicular_distance(&[0.5, 0.5], &direction) < 1e-12);
    assert!(perpendicular_distance(&[2.0, 2.0], &direction) < 1e-12);
    // Off-axis is strictly positive, and the value is the true perpendicular.
    let off = perpendicular_distance(&[1.0, 0.0], &[1.0, 0.0]);
    assert!(off < 1e-12);
    assert!(perpendicular_distance(&[0.0, 1.0], &[1.0, 0.0]) > 0.9);
}

#[test]
fn tst_metr_006_association_picks_the_nearest_direction_and_ties_go_low() {
    let directions = das_dennis(2, 2).unwrap(); // (0,1), (0.5,0.5), (1,0)
    let points = vec![vec![0.0, 1.0], vec![0.5, 0.5], vec![1.0, 0.0]];
    let associations = associate(&points, &directions);

    let references: Vec<usize> = associations.iter().map(|a| a.reference).collect();
    assert_eq!(references.len(), 3);
    assert_eq!(
        references
            .iter()
            .collect::<std::collections::HashSet<_>>()
            .len(),
        3,
        "three points on three distinct directions collapsed onto fewer"
    );

    // A point equidistant from two directions takes the lower index, so
    // association is a function of the point set. DEP-SRCH-009.
    let symmetric = associate(&[vec![1.0, 1.0]], &[vec![1.0, 0.0], vec![0.0, 1.0]]);
    assert_eq!(symmetric[0].reference, 0);
}

#[test]
fn tst_metr_006_the_least_occupied_niche_is_filled_first() {
    // Reference 0 already has two members, reference 1 has none. The single
    // free slot must go to reference 1's candidate even though reference 0's
    // is closer to its direction.
    let associations = vec![
        Association {
            reference: 0,
            distance: 0.0,
        },
        Association {
            reference: 0,
            distance: 0.0,
        },
        Association {
            reference: 0,
            distance: 0.01,
        },
        Association {
            reference: 1,
            distance: 0.9,
        },
    ];
    let chosen = niching_select(&associations, &[0, 1], &[2, 3], 2, 1);

    assert_eq!(chosen, vec![3], "selection filled an already-crowded niche");
}

#[test]
fn tst_metr_006_within_a_niche_the_closest_member_is_taken() {
    let associations = vec![
        Association {
            reference: 0,
            distance: 0.9,
        },
        Association {
            reference: 0,
            distance: 0.1,
        },
        Association {
            reference: 0,
            distance: 0.5,
        },
    ];
    assert_eq!(
        niching_select(&associations, &[], &[0, 1, 2], 1, 1),
        vec![1]
    );
}

#[test]
fn tst_metr_006_every_tie_is_broken_deterministically_not_by_insertion_order() {
    // Three candidates, identical distances, identical niches. Whichever is
    // chosen, it must be the same one however the candidate list is ordered —
    // a tie-break that fell through to insertion order would make selection a
    // function of how the population happened to be built.
    let associations = vec![
        Association {
            reference: 0,
            distance: 0.5,
        },
        Association {
            reference: 0,
            distance: 0.5,
        },
        Association {
            reference: 0,
            distance: 0.5,
        },
    ];
    let orders: [&[usize]; 3] = [&[0, 1, 2], &[2, 1, 0], &[1, 2, 0]];
    let chosen: Vec<Vec<usize>> = orders
        .iter()
        .map(|order| niching_select(&associations, &[], order, 1, 1))
        .collect();

    assert!(
        chosen.windows(2).all(|w| w[0] == w[1]),
        "a tie resolved by insertion order: {chosen:?}"
    );
    assert_eq!(chosen[0], vec![0]);
}

#[test]
fn tst_metr_006_more_members_than_slots_in_one_niche_is_handled() {
    // The boundary case: four candidates all associated with one reference,
    // two slots. Selection must take two, not spin and not overfill.
    let associations = vec![
        Association {
            reference: 0,
            distance: 0.4,
        },
        Association {
            reference: 0,
            distance: 0.1,
        },
        Association {
            reference: 0,
            distance: 0.3,
        },
        Association {
            reference: 0,
            distance: 0.2,
        },
    ];
    let chosen = niching_select(&associations, &[], &[0, 1, 2, 3], 3, 2);

    assert_eq!(chosen.len(), 2);
    assert_eq!(chosen, vec![1, 3], "the two closest were not taken");
}

#[test]
fn tst_metr_006_an_empty_reference_point_does_not_stall_selection() {
    // Reference 2 has no candidate at all. Without excluding it, the
    // least-occupied rule keeps choosing it and the loop never terminates.
    let associations = vec![
        Association {
            reference: 0,
            distance: 0.2,
        },
        Association {
            reference: 1,
            distance: 0.3,
        },
    ];
    let chosen = niching_select(&associations, &[], &[0, 1], 3, 2);
    assert_eq!(chosen, vec![0, 1]);
}

#[test]
fn tst_metr_006_selection_returns_at_most_the_slots_available() {
    let associations = vec![
        Association {
            reference: 0,
            distance: 0.2,
        },
        Association {
            reference: 1,
            distance: 0.3,
        },
        Association {
            reference: 2,
            distance: 0.4,
        },
    ];
    assert_eq!(
        niching_select(&associations, &[], &[0, 1, 2], 3, 0).len(),
        0
    );
    assert_eq!(
        niching_select(&associations, &[], &[0, 1, 2], 3, 5).len(),
        3
    );
}

#[test]
fn tst_metr_006_the_chosen_set_is_returned_in_index_order() {
    let associations = vec![
        Association {
            reference: 2,
            distance: 0.1,
        },
        Association {
            reference: 1,
            distance: 0.1,
        },
        Association {
            reference: 0,
            distance: 0.1,
        },
    ];
    let chosen = niching_select(&associations, &[], &[0, 1, 2], 3, 3);
    assert_eq!(
        chosen,
        vec![0, 1, 2],
        "the caller could depend on selection order"
    );
}

// --------------------------------------------------------------------------- //
// UNKNOWN through the NSGA-III path specifically
// --------------------------------------------------------------------------- //

#[test]
fn tst_metr_006_an_unknown_member_is_not_placed_in_objective_space() {
    // An individual with an unknown term has no coordinate, so associating it
    // with a direction would be inventing one. It is reported separately for
    // the caller to place at the back of its feasibility class, per §4 — not
    // silently given a position that would let niching select on a number
    // nobody measured.
    let population = vec![
        Ranked::feasible(vec![Some(1.0), Some(1.0)]),
        Ranked::feasible(vec![Some(2.0), None]),
        Ranked::feasible(vec![Some(3.0), Some(0.5)]),
    ];
    let (known, unknown) = placeable(&population, &[0, 1, 2]);

    assert_eq!(known, vec![0, 2]);
    assert_eq!(unknown, vec![1], "an unknown member was given a coordinate");
}

#[test]
fn tst_metr_006_an_infeasible_member_is_still_placeable() {
    // Feasibility is the caller's ordering concern, not a coordinate problem:
    // an infeasible individual has real objective values and can be placed.
    // Conflating the two would drop infeasible individuals from selection,
    // which §4 forbids.
    let population = vec![Ranked::infeasible(vec![Some(1.0), Some(1.0)], 2.0)];
    let (known, unknown) = placeable(&population, &[0]);
    assert_eq!(known, vec![0]);
    assert!(unknown.is_empty());
}

#[test]
fn tst_metr_005_a_hand_built_zero_span_normalizes_to_zero_rather_than_nan() {
    // `normalize`'s own fallback never produces a zero span - it sets an
    // intercept of `ideal + 1` when the per-axis maximum is zero - so this
    // guard is unreachable through it and a control removing the guard does
    // not fire. It protects a `Normalization` a caller built directly, which
    // is what this constructs. Found by a negative control that stayed silent
    // for a reason worth writing down rather than annotating away.
    let normalization = Normalization {
        ideal: vec![5.0, 5.0],
        intercepts: vec![5.0, 9.0],
        degenerate: true,
    };
    let applied = normalization.apply(&[5.0, 7.0]);

    assert!(
        applied.iter().all(|v| !v.is_nan()),
        "a zero span produced NaN"
    );
    assert_eq!(applied[0], 0.0);
    assert!((applied[1] - 0.5).abs() < 1e-12);
}
