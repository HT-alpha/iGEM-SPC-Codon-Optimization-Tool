"""TST-DEC-003. Front reduction keeps every per-objective extreme.

A front of one point passes this vacuously, and so does a cap large enough to
keep everyone. The fixture has three axis extremes that lose on the other
axes, plus interior points, and a cap of two. Dropping the extreme-point
union leaves at least one extreme out.
"""

from __future__ import annotations

from dna_design_suite.codon.decision.reduce import reduce_candidates


def _candidate(cid: str, vector: list[float]) -> dict:
    return {"candidate_id": cid, "objective_vector": vector}


def test_reduction_keeps_every_axis_extreme_past_the_cap() -> None:
    senses = ["minimize", "minimize", "minimize"]
    candidates = [
        _candidate("axis0", [0.0, 10.0, 10.0]),
        _candidate("axis1", [10.0, 0.0, 10.0]),
        _candidate("axis2", [10.0, 10.0, 0.0]),
        _candidate("mid_a", [4.0, 4.0, 4.0]),
        _candidate("mid_b", [3.0, 5.0, 4.0]),
        _candidate("mid_c", [5.0, 3.0, 4.0]),
        _candidate("mid_d", [4.0, 4.0, 3.0]),
    ]
    presented, reduction = reduce_candidates(
        candidates,
        senses,
        "reference_direction_representatives",
        max_candidates=2,
    )
    ids = {item["candidate_id"] for item in presented}
    assert {"axis0", "axis1", "axis2"} <= ids
    assert reduction["raw_nondominated_count"] == len(candidates)
    assert reduction["presented_count"] == len(presented)
    assert reduction["method"] == "reference_direction_representatives"
    assert len(presented) >= 3


def test_clustering_keeps_an_unknown_extreme_and_still_returns() -> None:
    """An extreme with an UNKNOWN axis is protected but not a distance anchor."""
    senses = ["minimize", "minimize"]
    candidates = [
        _candidate("axis0", [0.0, None]),
        _candidate("axis1", [None, 0.0]),
        _candidate("mid", [1.0, 1.0]),
        _candidate("far", [2.0, 3.0]),
    ]
    presented, reduction = reduce_candidates(candidates, senses, "clustering", max_candidates=3)
    ids = {item["candidate_id"] for item in presented}
    assert {"axis0", "axis1"} <= ids
    assert reduction["presented_count"] == len(presented)
