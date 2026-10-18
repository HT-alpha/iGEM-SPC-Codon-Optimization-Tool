"""Front reduction. `P2-15`, `AMD-008`, `TST-DEC-003`.

The stored front may be the full non-dominated set (`method: "none"`).
Presentation picks a smaller set and records how. Extreme points are kept
even when that exceeds `max_candidates`: a reduced front that hides the best
value on an axis is the failure this module exists to make structural.

Directions come from `cf-metrics` (`indicator_das_dennis`). A second generator
here would be another implementation of the same simplex.
"""

from __future__ import annotations

import json
from typing import Any

from dna_design_suite import _core

_METHODS = {
    "none",
    "reference_direction_representatives",
    "knee_selection",
    "clustering",
    "epsilon_grid",
}


def reduce_candidates(
    candidates: list[dict[str, Any]],
    senses: list[str],
    method: str,
    max_candidates: int = 20,
    recommended_id: str | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Return the presented candidates and the `reduction` object.

    `recommended_id`, when it is in the set, stays at index 0 so a reader of
    `candidates[0]` and a reader of `recommended` still name one design.
    """
    if method not in _METHODS:
        raise ValueError(f"reduction method {method!r} is not one of {sorted(_METHODS)}")
    if max_candidates < 1:
        raise ValueError("max_candidates must be at least 1")

    raw = len(candidates)
    if method == "none" or raw == 0:
        return list(candidates), _record(raw, len(candidates), method, max_candidates)

    protected = _ids(extreme_ids(candidates, senses))
    if recommended_id is not None and any(c["candidate_id"] == recommended_id for c in candidates):
        protected.add(recommended_id)

    if method == "reference_direction_representatives":
        chosen = _by_directions(candidates, senses, max_candidates, protected)
    elif method == "knee_selection":
        chosen = _knee(candidates, senses, max_candidates, protected)
    elif method == "clustering":
        chosen = _cluster(candidates, senses, max_candidates, protected)
    else:
        chosen = _epsilon_grid(candidates, senses, max_candidates, protected)

    ordered = _order(candidates, chosen, recommended_id)
    return ordered, _record(raw, len(ordered), method, max_candidates)


def present(
    document: dict[str, Any],
    method: str,
    max_candidates: int = 20,
) -> dict[str, Any]:
    """A copy of a front document with reduction applied and the digest redone."""
    recommended = document.get("recommended") or {}
    recommended_id = recommended.get("candidate_id") if isinstance(recommended, dict) else None
    presented, reduction = reduce_candidates(
        list(document["candidates"]),
        list(document["objective_sense"]),
        method,
        max_candidates,
        recommended_id,
    )
    copy = dict(document)
    copy["candidates"] = presented
    copy["reduction"] = reduction
    copy["digest"] = "sha256:" + "0" * 64
    copy["digest"] = _core.digest_json(json.dumps(copy), "digest")
    return copy


def extreme_ids(candidates: list[dict[str, Any]], senses: list[str]) -> list[str]:
    """One id per axis: the best numeric value on that axis. Ties take the smaller id."""
    found: list[str] = []
    for axis, sense in enumerate(senses):
        best_id: str | None = None
        best_value: float | None = None
        for candidate in candidates:
            vector = candidate["objective_vector"]
            if axis >= len(vector) or vector[axis] is None:
                continue
            value = _minimized(float(vector[axis]), sense)
            cid = candidate["candidate_id"]
            if best_value is None or value < best_value or (value == best_value and cid < best_id):
                best_value = value
                best_id = cid
        if best_id is not None:
            found.append(best_id)
    return found


def _by_directions(
    candidates: list[dict[str, Any]],
    senses: list[str],
    max_candidates: int,
    protected: set[str],
) -> set[str]:
    usable = [c for c in candidates if _complete(c, senses)]
    chosen = set(protected)
    if not usable:
        return chosen
    rows = _normalized([_row(c, senses) for c in usable])
    objectives = len(senses)
    if objectives == 1:
        return chosen
    divisions = 1
    directions: list[list[float]] = []
    while divisions <= 8:
        directions = _core.indicator_das_dennis(objectives, divisions)
        if len(directions) >= max_candidates:
            break
        divisions += 1
    buckets: dict[int, list[int]] = {}
    for index, row in enumerate(rows):
        nearest = min(range(len(directions)), key=lambda d: _perpendicular(row, directions[d]))
        buckets.setdefault(nearest, []).append(index)
    for direction in sorted(buckets):
        if len(chosen) >= max_candidates:
            break
        members = buckets[direction]
        winner = min(members, key=lambda i: (_magnitude(rows[i]), usable[i]["candidate_id"]))
        chosen.add(usable[winner]["candidate_id"])
    return chosen


def _knee(
    candidates: list[dict[str, Any]],
    senses: list[str],
    max_candidates: int,
    protected: set[str],
) -> set[str]:
    usable = [c for c in candidates if _complete(c, senses)]
    chosen = set(protected)
    if not usable:
        return chosen
    rows = _normalized([_row(c, senses) for c in usable])
    scores: list[tuple[float, str]] = []
    if len(senses) == 2:
        ends = [i for i, c in enumerate(usable) if c["candidate_id"] in protected]
        if len(ends) >= 2:
            a, b = rows[ends[0]], rows[ends[1]]
            for index, row in enumerate(rows):
                scores.append((_segment_distance(row, a, b), usable[index]["candidate_id"]))
    if not scores:
        for index, row in enumerate(rows):
            # Product of the remaining gap to the nadir corner (1, 1, ...).
            # Larger means closer to the ideal corner on every axis at once.
            gap = 1.0
            for value in row:
                gap *= max(0.0, 1.0 - value)
            scores.append((gap, usable[index]["candidate_id"]))
    scores.sort(key=lambda item: (-item[0], item[1]))
    for _score, cid in scores:
        if len(chosen) >= max_candidates:
            break
        chosen.add(cid)
    return chosen


def _cluster(
    candidates: list[dict[str, Any]],
    senses: list[str],
    max_candidates: int,
    protected: set[str],
) -> set[str]:
    usable = [c for c in candidates if _complete(c, senses)]
    chosen = set(protected)
    if not usable:
        return chosen
    rows = {
        c["candidate_id"]: row
        for c, row in zip(usable, _normalized([_row(c, senses) for c in usable]), strict=True)
    }
    while len(chosen) < max_candidates:
        remaining = [c["candidate_id"] for c in usable if c["candidate_id"] not in chosen]
        if not remaining:
            break
        anchors = [cid for cid in chosen if cid in rows]
        # Protected extremes can be incomplete (an UNKNOWN axis). They stay in
        # the set, but they are not points the distance can be measured from.
        if not chosen or not anchors:
            chosen.add(min(remaining))
            continue
        nxt = max(
            remaining,
            key=lambda cid: (
                min(_euclidean(rows[cid], rows[other]) for other in anchors),
                cid,
            ),
        )
        chosen.add(nxt)
    return chosen


def _epsilon_grid(
    candidates: list[dict[str, Any]],
    senses: list[str],
    max_candidates: int,
    protected: set[str],
) -> set[str]:
    usable = [c for c in candidates if _complete(c, senses)]
    chosen = set(protected)
    if not usable:
        return chosen
    rows = _normalized([_row(c, senses) for c in usable])
    cells = max(1, max_candidates)
    buckets: dict[tuple[int, ...], list[int]] = {}
    for index, row in enumerate(rows):
        cell = tuple(min(cells - 1, int(value * cells)) for value in row)
        buckets.setdefault(cell, []).append(index)
    for cell in sorted(buckets):
        if len(chosen) >= max_candidates:
            break
        members = buckets[cell]
        winner = min(members, key=lambda i: (_magnitude(rows[i]), usable[i]["candidate_id"]))
        chosen.add(usable[winner]["candidate_id"])
    return chosen


def _order(
    candidates: list[dict[str, Any]],
    chosen: set[str],
    recommended_id: str | None,
) -> list[dict[str, Any]]:
    by_id = {c["candidate_id"]: c for c in candidates}
    ordered: list[dict[str, Any]] = []
    if recommended_id in chosen and recommended_id in by_id:
        ordered.append(by_id[recommended_id])
    for candidate in candidates:
        cid = candidate["candidate_id"]
        if cid in chosen and cid != recommended_id:
            ordered.append(candidate)
    return ordered


def _record(raw: int, presented: int, method: str, max_candidates: int) -> dict[str, Any]:
    return {
        "raw_nondominated_count": raw,
        "presented_count": presented,
        "method": method,
        "parameters": {"max_candidates": max_candidates},
    }


def _ids(values: list[str]) -> set[str]:
    return set(values)


def _minimized(value: float, sense: str) -> float:
    if sense == "minimize":
        return value
    if sense == "maximize":
        return -value
    raise ValueError(f"objective sense {sense!r} is not minimize or maximize")


def _complete(candidate: dict[str, Any], senses: list[str]) -> bool:
    vector = candidate["objective_vector"]
    return len(vector) == len(senses) and all(value is not None for value in vector)


def _row(candidate: dict[str, Any], senses: list[str]) -> list[float]:
    return [
        _minimized(float(value), sense)
        for value, sense in zip(candidate["objective_vector"], senses, strict=True)
    ]


def _normalized(rows: list[list[float]]) -> list[list[float]]:
    width = len(rows[0])
    lows = [min(row[axis] for row in rows) for axis in range(width)]
    highs = [max(row[axis] for row in rows) for axis in range(width)]
    out = []
    for row in rows:
        scaled = []
        for axis, value in enumerate(row):
            span = highs[axis] - lows[axis]
            scaled.append(0.0 if span == 0.0 else (value - lows[axis]) / span)
        out.append(scaled)
    return out


def _perpendicular(point: list[float], direction: list[float]) -> float:
    denom = sum(component * component for component in direction) or 1.0
    scale = sum(p * d for p, d in zip(point, direction, strict=True)) / denom
    return sum((p - scale * d) ** 2 for p, d in zip(point, direction, strict=True)) ** 0.5


def _magnitude(row: list[float]) -> float:
    return sum(value * value for value in row) ** 0.5


def _euclidean(left: list[float], right: list[float]) -> float:
    return sum((a - b) ** 2 for a, b in zip(left, right, strict=True)) ** 0.5


def _segment_distance(point: list[float], start: list[float], end: list[float]) -> float:
    span = [b - a for a, b in zip(start, end, strict=True)]
    length = sum(component * component for component in span)
    if length == 0.0:
        return _euclidean(point, start)
    scale = (
        sum((p - a) * component for p, a, component in zip(point, start, span, strict=True))
        / length
    )
    scale = min(1.0, max(0.0, scale))
    projection = [a + scale * component for a, component in zip(start, span, strict=True)]
    return _euclidean(point, projection)
