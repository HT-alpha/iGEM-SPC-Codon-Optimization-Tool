"""A secondary model of a numeric readout, or an explicit refusal to fit one.

The only feature used is the GC fraction of an explicit CDS. It is counted
from that sequence. No codon weight is invented, and a design reference that
does not carry the CDS is not scored. Censored readouts are not treated as
points. Rank-deficient or too-small sets stay `not_fitted`.
"""

from __future__ import annotations

from typing import Any


def _gc_fraction(cds: str) -> float | None:
    text = cds.upper().replace("U", "T")
    if not text or any(base not in "ACGT" for base in text):
        return None
    return sum(1 for base in text if base in "GC") / len(text)


def _point_readout(observation: dict[str, Any]) -> tuple[str, float] | None:
    """The single uncensored numeric readout, when the observation has exactly one."""
    readouts = observation.get("readouts")
    if not isinstance(readouts, list):
        return None
    points: list[tuple[str, float]] = []
    for readout in readouts:
        if not isinstance(readout, dict):
            return None
        censoring = (readout.get("dynamic_range") or {}).get("censoring")
        if censoring != "none":
            continue
        value = readout.get("value")
        measures = readout.get("measures")
        numeric = isinstance(value, (int, float)) and not isinstance(value, bool)
        if not isinstance(measures, str) or not numeric:
            return None
        points.append((measures, float(value)))
    if len(points) != 1:
        return None
    return points[0]


def fit_outcome(observations: list[dict[str, Any]]) -> dict[str, Any]:
    """Ordinary least squares of one latent quantity on GC fraction.

    Two or more explicit point observations of the same `measures`, with at
    least two distinct GC fractions, produce coefficients. Every other case
    returns `state: not_fitted` and no coefficients.
    """
    groups: dict[str, list[tuple[float, float]]] = {}
    explicit = 0
    for observation in observations:
        subject = observation.get("subject") or {}
        if subject.get("kind") != "explicit":
            continue
        explicit += 1
        point = _point_readout(observation)
        fraction = _gc_fraction(str(subject.get("cds", "")))
        if point is None or fraction is None:
            continue
        measures, value = point
        groups.setdefault(measures, []).append((fraction, value))

    eligible = {name: rows for name, rows in groups.items() if len(rows) >= 2}
    if len(eligible) == 0:
        return {
            "state": "not_fitted",
            "feature": "gc_fraction",
            "n": max((len(rows) for rows in groups.values()), default=0),
            "n_explicit": explicit,
            "reason": (
                "fewer than two explicit point observations share one measures value; "
                "no coefficient is reported"
            ),
        }
    if len(eligible) > 1:
        names = ", ".join(sorted(eligible))
        return {
            "state": "not_fitted",
            "feature": "gc_fraction",
            "n": 0,
            "n_explicit": explicit,
            "reason": (
                f"more than one latent quantity has at least two point rows ({names}); "
                "they are not pooled, and no quantity was chosen for them"
            ),
        }
    measures, rows = next(iter(eligible.items()))
    rows = sorted(rows, key=lambda item: (item[0], item[1]))
    gc = [item[0] for item in rows]
    values = [item[1] for item in rows]
    if len(set(gc)) < 2:
        return {
            "state": "not_fitted",
            "feature": "gc_fraction",
            "measures": measures,
            "n": len(rows),
            "n_explicit": explicit,
            "reason": (
                "gc_fraction does not vary across these observations, so the "
                "design matrix is rank-deficient and no coefficient is reported"
            ),
        }
    n = len(rows)
    sum_g = sum(gc)
    sum_gg = sum(item * item for item in gc)
    sum_y = sum(values)
    sum_gy = sum(g * y for g, y in zip(gc, values, strict=True))
    determinant = n * sum_gg - sum_g * sum_g
    if determinant == 0.0:
        return {
            "state": "not_fitted",
            "feature": "gc_fraction",
            "measures": measures,
            "n": n,
            "n_explicit": explicit,
            "reason": "the normal-equation matrix is singular and no coefficient is reported",
        }
    intercept = (sum_gg * sum_y - sum_g * sum_gy) / determinant
    slope = (n * sum_gy - sum_g * sum_y) / determinant
    residual = 0.0
    for g, y in zip(gc, values, strict=True):
        error = y - (intercept + slope * g)
        residual += error * error
    return {
        "state": "fitted",
        "feature": "gc_fraction",
        "measures": measures,
        "n": n,
        "n_explicit": explicit,
        "coefficients": {"intercept": intercept, "gc_fraction": slope},
        "residual_sum_of_squares": residual,
    }
