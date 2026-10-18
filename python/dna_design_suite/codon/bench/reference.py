"""The cross-run hypervolume reference point. `P2-19`, `DEP-BEN-001`.

`cf-search`'s per-generation hypervolume uses the population nadir plus one.
That number is a stagnation signal inside one run. Two runs with different
fronts get different nadirs, so comparing those hypervolumes does not compare
the fronts. This function takes only the published ranges. It has no front
argument, so a caller cannot pass the candidates in by accident.
"""

from __future__ import annotations

from dna_design_suite import _core


class BenchError(Exception):
    """A benchmark step that would produce a number nobody could compare."""

    def __init__(self, detail: str) -> None:
        self.detail = detail
        super().__init__(detail)


def fixed_reference_point(
    ranges: list[tuple[float, float]],
    senses: list[str],
) -> list[float]:
    """Nadir end of each published range, in native objective units.

    `minimize` uses the high end. `maximize` uses the low end. A reversed
    range (`lo > hi`) or a missing sense is refused: filling it from the
    front would make the point candidate-dependent again.
    """
    if len(ranges) != len(senses):
        raise BenchError(
            f"{len(ranges)} ranges for {len(senses)} objectives; a reference "
            "point has one coordinate per objective"
        )
    point: list[float] = []
    for (low, high), sense in zip(ranges, senses, strict=True):
        if low > high:
            raise BenchError(
                f"reference range {(low, high)} is reversed; the published interval is [low, high]"
            )
        if sense == "minimize":
            point.append(float(high))
        elif sense == "maximize":
            point.append(float(low))
        else:
            raise BenchError(f"objective sense {sense!r} is not minimize or maximize")
    return point


def hypervolume(
    points: list[list[float]],
    reference: list[float],
    senses: list[str],
) -> float:
    """Hypervolume in minimization space, via `cf-metrics`.

    Maximize axes are negated, point and reference together, so the Rust
    indicator's "all minimization" contract holds without a second formula.
    """
    if any(len(point) != len(reference) or len(point) != len(senses) for point in points):
        raise BenchError("a point, the reference and the senses do not share a width")
    minimized = [
        [_minimized(value, sense) for value, sense in zip(point, senses, strict=True)]
        for point in points
    ]
    minimized_reference = [
        _minimized(value, sense) for value, sense in zip(reference, senses, strict=True)
    ]
    return float(_core.indicator_hypervolume(minimized, minimized_reference))


def _minimized(value: float, sense: str) -> float:
    if sense == "minimize":
        return float(value)
    if sense == "maximize":
        return -float(value)
    raise BenchError(f"objective sense {sense!r} is not minimize or maximize")
