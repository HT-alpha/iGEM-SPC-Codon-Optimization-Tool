"""Evaluate `local_structure` rules with the installed folding engine.

The Rust checker still lists the rule as undetectable. That listing is what
keeps an in-core search from treating the rule as satisfied. This module is
the evaluation: when ViennaRNA is installed the orchestrator measures the
declared span and drops the `CF-FEAS-008` warning for those ids. A missing
span is `CF-STRU-004`. A missing engine is `CF-STRU-003`. `min_codon_frequency`
stays unsearched.

`DEP-CONS-006`.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

_METRICS = ("mean_unpaired_probability", "min_free_energy")
_REPAIR_CODONS = 8


def rules(constraints: dict[str, Any] | None) -> list[dict[str, Any]]:
    """Enabled `local_structure` rules, in document order."""
    found = []
    for item in (constraints or {}).get("constraints") or []:
        if item.get("enabled", True) is False:
            continue
        rule = item.get("rule") or {}
        if rule.get("kind") != "local_structure":
            continue
        found.append(
            {
                "id": str(item.get("id") or "local_structure"),
                "blocking": item.get("class", "blocking") == "blocking",
                "rule": rule,
            }
        )
    return found


def _refusal(code: str, detail: str) -> Exception:
    from dna_design_suite.codon.orchestrator.design import DesignError

    return DesignError(code, detail, exit_code=2)


def require_engine(
    constraints: dict[str, Any] | None,
    structure_requested: bool,
    structure_engine: str = "",
) -> None:
    """Refuse before a job when a fold was asked for and no engine is installed.

    `nussinov` may score the structure objective. It does not supply
    `mean_unpaired_probability` or `min_free_energy`, so a `local_structure`
    rule still needs ViennaRNA.
    """
    from dna_design_suite.codon.orchestrator import structure_diag

    if structure_engine == "nussinov" and rules(constraints):
        raise _refusal(
            "CF-STRU-003",
            "local_structure metrics are ViennaRNA measurements "
            "(mean unpaired probability or minimum free energy). "
            "The Nussinov score is a paired fraction, not those metrics.",
        )
    if not rules(constraints) and not structure_requested:
        return
    if structure_engine == "nussinov" and structure_requested:
        return
    if structure_diag.is_available():
        return
    raise _refusal(
        "CF-STRU-003",
        "No structure engine is installed and DNA Design Suite ships none",
    )


def validate(constraints: dict[str, Any] | None) -> None:
    """A span that is missing is not replaced with a default window."""
    for item in rules(constraints):
        rule = item["rule"]
        span = rule.get("span")
        if not isinstance(span, dict) or "start_nt" not in span or "end_nt" not in span:
            raise _refusal(
                "CF-STRU-004",
                f"{item['id']}: local_structure needs span.start_nt and span.end_nt; "
                "a default window would be an invented region",
            )
        metric = rule.get("metric")
        if metric not in _METRICS:
            raise _refusal(
                "CF-STRU-004",
                f"{item['id']}: local_structure metric {metric!r} is not "
                "mean_unpaired_probability or min_free_energy",
            )


def evaluated_ids(constraints: dict[str, Any] | None) -> set[str]:
    return {item["id"] for item in rules(constraints)}


def measure(
    construct: str,
    rule: dict[str, Any],
    cache: dict[tuple[Any, ...], float],
) -> float:
    """The metric on the declared span.

    `mean_unpaired_probability` is the mean of window length 1 accessibility
    on that span, folded in the construct so the flanks are context.
    `min_free_energy` is the fold of the declared span alone.
    """
    span = rule["span"]
    start = int(span["start_nt"])
    end = int(span["end_nt"])
    if start < 0 or end > len(construct) or start >= end:
        raise _refusal(
            "CF-STRU-004",
            f"local_structure span {start}:{end} does not fit a construct of "
            f"length {len(construct)}",
        )
    metric = rule["metric"]
    key = (metric, start, end, construct)
    cached = cache.get(key)
    if cached is not None:
        return cached
    if metric == "mean_unpaired_probability":
        value = _mean_unpaired(construct, start, end)
    else:
        value = _span_mfe(construct[start:end])
    cache[key] = value
    return value


def holds(value: float, rule: dict[str, Any]) -> bool:
    lower = rule.get("min")
    upper = rule.get("max")
    below = lower is not None and value < float(lower)
    above = upper is not None and value > float(upper)
    return not below and not above


def blocking_distance(construct: str, constraints: dict[str, Any] | None, cache: dict) -> float:
    total = 0.0
    for item in rules(constraints):
        if not item["blocking"]:
            continue
        value = measure(construct, item["rule"], cache)
        rule = item["rule"]
        lower = rule.get("min")
        upper = rule.get("max")
        if lower is not None:
            total += max(0.0, float(lower) - value)
        if upper is not None:
            total += max(0.0, value - float(upper))
    return total


def repair(
    cds: str,
    constraints: dict[str, Any] | None,
    flank5: str,
    flank3: str,
    translation_table: int,
    still_feasible: Callable[[str], bool],
    cache: dict[tuple[Any, ...], float],
) -> str:
    """Try a few synonymous edits inside the violated spans.

    An edit is kept only when it lowers the structure distance and the rest of
    the constraint set still accepts the sequence. Eight codon positions, one
    pass. A span that still fails stays a failed constraint, not a warning.
    """
    from dna_design_suite import _core

    if not rules(constraints):
        return cds
    construct = flank5 + cds + flank3
    if blocking_distance(construct, constraints, cache) <= 0.0:
        return cds
    positions = _positions(cds, constraints, flank5, flank3, cache)
    current = cds
    for index in positions[:_REPAIR_CODONS]:
        codon = current[index * 3 : (index + 1) * 3]
        best = current
        best_distance = blocking_distance(flank5 + current + flank3, constraints, cache)
        for synonym in _core.synonyms(codon, translation_table):
            if synonym == codon:
                continue
            trial = current[: index * 3] + synonym + current[(index + 1) * 3 :]
            if not still_feasible(trial):
                continue
            distance = blocking_distance(flank5 + trial + flank3, constraints, cache)
            if distance < best_distance:
                best = trial
                best_distance = distance
        current = best
        if best_distance <= 0.0:
            break
    return current


def _positions(
    cds: str,
    constraints: dict[str, Any] | None,
    flank5: str,
    flank3: str,
    cache: dict,
) -> list[int]:
    construct = flank5 + cds + flank3
    cds_lo = len(flank5)
    cds_hi = cds_lo + len(cds)
    chosen: list[int] = []
    for item in rules(constraints):
        if not item["blocking"]:
            continue
        if holds(measure(construct, item["rule"], cache), item["rule"]):
            continue
        span = item["rule"]["span"]
        start = max(int(span["start_nt"]), cds_lo)
        end = min(int(span["end_nt"]), cds_hi)
        if start >= end:
            continue
        first = (start - cds_lo) // 3
        last = (end - cds_lo + 2) // 3
        for index in range(first, min(last, len(cds) // 3)):
            if index not in chosen:
                chosen.append(index)
    chosen.sort()
    return chosen


def _mean_unpaired(construct: str, start: int, end: int) -> float:
    from dna_design_suite.codon.predictors.viennarna import ViennaRNAEngine

    values = ViennaRNAEngine().local_accessibility(
        construct, window_nt=1, span_nt=1, temperature_c=37.0
    )
    window = list(values)[start:end]
    if len(window) != end - start:
        raise _refusal(
            "CF-STRU-004",
            f"accessibility returned {len(window)} values for span {start}:{end}",
        )
    return sum(window) / len(window)


def _span_mfe(span: str) -> float:
    import RNA

    from dna_design_suite.codon.predictors.viennarna import FOLD_LOCK

    with FOLD_LOCK:
        RNA.cvar.temperature = 37.0
        _structure, energy = RNA.fold(span.upper().replace("T", "U"))
        return float(energy)
