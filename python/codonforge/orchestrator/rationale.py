"""Per-position rationale stored on the recommended candidate.

The numbers come from the same term scorer the front already uses. This module
does not recompute CAI by a second formula: a zero weight is floored inside a
copy of the host profile and `term_scores` is called again, and the Kazusa
table is scored the same way.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from codonforge import _core

_KAZUSA = Path(__file__).resolve().parents[3] / "examples" / "kazusa_e_coli_316407.json"
_KAZUSA_NAME = (
    "Kazusa E. coli genome-wide frequencies (taxid 316407; Nakamura, Gojobori "
    "and Ikemura, NAR 2000) converted to relative adaptiveness by dividing each "
    "codon by the maximum frequency in its family"
)


def build_rationale(
    request: Any,
    cds: str,
    tier_request: str,
    scores: dict[str, Any],
) -> dict[str, Any] | None:
    """Annotations for one candidate. Empty when CAI was not scored and nothing moved."""
    annotations: list[dict[str, Any]] = []
    codon_count = len(cds) // 3
    if codon_count < 1:
        return None
    span = {"start": 0, "end": codon_count}
    cai = scores.get("usage.cai") if "usage.cai" in scores else None
    if isinstance(cai, (int, float)):
        annotations.append(_cai_definition(request, float(cai), span))
        annotations.extend(_zero_weight(request, cds, tier_request, span))
        reference = _kazusa(request, cds, tier_request, span)
        if reference is not None:
            annotations.append(reference)
    annotations.extend(_deviations(tier_request, cds))
    if not annotations:
        return None
    return {"candidate_tags": [], "annotations": annotations}


def _cai_definition(request: Any, cai: float, span: dict[str, int]) -> dict[str, Any]:
    usage = request.host_profile.get("codon_usage", {})
    subset = usage.get("subset", "loaded")
    text = (
        f"CAI {cai:.6f} is the geometric mean of relative adaptiveness on the loaded "
        f"{subset} table for {request.host}. 1.0 means every counted codon is a top "
        "synonym of that table. Stops and single-codon families are excluded."
    )
    return _annotation("cai_table_relative", span, text)


def _zero_weight(
    request: Any, cds: str, tier_request: str, span: dict[str, int]
) -> list[dict[str, Any]]:
    table = request.host_profile.get("codon_usage", {}).get("table", {})
    weights = {
        str(codon).upper(): float(entry["relative_adaptiveness"])
        for codon, entry in table.items()
        if isinstance(entry, dict) and "relative_adaptiveness" in entry
    }
    positive = [weight for weight in weights.values() if weight > 0]
    if not positive:
        return []
    floor = min(positive)
    counted = _counted_cache(request.translation_table)
    zeros: list[tuple[int, str]] = []
    for index in range(0, len(cds) - len(cds) % 3, 3):
        codon = cds[index : index + 3]
        if not counted.get(codon, False):
            continue
        if weights.get(codon, 1.0) == 0.0:
            zeros.append((index // 3, codon))
    if not zeros:
        return []
    sensitivity = _rescore(tier_request, cds, _floored(request.host_profile, floor))
    shown = f"{sensitivity:.6f}" if sensitivity is not None else "unavailable"
    notes = []
    for position, codon in zeros:
        text = (
            f"Codon {codon} at position {position} has relative adaptiveness 0 on the "
            f"loaded table, so CAI is 0. Replacing every 0 in that table with the "
            f"smallest positive weight {floor:.6g} scores {shown}."
        )
        notes.append(
            _annotation("cai_zero_weight", {"start": position, "end": position + 1}, text)
        )
    return notes


def _kazusa(
    request: Any, cds: str, tier_request: str, span: dict[str, int]
) -> dict[str, Any] | None:
    if not _KAZUSA.is_file():
        return None
    profile = json.loads(_KAZUSA.read_text(encoding="utf-8"))
    loaded = request.host_profile.get("codon_usage", {}).get("table")
    if loaded == profile.get("codon_usage", {}).get("table"):
        return None
    value = _rescore(tier_request, cds, profile)
    if value is None:
        return None
    text = (
        f"The same CDS scores {value:.6f} on {_KAZUSA_NAME}. "
        "That table did not choose these codons."
    )
    return _annotation("cai_reference_table", span, text)


def _deviations(tier_request: str, cds: str) -> list[dict[str, Any]]:
    raw = json.loads(_core.explain_baseline(tier_request, cds))
    notes = []
    for item in raw:
        position = int(item["codon_index"])
        blocked = item.get("constraints_if_reverted") or []
        if blocked:
            reason = (
                "putting this codon back to "
                f"{item['baseline']} would raise {', '.join(blocked)}"
            )
        else:
            reason = (
                "putting this one codon back to "
                f"{item['baseline']} does not add a violation the checker reports; "
                "the choice can still be holding a constraint that only appears "
                "together with another position"
            )
        text = (
            f"Position {position} is {item['codon']} (weight {item['weight']:.6g}) "
            f"rather than the most frequent synonym {item['baseline']} "
            f"(weight {item['baseline_weight']:.6g}) because {reason}."
        )
        notes.append(
            _annotation(
                "baseline_deviation",
                {"start": position, "end": position + 1},
                text,
            )
        )
    return notes


def _floored(profile: dict[str, Any], floor: float) -> dict[str, Any]:
    usage = profile.get("codon_usage", {})
    table = {}
    for codon, entry in usage.get("table", {}).items():
        if not isinstance(entry, dict):
            table[codon] = entry
            continue
        copied = dict(entry)
        if float(copied.get("relative_adaptiveness", 1.0)) == 0.0:
            copied["relative_adaptiveness"] = floor
        table[codon] = copied
    return {"codon_usage": {**usage, "table": table}}


def _rescore(tier_request: str, cds: str, profile: dict[str, Any]) -> float | None:
    payload = json.loads(tier_request)
    payload["host_profile"] = profile
    payload["terms"] = ["usage.cai"]
    payload["weights"] = [1.0]
    try:
        scores = json.loads(_core.term_scores(json.dumps(payload), cds))
    except ValueError:
        return None
    value = scores.get("usage.cai")
    if not isinstance(value, (int, float)):
        return None
    return float(value)


def _counted_cache(translation_table: int) -> dict[str, bool]:
    cache: dict[str, bool] = {}
    bases = "ACGT"
    for a in bases:
        for b in bases:
            for c in bases:
                codon = a + b + c
                protein = _core.translate(codon, translation_table)
                family = _core.synonyms(codon, translation_table)
                cache[codon] = protein not in ("", "*") and len(family) > 1
    return cache


def _annotation(tag: str, span: dict[str, int], text: str) -> dict[str, Any]:
    return {
        "tag": tag,
        "span": span,
        "term_id": "usage.cai",
        "magnitude": None,
        "confidence": None,
        "text": text,
    }
