"""The structure term, wired in as a diagnostic.

`DEP-STRU-003`. Interim path, and deliberately labelled as one.

Why this is not a Pareto axis yet
---------------------------------
`structure.local_accessibility` was planned as an L1 term computed in-core. That
assumed the from-scratch Rust kernel `AMD-013` removed. With the engine out of
process, a per-move exact delta would cross the Python boundary once per
candidate — about 19 200 crossings for a 64/300 run against `ADR-0001`'s budget
of 300 — so the term cannot drive the search until tier 3 can hand out batches.
That is `P2-05a`, registered rather than absorbed into this round.

What ships instead is honest and useful on its own: the term is computed from a
real ViennaRNA fold, for every candidate on the final front, and reported. A
user can see which design is more accessible. Nothing pretends the search
optimized for it.

The visible-gap rule
--------------------
`CF-ORCH-006` exists because a tier that silently does not run makes a job look
like more search than it was. A term reported next to real objectives, in the
same document, with nothing saying it took no part in the optimization, is the
same failure in a different place. So `CF-LAND-013` says so in
`job_status.warnings`, with stage and subject provenance, on every job where the
term is computed.
"""

from __future__ import annotations

from typing import Any

from codonforge.orchestrator.budget import Warning_
from codonforge.predictors import StructureEngineError, available_engines, select_engine
from codonforge.predictors.cache import AccessibilityCache
from codonforge.predictors.viennarna import accessibility_penalty

TERM_ID = "structure.local_accessibility"

#: Windowed-accessibility geometry. `AMD-001`: windowed and local, never a
#: global MFE. These are defaults, overridable per job once the term is an axis.
WINDOW_NT = 6
SPAN_NT = 40
TEMPERATURE_C = 37.0


def is_available() -> bool:
    return bool(available_engines())


def evaluate(sequences: list[str], temperature_c: float = TEMPERATURE_C) -> list[float]:
    """Accessibility penalty per sequence, from a real fold.

    One engine and one cache across the batch: front members are usually near
    neighbours, so priming once and substituting would be tempting — but they
    differ at many positions rather than one, and a cache is only exact for the
    edits it was told about. Folding each independently is the correct thing,
    and pretending otherwise would be the stale-value bug `TST-STRU-006` exists
    to catch, introduced deliberately.
    """
    if not sequences:
        return []
    engine = select_engine()
    values = []
    for sequence in sequences:
        cache = AccessibilityCache(
            engine=engine,
            window_nt=WINDOW_NT,
            span_nt=SPAN_NT,
            temperature_c=temperature_c,
        )
        values.append(accessibility_penalty(cache.prime(sequence)))
    return values


def detail(value: float) -> dict[str, dict[str, Any]]:
    """The `objective_detail` entry for one candidate."""
    return {
        TERM_ID: {
            "value": value,
            "unit": "dimensionless",
            "percentile_vs_native": None,
        }
    }


def diagnostic_warning(engine_name: str) -> dict[str, Any]:
    """`CF-LAND-013` — the term is reported but did not drive the search."""
    return Warning_(
        code="CF-LAND-013",
        title="A term was computed for reporting only, not optimized against",
        detail=(
            f"{TERM_ID} was evaluated with {engine_name} for every candidate on the "
            "front and appears in objective_detail, but it is not one of the "
            "objectives the search ranked on and is absent from objective_vector. "
            "Tier 3 evaluates one candidate at a time, so an out-of-process engine "
            "cannot be a search objective until the generation loop hands out "
            "batches (P2-05a). Designs on this front were not selected for "
            "accessibility"
        ),
        stage="landscape",
        subject=TERM_ID,
    ).to_problem_details()


def unavailable_warning() -> dict[str, Any]:
    """No engine installed. Said out loud rather than omitted silently."""
    return Warning_(
        code="CF-LAND-013",
        title="A diagnostic term could not be computed",
        detail=(
            f"{TERM_ID} was not computed because no structure engine is installed. "
            "CodonForge ships none by design (AMD-013): install ViennaRNA to see "
            "accessibility reported alongside the front. The search is unaffected, "
            "since this term does not drive it"
        ),
        stage="landscape",
        subject=TERM_ID,
    ).to_problem_details()


def annotate(sequences: list[str]) -> tuple[list[dict[str, Any]] | None, dict[str, Any]]:
    """Details for each sequence, plus the warning that must accompany them.

    Returns `(None, warning)` when no engine is installed, so a caller cannot
    get details without also getting the warning that explains what they are.
    Pairing them in one return value is the cheapest way to stop the two
    drifting apart.
    """
    if not is_available():
        return None, unavailable_warning()
    try:
        engine_name = select_engine().identity().name
        values = evaluate(sequences)
    except StructureEngineError:
        return None, unavailable_warning()
    return [detail(value) for value in values], diagnostic_warning(engine_name)
