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


# A local fold costs about the same whatever the CDS length is, because the
# refolded region is the edit plus one span on each side. A full fold costs
# about 0.06 ms per nucleotide on the machine this was measured on (375 nt
# ~21 ms, 936 nt ~58 ms, 2463 nt ~157 ms; one codon edit ~8 ms). Past this
# many codon edits the full fold is the cheaper exact answer.
_INCREMENTAL_DIVISOR = 160


def _codon_edits(before: str, after: str) -> list[tuple[int, str]] | None:
    """Codon-aligned substitutions that turn `before` into `after`.

    `None` when the strings are not the same length, so the caller folds
    from scratch instead of inventing an alignment.
    """
    if len(before) != len(after) or len(before) % 3 != 0:
        return None
    edits = []
    for start in range(0, len(before), 3):
        codon = after[start : start + 3]
        if before[start : start + 3] != codon:
            edits.append((start, codon))
    return edits


class FoldSession:
    """Exact accessibility values for one job.

    An identical CDS is returned as stored. A sequence a few codons away from
    one already folded is edited through `AccessibilityCache.substitute`,
    which recomputes only the windows inside the folding span and reuses the
    rest. That reuse is exact: a window's value depends only on bases within
    `span_nt`, and `TST-STRU-002` checks the edited values against a full
    fold. Anything further away is folded from scratch. The reported penalty
    is the same number a fresh fold would return.
    """

    def __init__(self, temperature_c: float = TEMPERATURE_C) -> None:
        self.temperature_c = temperature_c
        self.hits = 0
        self.full_folds = 0
        self.local_edits = 0
        self._cache: AccessibilityCache | None = None
        self._penalty: dict[str, float] = {}
        self._windows: dict[str, list[float]] = {}
        self._recent: list[str] = []

    def evaluate(self, sequences: list[str]) -> list[float]:
        if not sequences:
            return []
        return [self._one(sequence) for sequence in sequences]

    def _cache_for(self) -> AccessibilityCache:
        if self._cache is None:
            self._cache = AccessibilityCache(
                engine=select_engine(),
                window_nt=WINDOW_NT,
                span_nt=SPAN_NT,
                temperature_c=self.temperature_c,
            )
        return self._cache

    def _closest(self, sequence: str) -> tuple[str, list[tuple[int, str]]] | None:
        limit = max(1, len(sequence) // _INCREMENTAL_DIVISOR)
        best: tuple[int, str, list[tuple[int, str]]] | None = None
        for parent in reversed(self._recent[-48:]):
            if len(parent) != len(sequence):
                continue
            edits = _codon_edits(parent, sequence)
            if edits is None:
                continue
            count = len(edits)
            if count == 0 or count > limit:
                continue
            if best is None or count < best[0]:
                best = (count, parent, edits)
                if count == 1:
                    break
        if best is None:
            return None
        return best[1], best[2]

    def _remember(self, sequence: str, values: list[float], penalty: float) -> float:
        self._penalty[sequence] = penalty
        self._windows[sequence] = list(values)
        self._recent.append(sequence)
        return penalty

    def _one(self, sequence: str) -> float:
        cached = self._penalty.get(sequence)
        if cached is not None:
            self.hits += 1
            return cached
        cache = self._cache_for()
        parent = self._closest(sequence)
        if parent is None:
            values = cache.prime(sequence)
            self.full_folds += 1
            return self._remember(sequence, values, accessibility_penalty(values))

        parent_sequence, edits = parent
        if cache.current_sequence() != parent_sequence:
            cache.hold(parent_sequence, self._windows[parent_sequence])
        values: list[float] = []
        for start, codon in edits:
            values = list(cache.substitute(start, codon))
            self.local_edits += 1
        if cache.current_sequence() != sequence:
            # The edit list did not reconstruct the sequence. Fold it rather
            # than store a penalty for a different string.
            values = cache.prime(sequence)
            self.full_folds += 1
        return self._remember(sequence, values, accessibility_penalty(values))


def evaluate(
    sequences: list[str],
    temperature_c: float = TEMPERATURE_C,
    session: FoldSession | None = None,
) -> list[float]:
    """Accessibility penalty per sequence, from a real fold.

    Repeated sequences and sequences a few codons apart share one
    `FoldSession`. The shared value is the full-fold value: identical
    strings are stored, and a short edit goes through `AccessibilityCache`,
    which is exact for the windows outside the edit (`TST-STRU-002`). A
    caller that passes `session` keeps that memory across calls, which is
    how one generation reuses the previous generation's folds.
    """
    if not sequences:
        return []
    if session is None or session.temperature_c != temperature_c:
        session = FoldSession(temperature_c)
    return session.evaluate(sequences)


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
