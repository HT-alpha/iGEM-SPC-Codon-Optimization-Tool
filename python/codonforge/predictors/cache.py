"""Incremental accessibility cache.

Ticket P2-10. `TST-STRU-002`, `TST-STRU-006`.

A folding call is the most expensive thing in a design job, and a synonymous
substitution changes three nucleotides. Refolding the whole construct for a
three-nucleotide edit is the waste this exists to remove.

What "incremental" can and cannot mean here
-------------------------------------------
`AMD-013` put the engine out of process and `DEP-STRU-002` made the seam
stateless, so this cache cannot reach inside the engine's partition function. It
works at the level the seam exposes: it re-folds a **bounded region** around the
edit and reuses cached values everywhere else.

That is sound only because accessibility is *local*. A base pair can only form
within `span_nt`, so an edit at position `p` can only change the accessibility of
windows whose folding neighbourhood contains `p` — windows ending within
`span_nt` of the edit. Outside that radius the value is unchanged, and reusing it
is exact rather than approximate.

The radius is the load-bearing claim. If it were too small the cache would serve
values that are stale in a way no agreement test on the happy path would notice,
which is exactly the failure mode the sixth fixture-design rule was written for.

What this deliberately does not do
----------------------------------
No eviction policy and no size cap. A design job folds one construct of known
length, so the cache is bounded by the sequence, not by the run. Adding an LRU
would be a second mechanism deciding what is valid, alongside the invalidation
radius, and the two would drift.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from codonforge.predictors.structure import StructureEngine, StructureEngineError


@dataclass
class CacheStats:
    """Counters, exposed because a cache that never caches is a real defect.

    `TST-STRU-002` — agreement with full recomputation — is satisfied perfectly
    by a cache that recomputes everything every time. These counters are what
    make the difference observable, so "it is a cache" is asserted rather than
    assumed.
    """

    folds: int = 0
    """Engine calls actually made."""

    windows_recomputed: int = 0
    """Window values obtained from the engine."""

    windows_reused: int = 0
    """Window values served from the previous state."""

    def reset(self) -> None:
        self.folds = 0
        self.windows_recomputed = 0
        self.windows_reused = 0


@dataclass
class AccessibilityCache:
    """Windowed accessibility for a sequence under edits.

    Holds one sequence at a time. `substitute` reports the new accessibility
    after an edit, folding only the affected region.
    """

    engine: StructureEngine
    window_nt: int
    span_nt: int
    temperature_c: float

    _sequence: str | None = None
    _values: list[float] = field(default_factory=list)
    stats: CacheStats = field(default_factory=CacheStats)

    def prime(self, sequence: str) -> list[float]:
        """Fold a sequence from scratch and hold it as the current state."""
        values = list(
            self.engine.local_accessibility(
                sequence, self.window_nt, self.span_nt, self.temperature_c
            )
        )
        self.stats.folds += 1
        self.stats.windows_recomputed += len(values)
        self._sequence = sequence
        self._values = values
        return list(values)

    def current(self) -> list[float]:
        if self._sequence is None:
            raise StructureEngineError(
                "CF-STRU-005", "the cache has not been primed with a sequence"
            )
        return list(self._values)

    def affected_window_range(self, start: int, stop: int) -> tuple[int, int]:
        """Which window indices an edit over `[start, stop)` can change.

        A window ending at sequence position `e` covers `[e - window + 1, e]` and
        folds against everything within `span_nt` of it. So an edit at `p`
        reaches windows whose end lies within `span_nt` of `p`, and nothing
        further.

        Returned as a half-open range of *window indices*, where index `k` is the
        window ending at sequence position `k + window_nt - 1` — the same
        convention the engine returns.
        """
        if self._sequence is None:
            raise StructureEngineError(
                "CF-STRU-005", "the cache has not been primed with a sequence"
            )
        first_end = max(self.window_nt - 1, start - self.span_nt)
        last_end = min(len(self._sequence) - 1, stop - 1 + self.span_nt)
        low = first_end - (self.window_nt - 1)
        high = last_end - (self.window_nt - 1) + 1
        return (max(0, low), max(0, min(len(self._values), high)))

    def substitute(self, start: int, replacement: str) -> list[float]:
        """Apply an edit and return the new accessibility.

        Exact, not approximate: values outside the affected range are unchanged
        by the edit, so reusing them is not a shortcut.
        """
        if self._sequence is None:
            raise StructureEngineError(
                "CF-STRU-005", "the cache has not been primed with a sequence"
            )
        stop = start + len(replacement)
        if start < 0 or stop > len(self._sequence):
            raise StructureEngineError(
                "CF-STRU-004",
                f"edit [{start}, {stop}) falls outside a sequence of {len(self._sequence)} nt",
            )

        edited = self._sequence[:start] + replacement + self._sequence[stop:]
        low, high = self.affected_window_range(start, stop)

        # Fold a region wide enough that every affected window's folding
        # neighbourhood is inside it. The region is the affected windows plus a
        # span of context on each side, because a window's value depends on what
        # it could pair with, not only on what it covers.
        region_start = max(0, low + self.window_nt - 1 - (self.window_nt - 1) - self.span_nt)
        region_stop = min(len(edited), high + self.window_nt - 1 + self.span_nt)
        region = edited[region_start:region_stop]

        if len(region) < self.window_nt:
            # Too short to fold meaningfully; fall back to a full recomputation
            # rather than reporting from a region that cannot hold a window.
            self._sequence = edited
            return self.prime(edited)

        folded = list(
            self.engine.local_accessibility(
                region, self.window_nt, self.span_nt, self.temperature_c
            )
        )
        self.stats.folds += 1

        values = list(self._values)
        recomputed = 0
        for index in range(low, high):
            end_in_sequence = index + self.window_nt - 1
            offset = end_in_sequence - region_start - (self.window_nt - 1)
            if 0 <= offset < len(folded):
                values[index] = folded[offset]
                recomputed += 1
        self.stats.windows_recomputed += recomputed
        self.stats.windows_reused += len(values) - recomputed

        self._sequence = edited
        self._values = values
        return list(values)
