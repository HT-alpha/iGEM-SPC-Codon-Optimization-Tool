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

from dna_design_suite.codon.predictors.structure import StructureEngine, StructureEngineError


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

    def current_sequence(self) -> str | None:
        return self._sequence

    def hold(self, sequence: str, values: list[float]) -> None:
        """Restore a sequence already folded, without calling the engine.

        `substitute` edits whatever this cache currently holds. A batch that
        branches from an earlier sequence has to put that sequence back first.
        """
        expected = len(sequence) - self.window_nt + 1
        if len(values) != expected:
            raise StructureEngineError(
                "CF-STRU-004",
                f"holding {len(values)} window values for a sequence of {len(sequence)} nt; "
                f"window {self.window_nt} needs {expected}",
            )
        self._sequence = sequence
        self._values = list(values)

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

    def _region_bounds(self, length: int, low: int, high: int) -> tuple[int, int]:
        """Nucleotide span whose fold covers windows `[low, high)`.

        The same span `substitute` folds for one edit. Windows outside it do
        not depend on bases inside the edit.
        """
        region_start = max(0, low - self.span_nt)
        region_stop = min(length, high + self.window_nt - 1 + self.span_nt)
        return region_start, region_stop

    def apply_edits(self, edits: list[tuple[int, str]]) -> tuple[list[float], str]:
        """Apply several same-length edits and fold each overlapping cluster once.

        A single edit goes through `substitute`, so that path is unchanged.
        Edits whose fold regions overlap or abut are one region fold: folding
        them separately would recompute the shared windows. Edits whose
        regions do not meet stay separate folds. One fold of the span between
        distant edits refolds bases that neither edit can change. When a
        cluster covers the whole sequence, the result is one full fold.
        The returned penalty matches a fold of the edited sequence.
        """
        if self._sequence is None:
            raise StructureEngineError(
                "CF-STRU-005", "the cache has not been primed with a sequence"
            )
        if len(edits) == 1:
            start, replacement = edits[0]
            return self.substitute(start, replacement), "substitute"
        sequence = self._sequence
        edited = sequence
        for start, replacement in sorted(edits, key=lambda item: item[0], reverse=True):
            stop = start + len(replacement)
            if start < 0 or stop > len(edited):
                raise StructureEngineError(
                    "CF-STRU-004",
                    f"edit [{start}, {stop}) falls outside a sequence of {len(edited)} nt",
                )
            edited = edited[:start] + replacement + edited[stop:]
        if len(edited) != len(sequence):
            return self.prime(edited), "prime"
        spans: list[tuple[int, int, int, int]] = []
        for start, replacement in edits:
            low, high = self.affected_window_range(start, start + len(replacement))
            region_start, region_stop = self._region_bounds(len(edited), low, high)
            spans.append((region_start, region_stop, low, high))
        spans.sort()
        clusters: list[tuple[int, int, int, int]] = []
        for region_start, region_stop, low, high in spans:
            if clusters and region_start <= clusters[-1][1]:
                prev_start, prev_stop, prev_low, prev_high = clusters[-1]
                clusters[-1] = (
                    prev_start,
                    max(prev_stop, region_stop),
                    min(prev_low, low),
                    max(prev_high, high),
                )
            else:
                clusters.append((region_start, region_stop, low, high))
        if any(
            stop - start >= len(edited) or stop <= start or stop - start < self.window_nt
            for start, stop, _, _ in clusters
        ):
            return self.prime(edited), "prime"
        values = list(self._values)
        recomputed_total = 0
        for region_start, region_stop, low, high in clusters:
            region = edited[region_start:region_stop]
            folded = list(
                self.engine.local_accessibility(
                    region, self.window_nt, self.span_nt, self.temperature_c
                )
            )
            self.stats.folds += 1
            for index in range(low, high):
                end_in_sequence = index + self.window_nt - 1
                offset = end_in_sequence - region_start - (self.window_nt - 1)
                if 0 <= offset < len(folded):
                    values[index] = folded[offset]
                    recomputed_total += 1
        self.stats.windows_recomputed += recomputed_total
        self.stats.windows_reused += len(values) - recomputed_total
        self._sequence = edited
        self._values = values
        return list(values), "region"
