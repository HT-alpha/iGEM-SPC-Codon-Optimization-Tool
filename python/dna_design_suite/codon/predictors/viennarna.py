"""ViennaRNA binding for the structure seam.

Ticket P2-08, re-scoped by `AMD-013`. Contract:
`python/dna_design_suite/codon/predictors/structure.py`.

DNA Design Suite does not vendor ViennaRNA and does not redistribute it. The user
installs it; this module detects and calls it. That arrangement is what makes
ViennaRNA's licence condition 1 — no redistribution for a fee beyond media costs
— satisfied by construction rather than by promise, and it is why the `.par`
parameter files' own licensing status never has to be resolved: ViennaRNA loads
its own tables and DNA Design Suite never handles them.

Attribution is a live obligation under condition 2. See `NOTICE`.

What this module is careful about
---------------------------------
Three things that would each produce plausible numbers while being wrong:

1. **1-based indexing.** `pfl_fold_up` returns a 1-based matrix. An off-by-one
   here shifts every accessibility by one nucleotide, which is a third of a
   codon, and nothing downstream would look unusual.
2. **Impossible windows read as certainty.** `up[i][u]` is exactly `0.0` when
   `u > i`, because no run of length `u` can end at position `i`. Passing that
   through would report the 5' edge as maximally *inaccessible* when the truth
   is that the question was not asked. This module returns one value per window
   *end* position instead, so the undefined positions have no entry rather than
   a confident wrong one.
3. **Temperature.** `RNA.cvar.temperature` is a process global. It is set on
   every call rather than once at import, so a caller cannot inherit a
   temperature some earlier call left behind.
"""

from __future__ import annotations

import threading
from collections.abc import Sequence

from dna_design_suite.codon.predictors.structure import EngineIdentity, StructureEngineError

#: `RNA.cvar.temperature` is process-global. One fold at a time keeps two
#: structure jobs from writing each other's temperature. CAI-only work does
#: not take this lock, because it does not fold.
FOLD_LOCK = threading.Lock()


class ViennaRNAEngine:
    """Windowed local accessibility via ViennaRNA's sliding-window partition
    function (`RNAplfold`'s algorithm, through `RNA.pfl_fold_up`).

    Deliberately exposes no global-MFE method. `AMD-001` keeps global MFE out of
    the search loop on cost grounds, and the cheapest way to keep it out is to
    give the loop no way to ask: a `global_mfe` here would sooner or later be
    called from a per-move path by someone who had a good reason at the time.

    **Not thread-safe within a process.** ViennaRNA's temperature is a library
    global, so two threads folding at different temperatures would race. The
    structure term runs as an L0 batch in a worker rather than inside the
    parallel block sweep, so this does not conflict with `GATE-DET-001` - but it
    is a real constraint on where this class may be called from.
    """

    def __init__(self) -> None:
        try:
            import RNA
        except ImportError as exc:  # pragma: no cover - guarded by select_engine
            raise StructureEngineError(
                "CF-STRU-003",
                "ViennaRNA is not installed; DNA Design Suite does not ship it "
                "(AMD-013: user-installed backend, never vendored)",
            ) from exc
        self._rna = RNA

    def identity(self) -> EngineIdentity:
        """What the lockfile records.

        `parameter_set` is ViennaRNA's own default table for the version
        installed. DNA Design Suite names it rather than reading it: the point is that
        two runs can be compared, not that DNA Design Suite holds the numbers.
        """
        return EngineIdentity(
            name="viennarna",
            version=str(self._rna.__version__),
            parameter_set=f"viennarna-default@{self._rna.__version__}",
        )

    def local_accessibility(
        self,
        sequence: str,
        window_nt: int,
        span_nt: int,
        temperature_c: float,
    ) -> Sequence[float]:
        """Probability of an unpaired run of `window_nt` ending at each position
        where such a run fits.

        Returns `len(sequence) - window_nt + 1` values. Index `k` is the run
        ending at 0-based sequence position `k + window_nt - 1`.
        """
        if window_nt < 1:
            raise StructureEngineError(
                "CF-STRU-004", f"window_nt must be at least 1, got {window_nt}"
            )
        if window_nt > len(sequence):
            raise StructureEngineError(
                "CF-STRU-004",
                f"window_nt {window_nt} exceeds the sequence length {len(sequence)}; "
                "no run of that length fits anywhere, so there is nothing to report",
            )
        if span_nt < window_nt:
            raise StructureEngineError(
                "CF-STRU-004",
                f"span_nt {span_nt} is smaller than window_nt {window_nt}; a base-pair "
                "span shorter than the window asks for the accessibility of a region "
                "larger than the structure being considered",
            )

        # ViennaRNA scores RNA. T and U happen to give identical results in this
        # build, but relying on that would make the binding depend on an
        # undocumented convenience.
        rna = sequence.upper().replace("T", "U")

        # 1-based. Position i holds the probabilities for runs of length 1..u
        # ending at i, so a full window first fits at i == window_nt.
        # The matrix is only valid while this thread owns the engine lock:
        # another fold would replace the temperature and the buffer.
        with FOLD_LOCK:
            self._rna.cvar.temperature = float(temperature_c)
            matrix = self._rna.pfl_fold_up(rna, window_nt, span_nt, span_nt)
            values = [float(matrix[end][window_nt]) for end in range(window_nt, len(rna) + 1)]
        return values


def accessibility_penalty(accessibility: Sequence[float]) -> float:
    """Collapse per-window accessibility into one objective value.

    Mean negative log accessibility, so the term is *minimized* and a strongly
    paired region costs more than a weakly paired one. The log matters: raw
    probabilities near zero differ by orders of magnitude that a linear mean
    would flatten into indistinguishable near-zeros, and those are exactly the
    regions the objective exists to find.
    """
    if not accessibility:
        raise StructureEngineError("CF-STRU-004", "no accessibility values to summarize")
    import math

    # Floor rather than skip. A window with probability 0 to within floating
    # point is the most inaccessible window there is, and dropping it would make
    # the worst case invisible; 1e-300 keeps the log finite without inventing a
    # different physical claim.
    total = sum(-math.log(max(value, 1e-300)) for value in accessibility)
    return total / len(accessibility)
