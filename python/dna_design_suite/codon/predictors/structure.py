"""The structure-engine seam.

Ticket P2-09. `AMD-013`, `AMD-001`, `DEP-STRU-001`, `DEP-STRU-002`.

DNA Design Suite does not fold RNA. It asks something else to, and records which
something. This module is that boundary.

Why the seam lives in Python rather than in `cf-biophys-structure`
------------------------------------------------------------------
The engines that exist are C libraries with Python bindings, and the hot path is
Rust. Calling a Python engine per move would cross the boundary once per
candidate, which `ADR-0001` caps at once per generation. So the structure term is
an **L0 term** evaluated in batches through `P2-05`'s batch protocol, on the
Python side, and the Rust search never calls an engine at all — it yields a batch
and absorbs scores.

That is not a workaround. It is what the batch protocol and `P2-06`'s rationing
were built for, and the structure term is the first term expensive enough to
need them.

What the protocol deliberately does not have
--------------------------------------------
`DEP-STRU-002`. The planning contract's `StructureEngine` had
`local_accessibility_incremental(changed_nt_range)` — a method taking a range and
no sequence, which only means anything if the engine holds state from a previous
call. That assumed DNA Design Suite's own from-scratch kernel, which would have owned
such a cache. ViennaRNA's windowed API is stateless: hand it a sequence, get a
matrix. A binding cannot implement the incremental method, and pretending
otherwise would have put a cache inside every future engine instead of once in
`P2-10`.

So engines here are **stateless and pure**. Incrementality is DNA Design Suite's
problem, solved once, outside the seam.
"""

from __future__ import annotations

import importlib.util
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Protocol, runtime_checkable


class StructureEngineError(RuntimeError):
    """A structure engine could not be used. Carries a DNA Design Suite error code."""

    def __init__(self, code: str, detail: str) -> None:
        super().__init__(f"{code}: {detail}")
        self.code = code
        self.detail = detail


@dataclass(frozen=True)
class EngineIdentity:
    """What went in the lockfile.

    `DEP-STRU-001`'s reproducibility half. Which engine, at which version, with
    which parameter set, decides every energy a run reported. Two runs that
    disagree here are not comparable, however similar their inputs look, and the
    lockfile is where that becomes checkable rather than assumed.

    `parameter_set` is the engine's own identifier for the thermodynamic table
    it used. For a user-installed ViennaRNA that is the table ViennaRNA loaded,
    which DNA Design Suite never sees the contents of and never redistributes.
    """

    name: str
    version: str
    parameter_set: str

    def as_document(self) -> dict[str, str]:
        return {
            "name": self.name,
            "version": self.version,
            "parameter_set": self.parameter_set,
        }


@runtime_checkable
class StructureEngine(Protocol):
    """A folding backend.

    Every method is pure: same arguments, same answer, no retained state. That
    is what lets a binding to a stateless C library satisfy the protocol, and
    what keeps `GATE-DET-001` checkable — a stateful engine would make the
    result depend on call order, which under a parallel sweep means depending on
    thread scheduling.
    """

    def identity(self) -> EngineIdentity: ...

    def local_accessibility(
        self,
        sequence: str,
        window_nt: int,
        span_nt: int,
        temperature_c: float,
    ) -> Sequence[float]:
        """Probability that each position is inside an unpaired run of
        `window_nt`, computed over a sliding span of `span_nt`.

        One value per position. `AMD-001` makes this the only structure quantity
        the search loop may ask for.
        """
        ...


def available_engines() -> list[str]:
    """Engines this environment could use, without importing them.

    `importlib.util.find_spec` rather than a `try: import`: `TST-STRU-003`
    requires that no import of an absent package is attempted, and an import
    inside a `try` still runs the package's top-level code when it happens to be
    present under a name we did not expect.
    """
    found = []
    for module, engine in (("RNA", "viennarna"), ("linearfold", "linearfold")):
        try:
            if importlib.util.find_spec(module) is not None:
                found.append(engine)
        except (ImportError, ValueError):
            continue
    return found


def select_engine(preferred: str | None = None) -> StructureEngine:
    """Pick a structure engine, or refuse.

    There is no silent fallback and no default engine, and that is the point.
    `AMD-013` re-scoped `P2-08` from a from-scratch kernel to a binding, which
    means a machine with nothing installed has nothing to fold with. The honest
    response is `CF-STRU-003` naming what to install.

    `preferred="nussinov"` is the one opt-in exception (`DEP-STRU-006`). It is
    a combinatorial matching, not a thermodynamic fold, and it is never
    selected because ViennaRNA is missing.

    The alternative — degrading to some cheap proxy for accessibility — would
    produce a number for every candidate, so every downstream check would pass
    and the objective would be silently meaningless. That is the failure
    `DEP-STRU-001` was written against, in a different costume.
    """
    from dna_design_suite.codon.predictors.viennarna import ViennaRNAEngine

    if preferred == "nussinov":
        from dna_design_suite.codon.predictors.nussinov import NussinovEngine

        return NussinovEngine()

    engines = available_engines()
    if preferred is not None:
        if preferred not in engines:
            raise StructureEngineError(
                "CF-STRU-003",
                f"structure engine {preferred!r} was requested but is not installed; "
                f"available: {engines or 'none'}",
            )
        chosen = preferred
    elif engines:
        chosen = engines[0]
    else:
        raise StructureEngineError(
            "CF-STRU-003",
            "no structure engine is installed, and DNA Design Suite ships none. Install "
            "ViennaRNA (`pip install ViennaRNA`) and re-run. DNA Design Suite deliberately "
            "does not fall back to an approximation here: a proxy would return a "
            "number for every candidate, which would make the structure objective "
            "silently meaningless rather than visibly absent",
        )

    if chosen == "viennarna":
        return ViennaRNAEngine()
    raise StructureEngineError("CF-STRU-003", f"no binding is implemented for engine {chosen!r}")
