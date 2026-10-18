"""Opt-in combinatorial structure score.

Nussinov and Jacobson, PNAS 1980, 77:6309–6313. Pairs are A–U, G–C and G–U,
each worth 1, with a hairpin of at least three unpaired bases. No stacking
and no Turner parameters. The objective value is the mean paired fraction
in [0, 1]. It is not a free energy and not an unpaired probability.

`DEP-STRU-006`. `select_engine()` does not choose this engine.
"""

from __future__ import annotations

import json
from collections.abc import Sequence

from dna_design_suite import _core
from dna_design_suite.codon.predictors.structure import EngineIdentity, StructureEngineError


class NussinovEngine:
    """Maximum matching over a sliding span. Temperature is ignored."""

    def identity(self) -> EngineIdentity:
        return EngineIdentity(
            name="nussinov",
            version="0.1.0",
            parameter_set="unit_pair_AU_GC_GU_min_hairpin_3",
        )

    def local_accessibility(
        self,
        sequence: str,
        window_nt: int,
        span_nt: int,
        temperature_c: float,
    ) -> Sequence[float]:
        del temperature_c
        if window_nt < 1 or span_nt < window_nt or window_nt > len(sequence):
            raise StructureEngineError(
                "CF-STRU-004",
                f"window_nt {window_nt} and span_nt {span_nt} do not fit a "
                f"sequence of {len(sequence)} nt",
            )
        try:
            values = _core.nussinov_accessibility(sequence, window_nt, span_nt)
        except ValueError as exc:
            raise StructureEngineError("CF-STRU-004", str(exc)) from exc
        if isinstance(values, str):
            values = json.loads(values)
        return [float(value) for value in values]


def paired_fraction_penalty(unpaired: Sequence[float]) -> float:
    """Mean paired fraction. Minimized. Not a log probability."""
    if not unpaired:
        raise StructureEngineError("CF-STRU-004", "no accessibility values to summarize")
    return sum(1.0 - float(value) for value in unpaired) / len(unpaired)
