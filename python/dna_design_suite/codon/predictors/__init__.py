"""Predictor bindings that sit behind DNA Design Suite's plugin seams.

Nothing here is vendored. Every engine is user-installed and detected at
runtime, per `AMD-013`.
"""

from dna_design_suite.codon.predictors.structure import (
    EngineIdentity,
    StructureEngine,
    StructureEngineError,
    available_engines,
    select_engine,
)

__all__ = [
    "EngineIdentity",
    "StructureEngine",
    "StructureEngineError",
    "available_engines",
    "select_engine",
]
