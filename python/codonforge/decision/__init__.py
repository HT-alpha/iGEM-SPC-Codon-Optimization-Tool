"""cf-decision — Pareto front post-processing and export.

Phase 1 delivers the export half (`P1-30`). Probabilistic dominance, robustness
scoring and knee detection are `P2-08`'s, and this module does not pretend to
them: a Phase 1 front is emitted with `reduction.method: "none"` because no
reduction has been applied, not because none was wanted.

Declared dependencies: cf-types, cf-metrics, cf-landscape-api.
"""

from .front import build, candidate

__all__ = ["build", "candidate"]
