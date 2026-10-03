"""cf-decision — Pareto front post-processing and export.

Phase 1 delivers the export half (`P1-30`). `P2-15` is the reduction: reference
directions, knee, clustering and an ε-grid, with extreme points retained
(`TST-DEC-003`). A document with `reduction.method: "none"` is the unreduced
set, which is what a job stores. Presentation calls `present`.

Probabilistic dominance is still not this module.

Declared dependencies: cf-types, cf-metrics, cf-landscape-api.
"""

from .front import build, candidate
from .reduce import extreme_ids, present, reduce_candidates

__all__ = ["build", "candidate", "extreme_ids", "present", "reduce_candidates"]
