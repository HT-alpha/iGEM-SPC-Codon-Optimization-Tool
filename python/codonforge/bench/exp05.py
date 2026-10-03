"""The `EXP-05` decision rule, applied to measurements the harness produced.

The rule is the one in the experiment register, unchanged: switch the default
to Baldwinian iff it wins on hypervolume by more than the Lamarckian
seed-to-seed standard deviation AND its repair-invocation count is under 2×.
Population diversity is recorded beside those two; it is part of what the
experiment reports and it is not a third vote that can overturn them.

Fewer than 20 seeds per arm is not a decision. The register requires 20, and
a shorter run would be a different experiment. This function does not change
the repair default; a caller that wants to change it has to do that in the
search configuration after reading the result.
"""

from __future__ import annotations

import statistics
from dataclasses import dataclass

from codonforge.bench.reference import BenchError

SEED_COUNT = 20


@dataclass(frozen=True)
class Trial:
    hypervolume: float
    diversity: float
    repair_invocations: int


@dataclass(frozen=True)
class Decision:
    switch_to_baldwinian: bool
    hv_margin: float
    hv_standard_deviation: float
    repair_ratio: float
    mean_diversity_lamarckian: float
    mean_diversity_baldwinian: float


def decide(lamarckian: list[Trial], baldwinian: list[Trial]) -> Decision:
    if len(lamarckian) < SEED_COUNT or len(baldwinian) < SEED_COUNT:
        raise BenchError(
            f"EXP-05 needs {SEED_COUNT} seeds on each arm; "
            f"got {len(lamarckian)} and {len(baldwinian)}"
        )
    left = [trial.hypervolume for trial in lamarckian]
    right = [trial.hypervolume for trial in baldwinian]
    margin = statistics.fmean(right) - statistics.fmean(left)
    deviation = statistics.pstdev(left)
    repairs_left = sum(trial.repair_invocations for trial in lamarckian)
    repairs_right = sum(trial.repair_invocations for trial in baldwinian)
    ratio = (repairs_right / repairs_left) if repairs_left else float("inf")
    switch = margin > deviation and ratio < 2.0
    return Decision(
        switch_to_baldwinian=switch,
        hv_margin=margin,
        hv_standard_deviation=deviation,
        repair_ratio=ratio,
        mean_diversity_lamarckian=statistics.fmean(t.diversity for t in lamarckian),
        mean_diversity_baldwinian=statistics.fmean(t.diversity for t in baldwinian),
    )
