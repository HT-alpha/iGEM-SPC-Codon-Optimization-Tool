"""When two benchmark records may be compared. `TST-BEN-002`.

Hypervolume is comparable only at one budget tier, one lockfile and one
reference point. A record that differs in a label (the protein, the arm, the
seed) is still comparable: those are the things a paired test is about.
"""

from __future__ import annotations

from dataclasses import dataclass

from codonforge.bench.reference import BenchError


@dataclass(frozen=True)
class RunRecord:
    arm: str
    protein: str
    seed: int
    budget_tier: str
    lockfile_digest: str
    reference_point: tuple[float, ...]
    hypervolume: float


def assert_comparable(left: RunRecord, right: RunRecord) -> None:
    """Refuse the comparison, or return.

    The three fields are checked independently so a failure names the one
    that differed. A helper that returned a boolean would let a caller
    compare anyway.
    """
    if left.budget_tier != right.budget_tier:
        raise BenchError(
            f"budget tier {left.budget_tier!r} is not {right.budget_tier!r}; "
            "hypervolume at two tiers is not one comparison"
        )
    if left.lockfile_digest != right.lockfile_digest:
        raise BenchError(
            "the two fronts were scored under different lockfiles; the numbers "
            "are not the same landscape"
        )
    if left.reference_point != right.reference_point:
        raise BenchError(
            "the two hypervolumes use different reference points and cannot be compared"
        )
