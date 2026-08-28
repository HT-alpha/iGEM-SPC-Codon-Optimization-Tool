"""The elite archive: what a job reports, and why it cannot regress.

Ticket P1-30, corrected. Code: `CF-ORCH-007`.

## The guarantee

A job reports the best candidate any stage produced. Not the last stage's
output, not the output of a stage privileged by name, and not a value that
depends on the order stages happened to run in. Any stage may regress — the
report cannot.

## Why the guarantee lives here rather than in a stage ordering

The obvious alternative is to forbid a non-monotone stage from running after an
optimizing one. That fails for a reason worth stating: **non-monotone stages are
legitimate.** Mutation and repair explore, crossover explores, and NSGA is
monotone in the Pareto sense rather than in any scalarization. A pipeline that
required each stage to improve on the last would forbid exploration, which is
most of what a search does. So the invariant is placed where it survives any
future stage ordering: at the point of reporting.

It also lives in `cf-orchestrator` rather than in the library surface that first
implemented it. A selection rule computed by whichever surface happens to render
the result is a rule each surface can drift from — precisely how tier 1 and tier
2 came to disagree about an L0 term before `TierPlan` made the decision once.

## The order is total, and feasibility comes first

Three keys, in order:

1. **Feasible beats infeasible.** A high-scoring sequence that violates a
   blocking constraint is not a better answer than a lower-scoring one that
   holds; it is not an answer. Scoring alone would have let the *input* — which
   is recorded like any other candidate and may carry actionable violations —
   win on score.
2. **Higher scalarized score**, under the same scalarization the optimizing
   tiers maximize, so "better" means here what it meant during the search.
3. **Lexicographically smallest sequence digest.** Genuine ties happen, and a
   tie broken by list position makes the reported answer depend on the order
   stages were recorded in — reintroducing the order-dependence the whole
   mechanism exists to remove. A content digest is stable, has nothing to do
   with insertion order, and is already computed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

#: The `provenance.origin` vocabulary from `pareto_front.schema.json`.
ORIGINS = (
    "input",
    "random_init",
    "usage_sampled_init",
    "seed_import",
    "tier1_dp",
    "tier2_block_ca",
    "mutation",
    "crossover",
    "repair",
    "joint_reconciliation",
)


@dataclass(frozen=True)
class Candidate:
    """One sequence a stage produced, with everything selection needs."""

    #: A `provenance.origin` value, not a tier name: the front records where a
    #: sequence came from, and a stage may produce candidates of several
    #: origins.
    origin: str
    cds: str
    score: float
    feasible: bool
    #: Digest of the `sequence.schema.json` document, used as the tie-break.
    digest: str
    #: Position in the recorded order. Reporting only — never a selection key.
    recorded: int

    def to_json(self) -> dict[str, Any]:
        return {
            "origin": self.origin,
            "cds": self.cds,
            "score": self.score,
            "feasible": self.feasible,
            "digest": self.digest,
            "recorded": self.recorded,
        }

    @classmethod
    def from_json(cls, document: dict[str, Any]) -> Candidate:
        return cls(
            origin=document["origin"],
            cds=document["cds"],
            score=float(document["score"]),
            feasible=bool(document["feasible"]),
            digest=document["digest"],
            recorded=int(document["recorded"]),
        )

    def key(self) -> tuple[int, float, str]:
        """The total order. Feasibility, then score, then digest.

        Returned as something sortable rather than compared pairwise so the
        rule has exactly one statement and a caller cannot apply two of the
        three keys.
        """
        # The digest is negated in effect by sorting ascending on it while
        # sorting descending on the first two; `best` handles the direction so
        # this stays a plain tuple.
        return (int(self.feasible), self.score, self.digest)


class CandidateArchive:
    """Every candidate a job produced, and the best of them.

    Append-only. A stage records what it produced and never removes anything,
    because "the best across every stage that ran" is not a claim that can be
    made by a structure that forgets.
    """

    def __init__(self, candidates: list[Candidate] | None = None) -> None:
        self._candidates: list[Candidate] = list(candidates or [])

    def record(self, origin: str, cds: str, score: float, feasible: bool, digest: str) -> Candidate:
        if origin not in ORIGINS:
            raise ValueError(
                f"{origin!r} is not a pareto_front provenance origin; expected one of "
                f"{', '.join(ORIGINS)}"
            )
        candidate = Candidate(
            origin=origin,
            cds=cds,
            score=score,
            feasible=feasible,
            digest=digest,
            recorded=len(self._candidates),
        )
        self._candidates.append(candidate)
        return candidate

    def __len__(self) -> int:
        return len(self._candidates)

    def __iter__(self):
        return iter(self._candidates)

    @property
    def candidates(self) -> list[Candidate]:
        return list(self._candidates)

    def best(self) -> Candidate | None:
        """The reported candidate.

        Independent of the order candidates were recorded in: the same set in
        any order yields the same answer, which is what makes the guarantee
        hold for stage orderings that do not exist yet.
        """
        if not self._candidates:
            return None
        return max(
            self._candidates,
            # Ascending on feasibility and score, descending on digest, so that
            # `max` takes the feasible high scorer with the smallest digest.
            key=lambda c: (int(c.feasible), c.score, _inverted(c.digest)),
        )

    def last_recorded(self) -> Candidate | None:
        """The most recently recorded candidate, whatever it was.

        For reporting a regression: naming the stage that ran last requires
        knowing which one it was, and that is a fact about this run rather than
        a constant. The previous implementation hardcoded `tier3`.
        """
        if not self._candidates:
            return None
        return max(self._candidates, key=lambda c: c.recorded)

    def regressed(self) -> bool:
        """Whether the last stage to record is not the one being reported."""
        best, last = self.best(), self.last_recorded()
        if best is None or last is None:
            return False
        return best.digest != last.digest

    def to_json(self) -> list[dict[str, Any]]:
        return [c.to_json() for c in self._candidates]

    @classmethod
    def from_json(cls, documents: list[dict[str, Any]]) -> CandidateArchive:
        return cls([Candidate.from_json(d) for d in documents])


def _inverted(digest: str) -> tuple[int, ...]:
    """Order a digest so that `max` prefers the lexicographically smallest.

    Inverting the code points rather than negating a hash keeps the comparison
    a pure function of the string, so two processes agree — the same reason the
    tabu snapshot is sorted rather than left to `HashMap` order.
    """
    return tuple(-ord(ch) for ch in digest)
