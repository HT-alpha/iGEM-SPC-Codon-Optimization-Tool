"""Compute budget and tier plan.

Ticket P1-29. Specification: ADR-0008 (the degradation ladder) and
`04-model-resolution.md` §2.1 (what a planner does with a term it cannot afford
to score).

WHY THE PLAN IS AN OBJECT AND NOT A FUNCTION CALL. §2.1 says "the planner"
decides whether an unaffordable term is dropped with a warning or made fatal.
Until now there was no planner: tier 1 refused an L0 term outright with
`CF-LAND-002` and tier 2 dropped it with `CF-SRCH-010`, so the same term in the
same job got two different answers depending on which tier looked at it first.
`TierPlan` makes the decision once, records why, and hands each tier a term list
it can actually score.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

#: Named budget tiers, from `common.schema.json#/$defs/BudgetTier`.
TIERS = ("fast", "standard", "exact")

#: Which search tiers each budget tier runs.
#:
#: Read from ADR-0008's ladder rather than invented here: `fast` seeds and
#: polishes, `standard` adds the population loop, `exact` additionally
#: re-ranks on the true objective. A tier absent from this map is not run, and
#: a caller asking for one gets told rather than silently given less.
TIER_PLAN: dict[str, tuple[str, ...]] = {
    "fast": ("tier1_seed", "tier2_block_ca"),
    "standard": ("tier1_seed", "tier2_block_ca", "tier3_moea"),
    "exact": ("tier1_seed", "tier2_block_ca", "tier3_moea", "decision"),
}

#: Tiers the job runner actually executes.
#:
#: `P1-29` delivered the lifecycle and drove the generation loop alone, so tier
#: 1 and tier 2 were reported `skipped` with a `CF-ORCH-006` warning naming what
#: was missing rather than reporting success over less search than the plan
#: described. `P1-30` assembles them, and those warnings disappearing for a job
#: that runs them is the evidence the wiring is real rather than merely present
#: - which is what `TST-CLI-003` asserts.
#:
#: `decision` remains unassembled: the reduction and knee-selection it names are
#: `cf-decision`'s and belong to `P2-08`. An `exact` job still warns about it.
EXECUTED_TIERS: tuple[str, ...] = ("tier1_seed", "tier2_block_ca", "tier3_moea")

#: The only tier whose progress is denominated in generations.
#:
#: Tier 1 runs once and tier 2 runs sweeps, so reporting the loop's
#: generation count against them would state something untrue about work
#: they did not do.
GENERATIONAL_TIER = "tier3_moea"


@dataclass
class Warning_:
    """One `job_status.warnings` entry, with provenance.

    `code` separates classes of problem and `stage` separates producers of the
    same class — the tier-1 and tier-2 capability prechecks emit the same codes
    for the same reasons, so a status document carrying several is unreadable
    without both. `DEP-ORCH-001`.
    """

    code: str
    title: str
    detail: str
    stage: str
    subject: str
    status: int = 200

    def to_problem_details(self) -> dict[str, Any]:
        return {
            "type": f"https://codonforge.org/problems/{self.code.lower()}",
            "title": self.title,
            "status": self.status,
            "code": self.code,
            "detail": self.detail,
            "stage": self.stage,
            "subject": self.subject,
        }


@dataclass
class ComputeBudget:
    """How much work a job may do, and how much is left.

    Generations are the unit because they are the only unit the search engine
    can be interrupted at: ADR-0001 caps boundary crossings at one per
    generation, so a budget denominated in anything finer could not be enforced
    without violating it. Wall-clock and evaluation budgets belong on top of
    this, not instead of it.
    """

    tier: str
    total_generations: int
    #: Generations per segment between checkpoints.
    segment: int = 8
    consumed: int = 0

    def __post_init__(self) -> None:
        if self.tier not in TIERS:
            raise ValueError(f"unknown budget tier {self.tier!r}; expected one of {TIERS}")
        if self.total_generations < 0:
            raise ValueError("total_generations must not be negative")
        if self.segment < 1:
            raise ValueError("segment must be at least one generation")

    @property
    def remaining(self) -> int:
        return max(0, self.total_generations - self.consumed)

    @property
    def exhausted(self) -> bool:
        return self.remaining == 0

    def next_segment(self) -> int:
        """Generations to run next: the segment size, or what is left."""
        return min(self.segment, self.remaining)

    def spend(self, generations: int) -> None:
        if generations > self.remaining:
            raise ValueError(
                f"cannot spend {generations} generations against {self.remaining} remaining"
            )
        self.consumed += generations

    def stages(self) -> tuple[str, ...]:
        return TIER_PLAN[self.tier]


@dataclass
class TierPlan:
    """The capability decision, taken once for the whole job.

    Splitting it per tier is what produced the inconsistency this replaces. The
    same term, in the same job, now gets the same answer whichever tier asks.
    """

    budget: ComputeBudget
    scored: list[str] = field(default_factory=list)
    #: Terms tier 1's DP scores. In-core only: the DP needs a delta table.
    tier1: list[str] = field(default_factory=list)
    #: Terms tier 2's block ascent scores. In-core only, for the same reason.
    tier2: list[str] = field(default_factory=list)
    #: Requested terms that are computed and reported but never ranked on.
    diagnostics: list[str] = field(default_factory=list)
    dropped: list[str] = field(default_factory=list)
    warnings: list[Warning_] = field(default_factory=list)
    #: Optimization direction per scored term, as its provider declared it.
    #:
    #: Carried rather than derived. `job_status.resolved_objectives` requires a
    #: direction, and there is no defensible default: guessing `maximize` for a
    #: term the provider declared `minimize` would invert the objective in the
    #: document a user reads while the search did the opposite.
    directions: dict[str, str] = field(default_factory=dict)

    @classmethod
    def plan(
        cls,
        budget: ComputeBudget,
        terms: Iterable[str],
        capabilities: dict[str, str],
        radii: dict[str, int | None] | None = None,
        block_size: int | None = None,
        directions: dict[str, str] | None = None,
        roles: dict[str, str] | None = None,
    ) -> TierPlan:
        """Decide which terms each tier scores.

        `capabilities` maps a term id to `L0`, `L1` or `L2`; `radii` maps it to
        an interaction radius in codons, or `None` for unbounded. Both come from
        the registry rather than from configuration, because a capability a
        caller asserts is not a capability.
        """
        radii = radii or {}
        directions = directions or {}
        roles = roles or {}
        plan = cls(budget=budget)

        for term in terms:
            # A term may only be ranked on if its registration says so.
            #
            # `CF-LAND-014`. The registry already refuses to *register* a
            # provider in a role above its declared ceiling; nothing stopped the
            # orchestrator from then treating a diagnostic registration as a
            # Pareto axis. `usage.gc_content` is the case that exposed it: GC
            # content has no direction - too high folds the transcript, too low
            # does the same for other reasons - so its provider deliberately
            # caps itself at `diagnostic` and says in as many words that this
            # "keeps it out of any ranking". The search ranked on it anyway,
            # minimizing a quantity whose direction the biology does not have.
            #
            # Demoted rather than refused: the number is still worth reporting,
            # and refusing would reject jobs that have been accepted for
            # several rounds. But it stops being an axis, and the warning says
            # so with the same (code, stage, subject) provenance CF-ORCH-006
            # uses for a skipped tier.
            role = roles.get(term, "pareto_axis")
            if role != "pareto_axis":
                plan.warnings.append(
                    Warning_(
                        code="CF-LAND-014",
                        title="A term was requested as an objective but cannot be ranked on",
                        detail=(
                            f"{term} is registered as {role!r}, so it is computed and "
                            f"reported but takes no part in the search. Ranking on it "
                            f"would order designs by a quantity its provider declares "
                            f"has no meaningful direction"
                        ),
                        stage="tier_plan",
                        subject=term,
                    )
                )
                plan.diagnostics.append(term)
                continue

            capability = capabilities.get(term)
            if capability is None:
                raise KeyError(f"no provider is registered for term {term!r}")

            if capability == "L0":
                # Kept, but only for tier 3. `P2-11`.
                #
                # This used to drop the term outright with CF-SRCH-010, and
                # that was right while nothing could score an out-of-core term
                # at all. `P2-05a` changed the fact the rule rested on: tier 3
                # now evaluates such a term once per generation on the whole
                # population, which is exactly the amortization §2.1 says makes
                # it affordable. Tiers 1 and 2 still cannot take it - both need
                # a delta table, and an L0 provider has none - so it is absent
                # from their lists rather than from the search.
                plan.scored.append(term)
                plan.directions[term] = directions.get(term, "minimize")
                plan.warnings.append(
                    Warning_(
                        code="CF-SRCH-010",
                        title="A term is scored by the population loop only",
                        detail=(
                            f"{term}'s provider is L0, so it has no exact delta and cannot "
                            "be scored by the seeding DP or the block ascent. Tier 3 scores "
                            "it once per generation on the whole population; the seeds those "
                            "earlier tiers produce were chosen without it"
                        ),
                        stage="tier_plan",
                        subject=term,
                    )
                )
                continue

            radius = radii.get(term, 0)
            if block_size is not None and (radius is None or radius > block_size):
                # Not a warning: tier 2's partition is unsound under such a
                # term, and running anyway would produce a plausible answer
                # that is not the coordinate-ascent step it claims to be.
                raise TierPlanError(
                    code="CF-SRCH-008",
                    detail=(
                        f"{term} couples positions "
                        + ("without bound" if radius is None else f"up to {radius} codons apart")
                        + f", against a tier 2 block size of {block_size}. Raise the block "
                        "size or drop the term; the overlapping-block reconciliation that "
                        "answers this case arrives with P1-26"
                    ),
                    subject=term,
                )
            direction = directions.get(term)
            if direction not in ("minimize", "maximize"):
                raise ValueError(
                    f"term {term!r} has no declared optimization direction; the registry "
                    "must supply one, because job_status.resolved_objectives requires it "
                    "and there is no defensible default"
                )
            plan.scored.append(term)
            plan.tier1.append(term)
            plan.tier2.append(term)
            plan.directions[term] = direction

        return plan


class TierPlanError(Exception):
    """A planning failure that refuses the job rather than warning about it.

    Kept distinct from `Warning_` on purpose. `CF-SRCH-008` is sometimes
    described alongside the warning codes, but a term wider than the block does
    not degrade the search — it invalidates the partition, so it belongs in
    `job_status.problem` and stops the job.
    """

    def __init__(self, code: str, detail: str, subject: str) -> None:
        self.code = code
        self.detail = detail
        self.subject = subject
        super().__init__(f"{code}: {detail}")

    def to_problem_details(self) -> dict[str, Any]:
        return {
            "type": f"https://codonforge.org/problems/{self.code.lower()}",
            "title": "A term couples positions more widely than the tier 2 block size",
            "status": 422,
            "code": self.code,
            "detail": self.detail,
            "stage": "tier_plan",
            "subject": self.subject,
        }
