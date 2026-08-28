"""Job lifecycle: seal, plan, run in segments, checkpoint, cancel, report.

Ticket P1-29. Contract: `schemas/v1/job_status.schema.json`.

## The order is the guarantee

    seal lockfile → plan tiers → run segments → emit status

`TST-ORCH-001` is about the first arrow: no search work may begin before the
lockfile is sealed, because a job that resolved a model *during* the search
could observe a promotion mid-flight, which is the thing `DEP-ACT-001` exists to
prevent. Sealing first is not a convenience; it is the mechanism.

## Why the run is a loop over segments

The search engine is Rust and the lifecycle is Python. ADR-0001 caps boundary
crossings at one per generation, so this loop calls into Rust once per
*segment* — a whole number of generations — and decides between segments
whether to continue, checkpoint or stop. Cancellation and checkpointing
therefore happen at the same points, which is why there is no such thing here
as a checkpoint taken mid-cancellation: the persisted state is always a
completed generation.

## Three obligations discharged here

Three earlier tickets each wrote down something only a job could do:

* `P1-23` — `OperatorMix::warnings()` (`CF-SRCH-001`) had nowhere to go;
* `P1-18b` — the generation loop collected them into `GenerationReport::warnings`
  but could not deliver them;
* `P1-25` — `Tier2Plan::warnings()` (`CF-SRCH-010`) had the same problem.

All three now land in `job_status.warnings` through one path, each tagged with
the `stage` and `subject` of `DEP-ORCH-001` so a caller can attribute them.
"""

from __future__ import annotations

import datetime as _dt
import json
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from codonforge import _core
from codonforge.modelreg import ArtifactStore, Lockfile

from .budget import EXECUTED_TIERS, ComputeBudget, TierPlan, TierPlanError, Warning_

SCHEMA_VERSION = "1.0.0"

#: Job states, from `job_status.schema.json`.
STATES = ("queued", "running", "succeeded", "failed", "cancelled")


class Cancelled(Exception):
    """Raised when a run stops because cancellation was requested."""


@dataclass
class JobInputs:
    """The design problem, fixed at submission."""

    cds: str
    translation_table: int = 11
    five_prime_flank: str = ""
    three_prime_flank: str = ""
    constraint_set: dict[str, Any] = field(default_factory=dict)
    codon_weights: list[float] = field(default_factory=list)
    tabu_leak: float = 0.1


@dataclass
class Job:
    """One in-process design job.

    Holds nothing a resumed job could not reconstruct from its checkpoint and
    its sealed lockfile — which is the property that makes resume meaningful
    rather than merely possible.
    """

    job_id: str
    inputs: JobInputs
    budget: ComputeBudget
    lockfile: Lockfile
    plan: TierPlan
    store: ArtifactStore
    checkpoint_dir: Path | None = None
    state: str = "queued"
    created_at: str = field(
        default_factory=lambda: _dt.datetime.now(_dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    )
    started_at: str | None = None
    finished_at: str | None = None
    warnings: list[dict[str, Any]] = field(default_factory=list)
    problem: dict[str, Any] | None = None
    checkpoint: dict[str, Any] | None = None
    generations_run: int = 0
    #: Set to request a clean stop at the next segment boundary.
    cancel_requested: bool = False

    # ----------------------------------------------------------------- #
    # submission
    # ----------------------------------------------------------------- #

    @classmethod
    def submit(
        cls,
        job_id: str,
        inputs: JobInputs,
        budget: ComputeBudget,
        lockfile: Lockfile,
        store: ArtifactStore,
        terms: list[str],
        capabilities: dict[str, str],
        radii: dict[str, int | None] | None = None,
        block_size: int | None = None,
        directions: dict[str, str] | None = None,
        checkpoint_dir: Path | None = None,
    ) -> Job:
        """Seal, plan, and return a job ready to run.

        The lockfile must already be sealed: a job that could seal it later
        could also resolve a model later, and the whole point of pinning at
        submission is that there is no later.
        """
        if lockfile.digest is None:
            raise ValueError(
                "the landscape lockfile must be sealed before a job is created; "
                "an unsealed lockfile cannot pin anything"
            )

        plan = TierPlan.plan(budget, terms, capabilities, radii, block_size, directions)
        unwired = [tier for tier in budget.stages() if tier not in EXECUTED_TIERS]
        job = cls(
            job_id=job_id,
            inputs=inputs,
            budget=budget,
            lockfile=lockfile,
            plan=plan,
            store=store,
            checkpoint_dir=checkpoint_dir,
        )
        job.warnings.extend(w.to_problem_details() for w in plan.warnings)
        for tier in unwired:
            job.warnings.append(
                Warning_(
                    code="CF-ORCH-006",
                    title="A tier in the budget plan is not executed by this build",
                    detail=(
                        f"the {budget.tier} budget plans {tier}, which this build does not "
                        "run from a job: the tier exists in cf-search but the full pipeline "
                        "is assembled by P1-30. The job optimizes with the generation loop "
                        "alone, which is less search than the plan describes"
                    ),
                    stage="tier_plan",
                    subject=tier,
                ).to_problem_details()
            )
        return job

    # ----------------------------------------------------------------- #
    # the RNG plan
    # ----------------------------------------------------------------- #

    def root_seed(self) -> int:
        """The seed every stream in this job derives from.

        Taken from the sealed lockfile and from nowhere else. A job layer that
        opened its own seed source would make the lockfile's `rng` block a
        description of something that did not happen, and two jobs replayed
        from the same lockfile would diverge while both looked reproducible.
        """
        return int(self.lockfile.rng["root_seed"])

    def rng_plan(self) -> dict[str, Any]:
        """The `rng` block this job is running under.

        Returned rather than reconstructed so a caller can check that the
        algorithm named in the lockfile is the one the engine implements —
        `DEP-DET-001` exists because it once was not.
        """
        return dict(self.lockfile.rng)

    # ----------------------------------------------------------------- #
    # execution
    # ----------------------------------------------------------------- #

    def resolved_models(self) -> dict[str, Any]:
        """Replay the sealed lockfile against the store.

        By digest, never by coordinate. This is the call that makes a resumed
        job blind to a promotion that happened while it was checkpointed: the
        digests were fixed at submission and `Lockfile.replay` will not
        substitute a newer artifact for a missing one, it fails.
        """
        return self.lockfile.replay(self.store)

    def run(self, on_segment: Callable[[Job], None] | None = None) -> Job:
        """Run to the budget, a segment at a time.

        `on_segment` is called after each segment with the job, which is where
        a test injects a cancellation or a mid-flight model publication.
        """
        if self.state == "cancelled":
            return self

        # TST-ORCH-001: nothing above this line touches the search, and the
        # lockfile is replayed here so a job whose pinned artifacts have gone
        # fails before doing any work rather than partway through.
        self.resolved_models()

        self.state = "running"
        self.started_at = _now()

        try:
            while not self.budget.exhausted:
                if self.cancel_requested:
                    raise Cancelled()

                generations = self.budget.next_segment()
                result = json.loads(
                    _core.search_run_segment(
                        json.dumps(self.checkpoint) if self.checkpoint else None,
                        self.inputs.cds,
                        self.inputs.translation_table,
                        self.inputs.five_prime_flank,
                        self.inputs.three_prime_flank,
                        json.dumps(self.inputs.constraint_set),
                        self.inputs.codon_weights,
                        self.inputs.tabu_leak,
                        self.root_seed(),
                        generations,
                    )
                )

                # Committed as one step: the checkpoint, the budget and the
                # counter move together or not at all. A cancellation between
                # any two of them would leave a job whose recorded progress and
                # recorded state disagreed.
                self.checkpoint = result["checkpoint"]
                self.budget.spend(result["generations_run"])
                self.generations_run += result["generations_run"]
                self._absorb_warnings(result.get("warnings", []))
                self.persist_checkpoint()

                if on_segment is not None:
                    on_segment(self)
        except Cancelled:
            # A clean stop: the last completed segment is the state, and it was
            # already persisted. Nothing partial exists to roll back, because
            # nothing partial was ever written.
            self.state = "cancelled"
            self.finished_at = _now()
            self.problem = {
                "type": "https://codonforge.org/problems/cf-orch-003",
                "title": "Job was cancelled",
                "status": 200,
                "code": "CF-ORCH-003",
                "detail": (
                    f"cancelled after {self.generations_run} generations; the checkpoint "
                    "is a completed generation and the job can be resumed from it"
                ),
                "stage": "tier3",
                "subject": self.job_id,
            }
            return self
        except TierPlanError as error:
            self.state = "failed"
            self.finished_at = _now()
            self.problem = error.to_problem_details()
            return self

        self.state = "succeeded"
        self.finished_at = _now()
        return self

    def _absorb_warnings(self, incoming: list[dict[str, Any]]) -> None:
        """Add warnings the engine produced, without duplicating them.

        A segment reports the same operator-availability warnings every time it
        runs, because they are properties of the mix rather than of the segment.
        Appending them per segment would make a long job's status a list of the
        same two sentences forty times over, which is how a warning stops being
        read.
        """
        seen = {(w.get("code"), w.get("subject")) for w in self.warnings}
        for warning in incoming:
            key = (warning.get("code"), warning.get("subject"))
            if key in seen:
                continue
            seen.add(key)
            self.warnings.append(warning)

    # ----------------------------------------------------------------- #
    # checkpoint and resume
    # ----------------------------------------------------------------- #

    def checkpoint_path(self) -> Path | None:
        if self.checkpoint_dir is None:
            return None
        return Path(self.checkpoint_dir) / f"{self.job_id}.checkpoint.json"

    def persist_checkpoint(self) -> None:
        path = self.checkpoint_path()
        if path is None or self.checkpoint is None:
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        # Written whole, then moved into place. A partially written checkpoint
        # is worse than no checkpoint: resume would read it and produce a run
        # that is neither the original nor a fresh one.
        temporary = path.with_suffix(".partial")
        temporary.write_text(
            json.dumps(
                {
                    "job_id": self.job_id,
                    "lockfile_digest": self.lockfile.digest,
                    "generations_run": self.generations_run,
                    "budget_consumed": self.budget.consumed,
                    "state": self.checkpoint,
                },
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )
        temporary.replace(path)

    def resume(self, checkpoint_document: dict[str, Any]) -> Job:
        """Restore from a checkpoint written by this job's lockfile.

        The lockfile digest is checked because a checkpoint carries no model
        pins of its own: resuming one job's state under another's lockfile
        would silently change which artifacts produced the numbers, which is
        the same substitution `CF-RESOLVE-006` refuses at replay.
        """
        if checkpoint_document.get("job_id") != self.job_id:
            raise ValueError(
                f"checkpoint belongs to job {checkpoint_document.get('job_id')!r}, "
                f"not {self.job_id!r}"
            )
        if checkpoint_document.get("lockfile_digest") != self.lockfile.digest:
            raise ValueError(
                "checkpoint was written under a different landscape lockfile; resuming "
                "would change which artifacts produced the numbers"
            )
        self.checkpoint = checkpoint_document["state"]
        self.generations_run = int(checkpoint_document["generations_run"])
        self.budget.consumed = int(checkpoint_document["budget_consumed"])
        return self

    @classmethod
    def load_checkpoint(cls, path: Path) -> dict[str, Any]:
        return json.loads(Path(path).read_text(encoding="utf-8"))

    # ----------------------------------------------------------------- #
    # status
    # ----------------------------------------------------------------- #

    def status(self) -> dict[str, Any]:
        """A `job_status.schema.json` document."""
        document: dict[str, Any] = {
            "schema_version": SCHEMA_VERSION,
            "job_id": self.job_id,
            "state": self.state,
            "created_at": self.created_at,
            "landscape_lockfile_digest": self.lockfile.digest,
            "resolved_objectives": [
                {
                    "term_id": term,
                    "role": "pareto_axis",
                    "direction": self.plan.directions[term],
                }
                for term in self.plan.scored
            ],
            "progress": self._progress(),
            "warnings": self.warnings,
            "links": {"self": f"/v1/jobs/{self.job_id}"},
        }
        if self.started_at:
            document["started_at"] = self.started_at
        if self.finished_at:
            document["finished_at"] = self.finished_at
        if self.problem:
            document["problem"] = self.problem
        return document

    def _progress(self) -> dict[str, Any]:
        """Per-tier progress, reporting unrun tiers as skipped rather than done."""
        running = self.state == "running"
        tiers = []
        for tier in self.budget.stages():
            if tier not in EXECUTED_TIERS:
                tiers.append({"tier": tier, "state": "skipped"})
                continue
            if self.state in ("succeeded",):
                state = "done"
            elif self.generations_run > 0 or running:
                state = "running" if running else "pending"
            else:
                state = "pending"
            tiers.append(
                {
                    "tier": tier,
                    "state": state,
                    "generation": self.generations_run,
                    "max_generations": self.budget.total_generations or None,
                    "fraction": (
                        self.generations_run / self.budget.total_generations
                        if self.budget.total_generations
                        else None
                    ),
                }
            )
        current = "none"
        if running:
            current = EXECUTED_TIERS[0]
        return {"current_tier": current, "tiers": tiers}

    def warnings_from(self, stage: str) -> list[dict[str, Any]]:
        """Warnings a given stage produced. `DEP-ORCH-001`'s whole point."""
        return [w for w in self.warnings if w.get("stage") == stage]


def _now() -> str:
    return _dt.datetime.now(_dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
