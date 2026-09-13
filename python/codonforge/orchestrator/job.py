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

from .archive import CandidateArchive
from .budget import (
    EXECUTED_TIERS,
    GENERATIONAL_TIER,
    ComputeBudget,
    TierPlan,
    TierPlanError,
    Warning_,
)

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
    #: `host_profile.schema.json` document, for the terms that need one.
    host_profile: dict[str, Any] = field(default_factory=dict)
    host_id: str = "escherichia_coli:511145/BL21"
    host_digest: str = "sha256:" + "0" * 64
    #: Scalarization coefficients, one per scored term.
    term_weights: list[float] = field(default_factory=list)
    #: Tier 3 population size and stagnation window.
    population: int = 12
    stagnation_window: int = 8
    #: `E_max` in §6 - exact out-of-core evaluations per generation.
    #:
    #: Defaults to the population, so a job that asks for a structure objective
    #: scores every candidate rather than silently searching on a partially
    #: UNKNOWN axis by default. Rationing below this is a deliberate choice a
    #: caller makes, and `CF-SRCH-016` reports what it cost.
    l0_budget: int = 12
    #: Tier 1 beam width and seed count; tier 2 block size and sweep cap.
    tier1_beam: int = 32
    tier1_seeds: int = 4
    tier2_block: int = 8
    tier2_sweeps: int = 4


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
    #: Per-stage outcome, keyed by tier name. Populated as stages complete.
    stage_reports: dict[str, Any] = field(default_factory=dict)
    #: The most recent tier-3 segment result, for progress reporting.
    tier3: dict[str, Any] | None = None
    #: The sequence tiers 1 and 2 produced, handed to the generation loop.
    seeded_cds: str | None = None
    #: Every candidate any stage produced. The job reports the best of them.
    #:
    #: The pipeline is not monotone and is not required to be: mutation and
    #: repair explore, crossover explores, and NSGA is monotone in the Pareto
    #: sense rather than in any scalarization. So the invariant is placed at the
    #: point of reporting rather than enforced on the stage order, where it
    #: holds for orderings that do not exist yet. See `archive.py`.
    archive: CandidateArchive = field(default_factory=CandidateArchive)
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
        roles: dict[str, str] | None = None,
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

        plan = TierPlan.plan(budget, terms, capabilities, radii, block_size, directions, roles)
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

    def tier3_seed(self) -> int:
        """The RNG root the generation loop and tier 3 derive from.

        `DEP-DET-003`. A method rather than an expression inlined at the call
        site, so a test can assert the job *uses* the namespaced root — a test
        that only checked the derivation function passes whichever root the
        job actually passes, which is how a negative control reverting this to
        the job seed went unnoticed.
        """
        return int(_core.tier3_root_seed(self.root_seed()))

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
            # Tiers run in plan order, each crossing the boundary once. They
            # are skipped on resume: a checkpoint is a state of the generation
            # loop, and re-seeding a resumed job would discard the very work
            # the checkpoint exists to preserve.
            if self.checkpoint is None:
                self._run_seeding_tiers()

            while not self.budget.exhausted:
                if self.cancel_requested:
                    raise Cancelled()

                generations = self.budget.next_segment()
                # P2-04: tier 3 is the job's population loop. The mutation-only
                # driver it replaces could not produce a front at all - it kept
                # one incumbent - so a multi-objective job reported a single
                # design however many objectives it had been given.
                scorer, axis = self._l0_scorer()
                result = json.loads(
                    _core.search_tier3_segment(
                        self._tier_request(self.seeded_cds or self.inputs.cds),
                        json.dumps(self.checkpoint) if self.checkpoint else None,
                        self.inputs.population,
                        generations,
                        self.inputs.stagnation_window,
                        scorer,
                        axis,
                        self.inputs.l0_budget,
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
                self._absorb_warnings(self._deferral_warnings(result.get("stranded", [])))
                self.tier3 = result
                for member in result["front"]:
                    # Every front member is recorded, not just one. The archive
                    # is what CF-ORCH-007 reports from, and handing it a single
                    # representative would discard the trade-off tier 3 exists
                    # to produce before anything downstream could see it.
                    self.record_candidate(
                        "mutation", member["cds"], objectives=tuple(member["objectives"])
                    )
                self.persist_checkpoint()

                if result["terminated"]:
                    # Real termination, not the end of a segment. Continuing
                    # would spend budget on a search that has stopped moving.
                    break

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

    def best_candidate(self):
        """The candidate this job reports, warning if a later stage regressed.

        The warning is derived rather than narrated: it names the origin that
        actually recorded last and the one actually being reported. The first
        version hardcoded the stage as `tier3` and explained the generation
        loop by name, which made it a statement about this round's pipeline
        rather than about whatever pipeline ran.
        """
        best = self.archive.best()
        if best is not None and self.archive.regressed():
            last = self.archive.last_recorded()
            assert last is not None
            self._add_warning(
                {
                    "type": "https://codonforge.org/problems/cf-orch-007",
                    "title": "A later stage did not improve on an earlier one",
                    "status": 200,
                    "code": "CF-ORCH-007",
                    "detail": (
                        f"the reported candidate came from {best.origin}, which scored "
                        f"{best.score:.6g}; the last stage to produce one was "
                        f"{last.origin}, at {last.score:.6g}. Stages are not required to "
                        "improve on each other - exploration is not monotone - so the "
                        "job reports the best candidate any stage produced rather than "
                        "the most recent one"
                    ),
                    "stage": "tier_plan",
                    "subject": last.origin,
                }
            )
        return best

    def _add_warning(self, warning: dict[str, Any]) -> None:
        """Append a warning unless one with the same code and subject is there."""
        key = (warning.get("code"), warning.get("subject"))
        if key not in {(w.get("code"), w.get("subject")) for w in self.warnings}:
            self.warnings.append(warning)

    def record_candidate(
        self, origin: str, cds: str, objectives: tuple[float | None, ...] = ()
    ) -> None:
        """Add a candidate to the archive, scored and feasibility-checked.

        Both are computed here rather than taken from the producing stage. A
        stage reporting its own score would be reporting the surrogate it
        optimized, which is not always the objective the front claims; and a
        stage asserting its own feasibility is the kind of unchecked claim
        `verify_capability` exists to refuse elsewhere.
        """
        probe = self._tier_request(cds)
        report = json.loads(_core.scalarized_score_report(probe, cds))
        self.archive.record(
            origin=origin,
            cds=cds,
            score=report["score"],
            feasible=self.is_feasible(cds),
            digest=self.sequence_digest(cds),
            objectives=objectives,
        )

    def _deferral_warnings(self, stranded: list[int]) -> list[dict[str, Any]]:
        """Report individuals the L0 evaluation budget never reached.

        `CF-SRCH-016`, P2-06. Rationing defers rather than discards, and the
        deferral is modelled correctly in `cf_search::rationing` — but modelling
        it correctly is not the same as making it visible. An individual that
        stays `UNKNOWN` to the end because `E_max` never reached it leaves a job
        that searched less than its budget plan describes, with nothing in the
        output saying so. That is structurally the same problem `CF-ORCH-006`
        exists to solve for skipped tiers, so it gets the same treatment: a
        warning carrying (code, stage, subject) provenance.

        Empty in this build, because no L0 term exists yet to ration - the
        structure kernel is blocked on `DEP-STRU-001` and the LM is Phase 3. The
        path is wired and tested now rather than after the first L0 term lands,
        since the failure it guards against is silence, and silence is exactly
        what an unwired path produces.
        """
        if not stranded:
            return []
        return [
            Warning_(
                code="CF-SRCH-016",
                title="Some individuals were never exactly evaluated",
                detail=(
                    f"{len(stranded)} individual(s) were deferred by the L0 evaluation "
                    f"budget (E_max) and reached the end of the run with UNKNOWN "
                    f"objectives. They were not discarded and not imputed, but they were "
                    f"also never scored, so the front was selected over less of the "
                    f"population than the plan describes"
                ),
                stage="tier3_moea",
                subject=",".join(str(index) for index in sorted(stranded)[:16]),
            ).to_problem_details()
        ]

    def is_feasible(self, cds: str) -> bool:
        """Whether a sequence has any *actionable* blocking violation.

        Actionable, per `DEP-CONS-002`: a violation confined to immutable
        sequence is a property of the input rather than a fault in the
        candidate, and every other layer draws the line the same way.
        """
        construct = self.inputs.five_prime_flank + cds + self.inputs.three_prime_flank
        start = len(self.inputs.five_prime_flank)
        report = json.loads(
            _core.constraint_violations(
                json.dumps(self.inputs.constraint_set),
                construct,
                [(start, start + len(cds))],
                "assembled",
            )
        )
        return not any(v["actionable"] for v in report["violations"])

    def sequence_digest(self, cds: str) -> str:
        return json.loads(
            _core.sequence_document(
                cds,
                self.inputs.translation_table,
                self.inputs.five_prime_flank,
                self.inputs.three_prime_flank,
            )
        )["digest"]

    def _l0_scorer(self):
        """The out-of-core scorer for this job, and which axis it fills.

        `P2-11`. Returns `(None, 0)` when every objective is in-core, which is
        every job that does not ask for a structure term.

        The axis is the count of in-core objectives, not a configured number:
        tier 3 appends the L0 axis after the evaluator's own (`DEP-SRCH-012`),
        so the index follows from the split rather than being maintained
        alongside it.
        """
        # Imported here rather than at module scope: `structure_diag` reaches
        # the predictors package, which reaches an optional engine, and a job
        # that asks for no structure objective should not pay for that import
        # or fail on it.
        from codonforge.orchestrator import structure_diag

        out_of_core = [t for t in self.plan.scored if t == structure_diag.TERM_ID]
        if not out_of_core:
            return None, 0
        axis = len([t for t in self.plan.scored if t != structure_diag.TERM_ID])

        def score(sequences: list[str]) -> list[float | None]:
            # One call per generation. The engine is stateless, so this holds no
            # cache across calls - front members differ at many positions and a
            # cache is exact only for edits it was told about.
            return structure_diag.evaluate(sequences)

        return score, axis

    def _tier_request(self, cds: str, terms: list[str] | None = None) -> str:
        """The shared document tiers 1 and 2 both take."""
        return json.dumps(
            {
                "cds": cds,
                "translation_table": self.inputs.translation_table,
                "five_prime_flank": self.inputs.five_prime_flank,
                "three_prime_flank": self.inputs.three_prime_flank,
                "constraint_set": self.inputs.constraint_set,
                "host_profile": self.inputs.host_profile,
                "host_id": self.inputs.host_id,
                "host_digest": self.inputs.host_digest,
                # Defaults to every scored term. Tiers 1 and 2 pass their own
                # narrower list: an out-of-core term has no delta table, so
                # handing it to tier 1 is a CF-LAND-002 rather than a slow
                # path. The plan already knows which tier scores what; this
                # just stops the request from over-reporting it.
                "terms": terms if terms is not None else self.plan.scored,
                "weights": self._weights(terms),
                "root_seed": self.root_seed(),
            }
        )

    def _weights(self, terms: list[str] | None = None) -> list[float]:
        """Scalarization coefficients. Magnitudes, deliberately unsigned.

        Tier 1 and tier 2 each read `provider.direction()` and apply the sense
        multiplier themselves, so the caller supplies magnitudes only. Signing
        them here as well would apply the direction twice and invert the
        objective — the search would minimize what the documents said it
        maximized, and every document would still read correctly, which is
        exactly the silent inversion `direction` is carried to prevent.

        This is not hypothetical: this method did sign them, and a
        `minimize` term came out maximized end to end. `TST-CLI-004` pins it.
        """
        supplied = self.inputs.term_weights or [1.0] * len(self.plan.scored)
        if any(weight < 0 for weight in supplied):
            raise ValueError(
                "scalarization weights are magnitudes; the tiers apply each term's "
                "declared direction themselves, so a negative weight applies it twice"
            )
        if terms is None:
            return list(supplied)
        # A tier scoring a subset needs the weights for that subset, aligned by
        # position. Zipping a full weight list against a shorter term list would
        # silently give each term its neighbour's weight.
        by_term = dict(zip(self.plan.scored, supplied, strict=False))
        return [by_term.get(term, 1.0) for term in terms]

    def _run_seeding_tiers(self) -> None:
        """Run tier 1 and tier 2, in plan order, if the budget includes them."""
        stages = self.budget.stages()
        cds = self.inputs.cds
        self.record_candidate("input", cds)

        if "tier1_seed" in stages:
            report = json.loads(
                _core.search_tier1_seed(
                    self._tier_request(cds, self.plan.tier1),
                    self.inputs.tier1_beam,
                    self.inputs.tier1_seeds,
                )
            )
            if report["seeds"]:
                # Best surrogate first, and tier 1 already returns them in that
                # order. Taking [0] rather than re-sorting avoids introducing a
                # tie-break this layer would have to define.
                cds = _cds_of(report["seeds"][0])
                self.record_candidate("tier1_dp", cds)
            self.stage_reports["tier1_seed"] = report
            self._absorb_warnings(report.get("warnings", []))

        if "tier2_block_ca" in stages:
            report = json.loads(
                _core.search_tier2_polish(
                    self._tier_request(cds, self.plan.tier2),
                    self.inputs.tier2_block,
                    self.inputs.tier2_sweeps,
                )
            )
            cds = _cds_of(report["incumbent"])
            self.record_candidate("tier2_block_ca", cds)
            self.stage_reports["tier2_block_ca"] = report
            self._absorb_warnings(report.get("warnings", []))

        self.seeded_cds = cds

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
                    # The lockfile itself, not just its digest. A resumed job
                    # that re-resolved would pick up any artifact promoted
                    # while it was stopped - which is the exact thing being
                    # checkpointed makes possible, and the exact thing pinning
                    # at submission exists to prevent. Carrying the sealed
                    # document means resume replays digests instead.
                    "lockfile": self.lockfile.to_json(),
                    "generations_run": self.generations_run,
                    "budget_consumed": self.budget.consumed,
                    # What the seeding tiers did, and what they produced. A
                    # resumed job does not re-run them, so recomputing
                    # provenance from what this invocation happened to execute
                    # would report a candidate tier 2 produced as having come
                    # straight from the input.
                    "stage_reports": self.stage_reports,
                    "archive": self.archive.to_json(),
                    "seeded_cds": self.seeded_cds,
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
        self.stage_reports = checkpoint_document.get("stage_reports", {})
        self.archive = CandidateArchive.from_json(checkpoint_document.get("archive", []))
        self.seeded_cds = checkpoint_document.get("seeded_cds")
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
        summary = self.front_summary()
        if summary is not None:
            document["current_front_summary"] = summary
        return document

    def front_summary(self) -> dict[str, Any] | None:
        """`job_status.current_front_summary`, from tier 3's own report.

        Filled from the population loop rather than recomputed, so the number a
        UI plots is the number termination was decided on. `best_per_objective`
        is the best *value* attained on each axis — a property of the front,
        not a selected candidate. It is deliberately not a third notion of
        "best" alongside `front()` and `best()`: those two answer "which
        designs" and this answers "how far did each axis get".
        """
        if self.tier3 is None:
            return None
        front = self.tier3["front"]
        summary: dict[str, Any] = {
            "hypervolume_history": list(self.tier3["hypervolume_history"]),
            "n_nondominated": len(front),
            "n_feasible": sum(1 for m in front if m["feasible"]),
        }
        if summary["hypervolume_history"]:
            summary["hypervolume"] = summary["hypervolume_history"][-1]
        if front:
            summary["best_per_objective"] = {
                term: min(
                    m["objectives"][index] for m in front if m["objectives"][index] is not None
                )
                for index, term in enumerate(self.plan.scored)
                if any(m["objectives"][index] is not None for m in front)
            }
        return summary

    def _progress(self) -> dict[str, Any]:
        """Per-tier progress: one entry per planned tier, each reporting itself.

        Structured rather than a flat pair of counters because the tiers are
        not interchangeable. A job that ran tier 1 and skipped tier 2 and one
        that ran both are different outcomes, and a single percentage cannot
        say which happened.

        Only `tier3_moea` carries generation counters, because only it counts
        in generations. Tier 1 runs once and tier 2 runs sweeps, and reporting
        the loop's generation count against them would state something untrue
        about work they did not do — the mirror image of flattening.
        """
        running = self.state == "running"
        tiers = []
        for tier in self.budget.stages():
            if tier not in EXECUTED_TIERS:
                tiers.append({"tier": tier, "state": "skipped"})
                continue

            if tier in self.stage_reports:
                tiers.append({"tier": tier, "state": "done", "fraction": 1.0})
                continue

            if tier != GENERATIONAL_TIER:
                # Planned, executable, and not yet reached.
                tiers.append({"tier": tier, "state": "running" if running else "pending"})
                continue

            if self.state == "succeeded":
                state = "done"
            elif running:
                state = "running"
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
            unfinished = [t["tier"] for t in tiers if t["state"] in ("running", "pending")]
            current = unfinished[0] if unfinished else "none"
        return {"current_tier": current, "tiers": tiers}

    def warnings_from(self, stage: str) -> list[dict[str, Any]]:
        """Warnings a given stage produced. `DEP-ORCH-001`'s whole point."""
        return [w for w in self.warnings if w.get("stage") == stage]


def _cds_of(document: dict[str, Any]) -> str:
    """The coding sequence out of a `sequence.schema.json` document."""
    return str(document["cds"])


def _now() -> str:
    return _dt.datetime.now(_dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
