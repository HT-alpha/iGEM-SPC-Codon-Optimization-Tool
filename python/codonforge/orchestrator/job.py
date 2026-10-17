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
    #: Digest of the Potts artifact, when the context term is in the job.
    context_digest: str | None = None


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
    #: Which engine call `run` is inside. The single `ValueError` handler
    #: reads it so a problem document is attributed without a per-tier copy
    #: of the lifting logic. `DEP-ORCH-004`.
    _engine_stage: str = "job"
    #: Folds already measured for this job. The engine is deterministic, so a
    #: repeated sequence is the same number. Not checkpointed: the archive
    #: stores the objective values the fold produced.
    _structure_cache: dict[tuple[Any, ...], float] = field(default_factory=dict)

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

        plan = TierPlan.plan(budget, terms, capabilities, directions, roles)
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

            self._engine_stage = "tier3"
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
                    # The engine stores in-core axes and then the appended L0
                    # axis (`DEP-SRCH-012`). The archive, the knee and
                    # `best_per_objective` read `plan.scored`. Rewriting the
                    # front here, and not the checkpoint, keeps a resume on
                    # the engine's own order.
                    member["objectives"] = list(
                        self._to_plan_order(tuple(member["objectives"]))
                    )
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
        except ValueError as exc:
            # `DEP-ORCH-004`: one lift for every engine call inside this try —
            # tier 1, tier 2, and the tier-3 L0 scorer. A second catch beside
            # `search_tier3_segment` would be another copy of `DEP-ORCH-003`.
            # A `ValueError` that is not a `CF-*` problem document is the
            # original exception, and it still crashes.
            refusal = _engine_refusal(exc, self._engine_stage)
            if refusal is exc:
                raise
            self.state = "failed"
            self.finished_at = _now()
            self.problem = refusal.to_problem_details()
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
            objectives=self._objectives_for(cds, objectives),
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
        if any(v["actionable"] for v in report["violations"]):
            return False
        from codonforge.orchestrator import local_structure

        if not local_structure.rules(self.inputs.constraint_set):
            return True
        construct = self.inputs.five_prime_flank + cds + self.inputs.three_prime_flank
        return (
            local_structure.blocking_distance(
                construct, self.inputs.constraint_set, self._structure_cache
            )
            <= 0.0
        )

    def _objectives_for(
        self, cds: str, supplied: tuple[float | None, ...]
    ) -> tuple[float | None, ...]:
        """Minimized objectives, one slot per scored term.

        Tier 3 already returns this width, including an L0 axis it filled.
        Tiers 1 and 2 record a scalar. Leaving those vectors empty makes a
        length mismatch, and two candidates of different lengths never
        dominate each other, so the front mixes an unscored point with the
        population. Filling the missing slots here is the same sign convention
        tier 3 uses: minimize stays, maximize is negated.
        """
        width = len(self.plan.scored)
        if len(supplied) == width and all(value is not None for value in supplied):
            return tuple(supplied)
        measured = self._minimized_objectives(cds)
        if len(supplied) != width:
            return measured
        return tuple(
            old if old is not None else new for old, new in zip(supplied, measured, strict=True)
        )

    def _minimized_objectives(self, cds: str) -> tuple[float | None, ...]:
        scores = json.loads(_core.term_scores(self._tier_request(cds), cds))
        out: list[float | None] = []
        for term in self.plan.scored:
            raw = scores.get(term)
            # An L0 term has no in-core value. Tier 3 supplies it on the
            # candidates it scores; `propagate_objectives` copies that slot
            # onto the earlier record of the same sequence. Folding here would
            # score outside the batch path the L0 refusals are tested on.
            if raw is None:
                out.append(None)
                continue
            sign = 1.0 if self.plan.directions.get(term) == "minimize" else -1.0
            out.append(sign * float(raw))
        return tuple(out)

    def _repair_local_structure(self, cds: str) -> str:
        from codonforge.orchestrator import local_structure

        if not local_structure.rules(self.inputs.constraint_set):
            return cds

        def still_feasible(trial: str) -> bool:
            construct = self.inputs.five_prime_flank + trial + self.inputs.three_prime_flank
            start = len(self.inputs.five_prime_flank)
            report = json.loads(
                _core.constraint_violations(
                    json.dumps(self.inputs.constraint_set),
                    construct,
                    [(start, start + len(trial))],
                    "assembled",
                )
            )
            return not any(item["actionable"] for item in report["violations"])

        return local_structure.repair(
            cds,
            self.inputs.constraint_set,
            self.inputs.five_prime_flank,
            self.inputs.three_prime_flank,
            self.inputs.translation_table,
            still_feasible,
            self._structure_cache,
        )

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

    def _to_plan_order(self, supplied: tuple[float | None, ...]) -> tuple[float | None, ...]:
        """Map a tier-3 vector onto `plan.scored`.

        Tier 3's evaluator emits the in-core terms in request order and then
        appends the structure axis. `plan.scored` keeps the request order, so
        a structure term that is not last would otherwise be read as whichever
        in-core term followed it.
        """
        from codonforge.orchestrator import structure_diag

        engine = [term for term in self.plan.scored if term != structure_diag.TERM_ID]
        engine.extend(term for term in self.plan.scored if term == structure_diag.TERM_ID)
        if len(supplied) != len(engine) or engine == list(self.plan.scored):
            return supplied
        by_term = dict(zip(engine, supplied, strict=True))
        return tuple(by_term[term] for term in self.plan.scored)

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
                "context_digest": self.inputs.context_digest,
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
        self._absorb_warnings(self._unsearched_warnings(cds))

        if "tier1_seed" in stages:
            self._engine_stage = "tier1"
            try:
                report = json.loads(
                    _core.search_tier1_seed(
                        self._tier_request(cds, self.plan.tier1),
                        self.inputs.tier1_beam,
                        self.inputs.tier1_seeds,
                    )
                )
            except ValueError as exc:
                # Tier 1 ran. Recording that before the refusal leaves the
                # status is what stops a failed job from showing the tier as
                # still pending.
                refusal = _engine_refusal(exc, "tier1")
                self.stage_reports["tier1_seed"] = {
                    "seeds": [],
                    "refusal": refusal.code if isinstance(refusal, TierPlanError) else "error",
                }
                if not isinstance(refusal, TierPlanError) or refusal.code != "CF-SRCH-006":
                    raise
                # DEP-SRCH-015. The DP cannot carry GC windows or scattered
                # repeats. That refusal is the start of seed repair, not the
                # end of the job, when a synonymous sequence still exists.
                cds = self._repair_seed(cds)
            else:
                if report["seeds"]:
                    # Best surrogate first, and tier 1 already returns them in that
                    # order. Taking [0] rather than re-sorting avoids introducing a
                    # tie-break this layer would have to define.
                    cds = self._repair_local_structure(_cds_of(report["seeds"][0]))
                    self.record_candidate("tier1_dp", cds)
                self.stage_reports["tier1_seed"] = report
                self._absorb_warnings(report.get("warnings", []))

        if "tier2_block_ca" in stages:
            # No catch here. `run` lifts every `CF-*` problem document once.
            # A tier-2-only handler was how tier 3's `CF-LAND-012` stayed a
            # raw `ValueError`. `DEP-ORCH-004`.
            self._engine_stage = "tier2"
            report = json.loads(
                _core.search_tier2_polish(
                    self._tier_request(cds, self.plan.tier2),
                    self.inputs.tier2_block,
                    self.inputs.tier2_sweeps,
                )
            )
            polished = _cds_of(report["incumbent"])
            cds = self._repair_local_structure(polished)
            self.record_candidate(
                "repair" if cds != polished else "tier2_block_ca", cds
            )
            self.stage_reports["tier2_block_ca"] = report
            self._absorb_warnings(report.get("warnings", []))

        self.seeded_cds = cds

    def _unsearched_warnings(self, cds: str) -> list[dict[str, Any]]:
        """`CF-FEAS-008` for rules this version does not evaluate.

        A warning, not a failure. The empty demonstration constraint set has
        none of these rules, and a job that included `local_structure` must
        say so rather than report the rule as satisfied.
        """
        report = json.loads(_core.constraint_unsearched(self._tier_request(cds, self.plan.tier1)))
        warnings = report.get("warnings", [])
        if not isinstance(warnings, list):
            return []
        from codonforge.orchestrator import local_structure

        evaluated = local_structure.evaluated_ids(self.inputs.constraint_set)
        return [
            warning
            for warning in warnings
            if not (
                warning.get("code") == "CF-FEAS-008" and warning.get("subject") in evaluated
            )
        ]

    def _repair_seed(self, cds: str) -> str:
        """Continue from a repaired seed, or refuse with the rules still broken."""
        report = json.loads(_core.resolve_feasible_seed(self._tier_request(cds, self.plan.tier1)))
        self._absorb_warnings(report.get("unsearched", []))
        if not report.get("feasible"):
            problem = report.get("problem")
            if not isinstance(problem, dict):
                raise ValueError(
                    json.dumps({"code": "CF-SRCH-006", "detail": "seed repair failed"})
                )
            raise ValueError(json.dumps(problem))
        repaired = self._repair_local_structure(str(report["cds"]))
        self.record_candidate("repair", repaired)
        warning = report.get("warning")
        if isinstance(warning, dict):
            self._add_warning(warning)
        return repaired

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
                # Planned, executable, and not yet reached. A failed job did
                # not leave the tier waiting; it never got there.
                if self.state == "failed":
                    tiers.append({"tier": tier, "state": "skipped"})
                else:
                    tiers.append({"tier": tier, "state": "running" if running else "pending"})
                continue

            if self.state == "succeeded":
                state = "done"
            elif running:
                state = "running"
            elif self.state == "failed":
                state = "skipped"
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


def _engine_refusal(exc: ValueError, stage: str) -> Exception:
    """Lift a CF-* problem detail from the extension into `TierPlanError`.

    Programming mistakes that are not RFC 9457 documents still surface as
    crashes rather than as a forged refusal. One function serves every tier:
    copying it beside each engine call is how tier 3's `CF-LAND-012` stayed
    unstructured after tier 2 was fixed. `DEP-ORCH-003`, generalized by
    `DEP-ORCH-004`. `stage` is whichever call `run` was inside.
    """
    try:
        payload = json.loads(str(exc))
    except json.JSONDecodeError:
        return exc
    if not isinstance(payload, dict):
        return exc
    code = payload.get("code")
    if not isinstance(code, str) or not code.startswith("CF-"):
        return exc
    detail = str(payload.get("detail") or exc)
    title = str(payload.get("title") or "")
    status = payload.get("status")
    if not isinstance(status, int):
        status = 422
    return TierPlanError(
        code=code,
        detail=detail,
        subject=_subject_of(code, detail),
        stage=stage,
        title=title,
        status=status,
    )


def _subject_of(code: str, detail: str) -> str:
    """A machine-comparable subject for `DEP-ORCH-001`.

    `CF-SRCH-008` names a term; the term id is the first token of the detail.
    `CF-SRCH-009` names a configuration field.
    """
    if code == "CF-SRCH-009":
        if detail.startswith("block size"):
            return "tier2_block"
        if "weights" in detail:
            return "weights"
        return "tier2"
    if code in ("CF-LAND-011", "CF-LAND-012"):
        return "batch"
    if code.startswith("CF-STRU-"):
        return "structure.local_accessibility"
    token = detail.split(" ", 1)[0]
    return token


def _cds_of(document: dict[str, Any]) -> str:
    """The coding sequence out of a `sequence.schema.json` document."""
    return str(document["cds"])


def _now() -> str:
    return _dt.datetime.now(_dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
