"""The library surface: `cf.design(...)`.

Ticket P1-30. Specification: `planning/04-api/02-cli-and-library.md` §2.

Lives in `cf-orchestrator` rather than `cf-cli` because it needs `cf-decision`
to export a front, and `cf-cli` is declared as "a thin argument-parsing shell
over cf-orchestrator" that does not depend on `cf-decision`. Putting it here
keeps the CLI thin in fact and not just in description, and `codonforge.design`
re-exports it so the documented `cf.design(...)` spelling is the real one.

One implementation, three faces. The CLI builds a `DesignRequest` from its
flags and calls straight into here; there is no second validation path and no
second set of defaults. `TST-CLI-001` asserts the equivalence, but equivalence
between surfaces that share an implementation is nearly free — the load-bearing
test is `TST-CLI-003`, which asserts the pipeline *ran* what the budget planned.
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any

from codonforge import _core
from codonforge.decision import front as front_module
from codonforge.modelreg import ArtifactStore, Lockfile, ResolutionRequest, resolve, seal

#: Terms this build ships providers for.
#:
#: Imported rather than restated. There were three copies of this list - here,
#: in `modelreg.resolve`, and in the Rust registry - and the round that added
#: `structure.local_accessibility` updated one, which is how a term ended up
#: advertised by the CLI and rejected by the resolver. Two remain, because the
#: Rust one cannot be imported from Python; `TST-CLI-010` holds those together.
from codonforge.modelreg.resolve import BUILTIN_TERMS
from codonforge.orchestrator import structure_diag

from .budget import ComputeBudget
from .job import Job, JobInputs


class DesignError(Exception):
    """A design request that cannot be honoured, carrying its `CF-*` code."""

    def __init__(self, code: str, detail: str, exit_code: int = 1) -> None:
        self.code = code
        self.detail = detail
        self.exit_code = exit_code
        super().__init__(f"{code}: {detail}")


@dataclass
class DesignRequest:
    """The one request object all three surfaces build.

    The same shape the REST body validates against and the CLI builds from its
    flags. When surfaces drift it is always because someone added a code path;
    a single request type makes that structurally hard rather than discouraged.
    """

    host: str
    cds: str
    host_profile: dict[str, Any]
    objectives: list[str] = field(default_factory=lambda: ["usage.cai"])
    constraints: dict[str, Any] = field(default_factory=dict)
    five_prime_flank: str = ""
    three_prime_flank: str = ""
    translation_table: int = 11
    budget: str = "standard"
    seed: int = 0
    generations: int = 12
    segment: int = 4
    weights: list[float] = field(default_factory=list)
    host_digest: str = "sha256:" + "0" * 64
    job_id: str | None = None
    checkpoint_dir: str | None = None
    resume: bool = False

    def to_json(self) -> str:
        return json.dumps(asdict(self), sort_keys=True)

    #: Fields that describe *where output goes*, not what is computed.
    #:
    #: Excluded from the request's identity. A job told where to write its
    #: checkpoints resolves the same landscape as one that was not, and
    #: `--resume` must land on the same job id as the run it resumes or it
    #: cannot find the checkpoint. Deriving the id from the whole request made
    #: both false: two identical designs got different lockfiles because one
    #: had `--checkpoint-dir`.
    PLUMBING = ("job_id", "checkpoint_dir", "resume")

    def identity(self) -> str:
        """The canonical form of everything that determines the result."""
        fields = {k: v for k, v in asdict(self).items() if k not in self.PLUMBING}
        return json.dumps(fields, sort_keys=True)

    @classmethod
    def from_json(cls, text: str) -> DesignRequest:
        return cls(**json.loads(text))


@dataclass
class DesignResult:
    """What a completed (or cancelled) design returns."""

    job: Job
    front: dict[str, Any]

    @property
    def status(self) -> dict[str, Any]:
        return self.job.status()

    @property
    def lockfile(self) -> Lockfile:
        return self.job.lockfile

    @property
    def warnings(self) -> list[dict[str, Any]]:
        return self.job.status()["warnings"]

    def sequence_document(self) -> dict[str, Any]:
        """The `sequence.schema.json` document for the reported candidate.

        Built by the job, so the digest a writer round-trips against is the
        same one the archive selected on.
        """
        cds = self.front["candidates"][0]["sequences"][0]["cds"]
        return json.loads(
            _core.sequence_document(
                cds,
                self.job.inputs.translation_table,
                self.job.inputs.five_prime_flank,
                self.job.inputs.three_prime_flank,
            )
        )

    def best(self) -> dict[str, Any]:
        if not self.front["candidates"]:
            raise DesignError("CF-CLI-002", "the front is empty", exit_code=1)
        return self.front["candidates"][0]

    def to_json(self, indent: int | None = 2) -> str:
        return json.dumps(
            {"status": self.status, "front": self.front}, indent=indent, sort_keys=True
        )


def score_sequences(request: DesignRequest, sequences: list[str]) -> list[dict[str, Any]]:
    """Score existing CDS strings on this request's landscape.

    The same `_directions_request` a job builds, with only the CDS swapped.
    `/v1/score` and the benchmark harness both call this, so a competitor
    sequence and a user's sequence are not scored by two implementations.
    """
    scored: list[dict[str, Any]] = []
    for cds in sequences:
        document = json.loads(
            _core.term_scores(_directions_request(replace(request, cds=cds)), cds)
        )
        scored.append(document)
    return scored


def assigned_job_id(request: DesignRequest) -> str:
    """The id `design` will use: the caller's, or the one derived from the request."""
    return request.job_id or _derived_id("job", request.identity())


def design(
    request: DesignRequest,
    store: ArtifactStore,
    *,
    background: bool = False,
    on_segment: Any = None,
) -> DesignResult:
    """Run a design job to completion, in process.

    Synchronous, because a researcher in a notebook wants a blocking call and
    because `P1-29` executes in process with no daemon.
    """
    if background:
        # DEP-CLI-004. Returning a fake handle that had already finished, or
        # running synchronously while claiming otherwise, would make every
        # timing assumption in the calling code wrong.
        raise DesignError(
            "CF-CLI-003",
            "background=True needs the worker pool cf-orchestrator is declared to own "
            "and P1-29 did not build; Phase 1 runs synchronously only",
            exit_code=2,
        )

    unknown = [t for t in request.objectives if t not in BUILTIN_TERMS]
    if unknown:
        raise DesignError(
            "CF-LAND-001",
            f"no built-in provider supplies {', '.join(unknown)}; this build ships "
            f"{', '.join(BUILTIN_TERMS)}",
            exit_code=2,
        )

    job_id = request.job_id or _derived_id("job", request.identity())

    # On resume the lockfile comes out of the checkpoint, never from a fresh
    # resolution. Re-sealing would re-resolve, and a resumed job that
    # re-resolved would observe any artifact promoted while it was stopped -
    # reintroducing at this layer the hazard P1-29 closed at the one below.
    resumed_document = _resume_document(request, job_id) if request.resume else None
    if resumed_document is not None and "lockfile" in resumed_document:
        lockfile = Lockfile.from_json(resumed_document["lockfile"])
    else:
        lockfile = _seal(request, store)

    # Directions come from the providers, through the same registry the tiers
    # build. Never defaulted: a term declared `minimize` and reported as
    # maximized would invert the meaning of the document a user reads.
    directions = json.loads(_core.term_directions(_directions_request(request)))

    inputs = JobInputs(
        cds=request.cds,
        translation_table=request.translation_table,
        five_prime_flank=request.five_prime_flank,
        three_prime_flank=request.three_prime_flank,
        constraint_set=request.constraints,
        tabu_leak=0.0,
        host_profile=request.host_profile,
        host_id=request.host,
        host_digest=request.host_digest,
        term_weights=list(request.weights) or [1.0] * len(request.objectives),
    )

    job = Job.submit(
        job_id=job_id,
        inputs=inputs,
        budget=ComputeBudget(
            tier=request.budget,
            total_generations=request.generations,
            segment=request.segment,
        ),
        lockfile=lockfile,
        store=store,
        terms=list(request.objectives),
        # From the registry, not asserted. TierPlan's own docstring says a
        # capability a caller asserts is not a capability; this call site
        # asserted L1 for everything, which was invisible while every built-in
        # was L1 and wrong the moment one was not.
        capabilities=json.loads(_core.term_capabilities(_directions_request(request))),
        directions=directions,
        # From the registry's own registrations, so a term that cannot be
        # ranked on is demoted here rather than silently optimized.
        roles=json.loads(_core.term_roles(_directions_request(request))),
        checkpoint_dir=Path(request.checkpoint_dir) if request.checkpoint_dir else None,
    )

    if resumed_document is not None:
        job.resume(resumed_document)

    job.run(on_segment=on_segment)
    return DesignResult(job=job, front=_front(job, request, directions))


def _resume_document(request: DesignRequest, job_id: str) -> dict[str, Any]:
    """Load the checkpoint `--resume` names, or say why it cannot."""
    if not request.checkpoint_dir:
        raise DesignError(
            "CF-CLI-001",
            "--resume needs --checkpoint-dir; there is nowhere to resume from",
            exit_code=2,
        )
    path = Path(request.checkpoint_dir) / f"{job_id}.checkpoint.json"
    if not path.exists():
        raise DesignError(
            "CF-ORCH-001",
            f"--resume was given but no checkpoint exists at {path}",
            exit_code=2,
        )
    return Job.load_checkpoint(path)


def _directions_request(request: DesignRequest) -> str:
    return json.dumps(
        {
            "cds": request.cds,
            "translation_table": request.translation_table,
            "five_prime_flank": request.five_prime_flank,
            "three_prime_flank": request.three_prime_flank,
            "constraint_set": request.constraints,
            "host_profile": request.host_profile,
            "host_id": request.host,
            "host_digest": request.host_digest,
            "terms": list(request.objectives),
            "weights": [1.0] * len(request.objectives),
            "root_seed": request.seed,
        }
    )


#: Crockford base32, the alphabet ULIDs use. I, L, O and U are excluded.
#: Y is included: the schema class `[0-9A-HJKMNP-TV-Z]` has `V-Z`. A previous
#: reading treated Y as out of alphabet; the official examples already contain
#: it (`job_01JBQ8Z4K7M2N3P5R6T7V8W9XY`). This string is the only encoder;
#: `_ID_BODY` is the schema pattern applied to the 26-character body so an
#: alphabet edit cannot emit a document the schema would reject.
_CROCKFORD = "0123456789ABCDEFGHJKMNPQRSTVWXYZ"
_ID_BODY = re.compile(r"^[0-9A-HJKMNP-TV-Z]{26}$")


def _derived_id(prefix: str, material: str) -> str:
    """A ULID-shaped id that is a function of the request, not of the clock.

    Random ids made `--resume` impossible and `TST-CLI-001` unsatisfiable at
    once: a resumed job could never match the lockfile digest it checkpointed
    under, and three surfaces running the same request would produce three
    different lockfiles by construction. Deriving them from the request means
    the same request names the same job, which is what both features assume.

    The trade is that two genuinely separate submissions of an identical
    request share an id. That is the correct reading here — they are the same
    computation, they resolve the same artifacts, and they produce the same
    answer. When jobs need to be distinguished for reasons other than their
    content, the caller passes `job_id`.
    """
    digest = _core.digest_json(json.dumps({"id": material}), None)
    value = int(digest.split(":")[1], 16)
    out = []
    for _ in range(26):
        out.append(_CROCKFORD[value % 32])
        value //= 32
    body = "".join(reversed(out))
    if not _ID_BODY.match(body):
        raise RuntimeError(
            f"derived {prefix} id body {body!r} is outside the schema Crockford "
            "class [0-9A-HJKMNP-TV-Z]{{26}}"
        )
    return f"{prefix}_{body}"


def _seal(request: DesignRequest, store: ArtifactStore) -> Lockfile:
    """Pin the landscape before any search work. `TST-ORCH-001`."""
    entries = resolve(
        store,
        ResolutionRequest(
            terms=list(request.objectives),
            host_id=request.host,
            contract_major=1,
        ),
    )
    return seal(
        entries,
        lockfile_id=_derived_id("lock", request.identity()),
        contract_major=1,
        host_profile={
            "host_id": request.host,
            "version": "1.0.0",
            "digest": request.host_digest,
        },
        inputs={
            "construct_digest": _core.digest_json(request.identity(), None),
            "constraint_set_digest": _core.digest_json(json.dumps(request.constraints), None),
            "input_sequence_digests": [_core.digest_json(json.dumps(request.cds), None)],
        },
        rng={
            "algorithm": "xoshiro256pp",
            "root_seed": request.seed,
            "stream_assignment": "per_block_index",
            "deterministic_reduction": True,
        },
        environment={
            "codonforge_version": _core.version(),
            "git_sha": "0" * 40,
            "platform": "linux-x86_64",
            "thread_count": 1,
            "budget_tier": request.budget,
        },
    )


def _front(job: Job, request: DesignRequest, directions: dict[str, str]) -> dict[str, Any]:
    """The front document for a finished job.

    Reports whatever `job.archive` selects. The selection rule — feasible over
    infeasible, then higher scalarized score, then smallest digest — lives in
    `cf-orchestrator` so that this surface, the REST surface and the output
    writers all report the same candidate rather than each deriving one.
    """
    best = job.best_candidate()
    if best is None:
        raise DesignError("CF-CLI-002", "the job produced no candidates", exit_code=1)
    cds = best.cds

    document = json.loads(
        _core.sequence_document(
            cds,
            request.translation_table,
            request.five_prime_flank,
            request.three_prime_flank,
        )
    )
    _ = document
    reason = _termination_reason(job, best)

    # Axes `best()` considered. The lists are a byproduct of the same
    # scalarization loop that produced `best.score` — not a second pass over
    # capabilities that could drift from it. `TST-CLI-012`.
    scalarization = json.loads(_core.scalarized_score_report(job._tier_request(best.cds), best.cds))

    # Every non-dominated candidate, not just the best-scored one. Emitting a
    # single design from a document named `pareto_front` was defensible while
    # the pipeline was single-objective and produced one; now that tier 3
    # produces a trade-off set, reporting one point of it would discard the
    # thing the whole tier exists to compute - and would do it in the document
    # a user actually reads. `DEP-CLI-003`, updated.
    #
    # The best-scored candidate is listed first so a caller reading
    # `candidates[0]` gets what `best()` reports and the two answers agree on
    # a single-objective job, where the front has one member anyway.
    members = job.archive.front()
    ordered = [c for c in members if c.digest == best.digest]
    ordered += [c for c in members if c.digest != best.digest]

    # DEP-STRU-003: the structure term is computed for reporting, not ranked on.
    # Details and the warning that explains them are produced together, so a
    # front cannot carry the numbers without carrying the caveat.
    # Roles come from the providers, through the same accessor `directions`
    # uses. A mapping maintained here would be a second copy of what
    # TermProvider::max_role already says.
    roles = json.loads(_core.term_roles(_directions_request(request)))

    # Only when the term is *not* already an objective. Once P2-11 made it a
    # real axis, computing it a second time here would fold the same construct
    # twice per candidate and attach CF-LAND-013 - a warning saying the term
    # took no part in a search that had just ranked on it.
    if structure_diag.TERM_ID in job.plan.scored:
        structure_details, structure_warning = None, None
    else:
        structure_details, structure_warning = structure_diag.annotate(
            [member.cds for member in ordered]
        )
        job._absorb_warnings([structure_warning])

    candidates = []
    for position, member in enumerate(ordered):
        member_document = json.loads(
            _core.sequence_document(
                member.cds,
                request.translation_table,
                request.five_prime_flank,
                request.three_prime_flank,
            )
        )
        member_scores = json.loads(_core.term_scores(_directions_request(request), member.cds))
        vector = [
            _axis_value(member.objectives, index, member_scores.get(term))
            for index, term in enumerate(job.plan.scored)
        ]
        candidates.append(
            front_module.candidate(
                candidate_id=f"cand_{member_document['digest'].split(':')[1][:16]}",
                cds=member.cds,
                sequence_digest=member_document["digest"],
                # Same `_axis_value` as `unknown_terms` and `objective_detail`.
                # Building the vector from in-core `term_scores` alone left a
                # scored L0 axis as `null` in the vector while `unknown_terms`
                # was empty and `objective_detail` carried the number — three
                # fields describing one fact, two of them lying. The schema
                # says a term in `unknown_terms` is UNKNOWN in the vector;
                # the converse has to hold too.
                objective_vector=vector,
                origin=member.origin,
                feasible=member.feasible,
                generation_found=job.generations_run,
                objective_detail=_with_roles(
                    structure_details[position] if structure_details else None,
                    member_scores,
                    job.plan.scored,
                    roles,
                    member.objectives,
                ),
                unknown_terms=[
                    term for index, term in enumerate(job.plan.scored) if vector[index] is None
                ],
            )
        )

    return front_module.build(
        job_id=job.job_id,
        front_id=_derived_id("front", f"{job.job_id}:{cds}"),
        lockfile_digest=job.lockfile.digest or "",
        objective_order=list(job.plan.scored),
        objective_sense=[directions[t] for t in job.plan.scored],
        candidates=candidates,
        budget_tier=request.budget,
        termination_reason=reason,
        n_evaluations=job.generations_run,
        input_sequence=_baseline(job, request),
        recommended={
            "candidate_id": f"cand_{best.digest.split(':')[1][:16]}",
            "rule": "feasible_then_score_then_digest",
            # CF-ORCH-007's guarantee, stated rather than implied: the reported
            # answer is the best across every stage, which is only interesting
            # to a reader when it is *not* the last stage's output.
            "differs_from_last_stage": (
                job.archive.last_recorded() is not None
                and job.archive.last_recorded().digest != best.digest
            ),
            # Same discipline CF-ORCH-006 applies to a skipped tier: when the
            # scalarization excludes an axis that is on the front, say so.
            "informed_by": list(scalarization["informed_by"]),
            "not_informed_by": list(scalarization["not_informed_by"]),
        },
    )


def _termination_reason(job: Job, best: Any) -> str:
    """The `pareto_front.summary.termination_reason`, from the run.

    Read from tier 3's own report where there is one, so the document says why
    the search stopped rather than restating why the *budget* ended.
    """
    if job.state == "cancelled":
        return "user_cancel"
    if not best.feasible:
        return "no_feasible_solution"
    if job.tier3 is not None:
        return str(job.tier3["termination"])
    return "generation_limit"


def _baseline(job: Job, request: DesignRequest) -> dict[str, Any] | None:
    """The input's own scores, from the archive rather than re-evaluated."""
    recorded = next((c for c in job.archive if c.origin == "input"), None)
    if recorded is None:
        return None
    scores = json.loads(_core.term_scores(_directions_request(request), recorded.cds))
    return {
        "digest": recorded.digest,
        "objective_vector": [scores[t] for t in job.plan.scored],
        "feasible": recorded.feasible,
    }


def _axis_value(
    objectives: tuple[float | None, ...] | None, index: int, in_core: float | None
) -> float | None:
    """The value for one axis, preferring the in-core score.

    An out-of-core term has no in-core value, so it falls back to the objective
    vector tier 3 recorded - which is where a batched score actually lands.
    Returning `None` when neither exists is the honest answer, and is what
    `unknown_terms` reports.
    """
    if in_core is not None:
        return in_core
    if objectives is not None and index < len(objectives):
        return objectives[index]
    return None


def _with_roles(
    diagnostic: dict[str, Any] | None,
    scores: dict[str, float],
    scored: list[str],
    roles: dict[str, str],
    objectives: tuple[float | None, ...] | None = None,
) -> dict[str, Any] | None:
    """Assemble `objective_detail`, every entry carrying its declared role.

    The optimized terms and the diagnostic ones end up in the same map, which is
    the point: a caller reads one place and each entry says what it was. The
    alternative the schema had until now - infer `diagnostic` from absence in
    `objective_order` - works only while absence has exactly one cause, and a
    term dropped for an unsupported host would have been indistinguishable.
    """
    detail: dict[str, Any] = {}
    for index, term in enumerate(scored):
        value = _axis_value(objectives, index, scores.get(term))
        if value is None:
            continue
        detail[term] = {
            "value": value,
            "unit": "dimensionless",
            "role": roles.get(term, "pareto_axis"),
        }
    if diagnostic:
        for term, entry in diagnostic.items():
            detail[term] = {**entry, "role": roles.get(term, "diagnostic")}
    return detail or None
