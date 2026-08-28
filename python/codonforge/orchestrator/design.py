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
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from codonforge import _core
from codonforge.decision import front as front_module
from codonforge.modelreg import ArtifactStore, Lockfile, ResolutionRequest, resolve, seal

from .budget import ComputeBudget
from .job import Job, JobInputs

#: Terms Phase 1 ships providers for.
BUILTIN_TERMS = ("usage.cai", "usage.gc_content", "constraint.violation_count")


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

    def best(self) -> dict[str, Any]:
        if not self.front["candidates"]:
            raise DesignError("CF-CLI-002", "the front is empty", exit_code=1)
        return self.front["candidates"][0]

    def to_json(self, indent: int | None = 2) -> str:
        return json.dumps(
            {"status": self.status, "front": self.front}, indent=indent, sort_keys=True
        )


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
            f"no built-in provider supplies {', '.join(unknown)}; Phase 1 ships "
            f"{', '.join(BUILTIN_TERMS)}",
            exit_code=2,
        )

    job_id = request.job_id or _derived_id("job", request.to_json())

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
        capabilities={term: "L1" for term in request.objectives},
        radii={term: 0 for term in request.objectives},
        directions=directions,
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


#: Crockford base32, the alphabet ULIDs use: no I, L, O or U.
_CROCKFORD = "0123456789ABCDEFGHJKMNPQRSTVWXYZ"


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
    return f"{prefix}_" + "".join(reversed(out))


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
        lockfile_id=_derived_id("lock", request.to_json()),
        contract_major=1,
        host_profile={
            "host_id": request.host,
            "version": "1.0.0",
            "digest": request.host_digest,
        },
        inputs={
            "construct_digest": _core.digest_json(request.to_json(), None),
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
    """The front document for a finished job."""
    cds = job.seeded_cds or request.cds
    if job.checkpoint is not None:
        cds = job.checkpoint["incumbent"]["cds"]

    document = json.loads(
        _core.sequence_document(
            cds,
            request.translation_table,
            request.five_prime_flank,
            request.three_prime_flank,
        )
    )
    scores = json.loads(_core.term_scores(_directions_request(request), cds))

    origin = "tier2_block_ca" if "tier2_block_ca" in job.stage_reports else "input"
    reason = "user_cancel" if job.state == "cancelled" else "generation_limit"

    return front_module.build(
        job_id=job.job_id,
        front_id=_derived_id("front", f"{job.job_id}:{cds}"),
        lockfile_digest=job.lockfile.digest or "",
        objective_order=list(job.plan.scored),
        objective_sense=[directions[t] for t in job.plan.scored],
        candidates=[
            front_module.candidate(
                candidate_id=f"cand_{document['digest'].split(':')[1][:16]}",
                cds=cds,
                sequence_digest=document["digest"],
                objective_vector=[scores[t] for t in job.plan.scored],
                origin=origin,
                generation_found=job.generations_run,
            )
        ],
        budget_tier=request.budget,
        termination_reason=reason,
        n_evaluations=job.generations_run,
    )
