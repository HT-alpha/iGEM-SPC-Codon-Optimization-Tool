"""The library surface: `cf.design(...)`.

Ticket P1-30. Specification: `planning/04-api/02-cli-and-library.md` §2.

Lives in `cf-orchestrator` rather than `cf-cli` because it needs `cf-decision`
to export a front, and `cf-cli` is declared as "a thin argument-parsing shell
over cf-orchestrator" that does not depend on `cf-decision`. Putting it here
keeps the CLI thin in fact and not just in description, and `dna_design_suite.design`
re-exports it so the documented `cf.design(...)` spelling is the real one.

One implementation, three faces. The CLI builds a `DesignRequest` from its
flags and calls straight into here; there is no second validation path and no
second set of defaults. `TST-CLI-001` asserts the equivalence, but equivalence
between surfaces that share an implementation is nearly free — the load-bearing
test is `TST-CLI-003`, which asserts the pipeline *ran* what the budget planned.
"""

from __future__ import annotations

import json
import os
import re
import threading
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any

from dna_design_suite import _core
from dna_design_suite.codon.decision import front as front_module
from dna_design_suite.codon.orchestrator import structure_diag
from dna_design_suite.modelreg import ArtifactStore, Lockfile, ResolutionRequest, resolve, seal

#: Terms this build ships providers for.
#:
#: Imported rather than restated. There were three copies of this list - here,
#: in `modelreg.resolve`, and in the Rust registry - and the round that added
#: `structure.local_accessibility` updated one, which is how a term ended up
#: advertised by the CLI and rejected by the resolver. Two remain, because the
#: Rust one cannot be imported from Python; `TST-CLI-010` holds those together.
from dna_design_suite.modelreg.resolve import BUILTIN_TERMS, Pin

from .budget import ComputeBudget
from .job import Job, JobInputs
from .rationale import build_rationale


class DesignError(Exception):
    """A design request that cannot be honoured, carrying its `CF-*` code."""

    def __init__(self, code: str, detail: str, exit_code: int = 1) -> None:
        self.code = code
        self.detail = detail
        self.exit_code = exit_code
        super().__init__(f"{code}: {detail}")


class BackgroundDesign:
    """A handle for a design that is running on the worker pool.

    `DEP-CLI-006`. `done` is false until the pool thread finishes. `result`
    waits for that thread. The handle is not a finished job.
    """

    def __init__(self, job_id: str, future: Future) -> None:
        self.job_id = job_id
        self._future = future

    def done(self) -> bool:
        return self._future.done()

    def result(self, timeout: float | None = None) -> DesignResult:
        return self._future.result(timeout=timeout)


_POOL: ThreadPoolExecutor | None = None
_POOL_LOCK = threading.Lock()


def _worker_pool() -> ThreadPoolExecutor:
    """One pool for the process. Independent designs overlap on it."""
    global _POOL
    with _POOL_LOCK:
        if _POOL is None:
            workers = min(8, os.cpu_count() or 2)
            _POOL = ThreadPoolExecutor(max_workers=workers, thread_name_prefix="dna-design")
        return _POOL


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
    #: Hypervolume stagnation window. The engine default is 8. A campaign
    #: that is testing a generation budget sets this above that budget so
    #: the budget, not the window, is what stops the run.
    stagnation_window: int = 8
    weights: list[float] = field(default_factory=list)
    host_digest: str = "sha256:" + "0" * 64
    #: Set by `design()` when the Potts artifact is installed. Excluded from
    #: the request identity: the lockfile resolution entry is the record.
    context_digest: str | None = None
    job_id: str | None = None
    checkpoint_dir: str | None = None
    resume: bool = False
    #: Wall cap in milliseconds. `None` uses the ADR-0008 length-scaled target.
    wall_budget_ms: int | None = None
    #: Resident-set cap in bytes. `None` uses the laptop profile's 8 GiB.
    max_rss_bytes: int | None = None
    #: Synonymous CDS to search from. Empty means the supplied `cds`.
    #: Omitted from the identity when empty, so a job that does not set it
    #: keeps the lockfile it had before this field existed. `DEP-API-003`.
    starting_cds: str = ""
    #: `viennarna` (the default fold) or `nussinov` (opt-in approximate
    #: matching). Empty and `viennarna` are omitted from the identity.
    #: `DEP-STRU-006`.
    structure_engine: str = ""

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
    PLUMBING = ("job_id", "checkpoint_dir", "resume", "context_digest")

    def identity(self) -> str:
        """The canonical form of everything that determines the result."""
        fields = {k: v for k, v in asdict(self).items() if k not in self.PLUMBING}
        if not fields.get("starting_cds"):
            fields.pop("starting_cds", None)
        if fields.get("structure_engine", "") in ("", "viennarna"):
            fields.pop("structure_engine", None)
        return json.dumps(fields, sort_keys=True)

    @classmethod
    def from_json(cls, text: str) -> DesignRequest:
        payload = json.loads(text)
        payload.setdefault("starting_cds", "")
        payload.setdefault("structure_engine", "")
        return cls(**payload)


@dataclass
class DesignResult:
    """What a completed (or cancelled) design returns."""

    job: Job
    front: dict[str, Any] | None

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
        if self.front is None:
            raise DesignError(
                "CF-CLI-002", "a failed design has no sequence to report", exit_code=1
            )
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
        if self.front is None or not self.front["candidates"]:
            raise DesignError("CF-CLI-002", "the front is empty", exit_code=1)
        return self.front["candidates"][0]

    def to_json(self, indent: int | None = 2) -> str:
        return json.dumps(
            {"status": self.status, "front": self.front}, indent=indent, sort_keys=True
        )

    def to_sbol(self) -> str:
        """SBOL 3 for the recommended candidate and the request documents.

        The writer renders ``front``, ``sequence_document()``, the host
        profile, and the constraint set. It does not score the sequence again.
        A missing ``sbol3`` install is ``CF-SBOL-001``.
        """
        from dna_design_suite.export.sbol import to_sbol

        return to_sbol(
            self.front or {},
            self.sequence_document(),
            self.job.inputs.host_profile,
            self.job.inputs.constraint_set,
            design_host=self.job.inputs.host_id,
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


def _apply_starting_cds(request: DesignRequest) -> DesignRequest:
    """Use a validated synonymous CDS as the sequence the search starts from.

    The protein must match the request CDS. A difference is `CF-SEQ-004`
    with the residue positions, the same report seed import already makes.
    """
    starting = "".join((request.starting_cds or "").split()).upper().replace("U", "T")
    if not starting:
        return request
    try:
        target = json.loads(
            _core.sequence_document(request.cds, request.translation_table, None, None, None)
        )
        found = json.loads(
            _core.sequence_document(starting, request.translation_table, None, None, None)
        )
    except ValueError as exc:
        raise DesignError("CF-SEQ-004", str(exc)) from exc
    expected = str(target["protein"])
    actual = str(found["protein"])
    if expected != actual:
        diffs = []
        for index in range(max(len(expected), len(actual))):
            left = expected[index] if index < len(expected) else "-"
            right = actual[index] if index < len(actual) else "-"
            if left != right:
                diffs.append(f"{index}:{left}->{right}")
                if len(diffs) == 8:
                    break
        raise DesignError(
            "CF-SEQ-004",
            "starting CDS translates to a different protein at " + ", ".join(diffs),
        )
    return replace(request, cds=starting, starting_cds=starting)


def assigned_job_id(request: DesignRequest) -> str:
    """The id `design` will use: the caller's, or the one derived from the request."""
    return request.job_id or _derived_id("job", request.identity())


def design(
    request: DesignRequest,
    store: ArtifactStore,
    *,
    background: bool = False,
    on_segment: Any = None,
) -> DesignResult | BackgroundDesign:
    """Run a design job.

    The default call blocks until the job finishes. `background=True` returns
    a handle immediately and runs the same function on the worker pool
    (`DEP-CLI-006`). The handle's `result` is that job. It is not a finished
    document returned early.
    """
    _stamp_harmonize_native(request)
    if background:
        if request.job_id is None:
            request = replace(request, job_id=assigned_job_id(request))
        future = _worker_pool().submit(_run_design, request, store, on_segment)
        return BackgroundDesign(request.job_id or "", future)
    return _run_design(request, store, on_segment)


def _run_design(
    request: DesignRequest,
    store: ArtifactStore,
    on_segment: Any = None,
) -> DesignResult:
    """The body of `design`, on the caller's thread or a pool thread."""
    _stamp_harmonize_native(request)
    from dna_design_suite.codon.epistasis.dwell import DWELL_TERM, dwell_status

    if DWELL_TERM in request.objectives:
        dwell = dwell_status()
        # Not fitted in this build. Do not fall through into a job that
        # would omit the axis or report it as zero.
        if dwell["state"] != "fitted":
            raise DesignError("CF-LAND-001", dwell["reason"])
        raise DesignError(
            "CF-LAND-001",
            "epistasis.dwell_deviation is not scored in this build; "
            "a fitted GSE53767 head is not installed",
        )

    request = _apply_starting_cds(request)
    from dna_design_suite.codon.orchestrator import local_structure

    local_structure.validate(request.constraints)
    local_structure.require_engine(
        request.constraints,
        structure_diag.TERM_ID in request.objectives,
        request.structure_engine,
    )

    from dna_design_suite.codon.epistasis.context import TERM_ID as CONTEXT_TERM

    unknown = [t for t in request.objectives if t not in BUILTIN_TERMS and t != CONTEXT_TERM]
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
    from dna_design_suite.codon.epistasis.context import TERM_ID as CONTEXT_TERM
    from dna_design_suite.codon.epistasis.context import ContextError, install_pinned, prepare

    resumed_document = _resume_document(request, job_id) if request.resume else None
    context_window = None
    if resumed_document is not None and "lockfile" in resumed_document:
        lockfile = Lockfile.from_json(resumed_document["lockfile"])
        for entry in lockfile.resolution:
            if entry.term_id != CONTEXT_TERM:
                continue
            try:
                installed = install_pinned(store, entry.digest, entry.artifact_id)
            except ContextError as exc:
                raise DesignError(exc.code, exc.detail, exit_code=2) from exc
            request.context_digest = installed["digest"]
            context_window = int(installed["window"])
            break
    else:
        context_pin = None
        if CONTEXT_TERM in request.objectives:
            try:
                installed = prepare(store, request.host)
            except ContextError as exc:
                raise DesignError(exc.code, exc.detail, exit_code=2) from exc
            request.context_digest = installed["digest"]
            context_window = int(installed["window"])
            context_pin = Pin(installed["artifact_id"], installed["digest"])
        lockfile = _seal(request, store, context_pin)

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
        stagnation_window=request.stagnation_window,
        tier2_block=max(8, context_window or 0),
        context_digest=request.context_digest,
        wall_budget_ms=request.wall_budget_ms,
        max_rss_bytes=request.max_rss_bytes,
        structure_engine=request.structure_engine,
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
        radii=json.loads(_core.term_radii(_directions_request(request))),
        checkpoint_dir=Path(request.checkpoint_dir) if request.checkpoint_dir else None,
    )

    if resumed_document is not None:
        job.resume(resumed_document)

    job.run(on_segment=on_segment)
    # A failed job's archive holds the input, recorded before the search
    # refused. Building a front from it would present that input as the design.
    if job.state == "failed":
        return DesignResult(job=job, front=None)
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


def _stamp_harmonize_native(request: DesignRequest) -> None:
    """Pin the input CDS as the harmonization source.

    `DEP-LAND-006`. The score compares each candidate codon with the codon
    that was in the input, under the source organism's relative adaptiveness.
    Rebuilding that vector from the candidate would score every sequence as
    zero whenever the source table and the host table are the same object.
    The stamp is part of the host profile, so it is part of the job identity.
    A job that did not ask for the term is left untouched.
    """
    if "usage.harmonize_rca" not in request.objectives:
        return
    profile = dict(request.host_profile)
    source = profile.get("source_codon_usage")
    if not isinstance(source, dict) or not isinstance(source.get("table"), dict):
        raise DesignError(
            "CF-HOST-003",
            "usage.harmonize_rca needs host_profile.source_codon_usage with a table "
            "of relative adaptiveness for the organism the input CDS comes from. "
            "codon_usage on the same profile is the expression host. The method is "
            "the mean absolute difference of relative codon adaptiveness (Claassens "
            "et al., PLoS ONE 2017, doi:10.1371/journal.pone.0184355). This build "
            "does not invent a source table.",
            exit_code=2,
        )
    native = "".join(request.cds.upper().split()).replace("U", "T")
    profile["harmonize_native_cds"] = native
    request.host_profile = profile


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
            "context_digest": request.context_digest,
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


def _seal(request: DesignRequest, store: ArtifactStore, context_pin: Pin | None = None) -> Lockfile:
    """Pin the landscape before any search work. `TST-ORCH-001`."""
    from dna_design_suite.codon.epistasis.context import TERM_ID as CONTEXT_TERM

    pins = {CONTEXT_TERM: context_pin} if context_pin is not None else None
    entries = resolve(
        store,
        ResolutionRequest(
            terms=list(request.objectives),
            host_id=request.host,
            contract_major=1,
            pins=pins,
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
    members = job.archive.front()
    chosen, rule, informed, excluded = _recommendation(job, best, members)
    cds = chosen.cds

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

    # Every non-dominated candidate, not just the best-scored one. Emitting a
    # single design from a document named `pareto_front` was defensible while
    # the pipeline was single-objective and produced one; now that tier 3
    # produces a trade-off set, reporting one point of it would discard the
    # thing the whole tier exists to compute - and would do it in the document
    # a user actually reads. `DEP-CLI-003`, updated.
    #
    # `candidates[0]` is the recommended candidate. On one objective that is
    # `best()`. On two or more complete axes it is the knee (`DEP-CLI-005`).
    ordered = [c for c in members if c.digest == chosen.digest]
    ordered += [c for c in members if c.digest != chosen.digest]
    if not ordered:
        ordered = [chosen]

    # DEP-STRU-003: the structure term is computed for reporting, not ranked on.
    # Details and the warning that explains them are produced together, so a
    # front cannot carry the numbers without carrying the caveat.
    # Roles come from the providers, through the same accessor `directions`
    # uses. A mapping maintained here would be a second copy of what
    # TermProvider::max_role already says.
    roles = json.loads(_core.term_roles(_directions_request(request)))

    # CAI-only fronts do not fold. The old diagnostic called ViennaRNA once
    # per candidate so the report could show accessibility the search had
    # not ranked on. That fold does not change the CDS, and a batch of
    # CAI-only designs was waiting on it. When structure is an objective it
    # is scored by the generation loop, or by the short-construct polish in
    # tiers 1 and 2. `DEP-STRU-005`.
    structure_details = None

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
                rationale=(
                    build_rationale(
                        request,
                        member.cds,
                        _directions_request(request),
                        member_scores,
                    )
                    if position == 0
                    else None
                ),
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
            "candidate_id": _candidate_id(chosen.digest),
            "rule": rule,
            # True when the named answer is not the last stage's output.
            "differs_from_last_stage": (
                job.archive.last_recorded() is not None
                and job.archive.last_recorded().digest != chosen.digest
            ),
            "informed_by": informed,
            "not_informed_by": excluded,
        },
    )


def _candidate_id(digest: str) -> str:
    return f"cand_{digest.split(':')[1][:16]}"


def _recommendation(
    job: Job, best: Any, members: list[Any]
) -> tuple[Any, str, list[str], list[str]]:
    """Who the front names, and which axes that choice used.

    One objective keeps `feasible_then_score_then_digest`. `archive.best()` is
    that rule, and `CF-ORCH-007` still compares stages with it. Two or more
    objectives with at least one complete feasible vector use the knee: the
    product of the remaining gap to the worst value on each minimized axis,
    ties broken by the smaller candidate id. The scalar `best()` skips an
    out-of-core axis, so using it here would put structure in
    `not_informed_by` on a job that ranked on structure. `DEP-CLI-005`.
    """
    scored = list(job.plan.scored)
    complete = [
        candidate
        for candidate in members
        if candidate.feasible
        and len(candidate.objectives) == len(scored)
        and all(value is not None for value in candidate.objectives)
    ]
    if len(scored) >= 2 and complete:
        chosen = _knee(complete)
        informed = [
            term for index, term in enumerate(scored) if chosen.objectives[index] is not None
        ]
        excluded = [term for term in scored if term not in informed]
        return chosen, "knee_of_complete_front", informed, excluded

    report = json.loads(_core.scalarized_score_report(job._tier_request(best.cds), best.cds))
    return (
        best,
        "feasible_then_score_then_digest",
        list(report["informed_by"]),
        list(report["not_informed_by"]),
    )


def _knee(candidates: list[Any]) -> Any:
    width = len(candidates[0].objectives)
    bounds = []
    for axis in range(width):
        column = [candidate.objectives[axis] for candidate in candidates]
        bounds.append((min(column), max(column)))
    winner = None
    winner_gap = None
    winner_id = None
    for candidate in candidates:
        gap = 1.0
        for axis, (low, high) in enumerate(bounds):
            span = high - low
            normalized = 0.0 if span == 0.0 else (candidate.objectives[axis] - low) / span
            gap *= max(0.0, 1.0 - normalized)
        cid = _candidate_id(candidate.digest)
        if winner is None or gap > winner_gap or (gap == winner_gap and cid < winner_id):
            winner = candidate
            winner_gap = gap
            winner_id = cid
    return winner


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
