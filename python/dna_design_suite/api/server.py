"""HTTP surface for the library. `P2-18`, `DEP-API-001`.

The routes are the OpenAPI paths a caller can exercise today. The body of
`POST /v1/jobs` is the library `DesignRequest`, the same object the CLI runs.
Routes whose backing module is not in this build answer `CF-CLI-003` rather
than an empty document that looks like success.

Jobs run on a thread via `design(background=False)`. `background=True` on
the library returns a worker-pool handle (`DEP-CLI-006`). Seed import stays
`CF-CLI-003`.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import threading
from dataclasses import dataclass, field
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse

from dna_design_suite import __version__
from dna_design_suite.codon.decision.reduce import present
from dna_design_suite.codon.orchestrator import (
    DesignError,
    DesignRequest,
    assigned_job_id,
    design,
    score_sequences,
)
from dna_design_suite.modelreg import ArtifactStore

_SCHEMA = "1.0.0"
# Reasons are the product. An empty 200 would look implemented.
_UNBUILT = {
    "/v1/seeds:import": "Importing an external tool's output is not wired on this route.",
}
_PRESETS = {
    "cai": ["usage.cai"],
    "cai_burden": ["usage.cai", "burden.growth_cost"],
    "cai_structure": ["usage.cai", "structure.local_accessibility"],
    "cai_context": ["usage.cai", "epistasis.context_naturalness"],
}
_HOST_FILES = (
    "mg1655_ribosomal_usage.json",
    "kazusa_e_coli_316407.json",
)
_STATIC = Path(__file__).resolve().parent / "static"
_STATIC_FILES = {
    "/": ("index.html", "text/html; charset=utf-8"),
    "/index.html": ("index.html", "text/html; charset=utf-8"),
    "/app.js": ("app.js", "text/javascript; charset=utf-8"),
    "/formats.js": ("formats.js", "text/javascript; charset=utf-8"),
    "/helix.js": ("helix.js", "text/javascript; charset=utf-8"),
    "/viewer.js": ("viewer.js", "text/javascript; charset=utf-8"),
    "/codonforge-ui.css": ("codonforge-ui.css", "text/css; charset=utf-8"),
}
# The guided page offers these as one-click presets. They are the files in
# examples/, not a second copy, so the button and the README read one table.
_EXAMPLE_FILES = {
    "/examples/mg1655_ribosomal_usage.json": (
        "mg1655_ribosomal_usage.json",
        "application/json",
    ),
    "/examples/rpsL.cds": ("rpsL.cds", "text/plain; charset=utf-8"),
    "/examples/rpsL_constraints.json": (
        "rpsL_constraints.json",
        "application/json",
    ),
    "/examples/kazusa_e_coli_316407.json": (
        "kazusa_e_coli_316407.json",
        "application/json",
    ),
    "/examples/biobrick_rfc10.json": (
        "biobrick_rfc10.json",
        "application/json",
    ),
    "/examples/golden_gate.json": (
        "golden_gate.json",
        "application/json",
    ),
    "/examples/synthesis.json": (
        "synthesis.json",
        "application/json",
    ),
    "/examples/BBa_E0040.fa": (
        "BBa_E0040.fa",
        "text/plain; charset=utf-8",
    ),
}


def problem(code: str, title: str, status: int, detail: str) -> dict[str, Any]:
    return {
        "type": f"https://dna_design_suite.org/problems/{code.lower()}",
        "title": title,
        "status": status,
        "code": code,
        "detail": detail,
    }


def _git_sha() -> str:
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "--short=12", "HEAD"],
            check=False,
            capture_output=True,
            text=True,
            timeout=2,
        )
    except (OSError, subprocess.TimeoutExpired):
        return "unavailable"
    sha = completed.stdout.strip()
    return sha or "unavailable"


@dataclass
class Record:
    tenant: str
    body_hash: str
    status: dict[str, Any]
    front: dict[str, Any] | None = None
    lockfile: dict[str, Any] | None = None
    cancel: threading.Event = field(default_factory=threading.Event)


class Api:
    def __init__(self, store: ArtifactStore) -> None:
        self.store = store
        self.records: dict[str, Record] = {}
        self.primer_records: dict[str, Any] = {}
        self.idempotency: dict[tuple[str, str], tuple[str, str]] = {}
        self.lock = threading.Lock()
        self.run_lock = threading.Lock()
        self.git_sha = _git_sha()

    def version(self) -> dict[str, Any]:
        return {
            "version": __version__,
            "git_sha": self.git_sha,
            "contract_majors": [1],
            "schema_versions": {"job_status": _SCHEMA, "pareto_front": _SCHEMA},
        }

    def health(self) -> dict[str, Any]:
        with self.lock:
            depth = sum(
                1
                for record in self.records.values()
                if record.status["state"] in {"queued", "running"}
            )
        from dna_design_suite.codon.predictors.structure import available_engines

        return {
            "status": "ok",
            "workers_ready": 0,
            "gpu_available": False,
            "queue_depth": depth,
            # Names `available_engines` actually found. An empty list means
            # ViennaRNA is not installed; the page must not invent a fold.
            "structure_engines": available_engines(),
        }

    def submit(
        self, tenant: str, body: dict[str, Any], idempotency_key: str | None
    ) -> tuple[int, dict[str, Any]]:
        prepared = _prepare_job_body(body)
        if isinstance(prepared, tuple):
            return prepared
        body = prepared
        try:
            request = DesignRequest.from_json(json.dumps(body))
        except (TypeError, json.JSONDecodeError) as exc:
            return 400, problem(
                "CF-CLI-001",
                "Command-line arguments are missing or malformed",
                400,
                "POST /v1/jobs accepts the library DesignRequest "
                f"(host, cds, host_profile, objectives). {exc}",
            )
        digest = hashlib.sha256(json.dumps(body, sort_keys=True).encode()).hexdigest()
        if idempotency_key:
            with self.lock:
                previous = self.idempotency.get((tenant, idempotency_key))
                if previous is not None:
                    previous_hash, job_id = previous
                    if previous_hash != digest:
                        return 409, problem(
                            "CF-IDEM-001",
                            "Idempotency key was reused with a different body",
                            409,
                            "the same Idempotency-Key was sent with a different request body",
                        )
                    return 200, self.records[job_id].status
        request.job_id = assigned_job_id(request)
        accepted = _queued(request.job_id)
        record = Record(
            tenant=tenant,
            body_hash=digest,
            status=dict(accepted),
        )
        with self.lock:
            self.records[request.job_id] = record
            if idempotency_key:
                self.idempotency[(tenant, idempotency_key)] = (digest, request.job_id)
        thread = threading.Thread(target=self._run, args=(record, request), daemon=True)
        thread.start()
        # The worker replaces `record.status`. The 202 body is the document
        # at accept time, so a fast start cannot change what the client was told.
        return 202, accepted

    def _run(self, record: Record, request: DesignRequest) -> None:
        with self.lock:
            record.status = {**record.status, "state": "running"}

        def on_segment(job: Any) -> None:
            if record.cancel.is_set():
                job.cancel_requested = True
            with self.lock:
                record.status = job.status()

        try:
            with self.run_lock:
                if record.cancel.is_set():
                    with self.lock:
                        record.status = {**record.status, "state": "cancelled"}
                    return
                result = design(request, self.store, on_segment=on_segment)
        except DesignError as exc:
            failed = problem(exc.code, exc.code, exc.exit_code, exc.detail)
            # The registry title is not on DesignError. Keep the code and the
            # detail; the HTTP status for a design refusal is 422 when the
            # code's own status is not known here. exit_code 2 is invalid input.
            http = 422 if exc.exit_code != 2 else 400
            failed["status"] = http
            failed["title"] = exc.code
            with self.lock:
                record.status = {
                    **_queued(request.job_id or record.status["job_id"]),
                    "state": "failed",
                    "problem": failed,
                }
            return
        with self.lock:
            record.status = result.status
            # A failed job's archive still holds the unchanged input. Publishing
            # that document would present the input as the design.
            record.front = result.front if result.job.state == "succeeded" else None
            record.lockfile = result.lockfile.to_json()

    def get(self, tenant: str, job_id: str) -> tuple[int, dict[str, Any]]:
        record = self._visible(tenant, job_id)
        if isinstance(record, dict):
            return 403, record
        return 200, record.status

    def front(
        self, tenant: str, job_id: str, method: str, max_candidates: int
    ) -> tuple[int, dict[str, Any]]:
        record = self._visible(tenant, job_id)
        if isinstance(record, dict):
            return 403, record
        if record.front is None:
            if record.status.get("state") == "failed" and record.status.get("problem"):
                return 409, record.status["problem"]
            return 409, problem(
                "CF-CLI-003",
                "A requested surface needs a capability this build does not have",
                409,
                "the job has not produced a front yet",
            )
        if method == "none":
            return 200, record.front
        resolved = "reference_direction_representatives" if method == "default" else method
        names = {
            "knee": "knee_selection",
            "clustering": "clustering",
            "epsilon_grid": "epsilon_grid",
            "reference_direction_representatives": "reference_direction_representatives",
        }
        if resolved not in names:
            return 400, problem(
                "CF-CLI-001",
                "Command-line arguments are missing or malformed",
                400,
                f"reduction {method!r} is not a presentation method",
            )
        return 200, present(record.front, names[resolved], max_candidates)

    def lockfile(self, tenant: str, job_id: str) -> tuple[int, dict[str, Any]]:
        record = self._visible(tenant, job_id)
        if isinstance(record, dict):
            return 403, record
        if record.lockfile is None:
            return 409, problem(
                "CF-CLI-003",
                "A requested surface needs a capability this build does not have",
                409,
                "the job has not sealed a lockfile yet",
            )
        return 200, record.lockfile

    def cancel(self, tenant: str, job_id: str) -> tuple[int, dict[str, Any]]:
        record = self._visible(tenant, job_id)
        if isinstance(record, dict):
            return 403, record
        if record.status["state"] in {"succeeded", "failed", "cancelled"}:
            return 409, problem(
                "CF-CLI-001",
                "Command-line arguments are missing or malformed",
                409,
                f"job is already {record.status['state']}",
            )
        record.cancel.set()
        return 202, record.status

    def list_jobs(self, tenant: str) -> dict[str, Any]:
        with self.lock:
            items = [record.status for record in self.records.values() if record.tenant == tenant]
        return {"items": items, "next_cursor": None}

    def score(self, body: dict[str, Any]) -> tuple[int, dict[str, Any] | list[Any]]:
        sequences = body.get("sequences")
        if not isinstance(sequences, list) or not sequences:
            return 400, problem(
                "CF-CLI-001",
                "Command-line arguments are missing or malformed",
                400,
                "sequences must be a non-empty list of CDS strings",
            )
        fields = {key: value for key, value in body.items() if key != "sequences"}
        try:
            request = DesignRequest.from_json(json.dumps(fields))
        except TypeError as exc:
            return 400, problem(
                "CF-CLI-001",
                "Command-line arguments are missing or malformed",
                400,
                f"score needs the library DesignRequest fields plus sequences. {exc}",
            )
        try:
            results = score_sequences(request, [str(item) for item in sequences])
        except DesignError as exc:
            return 422, problem(exc.code, exc.code, 422, exc.detail)
        return 200, {"results": results}

    def _candidate_index(self) -> dict[tuple[str, str], str]:
        found: dict[tuple[str, str], str] = {}
        with self.lock:
            records = list(self.records.values())
        for record in records:
            front = record.front
            if not front:
                continue
            for candidate in front.get("candidates") or []:
                sequences = candidate.get("sequences") or []
                if not sequences:
                    continue
                found[(record.status["job_id"], str(candidate.get("candidate_id")))] = sequences[0][
                    "digest"
                ]
        return found

    def observations(self, tenant: str) -> dict[str, Any]:
        from dna_design_suite.active.observations import ObservationLog
        from dna_design_suite.active.outcome import fit_outcome

        documents = ObservationLog(self.store.root).list(tenant)
        return {"observations": documents, "outcome_model": fit_outcome(documents)}

    def submit_observation(self, tenant: str, body: dict[str, Any]) -> tuple[int, dict[str, Any]]:
        from dna_design_suite.active.observations import (
            TITLES,
            ObservationError,
            ObservationLog,
            accept,
        )
        from dna_design_suite.active.outcome import fit_outcome

        declared = body.get("tenant_id")
        if isinstance(declared, str) and declared != tenant:
            return 403, problem(
                "CF-AUTH-003",
                "Artifact belongs to another tenant",
                403,
                "observation tenant_id does not match the caller",
            )
        try:
            stored = accept(body, candidates=self._candidate_index())
        except ObservationError as exc:
            return 422, problem(exc.code, TITLES.get(exc.code, exc.code), 422, exc.detail)
        if stored.get("tenant_id") != tenant:
            return 403, problem(
                "CF-AUTH-003",
                "Artifact belongs to another tenant",
                403,
                "observation tenant_id does not match the caller",
            )
        log = ObservationLog(self.store.root)
        try:
            stored = log.append(tenant, stored)
        except ObservationError as exc:
            return 422, problem(exc.code, TITLES.get(exc.code, exc.code), 422, exc.detail)
        documents = log.list(tenant)
        outcome = fit_outcome(documents)
        return 202, {
            "observation": stored,
            "outcome_model": outcome,
            "recalibration": _schedule_recalibration(self.store, tenant, outcome),
        }

    def _visible(self, tenant: str, job_id: str) -> Record | dict[str, Any]:
        with self.lock:
            record = self.records.get(job_id)
        if record is None or record.tenant != tenant:
            return problem(
                "CF-AUTH-003",
                "Artifact belongs to another tenant",
                403,
                "no job with this id is visible to this tenant",
            )
        return record


def _unbuilt_reason(path: str) -> str | None:
    if path in _UNBUILT:
        return f"{path}: {_UNBUILT[path]}"
    return None


def _prepare_job_body(body: dict[str, Any]) -> dict[str, Any] | tuple[int, dict[str, Any]]:
    """Turn a shipped input shape into a library DesignRequest body.

    `DEP-API-002`. Protein uses the caller's usage table. A construct uses
    the GenBank CDS features already on disk in the body. An unknown preset
    is a malformed request, not an unimplemented route.
    """
    prepared = dict(body)
    preset = prepared.get("objectives_preset")
    if preset not in (None, "none"):
        terms = _PRESETS.get(preset) if isinstance(preset, str) else None
        if terms is None:
            return 422, problem(
                "CF-CLI-001",
                "Command-line arguments are missing or malformed",
                422,
                f"objectives_preset {preset!r} is not one of {', '.join(sorted(_PRESETS))}",
            )
        prepared["objectives"] = list(terms)
    prepared.pop("objectives_preset", None)
    incoming = prepared.get("input")
    if not isinstance(incoming, dict):
        return prepared
    kind = incoming.get("kind")
    if kind == "protein":
        if not isinstance(prepared.get("host_profile"), dict):
            return 400, problem(
                "CF-CLI-001",
                "Command-line arguments are missing or malformed",
                400,
                "input.kind 'protein' needs the caller's host_profile; "
                "a missing table is not filled in",
            )
        protein = incoming.get("protein")
        if not isinstance(protein, str) or not protein:
            return 400, problem(
                "CF-CLI-001",
                "Command-line arguments are missing or malformed",
                400,
                "input.kind 'protein' needs a protein string",
            )
        table = int(prepared.get("translation_table") or 11)
        try:
            prepared["cds"] = _reverse_translate(protein, table, prepared["host_profile"])
        except ValueError as exc:
            return 422, problem(
                "CF-CLI-001",
                "Command-line arguments are missing or malformed",
                422,
                str(exc),
            )
    elif kind == "single_cds":
        sequence = incoming.get("sequence") or incoming.get("cds")
        if not isinstance(sequence, str) or not sequence:
            return 400, problem(
                "CF-CLI-001",
                "Command-line arguments are missing or malformed",
                400,
                "input.kind 'single_cds' needs the caller's sequence",
            )
        prepared["cds"] = sequence
    elif kind == "construct":
        text = incoming.get("genbank") or incoming.get("text")
        if not isinstance(text, str) or not text:
            return 400, problem(
                "CF-CLI-001",
                "Command-line arguments are missing or malformed",
                400,
                "input.kind 'construct' needs a GenBank record",
            )
        from dna_design_suite.ingest.genbank import IngestError, extract_cds_sequences

        try:
            pairs = extract_cds_sequences(text)
        except IngestError as exc:
            return 422, problem(exc.code, "Construct input has no CDS", 422, exc.detail)
        locus = incoming.get("locus")
        if len(pairs) > 1 and not locus:
            return 422, problem(
                "CF-CLI-001",
                "Command-line arguments are missing or malformed",
                422,
                "a construct with more than one CDS needs input.locus",
            )
        if locus:
            matched = [sequence for name, sequence in pairs if name == locus]
            if len(matched) != 1:
                return 422, problem(
                    "CF-CLI-001",
                    "Command-line arguments are missing or malformed",
                    422,
                    f"input.locus {locus!r} does not name exactly one CDS",
                )
            prepared["cds"] = matched[0]
        else:
            prepared["cds"] = pairs[0][1]
    elif kind is not None:
        return 422, problem(
            "CF-CLI-001",
            "Command-line arguments are missing or malformed",
            422,
            f"input.kind {kind!r} is not protein, construct, or single_cds",
        )
    prepared.pop("input", None)
    return prepared


def _usage_weight(profile: dict[str, Any], codon: str) -> float:
    table = (profile.get("codon_usage") or {}).get("table") or {}
    entry = table.get(codon) or {}
    weight = entry.get("relative_adaptiveness", 0.0)
    try:
        return float(weight)
    except (TypeError, ValueError):
        return 0.0


def _reverse_translate(protein: str, translation_table: int, profile: dict[str, Any]) -> str:
    """One codon per residue, from the caller's table.

    The same rule as `cai_greedy`: highest relative adaptiveness, and the
    lexicographically smaller codon on a tie. A codon absent from the table
    has weight 0. Stops are whatever the translation table calls a stop.
    """
    from dna_design_suite import _core

    bases = "ACGT"
    by_amino: dict[str, list[str]] = {}
    for first in bases:
        for second in bases:
            for third in bases:
                codon = first + second + third
                amino = _core.translate(codon, translation_table)
                by_amino.setdefault(amino, []).append(codon)
    sequence = protein if protein.endswith("*") else protein + "*"
    chosen: list[str] = []
    for amino in sequence:
        options = by_amino.get(amino)
        if not options:
            raise ValueError(f"{amino!r} is not a residue in translation table {translation_table}")
        chosen.append(min(options, key=lambda codon: (-_usage_weight(profile, codon), codon)))
    return "".join(chosen)


def _schedule_recalibration(
    store: ArtifactStore, tenant: str, outcome: dict[str, Any]
) -> dict[str, Any]:
    """Store the Phase 5 OLS fit when it is fitted. Otherwise do not schedule.

    The stored document is the coefficients `fit_outcome` already computed.
    It does not move a model channel and it does not invent a second model.
    """
    if outcome.get("state") != "fitted":
        return {
            "state": "not_scheduled",
            "reason": "fit_outcome is not fitted, so no recalibration is stored",
        }
    record = {
        "tenant_id": tenant,
        "model": "phase5_ols",
        "feature": outcome.get("feature"),
        "measures": outcome.get("measures"),
        "n": outcome.get("n"),
        "coefficients": outcome.get("coefficients"),
        "residual_sum_of_squares": outcome.get("residual_sum_of_squares"),
    }
    canonical = json.dumps(record, sort_keys=True, separators=(",", ":"))
    digest = "sha256:" + hashlib.sha256(canonical.encode()).hexdigest()
    record["digest"] = digest
    directory = Path(store.root) / "recalibrations" / tenant
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{digest.removeprefix('sha256:')}.json"
    if not path.exists():
        path.write_text(json.dumps(record, indent=2, sort_keys=True), encoding="utf-8")
    return {"state": "fitted", "digest": digest}


def list_recalibrations(store: ArtifactStore, tenant: str) -> dict[str, Any]:
    directory = Path(store.root) / "recalibrations" / tenant
    records: list[dict[str, Any]] = []
    if directory.is_dir():
        for path in sorted(directory.glob("*.json")):
            document = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(document, dict):
                records.append(document)
    return {"recalibrations": records}


def list_hosts() -> dict[str, Any]:
    """Example profiles that are actually on disk. No catalogue is invented."""
    directory = _examples_dir()
    hosts = []
    for name in _HOST_FILES:
        path = directory / name
        if not path.is_file():
            continue
        hosts.append(
            {
                "id": path.stem,
                "filename": name,
                "installed": True,
                "profile": json.loads(path.read_text(encoding="utf-8")),
            }
        )
    return {"hosts": hosts}


def list_models() -> dict[str, Any]:
    """The verified Potts file, or an empty list when it is not installed."""
    from dna_design_suite.codon.epistasis.artifact import (
        ArtifactError,
        artifact_path,
        verify_artifact,
    )

    path = artifact_path()
    if not path.is_file():
        return {"models": [], "installed": False}
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
        verify_artifact(document)
    except (ArtifactError, json.JSONDecodeError, OSError, KeyError, TypeError):
        return {"models": [], "installed": False}
    return {
        "models": [
            {
                "id": "epistasis.context_naturalness",
                "filename": path.name,
                "installed": True,
                "digest": document.get("digest"),
            }
        ],
        "installed": True,
    }


def analyse_feasibility(body: dict[str, Any]) -> tuple[int, dict[str, Any]]:
    """The five checks in `cf-constraints`, on the construct the caller sent."""
    from dna_design_suite import _core

    cds = body.get("cds")
    if not isinstance(cds, str) or not cds:
        return 400, problem(
            "CF-CLI-001",
            "Command-line arguments are missing or malformed",
            400,
            "POST /v1/feasibility needs cds",
        )
    constraints = body.get("constraints")
    if constraints is None:
        constraints = body.get("constraint_set")
    if not isinstance(constraints, dict):
        constraints = {
            "schema_version": "1.0.0",
            "constraint_set_id": "cst_none",
            "constraints": [],
        }
    flank5 = body.get("five_prime_flank") or ""
    flank3 = body.get("three_prime_flank") or ""
    if not isinstance(flank5, str) or not isinstance(flank3, str):
        return 400, problem(
            "CF-CLI-001",
            "Command-line arguments are missing or malformed",
            400,
            "flanks must be strings",
        )
    try:
        table = int(body.get("translation_table") or 11)
    except (TypeError, ValueError):
        return 400, problem(
            "CF-CLI-001",
            "Command-line arguments are missing or malformed",
            400,
            "translation_table must be an integer",
        )
    construct = flank5 + cds + flank3
    try:
        report = json.loads(
            _core.feasibility_analyse(
                json.dumps(constraints),
                construct,
                cds,
                table,
                len(flank5),
            )
        )
    except ValueError as exc:
        return 422, problem(
            "CF-FEAS-005",
            "The constraint set is not feasible to analyse",
            422,
            str(exc),
        )
    return 200, report


def _queued(job_id: str) -> dict[str, Any]:
    from datetime import UTC, datetime

    return {
        "schema_version": _SCHEMA,
        "job_id": job_id,
        "state": "queued",
        "created_at": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "progress": {"current_tier": "none", "tiers": []},
        "links": {"self": f"/v1/jobs/{job_id}"},
        "warnings": [],
    }


def _tenant(header: str | None) -> str | None:
    if not header or not header.startswith("Bearer "):
        return None
    token = header.removeprefix("Bearer ").strip()
    if len(token) < 10 or not token.startswith("t_"):
        return None
    body = token[2:]
    if not 8 <= len(body) <= 32 or not body.isalnum():
        return None
    return token


class Handler(BaseHTTPRequestHandler):
    server: Server

    def log_message(self, fmt: str, *args: Any) -> None:
        return

    def do_GET(self) -> None:  # noqa: N802
        self._dispatch("GET")

    def do_POST(self) -> None:  # noqa: N802
        self._dispatch("POST")

    def _dispatch(self, method: str) -> None:
        parsed = urlparse(self.path)
        path = parsed.path.rstrip("/") or "/"
        if method == "GET" and path == "/favicon.ico":
            self.send_response(204)
            self.end_headers()
            return
        if method == "GET" and path in _STATIC_FILES:
            name, content_type = _STATIC_FILES[path]
            payload = (_STATIC / name).read_bytes()
            self.send_response(200)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
            return
        if method == "GET" and path in _EXAMPLE_FILES:
            name, content_type = _EXAMPLE_FILES[path]
            file = _examples_dir() / name
            if not file.is_file():
                message = (
                    f"{name} is not in this installation. Upload a codon-usage "
                    "table from the page instead.\n"
                ).encode()
                self.send_response(404)
                self.send_header("Content-Type", "text/plain; charset=utf-8")
                self.send_header("Content-Length", str(len(message)))
                self.end_headers()
                self.wfile.write(message)
                return
            payload = file.read_bytes()
            self.send_response(200)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
            return
        reason = _unbuilt_reason(path)
        if reason:
            self._send(
                501,
                problem(
                    "CF-CLI-003",
                    "A requested surface needs a capability this build does not have",
                    501,
                    reason,
                ),
                problem_body=True,
            )
            return
        if path == "/v1/hosts" and method == "GET":
            tenant = _tenant(self.headers.get("Authorization"))
            if tenant is None:
                self._unauthorized()
                return
            self._send(200, list_hosts())
            return
        if path == "/v1/models" and method == "GET":
            tenant = _tenant(self.headers.get("Authorization"))
            if tenant is None:
                self._unauthorized()
                return
            self._send(200, list_models())
            return
        if path == "/v1/recalibrations" and method == "GET":
            tenant = _tenant(self.headers.get("Authorization"))
            if tenant is None:
                self._unauthorized()
                return
            self._send(200, list_recalibrations(self.server.api.store, tenant))
            return
        if path == "/v1/feasibility" and method == "POST":
            tenant = _tenant(self.headers.get("Authorization"))
            if tenant is None:
                self._unauthorized()
                return
            status, payload = analyse_feasibility(self._json())
            self._send(status, payload, problem_body=status >= 400)
            return
        if path == "/v1/version" and method == "GET":
            self._send(200, self.server.api.version())
            return
        if path == "/v1/healthz" and method == "GET":
            self._send(200, self.server.api.health())
            return
        if path == "/v1/calibration" and method == "GET":
            from dna_design_suite.codon.epistasis.artifact import ArtifactError, load_calibration

            try:
                from dna_design_suite.codon.epistasis.dwell import dwell_status

                payload = load_calibration()
                payload["dwell"] = dwell_status()
            except (ArtifactError, json.JSONDecodeError, OSError) as exc:
                self._send(
                    500,
                    problem(
                        "CF-RESOLVE-006",
                        "Artifact content does not match its digest",
                        500,
                        str(exc),
                    ),
                    problem_body=True,
                )
                return
            self._send(200, payload)
            return
        tenant = _tenant(self.headers.get("Authorization"))
        if tenant is None:
            self._send(
                401,
                problem(
                    "CF-AUTH-001",
                    "Missing or malformed credential",
                    401,
                    "Authorization must be Bearer t_<8-32 alphanumerics>",
                ),
                problem_body=True,
            )
            return
        if path == "/v1/jobs" and method == "GET":
            self._send(200, self.server.api.list_jobs(tenant))
            return
        if path == "/v1/jobs" and method == "POST":
            status, payload = self.server.api.submit(
                tenant, self._json(), self.headers.get("Idempotency-Key")
            )
            self._send(status, payload, problem_body=status >= 400)
            return
        if path == "/v1/primers/candidates" and method == "POST":
            from dna_design_suite.primer.service import submit as submit_primers

            status, payload = submit_primers(self.server.api, tenant, self._json())
            self._send(status, payload, problem_body=status >= 400)
            return
        if path == "/v1/primers/analyze" and method == "POST":
            from dna_design_suite.primer.service import analyze as analyze_primers

            body = self._json()
            forward = body.get("forward")
            reverse = body.get("reverse")
            if not isinstance(forward, str) or not isinstance(reverse, str):
                self._send(
                    422,
                    problem(
                        "CF-PRIM-006",
                        "Primer search configuration is malformed",
                        422,
                        "POST /v1/primers/analyze requires forward and reverse strings",
                    ),
                    problem_body=True,
                )
                return
            try:
                self._send(200, analyze_primers(forward, reverse))
            except ValueError as exc:
                message = str(exc)
                code = message.split()[0] if message.startswith("CF-") else "CF-PRIM-006"
                self._send(
                    422,
                    problem(code, "Primer analysis failed", 422, message),
                    problem_body=True,
                )
            return
        if path == "/v1/score" and method == "POST":
            status, payload = self.server.api.score(self._json())
            self._send(status, payload, problem_body=status >= 400)
            return
        if path == "/v1/observations" and method == "GET":
            self._send(200, self.server.api.observations(tenant))
            return
        if path == "/v1/observations" and method == "POST":
            status, payload = self.server.api.submit_observation(tenant, self._json())
            self._send(status, payload, problem_body=status >= 400)
            return
        parts = [part for part in path.split("/") if part]
        if len(parts) >= 3 and parts[0] == "v1" and parts[1] == "primers":
            job_id = parts[2]
            record = self.server.api.primer_records.get(job_id)
            if record is None or record.tenant != tenant:
                self._send(
                    404,
                    problem(
                        "CF-AUTH-003",
                        "The job is not visible to this caller",
                        404,
                        "No primer job with that id is visible.",
                    ),
                    problem_body=True,
                )
                return
            if method == "GET" and len(parts) == 3:
                self._send(200, _primer_payload(record.status, parsed.query))
                return
            if method == "GET" and len(parts) == 5 and parts[3] == "pair":
                pair = record.pairs.get(parts[4])
                if pair is None:
                    self._send(
                        404,
                        problem(
                            "CF-PRIM-003",
                            "The amplification region is outside the sequence",
                            404,
                            f"No pair {parts[4]} on this job.",
                        ),
                        problem_body=True,
                    )
                    return
                self._send(200, _primer_payload(pair, parsed.query))
                return
        if len(parts) >= 3 and parts[0] == "v1" and parts[1] == "jobs":
            job_id = parts[2]
            tail = parts[3] if len(parts) > 3 else ""
            if method == "GET" and tail == "":
                status, payload = self.server.api.get(tenant, job_id)
                self._send(status, payload, problem_body=status >= 400)
                return
            if method == "GET" and tail == "front":
                query = parse_qs(parsed.query)
                reduction = query.get("reduction", ["default"])[0]
                try:
                    limit = int(query.get("max_candidates", ["20"])[0])
                except ValueError:
                    limit = 0
                status, payload = self.server.api.front(tenant, job_id, reduction, limit)
                self._send(status, payload, problem_body=status >= 400)
                return
            if method == "GET" and tail == "lockfile":
                status, payload = self.server.api.lockfile(tenant, job_id)
                self._send(status, payload, problem_body=status >= 400)
                return
            if method == "GET" and tail == "events":
                status, payload = self.server.api.get(tenant, job_id)
                if status != 200:
                    self._send(status, payload, problem_body=True)
                    return
                data = json.dumps(payload).encode()
                body = b"event: status\ndata: " + data + b"\n\n"
                self.send_response(200)
                self.send_header("Content-Type", "text/event-stream")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return
            if method == "POST" and tail == "cancel":
                status, payload = self.server.api.cancel(tenant, job_id)
                self._send(status, payload, problem_body=status >= 400)
                return
        self._send(
            501,
            problem(
                "CF-CLI-003",
                "A requested surface needs a capability this build does not have",
                501,
                f"{method} {path} is not implemented in this build",
            ),
            problem_body=True,
        )

    def _unauthorized(self) -> None:
        self._send(
            401,
            problem(
                "CF-AUTH-001",
                "Missing or malformed credential",
                401,
                "Authorization must be Bearer t_<8-32 alphanumerics>",
            ),
            problem_body=True,
        )

    def _json(self) -> dict[str, Any]:
        length = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(length) if length else b"{}"
        try:
            value = json.loads(raw.decode())
        except json.JSONDecodeError:
            return {}
        return value if isinstance(value, dict) else {}

    def _send(self, status: int, payload: Any, problem_body: bool = False) -> None:
        body = json.dumps(payload).encode()
        self.send_response(status)
        kind = "application/problem+json" if problem_body else "application/json"
        self.send_header("Content-Type", kind)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


class Server(ThreadingHTTPServer):
    def __init__(self, address: tuple[str, int], api: Api) -> None:
        super().__init__(address, Handler)
        self.api = api


def _primer_payload(payload: dict[str, Any], query: str) -> dict[str, Any]:
    """Primer JSON for one response.

    `product_sequence` is omitted unless `?product_sequence=1`. The default
    is documented on `dna_design_suite.primer.service.wants_product_sequence`.
    """
    from dna_design_suite.primer.service import project_primer_payload, wants_product_sequence

    return project_primer_payload(payload, wants_product_sequence(query))


def _examples_dir() -> Path:
    """Demo files the page fetches.

    A source checkout uses ``examples/`` at the repo root, so the page and
    the README read one copy. An installed wheel has no checkout; it serves
    the copy staged under ``dna_design_suite/api/examples``, which
    ``tools/stage_wheel_examples.py`` keeps byte-identical.
    """
    here = Path(__file__).resolve()
    packaged = here.parent / "examples"
    candidates = [here.parents[3] / "examples", packaged, Path.cwd() / "examples"]
    for path in candidates:
        if (path / "mg1655_ribosomal_usage.json").is_file():
            return path
    return candidates[0]


def serve(store: ArtifactStore, host: str = "127.0.0.1", port: int = 0) -> Server:
    server = Server((host, port), Api(store))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server


def run_page(
    store: ArtifactStore | None = None,
    host: str = "127.0.0.1",
    port: int = 8000,
    *,
    open_browser: bool = True,
    opener: Any = None,
    block: bool = True,
) -> Server:
    """Serve the page for a person at the keyboard.

    ``serve`` starts a daemon thread and returns, which is what tests use.
    This listens on port 8000, opens a browser, and blocks until interrupted.
    """
    import webbrowser

    if store is None:
        store = ArtifactStore(Path(".dna-design-suite/artifacts"))
    server = Server((host, port), Api(store))
    bound_host, bound_port = server.server_address[:2]
    url = f"http://{bound_host}:{bound_port}/"
    if open_browser:
        (opener or webbrowser.open)(url)
    if not block:
        threading.Thread(target=server.serve_forever, daemon=True).start()
        return server
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        server.shutdown()
    return server
