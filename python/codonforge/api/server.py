"""HTTP surface for the library. `P2-18`, `DEP-API-001`.

The routes are the OpenAPI paths a caller can exercise today. The body of
`POST /v1/jobs` is the library `DesignRequest`, the same object the CLI runs.
Routes whose backing module is not in this build answer `CF-CLI-003` rather
than an empty document that looks like success.

Jobs run on a thread via `design(background=False)`. `background=True` stays
the `CF-CLI-003` refusal it already is (`DEP-CLI-004`); the thread is a real
call, not a handle that pretends the work finished.
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

from codonforge import __version__
from codonforge.decision.reduce import present
from codonforge.modelreg import ArtifactStore
from codonforge.orchestrator import (
    DesignError,
    DesignRequest,
    assigned_job_id,
    design,
    score_sequences,
)

_SCHEMA = "1.0.0"
# Reasons are the product. An empty 200 would look implemented.
_UNBUILT = {
    "/v1/hosts": (
        "No host catalog is installed. A job carries the caller's host_profile; "
        "an empty list would look like a catalog."
    ),
    "/v1/models": (
        "Public model channels are not this build (P3-15). This route does not "
        "list the local store as if it were a catalog."
    ),
    "/v1/observations": "Wet-lab observation ingestion is not implemented.",
    "/v1/recalibrations": "Recalibration jobs are not implemented.",
    "/v1/seeds:import": "Importing an external tool's output is not wired on this route.",
    "/v1/feasibility": (
        "The five-check analyser is not this build. Feasibility of a finished "
        "candidate is constraint_report.feasible on the front, a boolean, with "
        "unknown_terms beside it."
    ),
}
_STATIC = Path(__file__).resolve().parent / "static"
_STATIC_FILES = {
    "/": ("index.html", "text/html; charset=utf-8"),
    "/index.html": ("index.html", "text/html; charset=utf-8"),
    "/app.js": ("app.js", "text/javascript; charset=utf-8"),
    "/codonforge-ui.css": ("codonforge-ui.css", "text/css; charset=utf-8"),
}


def problem(code: str, title: str, status: int, detail: str) -> dict[str, Any]:
    return {
        "type": f"https://codonforge.org/problems/{code.lower()}",
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
        return {
            "status": "ok",
            "workers_ready": 0,
            "gpu_available": False,
            "queue_depth": depth,
        }

    def submit(
        self, tenant: str, body: dict[str, Any], idempotency_key: str | None
    ) -> tuple[int, dict[str, Any]]:
        refused = _unbuilt_job_body(body)
        if refused:
            return 501, problem(
                "CF-CLI-003",
                "A requested surface needs a capability this build does not have",
                501,
                refused,
            )
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
            record.front = result.front
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
    for prefix, key in (
        ("/v1/hosts/", "/v1/hosts"),
        ("/v1/models/", "/v1/models"),
        ("/v1/observations/", "/v1/observations"),
    ):
        if path.startswith(prefix):
            return f"{path}: {_UNBUILT[key]}"
    return None


def _unbuilt_job_body(body: dict[str, Any]) -> str | None:
    """OpenAPI shapes this build does not translate into a DesignRequest.

    Returning 400 would call a documented shape a typo. Returning a stub
    document would look like the shape had been designed.
    """
    incoming = body.get("input")
    if isinstance(incoming, dict):
        kind = incoming.get("kind")
        if kind in {"protein", "construct", "single_cds"}:
            return (
                f"input.kind {kind!r} is not accepted. This build designs a "
                "library DesignRequest (host, cds, host_profile). A protein is "
                "not reverse-translated, and a construct is not expanded."
            )
    preset = body.get("objectives_preset")
    if preset not in (None, "none"):
        return (
            f"objectives_preset {preset!r} is not expanded in this build. "
            "Send objectives as a list of term ids."
        )
    return None


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
        if path == "/v1/version" and method == "GET":
            self._send(200, self.server.api.version())
            return
        if path == "/v1/healthz" and method == "GET":
            self._send(200, self.server.api.health())
            return
        if path == "/v1/calibration" and method == "GET":
            from codonforge.epistasis.artifact import ArtifactError, load_calibration

            try:
                payload = load_calibration()
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
            from codonforge.primers.service import submit as submit_primers

            status, payload = submit_primers(self.server.api, tenant, self._json())
            self._send(status, payload, problem_body=status >= 400)
            return
        if path == "/v1/primers/analyze" and method == "POST":
            from codonforge.primers.service import analyze as analyze_primers

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
                self._send(200, record.status)
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
                self._send(200, pair)
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


def serve(store: ArtifactStore, host: str = "127.0.0.1", port: int = 0) -> Server:
    server = Server((host, port), Api(store))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server
