"""TST-API-001, TST-API-002, TST-API-004 for the routes this build serves.

`DEP-API-001`. A job is the library `DesignRequest`. Idempotency replays the
same body and refuses a different one. A second tenant does not see the job.
Unimplemented routes are `CF-CLI-003` with `application/problem+json`.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from pathlib import Path

import jsonschema
import pytest
from referencing import Registry, Resource
from test_orchestrator import CDS, HOST, constraint_set, host_profile

from codonforge.api import serve
from codonforge.modelreg import ArtifactStore
from codonforge.orchestrator import score_sequences
from codonforge.orchestrator.design import DesignRequest

TENANT = "t_localdev1"
OTHER = "t_otherdev1"
SCHEMA = Path(__file__).resolve().parents[2] / "schemas" / "v1" / "job_status.schema.json"


def _call(method: str, url: str, headers: dict[str, str], body: dict | None = None):
    data = None if body is None else json.dumps(body).encode()
    request = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(request) as response:
            raw = response.read()
            payload = json.loads(raw) if raw else None
            return response.status, response.headers.get("Content-Type"), payload
    except urllib.error.HTTPError as exc:
        return exc.code, exc.headers.get("Content-Type"), json.loads(exc.read())


def _headers(tenant: str = TENANT, key: str | None = None) -> dict[str, str]:
    headers = {"Authorization": f"Bearer {tenant}", "Content-Type": "application/json"}
    if key:
        headers["Idempotency-Key"] = key
    return headers


def _body(**overrides: object) -> dict:
    document = {
        "host": HOST,
        "cds": CDS,
        "host_profile": host_profile(),
        "objectives": ["usage.cai"],
        "constraints": constraint_set(),
        "translation_table": 1,
        "budget": "fast",
        "seed": 7,
        "generations": 2,
        "segment": 2,
        "host_digest": "sha256:" + "1" * 64,
    }
    document.update(overrides)
    return document


@pytest.fixture()
def base(tmp_path: Path):
    server = serve(ArtifactStore(tmp_path / "artifacts"))
    url = f"http://127.0.0.1:{server.server_address[1]}"
    yield url
    server.shutdown()


def test_version_and_health_need_no_credential(base: str) -> None:
    status, _kind, payload = _call("GET", f"{base}/v1/version", {})
    assert status == 200
    assert payload["contract_majors"] == [1]
    assert payload["git_sha"]
    status, _kind, health = _call("GET", f"{base}/v1/healthz", {})
    assert status == 200
    assert health["status"] == "ok"
    assert health["workers_ready"] == 0


def test_a_missing_credential_is_a_problem(base: str) -> None:
    status, kind, payload = _call("GET", f"{base}/v1/jobs", {})
    assert status == 401
    assert kind.startswith("application/problem+json")
    assert payload["code"] == "CF-AUTH-001"


def test_an_unimplemented_route_is_cf_cli_003(base: str) -> None:
    status, kind, payload = _call("GET", f"{base}/v1/hosts", _headers())
    assert status == 501
    assert kind.startswith("application/problem+json")
    assert payload["code"] == "CF-CLI-003"
    assert payload["status"] == 501


def test_idempotency_replays_and_a_different_body_conflicts(base: str) -> None:
    schema_dir = SCHEMA.parent
    registry = Registry().with_resources(
        [
            (path.name, Resource.from_contents(json.loads(path.read_text())))
            for path in schema_dir.glob("*.schema.json")
        ]
    )
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    validator = jsonschema.Draft202012Validator(schema, registry=registry)
    # The schema references siblings. Validate the structural required keys
    # the handler itself emits; the finished job document is the library's,
    # which the orchestrator tests already validate end to end.
    first, kind, status_doc = _call("POST", f"{base}/v1/jobs", _headers(key="same"), _body())
    assert first == 202
    assert kind.startswith("application/json")
    assert status_doc["state"] == "queued"
    for key in ("schema_version", "job_id", "state", "created_at", "progress", "links"):
        assert key in status_doc
    validator.validate(status_doc)

    again, _kind, replayed = _call("POST", f"{base}/v1/jobs", _headers(key="same"), _body())
    assert again == 200
    assert replayed["job_id"] == status_doc["job_id"]

    conflict, kind, problem = _call("POST", f"{base}/v1/jobs", _headers(key="same"), _body(seed=8))
    assert conflict == 409
    assert kind.startswith("application/problem+json")
    assert problem["code"] == "CF-IDEM-001"

    other, _kind, hidden = _call("GET", f"{base}/v1/jobs/{status_doc['job_id']}", _headers(OTHER))
    assert other == 403
    assert hidden["code"] == "CF-AUTH-003"

    deadline = time.time() + 30
    finished = status_doc
    while finished["state"] in {"queued", "running"} and time.time() < deadline:
        time.sleep(0.05)
        _status, _kind, finished = _call(
            "GET", f"{base}/v1/jobs/{status_doc['job_id']}", _headers()
        )
    assert finished["state"] == "succeeded"
    validator.validate(finished)
    front_status, _kind, front = _call(
        "GET",
        f"{base}/v1/jobs/{status_doc['job_id']}/front?reduction=default&max_candidates=5",
        _headers(),
    )
    assert front_status == 200
    assert front["reduction"]["method"] == "reference_direction_representatives"
    assert front["candidates"]


def test_score_matches_the_library_helper(base: str) -> None:
    body = _body()
    body["sequences"] = [CDS, CDS]
    request = DesignRequest.from_json(
        json.dumps({k: v for k, v in body.items() if k != "sequences"})
    )
    expected = score_sequences(request, [CDS, CDS])
    status, _kind, payload = _call("POST", f"{base}/v1/score", _headers(), body)
    assert status == 200
    assert payload["results"] == expected
