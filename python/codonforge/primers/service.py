"""HTTP-facing primer jobs.

The search itself is `cf-primer-search`. This module only keeps the job
document in the shape `job_status.schema.json` already uses, so the page can
poll it the same way it polls a codon job.
"""

from __future__ import annotations

import json
import threading
from collections.abc import Callable
from datetime import UTC, datetime
from typing import Any

from codonforge import _core
from codonforge.orchestrator.design import _derived_id
from codonforge.primers.sequence import clean_template

SCHEMA = "1.0.0"
# Tests set this to a callable. Production leaves it as a no-op so a search
# is not delayed for the sake of a poll.
PROGRESS_GATE: Callable[[], None] | None = None


def _now() -> str:
    return datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _tier(name: str, state: str, fraction: float) -> dict[str, Any]:
    return {"tier": name, "state": state, "fraction": fraction}


def _status(job_id: str, state: str, fraction: float, tier: str) -> dict[str, Any]:
    generate = "done" if fraction > 0.02 else "running" if state == "running" else "pending"
    prune = "done" if fraction >= 0.25 else "running" if fraction > 0.02 else "pending"
    evaluate = "done" if fraction >= 1.0 else "running" if fraction >= 0.25 else "pending"
    evaluated = 0.0
    if fraction >= 0.25:
        evaluated = min(1.0, max(0.0, (fraction - 0.25) / 0.75))
    return {
        "schema_version": SCHEMA,
        "job_id": job_id,
        "state": state,
        "created_at": _now(),
        "progress": {
            "current_tier": tier if state == "running" else "none",
            "tiers": [
                _tier("primer_generate", generate, 1.0 if fraction > 0.02 else fraction),
                _tier("primer_prune", prune, 1.0 if fraction >= 0.25 else 0.0),
                _tier("primer_evaluate", evaluate, evaluated),
                _tier(
                    "primer_rank",
                    "done" if state == "succeeded" else "pending",
                    1.0 if state == "succeeded" else 0.0,
                ),
            ],
        },
        "links": {"self": f"/v1/primers/{job_id}"},
        "warnings": [],
    }


def _problem(message: str) -> dict[str, Any]:
    code = message.split()[0] if message.startswith("CF-") else "CF-PRIM-006"
    return {
        "type": f"https://codonforge.org/problems/{code.lower()}",
        "title": message,
        "status": 422,
        "code": code,
        "detail": message,
        "stage": "primer_evaluate",
        "subject": "primers",
    }


def analyze(forward: str, reverse: str) -> dict[str, Any]:
    forward_clean, forward_problem = clean_template(forward)
    reverse_clean, reverse_problem = clean_template(reverse)
    problem_text = forward_problem or reverse_problem
    if problem_text:
        raise ValueError(problem_text)
    return json.loads(_core.primer_analyze(forward_clean, reverse_clean))


def submit(api: Any, tenant: str, body: dict[str, Any]) -> tuple[int, dict[str, Any]]:
    try:
        job_id = _derived_id("job", json.dumps(body, sort_keys=True))
    except (TypeError, ValueError) as exc:
        return 422, _problem(f"CF-PRIM-006 {exc}")
    document = _status(job_id, "queued", 0.0, "primer_generate")
    with api.lock:
        api.primer_records[job_id] = _Record(tenant, document, {})
    thread = threading.Thread(target=_run, args=(api, job_id, body), daemon=True)
    thread.start()
    return 202, document


def _run(api: Any, job_id: str, body: dict[str, Any]) -> None:
    def publish(
        state: str,
        fraction: float,
        tier: str,
        result: dict[str, Any] | None = None,
        problem: dict[str, Any] | None = None,
    ) -> None:
        with api.lock:
            record = api.primer_records[job_id]
            created = record.status["created_at"]
            record.status = _status(job_id, state, fraction, tier)
            record.status["created_at"] = created
            if result is not None:
                record.status["primer_result"] = result
                record.pairs = {pair["pair_id"]: pair for pair in result.get("pairs", [])}
            if problem is not None:
                record.status["problem"] = problem
                record.status["state"] = "failed"

    publish("running", 0.02, "primer_generate")
    if PROGRESS_GATE is not None:
        PROGRESS_GATE()

    def on_progress(fraction: float) -> None:
        publish("running", float(fraction), "primer_evaluate")

    sequence = body.get("sequence")
    if not isinstance(sequence, str):
        publish("failed", 1.0, "primer_rank", problem=_problem("CF-PRIM-006 sequence is required"))
        return
    cleaned, problem_text = clean_template(sequence)
    if problem_text:
        publish("failed", 1.0, "primer_rank", problem=_problem(problem_text))
        return
    body = {**body, "sequence": cleaned}
    try:
        raw = _core.primer_search(json.dumps(body), on_progress)
    except ValueError as exc:
        publish("failed", 1.0, "primer_rank", problem=_problem(str(exc)))
        return
    result = json.loads(raw)
    publish("succeeded", 1.0, "primer_rank", result=result)
    if result.get("truncated"):
        with api.lock:
            api.primer_records[job_id].status["warnings"].append(
                {
                    "type": "https://codonforge.org/problems/cf-prim-007",
                    "title": "Primer candidates were limited to those closest to the Tm target",
                    "status": 200,
                    "code": "CF-PRIM-007",
                    "detail": "More than 40 primers on a flank survived the Tm pre-filter.",
                    "stage": "primer_prune",
                    "subject": "flank",
                }
            )


class _Record:
    def __init__(self, tenant: str, status: dict[str, Any], pairs: dict[str, Any]) -> None:
        self.tenant = tenant
        self.status = status
        self.pairs = pairs
