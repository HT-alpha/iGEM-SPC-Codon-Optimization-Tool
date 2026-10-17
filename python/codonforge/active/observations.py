"""Validate and append a wet-lab observation.

The eight checks are the ones in
`planning/04-api/03-active-learning-ingestion.md`. A failure appends nothing.
Recalibration is a different route and is not scheduled here.
"""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any

from jsonschema import Draft202012Validator, FormatChecker
from referencing import Registry, Resource

from codonforge import _core
from codonforge.modelreg.manifest import digest_of

_HOST = re.compile(r"^[a-z0-9]+(?:_[a-z0-9]+)*:[0-9]+(?:/[A-Za-z0-9._-]+)?$")
_CONDITION = re.compile(r"^cond_[A-Za-z0-9_.-]{1,48}$")
_CENSORING = {"none", "left", "right", "interval"}
_MEASURES = {
    "protein_abundance",
    "mrna_abundance",
    "translation_efficiency",
    "host_fitness",
    "solubility",
    "activity",
}
_SELECTION = {
    "top_ranked",
    "uniform_random_from_front",
    "max_uncertainty",
    "random_feasible",
    "user_chosen",
    "unknown",
}

TITLES = {
    "CF-OBS-001": "Observation failed schema validation",
    "CF-OBS-002": "As-built digest does not match the referenced candidate",
    "CF-OBS-003": "Submission contains no control",
    "CF-OBS-004": "Assay dynamic range and censoring must be declared",
    "CF-OBS-005": "Latent measured quantity must be declared",
    "CF-OBS-006": "Exploration-arm assignment must be declared",
    "CF-OBS-007": "Host or growth condition could not be resolved",
    "CF-OBS-008": "Replicate count is inconsistent with the reported values",
}


class ObservationError(Exception):
    """One observation check failed. Nothing was stored."""

    def __init__(self, code: str, detail: str) -> None:
        self.code = code
        self.detail = detail
        super().__init__(f"{code}: {detail}")


def _schema_dir() -> Path:
    return Path(__file__).resolve().parents[3] / "schemas" / "v1"


@lru_cache(maxsize=1)
def _validator() -> Draft202012Validator:
    schema_dir = _schema_dir()
    resources = [
        (path.name, Resource.from_contents(json.loads(path.read_text(encoding="utf-8"))))
        for path in sorted(schema_dir.glob("*.schema.json"))
    ]
    registry = Registry().with_resources(resources)
    schema = json.loads((schema_dir / "observation.schema.json").read_text(encoding="utf-8"))
    return Draft202012Validator(schema, registry=registry, format_checker=FormatChecker())


def _fail(code: str, detail: str) -> None:
    raise ObservationError(code, detail)


def _semantic(
    document: dict[str, Any],
    *,
    candidates: dict[tuple[str, str], str] | None,
    translation_table: int,
) -> None:
    controls = document.get("controls")
    if not isinstance(controls, list) or len(controls) < 1:
        _fail("CF-OBS-003", "at least one control is required")
    exploration = document.get("exploration")
    if (
        not isinstance(exploration, dict)
        or "is_exploration_arm" not in exploration
        or exploration.get("selection_rule") not in _SELECTION
        or not isinstance(exploration.get("is_exploration_arm"), bool)
    ):
        _fail("CF-OBS-006", "exploration.is_exploration_arm and selection_rule are required")
    readouts = document.get("readouts")
    if not isinstance(readouts, list) or len(readouts) < 1:
        _fail("CF-OBS-005", "readouts must declare the latent quantity they measure")
    for readout in readouts:
        if not isinstance(readout, dict):
            _fail("CF-OBS-004", "a readout must declare dynamic_range.censoring")
        dynamic = readout.get("dynamic_range")
        if not isinstance(dynamic, dict) or dynamic.get("censoring") not in _CENSORING:
            _fail("CF-OBS-004", "dynamic_range.censoring must be declared")
        if readout.get("measures") not in _MEASURES:
            _fail("CF-OBS-005", "readout.measures must name the latent quantity")
        values = readout.get("replicate_values", [])
        count = readout.get("n_replicates")
        if isinstance(values, list) and len(values) > 0 and count != len(values):
            _fail(
                "CF-OBS-008",
                f"n_replicates is {count} and replicate_values has {len(values)} entries",
            )
    host = document.get("host_id")
    condition = document.get("condition")
    condition_id = condition.get("condition_id") if isinstance(condition, dict) else None
    if not isinstance(host, str) or _HOST.match(host) is None:
        _fail("CF-OBS-007", f"host_id {host!r} is not a resolvable host identifier")
    if not isinstance(condition_id, str) or _CONDITION.match(condition_id) is None:
        _fail("CF-OBS-007", f"condition_id {condition_id!r} is not resolvable")
    subject = document.get("subject")
    if not isinstance(subject, dict):
        _fail("CF-OBS-002", "subject is missing, so the as-built sequence cannot be checked")
    kind = subject.get("kind")
    if kind == "explicit":
        cds = subject.get("cds")
        if not isinstance(cds, str):
            _fail("CF-OBS-002", "explicit subject has no cds to hash")
        try:
            sequence = json.loads(_core.sequence_document(cds, translation_table, None, None, None))
        except Exception as exc:
            _fail(
                "CF-OBS-002",
                "explicit cds is not a coding sequence under translation table "
                f"{translation_table}, which is the table this digest uses. {exc}",
            )
        if subject.get("digest") != sequence["digest"]:
            _fail(
                "CF-OBS-002",
                "subject.digest does not match the sequence document of the explicit cds",
            )
    elif kind == "design":
        key = (str(subject.get("job_id")), str(subject.get("candidate_id")))
        known = None if candidates is None else candidates.get(key)
        if known is None:
            _fail(
                "CF-OBS-002",
                "referenced candidate is not available to check as_built_digest against",
            )
        if subject.get("as_built_digest") != known:
            _fail("CF-OBS-002", "as_built_digest does not match the referenced candidate")
    else:
        _fail("CF-OBS-002", f"subject.kind {kind!r} is not design or explicit")


def _schema(document: dict[str, Any]) -> None:
    errors = sorted(_validator().iter_errors(document), key=lambda item: list(item.path))
    if errors:
        first = errors[0]
        where = "/".join(str(part) for part in first.path)
        detail = first.message if not where else f"{where}: {first.message}"
        _fail("CF-OBS-001", detail)


def accept(
    document: dict[str, Any],
    *,
    candidates: dict[tuple[str, str], str] | None = None,
    translation_table: int = 11,
) -> dict[str, Any]:
    """Return the document to store, or raise `ObservationError`.

    `ingest_status` is replaced with `validated` and the stored digest is
    recomputed. A caller-supplied digest that disagrees with the submitted
    document is `CF-OBS-001`. The caller-supplied body is not mutated.
    """
    if not isinstance(document, dict):
        _fail("CF-OBS-001", "observation must be a JSON object")
    _semantic(document, candidates=candidates, translation_table=translation_table)
    supplied = document.get("digest")
    if supplied is not None and supplied != digest_of(document):
        _fail("CF-OBS-001", "digest does not match the submitted document")
    stored = {
        key: value for key, value in document.items() if key not in {"digest", "ingest_status"}
    }
    stored["ingest_status"] = {
        "state": "validated",
        "validation_errors": [],
        "incorporated_into": [],
    }
    stored["digest"] = digest_of(stored)
    _schema(stored)
    return stored


class ObservationLog:
    """Append-only files under the artifact store, one directory per tenant."""

    def __init__(self, root: Path) -> None:
        self.root = Path(root) / "observations"

    def append(self, tenant: str, document: dict[str, Any]) -> dict[str, Any]:
        digest = str(document["digest"])
        body = digest.removeprefix("sha256:")
        path = self.root / tenant / f"{body}.json"
        if path.is_file():
            existing = json.loads(path.read_text(encoding="utf-8"))
            if existing != document:
                _fail("CF-OBS-001", f"digest {digest} is already stored with different content")
            return existing
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(document, sort_keys=True), encoding="utf-8")
        return document

    def list(self, tenant: str) -> list[dict[str, Any]]:
        directory = self.root / tenant
        if not directory.is_dir():
            return []
        documents = [
            json.loads(path.read_text(encoding="utf-8"))
            for path in sorted(directory.glob("*.json"))
        ]
        documents.sort(key=lambda item: str(item.get("observation_id", "")))
        return documents
