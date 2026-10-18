"""TST-CLI-013: derived ids stay inside the schema's Crockford class.

`_derived_id` is the only encoder for `job_`, `lock_` and `front_` ids. The
schema class is `[0-9A-HJKMNP-TV-Z]` — Crockford, which includes Y (`V-Z`) and
excludes I, L, O and U. A batch is generated and each id is matched against the
schema pattern loaded from disk, not against a restated copy of the regex.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from dna_design_suite.codon.orchestrator.design import _CROCKFORD, _derived_id

SCHEMAS = Path(__file__).resolve().parents[2] / "schemas" / "v1"
BATCH = 256


def _pattern(schema_name: str, *keys: str) -> re.Pattern[str]:
    document = json.loads((SCHEMAS / schema_name).read_text())
    node: object = document
    for key in keys:
        assert isinstance(node, dict), schema_name
        node = node[key]
    assert isinstance(node, str), f"{schema_name} {'/'.join(keys)} is not a pattern"
    return re.compile(node)


def test_tst_cli_013_derived_ids_match_the_schema_alphabet() -> None:
    front_pat = _pattern("pareto_front.schema.json", "properties", "front_id", "pattern")
    lock_pat = _pattern("landscape_lockfile.schema.json", "properties", "lockfile_id", "pattern")
    job_pat = _pattern("common.schema.json", "$defs", "JobId", "pattern")

    # The encoder's alphabet is the schema class. Checking this first is how an
    # injected I/L/O/U is caught even if a particular hash never draws it.
    schema_class = re.compile(r"[0-9A-HJKMNP-TV-Z]")
    assert len(_CROCKFORD) == 32
    assert len(set(_CROCKFORD)) == 32
    illegal = [c for c in _CROCKFORD if schema_class.fullmatch(c) is None]
    assert illegal == [], f"_CROCKFORD contains {illegal!r}, outside the schema class"

    saw_y = False
    for index in range(BATCH):
        front = _derived_id("front", f"front-material-{index}")
        lock = _derived_id("lock", f"lock-material-{index}")
        job = _derived_id("job", f"job-material-{index}")
        assert front_pat.fullmatch(front), front
        assert lock_pat.fullmatch(lock), lock
        assert job_pat.fullmatch(job), job
        for ident in (front, lock, job):
            body = ident.split("_", 1)[1]
            assert len(body) == 26
            for forbidden in "ILOU":
                assert forbidden not in body, ident
            saw_y = saw_y or ("Y" in body)

    # Y is in Crockford (`V-Z`). A "fix" that dropped it would still pass the
    # forbidden-character loop and fail here.
    assert saw_y, "a 256-id batch never emitted Y; the encoder may have dropped a legal symbol"
