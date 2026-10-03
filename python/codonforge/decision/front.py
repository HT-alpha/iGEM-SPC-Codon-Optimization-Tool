"""Pareto front emission.

Ticket P1-30, `cf-decision`'s half — "export to the ParetoFront schema" is this
module's declared responsibility, so the document is assembled here rather than
in whichever surface happens to want one.

Contract: `schemas/v1/pareto_front.schema.json`.

`DEP-CLI-003` records what a Phase 1 front was. A job still stores the full
non-dominated set with `reduction.method: "none"`. `P2-15` reduces a copy for
presentation (`reduce.present`); the stored document stays the one replay
compares.
"""

from __future__ import annotations

import datetime as _dt
import json
from typing import Any

from codonforge import _core

SCHEMA_VERSION = "1.0.0"


def build(
    job_id: str,
    front_id: str,
    lockfile_digest: str,
    objective_order: list[str],
    objective_sense: list[str],
    candidates: list[dict[str, Any]],
    budget_tier: str,
    termination_reason: str,
    n_evaluations: int = 0,
    input_sequence: dict[str, Any] | None = None,
    recommended: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Assemble and self-digest a front document.

    `objective_sense` is passed in rather than derived here, and its caller
    reads it from the providers. A sense guessed at this layer would disagree
    with the search that produced the numbers, and nothing in the document
    would reveal the disagreement.
    """
    if len(objective_sense) != len(objective_order):
        raise ValueError(
            f"{len(objective_order)} objectives but {len(objective_sense)} senses; "
            "a front whose senses do not line up with its objectives is unreadable"
        )

    document: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "front_id": front_id,
        "job_id": job_id,
        "created_at": _dt.datetime.now(_dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "landscape_lockfile_digest": lockfile_digest,
        "objective_order": objective_order,
        "objective_sense": objective_sense,
        "candidates": candidates,
        "reduction": {
            "raw_nondominated_count": len(candidates),
            "presented_count": len(candidates),
            "method": "none",
        },
        "summary": {
            "hypervolume": None,
            "hypervolume_normalized": None,
            "n_feasible": sum(1 for c in candidates if c["constraint_report"]["feasible"]),
            "n_evaluations": n_evaluations,
            "budget_tier": budget_tier,
            "termination_reason": termination_reason,
        },
        "digest": "sha256:" + "0" * 64,
    }
    # Which candidate is the answer, and which axes informed that choice.
    # Not recomputed here: the caller passes the outcome of
    # `CandidateArchive.best()` and the axis lists from the same
    # `scalarized_score_report` loop that produced its score. Deriving either
    # a second time at the serialization boundary is the pattern this project
    # has found five times.
    document["recommended"] = recommended
    if input_sequence is not None:
        # "Present so every improvement claim has a baseline in the same
        # document", per the schema. The baseline is the input candidate the
        # archive already scored and feasibility-checked, not a second
        # evaluation that could disagree with it.
        document["input_sequence"] = input_sequence
    # `created_at` is outside the digest by DEP-VER-003: a front holding the
    # same candidates from the same lockfile is the same front, and
    # GATE-REPRO-001 asks replay to reproduce one *exactly*.
    document["digest"] = _core.digest_json(json.dumps(document), "digest")
    return document


def candidate(
    candidate_id: str,
    cds: str,
    sequence_digest: str,
    objective_vector: list[float],
    origin: str,
    feasible: bool = True,
    generation_found: int | None = None,
    objective_detail: dict[str, dict[str, Any]] | None = None,
    unknown_terms: list[str] | None = None,
) -> dict[str, Any]:
    """One front entry.

    `sequence_digest` is the digest of the `sequence.schema.json` *document*,
    not of the bare `cds` string — the schema says so because
    `observation.subject.as_built_digest` is checked against it, and two
    conventions would make that check silently wrong for any construct with
    flanks.
    """
    entry = {
        "candidate_id": candidate_id,
        "cistron_id": None,
        "sequences": [{"cistron_id": None, "cds": cds, "digest": sequence_digest}],
        "nondomination_rank": 0,
        "objective_vector": objective_vector,
        "constraint_report": {
            "feasible": feasible,
            "n_blocking_violations": 0 if feasible else 1,
            "violations": [],
        },
        "provenance": {
            "origin": origin,
            "seed_source": None,
            "generation_found": generation_found,
            "parents": [],
        },
    }
    # Diagnostic terms live here rather than in `objective_vector`, which holds
    # the ranked axes. A term that did not drive the search must not appear in
    # the vector the front is ordered by, or a reader comparing candidates would
    # be reading a dimension nothing optimized. DEP-STRU-003.
    if objective_detail:
        entry["objective_detail"] = objective_detail
    # Always emitted, even when empty. An absent key would make "fully scored"
    # and "the writer forgot" indistinguishable, which is the ambiguity this
    # field exists to remove rather than relocate.
    entry["unknown_terms"] = list(unknown_terms or [])
    return entry
