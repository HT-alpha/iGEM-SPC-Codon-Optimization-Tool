"""TST-LAPI-001, Python leg: the two runtimes agree on the landscape contract.

Ticket P1-20.

ADR-0001 accepts duplicated types across the Rust/Python boundary as the cost of
a Rust hot path, and names contract tests as the mitigation. The Rust half of
this test (`crates/cf-landscape-api/tests/contract.rs`) shows a document
surviving Rust -> JSON -> Rust -> JSON. This half closes the loop: the same
documents validate against the shipped schema, and canonicalizing them in Python
reproduces the bytes Rust produced.

A drift that this catches and a type check would not: a field one side treats as
optional-with-a-default and the other treats as absent. Both sides parse; the
canonical forms differ by one key.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from codonforge import _core

SCHEMA_DIR = Path(__file__).resolve().parents[2] / "schemas" / "v1"

HOST_PROFILE_DIGEST = "sha256:" + "1" * 64
BASE_DIGEST = "sha256:" + "2" * 64
CONSTRUCT_DIGEST = "sha256:" + "3" * 64
ARTIFACT_DIGEST = "sha256:" + "4" * 64

HOST_CONTEXT = {
    "host_id": "escherichia_coli:511145/BL21",
    "host_profile_digest": HOST_PROFILE_DIGEST,
    "condition": {"condition_id": "cond_m9_glucose_37c"},
}


def _query(**overrides: object) -> dict:
    document = {
        "schema_version": "1.0.0",
        "contract_major": 1,
        "query_id": "q_01ARZ3NDEKTSV4RRFFQ69G5FAV",
        "target": {
            "kind": "full",
            "sequence": json.loads(
                _core.sequence_document("ATGCTGCGTAGCGGCGTGTAA", 1, None, None, None)
            ),
        },
        "host_context": HOST_CONTEXT,
        "requested_terms": ["usage.cai", "structure.local_accessibility"],
        "budget": {"tier": "standard", "max_latency_ms": 250},
    }
    document.update(overrides)
    return document


def _response(**overrides: object) -> dict:
    document = {
        "schema_version": "1.0.0",
        "contract_major": 1,
        "query_id": "q_01ARZ3NDEKTSV4RRFFQ69G5FAV",
        "status": "ok",
        "terms": {
            "usage.cai": {
                "value": 0.7421,
                "unit": "dimensionless",
                "direction": "maximize",
                "capability": "L1",
                "provenance": {
                    "resolved_by": "builtin",
                    "engine": "cai_native",
                    "engine_version": "0.1.0",
                },
                "uncertainty": {
                    "method": "none",
                    "calibrated": True,
                    "sigma": 0.0,
                    "shared_fraction": 0.0,
                },
            }
        },
        "objective_vector": {
            "order": ["usage.cai"],
            "values": [0.7421],
            "sense": ["maximize"],
        },
        "diagnostics": {"cache_hits": 1, "cache_misses": 0, "evaluations": 1},
    }
    document.update(overrides)
    return document


QUERIES = {
    "full": _query(),
    "delta": _query(
        query_id="q_01ARZ3NDEKTSV4RRFFQ69G5FAW",
        target={
            "kind": "delta",
            "base_digest": BASE_DIGEST,
            "edits": [{"codon_index": 1, "codon": "CTC"}],
        },
        requested_terms=["usage.cai"],
        budget={"tier": "fast"},
    ),
    "delta_table": _query(
        query_id="q_01ARZ3NDEKTSV4RRFFQ69G5FAX",
        target={
            "kind": "delta_table",
            "base_digest": BASE_DIGEST,
            "positions": [1, 2, 3, 4],
        },
        outputs={"scalars": True, "deltas": True},
        construct_context={"kind": "ref", "construct_digest": CONSTRUCT_DIGEST},
        model_pins={
            "epistasis.dwell_deviation": {
                "artifact_id": "cf:epistasis_dwell:escherichia_coli:511145/BL21:1.2.0",
                "digest": ARTIFACT_DIGEST,
            }
        },
        requested_terms=["epistasis.dwell_deviation", "structure.local_accessibility"],
        budget={"tier": "exact"},
        seed=20260826,
    ),
}

RESPONSES = {
    "ok": _response(),
    "partial": _response(
        status="partial",
        objective_vector={
            "order": ["usage.cai", "structure.local_accessibility"],
            "values": [0.7421, None],
            "sense": ["maximize", "minimize"],
        },
        diagnostics={
            "cache_hits": 0,
            "cache_misses": 2,
            "evaluations": 1,
            "terms_deferred": ["structure.local_accessibility"],
            "engine_notes": ["deferred by latency budget"],
        },
    ),
}


def _validate(document: dict, schema_name: str) -> None:
    jsonschema = pytest.importorskip("jsonschema")
    from referencing import Registry, Resource

    resources = [
        (path.name, Resource.from_contents(json.loads(path.read_text())))
        for path in SCHEMA_DIR.glob("*.schema.json")
    ]
    registry = Registry().with_resources(resources)
    schema = json.loads((SCHEMA_DIR / schema_name).read_text())
    jsonschema.Draft202012Validator(schema, registry=registry).validate(document)


@pytest.mark.parametrize("name", sorted(QUERIES))
def test_tst_lapi_001_query_round_trips_across_the_boundary(name: str) -> None:
    document = QUERIES[name]
    _validate(document, "landscape_query.schema.json")

    once = _core.landscape_query_round_trip(json.dumps(document))
    twice = _core.landscape_query_round_trip(once)
    assert once == twice, f"{name}: the round trip is not a fixed point"

    # The canonical form Python computes over the parsed document must be the
    # canonical form Rust emitted. If Rust dropped or invented a field, these
    # differ even though both sides parsed without complaint.
    assert _core.canonical_json(once) == once


@pytest.mark.parametrize("name", sorted(RESPONSES))
def test_tst_lapi_001_response_round_trips_across_the_boundary(name: str) -> None:
    document = RESPONSES[name]
    _validate(document, "landscape_response.schema.json")

    once = _core.landscape_response_round_trip(json.dumps(document))
    twice = _core.landscape_response_round_trip(once)
    assert once == twice
    assert _core.canonical_json(once) == once


def test_tst_lapi_001_construct_context_survives_without_being_interpreted() -> None:
    """A block no Phase 1 term reads is preserved, not dropped."""
    document = QUERIES["delta_table"]
    emitted = json.loads(_core.landscape_query_round_trip(json.dumps(document)))
    assert emitted["construct_context"] == {
        "kind": "ref",
        "construct_digest": CONSTRUCT_DIGEST,
    }


def test_tst_lapi_001_gradient_request_is_refused_at_the_boundary() -> None:
    """ADR-0002 closed the L2 path; asking for it fails rather than being ignored."""
    document = _query(
        outputs={"scalars": True, "gradients": True},
    )
    # Schema-valid: the contract admits the field so a non-CodonForge consumer
    # can use it. The Rust contract still refuses, which is the point.
    _validate(document, "landscape_query.schema.json")

    with pytest.raises(ValueError) as excinfo:
        _core.landscape_query_round_trip(json.dumps(document))
    problem = json.loads(str(excinfo.value))
    assert problem["code"] == "CF-LAND-002"


def test_tst_lapi_001_contract_major_mismatch_is_refused() -> None:
    document = _query(contract_major=2)
    with pytest.raises(ValueError) as excinfo:
        _core.landscape_query_round_trip(json.dumps(document))
    problem = json.loads(str(excinfo.value))
    assert problem["code"] == "CF-RESOLVE-002"


def test_tst_lapi_001_unknown_shared_fraction_is_not_defaulted() -> None:
    """An omitted shared_fraction stays omitted; it does not become 0.

    ADR-0005 rejects the independence assumption explicitly, and 0 is that
    assumption. A round trip that materialized it would make an unmeasured
    decomposition indistinguishable from a measured one.
    """
    document = _response(
        terms={
            "lm.naturalness": {
                "value": -2.7,
                "unit": "nats/codon",
                "direction": "maximize",
                "capability": "L0",
                "provenance": {"resolved_by": "base"},
                "uncertainty": {
                    "method": "deep_ensemble",
                    "calibrated": True,
                    "sigma": 0.31,
                    "n_members": 5,
                },
            }
        },
        objective_vector={
            "order": ["lm.naturalness"],
            "values": [-2.7],
            "sense": ["maximize"],
        },
    )
    emitted = json.loads(_core.landscape_response_round_trip(json.dumps(document)))
    assert "shared_fraction" not in emitted["terms"]["lm.naturalness"]["uncertainty"]
