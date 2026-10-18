"""Limits that now run, and the ones that stay refused.

Each test is aimed at one failure mode. A short construct that is scored
without ViennaRNA would still satisfy a "the job succeeded" check, so the
structure test compares the accessibility penalty. A feasibility route that
returns five empty names would fail the infeasible fixture.
"""

from __future__ import annotations

import json

from test_cli import request_for
from test_ingest import GENBANK
from test_orchestrator import HOST, constraint_set, host_profile

from dna_design_suite import _core
from dna_design_suite.api.server import (
    _prepare_job_body,
    _reverse_translate,
    analyse_feasibility,
    list_models,
)
from dna_design_suite.codon.orchestrator import DesignError, DesignRequest, design
from dna_design_suite.codon.orchestrator.structure_diag import TERM_ID, FoldSession
from dna_design_suite.modelreg import ArtifactStore


def test_a_short_construct_is_polished_for_accessibility(tmp_path) -> None:
    """Longer than one span keeps CF-SRCH-010. A construct inside the span does not."""
    store = ArtifactStore(tmp_path / "store")
    short = "ATG" + "CTG" * 8 + "TAA"
    assert len(short) <= 40
    structured = design(
        DesignRequest(
            host=HOST,
            cds=short,
            host_profile=host_profile(),
            objectives=[TERM_ID],
            constraints={
                "schema_version": "1.0.0",
                "constraint_set_id": "cst_none",
                "constraints": [],
            },
            translation_table=11,
            budget="fast",
            seed=7,
            generations=0,
            host_digest="sha256:" + "1" * 64,
        ),
        store,
    )
    assert structured.job.state == "succeeded"
    assert not [w for w in structured.status["warnings"] if w["code"] == "CF-SRCH-010"]
    cai = design(
        DesignRequest(
            host=HOST,
            cds=short,
            host_profile=host_profile(),
            objectives=["usage.cai"],
            constraints={
                "schema_version": "1.0.0",
                "constraint_set_id": "cst_none",
                "constraints": [],
            },
            translation_table=11,
            budget="fast",
            seed=7,
            generations=0,
            host_digest="sha256:" + "1" * 64,
        ),
        ArtifactStore(tmp_path / "cai"),
    )
    session = FoldSession()
    structured_cds = structured.front["candidates"][0]["sequences"][0]["cds"]
    cai_cds = cai.front["candidates"][0]["sequences"][0]["cds"]
    input_penalty, structured_penalty, cai_penalty = session.evaluate(
        [short, structured_cds, cai_cds]
    )
    assert structured_penalty < input_penalty
    assert structured_penalty <= cai_penalty
    assert not [w for w in cai.status["warnings"] if w["code"] == "CF-LAND-013"]

    long = design(request_for(budget="fast", generations=0), store)
    assert any(
        w["code"] == "CF-SRCH-010" and w["subject"] == TERM_ID for w in long.status["warnings"]
    )


def test_feasibility_runs_five_checks_and_an_impossible_band_fires() -> None:
    status, empty = analyse_feasibility({})
    assert status == 400
    assert empty["code"] == "CF-CLI-001"
    status, report = analyse_feasibility(
        {
            "cds": "ATGAAATAA",
            "constraints": constraint_set(),
            "translation_table": 11,
        }
    )
    assert status == 200
    assert report["checks_run"] == [
        "locked_span",
        "counting_bound",
        "per_position",
        "automaton_reachability",
        "gc_attainability",
    ]
    status, refused = analyse_feasibility(
        {
            "cds": "ATGAAATAA",
            "translation_table": 11,
            "constraints": {
                "schema_version": "1.0.0",
                "constraint_set_id": "cst_impossible_gc",
                "constraints": [
                    {
                        "id": "c_gc",
                        "class": "blocking",
                        "scope": "cds",
                        "rule": {"kind": "gc_window", "window_nt": 9, "min": 0.9, "max": 1.0},
                    }
                ],
            },
        }
    )
    assert status == 200
    assert refused["feasible"] is False
    assert any(item["check"] == "gc_attainability" for item in refused["findings"])


def test_protein_and_construct_inputs_use_the_caller_data() -> None:
    profile = {
        "codon_usage": {
            "subset": "user_supplied",
            "table": {
                "ATG": {"relative_adaptiveness": 1.0},
                "AAA": {"relative_adaptiveness": 1.0},
                "AAG": {"relative_adaptiveness": 0.0},
                "ACT": {"relative_adaptiveness": 1.0},
                "TAA": {"relative_adaptiveness": 1.0},
            },
        }
    }
    cds = _reverse_translate("MKT", 11, profile)
    assert cds == "ATGAAAACTTAA"
    assert _core.translate(cds, 11) == "MKT*"
    prepared = _prepare_job_body(
        {
            "host": HOST,
            "host_profile": profile,
            "input": {"kind": "protein", "protein": "MKT"},
            "objectives_preset": "cai",
            "translation_table": 11,
            "generations": 0,
            "budget": "fast",
        }
    )
    assert prepared["cds"] == cds
    assert prepared["objectives"] == ["usage.cai"]
    unknown = _prepare_job_body({"objectives_preset": "invented", "host": HOST, "cds": "ATGAAATAA"})
    assert unknown[0] == 422
    assert unknown[1]["code"] == "CF-CLI-001"
    construct = _prepare_job_body(
        {
            "host": HOST,
            "host_profile": profile,
            "input": {"kind": "construct", "genbank": GENBANK, "locus": "reporter"},
        }
    )
    assert construct["cds"] == "ATGAAACGTATTCTGTAA"
    single = _prepare_job_body(
        {
            "host": HOST,
            "host_profile": profile,
            "input": {"kind": "single_cds", "sequence": "ATGAAATAA"},
        }
    )
    assert single["cds"] == "ATGAAATAA"


def test_dwell_stays_refused_and_models_are_the_verified_file(tmp_path) -> None:
    store = ArtifactStore(tmp_path / "store")
    try:
        design(
            request_for(objectives=["usage.cai", "epistasis.dwell_deviation"], generations=0),
            store,
        )
    except DesignError as exc:
        assert exc.code == "CF-LAND-001"
    else:
        raise AssertionError("dwell was scored")
    models = list_models()
    if models["installed"]:
        assert models["models"][0]["id"] == "epistasis.context_naturalness"
        assert models["models"][0]["digest"]
    else:
        assert models["models"] == []


def test_a_preset_names_only_shipped_terms() -> None:
    body = _prepare_job_body(
        {
            "host": HOST,
            "cds": "ATGAAATAA",
            "host_profile": host_profile(),
            "objectives_preset": "cai_burden",
        }
    )
    assert body["objectives"] == ["usage.cai", "burden.growth_cost"]
    raw = json.dumps(body["objectives"])
    assert "s_value" not in raw
