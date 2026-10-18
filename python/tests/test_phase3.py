"""Phase 3: context naturalness from the verified Potts artifact."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from test_orchestrator import CDS, FIVE_PRIME, HOST, THREE_PRIME, host_profile

from dna_design_suite import _core
from dna_design_suite.codon.epistasis.artifact import artifact_path
from dna_design_suite.codon.epistasis.context import TERM_ID
from dna_design_suite.codon.orchestrator import DesignError, DesignRequest, design
from dna_design_suite.modelreg import ArtifactStore

ROOT = Path(__file__).resolve().parents[2]


def _request(**overrides: object) -> DesignRequest:
    base = dict(
        host=HOST,
        cds=CDS,
        host_profile=host_profile(),
        objectives=[TERM_ID],
        constraints={
            "schema_version": "1.0.0",
            "constraint_set_id": "cst_phase3",
            "constraints": [],
        },
        five_prime_flank=FIVE_PRIME,
        three_prime_flank=THREE_PRIME,
        translation_table=1,
        budget="fast",
        seed=7,
        generations=0,
        host_digest="sha256:" + "1" * 64,
        job_id="job_01JBQ8Z4K7M2N3P5R6T7V8W9XY",
    )
    base.update(overrides)
    return DesignRequest(**base)  # type: ignore[arg-type]


def test_a_missing_potts_artifact_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "dna_design_suite.codon.epistasis.context.artifact_path",
        lambda: tmp_path / "missing-manifest.json",
    )
    with pytest.raises(DesignError) as excinfo:
        design(_request(), ArtifactStore(tmp_path / "artifacts"))
    assert excinfo.value.code == "CF-LAND-001"
    assert "zero" in excinfo.value.detail


def test_a_non_coli_host_does_not_borrow_the_coli_artifact(tmp_path: Path) -> None:
    with pytest.raises(DesignError) as excinfo:
        design(
            _request(host="saccharomyces_cerevisiae:4932/S288C"),
            ArtifactStore(tmp_path / "artifacts"),
        )
    assert excinfo.value.code == "CF-LAND-001"


def test_dwell_stays_refused(tmp_path: Path) -> None:
    with pytest.raises(DesignError) as excinfo:
        design(
            _request(objectives=["epistasis.dwell_deviation"]),
            ArtifactStore(tmp_path / "artifacts"),
        )
    assert excinfo.value.code == "CF-LAND-001"


def test_the_installed_artifact_is_the_objective(tmp_path: Path) -> None:
    result = design(_request(), ArtifactStore(tmp_path / "artifacts"))
    assert result.job.state == "succeeded", result.status.get("problem")
    front = result.front
    assert front["objective_sense"] == ["maximize"]
    assert front["objective_order"] == [TERM_ID]
    vector = front["candidates"][0]["objective_vector"]
    assert vector[0] is not None
    designed = front["candidates"][0]["sequences"][0]["cds"]
    document = json.loads(artifact_path().read_text(encoding="utf-8"))
    window = int(document["parameters"]["selected_window"])
    chosen = next(row for row in document["parameters"]["models"] if int(row["window"]) == window)
    scored = json.loads(
        _core.score_potts_context(
            json.dumps(
                {
                    "sequences": [designed],
                    "window": window,
                    "h": chosen["h"],
                    "couplings": chosen["couplings"],
                }
            )
        )
    )
    assert vector[0] == pytest.approx(scored["nats_per_codon"])
    pinned = next(
        entry for entry in result.job.lockfile.resolution if entry.term_id == TERM_ID
    )
    assert pinned.digest == document["digest"]
    assert pinned.resolved_by == "pin"
    assert result.job.inputs.tier2_block >= window
