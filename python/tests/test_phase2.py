"""Phase 2: burden, structure on the recommendation, and local_structure."""

from __future__ import annotations

from pathlib import Path

import pytest
from test_orchestrator import CDS, FIVE_PRIME, HOST, THREE_PRIME, host_profile

from codonforge.modelreg import ArtifactStore
from codonforge.orchestrator import (
    DesignError,
    DesignRequest,
    design,
    local_structure,
    structure_diag,
)

ROOT = Path(__file__).resolve().parents[2]


def _request(**overrides: object) -> DesignRequest:
    base = dict(
        host=HOST,
        cds=CDS,
        host_profile=host_profile(),
        objectives=["usage.cai"],
        constraints={
            "schema_version": "1.0.0",
            "constraint_set_id": "cst_phase2",
            "constraints": [],
        },
        five_prime_flank=FIVE_PRIME,
        three_prime_flank=THREE_PRIME,
        translation_table=1,
        budget="fast",
        seed=7,
        generations=2,
        host_digest="sha256:" + "1" * 64,
        job_id="job_01JBQ8Z4K7M2N3P5R6T7V8W9XY",
    )
    base.update(overrides)
    return DesignRequest(**base)  # type: ignore[arg-type]


def test_structure_without_an_engine_is_refused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(structure_diag, "is_available", lambda: False)
    with pytest.raises(DesignError) as excinfo:
        design(
            _request(objectives=["usage.cai", "structure.local_accessibility"]),
            ArtifactStore(tmp_path / "artifacts"),
        )
    assert excinfo.value.code == "CF-STRU-003"


def test_cai_only_still_runs_when_structure_is_requested_nowhere(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(structure_diag, "is_available", lambda: False)
    result = design(_request(), ArtifactStore(tmp_path / "artifacts"))
    assert result.job.state == "succeeded", result.status.get("problem")
    assert result.front["recommended"]["rule"] == "feasible_then_score_then_digest"
    assert result.front["recommended"]["informed_by"] == ["usage.cai"]


def test_local_structure_is_measured_and_min_codon_frequency_is_not(tmp_path: Path) -> None:
    cds = "ATGAAACCCGGGTTTAAATAA"
    constraints = {
        "schema_version": "1.0.0",
        "constraint_set_id": "cst_local",
        "constraints": [
            {
                "id": "c_struct",
                "class": "blocking",
                "scope": "cds",
                "rule": {
                    "kind": "local_structure",
                    "region": "custom",
                    "span": {"start_nt": 0, "end_nt": 12, "strand": "+"},
                    "metric": "mean_unpaired_probability",
                    "min": 0.0,
                    "max": 1.0,
                    "engine": "native_window",
                },
            },
            {
                "id": "c_freq",
                "class": "blocking",
                "scope": "cds",
                "rule": {"kind": "min_codon_frequency", "min_relative_adaptiveness": 0.1},
            },
        ],
    }
    result = design(
        _request(cds=cds, constraints=constraints, five_prime_flank="", three_prime_flank=""),
        ArtifactStore(tmp_path / "artifacts"),
    )
    assert result.job.state == "succeeded", result.status.get("problem")
    warnings = result.status["warnings"]
    assert not any(
        warning["code"] == "CF-FEAS-008" and warning.get("subject") == "c_struct"
        for warning in warnings
    )
    assert any(
        warning["code"] == "CF-FEAS-008" and warning.get("subject") == "c_freq"
        for warning in warnings
    )
    designed = result.front["candidates"][0]["sequences"][0]["cds"]
    rule = constraints["constraints"][0]["rule"]
    measured = local_structure.measure(designed, rule, {})
    from codonforge.predictors.viennarna import ViennaRNAEngine

    direct = ViennaRNAEngine().local_accessibility(
        designed, window_nt=1, span_nt=1, temperature_c=37.0
    )
    direct_mean = sum(direct[:12]) / 12
    assert measured == pytest.approx(direct_mean)


def test_a_local_structure_rule_without_a_span_is_refused(tmp_path: Path) -> None:
    constraints = {
        "schema_version": "1.0.0",
        "constraint_set_id": "cst_nospan",
        "constraints": [
            {
                "id": "c_struct",
                "class": "blocking",
                "rule": {
                    "kind": "local_structure",
                    "region": "initiation_window",
                    "metric": "mean_unpaired_probability",
                    "min": 0.35,
                },
            }
        ],
    }
    with pytest.raises(DesignError) as excinfo:
        design(_request(constraints=constraints), ArtifactStore(tmp_path / "artifacts"))
    assert excinfo.value.code == "CF-STRU-004"


def test_cai_and_burden_produce_a_complete_front(tmp_path: Path) -> None:
    result = design(
        _request(
            objectives=["usage.cai", "burden.growth_cost"],
            budget="standard",
            generations=4,
        ),
        ArtifactStore(tmp_path / "artifacts"),
    )
    assert result.job.state == "succeeded", result.status.get("problem")
    front = result.front
    assert front["objective_sense"] == ["maximize", "minimize"]
    complete = [
        candidate["objective_vector"]
        for candidate in front["candidates"]
        if all(value is not None for value in candidate["objective_vector"])
    ]
    assert complete
    assert front["recommended"]["rule"] == "knee_of_complete_front"
    assert front["recommended"]["informed_by"] == ["usage.cai", "burden.growth_cost"]
    assert front["recommended"]["not_informed_by"] == []
    assert front["candidates"][0]["candidate_id"] == front["recommended"]["candidate_id"]
    for left in complete:
        for right in complete:
            if left is right:
                continue
            assert not (left[0] >= right[0] and left[1] <= right[1] and left != right)


def test_structure_stays_on_its_own_axis_when_burden_is_also_ranked(tmp_path: Path) -> None:
    """Tier 3 appends the L0 axis. The document still reads request order.

    With CAI, structure, then burden, the engine vector is
    [CAI, burden, structure]. Reading it as the request order reports the
    burden value twice and drops the fold.
    """
    result = design(
        _request(
            objectives=[
                "usage.cai",
                "structure.local_accessibility",
                "burden.growth_cost",
            ],
            budget="standard",
            generations=2,
        ),
        ArtifactStore(tmp_path / "artifacts"),
    )
    assert result.job.state == "succeeded", result.status.get("problem")
    front = result.front
    assert front["objective_order"] == [
        "usage.cai",
        "structure.local_accessibility",
        "burden.growth_cost",
    ]
    folded = {
        candidate["sequences"][0]["cds"]: value
        for candidate, value in zip(
            front["candidates"],
            structure_diag.evaluate(
                [candidate["sequences"][0]["cds"] for candidate in front["candidates"]]
            ),
            strict=True,
        )
    }
    filled = []
    for candidate in front["candidates"]:
        vector = candidate["objective_vector"]
        fold = folded[candidate["sequences"][0]["cds"]]
        # The batch budget can leave an individual unknown on this axis.
        if vector[1] is None:
            continue
        assert vector[1] == pytest.approx(fold)
        assert vector[2] != pytest.approx(vector[1])
        filled.append(vector[1])
    assert filled
    summary = result.status["current_front_summary"]["best_per_objective"]
    assert summary["structure.local_accessibility"] == pytest.approx(min(filled))


def test_explain_prints_a_tradeoff_from_the_front(tmp_path: Path) -> None:
    from codonforge.cli import writers

    result = design(
        _request(
            objectives=["usage.cai", "burden.growth_cost"],
            budget="standard",
            generations=4,
        ),
        ArtifactStore(tmp_path / "artifacts"),
    )
    rendered = writers.explain(result.front, result.status)
    if len(result.front["candidates"]) >= 2:
        assert "trade-off" in rendered
        assert "burden.growth_cost" in rendered
