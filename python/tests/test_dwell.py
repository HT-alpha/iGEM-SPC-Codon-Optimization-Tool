"""The dwell head is optional, and a missing fit is not a zero."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from codonforge.epistasis.artifact import ArtifactError
from codonforge.epistasis.dwell import dwell_status, verify_dwell
from codonforge.modelreg import ArtifactStore
from codonforge.orchestrator.design import DesignError, DesignRequest, design


def test_requesting_dwell_without_a_fitted_head_is_cf_land_001(tmp_path: Path) -> None:
    request = DesignRequest(
        host="escherichia_coli:511145/MG1655",
        cds="ATGGGGTTT",
        host_profile={},
        objectives=["usage.cai", "epistasis.dwell_deviation"],
    )
    with pytest.raises(DesignError) as caught:
        design(request, ArtifactStore(tmp_path))
    assert caught.value.code == "CF-LAND-001"
    assert "dwell" in caught.value.detail
    assert "zero" in caught.value.detail


def test_a_replaced_dwell_artifact_fails_the_digest() -> None:
    document = json.loads(Path("benchmarks/epistasis/dwell.json").read_text(encoding="utf-8"))
    verify_dwell(document)
    document["a_site_offset_nt"] = 15
    with pytest.raises(ArtifactError, match="digest"):
        verify_dwell(document)


def test_the_shipped_dwell_head_is_not_fitted() -> None:
    status = dwell_status()
    assert status["state"] == "not_fitted"
    assert status["accession"] == "GSE53767"
    assert status["a_site_offset_nt"] is None
    assert status["normalization"] is None
    titles = " ".join(item["paper"] for item in status["checked"])
    assert "Cell 2014" in titles
    assert "Nature 2012" in titles
    assert "Science 2009" in titles
    assert "eLife 2019" in titles
    assert "RiboToolkit" in titles
    assert "RiboCode" in titles
    assert "10 nt" in status["reason"]
    assert status["what_would_resolve"]
    assert dwell_status(Path("benchmarks/epistasis/does-not-exist.json"))["state"] == "absent"
