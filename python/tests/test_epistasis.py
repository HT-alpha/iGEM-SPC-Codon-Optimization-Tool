"""Held-out split, artifact digest, and the calibration note."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dna_design_suite.codon.epistasis.artifact import (
    ArtifactError,
    load_calibration,
    verify_artifact,
)
from dna_design_suite.codon.epistasis.train import held_out_split
from dna_design_suite.modelreg.manifest import digest_of


def test_the_panel_genes_are_removed_before_the_fit() -> None:
    panel = json.loads((Path("benchmarks/corpus/panel.json")).read_text(encoding="utf-8"))
    held_out = list(panel["potts_held_out"])
    train, held = held_out_split(held_out)
    assert len(held_out) == 13
    assert len(held) == 13
    assert len(train) == 4127
    assert {row["gene"] for row in held} == set(held_out)
    assert not ({row["gene"] for row in train} & set(held_out))


def test_the_shipped_artifact_rejects_a_replaced_coupling() -> None:
    """The digest guard on the fitted manifest, not only on a toy document."""
    path = Path("benchmarks/epistasis/manifest.json")
    document = json.loads(path.read_text(encoding="utf-8"))
    verify_artifact(document)
    document["parameters"]["models"][1]["couplings"][0] += 1.0
    with pytest.raises(ArtifactError, match="digest"):
        verify_artifact(document)


def test_the_artifact_digest_covers_a_coupling(tmp_path: Path) -> None:
    document = {
        "schema_version": "1.0.0",
        "parameters": {"models": [{"window": 5, "couplings": [0.0, 1.25]}]},
    }
    document["digest"] = digest_of(document)
    verify_artifact(document)
    document["parameters"]["models"][0]["couplings"][1] = 9.0
    with pytest.raises(ArtifactError, match="digest"):
        verify_artifact(document)
    assert load_calibration(tmp_path / "missing.json") == {"state": "absent"}
