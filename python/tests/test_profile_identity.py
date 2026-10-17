"""Pinned outputs for the v2.9.2 hot-path edits.

Seed repair stopped rebuilding a k-mer index on every pattern trial, and the
homopolymer excess became a local delta. Primer search stopped repeating
hairpin, self-dimer and Tm for every partner. Non-dominated sorting stopped
allocating a copy of each objective vector per comparison. None of those
edits may change a documented sequence, score, or primer document.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pytest

from codonforge import _core
from codonforge.modelreg import ArtifactStore
from codonforge.orchestrator.design import DesignRequest, design
from codonforge.predictors import available_engines

ROOT = Path(__file__).resolve().parents[2]
HOST = "escherichia_coli:511145/BL21"
DIGEST = "sha256:" + "ab" * 32


def _fasta(name: str) -> str:
    text = (ROOT / "benchmarks" / "validation" / "sequences" / name).read_text(encoding="utf-8")
    return "".join(
        line.strip() for line in text.splitlines() if line and not line.startswith(">")
    ).upper()


def _assembly() -> dict:
    return json.loads(
        (ROOT / "benchmarks" / "validation" / "constraints_assembly.json").read_text(
            encoding="utf-8"
        )
    )


def _profile() -> dict:
    return json.loads(
        (ROOT / "examples" / "mg1655_ribosomal_usage.json").read_text(encoding="utf-8")
    )


def test_waas_assembly_cai_vector_is_unchanged(tmp_path: Path) -> None:
    cds = _fasta("waaS.fa")
    result = design(
        DesignRequest(
            host=HOST,
            cds=cds,
            host_profile=_profile(),
            objectives=["usage.cai"],
            constraints=_assembly(),
            translation_table=11,
            budget="fast",
            seed=7,
            generations=0,
            stagnation_window=21,
            segment=2,
            host_digest=DIGEST,
        ),
        ArtifactStore(tmp_path / "artifacts"),
    )
    assert result.job.state == "succeeded", result.status.get("problem")
    assert result.front is not None
    vector = result.front["candidates"][0]["objective_vector"]
    assert vector == [0.9173459052341191]
    designed = result.front["candidates"][0]["sequences"][0]["cds"]
    assert _core.translate(designed, 11) == _core.translate(cds, 11)


def test_yeej_primer_document_is_byte_identical() -> None:
    sequence = _fasta("yeeJ.fa")
    flank = 300
    start = flank
    end = len(sequence) - flank
    region = end - start
    body = {
        "sequence": sequence,
        "region": {"start": start, "end": end},
        "config": {
            "tm_target": 60.0,
            "product_size_range": {"min": region, "max": region + 2 * flank},
            "flank_length": flank,
            "primer_length_range": {"min": 18, "max": 25},
        },
    }
    raw = _core.primer_search(json.dumps(body), lambda _fraction: None)
    assert (
        hashlib.sha256(raw.encode()).hexdigest()
        == "5f1fa115bcefde58d0c3d92a35e7adbc5fcbf04bf78b6efe43bf75f946aa149f"
    )
    recommended = json.loads(raw)["pairs"][0]
    assert recommended["product_size"] == 6843
    assert recommended["hetero_dimer_dG"] == -1.283850000000001
    assert recommended["refused"] is False
    assert recommended["forward_sequence"] == "TGACCGACAAATATTATG"
    assert recommended["reverse_sequence"] == "GAACTACAGTGGTTTAAA"


def test_yeej_structure_vector_and_fold_counts_are_unchanged(tmp_path: Path) -> None:
    """Tier 3 k-mer reuse must not move the documented structure front.

    `docs/benchmarks/v2.9.2-profiles.md`: yeeJ CAI + structure, assembly
    constraints, budget fast, generations 2, seed 7. The fold cache stays at
    1 full fold, 21 local edits, 12 hits.
    """
    if "viennarna" not in available_engines():
        if os.environ.get("CF_REQUIRE_STRUCTURE_ENGINE") == "1":
            pytest.fail("CF_REQUIRE_STRUCTURE_ENGINE=1 but ViennaRNA is not installed")
        pytest.skip("ViennaRNA is not installed; it is a user-installed backend")
    cds = _fasta("yeeJ.fa")
    result = design(
        DesignRequest(
            host=HOST,
            cds=cds,
            host_profile=_profile(),
            objectives=["usage.cai", "structure.local_accessibility"],
            constraints=_assembly(),
            translation_table=11,
            budget="fast",
            seed=7,
            generations=2,
            stagnation_window=21,
            segment=2,
            host_digest=DIGEST,
        ),
        ArtifactStore(tmp_path / "artifacts"),
    )
    assert result.job.state == "succeeded", result.status.get("problem")
    assert result.front is not None
    vector = result.front["candidates"][0]["objective_vector"]
    assert vector == [0.9575892496987667, 2.408954059933615]
    designed = result.front["candidates"][0]["sequences"][0]["cds"]
    assert _core.translate(designed, 11) == _core.translate(cds, 11)
    session = result.job._fold_session
    assert session is not None
    assert session.hits == 12
    assert session.full_folds == 1
    assert session.local_edits == 21
