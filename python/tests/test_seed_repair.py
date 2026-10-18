"""DEP-SRCH-015: a tier-1 refusal is repaired, or it fails honestly.

The elastin-like CDS trips GC windows and repeated 12-mers, which the tier-1
DP cannot carry. Returning that CDS with exit 0 is the v2.5.0 failure. An
impossible GC band must exit 3 and must not print the input as the design.
"""

from __future__ import annotations

import json
import subprocess
import sys
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

from test_orchestrator import HOST, host_profile

from dna_design_suite.cli import main as cli_main
from dna_design_suite.cli import writers
from dna_design_suite.codon.orchestrator import DesignRequest, design
from dna_design_suite.codon.orchestrator.rationale import build_rationale
from dna_design_suite.modelreg import ArtifactStore

ROOT = Path(__file__).resolve().parents[2]


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


def test_seed_repair_returns_a_feasible_elastin_design(tmp_path: Path) -> None:
    import importlib

    sys.path.insert(0, str(ROOT / "benchmarks" / "validation"))
    cai_independent = importlib.import_module("cai_independent")
    ribosomal = json.loads(
        (ROOT / "examples" / "mg1655_ribosomal_usage.json").read_text(encoding="utf-8")
    )
    weights = cai_independent.load_ribosomal_weights(
        ROOT / "examples" / "mg1655_ribosomal_usage.json"
    )
    cds = _fasta("elp_vpgvg.fa")
    result = design(
        DesignRequest(
            host=HOST,
            cds=cds,
            host_profile=ribosomal,
            objectives=["usage.cai"],
            constraints=_assembly(),
            translation_table=11,
            budget="fast",
            seed=7,
            generations=0,
            host_digest="sha256:" + "1" * 64,
            job_id="job_01JBQ8Z4K7M2N3P5R6T7V8W9XY",
        ),
        ArtifactStore(tmp_path / "artifacts"),
    )
    assert result.job.state == "succeeded", result.status.get("problem")
    assert result.front is not None
    designed = result.front["candidates"][0]["sequences"][0]["cds"]
    assert designed != cds
    assert result.front["candidates"][0]["constraint_report"]["feasible"]
    assert cai_independent.translate(designed) == cai_independent.translate(cds)
    measured = float(cai_independent.score(designed, weights)["cai"])
    assert measured > 0.369953, measured
    assert any(warning["code"] == "CF-SRCH-017" for warning in result.status["warnings"])
    states = {tier["tier"]: tier["state"] for tier in result.status["progress"]["tiers"]}
    assert states["tier1_seed"] == "done"
    assert states["tier2_block_ca"] == "done"


def test_an_impossible_gc_band_exits_3_and_does_not_print_the_input(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    cds = "ATG" + "GGG" * 12 + "TAA"
    constraints = {
        "schema_version": "1.0.0",
        "constraint_set_id": "cst_impossible_gc",
        "constraints": [
            {
                "id": "c_gc",
                "class": "blocking",
                "scope": "cds",
                "rule": {"kind": "gc_window", "window_nt": 12, "min": 0.0, "max": 0.1},
            }
        ],
    }
    profile = tmp_path / "host.json"
    profile.write_text(json.dumps(host_profile()), encoding="utf-8")
    rules = tmp_path / "constraints.json"
    rules.write_text(json.dumps(constraints), encoding="utf-8")
    code = cli_main(
        [
            "design",
            "--host",
            HOST,
            "--cds",
            cds,
            "--host-profile",
            str(profile),
            "--constraints",
            str(rules),
            "--budget",
            "fast",
            "--generations",
            "1",
            "--store",
            str(tmp_path / "artifacts"),
        ]
    )
    captured = capsys.readouterr()
    assert code == 3
    assert "CF-SRCH-006" in captured.err
    assert cds not in captured.out
    assert "tier1_seed" in captured.out
    assert "pending" not in captured.out


def _empty_constraints() -> dict:
    return {
        "schema_version": "1.0.0",
        "constraint_set_id": "cst_none",
        "constraints": [],
    }


def _design(tmp_path: Path, cds: str, profile: dict, constraints: dict):
    return design(
        DesignRequest(
            host=HOST,
            cds=cds,
            host_profile=profile,
            objectives=["usage.cai"],
            constraints=constraints,
            translation_table=11,
            budget="fast",
            seed=7,
            generations=2,
            host_digest="sha256:" + "1" * 64,
            job_id="job_01JBQ8Z4K7M2N3P5R6T7V8W9XY",
        ),
        ArtifactStore(tmp_path / "artifacts"),
    )


def test_explain_states_what_cai_means_and_names_a_zero_weight(tmp_path: Path) -> None:
    profile = json.loads(
        (ROOT / "examples" / "mg1655_ribosomal_usage.json").read_text(encoding="utf-8")
    )
    cds = "ATGAGGCCCTAA"
    result = _design(tmp_path, cds, profile, _empty_constraints())
    assert result.job.state == "succeeded", result.status.get("problem")
    assert result.front is not None
    explained = writers.explain(result.front, result.status)
    assert "1.0 means every counted codon is a top synonym of that table" in explained
    assert "Stops and single-codon families are excluded" in explained
    assert "Kazusa" in explained
    reading = build_rationale(
        _request_from(result, profile, cds),
        cds,
        result.job._tier_request(cds),
        {"usage.cai": 0.0},
    )
    assert reading is not None
    by_tag = {item["tag"]: item["text"] for item in reading["annotations"]}
    assert (
        "1.0 means every counted codon is a top synonym of that table"
        in by_tag["cai_table_relative"]
    )
    assert "Stops and single-codon families are excluded" in by_tag["cai_table_relative"]
    assert "Codon AGG at position 1" in by_tag["cai_zero_weight"]
    assert "0.002392" in by_tag["cai_zero_weight"]
    assert "Kazusa" in by_tag["cai_reference_table"]


def _request_from(result, profile: dict, cds: str) -> DesignRequest:
    """The request the rationale reader expects. The job already ran it."""
    return DesignRequest(
        host=HOST,
        cds=cds,
        host_profile=profile,
        objectives=["usage.cai"],
        translation_table=result.job.inputs.translation_table,
    )


def test_biobrick_preset_removes_an_ecori_site(tmp_path: Path) -> None:
    preset = json.loads((ROOT / "examples" / "biobrick_rfc10.json").read_text(encoding="utf-8"))
    cds = "ATGGAATTCAAATAA"
    result = _design(tmp_path, cds, host_profile(), preset)
    assert result.job.state == "succeeded", result.status.get("problem")
    assert result.front is not None
    designed = result.front["candidates"][0]["sequences"][0]["cds"]
    assert designed != cds
    assert result.front["candidates"][0]["constraint_report"]["feasible"]
    assert "GAATTC" not in designed


def test_cli_module_execution_does_not_warn() -> None:
    completed = subprocess.run(
        [sys.executable, "-W", "error", "-m", "dna_design_suite.cli.main", "--help"],
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert "RuntimeWarning" not in completed.stderr


def test_importing_main_from_the_package_still_works() -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        from dna_design_suite.cli import main

        assert callable(main)


def _corpus_cds(gene: str) -> str:
    import gzip

    found = ""
    current = None
    buf: list[str] = []
    path = ROOT / "benchmarks" / "corpus" / "NC_000913.3_cds.fna.gz"
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if current == gene:
                    found = "".join(buf)
                    break
                current = None
                if "[gene=" in line:
                    current = line.split("[gene=", 1)[1].split("]", 1)[0]
                buf = []
            else:
                buf.append(line.strip())
        else:
            if current == gene:
                found = "".join(buf)
    if not found:
        raise AssertionError(f"{gene} missing from the corpus CDS")
    return found.upper()


def test_rpoc_assembly_seed_is_feasible(tmp_path: Path) -> None:
    """rpoC assembly used to stop with GC windows still open.

    Single-codon sweeps left five 50 nt windows outside the band, and k-mer
    repair refuses to edit while that penalty remains, so the job was
    ``CF-SRCH-006``. A synonymous CDS exists. ``DEP-SRCH-023`` tries one pair
    inside the stuck window, at most four times, and only after a cycle that
    made no progress. Skipping that branch fails this test: the state goes
    back to failed and the problem code is ``CF-SRCH-006``.
    """
    import importlib

    sys.path.insert(0, str(ROOT / "benchmarks" / "validation"))
    cai_independent = importlib.import_module("cai_independent")
    cds = _corpus_cds("rpoC")
    assert cds.startswith("GTG")
    ribosomal = json.loads(
        (ROOT / "examples" / "mg1655_ribosomal_usage.json").read_text(encoding="utf-8")
    )
    weights = cai_independent.load_ribosomal_weights(
        ROOT / "examples" / "mg1655_ribosomal_usage.json"
    )
    result = design(
        DesignRequest(
            host=HOST,
            cds=cds,
            host_profile=ribosomal,
            objectives=["usage.cai"],
            constraints=_assembly(),
            translation_table=11,
            budget="fast",
            seed=7,
            generations=0,
            stagnation_window=21,
            segment=2,
            host_digest="sha256:" + "ab" * 32,
        ),
        ArtifactStore(tmp_path / "artifacts"),
    )
    assert result.job.state == "succeeded", result.status.get("problem")
    assert result.front is not None
    designed = result.front["candidates"][0]["sequences"][0]["cds"]
    assert cai_independent.translate(designed) == cai_independent.translate(cds)
    assert result.front["candidates"][0]["constraint_report"]["feasible"]
    assert float(cai_independent.score(designed, weights)["cai"]) == 0.9593415607455033
