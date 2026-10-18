"""Harmonization, frequency matching, usage import, vendor screen, global fold.

The harmonization negative control is the one that has to be able to fail.
Scoring a synonym-changed CDS against the input CDS is above zero when the
source table and the host table are the same object. An implementation that
rebuilds the source vector from the candidate scores that edit as zero, and
this test is then red.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dna_design_suite import _core
from dna_design_suite.api.server import (
    _PRESETS,
    import_host_body,
    run_global_fold,
    run_vendor_screen,
)
from dna_design_suite.cli import main as cli_main
from dna_design_suite.codon.orchestrator import DesignError, DesignRequest, design
from dna_design_suite.codon.orchestrator.design import score_sequences
from dna_design_suite.codon.predictors.global_fold import global_fold_report, nussinov_pair_count
from dna_design_suite.ingest.usage_table import parse_usage_table
from dna_design_suite.modelreg import ArtifactStore
from dna_design_suite.vendor.screen import screen_sequence

ROOT = Path(__file__).resolve().parents[2]
HOST = "escherichia_coli:511145/MG1655"
CDS = "ATGCTGCGTAGCGGCGTGGCGCCGACCAAACTGCGTAGCTAA"
EMPTY = {
    "schema_version": "1.0.0",
    "constraint_set_id": "cst_none",
    "constraints": [],
}


def _profile() -> dict:
    path = ROOT / "examples" / "mg1655_ribosomal_usage.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _request(profile: dict, objectives: list[str], budget: str = "fast") -> DesignRequest:
    return DesignRequest(
        host=HOST,
        cds=CDS,
        host_profile=profile,
        objectives=objectives,
        constraints=EMPTY,
        translation_table=11,
        budget=budget,
        seed=7,
        generations=2,
        segment=2,
        stagnation_window=21,
        host_digest="sha256:" + "cd" * 32,
    )


def test_harmonize_native_is_the_input_and_a_synonym_scores_above_zero(tmp_path: Path) -> None:
    profile = _profile()
    profile["source_codon_usage"] = profile["codon_usage"]
    native_w = profile["codon_usage"]["table"]["CTG"]["relative_adaptiveness"]
    other_w = profile["codon_usage"]["table"]["TTA"]["relative_adaptiveness"]
    assert native_w != other_w, "CTG and TTA need different weights or the control is vacuous"
    profile["harmonize_native_cds"] = CDS
    request = _request(profile, ["usage.harmonize_rca"])
    native = score_sequences(request, [CDS])[0]["usage.harmonize_rca"]
    assert native == 0.0
    edited = "ATG" + "TTA" + CDS[6:]
    assert _core.translate(edited, 11) == _core.translate(CDS, 11)
    moved = score_sequences(request, [edited])[0]["usage.harmonize_rca"]
    assert moved > 0.0
    # The same edit scores 0 when the pinned native is the candidate. That is
    # the defect: rebuilding the source vector from the sequence being scored.
    broken = dict(profile)
    broken["harmonize_native_cds"] = edited
    self_score = score_sequences(_request(broken, ["usage.harmonize_rca"]), [edited])[0]
    assert self_score["usage.harmonize_rca"] == 0.0

    caller = _profile()
    caller["source_codon_usage"] = caller["codon_usage"]
    designed = design(_request(caller, ["usage.harmonize_rca"]), ArtifactStore(tmp_path / "store"))
    assert designed.job.state == "succeeded", designed.job.problem
    assert "harmonize_native_cds" not in caller
    score = designed.front["candidates"][0]["objective_vector"][0]
    assert score == 0.0
    assert designed.front["objective_sense"] == ["minimize"]


def test_harmonize_without_a_source_table_is_refused(tmp_path: Path) -> None:
    with pytest.raises(DesignError) as caught:
        design(_request(_profile(), ["usage.harmonize_rca"]), ArtifactStore(tmp_path / "store"))
    assert caught.value.code == "CF-HOST-003"


def test_cai_only_does_not_stamp_a_native_cds(tmp_path: Path) -> None:
    profile = _profile()
    result = design(_request(profile, ["usage.cai"]), ArtifactStore(tmp_path / "store"))
    assert result.job.state == "succeeded", result.job.problem
    assert "harmonize_native_cds" not in profile
    assert "source_codon_usage" not in profile


def test_match_frequency_matches_a_hand_sum_and_fast_warns(tmp_path: Path) -> None:
    profile = _profile()
    request = _request(profile, ["usage.match_frequency"], budget="fast")
    reported = score_sequences(request, [CDS])[0]["usage.match_frequency"]
    table = profile["codon_usage"]["table"]
    families: dict[str, list[str]] = {}
    for codon in table:
        try:
            amino = _core.translate(codon + "AAA", 11)[:1]
        except ValueError:
            continue
        if amino in {"*", "?"}:
            continue
        families.setdefault(amino, []).append(codon)
    counts: dict[str, int] = {}
    for index in range(0, len(CDS), 3):
        codon = CDS[index : index + 3]
        counts[codon] = counts.get(codon, 0) + 1
    hand = 0.0
    for members in families.values():
        if len(members) <= 1:
            continue
        weights = [float(table[codon]["relative_adaptiveness"]) for codon in members]
        total_w = sum(weights)
        if total_w == 0:
            continue
        n = sum(counts.get(codon, 0) for codon in members)
        if n == 0:
            continue
        for codon, weight in zip(members, weights, strict=True):
            observed = counts.get(codon, 0) / n
            diff = observed - weight / total_w
            hand += diff * diff
    assert reported == pytest.approx(hand, rel=1e-9, abs=1e-12)
    assert hand > 0.0

    result = design(request, ArtifactStore(tmp_path / "store"))
    assert result.job.state == "succeeded", result.job.problem
    warnings = [item for item in result.warnings if item["code"] == "CF-SRCH-010"]
    assert warnings
    assert "L0" not in warnings[0]["detail"]
    assert isinstance(result.front["candidates"][0]["objective_vector"][0], float)
    radii = json.loads(
        _core.term_radii(
            json.dumps(
                {
                    "cds": CDS,
                    "translation_table": 11,
                    "five_prime_flank": "",
                    "three_prime_flank": "",
                    "constraint_set": EMPTY,
                    "terms": ["usage.match_frequency"],
                    "weights": [1.0],
                    "host_profile": profile,
                    "host_id": HOST,
                    "host_digest": "sha256:" + "cd" * 32,
                    "root_seed": 7,
                }
            )
        )
    )
    assert radii["usage.match_frequency"] == "unbounded"


def test_kazusa_csv_import_keeps_user_supplied_and_normalizes() -> None:
    text = (ROOT / "benchmarks" / "validation" / "kazusa" / "e_coli_316407.csv").read_text(
        encoding="utf-8"
    )
    profile = parse_usage_table(text, organism="Kazusa taxid 316407, file in this repo")
    usage = profile["codon_usage"]
    assert usage["subset"] == "user_supplied"
    alanine = [
        usage["table"][codon]["relative_adaptiveness"] for codon in ("GCA", "GCT", "GCC", "GCG")
    ]
    assert max(alanine) == 1.0
    assert usage["table"]["GCG"]["relative_adaptiveness"] == 1.0
    status, body = import_host_body({"text": text, "subset": "user_supplied"})
    assert status == 200
    assert body["codon_usage"]["table"]["GCG"]["relative_adaptiveness"] == 1.0


def test_vendor_screen_flags_a_homopolymer_and_is_not_a_vendor_score() -> None:
    report = screen_sequence("AAAAAAAATGAAATAA")
    assert report["vendor_score"] is False
    assert report["wet_lab_expression"] == "out_of_scope"
    assert report["hairpin"]["scored"] is False
    flagged = [row for row in report["checks"] if row["kind"] == "homopolymer"]
    assert flagged and all(row["passed"] is False for row in flagged)
    assert all(row["vendor_cutoff_transcribed"] is False for row in report["checks"])
    status, body = run_vendor_screen({"cds": "AAAAAAAATGAAATAA"})
    assert status == 200
    assert body["n_flags"] >= 1


def test_global_fold_labels_the_engine(monkeypatch: pytest.MonkeyPatch) -> None:
    assert nussinov_pair_count("AAAA") == 0
    assert nussinov_pair_count("GAAAAAC") == 1
    live = global_fold_report("GAAAAAC")
    assert live["not_windowed_accessibility"] is True
    assert live["wet_lab_expression"] == "out_of_scope"
    if live["engine"] == "viennarna":
        assert live["unit"] == "kcal/mol"
        assert live["quantity"] == "global_mfe_kcal_per_mol"
        assert "not windowed" in live["label"]
    else:
        assert live["unit"] == "base_pairs"
        assert live["unit"] != "kcal/mol"

    monkeypatch.setattr(
        "dna_design_suite.codon.predictors.global_fold._vienna_mfe",
        lambda *_args, **_kwargs: None,
    )
    fallback = global_fold_report("GAAAAAC")
    assert fallback["engine"] == "nussinov_fallback"
    assert fallback["unit"] == "base_pairs"
    assert fallback["value"] == 1
    assert "kcal/mol" not in fallback["label"] or "not a kcal/mol" in fallback["label"]
    long = "A" * 401
    refused = global_fold_report(long)
    assert refused["engine"] == "unavailable"
    assert refused["value"] is None
    status, body = run_global_fold({"cds": "GAAAAAC"})
    assert status == 200
    assert body["kind"] == "global_fold"


def test_presets_and_cli_import(capsys: pytest.CaptureFixture[str]) -> None:
    assert _PRESETS["harmonize"] == ["usage.harmonize_rca"]
    assert _PRESETS["match_frequency"] == ["usage.match_frequency"]
    code = cli_main(
        [
            "hosts",
            "import",
            "--input",
            str(ROOT / "benchmarks" / "validation" / "kazusa" / "e_coli_316407.csv"),
        ]
    )
    assert code == 0
    printed = json.loads(capsys.readouterr().out)
    assert printed["codon_usage"]["subset"] == "user_supplied"
    page = (ROOT / "python" / "dna_design_suite" / "api" / "static" / "app.js").read_text(
        encoding="utf-8"
    )
    html = (ROOT / "python" / "dna_design_suite" / "api" / "static" / "index.html").read_text(
        encoding="utf-8"
    )
    needles = (
        "/v1/hosts:import",
        "/v1/vendor-screen",
        "/v1/structure:global",
        "readUsageFile",
    )
    for needle in needles:
        assert needle in page
    assert "usage.harmonize_rca" in html
    assert "usage.match_frequency" in html
