"""Gibson context: flanks participate in the 12-mer scan, arms are copies.

The context case is the one that has to fail if the flank is dropped. A CDS
that contains ``ACGTACGTACGT`` once is legal on its own. The same CDS in front
of a promoter that also contains that 12-mer is a repeat the codon search can
edit. A repeat that never touches the CDS is reported and is not actionable.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dna_design_suite import _core
from dna_design_suite.api import serve
from dna_design_suite.api.server import run_gibson
from dna_design_suite.cli import main as cli_main
from dna_design_suite.codon.gibson import (
    ARM_LENGTH_DEFAULT,
    GibsonError,
    assess,
    constraint_set,
    format_text,
)
from dna_design_suite.codon.orchestrator import DesignRequest, design
from dna_design_suite.modelreg import ArtifactStore

ROOT = Path(__file__).resolve().parents[2]
SHARED = "ACGTACGTACGT"
CDS = "ATG" + SHARED + "GCTTAA"
PROMOTER = "GGGGGGGG" + SHARED
HOST = "escherichia_coli:511145/BL21"


def _profile() -> dict:
    from test_orchestrator import host_profile

    return host_profile()


def test_a_promoter_repeat_is_visible_only_when_the_flank_is_supplied() -> None:
    alone = assess(CDS, "", "", arm_length=20)
    assert all(site["kmer"] != SHARED for site in alone["recombination"]["sites"])
    report = assess(CDS, PROMOTER, "", arm_length=20)
    matches = [site for site in report["recombination"]["sites"] if site["kmer"] == SHARED]
    assert len(matches) == 1
    site = matches[0]
    assert site["class"] == "cds_context"
    assert site["actionable"] is True
    assert site["count"] == 2
    assert site["in_five_prime_arm"] is True
    places = {item["where"] for item in site["occurrences"]}
    assert places == {"five_prime_flank", "cds"}
    assert report["guarantees_assembly"] is False
    assert report["homology_arms"]["tm_celsius"] is None
    assert report["repeat_rule"]["k"] == 12
    assert report["repeat_rule"]["canonical_strand"] is False
    assert report["repeat_rule"]["source"] == "examples/gibson.json"


def test_a_flank_only_repeat_is_not_actionable() -> None:
    flank = SHARED + "AAAA" + SHARED
    coding = "ATGAAAGGGCCCTTTTAA"
    assert SHARED not in coding
    report = assess(coding, flank, "", arm_length=20)
    matches = [site for site in report["recombination"]["sites"] if site["kmer"] == SHARED]
    assert len(matches) == 1
    assert matches[0]["class"] == "flank_only"
    assert matches[0]["actionable"] is False
    assert report["recombination"]["n_flank_only"] >= 1


def test_two_copies_inside_the_cds_are_an_internal_repeat() -> None:
    coding = "ATG" + SHARED + "GGG" + SHARED + "TAA"
    report = assess(coding, "TTTTTTTTTTTTTTTTTTTT", "CCCCCCCCCCCCCCCCCCCC", arm_length=20)
    matches = [site for site in report["recombination"]["sites"] if site["kmer"] == SHARED]
    assert len(matches) == 1
    assert matches[0]["class"] == "cds_internal"
    assert matches[0]["actionable"] is True
    assert all(item["where"] == "cds" for item in matches[0]["occurrences"])


def test_an_arm_is_a_copy_and_a_short_flank_is_not_padded() -> None:
    five = "TTGACAATTAATCATCGGCTCGTATAATGTGTGGA"
    three = "AAAGCCCGCTCTAGAGGATCCCCG"
    report = assess("ATGTAA", five, three, arm_length=20)
    assert report["homology_arms"]["five_prime"]["sequence"] == five[-20:]
    assert report["homology_arms"]["three_prime"]["sequence"] == three[:20]
    assert report["homology_arms"]["fragment"] == five[-20:] + "ATGTAA" + three[:20]
    short = assess("ATGTAA", "ACGTACGTACGTACGT", three, arm_length=20)
    assert short["homology_arms"]["five_prime"]["designed"] is False
    assert short["homology_arms"]["five_prime"]["sequence"] is None
    assert short["homology_arms"]["five_prime"]["code"] == "CF-GIB-001"
    assert short["homology_arms"]["fragment"] is None
    assert "invented" in short["homology_arms"]["fragment_note"]
    assert any(item["code"] == "CF-GIB-001" for item in short["problems"])
    assert short["homology_arms"]["three_prime"]["designed"] is True


def test_malformed_input_is_cf_gib_003() -> None:
    with pytest.raises(GibsonError) as empty:
        assess("", "ACGT" * 10, "ACGT" * 10)
    assert empty.value.code == "CF-GIB-003"
    with pytest.raises(GibsonError) as base:
        assess("ATGNAA", "ACGT" * 10, "")
    assert "N" in base.value.detail
    with pytest.raises(GibsonError) as uracil:
        assess("ATGUTAA", "ACGT" * 10, "")
    assert uracil.value.code == "CF-GIB-003"
    with pytest.raises(GibsonError) as low:
        assess("ATGTAA", "ACGT" * 10, "", arm_length=19)
    assert low.value.code == "CF-GIB-003"
    with pytest.raises(GibsonError) as high:
        assess("ATGTAA", "ACGT" * 10, "", arm_length=151)
    assert high.value.code == "CF-GIB-003"
    missing = assess("ATGTAA", "", "", arm_length=ARM_LENGTH_DEFAULT)
    assert {item["code"] for item in missing["problems"]} == {"CF-GIB-002"}
    assert missing["homology_arms"]["fragment"] is None


def test_the_report_states_what_it_does_not_scan() -> None:
    report = assess("ATGTAA", "ACGT" * 10, "TGCA" * 10, arm_length=20)
    text = " ".join(report["limitations"])
    assert "circular" in text
    assert "reverse complement" in text
    assert "not supplied" in text
    assert "Melting temperature is not computed" in report["homology_arms"]["method"]
    assert "every gene" in text
    file_rule = next(
        item["rule"]
        for item in constraint_set()["constraints"]
        if item["rule"]["kind"] == "kmer_repeat"
    )
    assert report["repeat_rule"]["k"] == file_rule["k"]
    assert report["repeat_rule"]["max_occurrences"] == file_rule["max_occurrences"]


def test_registered_titles_match_the_rust_registry() -> None:
    registry = (ROOT / "crates" / "cf-types" / "src" / "error.rs").read_text(encoding="utf-8")
    from dna_design_suite.codon.gibson import TITLES

    for code, title in TITLES.items():
        assert f'"{code}"' in registry
        assert title in registry
    documented = (ROOT / "docs" / "gibson-assembly.md").read_text(encoding="utf-8")
    for code in ("CF-GIB-001", "CF-GIB-002", "CF-GIB-003"):
        assert code in documented
        assert f'"{code}"' in registry


def test_cli_feasibility_sees_the_promoter_and_misses_it_without_the_flank(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    path = tmp_path / "gibson.json"
    path.write_text(json.dumps(constraint_set()), encoding="utf-8")
    bare = cli_main(["feasibility", "--cds", CDS, "--constraints", str(path), "--output", "json"])
    assert bare == 0
    capsys.readouterr()
    flanked = cli_main(
        [
            "feasibility",
            "--cds",
            CDS,
            "--constraints",
            str(path),
            "--five-prime-flank",
            PROMOTER,
            "--output",
            "json",
        ]
    )
    assert flanked == 3
    payload = json.loads(capsys.readouterr().out)
    assert any(SHARED in item["detail"] and item["actionable"] for item in payload["violations"])


def test_cli_gibson_exits_3_when_the_arm_cannot_be_copied(
    capsys: pytest.CaptureFixture[str],
) -> None:
    code = cli_main(
        [
            "gibson",
            "--cds",
            "ATGTAA",
            "--five-prime-flank",
            "ACGT",
            "--three-prime-flank",
            "ACGT" * 10,
            "--arm-length",
            "20",
        ]
    )
    assert code == 3
    report = json.loads(capsys.readouterr().out)
    assert report["homology_arms"]["five_prime"]["sequence"] is None
    text = format_text(report)
    assert "CF-GIB-001" in text
    assert "not reported" in text


def test_post_v1_gibson_returns_the_report(tmp_path: Path) -> None:
    status, body = run_gibson(
        {
            "cds": "ATGTAA",
            "five_prime_flank": "TTGACAATTAATCATCGGCTCGTATAATGTGTGGA",
            "three_prime_flank": "AAAGCCCGCTCTAGAGGATCCCCG",
            "arm_length": 20,
        }
    )
    assert status == 200
    assert body["kind"] == "gibson_assembly"
    assert body["guarantees_assembly"] is False
    refused, problem = run_gibson({"cds": "ATGUAA", "arm_length": 40})
    assert refused == 422
    assert problem["code"] == "CF-GIB-003"
    server = serve(ArtifactStore(tmp_path / "artifacts"))
    url = f"http://127.0.0.1:{server.server_address[1]}/v1/gibson"
    import urllib.request

    posted = {"cds": "ATGTAA", "five_prime_flank": "ACGT" * 10, "arm_length": 20}
    request = urllib.request.Request(
        url,
        data=json.dumps(posted).encode(),
        headers={"Authorization": "Bearer t_localdev1", "Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request) as response:
        payload = json.loads(response.read())
    server.shutdown()
    assert payload["homology_arms"]["five_prime"]["designed"] is True


def test_the_page_offers_the_preset_and_posts_the_report() -> None:
    html = (ROOT / "python" / "dna_design_suite" / "api" / "static" / "index.html").read_text(
        encoding="utf-8"
    )
    script = (ROOT / "python" / "dna_design_suite" / "api" / "static" / "app.js").read_text(
        encoding="utf-8"
    )
    assert 'data-preset="gibson.json"' in html
    assert 'id="gibson-report"' in html
    assert "/v1/gibson" in script
    assert "renderGibsonReport" in script
    assert "guarantee_statement" in script


def test_design_removes_the_promoter_12mer_and_keeps_the_protein(tmp_path: Path) -> None:
    """The input is infeasible. A search that ignores the flank still returns it."""
    protein = _core.translate(CDS, 11)
    result = design(
        DesignRequest(
            host=HOST,
            cds=CDS,
            host_profile=_profile(),
            objectives=["usage.cai"],
            constraints=constraint_set(),
            five_prime_flank=PROMOTER,
            translation_table=11,
            budget="fast",
            seed=7,
            generations=1,
            host_digest="sha256:" + "1" * 64,
        ),
        ArtifactStore(tmp_path / "artifacts"),
    )
    assert result.job.state == "succeeded", result.status.get("problem")
    best = result.best()["sequences"][0]["cds"]
    assert _core.translate(best, 11) == protein
    assert SHARED not in best
    report = assess(best, PROMOTER, "", arm_length=20)
    assert report["recombination"]["n_actionable"] == 0
    assert result.best()["constraint_report"]["feasible"] is True
