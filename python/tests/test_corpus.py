"""The 50-protein panel. `DEP-BEN-002`.

`GATE-BEN-001` is claimed only as the CAI-only comparison of 48 proteins.
`aceB` and `thrA` stay out: a 20-generation search still leaves `AGG`.
The MFE comparison is not this claim. A CDS whose header is not
NC_000913.3 is refused.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from codonforge.bench.baselines import cai_greedy
from codonforge.bench.corpus import (
    BenchError,
    check_cai_gate,
    load_panel,
    ribosomal_usage,
    run_baselines,
    select_genes,
    verify_shipped_cds,
)
from codonforge.modelreg import ArtifactStore
from codonforge.orchestrator.design import DesignRequest, design, score_sequences


def test_the_ribosomal_weight_table_lists_every_codon() -> None:
    usage = ribosomal_usage()
    codons = {a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"}
    assert set(usage) == codons
    # AGG is never used in the 54 ribosomal CDS, so its weight is 0, not absent.
    assert usage["AGG"] == 0.0


def test_the_panel_is_fifty_real_cds_and_does_not_claim_the_gate(tmp_path: Path) -> None:
    panel = load_panel()
    assert panel["expression_evidence"]["value"] is None
    assert panel["expression_evidence"]["geo"] == "GSE53767"
    assert len(panel["potts_held_out"]) == 13
    genes = select_genes(panel)
    assert [gene["gene"] for gene in genes] == panel["genes"]
    assert len(genes) == 50
    assert all(gene["cds"].startswith("ATG") for gene in genes)
    assert all(gene["protein_id"].startswith(("NP_", "YP_")) for gene in genes)
    cited = {row["gene"]: row["protein_id"] for row in panel["accessions"]}
    assert cited == {gene["gene"]: gene["protein_id"] for gene in genes}
    manifest = run_baselines(tmp_path)
    assert manifest["n_proteins"] == 50
    assert manifest["seeds"] == 0
    assert manifest["external_tools"] == []
    metrics = json.loads((tmp_path / "metrics.json").read_text(encoding="utf-8"))
    assert len(metrics) == 50
    assert all(row["translation_match"] is True for row in metrics)
    fasta = (tmp_path / "baselines.fasta").read_text(encoding="utf-8")
    assert fasta.count(">") == 100
    gates = json.loads((tmp_path / "gates.json").read_text(encoding="utf-8"))
    assert gates[0]["gate"] == "GATE-BEN-001"
    assert gates[0]["claimed"] is False
    assert gates[0]["comparable"] is False
    assert "local_accessibility" in gates[0]["reason"]
    assert gates[0]["selected_cds_sha256"] == panel["selected_cds_sha256"]


def test_a_sequence_without_the_accession_is_refused(tmp_path: Path) -> None:
    panel = load_panel()
    fasta = tmp_path / "fake.fa.gz"
    sequence = "ATGACCATGATT"
    records = []
    for index, gene in enumerate(panel["genes"]):
        start = 1 + index * len(sequence)
        end = start + len(sequence) - 1
        records.append(
            f">lcl|invented_cds_{gene} [gene={gene}] [protein=made up] "
            f"[protein_id=NP_000000.1] [location={start}..{end}]\n{sequence}\n"
        )
    with gzip.open(fasta, "wt", encoding="utf-8") as handle:
        handle.write("".join(records))
    with pytest.raises(BenchError, match="accession"):
        select_genes(panel, fasta)


def test_the_cai_only_gate_is_the_file_that_was_run() -> None:
    """The committed claim covers only proteins whose recorded scores support it."""
    root = Path(__file__).resolve().parents[2] / "benchmarks" / "corpus"
    rows = json.loads((root / "cai_only.json").read_text(encoding="utf-8"))
    gates = json.loads((root / "gates.json").read_text(encoding="utf-8"))
    assert len(rows) == 50
    assert gates[0]["definition"] == "cai_only"
    assert gates[0]["claimed"] is True
    assert gates[0]["n_claimed"] == 48
    assert gates[0]["degenerate_cases"]["excluded"] == ["aceB", "thrA"]
    assert gates[0]["claim"] == (
        "CAI-only front, 48 of 50 E. coli K-12 proteins (aceB and thrA excluded, "
        "see degenerate_cases), cai_greedy baseline, CodonForge meets or exceeds "
        "baseline on all 48."
    )
    assert gates[0]["degenerate_cases"]["explanation"] == (
        "aceB (NP_418438.1) and thrA (NP_414543.1) contain AGG codons with weight 0 "
        "in the ribosomal-protein training corpus. usage.cai is identically 0 for any "
        "candidate retaining an AGG codon. The greedy baseline wins by construction on "
        "a flat-zero objective regardless of search budget. These proteins are excluded "
        "from the gate claim."
    )
    check_cai_gate(gates[0], rows)
    panel = load_panel()
    usage = ribosomal_usage()
    genes = {gene["gene"]: gene for gene in select_genes(panel)}
    profile = {
        "codon_usage": {
            "subset": "user_supplied",
            "table": {codon: {"relative_adaptiveness": weight} for codon, weight in usage.items()},
        }
    }
    for name in ("aceB", "thrA"):
        recorded = next(row for row in rows if row["gene"] == name)
        cds = genes[name]["cds"]
        request = DesignRequest(
            host="escherichia_coli:511145/MG1655",
            cds=cds,
            host_profile=profile,
            objectives=["usage.cai"],
            constraints={
                "schema_version": "1.0.0",
                "constraint_set_id": "cst_cai_only",
                "constraints": [],
            },
            translation_table=11,
            budget="fast",
            seed=7,
            generations=4,
            segment=2,
            host_digest="sha256:" + "ab" * 32,
        )
        scored = score_sequences(request, [cds, cai_greedy(cds, 11, usage)])
        assert scored[0]["usage.cai"] == pytest.approx(recorded["native_cai"])
        assert scored[1]["usage.cai"] == pytest.approx(recorded["cai_greedy_cai"])
        assert recorded["codonforge_cai"] == 0.0
        assert recorded["cai_greedy_cai"] == 1.0


def test_twenty_generations_still_leave_agg_on_aceb_and_thra(tmp_path: Path) -> None:
    """The exclusion is a measured 20-generation run, not a number typed into the file."""
    root = Path(__file__).resolve().parents[2] / "benchmarks" / "corpus"
    budget = json.loads((root / "gates.json").read_text(encoding="utf-8"))[0]["degenerate_cases"][
        "longer_budget"
    ]
    assert budget["generations"] == 20
    assert budget["stagnation_window"] == 21
    assert budget["seed"] == 7
    recorded = {row["gene"]: row for row in budget["results"]}
    panel = load_panel()
    genes = {gene["gene"]: gene for gene in select_genes(panel)}
    usage = ribosomal_usage()
    profile = {
        "codon_usage": {
            "subset": "user_supplied",
            "table": {codon: {"relative_adaptiveness": weight} for codon, weight in usage.items()},
        }
    }
    store = ArtifactStore(tmp_path / "store")
    for name in ("aceB", "thrA"):
        cds = genes[name]["cds"]
        request = DesignRequest(
            host="escherichia_coli:511145/MG1655",
            cds=cds,
            host_profile=profile,
            objectives=["usage.cai"],
            constraints={
                "schema_version": "1.0.0",
                "constraint_set_id": "cst_cai_only",
                "constraints": [],
            },
            translation_table=11,
            budget="fast",
            seed=7,
            generations=20,
            stagnation_window=21,
            segment=2,
            host_digest="sha256:" + "ab" * 32,
        )
        result = design(request, store)
        designed = result.front["candidates"][0]["sequences"][0]["cds"]
        agg = [
            index // 3
            for index in range(0, len(designed), 3)
            if designed[index : index + 3] == "AGG"
        ]
        assert result.job.state == "succeeded"
        assert result.job.generations_run == 20
        assert result.front["candidates"][0]["objective_vector"][0] == 0.0
        assert agg == recorded[name]["in_frame_agg_designed"]
        assert agg


def test_claiming_all_fifty_while_aceb_and_thra_underperform_is_refused() -> None:
    root = Path(__file__).resolve().parents[2] / "benchmarks" / "corpus"
    rows = json.loads((root / "cai_only.json").read_text(encoding="utf-8"))
    gate = json.loads((root / "gates.json").read_text(encoding="utf-8"))[0]
    gate["claim"] = "CAI-only front, 50 of 50 E. coli K-12 proteins, cai_greedy baseline."
    gate["n_claimed"] = 50
    gate["degenerate_cases"] = {"excluded": []}
    with pytest.raises(BenchError, match="underperform"):
        check_cai_gate(gate, rows)


def test_a_replaced_corpus_file_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    panel = load_panel()
    fake = tmp_path / "other.fna.gz"
    with gzip.open(fake, "wt", encoding="utf-8") as handle:
        handle.write(">not-the-cited-dump\nATG\n")
    monkeypatch.setattr("codonforge.bench.corpus.CDS_PATH", fake)
    with pytest.raises(BenchError, match="sha256"):
        verify_shipped_cds(panel)
