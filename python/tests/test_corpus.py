"""The 50-protein panel. `DEP-BEN-002`.

`GATE-BEN-001` is claimed as the CAI-only comparison of all 50 proteins.
`aceB` and `thrA` meet `cai_greedy` at the gate budget (`DEP-SRCH-020`).
The MFE comparison is not this claim. A CDS whose header is not
NC_000913.3 is refused.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from dna_design_suite import _core
from dna_design_suite.codon.bench.baselines import cai_greedy
from dna_design_suite.codon.bench.corpus import (
    BenchError,
    cai_gate_document,
    check_cai_gate,
    load_panel,
    ribosomal_usage,
    run_baselines,
    select_genes,
    verify_shipped_cds,
)
from dna_design_suite.codon.orchestrator.design import DesignRequest, design, score_sequences
from dna_design_suite.modelreg import ArtifactStore


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
    """The committed claim is the document the rows support, and only that."""
    root = Path(__file__).resolve().parents[2] / "benchmarks" / "corpus"
    rows = json.loads((root / "cai_only.json").read_text(encoding="utf-8"))
    gates = json.loads((root / "gates.json").read_text(encoding="utf-8"))
    panel = load_panel()
    assert len(rows) == 50
    assert gates == [cai_gate_document(rows, selected_cds_sha256=panel["selected_cds_sha256"])]
    assert gates[0]["claimed"] is True
    assert gates[0]["n_claimed"] == 50
    assert gates[0]["degenerate_cases"]["excluded"] == []
    check_cai_gate(gates[0], rows)
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
        assert recorded["codonforge_cai"] == 1.0
        assert recorded["cai_greedy_cai"] == 1.0


def test_aceb_and_thra_meet_cai_greedy_at_the_gate_budget(tmp_path: Path) -> None:
    """The flat arithmetic delta used to leave in-frame AGG and CAI 0.

    The gate budget is fast, 4 generations, seed 7. A 20-generation run of the
    old linearization still scored 0. This is the same two proteins, and the
    assertion is the one the gate uses: DNA Design Suite at least `cai_greedy`, and
    the protein unchanged.
    """
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
        greedy = cai_greedy(cds, 11, usage)
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
            stagnation_window=21,
            segment=2,
            host_digest="sha256:" + "ab" * 32,
        )
        result = design(request, store)
        designed = result.front["candidates"][0]["sequences"][0]["cds"]
        scored = score_sequences(request, [designed, greedy])
        agg = [
            index // 3
            for index in range(0, len(designed), 3)
            if designed[index : index + 3] == "AGG"
        ]
        assert result.job.state == "succeeded"
        assert agg == []
        assert _core.translate(designed, 11) == _core.translate(cds, 11)
        assert scored[0]["usage.cai"] + 1e-12 >= scored[1]["usage.cai"]
        assert result.front["candidates"][0]["objective_vector"][0] == scored[0]["usage.cai"]


def test_a_below_greedy_score_on_a_claimed_protein_is_refused() -> None:
    """Injecting the old aceB score of 0 must make the 50-protein claim fail."""
    root = Path(__file__).resolve().parents[2] / "benchmarks" / "corpus"
    rows = json.loads((root / "cai_only.json").read_text(encoding="utf-8"))
    gate = json.loads((root / "gates.json").read_text(encoding="utf-8"))[0]
    for row in rows:
        if row["gene"] == "aceB":
            row["codonforge_cai"] = 0.0
    with pytest.raises(BenchError, match="underperform"):
        check_cai_gate(gate, rows)


def test_excluding_a_protein_that_meets_cai_greedy_is_refused() -> None:
    root = Path(__file__).resolve().parents[2] / "benchmarks" / "corpus"
    rows = json.loads((root / "cai_only.json").read_text(encoding="utf-8"))
    gate = json.loads((root / "gates.json").read_text(encoding="utf-8"))[0]
    gate["n_claimed"] = 49
    gate["claim"] = (
        "CAI-only front, 49 of 50 E. coli K-12 proteins (aceB excluded, "
        "see degenerate_cases), cai_greedy baseline, DNA Design Suite meets or exceeds "
        "baseline on all 49."
    )
    gate["degenerate_cases"] = {"excluded": ["aceB"]}
    with pytest.raises(BenchError, match="excluded proteins that meet"):
        check_cai_gate(gate, rows)


def test_a_replaced_corpus_file_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    panel = load_panel()
    fake = tmp_path / "other.fna.gz"
    with gzip.open(fake, "wt", encoding="utf-8") as handle:
        handle.write(">not-the-cited-dump\nATG\n")
    monkeypatch.setattr("dna_design_suite.codon.bench.corpus.CDS_PATH", fake)
    with pytest.raises(BenchError, match="sha256"):
        verify_shipped_cds(panel)
