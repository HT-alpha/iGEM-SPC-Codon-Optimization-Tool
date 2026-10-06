"""The 50-protein panel. `DEP-BEN-002`.

`GATE-BEN-001` is a non-claim: LinearFold's global MFE is not
`structure.local_accessibility`. A CDS whose header is not NC_000913.3
is refused.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from codonforge.bench.corpus import (
    BenchError,
    load_panel,
    run_baselines,
    select_genes,
    verify_shipped_cds,
)


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


def test_a_replaced_corpus_file_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    panel = load_panel()
    fake = tmp_path / "other.fna.gz"
    with gzip.open(fake, "wt", encoding="utf-8") as handle:
        handle.write(">not-the-cited-dump\nATG\n")
    monkeypatch.setattr("codonforge.bench.corpus.CDS_PATH", fake)
    with pytest.raises(BenchError, match="sha256"):
        verify_shipped_cds(panel)
