"""The 13-protein panel. `DEP-BEN-002`.

`GATE-BEN-001` is a non-claim: fewer than 50 proteins, no competitor arm.
A CDS whose header is not NC_000913.3 is refused.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from codonforge.bench.corpus import BenchError, load_panel, run_baselines, select_genes


def test_the_panel_is_thirteen_real_cds_and_does_not_claim_the_gate(tmp_path: Path) -> None:
    panel = load_panel()
    assert panel["expression_evidence"]["value"] is None
    assert panel["expression_evidence"]["geo"] == "GSE53767"
    genes = select_genes(panel)
    assert [gene["gene"] for gene in genes] == panel["genes"]
    assert len(genes) == 13
    assert all(gene["cds"].startswith("ATG") for gene in genes)
    assert all(gene["protein_id"].startswith(("NP_", "YP_")) for gene in genes)
    manifest = run_baselines(tmp_path)
    assert manifest["n_proteins"] == 13
    assert manifest["seeds"] == 0
    assert manifest["external_tools"] == []
    metrics = json.loads((tmp_path / "metrics.json").read_text(encoding="utf-8"))
    assert len(metrics) == 13
    assert all(row["translation_match"] is True for row in metrics)
    gates = json.loads((tmp_path / "gates.json").read_text(encoding="utf-8"))
    assert gates[0]["gate"] == "GATE-BEN-001"
    assert gates[0]["claimed"] is False


def test_a_sequence_without_the_accession_is_refused(tmp_path: Path) -> None:
    fasta = tmp_path / "fake.fa.gz"
    record = ">lcl|invented_cds_none [gene=lacZ] [protein=made up] [location=1..12]\nATGACCATGATT\n"
    with gzip.open(fasta, "wt", encoding="utf-8") as handle:
        handle.write(record)
    with pytest.raises(BenchError, match="accession"):
        select_genes(load_panel(), fasta)
