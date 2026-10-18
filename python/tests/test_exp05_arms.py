"""The two EXP-05 arms must not be the same run.

Baldwinian scores the repaired sequence and breeds from the unrepaired
genome, then repairs that genome again on the next generation. The repair
counts are the evidence. Identical counts mean the flag was not read.
"""

from __future__ import annotations

from pathlib import Path

from dna_design_suite.codon.bench.corpus import load_panel, select_genes
from dna_design_suite.codon.bench.exp05_campaign import host_profile, run_one
from dna_design_suite.modelreg import ArtifactStore


def test_baldwinian_repair_counts_differ_from_lamarckian(tmp_path: Path) -> None:
    panel = load_panel()
    gene = {row["gene"]: row for row in select_genes(panel)}["rpsJ"]
    profile = host_profile()
    store = ArtifactStore(tmp_path / "store")
    lamarckian = 0
    baldwinian = 0
    for seed in range(4):
        lamarckian += run_one(gene, profile, seed, True, store).repair_invocations
        baldwinian += run_one(gene, profile, seed, False, store).repair_invocations
    assert baldwinian != lamarckian, (
        f"repair counts match ({lamarckian}); Baldwinian did not reach the search"
    )
