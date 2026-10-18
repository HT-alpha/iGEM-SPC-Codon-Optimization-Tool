"""A run of zero-weight codons must not flatten the CAI search.

`usage.cai` is 0 while any counted codon has weight 0. The arithmetic delta
of removing one of them is therefore also 0, and so is the delta of every
other codon. Tiers 1 and 2 sum those deltas. `DEP-SRCH-020` linearizes CAI as
a sum of log-weights in that case. The reported score stays the geometric
mean.

Seed 11 on twelve `AGG` codons is the case that showed the flat sum: the
search removed the zeros and still stopped at CAI 0.399, below `cai_greedy`.
"""

from __future__ import annotations

import json
from pathlib import Path

from dna_design_suite import _core
from dna_design_suite.codon.bench.baselines import cai_greedy
from dna_design_suite.codon.bench.corpus import ribosomal_usage
from dna_design_suite.codon.orchestrator.design import DesignRequest, design, score_sequences
from dna_design_suite.modelreg import ArtifactStore

ROOT = Path(__file__).resolve().parents[2]


def _profile() -> tuple[dict[str, float], dict]:
    weights = ribosomal_usage()
    profile = {
        "codon_usage": {
            "subset": "user_supplied",
            "table": {
                codon: {"relative_adaptiveness": weight} for codon, weight in weights.items()
            },
        }
    }
    return weights, profile


def test_twelve_zero_weight_arginines_match_cai_greedy(tmp_path: Path) -> None:
    weights, profile = _profile()
    assert weights["AGG"] == 0.0
    cds = "ATG" + "AGG" * 12 + "CTG" * 4 + "AAATAA"
    greedy = cai_greedy(cds, 11, weights)
    protein = _core.translate(cds, 11)
    store = ArtifactStore(tmp_path / "store")
    for seed in (7, 11):
        request = DesignRequest(
            host="escherichia_coli:511145/MG1655",
            cds=cds,
            host_profile=profile,
            objectives=["usage.cai"],
            constraints={
                "schema_version": "1.0.0",
                "constraint_set_id": "cst_zero",
                "constraints": [],
            },
            translation_table=11,
            budget="fast",
            seed=seed,
            generations=4,
            segment=2,
            stagnation_window=21,
            host_digest="sha256:" + "ab" * 32,
        )
        native = score_sequences(request, [cds])[0]["usage.cai"]
        assert native == 0.0
        result = design(request, store)
        assert result.job.state == "succeeded", result.job.problem
        designed = result.front["candidates"][0]["sequences"][0]["cds"]
        cai = result.front["candidates"][0]["objective_vector"][0]
        greedy_cai = score_sequences(request, [greedy])[0]["usage.cai"]
        codons = [designed[index : index + 3] for index in range(0, len(designed), 3)]
        assert "AGG" not in codons
        assert _core.translate(designed, 11) == protein
        assert cai == greedy_cai
        assert cai == 1.0


def test_the_ribosomal_table_on_disk_still_lists_agg_at_zero() -> None:
    document = json.loads((ROOT / "examples" / "mg1655_ribosomal_usage.json").read_text())
    assert document["codon_usage"]["table"]["AGG"]["relative_adaptiveness"] == 0.0
