"""Two cistrons on one flux-weighted demand ledger.

The joint growth cost is `burden.growth_cost` of the summed counts. It is not
the sum of the two solo scores. A pass that ignores the other cistron leaves
the piled pair in place, and the true joint cost of that pair is worse.
"""

from __future__ import annotations

import json

from test_orchestrator import HOST

from codonforge import _core
from codonforge.modelreg import ArtifactStore
from codonforge.orchestrator.design import DesignRequest
from codonforge.orchestrator.operon import (
    Cistron,
    _growth,
    _joint_cost,
    _weighted_counts,
    codon_counts,
    design_operon,
)

# M L *, both starting on the host's preferred leucine codon.
CDS = "ATGCTGTAA"
TABLE = 11


def leucine_profile() -> dict:
    weights = {
        "CTG": 1.0,
        "TTA": 0.8,
        "TTG": 0.6,
        "CTT": 0.4,
        "CTC": 0.2,
        "CTA": 0.1,
        "ATG": 1.0,
    }
    return {
        "codon_usage": {
            "subset": "highly_expressed",
            "table": {
                codon: {"relative_adaptiveness": weight} for codon, weight in weights.items()
            },
        }
    }


def test_integer_counts_match_the_sequence_term() -> None:
    profile = leucine_profile()
    request = DesignRequest(
        host=HOST,
        cds=CDS,
        host_profile=profile,
        objectives=["burden.growth_cost"],
        translation_table=TABLE,
        host_digest="sha256:" + "1" * 64,
    )
    scored = json.loads(_core.term_scores(_directions(request), CDS))
    direct = _growth(profile, TABLE, codon_counts(CDS, TABLE))
    assert direct == scored["burden.growth_cost"]


def test_joint_cost_is_the_shared_ledger(tmp_path) -> None:
    profile = leucine_profile()
    result = design_operon(
        [Cistron(CDS, 1.0, "a"), Cistron(CDS, 1.0, "b")],
        ArtifactStore(tmp_path / "store"),
        host=HOST,
        host_profile=profile,
        translation_table=TABLE,
        budget="fast",
        seed=7,
        generations=0,
        passes=2,
        host_digest="sha256:" + "1" * 64,
    )
    assert result.proteins == ["ML*", "ML*"]
    direct = _growth(profile, TABLE, _weighted_counts(result.sequences, result.fluxes, TABLE))
    assert result.joint_growth_cost == direct
    assert result.joint_growth_cost != sum(result.solo_growth_costs)
    assert result.joint_growth_cost < result.independent_joint_growth_cost
    assert result.swaps_per_pass[0] >= 1
    assert result.swaps_per_pass[1] == 0
    piled = _joint_cost([CDS, CDS], [1.0, 1.0], TABLE, profile)
    assert result.joint_growth_cost < piled


def test_ignoring_the_other_cistron_leaves_the_worse_joint(tmp_path, monkeypatch) -> None:
    """A scorer that drops background counts cannot find the joint improvement."""
    profile = leucine_profile()
    real = _weighted_counts

    def only_the_first(sequences, fluxes, translation_table):
        return real(sequences[:1], fluxes[:1], translation_table)

    monkeypatch.setattr("codonforge.orchestrator.operon._weighted_counts", only_the_first)
    broken = design_operon(
        [Cistron(CDS, 1.0, "a"), Cistron(CDS, 1.0, "b")],
        ArtifactStore(tmp_path / "broken"),
        host=HOST,
        host_profile=profile,
        translation_table=TABLE,
        budget="fast",
        seed=7,
        generations=0,
        passes=2,
        host_digest="sha256:" + "1" * 64,
    )
    monkeypatch.undo()
    assert broken.sequences == broken.independent_sequences
    ignored = _joint_cost(broken.sequences, [1.0, 1.0], TABLE, profile)
    honest = design_operon(
        [Cistron(CDS, 1.0, "a"), Cistron(CDS, 1.0, "b")],
        ArtifactStore(tmp_path / "honest"),
        host=HOST,
        host_profile=profile,
        translation_table=TABLE,
        budget="fast",
        seed=7,
        generations=0,
        passes=2,
        host_digest="sha256:" + "1" * 64,
    )
    assert honest.joint_growth_cost < ignored


def _directions(request: DesignRequest) -> str:
    return json.dumps(
        {
            "cds": request.cds,
            "translation_table": request.translation_table,
            "five_prime_flank": "",
            "three_prime_flank": "",
            "constraint_set": {
                "schema_version": "1.0.0",
                "constraint_set_id": "cst_none",
                "constraints": [],
            },
            "host_profile": request.host_profile,
            "host_id": request.host,
            "host_digest": request.host_digest,
            "context_digest": None,
            "terms": list(request.objectives),
            "weights": [1.0],
            "root_seed": request.seed,
        }
    )
