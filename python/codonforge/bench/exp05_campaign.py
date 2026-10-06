"""EXP-05 on three panel proteins.

Both arms use the same CDS, the same host profile, the same generation
budget and the same fixed hypervolume reference. The repair flag is the
only intentional difference. `decide` is applied as written.

The native and CAI-greedy baselines do not run the population loop, so
they have no generation budget to copy. The budget below is the one both
arms share.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Any

from codonforge.bench.corpus import load_panel, ribosomal_usage, select_genes
from codonforge.bench.exp05 import Trial, decide
from codonforge.bench.reference import hypervolume
from codonforge.modelreg import ArtifactStore
from codonforge.orchestrator.design import DesignRequest, design

# Three shortest panel CDS. Twenty seeds is 120 jobs at this budget.
PROTEINS = ("rpsJ", "rplL", "rpsL")
SEEDS = 20
GENERATIONS = 20
SEGMENT = 4
# The engine stops on hypervolume stagnation after 8 generations. The
# 4-generation campaign never reached that window. This one sets the
# window above the budget so both arms actually run 20 generations.
STAGNATION_WINDOW = 21
BUDGET = "standard"
# Low end of the defined interval of a geometric mean of weights in (0, 1].
CAI_REFERENCE = 0.0
REFERENCE = (CAI_REFERENCE,)
SENSES = ["maximize"]


def host_profile() -> dict[str, Any]:
    """Ribosomal-protein relative adaptiveness, labelled user_supplied.

    The counts are the corpus rule (gene names rpl/rps/rpm). The subset
    name is user_supplied so this is not presented as the engine's
    highly_expressed table.
    """
    usage = ribosomal_usage()
    table = {codon: {"relative_adaptiveness": weight} for codon, weight in usage.items()}
    return {"codon_usage": {"subset": "user_supplied", "table": table}}


def constraints_for(cds: str, lamarckian: bool) -> dict[str, Any]:
    """A motif that occurs in this CDS, so repair has something to do."""
    motif = cds[30:34]
    if motif not in cds:
        raise ValueError("interior motif is not in the CDS")
    return {
        "schema_version": "1.0.0",
        "constraint_set_id": "cst_exp05",
        "constraints": [
            {
                "id": "interior_motif",
                "class": "blocking",
                "scope": "construct",
                "rule": {
                    "kind": "forbidden_motif",
                    "motif": motif,
                    "strand": "both",
                    "max_occurrences": 0,
                },
            }
        ],
        "repair_policy": {"lamarckian": lamarckian},
    }


def _diversity(members: list[dict[str, Any]]) -> float:
    sequences = []
    for member in members:
        document = member.get("inherited") or member["sequence"]
        sequences.append(str(document["cds"]))
    if len(sequences) < 2:
        return 0.0
    total = 0.0
    pairs = 0
    width = len(sequences[0]) // 3
    for left in range(len(sequences)):
        for right in range(left + 1, len(sequences)):
            mismatches = sum(
                sequences[left][index : index + 3] != sequences[right][index : index + 3]
                for index in range(0, len(sequences[left]), 3)
            )
            total += mismatches / width
            pairs += 1
    return total / pairs


def _hypervolume(front: dict[str, Any]) -> float:
    points = []
    for candidate in front["candidates"]:
        vector = candidate.get("objective_vector") or []
        if not vector or vector[0] is None:
            continue
        if candidate.get("constraint_report", {}).get("feasible") is False:
            continue
        points.append([float(vector[0])])
    if not points:
        return 0.0
    return hypervolume(points, list(REFERENCE), SENSES)


def run_one(
    gene: dict[str, Any],
    profile: dict[str, Any],
    seed: int,
    lamarckian: bool,
    store: ArtifactStore,
) -> Trial:
    request = DesignRequest(
        host="escherichia_coli:511145/MG1655",
        cds=gene["cds"],
        host_profile=profile,
        objectives=["usage.cai"],
        constraints=constraints_for(gene["cds"], lamarckian),
        translation_table=11,
        budget=BUDGET,
        seed=seed,
        generations=GENERATIONS,
        segment=SEGMENT,
        stagnation_window=STAGNATION_WINDOW,
        host_digest="sha256:" + "ab" * 32,
    )
    result = design(request, store)
    state = result.job.checkpoint or {}
    members = state.get("members") or []
    problem = ""
    if result.job.state == "failed" and result.job.problem:
        problem = str(result.job.problem.get("code", ""))
    return Trial(
        hypervolume=_hypervolume(result.front),
        diversity=_diversity(members),
        repair_invocations=int(state.get("repairs_attempted") or 0),
        generations_run=int(result.job.generations_run),
        problem_code=problem,
    )


def run_campaign(directory: Path) -> dict[str, Any]:
    """Twenty seeds, two arms, three proteins. Applies decide() unchanged."""
    directory.mkdir(parents=True, exist_ok=True)
    panel = load_panel()
    genes = {row["gene"]: row for row in select_genes(panel)}
    profile = host_profile()
    store = ArtifactStore(directory / "store")
    started = time.perf_counter()
    per_protein = []
    for name in PROTEINS:
        gene = genes[name]
        lamarckian = []
        baldwinian = []
        for seed in range(SEEDS):
            lamarckian.append(run_one(gene, profile, seed, True, store))
            baldwinian.append(run_one(gene, profile, seed, False, store))
        decision = decide(lamarckian, baldwinian)
        per_protein.append(
            {
                "gene": name,
                "codons": gene["codons"],
                "switch_to_baldwinian": decision.switch_to_baldwinian,
                "hv_margin": decision.hv_margin,
                "hv_standard_deviation": decision.hv_standard_deviation,
                "repair_ratio": (
                    None if math.isinf(decision.repair_ratio) else decision.repair_ratio
                ),
                "mean_diversity_lamarckian": decision.mean_diversity_lamarckian,
                "mean_diversity_baldwinian": decision.mean_diversity_baldwinian,
                "mean_repairs_lamarckian": sum(t.repair_invocations for t in lamarckian)
                / len(lamarckian),
                "mean_repairs_baldwinian": sum(t.repair_invocations for t in baldwinian)
                / len(baldwinian),
                "generations_run_min": min(t.generations_run for t in lamarckian + baldwinian),
                "generations_run_max": max(t.generations_run for t in lamarckian + baldwinian),
                "problem_codes": sorted(
                    {t.problem_code for t in lamarckian + baldwinian if t.problem_code}
                ),
            }
        )
    elapsed = time.perf_counter() - started
    # The shipped default changes only if every protein's decision says so.
    # A split result is inconclusive. decide() itself is not edited.
    switches = [row["switch_to_baldwinian"] for row in per_protein]
    if all(switches):
        outcome = "switch_to_baldwinian"
    elif not any(switches):
        outcome = "keep_lamarckian"
    else:
        outcome = "inconclusive"
    document = {
        "proteins": list(PROTEINS),
        "seeds_per_arm": SEEDS,
        "budget": {
            "tier": BUDGET,
            "generations": GENERATIONS,
            "segment": SEGMENT,
            "stagnation_window": STAGNATION_WINDOW,
            "generations_run_min": min(row["generations_run_min"] for row in per_protein),
            "generations_run_max": max(row["generations_run_max"] for row in per_protein),
        },
        "reference_point": list(REFERENCE),
        "senses": SENSES,
        "elapsed_seconds": elapsed,
        "outcome": outcome,
        "per_protein": per_protein,
        "default_changed": False,
    }
    (directory / "campaign.json").write_text(json.dumps(document, indent=2), encoding="utf-8")
    return document
