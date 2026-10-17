#!/usr/bin/env python3
"""Shared pieces for the v2.9.3 design-quality measurements.

The CAI number in the reports is ``benchmarks/validation/cai_independent.py``:
Sharp and Li, geometric mean of relative adaptiveness, stops and the
single-codon families Met and Trp excluded, and a counted weight of 0 making
the mean 0. ``usage.cai`` from the design engine is recorded beside it so a
disagreement is visible.

From the repository root, after ``make build``. The two study scripts import
this module; it is not a measurement by itself.
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmarks" / "validation"))

import cai_independent  # noqa: E402

from codonforge import __version__ as PACKAGE_VERSION  # noqa: E402
from codonforge import _core  # noqa: E402
from codonforge.modelreg import ArtifactStore  # noqa: E402
from codonforge.orchestrator.design import DesignRequest, design, score_sequences  # noqa: E402
from codonforge.orchestrator.structure_diag import (  # noqa: E402
    SPAN_NT,
    TEMPERATURE_C,
    WINDOW_NT,
)
from codonforge.orchestrator.structure_diag import (  # noqa: E402
    evaluate as accessibility_evaluate,
)

HOST = "escherichia_coli:511145/MG1655"
TRANSLATION_TABLE = 11
WEIGHTS_PATH = ROOT / "examples" / "mg1655_ribosomal_usage.json"
ASSEMBLY_PATH = ROOT / "benchmarks" / "validation" / "constraints_assembly.json"
GOLDEN_GATE_PATH = ROOT / "examples" / "golden_gate.json"
BIOBRICK_PATH = ROOT / "examples" / "biobrick_rfc10.json"

# Validation CDS used as the shared protein set. lacZ and yeeJ are included
# because both already appear in the published notes; they are long, and a
# caller that cannot finish them records a timeout rather than a blank.
GENES: list[tuple[str, str]] = [
    ("rpsL", "benchmarks/validation/sequences/rpsL.fa"),
    ("elastin_vpgvg20", "benchmarks/validation/sequences/elp_vpgvg.fa"),
    ("gfp_BBa_E0040", "benchmarks/validation/sequences/BBa_E0040_onestop.fa"),
    ("mrfp_BBa_E1010", "benchmarks/validation/sequences/BBa_E1010.fa"),
    ("mcherry_BBa_J06504", "benchmarks/validation/sequences/BBa_J06504.fa"),
    ("bla", "benchmarks/validation/sequences/bla.fa"),
    ("waaS", "benchmarks/validation/sequences/waaS.fa"),
    ("yagF", "benchmarks/validation/sequences/yagF.fa"),
    ("lacZ", "benchmarks/validation/sequences/lacZ.fa"),
    ("yeeJ", "benchmarks/validation/sequences/yeeJ.fa"),
]

CONSTRAINT_SETS: list[tuple[str, str]] = [
    ("assembly", "benchmarks/validation/constraints_assembly.json"),
    ("golden_gate", "examples/golden_gate.json"),
    ("biobrick_rfc10", "examples/biobrick_rfc10.json"),
]


def git_head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()


def load_weights() -> dict[str, float]:
    return cai_independent.load_ribosomal_weights(WEIGHTS_PATH)


def load_fasta(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    sequence = "".join(line.strip() for line in lines if line and not line.startswith(">"))
    sequence = sequence.upper().replace("U", "T")
    if not sequence or len(sequence) % 3 != 0:
        raise ValueError(f"{path} is not a whole number of codons")
    if any(base not in "ACGT" for base in sequence):
        raise ValueError(f"{path} contains a base other than ACGT")
    return sequence


def load_genes() -> list[dict[str, Any]]:
    genes = []
    for name, relative in GENES:
        cds = load_fasta(ROOT / relative)
        genes.append(
            {
                "gene": name,
                "path": relative,
                "cds": cds,
                "nt": len(cds),
                "codons": len(cds) // 3,
                "protein": cai_independent.translate(cds),
            }
        )
    return genes


def load_constraint_sets() -> list[dict[str, Any]]:
    loaded = []
    for name, relative in CONSTRAINT_SETS:
        document = json.loads((ROOT / relative).read_text(encoding="utf-8"))
        loaded.append({"name": name, "path": relative, "document": document})
    return loaded


def host_profile(weights: dict[str, float]) -> dict[str, Any]:
    return {
        "codon_usage": {
            "subset": "user_supplied",
            "table": {
                codon: {"relative_adaptiveness": weight}
                for codon, weight in sorted(weights.items())
            },
        }
    }


def profile_digest(profile: dict[str, Any]) -> str:
    payload = json.dumps(profile, sort_keys=True, separators=(",", ":")).encode()
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def gc_fraction(seq: str) -> float:
    return cai_independent.gc_fraction(seq)


def independent_cai(seq: str, weights: dict[str, float]) -> dict[str, Any]:
    scored = cai_independent.score(seq, weights)
    return {
        "cai": scored["cai"],
        "counted_codons": scored["counted_codons"],
        "non_top_codons": scored["non_top_codons"],
        "zero_weight_codons": scored["zero_weight_codons"],
        "all_counted_codons_are_top": scored["all_counted_codons_are_top"],
    }


def engine_cai(seq: str, profile: dict[str, Any], digest: str) -> float | None:
    """``usage.cai`` from the same scorer a design job uses."""
    request = DesignRequest(
        host=HOST,
        cds=seq,
        host_profile=profile,
        objectives=["usage.cai"],
        constraints={
            "schema_version": "1.0.0",
            "constraint_set_id": "cst_score_only",
            "constraints": [],
        },
        translation_table=TRANSLATION_TABLE,
        budget="fast",
        seed=0,
        generations=0,
        host_digest=digest,
        job_id="score-only",
    )
    scored = score_sequences(request, [seq])
    value = scored[0].get("usage.cai")
    if value is None:
        return None
    return float(value)


def violation_summary(seq: str, constraints: dict[str, Any]) -> dict[str, Any]:
    report = json.loads(
        _core.constraint_violations(
            json.dumps(constraints),
            seq,
            [(0, len(seq))],
            "assembled",
        )
    )
    counts: dict[str, int] = {}
    for item in report["violations"]:
        if not item["actionable"]:
            continue
        counts[item["constraint_id"]] = counts.get(item["constraint_id"], 0) + 1
    return {
        "feasible": not counts,
        "actionable_by_id": dict(sorted(counts.items())),
        "n_actionable": sum(counts.values()),
    }


def protein_identity(left: str, right: str) -> bool:
    return cai_independent.translate(left) == cai_independent.translate(right)


def accessibility(seq: str) -> float:
    """Post-hoc windowed accessibility penalty. Not a global MFE."""
    return float(accessibility_evaluate([seq])[0])


def empty_constraints(name: str = "cst_empty") -> dict[str, Any]:
    return {
        "schema_version": "1.0.0",
        "constraint_set_id": name,
        "constraints": [],
    }


def run_codonforge(
    cds: str,
    constraints: dict[str, Any],
    profile: dict[str, Any],
    digest: str,
    store: ArtifactStore,
    *,
    objectives: list[str] | None = None,
    budget: str = "standard",
    seed: int = 7,
    generations: int = 4,
    stagnation_window: int = 21,
    job_id: str,
) -> dict[str, Any]:
    """One design. Wall time is the ``design`` call, not the scoring around it."""
    request = DesignRequest(
        host=HOST,
        cds=cds,
        host_profile=profile,
        objectives=list(objectives or ["usage.cai"]),
        constraints=constraints,
        translation_table=TRANSLATION_TABLE,
        budget=budget,
        seed=seed,
        generations=generations,
        segment=max(1, generations) if generations else 1,
        stagnation_window=stagnation_window,
        host_digest=digest,
        job_id=job_id,
    )
    started = time.perf_counter()
    result = design(request, store)
    wall_s = time.perf_counter() - started
    status = result.status
    row: dict[str, Any] = {
        "state": result.job.state,
        "wall_s": wall_s,
        "budget": budget,
        "seed": seed,
        "generations_requested": generations,
        "stagnation_window": stagnation_window,
        "objectives": list(request.objectives),
        "tier1_beam": result.job.inputs.tier1_beam,
        "tier1_seeds": result.job.inputs.tier1_seeds,
        "tier2_block": result.job.inputs.tier2_block,
        "tier2_sweeps": result.job.inputs.tier2_sweeps,
        "population": result.job.inputs.population,
        "warnings": [
            {"code": item.get("code"), "subject": item.get("subject")}
            for item in status.get("warnings") or []
        ],
        "progress": status.get("progress"),
        "problem": status.get("problem"),
    }
    summary = status.get("current_front_summary")
    if summary is not None:
        row["hypervolume_history"] = summary.get("hypervolume_history")
        row["n_nondominated_engine"] = summary.get("n_nondominated")
        row["n_feasible_engine"] = summary.get("n_feasible")
    if result.front is None or not result.front.get("candidates"):
        row["cds"] = None
        return row
    front = result.front
    row["objective_sense"] = front.get("objective_sense")
    row["objective_order"] = front.get("objective_order")
    row["recommended_rule"] = (front.get("recommended") or {}).get("rule")
    candidates = []
    for candidate in front["candidates"]:
        vector = candidate.get("objective_vector")
        sequence = candidate["sequences"][0]["cds"]
        candidates.append(
            {
                "candidate_id": candidate["candidate_id"],
                "origin": (candidate.get("provenance") or {}).get("origin"),
                "generation_found": (candidate.get("provenance") or {}).get("generation_found"),
                "objective_vector": vector,
                "feasible_flag": (candidate.get("constraint_report") or {}).get("feasible"),
                "cds": sequence,
                "sha256": sha256_text(sequence),
            }
        )
    row["candidates"] = candidates
    row["cds"] = candidates[0]["cds"]
    archive = []
    for item in result.job.archive.to_json():
        archive.append(
            {
                "origin": item["origin"],
                "score": item["score"],
                "feasible": item["feasible"],
                "sha256": sha256_text(item["cds"]),
                "objectives": item["objectives"],
                "recorded": item["recorded"],
            }
        )
    row["archive"] = archive
    return row


def revcomp(seq: str) -> str:
    return cai_independent.revcomp(seq)


def local_patterns(constraints: dict[str, Any]) -> dict[str, Any]:
    """Patterns a codon DP can carry, and the rules it cannot.

    Restriction sites, homopolymers and fixed motifs are regular. A GC window
    and a scattered k-mer repeat are not, and they are named rather than
    dropped silently. ``canonical_strand`` false counts a k-mer and its
    reverse complement as different strings; that rule is still global.
    """
    patterns: set[str] = set()
    deferred: list[str] = []
    for constraint in constraints.get("constraints") or []:
        rule = constraint["rule"]
        kind = rule["kind"]
        if kind == "restriction_site":
            site = rule["site"].upper()
            patterns.add(site)
            if rule.get("strand", "both") == "both":
                patterns.add(revcomp(site))
        elif kind == "forbidden_motif":
            motif = rule["motif"].upper()
            patterns.add(motif)
            if rule.get("strand", "both") == "both":
                patterns.add(revcomp(motif))
        elif kind == "homopolymer":
            run = int(rule["max_run"]) + 1
            bases = "ACGT" if rule.get("base", "any") == "any" else rule["base"].upper()
            for base in bases:
                patterns.add(base * run)
        elif kind in {"gc_window", "kmer_repeat"}:
            deferred.append(constraint["id"])
        else:
            deferred.append(constraint["id"] + ":" + kind)
    return {
        "patterns": sorted(patterns),
        "deferred_constraint_ids": deferred,
        "exact_for_full_set": not deferred,
    }


def minimum_positive_weight(weights: dict[str, float]) -> float:
    positive = [weight for weight in weights.values() if weight > 0]
    if not positive:
        raise ValueError("the weight table has no positive codon")
    return min(positive)


def codon_usage_by_amino_acid(weights: dict[str, float]) -> dict[str, dict[str, float]]:
    """Relative adaptiveness grouped by amino acid, zeros left as zeros.

    DNA Chisel's ``MaximizeCAI`` takes ``log(frequency or 0.001)``. Passing
    relative adaptiveness unchanged keeps ``log(w) - log(w_max)`` equal to
    ``log(w)`` when the family maximum is 1, which these ribosomal weights are.
    Renormalising to frequencies would let the 0.001 floor outrank a rare codon.
    """
    table: dict[str, dict[str, float]] = {}
    for codon, amino in cai_independent.AA.items():
        table.setdefault(amino, {})[codon] = float(weights[codon])
    return table


def family_maxima(weights: dict[str, float]) -> dict[str, float]:
    maxima: dict[str, float] = {}
    for codon, amino in cai_independent.AA.items():
        maxima[amino] = max(maxima.get(amino, 0.0), float(weights[codon]))
    return maxima


def log_add(codon: str, weights: dict[str, float]) -> float | None:
    """Natural-log contribution of one codon to the counted CAI sum.

    ``None`` means the codon is a counted zero and cannot sit on a positive
    CAI path. Met, Trp and stops contribute 0 because the independent CAI
    excludes them.
    """
    amino = cai_independent.AA[codon]
    if amino in {"*", "M", "W"}:
        return 0.0
    weight = float(weights[codon])
    if weight == 0.0:
        return None
    if weight < 0:
        raise ValueError(f"negative weight on {codon}")
    return math.log(weight)


def tool_versions() -> dict[str, Any]:
    versions: dict[str, Any] = {
        "codonforge_package": PACKAGE_VERSION,
        "codonforge_core": _core.version(),
        "git_head": git_head(),
        "python": sys.version.split()[0],
        "viennarna_accessibility": {
            "window_nt": WINDOW_NT,
            "span_nt": SPAN_NT,
            "temperature_c": TEMPERATURE_C,
        },
    }
    try:
        import RNA

        versions["viennarna"] = RNA.__version__
    except ImportError:
        versions["viennarna"] = None
    return versions
