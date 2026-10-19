#!/usr/bin/env python3
"""Score DNA Design Suite and DNA Chisel sequences without importing either tool.

The interpreter for this file is Biopython 1.81. That is the last Biopython
release that still ships the corrected Sharp and Li calculator,
``Bio.SeqUtils.CodonAdaptationIndex``. Biopython 1.82 removed it, and the
project environment's Biopython 1.88 does not have it. Global GC and each GC
window use ``Bio.SeqUtils.gc_fraction``. Translation uses ``Bio.Seq.Seq`` with
NCBI table 11.

The reference weights are the relative-adaptiveness values in
``examples/mg1655_ribosomal_usage.json``. ``calculate`` skips ATG and TGG.
Stop codons are omitted from the index, so ``calculate`` skips those too. A
counted codon whose weight is 0 is reported as CAI 0: that is the geometric
mean, and ``calculate`` itself rejects ``log(0)``. This is not
``SharpEcoliIndex`` and not ``usage.cai``.

Hard constraints are read from the constraint JSON and counted with string
scans. A kind this checker does not implement raises, so a rule cannot pass
by being ignored.

From the repository root, with the scoring environment installed:

    benchmarks/gate_ben_001_chisel/.score-venv/bin/python \\
        benchmarks/gate_ben_001_chisel/independent_score.py --self-test
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

from Bio import __version__ as BIOPYTHON_VERSION
from Bio.Restriction import BsaI, BsmBI, EcoRI, NotI, PstI, SpeI, XbaI
from Bio.Seq import Seq
from Bio.SeqUtils import CodonAdaptationIndex, gc_fraction

ROOT = Path(__file__).resolve().parents[2]
WEIGHTS_PATH = ROOT / "examples" / "mg1655_ribosomal_usage.json"
CONSTRAINTS_PATH = ROOT / "benchmarks" / "validation" / "constraints_assembly.json"
STOPS = {"TAA", "TAG", "TGA"}
# Biopython's calculate() skips these two before the index lookup.
METHIONINE_AND_TRYPTOPHAN = {"ATG", "TGG"}
CAI_TIE_ABS = 1e-12
ENZYMES = {
    "BsaI": BsaI,
    "BsmBI": BsmBI,
    "EcoRI": EcoRI,
    "XbaI": XbaI,
    "SpeI": SpeI,
    "PstI": PstI,
    "NotI": NotI,
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_weights(path: Path) -> dict[str, float]:
    document = json.loads(path.read_text(encoding="utf-8"))
    table = document["codon_usage"]["table"]
    weights = {codon: float(row["relative_adaptiveness"]) for codon, row in table.items()}
    if len(weights) != 64 or any(len(codon) != 3 for codon in weights):
        raise SystemExit(f"{path} is not a 64-codon table")
    return weights


def load_constraints(path: Path) -> dict[str, Any]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if not document.get("constraints"):
        raise SystemExit(f"{path} has no constraints")
    return document


def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def translate(seq: str) -> str:
    if len(seq) % 3:
        raise ValueError(f"length {len(seq)} is not a multiple of 3")
    return str(Seq(seq).translate(table=11))


def _index(weights: dict[str, float]) -> CodonAdaptationIndex:
    """An index ``calculate`` can read. Stops are left out on purpose."""
    index = CodonAdaptationIndex(["ATGAAATAA"])
    index.clear()
    for codon, weight in weights.items():
        if codon in STOPS:
            continue
        index[codon] = float(weight)
    return index


def cai_for_cds(seq: str, weights: dict[str, float], index: CodonAdaptationIndex) -> float | None:
    """Sharp and Li geometric mean via Biopython 1.81 ``calculate``.

    Returns 0 when any counted codon has weight 0. Returns None when every
    codon is Met, Trp, or a stop, because ``calculate`` would divide by zero.
    """
    counted = 0
    for offset in range(0, len(seq), 3):
        codon = seq[offset : offset + 3]
        if codon in METHIONINE_AND_TRYPTOPHAN or codon in STOPS:
            continue
        if codon not in weights:
            raise ValueError(f"codon {codon} is not in the reference table")
        if weights[codon] == 0.0:
            return 0.0
        counted += 1
    if counted == 0:
        return None
    return float(index.calculate(seq))


def _overlapping(seq: str, pattern: str) -> int:
    count = 0
    start = 0
    while True:
        found = seq.find(pattern, start)
        if found < 0:
            return count
        count += 1
        start = found + 1


def _motif_hits(seq: str, motif: str, strand: str) -> int:
    if any(base not in "ACGT" for base in motif):
        raise ValueError(
            f"motif {motif} is not a literal ACGT string; this checker will not expand IUPAC"
        )
    hits = _overlapping(seq, motif)
    if strand == "both":
        complement = revcomp(motif)
        if complement != motif:
            hits += _overlapping(seq, complement)
    elif strand != "sense":
        raise ValueError(f"strand {strand} is not sense or both")
    return hits


def _gc_windows(seq: str, window: int, low: float, high: float) -> int:
    if window <= 0 or len(seq) < window:
        return 0
    violations = 0
    for start in range(0, len(seq) - window + 1):
        fraction = gc_fraction(seq[start : start + window])
        if fraction < low or fraction > high:
            violations += 1
    return violations


def _homopolymer_runs(seq: str, max_run: int, base: str) -> int:
    allowed = "ACGT" if base == "any" else base
    if any(item not in "ACGT" for item in allowed):
        raise ValueError(f"homopolymer base {base} is not ACGT or any")
    violations = 0
    index = 0
    while index < len(seq):
        end = index + 1
        while end < len(seq) and seq[end] == seq[index]:
            end += 1
        if seq[index] in allowed and (end - index) > max_run:
            violations += 1
        index = end
    return violations


def _kmer_offenders(seq: str, rule: dict[str, Any]) -> int:
    k = int(rule["k"])
    limit = int(rule["max_occurrences"])
    overlapping = bool(rule.get("count_overlapping", True))
    canonical = bool(rule.get("canonical_strand", False))
    if k <= 0 or len(seq) < k:
        return 0
    step = 1 if overlapping else k
    counts: dict[str, int] = {}
    for start in range(0, len(seq) - k + 1, step):
        kmer = seq[start : start + k]
        if canonical:
            complement = revcomp(kmer)
            if complement < kmer:
                kmer = complement
        counts[kmer] = counts.get(kmer, 0) + 1
    return sum(1 for count in counts.values() if count > limit)


def constraint_report(seq: str, document: dict[str, Any]) -> dict[str, Any]:
    """One plain count per constraint id. ``passes`` is zero violations."""
    by_id: dict[str, Any] = {}
    for constraint in document["constraints"]:
        rule = constraint["rule"]
        kind = rule["kind"]
        if kind == "gc_window":
            violations = _gc_windows(
                seq, int(rule["window_nt"]), float(rule["min"]), float(rule["max"])
            )
        elif kind in {"restriction_site", "forbidden_motif"}:
            motif = rule.get("site", rule.get("motif", "")).upper()
            hits = _motif_hits(seq, motif, rule.get("strand", "both"))
            violations = max(0, hits - int(rule.get("max_occurrences", 0)))
        elif kind == "homopolymer":
            violations = _homopolymer_runs(seq, int(rule["max_run"]), rule.get("base", "any"))
        elif kind == "kmer_repeat":
            violations = _kmer_offenders(seq, rule)
        else:
            raise ValueError(
                f"{constraint['id']}: kind {kind} has no plain checker; "
                "it is not treated as satisfied"
            )
        by_id[constraint["id"]] = {
            "kind": kind,
            "n_violations": violations,
            "passes": violations == 0,
        }
    total = sum(item["n_violations"] for item in by_id.values())
    return {
        "by_id": by_id,
        "n_violations": total,
        "passes": total == 0 and bool(by_id),
        "ids": [constraint["id"] for constraint in document["constraints"]],
    }


def score_sequence(
    seq: str | None,
    *,
    input_cds: str,
    weights: dict[str, float],
    index: CodonAdaptationIndex,
    constraints: dict[str, Any],
) -> dict[str, Any]:
    if not seq:
        return {
            "cai": None,
            "gc": None,
            "translation": None,
            "translation_match": False,
            "constraints": None,
            "feasible": False,
            "score_error": "no sequence",
        }
    sequence = seq.upper().replace("U", "T")
    if any(base not in "ACGT" for base in sequence):
        return {
            "cai": None,
            "gc": None,
            "translation": None,
            "translation_match": False,
            "constraints": None,
            "feasible": False,
            "score_error": "sequence is not ACGT",
        }
    try:
        protein = translate(sequence)
        reference = translate(input_cds)
    except ValueError as exc:
        return {
            "cai": None,
            "gc": float(gc_fraction(sequence)) if sequence else None,
            "translation": None,
            "translation_match": False,
            "constraints": constraint_report(sequence, constraints),
            "feasible": False,
            "score_error": str(exc),
        }
    try:
        value = cai_for_cds(sequence, weights, index)
    except (ValueError, ZeroDivisionError) as exc:
        value = None
        cai_error = str(exc)
    else:
        cai_error = None
    checked = constraint_report(sequence, constraints)
    match = protein == reference
    return {
        "cai": value,
        "gc": float(gc_fraction(sequence)),
        "translation": protein,
        "translation_match": match,
        "constraints": checked,
        "feasible": bool(match and checked["passes"] and value is not None),
        "score_error": cai_error,
    }


def _relation(left: float | None, right: float | None, *, tolerance: float = 0.0) -> str:
    if left is None or right is None:
        return "missing"
    if abs(left - right) <= tolerance:
        return "tie"
    if left > right:
        return "dna_design_suite_higher"
    return "dnachisel_higher"


def _names(genes: list[str]) -> str:
    if not genes:
        return "none"
    return ", ".join(genes)


def compare_rows(document: dict[str, Any]) -> dict[str, Any]:
    """Paired counts. Panel order, not sorted by a score."""
    order = document["meta"]["gene_order"]
    paired: dict[str, dict[str, dict[str, Any]]] = {}
    for row in document["rows"]:
        paired.setdefault(row["gene"], {})[row["tool"]] = row
    ids = document["meta"]["constraint_ids"]
    cai_both: list[dict[str, Any]] = []
    cai_feasible: list[dict[str, Any]] = []
    gc_pairs: list[dict[str, Any]] = []
    per_constraint = {
        constraint_id: {
            "suite_fail_chisel_pass": [],
            "chisel_fail_suite_pass": [],
            "both_fail": [],
            "both_pass": [],
            "suite_violations": 0,
            "chisel_violations": 0,
        }
        for constraint_id in ids
    }
    translation = {
        "suite_mismatch": [],
        "chisel_mismatch": [],
        "suite_returned": 0,
        "chisel_returned": 0,
        "suite_match": 0,
        "chisel_match": 0,
    }
    feasible = {
        "suite_only": [],
        "chisel_only": [],
        "both": [],
        "neither": [],
    }
    returned = {"suite_only": [], "chisel_only": [], "both": [], "neither": []}
    solver = {
        "suite_faster": [],
        "chisel_faster": [],
        "tie": [],
        "suite_sum": 0.0,
        "chisel_sum": 0.0,
    }
    wall = {"suite_faster": [], "chisel_faster": [], "tie": [], "suite_sum": 0.0, "chisel_sum": 0.0}
    product_disagree = []
    usage_deltas = []

    def _timed(bucket: dict[str, Any], gene: str, left: float | None, right: float | None) -> None:
        if left is None or right is None:
            return
        bucket["suite_sum"] += left
        bucket["chisel_sum"] += right
        if left < right:
            bucket["suite_faster"].append(gene)
        elif right < left:
            bucket["chisel_faster"].append(gene)
        else:
            bucket["tie"].append(gene)

    for gene in order:
        pair = paired[gene]
        suite = pair["dna_design_suite"]
        chisel = pair["dnachisel"]
        suite_seq = bool(suite.get("sequence"))
        chisel_seq = bool(chisel.get("sequence"))
        if suite_seq and chisel_seq:
            returned["both"].append(gene)
        elif suite_seq:
            returned["suite_only"].append(gene)
        elif chisel_seq:
            returned["chisel_only"].append(gene)
        else:
            returned["neither"].append(gene)
        if suite_seq:
            translation["suite_returned"] += 1
            if suite["translation_match"]:
                translation["suite_match"] += 1
            else:
                translation["suite_mismatch"].append(gene)
        if chisel_seq:
            translation["chisel_returned"] += 1
            if chisel["translation_match"]:
                translation["chisel_match"] += 1
            else:
                translation["chisel_mismatch"].append(gene)
        suite_ok = bool(suite.get("feasible"))
        chisel_ok = bool(chisel.get("feasible"))
        if suite_ok and chisel_ok:
            feasible["both"].append(gene)
        elif suite_ok:
            feasible["suite_only"].append(gene)
        elif chisel_ok:
            feasible["chisel_only"].append(gene)
        else:
            feasible["neither"].append(gene)
        if suite_seq and chisel_seq:
            cai_both.append(
                {
                    "gene": gene,
                    "relation": _relation(suite["cai"], chisel["cai"], tolerance=CAI_TIE_ABS),
                    "exact": suite["cai"] == chisel["cai"],
                    "dna_design_suite": suite["cai"],
                    "dnachisel": chisel["cai"],
                }
            )
            gc_pairs.append(
                {
                    "gene": gene,
                    "relation": _relation(suite["gc"], chisel["gc"]),
                    "dna_design_suite": suite["gc"],
                    "dnachisel": chisel["gc"],
                }
            )
        if suite_ok and chisel_ok:
            cai_feasible.append(
                {
                    "gene": gene,
                    "relation": _relation(suite["cai"], chisel["cai"], tolerance=CAI_TIE_ABS),
                    "exact": suite["cai"] == chisel["cai"],
                    "dna_design_suite": suite["cai"],
                    "dnachisel": chisel["cai"],
                }
            )
        for constraint_id in ids:
            slot = per_constraint[constraint_id]
            if suite_seq:
                slot["suite_violations"] += suite["constraints"]["by_id"][constraint_id][
                    "n_violations"
                ]
            if chisel_seq:
                slot["chisel_violations"] += chisel["constraints"]["by_id"][constraint_id][
                    "n_violations"
                ]
            if not suite_seq or not chisel_seq:
                continue
            suite_pass = suite["constraints"]["by_id"][constraint_id]["passes"]
            chisel_pass = chisel["constraints"]["by_id"][constraint_id]["passes"]
            if suite_pass and chisel_pass:
                slot["both_pass"].append(gene)
            elif not suite_pass and not chisel_pass:
                slot["both_fail"].append(gene)
            elif suite_pass:
                slot["chisel_fail_suite_pass"].append(gene)
            else:
                slot["suite_fail_chisel_pass"].append(gene)
        _timed(solver, gene, suite.get("solver_s"), chisel.get("solver_s"))
        _timed(wall, gene, suite.get("wall_s"), chisel.get("wall_s"))
        for row in (suite, chisel):
            diagnostic = row.get("diagnostic_product_checker")
            if row.get("sequence") and diagnostic is not None:
                plain = bool(row["constraints"]["passes"])
                if bool(diagnostic.get("feasible")) != plain:
                    product_disagree.append(
                        {
                            "gene": gene,
                            "tool": row["tool"],
                            "plain_passes": plain,
                            "product_feasible": diagnostic.get("feasible"),
                        }
                    )
            if row.get("cai") is not None and row.get("diagnostic_usage_cai") is not None:
                usage_deltas.append(abs(float(row["diagnostic_usage_cai"]) - float(row["cai"])))

    def _tally(pairs: list[dict[str, Any]]) -> dict[str, Any]:
        groups = {"dna_design_suite_higher": [], "dnachisel_higher": [], "tie": [], "missing": []}
        for item in pairs:
            groups[item["relation"]].append(item["gene"])
        return {
            "dna_design_suite_higher": groups["dna_design_suite_higher"],
            "dnachisel_higher": groups["dnachisel_higher"],
            "tie": groups["tie"],
            "n_dna_design_suite_higher": len(groups["dna_design_suite_higher"]),
            "n_dnachisel_higher": len(groups["dnachisel_higher"]),
            "n_tie": len(groups["tie"]),
            "n": len(pairs),
        }

    def _mean(pairs: list[dict[str, Any]], key: str) -> float | None:
        values = [item[key] for item in pairs if item[key] is not None]
        if not values:
            return None
        return float(sum(values) / len(values))

    return {
        "cai_tie_abs": CAI_TIE_ABS,
        "cai_both_returned": _tally(cai_both),
        "cai_both_feasible": _tally(cai_feasible),
        "cai_both_returned_mean": {
            "dna_design_suite": _mean(cai_both, "dna_design_suite"),
            "dnachisel": _mean(cai_both, "dnachisel"),
        },
        "cai_both_feasible_mean": {
            "dna_design_suite": _mean(cai_feasible, "dna_design_suite"),
            "dnachisel": _mean(cai_feasible, "dnachisel"),
        },
        "gc_both_returned": _tally(gc_pairs),
        "gc_note": "Global GC is not an objective. A higher GC is not a better design.",
        "translation": translation,
        "feasible": {key: value for key, value in feasible.items()},
        "n_feasible": {key: len(value) for key, value in feasible.items()},
        "returned": returned,
        "n_returned": {key: len(value) for key, value in returned.items()},
        "per_constraint": per_constraint,
        "solver_s": solver,
        "wall_s": wall,
        "product_checker_disagreements": product_disagree,
        "max_abs_usage_cai_minus_biopython": max(usage_deltas) if usage_deltas else None,
        "n_usage_cai_compared": len(usage_deltas),
    }


def _count_line(label: str, higher_us: int, higher_them: int, ties: int) -> str:
    return f"| {label} | {higher_us} | {higher_them} | {ties} |"


def render_report(document: dict[str, Any]) -> str:
    meta = document["meta"]
    comparison = document["comparison"]
    order = meta["gene_order"]
    paired: dict[str, dict[str, dict[str, Any]]] = {}
    for row in document["rows"]:
        paired.setdefault(row["gene"], {})[row["tool"]] = row
    versions = meta["versions"]
    lines: list[str] = []
    lines.append("# DNA Design Suite and DNA Chisel on the GATE-BEN-001 proteins")
    lines.append("")
    lines.append(
        "One head-to-head, measured at commit "
        f"`{meta['git_commit']}`. The sequences and the scores are in "
        "`docs/benchmarks/gate-ben-001-dnachisel-rows.json`, written by this "
        "run. CAI in this note is the Sharp and Li geometric mean of a fixed "
        "weight table. It is a sequence statistic. It is not a predicted "
        "expression level, and this note makes no wet-lab claim."
    )
    lines.append("")
    lines.append(
        "The published `GATE-BEN-001` claim is a different comparison: the "
        "same 50 proteins and the same ribosomal weights, an empty constraint "
        "list, and DNA Design Suite against `cai_greedy`. This note does not "
        "replace that claim and does not extend it. The hard constraints "
        "here are `benchmarks/validation/constraints_assembly.json`, given "
        "to both tools."
    )
    lines.append("")
    lines.append("## Versions, seeds, and the machine")
    lines.append("")
    lines.append("| piece | pin |")
    lines.append("| --- | --- |")
    lines.append(f"| git commit measured | `{meta['git_commit']}` |")
    lines.append(f"| DNA Design Suite package | {versions['dna_design_suite']} |")
    lines.append(f"| DNA Design Suite core | {versions['core']} |")
    lines.append(f"| DNA Chisel | {versions['dnachisel']} |")
    lines.append(f"| scoring Biopython | {versions['biopython_scorer']} |")
    lines.append(f"| project-environment Biopython | {versions['biopython_project']} |")
    lines.append(f"| numpy | {versions['numpy']} |")
    lines.append(f"| Python, design | {versions['python_design']} |")
    lines.append(f"| Python, scoring | {versions['python_scorer']} |")
    lines.append(f"| DNA Design Suite seed | {meta['seeds']['dna_design_suite']} |")
    lines.append(f"| DNA Chisel numpy seed | {meta['seeds']['dnachisel_numpy']} |")
    lines.append(f"| DNA Chisel `PYTHONHASHSEED` | {meta['seeds']['pythonhashseed']} |")
    lines.append("")
    machine = meta["machine"]
    lines.append(
        f"Machine: {machine['cpu']}, {machine['nproc']} CPUs, "
        f"{machine['mem_total_kib']} KiB RAM, {machine['platform']}. "
        "Genes ran one at a time. DNA Design Suite time is `time.perf_counter` "
        "around `design` in an already imported process, after one untimed "
        "warmup call. DNA Chisel solver time is `time.perf_counter` inside a "
        "fresh interpreter, around problem construction, "
        "`resolve_constraints`, and `optimize`. DNA Chisel wall time is the "
        "parent clock around that process, so it includes interpreter startup. "
        "Neither clock is a packaging claim about another machine."
    )
    lines.append("")
    lines.append("## What both tools were given")
    lines.append("")
    lines.append(
        f"Fifty CDS from `{meta['panel_id']}`, accession `{meta['accession']}`, "
        f"selected-CDS sha256 `{meta['selected_cds_sha256']}`. Host "
        f"`{meta['host']}`, translation table {meta['translation_table']}. "
        "Weights are `examples/mg1655_ribosomal_usage.json`, sha256 "
        f"`{meta['weights_sha256']}`. They match `ribosomal_usage()` on the "
        "shipped CDS gzip. Every amino-acid family in that table has maximum "
        f"1. The smallest positive weight is {meta['min_positive_weight']}. "
        f"Weight 0: {', '.join(meta['zero_weight_codons'])}."
    )
    lines.append("")
    lines.append(
        "DNA Design Suite settings are the `GATE-BEN-001` runner's settings, "
        f"with the assembly file in place of the empty list: budget "
        f"`{meta['suite_job']['budget']}`, generations "
        f"{meta['suite_job']['generations']}, seed {meta['seeds']['dna_design_suite']}, "
        f"segment {meta['suite_job']['segment']}, stagnation window "
        f"{meta['suite_job']['stagnation_window']}, objective `usage.cai` only. "
        "Budget `fast` runs tier 1 and tier 2. It does not run the "
        "generational tier, so the generation count is not extra search."
    )
    lines.append("")
    lines.append(
        "DNA Chisel settings are library defaults apart from the seed and "
        "a silent logger: `max_random_iters` "
        f"{meta['chisel_job']['max_random_iters']}, `mutations_per_iteration` "
        f"{meta['chisel_job']['mutations_per_iteration']}, "
        "`optimization_stagnation_tolerance` "
        f"{meta['chisel_job']['optimization_stagnation_tolerance']}, "
        "`randomization_threshold` "
        f"{meta['chisel_job']['randomization_threshold']}. The objective is "
        '`CodonOptimize(method="use_best_codon")` on this weight table. '
        "DNA Chisel's own score is `N * log` over every codon, including "
        "Met, Trp, and stops, and it replaces a zero with 0.001 before the "
        "logarithm. The smallest positive weight is above 0.001, so that "
        "floor does not reorder the positive codons. The CAI column below "
        "does not use that score."
    )
    lines.append("")
    lines.append("## Constraint translation")
    lines.append("")
    lines.append("| id | rule | DNA Chisel specification |")
    lines.append("| --- | --- | --- |")
    for item in meta["constraint_translation"]["translated"]:
        lines.append(f"| `{item['id']}` | {item['rule']} | {item['spec']} |")
    lines.append("")
    dropped = meta["constraint_translation"]["not_translated"]
    if dropped:
        lines.append(
            "These shared constraints were not sent, because no DNA Chisel "
            "specification is the same rule:"
        )
        lines.append("")
        for item in dropped:
            lines.append(f"- `{item['id']}`: {item['reason']}")
        lines.append("")
    else:
        lines.append(
            "Every constraint in the shared file was sent. None of those rules was dropped."
        )
        lines.append("")
    lines.append(
        "The file's `repair_policy` is empty. It is search configuration, "
        "not a sequence predicate, and it was not sent."
    )
    lines.append("")
    lines.append("Schema kinds that are not in the shared file were not given to either tool:")
    lines.append("")
    for item in meta["schema_kinds_not_in_file"]:
        lines.append(f"- `{item['kind']}`: {item['reason']}")
    lines.append("")
    lines.append("## How the sequences were scored")
    lines.append("")
    lines.append(meta["scoring"]["cai"])
    lines.append("")
    lines.append(meta["scoring"]["gc"])
    lines.append("")
    lines.append(meta["scoring"]["constraints"])
    lines.append("")
    lines.append(meta["scoring"]["translation"])
    lines.append("")
    if comparison["product_checker_disagreements"]:
        lines.append(
            "The product's constraint counter disagreed with the plain "
            f"checker on {len(comparison['product_checker_disagreements'])} "
            "sequences. The tables use the plain checker. The disagreements "
            "are in the JSON under `comparison.product_checker_disagreements`."
        )
    else:
        lines.append(
            "On every returned sequence the product's pass/fail agreed with "
            "the plain checker. The tables still use the plain checker. "
            f"Across {comparison['n_usage_cai_compared']} sequences the largest "
            "absolute gap between `usage.cai` and this Biopython CAI is "
            f"{comparison['max_abs_usage_cai_minus_biopython']}. That gap is "
            "recorded. It is not the headline."
        )
    lines.append("")
    lines.append("## Losses and ties")
    lines.append("")
    lines.append(
        "A loss is a metric on which DNA Chisel is strictly better on at "
        "least one gene, or a total on which DNA Chisel is strictly better. "
        "A tie is equality, with CAI equality at an absolute difference of "
        f"at most {CAI_TIE_ABS}. Global GC is reported as higher, lower, or "
        "equal, and is not called a win or a loss. Genes are named in panel "
        "order."
    )
    lines.append("")
    cai_any = comparison["cai_both_returned"]
    cai_ok = comparison["cai_both_feasible"]
    lines.append(
        "CAI where both tools returned a sequence, whether or not the "
        f"sequence passes the hard constraints: DNA Design Suite higher on "
        f"{cai_any['n_dna_design_suite_higher']}, DNA Chisel higher on "
        f"{cai_any['n_dnachisel_higher']}, tie on {cai_any['n_tie']} "
        f"(of {cai_any['n']})."
    )
    lines.append("")
    if cai_any["dnachisel_higher"]:
        lines.append(
            "DNA Chisel CAI higher, not conditioned on feasibility: "
            + _names(cai_any["dnachisel_higher"])
            + "."
        )
        lines.append("")
    if cai_any["tie"]:
        lines.append("CAI tie, not conditioned on feasibility: " + _names(cai_any["tie"]) + ".")
        lines.append("")
    lines.append(
        "CAI where both sequences are feasible (translation matches and "
        f"every hard constraint passes): DNA Design Suite higher on "
        f"{cai_ok['n_dna_design_suite_higher']}, DNA Chisel higher on "
        f"{cai_ok['n_dnachisel_higher']}, tie on {cai_ok['n_tie']} "
        f"(of {cai_ok['n']}). Paired means on that set: DNA Design Suite "
        f"{comparison['cai_both_feasible_mean']['dna_design_suite']}, "
        f"DNA Chisel {comparison['cai_both_feasible_mean']['dnachisel']}."
    )
    lines.append("")
    if cai_ok["dnachisel_higher"]:
        lines.append(
            "DNA Chisel CAI higher among feasible pairs: "
            + _names(cai_ok["dnachisel_higher"])
            + "."
        )
        lines.append("")
    if cai_ok["tie"]:
        lines.append("CAI tie among feasible pairs: " + _names(cai_ok["tie"]) + ".")
        lines.append("")
    feasible = comparison["feasible"]
    lines.append(
        "Feasible sequences: both "
        f"{len(feasible['both'])}, DNA Design Suite only "
        f"{len(feasible['suite_only'])}, DNA Chisel only "
        f"{len(feasible['chisel_only'])}, neither "
        f"{len(feasible['neither'])}."
    )
    lines.append("")
    if feasible["chisel_only"]:
        lines.append(
            "DNA Chisel feasible and DNA Design Suite not: " + _names(feasible["chisel_only"]) + "."
        )
        lines.append("")
    if not feasible["chisel_only"] and len(feasible["both"]) == len(order):
        lines.append("Feasibility is a tie on all 50: both tools are feasible on every gene.")
        lines.append("")
    translation = comparison["translation"]
    lines.append(
        "Translation back to the input protein: DNA Design Suite "
        f"{translation['suite_match']} of {translation['suite_returned']} "
        "returned sequences, DNA Chisel "
        f"{translation['chisel_match']} of {translation['chisel_returned']}."
    )
    lines.append("")
    if translation["suite_mismatch"]:
        lines.append(
            "DNA Design Suite translation mismatch: " + _names(translation["suite_mismatch"]) + "."
        )
        lines.append("")
    if translation["chisel_mismatch"]:
        lines.append(
            "DNA Chisel translation mismatch: " + _names(translation["chisel_mismatch"]) + "."
        )
        lines.append("")
    if not translation["suite_mismatch"] and not translation["chisel_mismatch"]:
        lines.append(
            "Translation is a tie on every returned sequence: each one "
            "translates back to its input protein."
        )
        lines.append("")
    returned = comparison["returned"]
    if returned["chisel_only"] or returned["suite_only"] or returned["neither"]:
        lines.append(
            "Returned a sequence: both "
            f"{len(returned['both'])}, DNA Design Suite only "
            f"{len(returned['suite_only'])}, DNA Chisel only "
            f"{len(returned['chisel_only'])}, neither {len(returned['neither'])}."
        )
        if returned["chisel_only"]:
            lines.append(
                "DNA Chisel returned a sequence and DNA Design Suite did not: "
                + _names(returned["chisel_only"])
                + "."
            )
        if returned["suite_only"]:
            lines.append(
                "DNA Design Suite returned a sequence and DNA Chisel did not: "
                + _names(returned["suite_only"])
                + "."
            )
        if returned["neither"]:
            lines.append("Neither tool returned a sequence: " + _names(returned["neither"]) + ".")
        lines.append("")
    lines.append("Per hard constraint, on genes where both tools returned a sequence:")
    lines.append("")
    lines.append(
        "| constraint | both pass | suite fails, Chisel passes | "
        "Chisel fails, suite passes | both fail | suite violation sum | "
        "Chisel violation sum |"
    )
    lines.append("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for constraint_id, slot in comparison["per_constraint"].items():
        lines.append(
            f"| `{constraint_id}` | {len(slot['both_pass'])} | "
            f"{len(slot['suite_fail_chisel_pass'])} | "
            f"{len(slot['chisel_fail_suite_pass'])} | {len(slot['both_fail'])} | "
            f"{slot['suite_violations']} | {slot['chisel_violations']} |"
        )
    lines.append("")
    for constraint_id, slot in comparison["per_constraint"].items():
        if slot["suite_fail_chisel_pass"]:
            lines.append(
                f"`{constraint_id}`, DNA Design Suite fails and DNA Chisel passes: "
                + _names(slot["suite_fail_chisel_pass"])
                + "."
            )
            lines.append("")
        if slot["suite_violations"] > slot["chisel_violations"]:
            lines.append(
                f"`{constraint_id}` violation sum is higher for DNA Design Suite "
                f"({slot['suite_violations']} against {slot['chisel_violations']})."
            )
            lines.append("")
        elif slot["suite_violations"] == slot["chisel_violations"]:
            lines.append(f"`{constraint_id}` violation sums tie at {slot['suite_violations']}.")
            lines.append("")
    gc = comparison["gc_both_returned"]
    lines.append(
        "Global GC, both returned a sequence, not an objective: DNA Design "
        f"Suite higher on {gc['n_dna_design_suite_higher']}, DNA Chisel higher "
        f"on {gc['n_dnachisel_higher']}, equal on {gc['n_tie']}."
    )
    lines.append("")
    if gc["tie"]:
        lines.append("Global GC equal: " + _names(gc["tie"]) + ".")
        lines.append("")
    solver = comparison["solver_s"]
    wall = comparison["wall_s"]
    lines.append(
        "Solver time, sum of the 50 genes: DNA Design Suite "
        f"{solver['suite_sum']} s, DNA Chisel {solver['chisel_sum']} s. "
        f"DNA Design Suite faster on {len(solver['suite_faster'])} genes, "
        f"DNA Chisel faster on {len(solver['chisel_faster'])}, "
        f"equal on {len(solver['tie'])}."
    )
    lines.append("")
    if solver["chisel_faster"]:
        lines.append("DNA Chisel solver time faster: " + _names(solver["chisel_faster"]) + ".")
        lines.append("")
    if solver["tie"]:
        lines.append("Solver time equal: " + _names(solver["tie"]) + ".")
        lines.append("")
    if solver["suite_sum"] > solver["chisel_sum"]:
        lines.append(
            "On total solver time, DNA Design Suite is slower "
            f"({solver['suite_sum']} s against {solver['chisel_sum']} s)."
        )
        lines.append("")
    elif solver["suite_sum"] == solver["chisel_sum"]:
        lines.append("Total solver time is a tie.")
        lines.append("")
    lines.append(
        "Wall time, which includes DNA Chisel's fresh interpreter and does "
        f"not include a fresh interpreter for DNA Design Suite: sums "
        f"{wall['suite_sum']} s and {wall['chisel_sum']} s. DNA Chisel's wall "
        f"is lower on {len(wall['chisel_faster'])} genes. That comparison is "
        "not the same boundary on both sides; it is reported anyway."
    )
    lines.append("")
    lines.append("## Per gene")
    lines.append("")
    lines.append(
        "CAI and GC are the independent scores. Feasible means the "
        "translation matches and every hard constraint passes. Solver "
        "seconds are the clocks defined above."
    )
    lines.append("")
    lines.append(
        "| gene | nt | suite CAI | Chisel CAI | CAI | suite GC | Chisel GC | "
        "suite feasible | Chisel feasible | suite violations | "
        "Chisel violations | suite solver s | Chisel solver s |"
    )
    lines.append(
        "| --- | ---: | ---: | ---: | --- | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: |"
    )
    for gene in order:
        suite = paired[gene]["dna_design_suite"]
        chisel = paired[gene]["dnachisel"]
        relation = "missing"
        if suite.get("sequence") and chisel.get("sequence"):
            relation = _relation(suite["cai"], chisel["cai"], tolerance=CAI_TIE_ABS)
            relation = {
                "dna_design_suite_higher": "suite higher",
                "dnachisel_higher": "Chisel higher",
                "tie": "tie",
            }[relation]
        lines.append(
            (
                "| {gene} | {nt} | {sc} | {cc} | {rel} | {sg} | {cg} | "
                "{sf} | {cf} | {sv} | {cv} | {ss} | {cs} |"
            ).format(
                gene=gene,
                nt=suite["nt"],
                sc=_cell(suite.get("cai")),
                cc=_cell(chisel.get("cai")),
                rel=relation,
                sg=_cell(suite.get("gc")),
                cg=_cell(chisel.get("gc")),
                sf=suite.get("feasible"),
                cf=chisel.get("feasible"),
                sv=_violations(suite),
                cv=_violations(chisel),
                ss=_cell(suite.get("solver_s")),
                cs=_cell(chisel.get("solver_s")),
            )
        )
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```sh")
    lines.append("make setup && make build")
    lines.append(".venv/bin/pip install 'dnachisel==3.2.16'")
    lines.append("python3 -m virtualenv benchmarks/gate_ben_001_chisel/.score-venv")
    lines.append("benchmarks/gate_ben_001_chisel/.score-venv/bin/pip install 'biopython==1.81'")
    lines.append("benchmarks/gate_ben_001_chisel/.score-venv/bin/python \\")
    lines.append("    benchmarks/gate_ben_001_chisel/independent_score.py --self-test")
    lines.append(".venv/bin/python benchmarks/gate_ben_001_chisel/run_benchmark.py \\")
    lines.append("    --out docs/benchmarks/gate-ben-001-dnachisel-rows.json \\")
    lines.append("    --doc docs/benchmarks/gate-ben-001-dnachisel.md")
    lines.append("```")
    lines.append("")
    lines.append(meta.get("negative_control", ""))
    lines.append("")
    return "\n".join(lines)


def _cell(value: Any) -> str:
    if value is None:
        return "—"
    return str(value)


def _violations(row: dict[str, Any]) -> str:
    checked = row.get("constraints")
    if not checked:
        return "—"
    return str(checked["n_violations"])


def self_test() -> None:
    weights = load_weights(WEIGHTS_PATH)
    constraints = load_constraints(CONSTRAINTS_PATH)
    index = _index(weights)
    if weights["AAA"] != 1.0 or weights["AAG"] <= 0 or weights["AGG"] != 0.0:
        raise SystemExit("the ribosomal table does not have the weights this self-test names")
    top_codons = ["ATG", "AAA", "AAC", "TTC", "TAC", "GTT", "TCT", "TAA"]
    for codon in top_codons:
        if codon in STOPS or codon in METHIONINE_AND_TRYPTOPHAN:
            continue
        if weights[codon] != 1.0:
            raise SystemExit(f"{codon} weight is {weights[codon]}, the all-top fixture needs 1")
    top = cai_for_cds("".join(top_codons), weights, index)
    if top != 1.0:
        raise SystemExit(f"all-top CAI is {top}, not 1")
    mixed = cai_for_cds("ATGAAAAAGTAA", weights, index)
    expected = math.exp((math.log(1.0) + math.log(weights["AAG"])) / 2.0)
    if mixed is None or abs(mixed - expected) > 1e-12:
        raise SystemExit(f"two-codon CAI is {mixed}, expected {expected}")
    if mixed is not None and abs(mixed - weights["AAG"]) <= 1e-9:
        raise SystemExit("CAI used the n-1 divisor from the deprecated CodonUsage class")
    tagged = cai_for_cds("ATGAAATAG", weights, index)
    if tagged != 1.0:
        raise SystemExit(f"a TAG stop changed CAI to {tagged}; stops are supposed to be excluded")
    zero = cai_for_cds("ATGAGGTAA", weights, index)
    if zero != 0.0:
        raise SystemExit(f"an AGG codon scored {zero}, not 0")
    if abs(gc_fraction("G" * 20 + "A" * 30) - 0.4) > 0:
        raise SystemExit("gc_fraction did not report 20/50 as 0.4")
    low = "G" * 19 + "A" * 31
    high = "G" * 31 + "A" * 19
    edge_low = "G" * 20 + "A" * 30
    edge_high = "G" * 30 + "A" * 20
    if _gc_windows(edge_low, 50, 0.4, 0.6) != 0 or _gc_windows(edge_high, 50, 0.4, 0.6) != 0:
        raise SystemExit("a window at exactly 0.4 or 0.6 was counted as a violation")
    if _gc_windows(low, 50, 0.4, 0.6) != 1 or _gc_windows(high, 50, 0.4, 0.6) != 1:
        raise SystemExit("a window outside 0.4-0.6 was accepted")
    for enzyme, site in {
        "BsaI": "GGTCTC",
        "BsmBI": "CGTCTC",
        "EcoRI": "GAATTC",
        "XbaI": "TCTAGA",
        "SpeI": "ACTAGT",
        "PstI": "CTGCAG",
        "NotI": "GCGGCCGC",
    }.items():
        if ENZYMES[enzyme].site != site:
            raise SystemExit(f"Biopython {enzyme} site is {ENZYMES[enzyme].site}, not {site}")
    forward = constraint_report("GGTCTC", constraints)
    antisense = constraint_report("GAGACC", constraints)
    if forward["by_id"]["c_bsai"]["n_violations"] < 1:
        raise SystemExit("forward BsaI GGTCTC was not counted")
    if antisense["by_id"]["c_bsai"]["n_violations"] < 1:
        raise SystemExit("BsaI reverse complement GAGACC was not counted")
    clean = constraint_report("AAAAAA", constraints)
    if clean["by_id"]["c_homopolymer"]["n_violations"] != 0:
        raise SystemExit("a run of 6 was counted against a limit of 6")
    if clean["by_id"]["c_bsai"]["n_violations"] != 0:
        raise SystemExit("AAAAAA was counted as BsaI")
    long_run = constraint_report("AAAAAAA", constraints)
    if long_run["by_id"]["c_homopolymer"]["n_violations"] < 1:
        raise SystemExit("a run of 7 was accepted")
    if constraint_report("AGGAGG", constraints)["by_id"]["c_sd"]["n_violations"] < 1:
        raise SystemExit("AGGAGG was accepted")
    if constraint_report("CCTCCT", constraints)["by_id"]["c_sd"]["n_violations"] < 1:
        raise SystemExit("Shine-Dalgarno reverse complement CCTCCT was not counted")
    repeated = "ACGTACGTACGT" * 2
    if constraint_report(repeated, constraints)["by_id"]["c_kmer12"]["n_violations"] < 1:
        raise SystemExit("a duplicated 12-mer was accepted")
    once = "ACGTACGTACGT" + "TTTTAAAA"
    if constraint_report(once, constraints)["by_id"]["c_kmer12"]["n_violations"] != 0:
        raise SystemExit("a single 12-mer was counted as a repeat")
    if constraint_report("A" * 13, constraints)["by_id"]["c_kmer12"]["n_violations"] < 1:
        raise SystemExit("overlapping 12-mers in a run of 13 A were not counted")
    original = "ATGAAATAA"
    mutant = "ATGTAATAA"
    same = score_sequence(
        original, input_cds=original, weights=weights, index=index, constraints=constraints
    )
    changed = score_sequence(
        mutant, input_cds=original, weights=weights, index=index, constraints=constraints
    )
    if not same["translation_match"]:
        raise SystemExit("an unchanged CDS did not translate back to itself")
    if changed["translation_match"]:
        raise SystemExit("ATGAAATAA versus ATGTAATAA was accepted as the same protein")
    if changed["feasible"]:
        raise SystemExit("a protein change was marked feasible")
    print("self-test passed", file=sys.stderr)


def score_document(
    document: dict[str, Any], weights: dict[str, float], constraints: dict[str, Any]
) -> dict[str, Any]:
    if sha256_file(WEIGHTS_PATH) != document["meta"]["weights_sha256"]:
        raise SystemExit("weights file sha256 does not match the measurement")
    if sha256_file(CONSTRAINTS_PATH) != document["meta"]["constraints_sha256"]:
        raise SystemExit("constraint file sha256 does not match the measurement")
    index = _index(weights)
    inputs = {gene["gene"]: gene["cds"] for gene in document["inputs"]}
    for row in document["rows"]:
        scored = score_sequence(
            row.get("sequence"),
            input_cds=inputs[row["gene"]],
            weights=weights,
            index=index,
            constraints=constraints,
        )
        row.update(scored)
    document["meta"]["versions"]["biopython_scorer"] = BIOPYTHON_VERSION
    document["meta"]["versions"]["python_scorer"] = sys.version.split()[0]
    document["meta"]["scoring"] = {
        "cai": (
            "CAI is `Bio.SeqUtils.CodonAdaptationIndex.calculate` from "
            f"Biopython {BIOPYTHON_VERSION}. The reference weights are the "
            "relative-adaptiveness values in `examples/mg1655_ribosomal_usage.json` "
            f"(sha256 {document['meta']['weights_sha256']}). `calculate` skips "
            "ATG and TGG. TAA, TAG, and TGA are omitted from the index, so "
            "they are skipped as well. A counted codon of weight 0 is reported "
            "as 0, because the geometric mean of a set that contains 0 is 0 "
            "and `calculate` raises on `log(0)`. The deprecated "
            "`Bio.SeqUtils.CodonUsage` class, which divides by `n - 1`, is "
            "not used. `SharpEcoliIndex` is not used."
        ),
        "gc": (
            "Global GC is `Bio.SeqUtils.gc_fraction` from the same Biopython. "
            "Each hard GC window is the same function on every 50 nt slice, "
            "step 1, and the window fails when that value is below 0.40 or "
            "above 0.60. The endpoints are included."
        ),
        "constraints": (
            "Each hard constraint is counted by a string scan in "
            "`independent_score.py`: overlapping literal matches for sites "
            "and motifs, the reverse complement as well when the rule says "
            "both strands (a palindrome is counted once), one violation per "
            "homopolymer run longer than 6, and one violation per distinct "
            "12-mer whose overlapping starts exceed 1. A kind with no scan "
            "raises. The product's constraint counter is stored beside these "
            "counts and is not the headline."
        ),
        "translation": (
            "Each CDS is translated with `Bio.Seq.Seq.translate(table=11)`, "
            "including the terminal stop as `*`. The design matches when that "
            "string equals the translation of the input CDS. A different stop "
            "codon still matches. An amino-acid change does not."
        ),
    }
    document["comparison"] = compare_rows(document)
    return document


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--score", type=Path, help="unscored measurement JSON")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--doc", type=Path)
    parser.add_argument("--doc-only", type=Path, help="rewrite the markdown from a scored JSON")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.doc_only:
        document = json.loads(args.doc_only.read_text(encoding="utf-8"))
        if args.doc is None:
            raise SystemExit("--doc-only needs --doc")
        args.doc.write_text(render_report(document), encoding="utf-8")
        return
    if args.score is None or args.out is None or args.doc is None:
        raise SystemExit("need --score, --out, and --doc, or --self-test")
    document = json.loads(args.score.read_text(encoding="utf-8"))
    if document.get("meta", {}).get("partial"):
        raise SystemExit("refusing to score a partial gene list into the report")
    scored = score_document(
        document, load_weights(WEIGHTS_PATH), load_constraints(CONSTRAINTS_PATH)
    )
    args.out.write_text(json.dumps(scored, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.doc.write_text(render_report(scored), encoding="utf-8")


if __name__ == "__main__":
    main()
