#!/usr/bin/env python3
"""CodonForge against DNA Chisel, cai_greedy, and LinearDesign if the binary runs.

Part A of the v2.9.3 design-quality study. Numbers come from this process.
A tool that does not run is recorded as a failure, not as a score.

From the repository root, after ``make build`` and ``pip install dnachisel==3.2.16``:

    .venv/bin/python tools/external_comparison.py \\
        --out docs/benchmarks/v2.9.3-external-rows.json

DNA Chisel is invoked in a child process (``--chisel-worker``) so a stall can
be killed. LinearDesign is the binary compiled from
``LinearDesignSoftware/LinearDesign``; this script does not vendor it.
``LINEAR_DESIGN_BIN`` overrides the default path. CDSfold is not invoked:
the clone's Makefile needs ViennaRNA 2.1.1 headers and ``libRNA``, and this
environment's ViennaRNA 2.7.2 wheel ships neither.

The shared CAI is the independent Sharp and Li mean on
``examples/mg1655_ribosomal_usage.json``. DNA Chisel is given that same table
through ``CodonOptimize(..., method="use_best_codon")``, which constructs
``MaximizeCAI``. Its score is ``N * log(CAI)`` with a floor of 0.001 on a
zero frequency. The smallest positive ribosomal weight is above that floor,
so the codon ranking matches the sum of log weights. Met, Trp and the stop
still enter DNA Chisel's sum and do not enter ours. The reported CAI column
is always ours.

LinearDesign's objective is folding free energy plus ``lambda`` times a codon
log-weight. ``lambda`` 0 drops the codon term. Neither setting is
``usage.cai``, and the tool does not take the assembly, BioBrick or Golden
Gate files. Its rows are a separate table.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import subprocess
import sys
import time
import traceback
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import benchmark_common as common  # noqa: E402
from codonforge.bench.baselines import cai_greedy  # noqa: E402
from codonforge.modelreg import ArtifactStore  # noqa: E402

LINEAR_DESIGN_BIN = Path(
    os.environ.get(
        "LINEAR_DESIGN_BIN",
        "/tmp/lineardesign-try/LinearDesign/bin/LinearDesign_2D",
    )
)
LINEAR_DESIGN_ROOT = LINEAR_DESIGN_BIN.resolve().parents[1]
CHISEL_TIMEOUT_S = 120
LINEAR_TIMEOUT_S = 120

# DNA Chisel floors a zero frequency at 0.001 inside MaximizeCAI.evaluate.
DNA_CHISEL_ZERO_FLOOR = 0.001


def chisel_spec(constraints: dict[str, Any]) -> dict[str, Any]:
    """Map one CodonForge constraint set onto DNA Chisel specifications.

    A rule that does not have the same accept/reject language is listed in
    ``not_comparable`` and is not sent. The assembly file's rules do map;
    the function still checks each one instead of assuming that.
    """
    specs: list[dict[str, Any]] = [{"kind": "translation"}]
    comparable: list[str] = []
    notes: list[str] = []
    for constraint in constraints.get("constraints") or []:
        rule = constraint["rule"]
        kind = rule["kind"]
        if kind == "gc_window":
            window_matches = (
                int(rule["window_nt"]) == 50
                and float(rule["min"]) == 0.4
                and float(rule["max"]) == 0.6
            )
            if not window_matches:
                notes.append(f"{constraint['id']}: gc window is not the 50 nt 0.4-0.6 band")
            specs.append(
                {
                    "kind": "gc",
                    "mini": float(rule["min"]),
                    "maxi": float(rule["max"]),
                    "window": int(rule["window_nt"]),
                }
            )
            comparable.append(constraint["id"])
            notes.append(
                f"{constraint['id']}: EnforceGCContent window steps by 1, "
                "the same starts the CodonForge checker uses"
            )
        elif kind == "restriction_site":
            site = rule["site"].upper()
            specs.append({"kind": "avoid", "pattern": site, "strand": "both"})
            comparable.append(constraint["id"])
            notes.append(
                f"{constraint['id']}: AvoidPattern({site}) with strand both, "
                "which also rejects the reverse complement. Biopython's "
                f"{rule.get('enzyme')} site matches this string."
            )
        elif kind == "homopolymer":
            run = int(rule["max_run"]) + 1
            bases = "ACGT" if rule.get("base", "any") == "any" else rule["base"].upper()
            for base in bases:
                specs.append({"kind": "homopolymer", "base": base, "n": run})
            comparable.append(constraint["id"])
            notes.append(
                f"{constraint['id']}: HomopolymerPattern length {run} rejects a run "
                f"longer than {rule['max_run']}"
            )
        elif kind == "forbidden_motif":
            motif = rule["motif"].upper()
            specs.append({"kind": "avoid", "pattern": motif, "strand": "both"})
            comparable.append(constraint["id"])
        elif kind == "kmer_repeat":
            if int(rule["max_occurrences"]) != 1 or not rule.get("count_overlapping", False):
                notes.append(f"{constraint['id']}: k-mer rule is not max 1 overlapping")
                continue
            include_rc = bool(rule.get("canonical_strand", False))
            specs.append(
                {
                    "kind": "kmers",
                    "k": int(rule["k"]),
                    "include_reverse_complement": include_rc,
                }
            )
            comparable.append(constraint["id"])
            notes.append(
                f"{constraint['id']}: UniquifyAllKmers(k={rule['k']}, "
                f"include_reverse_complement={include_rc}). "
                "canonical_strand false means a k-mer and its reverse complement "
                "are different, which is include_reverse_complement False. "
                "DNA Chisel's own docstring says the constraint may fail to "
                "resolve a k-mer that already occurs more than twice. "
                "Occurrences are counted at every start, including overlaps."
            )
        else:
            notes.append(f"{constraint['id']}: kind {kind} was not mapped")
    return {"specs": specs, "comparable_ids": comparable, "notes": notes}


def chisel_worker() -> None:
    payload = json.loads(sys.stdin.read())
    try:
        import numpy as np
        from dnachisel import (
            AvoidPattern,
            CodonOptimize,
            DnaOptimizationProblem,
            EnforceGCContent,
            EnforceTranslation,
            HomopolymerPattern,
            NoSolutionError,
            UniquifyAllKmers,
        )

        random.seed(int(payload["seed"]))
        np.random.seed(int(payload["seed"]))
        constraints = []
        if payload.get("enforce_translation", True):
            constraints.append(EnforceTranslation())
        for spec in payload["specs"]:
            kind = spec["kind"]
            if kind == "translation":
                continue
            if kind == "gc":
                constraints.append(
                    EnforceGCContent(mini=spec["mini"], maxi=spec["maxi"], window=spec["window"])
                )
            elif kind == "avoid":
                constraints.append(AvoidPattern(spec["pattern"]))
            elif kind == "homopolymer":
                constraints.append(AvoidPattern(HomopolymerPattern(spec["base"], spec["n"])))
            elif kind == "kmers":
                constraints.append(
                    UniquifyAllKmers(
                        spec["k"],
                        include_reverse_complement=spec["include_reverse_complement"],
                    )
                )
            else:
                raise ValueError(f"unknown spec {kind}")
        problem = DnaOptimizationProblem(
            sequence=payload["sequence"],
            constraints=constraints,
            objectives=[
                CodonOptimize(
                    method="use_best_codon",
                    codon_usage_table=payload["codon_usage"],
                )
            ],
            logger=None,
        )
        if payload.get("mode") == "evaluate":
            scores = []
            for constraint in problem.constraints:
                evaluation = constraint.evaluate(problem)
                scores.append(
                    {
                        "label": constraint.label(),
                        "score": float(evaluation.score),
                        "passes": bool(evaluation.passes),
                    }
                )
            print(
                json.dumps(
                    {
                        "ok": True,
                        "sequence": str(problem.sequence),
                        "all_pass": bool(problem.all_constraints_pass()),
                        "scores": scores,
                        "mode": "evaluate",
                    }
                )
            )
            return
        try:
            problem.resolve_constraints()
            problem.optimize()
        except NoSolutionError as exc:
            print(
                json.dumps({"ok": False, "error_type": "NoSolutionError", "error": str(exc)[:2000]})
            )
            return
        scores = []
        for constraint in problem.constraints:
            evaluation = constraint.evaluate(problem)
            scores.append(
                {
                    "label": constraint.label(),
                    "score": float(evaluation.score),
                    "passes": bool(evaluation.passes),
                }
            )
        objective_score = None
        if problem.objectives:
            objective_score = float(problem.objectives[0].evaluate(problem).score)
        print(
            json.dumps(
                {
                    "ok": True,
                    "sequence": str(problem.sequence),
                    "all_pass": bool(problem.all_constraints_pass()),
                    "scores": scores,
                    "objective_score": objective_score,
                }
            )
        )
    except Exception as exc:
        print(
            json.dumps(
                {
                    "ok": False,
                    "error_type": type(exc).__name__,
                    "error": str(exc)[:2000],
                    "trace": traceback.format_exc()[-2000:],
                }
            )
        )


def run_chisel(
    seq: str,
    specs: list[dict[str, Any]],
    codon_usage: dict[str, dict[str, float]],
    seed: int,
    timeout_s: int,
    *,
    mode: str = "optimize",
) -> dict[str, Any]:
    request = json.dumps(
        {
            "sequence": seq,
            "specs": specs,
            "codon_usage": codon_usage,
            "seed": seed,
            "mode": mode,
        }
    )
    started = time.perf_counter()
    try:
        proc = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--chisel-worker"],
            input=request,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            cwd=ROOT,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return {
            "ok": False,
            "error_type": "TimeoutExpired",
            "wall_s": time.perf_counter() - started,
        }
    wall_s = time.perf_counter() - started
    line = proc.stdout.strip().splitlines()
    if not line:
        return {
            "ok": False,
            "error_type": "EmptyStdout",
            "error": (proc.stderr or "")[-2000:],
            "returncode": proc.returncode,
            "wall_s": wall_s,
        }
    parsed = json.loads(line[-1])
    parsed["wall_s"] = wall_s
    parsed["returncode"] = proc.returncode
    return parsed


def write_lineardesign_table(path: Path, weights: dict[str, float]) -> dict[str, Any]:
    """Frequencies proportional to relative adaptiveness, per amino acid.

    A zero weight is written as 1e-12 before normalisation so ``log2(w)`` in
    LinearDesign's CAI is finite. That is their table, not ours. The study
    does not report their printed CAI as Sharp and Li on the ribosomal weights.
    """
    import cai_independent

    floor = 1e-12
    floored_codons = []
    by_aa: dict[str, list[tuple[str, float]]] = {}
    for codon, amino in cai_independent.AA.items():
        weight = float(weights[codon])
        used = weight if weight > 0 else floor
        by_aa.setdefault(amino, []).append((codon, used))
        if weight == 0:
            floored_codons.append(codon)
    lines = ["#codon,aa,freq"]
    for amino in sorted(by_aa):
        total = sum(weight for _, weight in by_aa[amino])
        for codon, weight in sorted(by_aa[amino]):
            rna = codon.replace("T", "U")
            lines.append(f"{rna},{amino},{weight / total:.12g}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {"zero_weight_codons_floored_at": floor, "codons": floored_codons}


def run_lineardesign(
    protein: str, table: Path, lambda_value: float, timeout_s: int
) -> dict[str, Any]:
    if not LINEAR_DESIGN_BIN.is_file():
        return {"ok": False, "error_type": "MissingBinary", "path": str(LINEAR_DESIGN_BIN)}
    started = time.perf_counter()
    try:
        proc = subprocess.run(
            [str(LINEAR_DESIGN_BIN), str(lambda_value), "0", str(table)],
            input=protein + "\n",
            capture_output=True,
            text=True,
            timeout=timeout_s,
            cwd=LINEAR_DESIGN_ROOT,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return {
            "ok": False,
            "error_type": "TimeoutExpired",
            "wall_s": time.perf_counter() - started,
        }
    wall_s = time.perf_counter() - started
    sequence = None
    reported_mfe = None
    reported_cai = None
    reported_runtime = None
    for raw in proc.stdout.replace("\r", "\n").splitlines():
        line = raw.strip()
        if line.startswith("mRNA sequence:"):
            sequence = line.split(":", 1)[1].strip().replace(" ", "").replace("U", "T")
        elif "folding free energy:" in line and "CAI:" in line:
            # "mRNA folding free energy: -230.50 kcal/mol; mRNA CAI: 0.963"
            left, right = line.split(";")
            reported_mfe = float(left.split(":")[1].split()[0])
            reported_cai = float(right.split(":")[1].strip())
        elif line.startswith("Runtime:"):
            reported_runtime = float(line.split()[1])
    if proc.returncode != 0 or not sequence:
        return {
            "ok": False,
            "error_type": "NonZero" if proc.returncode != 0 else "Unparsed",
            "returncode": proc.returncode,
            "stderr": (proc.stderr or "")[-1500:],
            "stdout_tail": proc.stdout[-1500:],
            "wall_s": wall_s,
        }
    return {
        "ok": True,
        "sequence": sequence,
        "reported_mfe_kcal_mol": reported_mfe,
        "reported_cai_their_definition": reported_cai,
        "reported_runtime_s": reported_runtime,
        "wall_s": wall_s,
        "returncode": proc.returncode,
    }


def score_row(
    seq: str | None,
    source: str,
    weights: dict[str, float],
    profile: dict[str, Any],
    digest: str,
    constraints: dict[str, Any] | None,
) -> dict[str, Any]:
    if not seq:
        return {}
    import cai_independent

    ours = common.independent_cai(seq, weights)
    engine = common.engine_cai(seq, profile, digest)
    delta = None
    if ours["cai"] is not None and engine is not None:
        delta = abs(float(ours["cai"]) - engine)
    produced_protein = cai_independent.translate(seq)
    source_protein = cai_independent.translate(source)
    row = {
        "cds": seq,
        "sha256": common.sha256_text(seq),
        "nt": len(seq),
        "gc": common.gc_fraction(seq),
        "independent_cai": ours,
        "usage_cai": engine,
        "cai_abs_delta": delta,
        "translation_equal_length": len(seq) == len(source),
    }
    if produced_protein == source_protein:
        row["protein_relation"] = "identical_including_stop"
    elif source_protein.endswith("*") and produced_protein == source_protein[:-1]:
        row["protein_relation"] = "matches_without_terminal_stop"
    elif produced_protein.endswith("*") and produced_protein[:-1] == source_protein.rstrip("*"):
        row["protein_relation"] = "same_sense_with_a_stop"
    else:
        row["protein_relation"] = "different_protein"
    if constraints is not None and len(seq) % 3 == 0:
        row["our_checker"] = common.violation_summary(seq, constraints)
    return row


def mapping_checks(
    weights: dict[str, float], codon_usage: dict[str, dict[str, float]]
) -> dict[str, Any]:
    """Fixtures that must fail both checkers, or the mapping is not the one described."""
    assembly = json.loads(common.ASSEMBLY_PATH.read_text(encoding="utf-8"))
    mapped = chisel_spec(assembly)
    # 54 nt so a 50 nt window exists. Evaluate the input; do not let Chisel repair it.
    cases = {
        "gc_window": "ATG" + "GGC" * 16 + "TAA",
        "bsai_reverse": "ATGGAGACC" + "AAC" * 14 + "TAA",
        "homopolymer_a": "ATG" + "AAA" * 4 + "GCT" * 12 + "TAA",
    }
    checks = []
    for name, seq in cases.items():
        ours = common.violation_summary(seq, assembly)
        chisel = run_chisel(
            seq, mapped["specs"], codon_usage, seed=7, timeout_s=60, mode="evaluate"
        )
        checks.append(
            {
                "name": name,
                "nt": len(seq),
                "our_feasible": ours["feasible"],
                "our_ids": ours["actionable_by_id"],
                "chisel_ok": chisel.get("ok"),
                "chisel_error_type": chisel.get("error_type"),
                "chisel_all_pass": chisel.get("all_pass"),
                "chisel_failing": [
                    item["label"]
                    for item in chisel.get("scores") or []
                    if item.get("passes") is False
                ],
            }
        )
    return {
        "mode": "evaluate_without_resolve",
        "cases": checks,
        "negative_cases_fail_both": all(
            item["our_feasible"] is False and item["chisel_all_pass"] is False for item in checks
        ),
    }


def finish_tool_row(
    gene: dict[str, Any],
    tool: str,
    constraint_name: str | None,
    raw: dict[str, Any],
    weights: dict[str, float],
    profile: dict[str, Any],
    digest: str,
    constraints: dict[str, Any] | None,
    objective_note: str,
) -> dict[str, Any]:
    seq = raw.get("sequence") or raw.get("cds")
    scored = score_row(seq, gene["cds"], weights, profile, digest, constraints) if seq else {}
    access = None
    if seq and scored.get("protein_relation") != "different_protein":
        try:
            access = common.accessibility(seq)
        except Exception as exc:  # a fold failure is data, not a reason to drop the row
            access = {"error": type(exc).__name__, "detail": str(exc)[:500]}
    return {
        "gene": gene["gene"],
        "nt": gene["nt"],
        "codons": gene["codons"],
        "tool": tool,
        "constraint_set": constraint_name,
        "objective_note": objective_note,
        "wall_s": raw.get("wall_s"),
        "tool_ok": raw.get("ok", raw.get("state") == "succeeded"),
        "tool_error_type": raw.get("error_type"),
        "tool_error": raw.get("error"),
        "our_checker": scored.get("our_checker"),
        "independent_cai": (scored.get("independent_cai") or {}).get("cai"),
        "independent_detail": scored.get("independent_cai"),
        "usage_cai": scored.get("usage_cai"),
        "cai_abs_delta": scored.get("cai_abs_delta"),
        "gc": scored.get("gc"),
        "protein_relation": scored.get("protein_relation"),
        "sha256": scored.get("sha256"),
        "accessibility_penalty": access,
        "cds": scored.get("cds"),
        "extra": {
            key: raw[key]
            for key in (
                "all_pass",
                "scores",
                "objective_score",
                "state",
                "warnings",
                "reported_mfe_kcal_mol",
                "reported_cai_their_definition",
                "reported_runtime_s",
                "lambda",
                "returncode",
                "stderr",
            )
            if key in raw
        },
    }


def cross_check_saved(path: Path) -> None:
    """Evaluate saved sequences with DNA Chisel. Does not re-optimise.

    The optimiser rows stay as they were written. This adds ``chisel_evaluate``
    on CodonForge and DNA Chisel rows so both checkers are on the same CDS.
    """
    document = json.loads(path.read_text(encoding="utf-8"))
    weights = common.load_weights()
    codon_usage = common.codon_usage_by_amino_acid(weights)
    specs = {item["name"]: chisel_spec(item["document"]) for item in common.load_constraint_sets()}
    evaluated = 0
    for row in document["rows"]:
        if row.get("tool") not in {"codonforge", "dnachisel"}:
            continue
        sequence = row.get("cds")
        constraint_name = row.get("constraint_set")
        if not sequence or constraint_name not in specs:
            continue
        print(f"evaluate {row['tool']} {row['gene']} {constraint_name}", flush=True)
        result = run_chisel(
            sequence,
            specs[constraint_name]["specs"],
            codon_usage,
            seed=7,
            timeout_s=60,
            mode="evaluate",
        )
        row["chisel_evaluate"] = {
            "ok": result.get("ok"),
            "all_pass": result.get("all_pass"),
            "error_type": result.get("error_type"),
            "failing": [
                item["label"] for item in result.get("scores") or [] if item.get("passes") is False
            ],
        }
        evaluated += 1
    document["chisel_evaluate_note"] = (
        "Follow-up evaluation of the sequences already in rows. "
        "DNA Chisel mode=evaluate, no resolve_constraints and no optimize. "
        "Same specification list as the optimisation runs."
    )
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"cross-check wrote {path} evaluations {evaluated}", flush=True)


def main() -> None:
    if "--chisel-worker" in sys.argv:
        chisel_worker()
        return
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--skip-lineardesign", action="store_true")
    parser.add_argument(
        "--cross-check-file",
        type=Path,
        help="Evaluate sequences already stored in this JSON. Does not re-optimise.",
    )
    parser.add_argument(
        "--genes",
        default="",
        help="Comma-separated gene names. Default: the full shared set.",
    )
    args = parser.parse_args()
    if args.cross_check_file is not None:
        cross_check_saved(args.cross_check_file)
        return

    weights = common.load_weights()
    profile = common.host_profile(weights)
    digest = common.profile_digest(profile)
    genes = common.load_genes()
    if args.genes:
        wanted = {item.strip() for item in args.genes.split(",") if item.strip()}
        genes = [gene for gene in genes if gene["gene"] in wanted]
    sets = common.load_constraint_sets()
    codon_usage = common.codon_usage_by_amino_acid(weights)
    minima = common.family_maxima(weights)
    floor_ok = common.minimum_positive_weight(weights) > DNA_CHISEL_ZERO_FLOOR
    store = ArtifactStore(Path("/tmp/cf-v293-external-store"))
    rows: list[dict[str, Any]] = []
    print("mapping checks", flush=True)
    checks = mapping_checks(weights, codon_usage)

    objective_chisel = (
        "DNA Chisel MaximizeCAI via CodonOptimize(method='use_best_codon') on the "
        "ribosomal relative-adaptiveness table, zeros left as zeros. "
        f"Zero floor inside DNA Chisel is {DNA_CHISEL_ZERO_FLOOR}; "
        f"minimum positive ribosomal weight is {common.minimum_positive_weight(weights)}. "
        "Ranking matches sum of log weights only while that minimum stays above the floor: "
        f"{floor_ok}. Reported CAI is the independent Sharp and Li mean, which excludes "
        "Met, Trp and stops. DNA Chisel's own objective score is N*log on every codon, "
        "with the floor. Seed is numpy and random seed 7, one run, unless a probe says otherwise."
    )
    objective_cf = (
        "CodonForge usage.cai, budget standard, generations 4, stagnation window 21, "
        "seed 7, host escherichia_coli:511145/MG1655, translation table 11. "
        "The generation cap is what stops tier 3; the window is longer than 4."
    )
    objective_greedy = (
        "cai_greedy picks the highest relative adaptiveness at each codon and breaks "
        "ties toward the lexicographically smaller codon. It does not see constraints. "
        "Under a constraint set its CAI is a floor only when the sequence is feasible; "
        "an infeasible CAI of 1 is not a feasible design."
    )

    for gene in genes:
        started = time.perf_counter()
        greedy = cai_greedy(gene["cds"], common.TRANSLATION_TABLE, weights)
        greedy_wall = time.perf_counter() - started
        for constraint in sets:
            print(f"codonforge {gene['gene']} {constraint['name']}", flush=True)
            designed = common.run_codonforge(
                gene["cds"],
                constraint["document"],
                profile,
                digest,
                store,
                budget="standard",
                seed=7,
                generations=4,
                stagnation_window=21,
                job_id=f"ext-{gene['gene']}-{constraint['name']}",
            )
            designed["sequence"] = designed.get("cds")
            designed["ok"] = designed["state"] == "succeeded" and bool(designed.get("cds"))
            rows.append(
                finish_tool_row(
                    gene,
                    "codonforge",
                    constraint["name"],
                    designed,
                    weights,
                    profile,
                    digest,
                    constraint["document"],
                    objective_cf,
                )
            )
            print(f"chisel {gene['gene']} {constraint['name']}", flush=True)
            mapped = chisel_spec(constraint["document"])
            chisel = run_chisel(
                gene["cds"], mapped["specs"], codon_usage, seed=7, timeout_s=CHISEL_TIMEOUT_S
            )
            row = finish_tool_row(
                gene,
                "dnachisel",
                constraint["name"],
                chisel,
                weights,
                profile,
                digest,
                constraint["document"],
                objective_chisel,
            )
            row["mapping"] = {
                "comparable_ids": mapped["comparable_ids"],
                "notes": mapped["notes"],
            }
            row["extra"]["chisel_seed"] = 7
            rows.append(row)
            greedy_raw = {
                "ok": True,
                "sequence": greedy,
                "wall_s": greedy_wall,
                "state": "succeeded",
            }
            rows.append(
                finish_tool_row(
                    gene,
                    "cai_greedy",
                    constraint["name"],
                    greedy_raw,
                    weights,
                    profile,
                    digest,
                    constraint["document"],
                    objective_greedy,
                )
            )
            native_raw = {"ok": True, "sequence": gene["cds"], "wall_s": 0.0, "state": "input"}
            rows.append(
                finish_tool_row(
                    gene,
                    "native_input",
                    constraint["name"],
                    native_raw,
                    weights,
                    profile,
                    digest,
                    constraint["document"],
                    "The input CDS. Not an optimiser.",
                )
            )

    # One DNA Chisel seed probe. The optimiser is stochastic; one seed is not the tool.
    probe_gene = next(gene for gene in genes if gene["gene"] == "rpsL")
    assembly = next(item for item in sets if item["name"] == "assembly")
    mapped = chisel_spec(assembly["document"])
    chisel_seeds = []
    for seed in (1, 7, 11):
        print(f"chisel seed probe rpsL assembly seed {seed}", flush=True)
        chisel = run_chisel(
            probe_gene["cds"], mapped["specs"], codon_usage, seed=seed, timeout_s=CHISEL_TIMEOUT_S
        )
        row = finish_tool_row(
            probe_gene,
            "dnachisel_seed_probe",
            "assembly",
            chisel,
            weights,
            profile,
            digest,
            assembly["document"],
            objective_chisel + f" Probe seed {seed}.",
        )
        row["extra"]["chisel_seed"] = seed
        chisel_seeds.append(row)
        rows.append(row)

    lineardesign_rows: list[dict[str, Any]] = []
    lineardesign_meta: dict[str, Any] = {
        "binary": str(LINEAR_DESIGN_BIN),
        "present": LINEAR_DESIGN_BIN.is_file(),
        "commit": "f0126ca89a8b853088b4bccfd2cc8c378d3678be",
        "source": "https://github.com/LinearDesignSoftware/LinearDesign",
        "compile": "make, g++ -std=c++11, linked to src/Utils/libraries/LinearDesign_linux64.so",
        "objective": (
            "CLI lambda is multiplied by 100 in src/linear_design.cpp and used as the "
            "coefficient on log codon weight in src/Utils/network.h. Lambda 0 removes "
            "that term. Folding free energy remains in the search at every lambda. "
            "This is not usage.cai and it is not given the constraint files. "
            "Their printed CAI is log2 geometric mean over every codon in the mRNA, "
            "including Met and Trp, on the frequency table written for this run."
        ),
        "comparable": False,
    }
    if not args.skip_lineardesign and LINEAR_DESIGN_BIN.is_file():
        table_path = Path("/tmp/cf-v293-lineardesign-codon-table.csv")
        lineardesign_meta["table"] = write_lineardesign_table(table_path, weights)
        for gene in genes:
            protein = gene["protein"]
            if protein.endswith("*"):
                protein = protein[:-1]
            if "*" in protein:
                lineardesign_rows.append(
                    {"gene": gene["gene"], "ok": False, "error_type": "InternalStop"}
                )
                continue
            for lambda_value in (0.0, 4.0):
                print(f"lineardesign {gene['gene']} lambda {lambda_value}", flush=True)
                raw = run_lineardesign(protein, table_path, lambda_value, LINEAR_TIMEOUT_S)
                raw["lambda"] = lambda_value
                row = finish_tool_row(
                    gene,
                    "lineardesign",
                    None,
                    raw,
                    weights,
                    profile,
                    digest,
                    None,
                    lineardesign_meta["objective"],
                )
                # Feasibility against each of our sets is an evaluation, not a request.
                if raw.get("sequence"):
                    row["our_checker_by_set"] = {
                        item["name"]: common.violation_summary(raw["sequence"], item["document"])
                        for item in sets
                    }
                row["extra"]["lambda"] = lambda_value
                lineardesign_rows.append(row)
                rows.append(row)
    elif not LINEAR_DESIGN_BIN.is_file():
        lineardesign_meta["error"] = f"{LINEAR_DESIGN_BIN} is not a file"

    try:
        import dnachisel

        chisel_version = dnachisel.__version__
    except ImportError:
        chisel_version = None

    document = {
        "study": "v2.9.3 external comparison",
        "versions": common.tool_versions(),
        "dnachisel": chisel_version,
        "biopython": __import__("Bio").__version__,
        "numpy": __import__("numpy").__version__,
        "weights": str(common.WEIGHTS_PATH.relative_to(common.ROOT)),
        "host": common.HOST,
        "translation_table": common.TRANSLATION_TABLE,
        "family_maxima_all_one": all(value == 1.0 for value in minima.values()),
        "dnachisel_zero_floor": DNA_CHISEL_ZERO_FLOOR,
        "minimum_positive_weight": common.minimum_positive_weight(weights),
        "zero_floor_does_not_invert_positive_codons": floor_ok,
        "commands": {
            "this_script": (
                ".venv/bin/python tools/external_comparison.py "
                "--out docs/benchmarks/v2.9.3-external-rows.json"
            ),
            "chisel_install": ".venv/bin/pip install 'dnachisel==3.2.16'",
            "lineardesign_compile": (
                "git clone https://github.com/LinearDesignSoftware/LinearDesign.git && make"
            ),
            "lineardesign_run": (
                "cd LinearDesign && ./bin/LinearDesign_2D <lambda> 0 <codon_table.csv>  "
                "with the amino-acid sequence on stdin, terminal stop removed"
            ),
        },
        "cdsfold": {
            "installed": False,
            "commit_cloned": "06f3ee8e227f39edb6ddb639082cc070c7df081a",
            "source": "https://github.com/gterai/CDSfold",
            "reason": (
                "Makefile sets VIENNA and compiles with -I$VIENNA/include/ViennaRNA "
                "and -lRNA. README.rd asks for ViennaRNA 2.1.1. This environment has "
                "ViennaRNA 2.7.2 as a Python wheel (RNA/__init__.py and the extension "
                "module only). There is no compile-time header and no libRNA to link. "
                "CDSfold was not patched onto 2.7.2. No CDSfold number is in this file."
            ),
        },
        "mapping_checks": checks,
        "lineardesign": lineardesign_meta,
        "chisel_seed_probe_sha256": [row.get("sha256") for row in chisel_seeds],
        "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.out} rows {len(rows)}", flush=True)


if __name__ == "__main__":
    main()
