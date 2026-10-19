#!/usr/bin/env python3
"""Head-to-head of DNA Design Suite and DNA Chisel on the GATE-BEN-001 proteins.

Both tools get the same 50 CDS, the same ribosomal relative-adaptiveness
table, and the same hard constraints
(``benchmarks/validation/constraints_assembly.json``). DNA Design Suite uses
the GATE-BEN-001 job settings (budget ``fast``, 4 generations, seed 7,
segment 2, objective ``usage.cai``). The published gate itself uses an empty
constraint list; this run does not, and it does not replace that claim.

Headline CAI, GC, translation, and constraint counts are computed by
``independent_score.py`` under Biopython 1.81. This process does not score
those headlines.

From the repository root, after ``make build`` and
``pip install 'dnachisel==3.2.16'``, and after the scoring environment in
the module docstring of ``independent_score.py``:

    .venv/bin/python benchmarks/gate_ben_001_chisel/run_benchmark.py \\
        --out docs/benchmarks/gate-ben-001-dnachisel-rows.json \\
        --doc docs/benchmarks/gate-ben-001-dnachisel.md
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import random
import subprocess
import sys
import time
import traceback
from importlib.metadata import version as package_version
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
WEIGHTS_PATH = ROOT / "examples" / "mg1655_ribosomal_usage.json"
CONSTRAINTS_PATH = ROOT / "benchmarks" / "validation" / "constraints_assembly.json"
SCORER = Path(__file__).with_name("independent_score.py")
SCORE_PYTHON = Path(__file__).with_name(".score-venv") / "bin" / "python"
HOST = "escherichia_coli:511145/MG1655"
SEED = 7
BUDGET = "fast"
GENERATIONS = 4
SEGMENT = 2
STAGNATION_WINDOW = 8
CHISEL_TIMEOUT_S = 180
CHISEL_ZERO_FLOOR = 0.001

# Kinds in the constraint schema that this shared file does not contain.
# They are recorded so a reader can see they were not silently dropped.
SCHEMA_KINDS_NOT_IN_FILE = [
    {
        "kind": "local_structure",
        "reason": (
            "No DNA Chisel specification is this predicate. AvoidHairpins is a "
            "hairpin rule, not a windowed unpaired probability or a pinned "
            "minimum free energy, and it does not name the folding engine. "
            "The shared file does not contain this kind, so neither tool was given it."
        ),
    },
    {
        "kind": "codon_ban",
        "reason": (
            "DNA Chisel has no hard specification that removes a listed codon "
            "and keeps the amino acid. A zero in the codon table changes the "
            "CAI objective and still leaves the codon legal under translation. "
            "The shared file does not contain this kind, so neither tool was given it."
        ),
    },
    {
        "kind": "min_codon_frequency",
        "reason": (
            "AvoidRareCodons compares a frequency in DNA Chisel's codon table. "
            "The shared numbers are relative adaptiveness, not frequencies, so "
            "the same cutoff would be a different rule. The shared file does "
            "not contain this kind, so neither tool was given it."
        ),
    },
    {
        "kind": "custom_regex",
        "reason": (
            "DNA Chisel SequencePattern can carry a regular expression, so a "
            "translation of this kind exists. The shared file does not contain "
            "one, so neither tool was given one."
        ),
    },
]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_commit() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()


def machine() -> dict[str, Any]:
    cpu = "unknown"
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text(encoding="utf-8").splitlines():
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    mem_kib = None
    meminfo = Path("/proc/meminfo")
    if meminfo.is_file():
        for line in meminfo.read_text(encoding="utf-8").splitlines():
            if line.startswith("MemTotal:"):
                mem_kib = int(line.split()[1])
                break
    return {
        "platform": platform.platform(),
        "cpu": cpu,
        "nproc": os.cpu_count(),
        "mem_total_kib": mem_kib,
    }


def translate_constraints(document: dict[str, Any]) -> dict[str, Any]:
    """Map each shared rule onto a DNA Chisel specification, or list it as untranslated.

    A rule that is not the same accept/reject predicate is not sent under a
    weaker specification.
    """
    specs: list[dict[str, Any]] = []
    translated: list[dict[str, str]] = []
    not_translated: list[dict[str, str]] = []
    notes: list[str] = []
    for constraint in document["constraints"]:
        rule = constraint["rule"]
        kind = rule["kind"]
        identity = constraint["id"]
        if kind == "gc_window":
            specs.append(
                {
                    "kind": "gc",
                    "mini": float(rule["min"]),
                    "maxi": float(rule["max"]),
                    "window": int(rule["window_nt"]),
                }
            )
            translated.append(
                {
                    "id": identity,
                    "rule": (
                        f"gc_window {rule['window_nt']} nt, "
                        f"{rule['min']} to {rule['max']}, every start"
                    ),
                    "spec": (
                        f"EnforceGCContent(mini={rule['min']}, maxi={rule['max']}, "
                        f"window={rule['window_nt']}). The implementation returns "
                        "len(sequence) - window + 1 windows at step 1. Bounds are inclusive."
                    ),
                }
            )
        elif kind == "restriction_site":
            if rule.get("strand", "both") != "both" or int(rule.get("max_occurrences", 0)) != 0:
                not_translated.append(
                    {
                        "id": identity,
                        "reason": "only a both-strands ban (max_occurrences 0) is mapped",
                    }
                )
                continue
            site = rule["site"].upper()
            specs.append({"kind": "avoid", "pattern": site})
            translated.append(
                {
                    "id": identity,
                    "rule": f"restriction_site {rule.get('enzyme')} {site}, both strands, max 0",
                    "spec": (
                        f"AvoidPattern({site}). With no location the strand is 0, "
                        "and initialization copies that strand onto the full sequence, "
                        "so the reverse complement is rejected. A palindrome is searched once."
                    ),
                }
            )
        elif kind == "forbidden_motif":
            if rule.get("strand", "both") != "both" or int(rule.get("max_occurrences", 0)) != 0:
                not_translated.append(
                    {
                        "id": identity,
                        "reason": "only a both-strands ban (max_occurrences 0) is mapped",
                    }
                )
                continue
            motif = rule["motif"].upper()
            specs.append({"kind": "avoid", "pattern": motif})
            translated.append(
                {
                    "id": identity,
                    "rule": f"forbidden_motif {motif}, both strands, max 0",
                    "spec": (
                        f"AvoidPattern({motif}) on both strands, same path as a restriction site."
                    ),
                }
            )
        elif kind == "homopolymer":
            run = int(rule["max_run"]) + 1
            bases = "ACGT" if rule.get("base", "any") == "any" else rule["base"].upper()
            for base in bases:
                specs.append({"kind": "homopolymer", "base": base, "n": run})
            translated.append(
                {
                    "id": identity,
                    "rule": (
                        f"homopolymer base {rule.get('base', 'any')}, max run {rule['max_run']}"
                    ),
                    "spec": (
                        f"AvoidPattern(HomopolymerPattern(base, {run})) for each of {bases}. "
                        "A match of that length is a run longer than the limit. Strand 0 also "
                        "sees the reverse complement, which is the complementary base."
                    ),
                }
            )
        elif kind == "kmer_repeat":
            overlapping = bool(rule.get("count_overlapping", True))
            if int(rule["max_occurrences"]) != 1 or not overlapping:
                not_translated.append(
                    {
                        "id": identity,
                        "reason": (
                            "UniquifyAllKmers requires every k-mer at most once and "
                            "counts every start. A different limit or a non-overlapping "
                            "count is not that specification."
                        ),
                    }
                )
                continue
            include_rc = bool(rule.get("canonical_strand", False))
            specs.append(
                {
                    "kind": "kmers",
                    "k": int(rule["k"]),
                    "include_reverse_complement": include_rc,
                }
            )
            translated.append(
                {
                    "id": identity,
                    "rule": (
                        f"kmer_repeat k={rule['k']}, max 1, overlapping "
                        f"{overlapping}, canonical_strand {include_rc}"
                    ),
                    "spec": (
                        f"UniquifyAllKmers(k={rule['k']}, "
                        f"include_reverse_complement={include_rc}). "
                        "Occurrences are every start in the reference span. "
                        "canonical_strand false keeps a k-mer distinct from its reverse "
                        "complement. DNA Chisel's docstring says this constraint may fail "
                        "to resolve a k-mer that already occurs more than twice; that "
                        "outcome is recorded as no sequence."
                    ),
                }
            )
            notes.append(
                "UniquifyAllKmers can fail to clear a k-mer that occurs more than twice. "
                "That is a solver limit of a specification that was sent, not an untranslated rule."
            )
        else:
            not_translated.append(
                {
                    "id": identity,
                    "reason": (
                        f"kind {kind} has no DNA Chisel specification that is the same predicate"
                    ),
                }
            )
    return {
        "specs": specs,
        "translated": translated,
        "not_translated": not_translated,
        "notes": notes,
    }


def codon_usage_table(weights: dict[str, float]) -> dict[str, dict[str, float]]:
    """Relative adaptiveness grouped by amino acid, the shape DNA Chisel reads.

    The values are not frequencies and do not sum to 1. Family maxima are 1,
    so MaximizeCAI's ``log(w) - log(w_max)`` equals ``log(w)``.
    """
    from Bio.Data import CodonTable

    table = CodonTable.unambiguous_dna_by_name["Bacterial"]
    grouped: dict[str, dict[str, float]] = {}
    for codon, amino in table.forward_table.items():
        grouped.setdefault(amino, {})[codon] = float(weights[codon])
    grouped["*"] = {codon: float(weights[codon]) for codon in table.stop_codons}
    return grouped


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
        constraints = [
            EnforceTranslation(genetic_table="Bacterial", start_codon="ATG"),
        ]
        for spec in payload["specs"]:
            kind = spec["kind"]
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
        started = time.perf_counter()
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
        initialized = time.perf_counter()
        error_type = None
        error = None
        try:
            problem.resolve_constraints()
            resolved = time.perf_counter()
            problem.optimize()
            optimized = time.perf_counter()
        except NoSolutionError as exc:
            resolved = time.perf_counter()
            optimized = resolved
            error_type = type(exc).__name__
            error = str(exc)[:500]
        checks = []
        for constraint in problem.constraints:
            evaluation = constraint.evaluate(problem)
            checks.append(
                {
                    "name": type(constraint).__name__,
                    "passes": bool(evaluation.passes),
                    "score": float(evaluation.score),
                }
            )
        result = {
            "ok": error_type is None,
            "sequence": None if error_type else str(problem.sequence).upper(),
            "error_type": error_type,
            "error": error,
            "solver_s": optimized - started,
            "initialize_s": initialized - started,
            "resolve_s": resolved - initialized,
            "optimize_s": optimized - resolved,
            "all_constraints_pass": bool(problem.all_constraints_pass())
            if error_type is None
            else False,
            "constraint_checks": checks,
            "numpy": np.__version__,
        }
    except Exception as exc:
        traceback.print_exc(file=sys.stderr)
        result = {
            "ok": False,
            "sequence": None,
            "error_type": type(exc).__name__,
            "error": str(exc)[:500],
            "solver_s": None,
            "all_constraints_pass": False,
            "constraint_checks": [],
        }
    json.dump(result, sys.stdout)
    sys.stdout.write("\n")


def run_chisel(sequence: str, specs: list[dict[str, Any]], usage: dict[str, Any]) -> dict[str, Any]:
    payload = json.dumps(
        {
            "sequence": sequence,
            "seed": SEED,
            "specs": specs,
            "codon_usage": usage,
        }
    )
    env = os.environ.copy()
    env["PYTHONHASHSEED"] = str(SEED)
    started = time.perf_counter()
    try:
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--chisel-worker"],
            input=payload,
            text=True,
            capture_output=True,
            timeout=CHISEL_TIMEOUT_S,
            cwd=ROOT,
            env=env,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        return {
            "ok": False,
            "sequence": None,
            "error_type": "TimeoutExpired",
            "error": f"no result after {CHISEL_TIMEOUT_S} s",
            "solver_s": None,
            "wall_s": time.perf_counter() - started,
            "stderr": (exc.stderr or "")[-500:] if isinstance(exc.stderr, str) else "",
            "all_constraints_pass": False,
        }
    wall_s = time.perf_counter() - started
    try:
        result = json.loads(completed.stdout)
    except json.JSONDecodeError:
        result = {
            "ok": False,
            "sequence": None,
            "error_type": "WorkerOutput",
            "error": completed.stdout[-500:],
            "solver_s": None,
            "all_constraints_pass": False,
        }
    result["wall_s"] = wall_s
    result["returncode"] = completed.returncode
    if completed.returncode != 0 and result.get("ok"):
        result["ok"] = False
        result["error_type"] = result.get("error_type") or "NonZeroExit"
    if completed.stderr:
        result["stderr"] = completed.stderr[-500:]
    return result


def suite_job(
    cds: str,
    profile: dict[str, Any],
    constraints: dict[str, Any],
    digest: str,
    store_dir: Path,
    gene: str,
) -> dict[str, Any]:
    from dna_design_suite.codon.orchestrator.design import DesignRequest, design, score_sequences
    from dna_design_suite.modelreg import ArtifactStore

    store = ArtifactStore(store_dir / gene)
    request = DesignRequest(
        host=HOST,
        cds=cds,
        host_profile=profile,
        objectives=["usage.cai"],
        constraints=constraints,
        translation_table=11,
        budget=BUDGET,
        seed=SEED,
        generations=GENERATIONS,
        segment=SEGMENT,
        stagnation_window=STAGNATION_WINDOW,
        host_digest=digest,
    )
    started = time.perf_counter()
    try:
        result = design(request, store)
    except Exception as exc:
        return {
            "ok": False,
            "sequence": None,
            "state": None,
            "error_type": type(exc).__name__,
            "error": str(exc)[:500],
            "solver_s": time.perf_counter() - started,
            "wall_s": time.perf_counter() - started,
            "warnings": [],
            "diagnostic_usage_cai": None,
            "diagnostic_product_checker": None,
        }
    solver_s = time.perf_counter() - started
    sequence = None
    n_candidates = 0
    if result.front and result.front.get("candidates"):
        n_candidates = len(result.front["candidates"])
        sequence = result.front["candidates"][0]["sequences"][0]["cds"].upper()
    warnings = []
    for item in result.job.warnings:
        warnings.append(
            {
                "code": item.get("code"),
                "stage": item.get("stage"),
                "title": item.get("title"),
            }
        )
    usage_cai = None
    product = None
    if sequence:
        try:
            usage_cai = score_sequences(request, [sequence])[0].get("usage.cai")
        except Exception as exc:
            usage_cai = None
            warnings.append({"code": "diagnostic", "stage": "score", "title": str(exc)[:200]})
        try:
            from dna_design_suite import _core

            report = json.loads(
                _core.constraint_violations(
                    json.dumps(constraints),
                    sequence,
                    [(0, len(sequence))],
                    "assembled",
                )
            )
            counts: dict[str, int] = {}
            for item in report["violations"]:
                if item["actionable"]:
                    counts[item["constraint_id"]] = counts.get(item["constraint_id"], 0) + 1
            product = {
                "feasible": not counts,
                "actionable_by_id": dict(sorted(counts.items())),
                "n_actionable": sum(counts.values()),
            }
        except Exception as exc:
            product = {"feasible": None, "error": str(exc)[:200]}
    return {
        "ok": result.job.state == "succeeded" and sequence is not None,
        "sequence": sequence,
        "state": result.job.state,
        "n_candidates": n_candidates,
        "error_type": None if sequence else "NoSequence",
        "error": None if sequence else f"state {result.job.state}",
        "solver_s": solver_s,
        "wall_s": solver_s,
        "warnings": warnings,
        "diagnostic_usage_cai": usage_cai,
        "diagnostic_product_checker": product,
    }


def chisel_diagnostics(
    sequence: str | None, constraints: dict[str, Any], profile: dict[str, Any], digest: str
) -> dict[str, Any]:
    """Product counters only. The headline scorer does not read these."""
    if not sequence:
        return {"diagnostic_usage_cai": None, "diagnostic_product_checker": None}
    from dna_design_suite import _core
    from dna_design_suite.codon.orchestrator.design import DesignRequest, score_sequences

    request = DesignRequest(
        host=HOST,
        cds=sequence,
        host_profile=profile,
        objectives=["usage.cai"],
        constraints=constraints,
        translation_table=11,
        budget=BUDGET,
        seed=SEED,
        generations=0,
        segment=SEGMENT,
        stagnation_window=STAGNATION_WINDOW,
        host_digest=digest,
        job_id="diagnostic",
    )
    usage_cai = None
    try:
        usage_cai = score_sequences(request, [sequence])[0].get("usage.cai")
    except Exception:
        usage_cai = None
    try:
        report = json.loads(
            _core.constraint_violations(
                json.dumps(constraints),
                sequence,
                [(0, len(sequence))],
                "assembled",
            )
        )
        counts: dict[str, int] = {}
        for item in report["violations"]:
            if item["actionable"]:
                counts[item["constraint_id"]] = counts.get(item["constraint_id"], 0) + 1
        product = {
            "feasible": not counts,
            "actionable_by_id": dict(sorted(counts.items())),
            "n_actionable": sum(counts.values()),
        }
    except Exception as exc:
        product = {"feasible": None, "error": str(exc)[:200]}
    return {"diagnostic_usage_cai": usage_cai, "diagnostic_product_checker": product}


def profile_for(weights: dict[str, float]) -> tuple[dict[str, Any], str]:
    profile = {
        "codon_usage": {
            "subset": "user_supplied",
            "table": {
                codon: {"relative_adaptiveness": weights[codon]} for codon in sorted(weights)
            },
        }
    }
    payload = json.dumps(profile, sort_keys=True, separators=(",", ":")).encode()
    return profile, "sha256:" + hashlib.sha256(payload).hexdigest()


def write_fasta(path: Path, document: dict[str, Any]) -> None:
    chunks = []
    inputs = {item["gene"]: item["cds"] for item in document["inputs"]}
    for gene in document["meta"]["gene_order"]:
        chunks.append(f">{gene}|input\n{inputs[gene]}\n")
    for row in document["rows"]:
        if row.get("sequence"):
            chunks.append(f">{row['gene']}|{row['tool']}\n{row['sequence']}\n")
    path.write_text("".join(chunks), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--doc", type=Path)
    parser.add_argument("--fasta", type=Path)
    parser.add_argument(
        "--max-genes", type=int, default=0, help="debug cap; refuses the official report"
    )
    parser.add_argument("--chisel-worker", action="store_true")
    args = parser.parse_args()
    if args.chisel_worker:
        chisel_worker()
        return
    if args.out is None or args.doc is None:
        raise SystemExit("need --out and --doc")
    if not SCORE_PYTHON.is_file():
        raise SystemExit(
            "missing scoring interpreter. Create it with: "
            "python3 -m virtualenv benchmarks/gate_ben_001_chisel/.score-venv && "
            "benchmarks/gate_ben_001_chisel/.score-venv/bin/pip install 'biopython==1.81'"
        )

    from dna_design_suite import __version__ as package
    from dna_design_suite import core_version
    from dna_design_suite.codon.bench.corpus import load_panel, ribosomal_usage, select_genes

    weights = ribosomal_usage()
    file_weights = {
        codon: float(row["relative_adaptiveness"])
        for codon, row in json.loads(WEIGHTS_PATH.read_text(encoding="utf-8"))["codon_usage"][
            "table"
        ].items()
    }
    if weights != file_weights:
        raise SystemExit("ribosomal_usage() does not match examples/mg1655_ribosomal_usage.json")
    constraints = json.loads(CONSTRAINTS_PATH.read_text(encoding="utf-8"))
    mapping = translate_constraints(constraints)
    panel = load_panel()
    genes = select_genes(panel)
    if args.max_genes:
        genes = genes[: args.max_genes]
    profile, digest = profile_for(weights)
    usage = codon_usage_table(weights)
    positives = [weight for weight in weights.values() if weight > 0]
    family_max_is_one = True
    for _amino, codons in usage.items():
        if abs(max(codons.values()) - 1.0) > 1e-12:
            family_max_is_one = False
    store_dir = Path("/tmp/dds-gate-ben-001-stores")
    store_dir.mkdir(parents=True, exist_ok=True)

    # Pay the first-import cost outside the timed genes.
    warmup_cds = genes[0]["cds"][:30]
    if len(warmup_cds) % 3:
        warmup_cds = warmup_cds[: len(warmup_cds) - (len(warmup_cds) % 3)]
    suite_job(warmup_cds, profile, constraints, digest, store_dir, "_warmup")

    rows: list[dict[str, Any]] = []
    for index, gene in enumerate(genes, start=1):
        cds = gene["cds"].upper()
        print(f"[{index}/{len(genes)}] {gene['gene']} {len(cds)} nt", file=sys.stderr, flush=True)
        suite = suite_job(cds, profile, constraints, digest, store_dir, gene["gene"])
        print(
            f"  suite ok={suite['ok']} solver_s={suite['solver_s']}",
            file=sys.stderr,
            flush=True,
        )
        chisel = run_chisel(cds, mapping["specs"], usage)
        extra = chisel_diagnostics(chisel.get("sequence"), constraints, profile, digest)
        print(
            "  chisel ok={} solver_s={} wall_s={}".format(
                chisel.get("ok"),
                chisel.get("solver_s"),
                chisel.get("wall_s"),
            ),
            file=sys.stderr,
            flush=True,
        )
        common = {
            "gene": gene["gene"],
            "protein_id": gene["protein_id"],
            "nt": len(cds),
            "codons": len(cds) // 3,
        }
        rows.append(
            {
                **common,
                "tool": "dna_design_suite",
                "sequence": suite.get("sequence"),
                "tool_ok": suite["ok"],
                "state": suite.get("state"),
                "n_candidates": suite.get("n_candidates"),
                "error_type": suite.get("error_type"),
                "error": suite.get("error"),
                "solver_s": suite.get("solver_s"),
                "wall_s": suite.get("wall_s"),
                "warnings": suite.get("warnings"),
                "diagnostic_usage_cai": suite.get("diagnostic_usage_cai"),
                "diagnostic_product_checker": suite.get("diagnostic_product_checker"),
            }
        )
        rows.append(
            {
                **common,
                "tool": "dnachisel",
                "sequence": chisel.get("sequence"),
                "tool_ok": bool(chisel.get("ok")),
                "state": None,
                "error_type": chisel.get("error_type"),
                "error": chisel.get("error"),
                "solver_s": chisel.get("solver_s"),
                "wall_s": chisel.get("wall_s"),
                "initialize_s": chisel.get("initialize_s"),
                "resolve_s": chisel.get("resolve_s"),
                "optimize_s": chisel.get("optimize_s"),
                "chisel_all_constraints_pass": chisel.get("all_constraints_pass"),
                "chisel_constraint_checks": chisel.get("constraint_checks"),
                "stderr": chisel.get("stderr"),
                "diagnostic_usage_cai": extra["diagnostic_usage_cai"],
                "diagnostic_product_checker": extra["diagnostic_product_checker"],
            }
        )

    import numpy

    document = {
        "meta": {
            "partial": bool(args.max_genes),
            "git_commit": git_commit(),
            "panel_id": panel["id"],
            "accession": panel["accession"],
            "selected_cds_sha256": panel["selected_cds_sha256"],
            "n_genes": len(genes),
            "gene_order": [gene["gene"] for gene in genes],
            "host": HOST,
            "translation_table": 11,
            "weights_path": "examples/mg1655_ribosomal_usage.json",
            "weights_sha256": sha256_file(WEIGHTS_PATH),
            "constraints_path": "benchmarks/validation/constraints_assembly.json",
            "constraints_sha256": sha256_file(CONSTRAINTS_PATH),
            "constraint_ids": [item["id"] for item in constraints["constraints"]],
            "constraint_translation": {
                "translated": mapping["translated"],
                "not_translated": mapping["not_translated"],
                "notes": mapping["notes"],
            },
            "schema_kinds_not_in_file": SCHEMA_KINDS_NOT_IN_FILE,
            "zero_weight_codons": sorted(
                codon for codon, weight in weights.items() if weight == 0.0
            ),
            "min_positive_weight": min(positives),
            "family_max_is_one": family_max_is_one,
            "zero_floor_does_not_invert_positive_codons": min(positives) > CHISEL_ZERO_FLOOR,
            "dnachisel_zero_floor": CHISEL_ZERO_FLOOR,
            "seeds": {
                "dna_design_suite": SEED,
                "dnachisel_numpy": SEED,
                "dnachisel_random": SEED,
                "pythonhashseed": SEED,
            },
            "suite_job": {
                "budget": BUDGET,
                "generations": GENERATIONS,
                "segment": SEGMENT,
                "stagnation_window": STAGNATION_WINDOW,
                "objectives": ["usage.cai"],
            },
            "chisel_job": {
                "method": "use_best_codon",
                "genetic_table": "Bacterial",
                "start_codon": "ATG",
                "max_random_iters": 1000,
                "mutations_per_iteration": 2,
                "optimization_stagnation_tolerance": 100,
                "randomization_threshold": 10000,
                "local_extensions": [0, 5],
                "logger": None,
                "timeout_s": CHISEL_TIMEOUT_S,
            },
            "versions": {
                "dna_design_suite": package,
                "core": core_version(),
                "dnachisel": package_version("dnachisel"),
                "numpy": numpy.__version__,
                "biopython_project": package_version("biopython"),
                "python_design": platform.python_version(),
                "python_codon_tables": package_version("python_codon_tables"),
            },
            "machine": machine(),
            "negative_control": (
                "Negative control, before this measurement. "
                "`independent_score.py --self-test` passed. The reverse-complement "
                "count in `_motif_hits` was then replaced with `hits += 0` and the "
                "self-test exited 1 with `BsaI reverse complement GAGACC was not "
                "counted`. That edit was reverted and the self-test passed again. "
                "The fixture is the antisense BsaI site, so a checker that only "
                "reads the forward strand does not pass."
            ),
        },
        "inputs": [
            {
                "gene": gene["gene"],
                "protein_id": gene["protein_id"],
                "codons": gene["codons"],
                "cds": gene["cds"].upper(),
            }
            for gene in genes
        ],
        "rows": rows,
    }
    raw_path = args.out.with_suffix(".unscored.json")
    raw_path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.max_genes:
        print(f"partial run written to {raw_path}; not scored into the report", file=sys.stderr)
        return
    completed = subprocess.run(
        [
            str(SCORE_PYTHON),
            str(SCORER),
            "--score",
            str(raw_path),
            "--out",
            str(args.out),
            "--doc",
            str(args.doc),
        ],
        cwd=ROOT,
        check=False,
    )
    if completed.returncode != 0:
        raise SystemExit(f"independent scorer exited {completed.returncode}")
    scored = json.loads(args.out.read_text(encoding="utf-8"))
    fasta = args.fasta or args.out.with_name("gate-ben-001-dnachisel-sequences.fasta")
    write_fasta(fasta, scored)
    raw_path.unlink(missing_ok=True)
    print(f"wrote {args.out}", file=sys.stderr)
    print(f"wrote {args.doc}", file=sys.stderr)
    print(f"wrote {fasta}", file=sys.stderr)


if __name__ == "__main__":
    main()
