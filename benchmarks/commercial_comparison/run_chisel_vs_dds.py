#!/usr/bin/env python3
"""DNA Design Suite against DNA Chisel on a fixed, predeclared gene list.

The protocol matches the study in ``docs/benchmarks/v2.9.3-external-comparison.md``
(pull request 31): host ``escherichia_coli:511145/MG1655``, ribosomal relative
adaptiveness, budget ``standard``, 4 generations, stagnation window 21, seed 7.
DNA Chisel receives that same table through ``CodonOptimize(method="use_best_codon")``
and the same constraint language, mapped by ``tools/external_comparison.py``.

Nothing in this script retunes a gene. The corpus rows are the first record
in each length bin of ``benchmarks/corpus/generality_sample.json`` with at
least 150 nt, excluding the four named genes, ordered by ``(nt, gene)``.

Commercial web tools are not called. From the repository root, after
``.venv/bin/maturin develop --release`` and
``.venv/bin/pip install 'dnachisel==3.2.16' 'ViennaRNA==2.7.2'``::

    .venv/bin/python benchmarks/commercial_comparison/run_chisel_vs_dds.py \\
        --out docs/benchmarks/commercial-comparison-rows.json
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import benchmark_common as common  # noqa: E402
import external_comparison as ext  # noqa: E402
from dna_design_suite import __version__ as PACKAGE_VERSION  # noqa: E402
from dna_design_suite import _core  # noqa: E402
from dna_design_suite.codon.bench.baselines import cai_greedy  # noqa: E402
from dna_design_suite.modelreg import ArtifactStore  # noqa: E402
from generality_corpus import load_sample  # noqa: E402

NAMED = ("rpsL", "waaS", "lacZ", "yeeJ")
CORPUS_MIN_NT = 150
LENGTH_BINS = ("short", "medium", "long", "very_long")
# Same job the v2.9.3 external comparison used. Not chosen per gene.
BUDGET = "standard"
GENERATIONS = 4
STAGNATION_WINDOW = 21
SEED = 7
CHISEL_TIMEOUT_S = 180


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
        "python": platform.python_version(),
        "cpu": cpu,
        "nproc": os.cpu_count(),
        "mem_total_kib": mem_kib,
    }


def corpus_genes(sample: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """One gene per length bin. The rule is the sort key, not a score."""
    chosen = []
    for name in LENGTH_BINS:
        eligible = [
            row
            for row in sample
            if row["length_bin"] == name
            and int(row["nt"]) >= CORPUS_MIN_NT
            and row["gene"] not in NAMED
        ]
        eligible.sort(key=lambda row: (int(row["nt"]), row["gene"]))
        if not eligible:
            raise SystemExit(f"no corpus gene in bin {name}")
        row = eligible[0]
        chosen.append(
            {
                "gene": row["gene"],
                "panel": "corpus",
                "length_bin": row["length_bin"],
                "protein_id": row["protein_id"],
                "source": "benchmarks/corpus/generality_sample.json",
                "cds": row["cds"],
                "nt": len(row["cds"]),
                "codons": len(row["cds"]) // 3,
            }
        )
    return chosen


def named_genes() -> list[dict[str, Any]]:
    loaded = {row["gene"]: row for row in common.load_genes()}
    genes = []
    for name in NAMED:
        row = loaded[name]
        genes.append(
            {
                "gene": name,
                "panel": "named",
                "length_bin": None,
                "protein_id": None,
                "source": row["path"],
                "cds": row["cds"],
                "nt": row["nt"],
                "codons": row["codons"],
            }
        )
    return genes


def flat_metrics(row: dict[str, Any]) -> dict[str, Any]:
    checker = row.get("our_checker") or {}
    by_id = checker.get("actionable_by_id") or {}
    return {
        "cai": row.get("independent_cai"),
        "usage_cai": row.get("usage_cai"),
        "cai_abs_delta": row.get("cai_abs_delta"),
        "gc": row.get("gc"),
        "protein_relation": row.get("protein_relation"),
        "feasible": checker.get("feasible"),
        "gc_window_violations": by_id.get("c_gc50", 0),
        "actionable_by_id": by_id,
        "n_actionable": checker.get("n_actionable"),
        "wall_s": row.get("wall_s"),
        "accessibility_penalty": row.get("accessibility_penalty"),
        "sha256": row.get("sha256"),
        "tool_ok": row.get("tool_ok"),
        "tool_error_type": row.get("tool_error_type"),
    }


def pair_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """DNA Design Suite versus DNA Chisel on rows where the comparison is the same CDS."""
    index: dict[tuple[str, str], dict[str, dict[str, Any]]] = {}
    for row in rows:
        if row["tool"] not in {"dna_design_suite", "dnachisel"}:
            continue
        if row.get("probe"):
            continue
        index.setdefault((row["constraint_set"], row["gene"]), {})[row["tool"]] = row
    summary = []
    for (constraint, gene), pair in sorted(index.items()):
        left = pair.get("dna_design_suite")
        right = pair.get("dnachisel")
        if left is None or right is None:
            continue
        left_cai = left.get("independent_cai")
        right_cai = right.get("independent_cai")
        relation = None
        if (
            left.get("tool_ok")
            and right.get("tool_ok")
            and left_cai is not None
            and right_cai is not None
        ):
            if left_cai > right_cai:
                relation = "dna_design_suite_higher"
            elif right_cai > left_cai:
                relation = "dnachisel_higher"
            else:
                relation = "tie"
        elif left.get("tool_ok") and not right.get("tool_ok"):
            relation = "dnachisel_no_sequence"
        elif right.get("tool_ok") and not left.get("tool_ok"):
            relation = "dna_design_suite_no_sequence"
        else:
            relation = "neither_sequence"
        summary.append(
            {
                "constraint_set": constraint,
                "gene": gene,
                "panel": left.get("panel"),
                "nt": left.get("nt"),
                "relation": relation,
                "dna_design_suite_cai": left_cai,
                "dnachisel_cai": right_cai,
                "dna_design_suite_feasible": (left.get("our_checker") or {}).get("feasible"),
                "dnachisel_feasible": (right.get("our_checker") or {}).get("feasible"),
                "dnachisel_self_pass": (right.get("extra") or {}).get("all_pass"),
                "dna_design_suite_chisel_pass": (left.get("chisel_evaluate") or {}).get("all_pass"),
                "protein_dna_design_suite": left.get("protein_relation"),
                "protein_dnachisel": right.get("protein_relation"),
            }
        )
    return summary


def fmt_num(value: Any) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        return json.dumps(value)
    return str(value)


def feas_text(row: dict[str, Any]) -> str:
    if not row.get("tool_ok"):
        return f"no sequence ({row.get('tool_error_type')})"
    checker = row.get("our_checker") or {}
    if checker.get("feasible"):
        text = "pass"
    else:
        ids = checker.get("actionable_by_id") or {}
        text = "fail " + ",".join(f"{key}:{ids[key]}" for key in ids)
    chisel = row.get("chisel_evaluate")
    if chisel is not None:
        text += ", DNA Chisel pass" if chisel.get("all_pass") else ", DNA Chisel fail"
    elif row["tool"] == "dnachisel":
        passed = (row.get("extra") or {}).get("all_pass")
        if passed is True:
            text += ", DNA Chisel pass"
        elif passed is False:
            text += ", DNA Chisel fail"
    return text


def access_text(value: Any) -> str:
    if isinstance(value, dict):
        return str(value.get("error", "—"))
    return fmt_num(value)


def md_row(cells: list[Any]) -> str:
    return "| " + " | ".join(str(cell) for cell in cells) + " |"


def render_tables(document: dict[str, Any]) -> str:
    lines = [
        "<!-- generated by benchmarks/commercial_comparison/run_chisel_vs_dds.py -->",
        "",
    ]
    rows = [row for row in document["rows"] if not row.get("probe")]
    constraints = []
    for row in rows:
        if row["constraint_set"] not in constraints:
            constraints.append(row["constraint_set"])
    for constraint in constraints:
        lines.append(f"## {constraint}")
        lines.append("")
        header = [
            "gene",
            "panel",
            "nt",
            "tool",
            "CAI",
            "GC",
            "protein",
            "feasibility",
            "wall s",
            "accessibility",
        ]
        lines.append(md_row(header))
        lines.append(
            md_row(["---", "---", "---:", "---", "---:", "---:", "---", "---", "---:", "---:"])
        )
        subset = [row for row in rows if row["constraint_set"] == constraint]
        for row in subset:
            lines.append(
                md_row(
                    [
                        row["gene"],
                        row["panel"],
                        row["nt"],
                        row["tool"],
                        fmt_num(row.get("independent_cai")),
                        fmt_num(row.get("gc")),
                        row.get("protein_relation") or "—",
                        feas_text(row).replace("|", "/"),
                        fmt_num(row.get("wall_s")),
                        access_text(row.get("accessibility_penalty")),
                    ]
                )
            )
        lines.append("")
    lines.append("## DNA Chisel seed probe, rpsL, assembly")
    lines.append("")
    lines.append("| seed | CAI | GC | sha256 | wall s | feasible |")
    lines.append("| ---: | ---: | ---: | --- | ---: | --- |")
    for row in document["rows"]:
        if not row.get("probe"):
            continue
        lines.append(
            "| {seed} | {cai} | {gc} | `{sha}` | {wall} | {feas} |".format(
                seed=row["extra"].get("chisel_seed"),
                cai=fmt_num(row.get("independent_cai")),
                gc=fmt_num(row.get("gc")),
                sha=row.get("sha256"),
                wall=fmt_num(row.get("wall_s")),
                feas=feas_text(row).replace("|", "/"),
            )
        )
    lines.append("")
    lines.append("## Pair relation")
    lines.append("")
    counts: dict[str, int] = {}
    for item in document["pair_summary"]:
        counts[item["relation"]] = counts.get(item["relation"], 0) + 1
    for key in sorted(counts):
        lines.append(f"- `{key}`: {counts[key]}")
    lines.append("")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--tables",
        type=Path,
        help="Also write the markdown tables. Default: <out stem>-tables.md next to --out.",
    )
    args = parser.parse_args()

    weights = common.load_weights()
    profile = common.host_profile(weights)
    digest = common.profile_digest(profile)
    codon_usage = common.codon_usage_by_amino_acid(weights)
    floor_ok = common.minimum_positive_weight(weights) > ext.DNA_CHISEL_ZERO_FLOOR
    sets = common.load_constraint_sets()
    sample = load_sample()
    genes = named_genes() + corpus_genes(sample)
    store = ArtifactStore(Path("/tmp/dds-commercial-comparison-store"))

    print("mapping checks", flush=True)
    checks = ext.mapping_checks(weights, codon_usage)

    objective_dds = (
        "DNA Design Suite usage.cai, budget standard, generations 4, "
        "stagnation window 21, seed 7, host escherichia_coli:511145/MG1655, "
        "translation table 11. Same job shape as docs/benchmarks/v2.9.3-external-comparison.md."
    )
    objective_chisel = (
        "DNA Chisel MaximizeCAI via CodonOptimize(method='use_best_codon') on the "
        "ribosomal relative-adaptiveness table. Reported CAI is the independent "
        "Sharp and Li mean, which excludes Met, Trp and stops. "
        f"Zero floor inside DNA Chisel is {ext.DNA_CHISEL_ZERO_FLOOR}; "
        f"minimum positive ribosomal weight is {common.minimum_positive_weight(weights)}; "
        f"floor does not reorder positive codons: {floor_ok}."
    )
    objective_greedy = (
        "cai_greedy picks the highest relative adaptiveness at each codon and breaks "
        "ties toward the lexicographically smaller codon. It does not see constraints."
    )

    rows: list[dict[str, Any]] = []
    specs = {item["name"]: ext.chisel_spec(item["document"]) for item in sets}

    def add_row(
        gene: dict[str, Any],
        tool: str,
        constraint_name: str,
        raw: dict[str, Any],
        note: str,
        constraints: dict[str, Any],
    ) -> dict[str, Any]:
        row = ext.finish_tool_row(
            gene, tool, constraint_name, raw, weights, profile, digest, constraints, note
        )
        row["panel"] = gene["panel"]
        row["length_bin"] = gene["length_bin"]
        row["protein_id"] = gene["protein_id"]
        row["source"] = gene["source"]
        row["metrics"] = flat_metrics(row)
        return row

    for gene in genes:
        started = time.perf_counter()
        greedy = cai_greedy(gene["cds"], common.TRANSLATION_TABLE, weights)
        greedy_wall = time.perf_counter() - started
        for constraint in sets:
            print(f"dna_design_suite {gene['gene']} {constraint['name']}", flush=True)
            designed = common.run_codonforge(
                gene["cds"],
                constraint["document"],
                profile,
                digest,
                store,
                budget=BUDGET,
                seed=SEED,
                generations=GENERATIONS,
                stagnation_window=STAGNATION_WINDOW,
                job_id=f"commercial-{gene['gene']}-{constraint['name']}",
            )
            designed["sequence"] = designed.get("cds")
            designed["ok"] = designed["state"] == "succeeded" and bool(designed.get("cds"))
            dds_row = add_row(
                gene,
                "dna_design_suite",
                constraint["name"],
                designed,
                objective_dds,
                constraint["document"],
            )
            if designed.get("sequence"):
                evaluated = ext.run_chisel(
                    designed["sequence"],
                    specs[constraint["name"]]["specs"],
                    codon_usage,
                    seed=SEED,
                    timeout_s=60,
                    mode="evaluate",
                )
                dds_row["chisel_evaluate"] = {
                    "ok": evaluated.get("ok"),
                    "all_pass": evaluated.get("all_pass"),
                    "error_type": evaluated.get("error_type"),
                    "failing": [
                        item["label"]
                        for item in evaluated.get("scores") or []
                        if item.get("passes") is False
                    ],
                }
            rows.append(dds_row)

            print(f"dnachisel {gene['gene']} {constraint['name']}", flush=True)
            chisel = ext.run_chisel(
                gene["cds"],
                specs[constraint["name"]]["specs"],
                codon_usage,
                seed=SEED,
                timeout_s=CHISEL_TIMEOUT_S,
            )
            chisel_row = add_row(
                gene,
                "dnachisel",
                constraint["name"],
                chisel,
                objective_chisel,
                constraint["document"],
            )
            chisel_row["mapping_comparable_ids"] = specs[constraint["name"]]["comparable_ids"]
            chisel_row["extra"]["chisel_seed"] = SEED
            rows.append(chisel_row)

            rows.append(
                add_row(
                    gene,
                    "cai_greedy",
                    constraint["name"],
                    {"ok": True, "sequence": greedy, "wall_s": greedy_wall, "state": "succeeded"},
                    objective_greedy,
                    constraint["document"],
                )
            )
            rows.append(
                add_row(
                    gene,
                    "native_input",
                    constraint["name"],
                    {"ok": True, "sequence": gene["cds"], "wall_s": 0.0, "state": "input"},
                    "The input CDS. Not an optimiser.",
                    constraint["document"],
                )
            )

    probe_gene = next(gene for gene in genes if gene["gene"] == "rpsL")
    assembly = next(item for item in sets if item["name"] == "assembly")
    for seed in (1, 7, 11):
        print(f"dnachisel seed probe rpsL assembly {seed}", flush=True)
        chisel = ext.run_chisel(
            probe_gene["cds"],
            specs["assembly"]["specs"],
            codon_usage,
            seed=seed,
            timeout_s=CHISEL_TIMEOUT_S,
        )
        row = add_row(
            probe_gene,
            "dnachisel_seed_probe",
            "assembly",
            chisel,
            objective_chisel + f" Probe seed {seed}.",
            assembly["document"],
        )
        row["probe"] = True
        row["extra"]["chisel_seed"] = seed
        rows.append(row)

    try:
        import dnachisel

        chisel_version = dnachisel.__version__
    except ImportError:
        chisel_version = None

    document: dict[str, Any] = {
        "study": "commercial codon-tool comparison, runnable arm",
        "package_version": PACKAGE_VERSION,
        "core_version": _core.version(),
        "optimizer_commit": "db867d4aa31a76229098393dbd0448b01a341032",
        "optimizer_commit_note": (
            "main at the time of the run. No file under crates/ or python/ "
            "was edited for this measurement."
        ),
        "git_head_at_runtime": common.git_head(),
        "machine": machine(),
        "build": ".venv/bin/maturin develop --release",
        "versions": {
            "package": PACKAGE_VERSION,
            "core": _core.version(),
            "python": platform.python_version(),
            "viennarna": common.tool_versions().get("viennarna"),
            "viennarna_accessibility": common.tool_versions()["viennarna_accessibility"],
        },
        "dnachisel": chisel_version,
        "biopython": __import__("Bio").__version__,
        "numpy": __import__("numpy").__version__,
        "weights": str(common.WEIGHTS_PATH.relative_to(common.ROOT)),
        "host": common.HOST,
        "translation_table": common.TRANSLATION_TABLE,
        "job": {
            "budget": BUDGET,
            "generations": GENERATIONS,
            "stagnation_window": STAGNATION_WINDOW,
            "seed": SEED,
            "objectives": ["usage.cai"],
        },
        "corpus_rule": (
            f"First gene in each of {list(LENGTH_BINS)} from generality_sample.json "
            f"with nt >= {CORPUS_MIN_NT}, gene not in {list(NAMED)}, sorted by (nt, gene)."
        ),
        "corpus_sample_sha256": json.loads(
            (ROOT / "benchmarks" / "corpus" / "generality_sample.json").read_text(encoding="utf-8")
        )["sample_sha256"],
        "genes": [
            {
                key: gene[key]
                for key in ("gene", "panel", "length_bin", "protein_id", "source", "nt", "codons")
            }
            for gene in genes
        ],
        "dnachisel_zero_floor": ext.DNA_CHISEL_ZERO_FLOOR,
        "minimum_positive_weight": common.minimum_positive_weight(weights),
        "zero_floor_does_not_invert_positive_codons": floor_ok,
        "wall_time_note": (
            "DNA Chisel wall_s is time.perf_counter around a fresh interpreter "
            "(tools/external_comparison.py run_chisel). It includes process startup. "
            "DNA Design Suite wall_s is the design() call inside the already loaded process."
        ),
        "accessibility_note": (
            "accessibility_penalty is dna_design_suite.codon.orchestrator.structure_diag.evaluate "
            "after the design: mean -log unpaired probability, window 6 nt, span 40 nt, 37 C. "
            "It was not an objective of either optimiser in this run."
        ),
        "mapping_checks": checks,
        "rows": rows,
    }
    document["pair_summary"] = pair_summary(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tables_path = args.tables or args.out.with_name("commercial-comparison-tables.md")
    tables_path.write_text(render_tables(document), encoding="utf-8")
    print(f"wrote {args.out} rows {len(rows)}", flush=True)
    print(f"wrote {tables_path}", flush=True)


if __name__ == "__main__":
    main()
