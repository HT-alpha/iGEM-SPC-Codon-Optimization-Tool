#!/usr/bin/env python3
"""Broad codon-design corpus. Reproducible from the shipped CDS gzip.

The sequences are the CDS records in ``benchmarks/corpus/NC_000913.3_cds.fna.gz``
(Escherichia coli K-12 MG1655, accession NC_000913.3). A record whose header
does not carry that accession, a protein_id, a gene name, and a location span
equal to the sequence length is refused. Nothing is reverse-translated and
nothing is written in by hand.

The sample is a stratified draw, seed 20261018, across length bins and GC
tertiles, plus the rare-codon-richest records under the ribosomal weight
table. The same list is reloaded from ``benchmarks/corpus/generality_sample.json``
when that file is present, and the script checks its digest against the gzip.

Presets, each a budget and a constraint set:

* ``empty`` — no hard constraints, ``usage.cai``, budget ``fast``, 0 generations
* ``assembly`` — the shipped assembly file, same budget
* ``harsh`` — assembly plus a 30 nt GC band, homopolymer 5, and a 10-mer cap of 2
* ``structure`` — assembly constraints, ``usage.cai`` and
  ``structure.local_accessibility``, budget ``fast``, 2 generations

``standard`` repeats empty and assembly with budget ``standard`` and 12
generations, so tier 3 actually runs. The Kazusa O157:H7 table is the other
shipped usage file; it is scored on empty and assembly.

From the repository root, after ``make build``::

    .venv/bin/python tools/generality_corpus.py sample
    .venv/bin/python tools/generality_corpus.py run --out /tmp/generality
    .venv/bin/python tools/generality_corpus.py summarize /tmp/generality/rows.jsonl
    .venv/bin/python tools/generality_corpus.py compare before.jsonl after.jsonl
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import multiprocessing as mp
import os
import random
import sys
import time
import traceback
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SAMPLE_PATH = ROOT / "benchmarks" / "corpus" / "generality_sample.json"
WEIGHTS_RIBO = ROOT / "examples" / "mg1655_ribosomal_usage.json"
WEIGHTS_KAZUSA = ROOT / "examples" / "kazusa_e_coli_316407.json"
ASSEMBLY_PATH = ROOT / "benchmarks" / "validation" / "constraints_assembly.json"

# Fixed draw. Not fitted to a gene, and not re-rolled when a row looks bad.
SAMPLE_SEED = 20261018
SAMPLE_TARGET = 220
HOST = "escherichia_coli:511145/MG1655"
TABLE = 11

# Length bins in nucleotides. Edges are round and cover the file; they are
# not the length of any one gene.
LENGTH_BINS: tuple[tuple[str, int, int], ...] = (
    ("short", 0, 300),
    ("medium", 300, 900),
    ("long", 900, 1800),
    ("very_long", 1800, 10**9),
)
BIN_QUOTA = {"short": 40, "medium": 70, "long": 70, "very_long": 40}


def _load_weights(path: Path) -> dict[str, float]:
    sys.path.insert(0, str(ROOT / "benchmarks" / "validation"))
    import cai_independent

    return cai_independent.load_ribosomal_weights(path)


def _rare_cutoff(weights: dict[str, float]) -> float:
    """10th percentile of the positive weights in the table.

    A counted codon at or below that weight is treated as rare for the
    stratification only. The cutoff is the table's own percentile, not a
    constant chosen for one gene.
    """
    positive = sorted(weight for weight in weights.values() if weight > 0.0)
    if not positive:
        raise SystemExit("usage table has no positive weight")
    index = max(0, math.ceil(0.10 * len(positive)) - 1)
    return positive[index]


def _gc(seq: str) -> float:
    return sum(base in "GC" for base in seq) / len(seq)


def _rare_fraction(seq: str, weights: dict[str, float], cutoff: float) -> float:
    sys.path.insert(0, str(ROOT / "benchmarks" / "validation"))
    import cai_independent

    scored = cai_independent.score(seq, weights)
    counted = int(scored["counted_codons"])
    if counted == 0:
        return 0.0
    rare = 0
    for index in range(0, len(seq), 3):
        codon = seq[index : index + 3]
        amino = cai_independent.AA[codon]
        if amino in {"*", "M", "W"}:
            continue
        # A zero weight is rarer than the percentile cutoff. Count it.
        if weights[codon] <= cutoff:
            rare += 1
    return rare / counted


def _length_bin(nt: int) -> str:
    for name, lo, hi in LENGTH_BINS:
        if lo <= nt < hi:
            return name
    return "very_long"


def _valid_cds(header: str, sequence: str) -> dict[str, Any] | None:
    from dna_design_suite.codon.bench.corpus import _gene, _header_problem

    problem = _header_problem(header, sequence)
    if problem:
        return None
    sequence = sequence.upper().replace("U", "T")
    if any(base not in "ACGT" for base in sequence):
        return None
    if not sequence.startswith(("ATG", "GTG", "TTG")):
        return None
    if sequence[-3:] not in {"TAA", "TAG", "TGA"}:
        return None
    from dna_design_suite import _core

    try:
        # from_cds refuses an internal stop. A record that fails that check
        # is left out; the sequence is not edited to remove the stop.
        _core.translate(sequence, TABLE)
    except ValueError:
        return None
    gene = _gene(header)
    if not gene:
        return None
    protein_id = header.split("[protein_id=", 1)[1].split("]", 1)[0]
    return {
        "gene": gene,
        "protein_id": protein_id,
        "nt": len(sequence),
        "cds": sequence,
    }


def load_records() -> list[dict[str, Any]]:
    from dna_design_suite.codon.bench.corpus import (
        CDS_PATH,
        iter_cds,
        load_panel,
        verify_shipped_cds,
    )

    panel = load_panel()
    verify_shipped_cds(panel)
    found: dict[str, dict[str, Any]] = {}
    order: list[str] = []
    for header, sequence in iter_cds(CDS_PATH):
        row = _valid_cds(header, sequence)
        if row is None or row["gene"] in found:
            continue
        found[row["gene"]] = row
        order.append(row["gene"])
    return [found[gene] for gene in order]


def _tertile_pick(
    rows: list[dict[str, Any]],
    quota: int,
    rng: random.Random,
) -> list[dict[str, Any]]:
    if len(rows) <= quota:
        return list(rows)
    ordered = sorted(rows, key=lambda row: (row["gc"], row["gene"]))
    thirds = [ordered[: len(ordered) // 3], ordered[len(ordered) // 3 : 2 * len(ordered) // 3]]
    thirds.append(ordered[2 * len(ordered) // 3 :])
    base, extra = divmod(quota, 3)
    picked: list[dict[str, Any]] = []
    for index, group in enumerate(thirds):
        take = base + (1 if index < extra else 0)
        take = min(take, len(group))
        picked.extend(rng.sample(group, take))
    if len(picked) < quota:
        chosen = {row["gene"] for row in picked}
        rest = [row for row in ordered if row["gene"] not in chosen]
        picked.extend(rng.sample(rest, min(quota - len(picked), len(rest))))
    return picked


def draw_sample(records: list[dict[str, Any]], weights: dict[str, float]) -> list[dict[str, Any]]:
    cutoff = _rare_cutoff(weights)
    enriched = []
    for row in records:
        item = {
            "gene": row["gene"],
            "protein_id": row["protein_id"],
            "nt": row["nt"],
            "gc": _gc(row["cds"]),
            "rare_fraction": _rare_fraction(row["cds"], weights, cutoff),
            "length_bin": _length_bin(row["nt"]),
            "cds": row["cds"],
        }
        enriched.append(item)
    rng = random.Random(SAMPLE_SEED)
    by_bin: dict[str, list[dict[str, Any]]] = {name: [] for name, _, _ in LENGTH_BINS}
    for row in enriched:
        by_bin[row["length_bin"]].append(row)
    selected: dict[str, dict[str, Any]] = {}
    for name, _, _ in LENGTH_BINS:
        for row in _tertile_pick(by_bin[name], BIN_QUOTA[name], rng):
            selected[row["gene"]] = row
    # The rare-codon tail is added whole when it is not already in the draw,
    # so the sample cannot miss that class by chance.
    rare_tail = sorted(enriched, key=lambda row: (-row["rare_fraction"], row["gene"]))[:25]
    for row in rare_tail:
        selected.setdefault(row["gene"], row)
    # Shortest and longest valid CDS, so the span is the file's span.
    by_length = sorted(enriched, key=lambda row: (row["nt"], row["gene"]))
    for row in by_length[:10] + by_length[-10:]:
        selected.setdefault(row["gene"], row)
    ordered = sorted(selected.values(), key=lambda row: (row["nt"], row["gene"]))
    if len(ordered) < 200:
        raise SystemExit(f"sample has {len(ordered)} genes; 200 were required")
    return ordered


def sample_digest(rows: list[dict[str, Any]]) -> str:
    payload = "".join(f"{row['gene']}\n{row['protein_id']}\n{row['cds']}\n" for row in rows)
    return hashlib.sha256(payload.encode()).hexdigest()


def write_sample(rows: list[dict[str, Any]]) -> None:
    public = []
    for row in rows:
        public.append({key: row[key] for key in row if key != "cds"})
    document = {
        "seed": SAMPLE_SEED,
        "accession": "NC_000913.3",
        "n": len(public),
        "sample_sha256": sample_digest(rows),
        "rare_cutoff_rule": (
            "10th percentile of positive relative-adaptiveness values in "
            "examples/mg1655_ribosomal_usage.json; a counted codon at or below "
            "that weight counts as rare for stratification"
        ),
        "length_bins_nt": {name: [lo, hi] for name, lo, hi in LENGTH_BINS},
        "genes": public,
    }
    SAMPLE_PATH.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


def load_sample() -> list[dict[str, Any]]:
    if not SAMPLE_PATH.is_file():
        raise SystemExit(f"missing {SAMPLE_PATH}; run the sample subcommand first")
    document = json.loads(SAMPLE_PATH.read_text(encoding="utf-8"))
    if document.get("seed") != SAMPLE_SEED:
        raise SystemExit("sample file seed does not match the script")
    wanted = {row["gene"]: row for row in document["genes"]}
    records = {row["gene"]: row for row in load_records()}
    loaded = []
    for gene, meta in wanted.items():
        if gene not in records:
            raise SystemExit(f"{gene} is in the sample file and not in the gzip")
        cds = records[gene]["cds"]
        if records[gene]["protein_id"] != meta["protein_id"]:
            raise SystemExit(f"{gene} protein_id changed relative to the sample file")
        loaded.append({**meta, "cds": cds})
    digest = sample_digest(loaded)
    if digest != document["sample_sha256"]:
        raise SystemExit("sample digest does not match the gzip; the draw is not the cited one")
    return loaded


def _constraints(name: str) -> dict[str, Any]:
    if name == "empty":
        return {
            "schema_version": "1.0.0",
            "constraint_set_id": "cst_empty",
            "constraints": [],
        }
    document = json.loads(ASSEMBLY_PATH.read_text(encoding="utf-8"))
    if name == "assembly":
        return document
    if name == "harsh":
        document = json.loads(json.dumps(document))
        document["constraint_set_id"] = "cst_harsh_assembly"
        document["name"] = "assembly plus a 30 nt GC window, homopolymer 5, and a 10-mer cap of 2"
        document["constraints"] = list(document["constraints"]) + [
            {
                "id": "c_gc30",
                "class": "blocking",
                "scope": "cds",
                "rule": {"kind": "gc_window", "window_nt": 30, "min": 0.35, "max": 0.65},
            },
            {
                "id": "c_homo5",
                "class": "blocking",
                "scope": "cds",
                "rule": {"kind": "homopolymer", "base": "any", "max_run": 5},
            },
            {
                "id": "c_kmer10",
                "class": "blocking",
                "scope": "cds",
                "rule": {
                    "kind": "kmer_repeat",
                    "k": 10,
                    "max_occurrences": 2,
                    "count_overlapping": True,
                    "canonical_strand": False,
                },
            },
        ]
        return document
    raise SystemExit(f"unknown constraint preset {name}")


def preset_table() -> list[dict[str, Any]]:
    """Every budget/constraint preset the report runs."""
    return [
        {
            "preset": "empty",
            "budget": "fast",
            "generations": 0,
            "objectives": ["usage.cai"],
            "constraints": "empty",
            "hosts": ("ribosomal", "kazusa"),
        },
        {
            "preset": "assembly",
            "budget": "fast",
            "generations": 0,
            "objectives": ["usage.cai"],
            "constraints": "assembly",
            "hosts": ("ribosomal", "kazusa"),
        },
        {
            "preset": "harsh",
            "budget": "fast",
            "generations": 0,
            "objectives": ["usage.cai"],
            "constraints": "harsh",
            "hosts": ("ribosomal",),
        },
        {
            "preset": "structure",
            "budget": "fast",
            "generations": 2,
            "objectives": ["usage.cai", "structure.local_accessibility"],
            "constraints": "assembly",
            "hosts": ("ribosomal",),
        },
        {
            "preset": "empty_standard",
            "budget": "standard",
            "generations": 12,
            "objectives": ["usage.cai"],
            "constraints": "empty",
            "hosts": ("ribosomal",),
        },
        {
            "preset": "assembly_standard",
            "budget": "standard",
            "generations": 12,
            "objectives": ["usage.cai"],
            "constraints": "assembly",
            "hosts": ("ribosomal",),
        },
    ]


def _job_timeout(nt: int, preset: str) -> float:
    # A backstop for the runner, not a product constant. Long CDS get more
    # wall time; a job that is still inside the Rust seed when this fires is
    # recorded as a timeout rather than allowed to stall the corpus.
    scale = 0.03 if preset in {"harsh", "structure", "assembly_standard"} else 0.015
    return max(25.0, scale * nt)


def _run_design(payload: dict[str, Any]) -> dict[str, Any]:
    os.chdir(ROOT)
    sys.path.insert(0, str(ROOT))
    sys.path.insert(0, str(ROOT / "benchmarks" / "validation"))
    import cai_independent

    from dna_design_suite.codon.orchestrator.design import DesignError, DesignRequest, design
    from dna_design_suite.modelreg import ArtifactStore

    weights = cai_independent.load_ribosomal_weights(Path(payload["weights_path"]))
    profile = {
        "codon_usage": {
            "subset": "user_supplied",
            "table": {
                codon: {"relative_adaptiveness": weight}
                for codon, weight in sorted(weights.items())
            },
        }
    }
    digest = (
        "sha256:"
        + hashlib.sha256(
            json.dumps(profile, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
    )
    started = time.perf_counter()
    try:
        result = design(
            DesignRequest(
                host=HOST,
                cds=payload["cds"],
                host_profile=profile,
                objectives=list(payload["objectives"]),
                constraints=payload["constraints"],
                translation_table=TABLE,
                budget=payload["budget"],
                seed=7,
                generations=int(payload["generations"]),
                stagnation_window=8,
                segment=4,
                host_digest=digest,
            ),
            ArtifactStore(Path(payload["store"])),
        )
    except DesignError as exc:
        wall = time.perf_counter() - started
        return {
            "ok": False,
            "wall_s": wall,
            "code": exc.code,
            "detail": exc.detail,
            "state": "failed",
        }
    wall = time.perf_counter() - started
    status = result.status
    problem = status.get("problem") or {}
    code = problem.get("code")
    if result.front is None or status.get("state") != "succeeded":
        return {
            "ok": False,
            "wall_s": wall,
            "code": code,
            "detail": problem.get("detail"),
            "state": status.get("state"),
            "generations_run": result.job.generations_run,
        }
    candidate = result.front["candidates"][0]
    cds = candidate["sequences"][0]["cds"]
    report = candidate.get("constraint_report") or {}
    vector = candidate.get("objective_vector") or []
    cai = cai_independent.score(cds, weights)
    session = getattr(result.job, "_fold_session", None)
    folds = None
    if session is not None:
        folds = {
            "hits": session.hits,
            "full_folds": session.full_folds,
            "local_edits": session.local_edits,
        }
    return {
        "ok": True,
        "wall_s": wall,
        "state": status.get("state"),
        "code": code,
        "cds": cds,
        "cai": cai["cai"],
        "feasible": bool(report.get("feasible")),
        "n_actionable": sum(1 for item in report.get("violations") or [] if item.get("actionable")),
        "vector": vector,
        "generations_run": result.job.generations_run,
        "termination": (result.job.tier3 or {}).get("termination") if result.job.tier3 else None,
        "folds": folds,
        "protein_same": cai_independent.translate(cds) == cai_independent.translate(payload["cds"]),
    }


def _serve(in_q: mp.Queue, out_q: mp.Queue) -> None:
    while True:
        payload = in_q.get()
        if payload is None:
            return
        try:
            out_q.put(_run_design(payload))
        except Exception:
            out_q.put(
                {
                    "ok": False,
                    "wall_s": None,
                    "code": "RUNNER",
                    "detail": traceback.format_exc(limit=8),
                    "state": "failed",
                }
            )


def _percentile(values: list[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    index = min(len(ordered) - 1, max(0, math.ceil(fraction * len(ordered)) - 1))
    return ordered[index]


def summarize(path: Path) -> dict[str, Any]:
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line]
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(f"{row['preset']}|{row['host']}", []).append(row)
    report: dict[str, Any] = {"n_rows": len(rows), "presets": {}}
    for key, items in sorted(groups.items()):
        walls = [
            float(item["wall_s"]) for item in items if isinstance(item.get("wall_s"), (int, float))
        ]
        ok = [
            item
            for item in items
            if item.get("ok") and item.get("feasible") and item.get("protein_same")
        ]
        codes: dict[str, int] = {}
        for item in items:
            if item.get("ok") and item.get("feasible"):
                continue
            code = item.get("code") or item.get("state") or "unknown"
            codes[str(code)] = codes.get(str(code), 0) + 1
        deltas_greedy = [
            item["cai"] - item["cai_greedy"]
            for item in ok
            if item.get("cai") is not None and item.get("cai_greedy") is not None
        ]
        deltas_opt = [
            item["cai"] - item["cai_optimum"]
            for item in ok
            if item.get("cai") is not None and item.get("cai_optimum") is not None
        ]
        by_bin: dict[str, list[float]] = {}
        per_kb: dict[str, list[float]] = {}
        for item in items:
            if not isinstance(item.get("wall_s"), (int, float)):
                continue
            bin_name = item.get("length_bin") or "?"
            by_bin.setdefault(bin_name, []).append(float(item["wall_s"]))
            nt = item.get("nt")
            if isinstance(nt, int) and nt > 0:
                per_kb.setdefault(bin_name, []).append(float(item["wall_s"]) * 1000.0 / nt)
        terminations: dict[str, int] = {}
        for item in items:
            label = item.get("termination") or "none"
            terminations[str(label)] = terminations.get(str(label), 0) + 1
        actionable = sum(1 for item in items if (item.get("n_actionable") or 0) > 0)
        report["presets"][key] = {
            "n": len(items),
            "success": len(ok),
            "success_rate": len(ok) / len(items) if items else None,
            "failures": codes,
            "actionable_rows": actionable,
            "terminations": terminations,
            "wall_median_s": _percentile(walls, 0.50),
            "wall_p95_s": _percentile(walls, 0.95),
            "wall_max_s": max(walls) if walls else None,
            "wall_by_length_bin_median_s": {
                name: _percentile(values, 0.50) for name, values in sorted(by_bin.items())
            },
            "sec_per_kb_by_length_bin_median": {
                name: _percentile(values, 0.50) for name, values in sorted(per_kb.items())
            },
            "cai_minus_greedy_median": _percentile(deltas_greedy, 0.50),
            "cai_minus_greedy_min": min(deltas_greedy) if deltas_greedy else None,
            "cai_below_greedy": sum(1 for delta in deltas_greedy if delta < -1e-12),
            "cai_minus_optimum_median": _percentile(deltas_opt, 0.50),
            "cai_minus_optimum_min": min(deltas_opt) if deltas_opt else None,
            "cai_below_optimum": sum(1 for delta in deltas_opt if delta < -1e-12),
        }
    return report


def _format_report(report: dict[str, Any]) -> str:
    lines = [f"rows {report['n_rows']}"]
    for key, item in report["presets"].items():
        lines.append(
            f"{key}: success {item['success']}/{item['n']} "
            f"({item['success_rate']:.3f}) failures {item['failures']} "
            f"wall median/p95/max "
            f"{item['wall_median_s']} / {item['wall_p95_s']} / {item['wall_max_s']} "
            f"dCAI greedy min/median {item['cai_minus_greedy_min']} / "
            f"{item['cai_minus_greedy_median']} "
            f"dCAI optimum min/median {item['cai_minus_optimum_min']} / "
            f"{item['cai_minus_optimum_median']} "
            f"below_greedy {item['cai_below_greedy']} "
            f"below_optimum {item['cai_below_optimum']} "
            f"actionable {item['actionable_rows']} "
            f"terminations {item['terminations']}"
        )
        lines.append(f"  wall by bin {item['wall_by_length_bin_median_s']}")
        lines.append(f"  sec/kb by bin {item['sec_per_kb_by_length_bin_median']}")
    return "\n".join(lines)


def compare_rows(before_path: Path, after_path: Path) -> dict[str, Any]:
    """CAI and wall of two corpus runs, keyed by preset, host, and gene."""

    def load(path: Path) -> dict[tuple[str, str, str], dict[str, Any]]:
        rows: dict[tuple[str, str, str], dict[str, Any]] = {}
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line:
                continue
            row = json.loads(line)
            rows[(row["preset"], row["host"], row["gene"])] = row
        return rows

    before = load(before_path)
    after = load(after_path)
    regressions = []
    improvements = []
    vector_changes = []
    hash_changes = 0
    missing = 0
    ratios: dict[str, list[float]] = {}
    for key, old in before.items():
        new = after.get(key)
        if new is None:
            missing += 1
            continue
        old_cai = old.get("cai")
        new_cai = new.get("cai")
        if isinstance(old_cai, float) and isinstance(new_cai, float):
            delta = new_cai - old_cai
            if delta < -1e-12:
                regressions.append(
                    {"preset": key[0], "host": key[1], "gene": key[2], "delta": delta}
                )
            elif delta > 1e-12:
                improvements.append(
                    {"preset": key[0], "host": key[1], "gene": key[2], "delta": delta}
                )
        if old.get("vector") != new.get("vector") and old.get("vector") and new.get("vector"):
            vector_changes.append({"preset": key[0], "host": key[1], "gene": key[2]})
        old_hash = old.get("cds_sha256")
        new_hash = new.get("cds_sha256")
        if old_hash and new_hash and old_hash != new_hash:
            hash_changes += 1
        old_wall = old.get("wall_s")
        new_wall = new.get("wall_s")
        if isinstance(old_wall, float) and isinstance(new_wall, float) and old_wall > 0:
            ratios.setdefault(f"{key[0]}|{key[1]}", []).append(new_wall / old_wall)
    wall = {
        key: {
            "n": len(values),
            "median_ratio": _percentile(values, 0.50),
            "p95_ratio": _percentile(values, 0.95),
        }
        for key, values in sorted(ratios.items())
    }
    return {
        "before_rows": len(before),
        "after_rows": len(after),
        "missing_after": missing,
        "cai_regressions": regressions,
        "cai_improvements": improvements,
        "vector_changes": vector_changes,
        "cds_hash_changes": hash_changes,
        "wall_ratio_after_over_before": wall,
    }


def run_corpus(out_dir: Path, presets: list[str] | None, limit: int | None) -> None:
    sys.path.insert(0, str(ROOT / "benchmarks" / "validation"))
    import cai_independent

    from dna_design_suite.codon.bench.baselines import cai_greedy

    out_dir.mkdir(parents=True, exist_ok=True)
    sample = load_sample()
    # Keep the length span: take every k-th record of the length-sorted sample.
    if limit is not None and limit < len(sample):
        step = len(sample) / limit
        picked = []
        for index in range(limit):
            picked.append(sample[min(len(sample) - 1, int(index * step))])
        sample = picked
    ribo = _load_weights(WEIGHTS_RIBO)
    kazusa = _load_weights(WEIGHTS_KAZUSA)
    tables = {"ribosomal": ribo, "kazusa": kazusa}
    paths = {"ribosomal": str(WEIGHTS_RIBO), "kazusa": str(WEIGHTS_KAZUSA)}
    optimum: dict[tuple[str, str], float | None] = {}
    greedy: dict[tuple[str, str], float | None] = {}
    for row in sample:
        for host, weights in tables.items():
            greedy[(row["gene"], host)] = float(
                cai_independent.score(cai_greedy(row["cds"], TABLE, weights), weights)["cai"]
            )
    rows_path = out_dir / "rows.jsonl"
    done: set[tuple[str, str, str]] = set()
    if rows_path.is_file():
        for line in rows_path.read_text(encoding="utf-8").splitlines():
            if not line:
                continue
            old = json.loads(line)
            done.add((old["preset"], old["host"], old["gene"]))
            if old["preset"] == "empty" and old.get("ok") and old.get("cai") is not None:
                optimum[(old["gene"], old["host"])] = old["cai"]
    ctx = mp.get_context("fork")
    in_q: mp.Queue = ctx.Queue()
    out_q: mp.Queue = ctx.Queue()
    worker = ctx.Process(target=_serve, args=(in_q, out_q))
    worker.start()
    store_root = out_dir / "stores"
    store_root.mkdir(exist_ok=True)
    planned = [item for item in preset_table() if presets is None or item["preset"] in presets]
    # Empty first, so the constraint-free optimum exists before later presets.
    planned.sort(key=lambda item: (item["constraints"] != "empty", item["preset"]))
    with rows_path.open("a", encoding="utf-8") as handle:
        for spec in planned:
            constraints = _constraints(spec["constraints"])
            for host in spec["hosts"]:
                for row in sample:
                    key = (spec["preset"], host, row["gene"])
                    if key in done:
                        continue
                    payload = {
                        "gene": row["gene"],
                        "cds": row["cds"],
                        "weights_path": paths[host],
                        "objectives": spec["objectives"],
                        "constraints": constraints,
                        "budget": spec["budget"],
                        "generations": spec["generations"],
                        "store": str(store_root / f"{spec['preset']}-{host}-{row['gene']}"),
                    }
                    timeout = _job_timeout(row["nt"], spec["preset"])
                    in_q.put(payload)
                    try:
                        outcome = out_q.get(timeout=timeout)
                    except Exception:
                        outcome = {
                            "ok": False,
                            "wall_s": timeout,
                            "code": "TIMEOUT",
                            "detail": f"no result within {timeout:.0f}s",
                            "state": "timeout",
                        }
                        worker.terminate()
                        worker.join(5)
                        in_q = ctx.Queue()
                        out_q = ctx.Queue()
                        worker = ctx.Process(target=_serve, args=(in_q, out_q))
                        worker.start()
                    empty_ok = spec["preset"] == "empty" and outcome.get("ok")
                    if empty_ok and outcome.get("cai") is not None:
                        optimum[(row["gene"], host)] = outcome["cai"]
                    record = {
                        "preset": spec["preset"],
                        "host": host,
                        "budget": spec["budget"],
                        "generations": spec["generations"],
                        "gene": row["gene"],
                        "protein_id": row["protein_id"],
                        "nt": row["nt"],
                        "gc": row["gc"],
                        "rare_fraction": row["rare_fraction"],
                        "length_bin": row["length_bin"],
                        "cai_greedy": greedy[(row["gene"], host)],
                        "cai_optimum": optimum.get((row["gene"], host)),
                        **{k: v for k, v in outcome.items() if k != "cds"},
                        "cds_sha256": hashlib.sha256(outcome["cds"].encode()).hexdigest()
                        if outcome.get("cds")
                        else None,
                    }
                    handle.write(json.dumps(record) + "\n")
                    handle.flush()
                    done.add(key)
    in_q.put(None)
    worker.join(30)
    if worker.is_alive():
        worker.terminate()
    report = summarize(rows_path)
    (out_dir / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(_format_report(report))


def main() -> None:
    parser = argparse.ArgumentParser(description="Broad codon-design corpus")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("sample")
    run = sub.add_parser("run")
    run.add_argument("--out", type=Path, required=True)
    run.add_argument("--presets", default="")
    run.add_argument("--limit", type=int, default=0)
    show = sub.add_parser("summarize")
    show.add_argument("rows", type=Path)
    diff = sub.add_parser("compare")
    diff.add_argument("before", type=Path)
    diff.add_argument("after", type=Path)
    args = parser.parse_args()
    if args.command == "sample":
        weights = _load_weights(WEIGHTS_RIBO)
        rows = draw_sample(load_records(), weights)
        write_sample(rows)
        bins: dict[str, int] = {}
        for row in rows:
            bins[row["length_bin"]] = bins.get(row["length_bin"], 0) + 1
        print(
            f"wrote {SAMPLE_PATH} n={len(rows)} "
            f"nt {rows[0]['nt']}..{rows[-1]['nt']} bins {bins} "
            f"sha256 {sample_digest(rows)}"
        )
        return
    if args.command == "summarize":
        print(_format_report(summarize(args.rows)))
        return
    if args.command == "compare":
        print(json.dumps(compare_rows(args.before, args.after), indent=2))
        return
    presets = [item for item in args.presets.split(",") if item] or None
    run_corpus(args.out, presets, args.limit or None)


if __name__ == "__main__":
    main()
