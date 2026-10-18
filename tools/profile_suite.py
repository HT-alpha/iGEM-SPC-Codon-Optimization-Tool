#!/usr/bin/env python3
"""Attribute DNA Design Suite and Primer design time, then record peak RSS.

Seed is 7. Host is ``escherichia_coli:511145/BL21``. Codon weights are
``examples/mg1655_ribosomal_usage.json``. This is the v2.9.1 challenging
sweep plus Primer design flank and region rows, stage timings, and the Rust
counters behind ``CF_SEED_PROFILE`` and ``CF_PRIMER_PROFILE``.

Each row is its own process. ``ru_maxrss`` is that process. Linux reports
kilobytes. cProfile rows are cumulative time inside that process, including
the one PyO3 call that owns a Rust stage.

From the repository root, after ``make build``:

    .venv/bin/python tools/profile_suite.py --suite stages --out /tmp/cf-stages.json
    .venv/bin/python tools/profile_suite.py --suite seed --out /tmp/cf-seed.json
    .venv/bin/python tools/profile_suite.py --suite primers --out /tmp/cf-primers.json
    .venv/bin/python tools/profile_suite.py --suite structure --out /tmp/cf-structure.json
    .venv/bin/python tools/profile_suite.py --suite operon --out /tmp/cf-operon.json
    .venv/bin/python tools/profile_suite.py --suite cprofile --out /tmp/cf-cprofile.json
    .venv/bin/python tools/profile_suite.py --suite potts --out /tmp/cf-potts.json

``--suite all`` runs those in order.
"""

from __future__ import annotations

import argparse
import cProfile
import io
import json
import os
import pstats
import resource
import subprocess
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SEED = 7
HOST = "escherichia_coli:511145/BL21"
DIGEST = "sha256:" + "ab" * 32
ASSEMBLY = "benchmarks/validation/constraints_assembly.json"


def _fasta(path: Path) -> str:
    return "".join(
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith(">")
    )


def _sequences() -> dict[str, str]:
    from dna_design_suite.codon.bench.corpus import load_panel, select_genes

    found = {gene["gene"]: gene["cds"] for gene in select_genes(load_panel())}
    validation = ROOT / "benchmarks" / "validation" / "sequences"
    for name, filename in (
        ("yagF", "yagF.fa"),
        ("waaS", "waaS.fa"),
        ("elp", "elp_vpgvg.fa"),
        ("yeeJ", "yeeJ.fa"),
        ("rpsL", "rpsL.fa"),
    ):
        found[name] = _fasta(validation / filename)
    return found


def _constraints(name: str | None) -> dict:
    if name is None:
        return {"schema_version": "1.0.0", "constraint_set_id": "cst_empty", "constraints": []}
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def _profile() -> dict:
    return json.loads(
        (ROOT / "examples" / "mg1655_ribosomal_usage.json").read_text(encoding="utf-8")
    )


def _rss_kb() -> int:
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def _wrap_core() -> dict[str, dict[str, float | int]]:
    from dna_design_suite import _core

    totals: dict[str, dict[str, float | int]] = {}

    def wrap(name: str, fn):
        def inner(*args, **kwargs):
            started = time.perf_counter()
            try:
                return fn(*args, **kwargs)
            finally:
                bucket = totals.setdefault(name, {"seconds": 0.0, "calls": 0})
                bucket["seconds"] += time.perf_counter() - started
                bucket["calls"] += 1

        return inner

    for name in (
        "search_tier1_seed",
        "resolve_feasible_seed",
        "search_tier2_polish",
        "search_tier3_segment",
        "explain_baseline",
        "term_scores",
        "constraint_violations",
        "constraint_unsearched",
        "primer_search",
    ):
        if hasattr(_core, name):
            setattr(_core, name, wrap(name, getattr(_core, name)))
    return totals


def _design(gene: str, objectives: list[str], constraints: str | None, generations: int) -> dict:
    from dna_design_suite.codon.orchestrator.design import DesignRequest, design
    from dna_design_suite.modelreg import ArtifactStore

    totals = _wrap_core()
    sequences = _sequences()
    cds = sequences[gene]
    started = time.perf_counter()
    result = design(
        DesignRequest(
            host=HOST,
            cds=cds,
            host_profile=_profile(),
            objectives=objectives,
            constraints=_constraints(constraints),
            translation_table=11,
            budget="fast",
            seed=SEED,
            generations=generations,
            stagnation_window=21,
            segment=2,
            host_digest=DIGEST,
        ),
        ArtifactStore(Path("/tmp") / f"cf-profile-{gene}-{generations}"),
    )
    wall = time.perf_counter() - started
    payload: dict = {
        "state": result.job.state,
        "wall_s": wall,
        "peak_rss_kb": _rss_kb(),
        "nt": len(cds),
        "stages": {
            name: {"seconds": round(item["seconds"], 6), "calls": item["calls"]}
            for name, item in totals.items()
        },
    }
    if result.front and result.front.get("candidates"):
        chosen = result.front["candidates"][0]
        payload["vector"] = chosen["objective_vector"]
        payload["front_bytes"] = len(json.dumps(result.front).encode())
        payload["status_bytes"] = len(json.dumps(result.status).encode())
        designed = chosen["sequences"][0]["cds"]
        from dna_design_suite import _core

        payload["protein_identical"] = _core.translate(designed, 11) == _core.translate(cds, 11)
    session = result.job._fold_session
    if session is not None:
        payload["folds"] = {
            "hits": session.hits,
            "full_folds": session.full_folds,
            "local_edits": session.local_edits,
        }
    return payload


def _seed_profile(gene: str) -> dict:
    from dna_design_suite import _core

    sequences = _sequences()
    cds = sequences[gene]
    request = {
        "cds": cds,
        "translation_table": 11,
        "five_prime_flank": "",
        "three_prime_flank": "",
        "constraint_set": _constraints(ASSEMBLY),
        "host_profile": _profile(),
        "host_id": HOST,
        "host_digest": DIGEST,
        "terms": ["usage.cai"],
        "weights": [1.0],
        "root_seed": SEED,
    }
    started = time.perf_counter()
    report = json.loads(_core.resolve_feasible_seed(json.dumps(request)))
    return {
        "wall_s": time.perf_counter() - started,
        "peak_rss_kb": _rss_kb(),
        "nt": len(cds),
        "feasible": report["feasible"],
        "edits": report["edits"],
        "cycles": report["cycles"],
        "cds_sha256": __import__("hashlib").sha256(report["cds"].encode()).hexdigest(),
    }


def _primer(
    sequence: str, start: int, end: int, flank: int, product_min: int | None = None
) -> dict:
    from dna_design_suite import _core

    region = end - start
    body = {
        "sequence": sequence,
        "region": {"start": start, "end": end},
        "config": {
            "tm_target": 60.0,
            "product_size_range": {
                "min": region if product_min is None else product_min,
                "max": region + 2 * flank,
            },
            "flank_length": flank,
            "primer_length_range": {"min": 18, "max": 25},
        },
    }
    started = time.perf_counter()
    raw = _core.primer_search(json.dumps(body), lambda _fraction: None)
    wall = time.perf_counter() - started
    result = json.loads(raw)
    recommended = (result.get("pairs") or [{}])[0]
    return {
        "wall_s": wall,
        "peak_rss_kb": _rss_kb(),
        "region_nt": region,
        "flank": flank,
        "response_bytes": len(raw.encode()),
        "survivors": result.get("survivors"),
        "truncated": result.get("truncated"),
        "n_pairs": len(result.get("pairs") or []),
        "product_size": recommended.get("product_size"),
        "hetero_dimer_dG": recommended.get("hetero_dimer_dG"),
        "refused": recommended.get("refused"),
        "sha256": __import__("hashlib").sha256(raw.encode()).hexdigest(),
    }


def _operon(genes: list[str], fluxes: list[float]) -> dict:
    from dna_design_suite.codon.orchestrator import Cistron, design_operon
    from dna_design_suite.modelreg import ArtifactStore

    sequences = _sequences()
    cistrons = [
        Cistron(sequences[gene], float(flux), gene)
        for gene, flux in zip(genes, fluxes, strict=True)
    ]
    started = time.perf_counter()
    result = design_operon(
        cistrons,
        ArtifactStore(Path("/tmp") / ("cf-profile-operon-" + "-".join(genes))),
        host=HOST,
        host_profile=_profile(),
        objectives=["burden.growth_cost"],
        translation_table=11,
        budget="fast",
        seed=SEED,
        generations=0,
        passes=2,
        host_digest=DIGEST,
    )
    return {
        "wall_s": time.perf_counter() - started,
        "peak_rss_kb": _rss_kb(),
        "nt": [len(sequence) for sequence in result.sequences],
        "joint_growth_cost": result.joint_growth_cost,
        "swaps_per_pass": result.swaps_per_pass,
    }


def _cprofile(kind: str) -> dict:
    profiler = cProfile.Profile()
    if kind == "yeeJ-assembly":
        profiler.enable()
        measured = _design("yeeJ", ["usage.cai"], ASSEMBLY, 0)
        profiler.disable()
    elif kind == "primer-yeeJ":
        sequence = _sequences()["yeeJ"]
        profiler.enable()
        measured = _primer(sequence, 300, len(sequence) - 300, 300)
        profiler.disable()
    else:
        raise ValueError(kind)
    stream = io.StringIO()
    stats = pstats.Stats(profiler, stream=stream).sort_stats("cumulative")
    stats.print_stats(25)
    rows = []
    for func, (calls, _reccalls, _excl, cum, _callers) in stats.stats.items():
        filename, line, name = func
        rows.append(
            {
                "function": name,
                "file": filename,
                "line": line,
                "calls": calls,
                "cumulative_s": cum,
            }
        )
    rows.sort(key=lambda item: item["cumulative_s"], reverse=True)
    return {
        "kind": kind,
        "measured": {key: value for key, value in measured.items() if key != "stages"},
        "top": rows[:25],
        "text": stream.getvalue(),
    }


def _potts(gene: str) -> dict:
    return _design(gene, ["usage.cai", "epistasis.context_naturalness"], None, 0)


def _run_job(job: dict) -> dict:
    kind = job["kind"]
    if kind == "design":
        measured = _design(job["gene"], job["objectives"], job["constraints"], job["generations"])
    elif kind == "seed":
        measured = _seed_profile(job["gene"])
    elif kind == "primer":
        measured = _primer(
            job["sequence"], job["start"], job["end"], job["flank"], job.get("product_min")
        )
    elif kind == "operon":
        measured = _operon(job["genes"], job["fluxes"])
    elif kind == "cprofile":
        measured = _cprofile(job["name"])
    elif kind == "potts":
        measured = _potts(job["gene"])
    else:
        raise ValueError(kind)
    return measured


def _worker() -> int:
    job = json.loads(sys.stdin.read())
    try:
        measured = _run_job(job)
        measured["state"] = measured.get("state", "succeeded")
    except Exception as exc:
        measured = {
            "state": "error",
            "problem": type(exc).__name__,
            "detail": str(exc),
            "trace": traceback.format_exc(limit=12),
            "peak_rss_kb": _rss_kb(),
        }
    json.dump(measured, sys.stdout)
    return 0


def _jobs(suite: str) -> list[dict]:
    sequences = None

    def sequence(gene: str) -> str:
        nonlocal sequences
        if sequences is None:
            sequences = _sequences()
        return sequences[gene]

    rows: list[dict] = []
    if suite in {"stages", "all"}:
        for gene in ("yeeJ", "rpoB", "lacZ", "waaS", "elp", "rpsL", "yagF", "aceB", "thrA"):
            rows.append(
                {
                    "suite": "stages",
                    "name": f"{gene} CAI empty",
                    "kind": "design",
                    "gene": gene,
                    "objectives": ["usage.cai"],
                    "constraints": None,
                    "generations": 0,
                }
            )
        for gene in ("yeeJ", "rpoB", "lacZ", "waaS", "elp", "rpsL", "yagF"):
            rows.append(
                {
                    "suite": "stages",
                    "name": f"{gene} CAI assembly",
                    "kind": "design",
                    "gene": gene,
                    "objectives": ["usage.cai"],
                    "constraints": ASSEMBLY,
                    "generations": 0,
                }
            )
        rows.append(
            {
                "suite": "stages",
                "name": "rpoB CAI+burden empty",
                "kind": "design",
                "gene": "rpoB",
                "objectives": ["usage.cai", "burden.growth_cost"],
                "constraints": None,
                "generations": 0,
            }
        )
        rows.append(
            {
                "suite": "stages",
                "name": "waaS CAI+burden assembly",
                "kind": "design",
                "gene": "waaS",
                "objectives": ["usage.cai", "burden.growth_cost"],
                "constraints": ASSEMBLY,
                "generations": 0,
            }
        )
    if suite in {"seed", "all"}:
        for gene in ("yeeJ", "rpoB", "waaS", "elp"):
            rows.append({"suite": "seed", "name": f"seed {gene}", "kind": "seed", "gene": gene})
    if suite in {"primers", "all"}:
        yeej = sequence("yeeJ")
        rpsl = sequence("rpsL")
        mid = len(yeej) // 2
        for flank in (50, 100, 200, 300):
            start, end = mid - 500, mid + 500
            rows.append(
                {
                    "suite": "primers",
                    "name": f"yeeJ region 1000 flank {flank}",
                    "kind": "primer",
                    "sequence": yeej,
                    "start": start,
                    "end": end,
                    "flank": flank,
                    "product_min": max(1, (end - start) - 2 * flank),
                }
            )
        for region in (300, 1000, len(yeej) - 600):
            start = 300
            end = start + region
            rows.append(
                {
                    "suite": "primers",
                    "name": f"yeeJ region {region} flank 300",
                    "kind": "primer",
                    "sequence": yeej,
                    "start": start,
                    "end": end,
                    "flank": 300,
                }
            )
        rows.append(
            {
                "suite": "primers",
                "name": "rpsL whole CDS flank 40",
                "kind": "primer",
                "sequence": rpsl,
                "start": 40,
                "end": len(rpsl) - 40,
                "flank": 40,
            }
        )
        synthetic = "ACGT" * 75
        template = synthetic + yeej + synthetic
        rows.append(
            {
                "suite": "primers",
                "name": "yeeJ synthetic flanks",
                "kind": "primer",
                "sequence": template,
                "start": len(synthetic),
                "end": len(synthetic) + len(yeej),
                "flank": 300,
            }
        )
    if suite in {"structure", "all"}:
        for gene, generations in (("rpsL", 2), ("waaS", 2), ("elp", 2), ("lacZ", 2), ("yeeJ", 2)):
            rows.append(
                {
                    "suite": "structure",
                    "name": f"{gene} CAI+structure assembly",
                    "kind": "design",
                    "gene": gene,
                    "objectives": ["usage.cai", "structure.local_accessibility"],
                    "constraints": ASSEMBLY,
                    "generations": generations,
                }
            )
    if suite in {"operon", "all"}:
        rows.append(
            {
                "suite": "operon",
                "name": "operon rpsL+waaS",
                "kind": "operon",
                "genes": ["rpsL", "waaS"],
                "fluxes": [1.0, 2.0],
            }
        )
        rows.append(
            {
                "suite": "operon",
                "name": "operon lacZ+thrA",
                "kind": "operon",
                "genes": ["lacZ", "thrA"],
                "fluxes": [1.0, 2.0],
            }
        )
    if suite in {"cprofile", "all"}:
        rows.append({"suite": "cprofile", "name": "yeeJ-assembly", "kind": "cprofile"})
        rows.append({"suite": "cprofile", "name": "primer-yeeJ", "kind": "cprofile"})
    if suite in {"potts", "all"}:
        for gene in ("rpsL", "waaS", "lacZ"):
            rows.append(
                {"suite": "potts", "name": f"{gene} CAI+Potts", "kind": "potts", "gene": gene}
            )
    return rows


def _child_env(job: dict) -> dict[str, str]:
    env = os.environ.copy()
    if job["kind"] == "seed":
        env["CF_SEED_PROFILE"] = "1"
    if job["kind"] == "primer":
        env["CF_PRIMER_PROFILE"] = "1"
    if job["kind"] == "design" and "structure.local_accessibility" in job.get("objectives", []):
        env["CF_TIER3_PROFILE"] = "1"
    return env


def _one(job: dict, timeout: int) -> dict:
    # The sequence is large. Pass the job on stdin.
    payload = dict(job)
    payload.pop("sequence", None)
    if job["kind"] == "primer":
        payload["sequence"] = job["sequence"]
    try:
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--worker"],
            input=json.dumps(payload),
            text=True,
            capture_output=True,
            timeout=timeout,
            cwd=ROOT,
            check=False,
            env=_child_env(job),
        )
    except subprocess.TimeoutExpired:
        return {"name": job["name"], "state": "timeout", "wall_s": timeout}
    stderr_profile = [
        line
        for line in completed.stderr.splitlines()
        if line.startswith(("CF_SEED_PROFILE ", "CF_PRIMER_PROFILE ", "CF_TIER3_PROFILE "))
    ]
    if completed.returncode != 0 or not completed.stdout.strip():
        return {
            "name": job["name"],
            "state": "error",
            "detail": (completed.stderr or completed.stdout)[-2000:],
            "returncode": completed.returncode,
            "profile_lines": stderr_profile,
        }
    measured = json.loads(completed.stdout)
    measured["name"] = job["name"]
    measured["suite"] = job["suite"]
    if stderr_profile:
        measured["rust_profile"] = stderr_profile[-1].split(" ", 1)[1]
    return measured


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument(
        "--suite",
        choices=("stages", "seed", "primers", "structure", "operon", "cprofile", "potts", "all"),
        default="all",
    )
    parser.add_argument("--out", type=Path)
    parser.add_argument("--timeout", type=int, default=240)
    args = parser.parse_args(argv)
    if args.worker:
        return _worker()
    if args.out is None:
        parser.error("--out is required")
    rows = []
    for job in _jobs(args.suite):
        print(job["name"], flush=True)
        row = _one(job, args.timeout)
        rows.append(row)
        print(f"  {row.get('state')} {row.get('wall_s')} s", flush=True)
    document = {"seed": SEED, "host": HOST, "suite": args.suite, "rows": rows}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
