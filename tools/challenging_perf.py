#!/usr/bin/env python3
"""Wall time and peak RSS on long, constrained, and multi-objective jobs.

Seed is 7 throughout.

Sweep (``--suite sweep``): budget ``fast``.
In-core jobs use generations 0, so only tier 1 and tier 2 run. That is the
fast plan. A structure objective is L0, so those jobs use generations 2 and
stagnation window 21. Otherwise the axis is never scored.

Standard rerun (``--suite standard``): the slowest design jobs from a sweep
file, budget ``standard``, generations 4, stagnation window 21.

Each case is a child process. ``ru_maxrss`` is that process, including the
interpreter and the loaded extension. Linux reports it in kilobytes.

From the repository root:

    .venv/bin/python tools/challenging_perf.py --suite sweep --out /tmp/cf-sweep.json
    .venv/bin/python tools/challenging_perf.py --suite standard \\
        --from /tmp/cf-sweep.json --out /tmp/cf-standard.json
"""

from __future__ import annotations

import argparse
import json
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
SYNTHETIC_FLANK = "ACGT" * 75


def _cases() -> list[dict]:
    """The sweep. Names are stable; ``--suite standard`` replays design rows."""
    assembly = "benchmarks/validation/constraints_assembly.json"
    golden = "examples/golden_gate.json"
    biobrick = "examples/biobrick_rfc10.json"
    rows: list[dict] = []

    def design(
        name: str,
        gene: str,
        objectives: list[str],
        constraints: str | None,
        generations: int,
    ) -> None:
        rows.append(
            {
                "kind": "design",
                "name": name,
                "gene": gene,
                "objectives": objectives,
                "constraints": constraints,
                "generations": generations,
                "budget": "fast",
            }
        )

    for gene in ("rpoB", "lacZ", "yeeJ"):
        design(f"{gene} CAI empty", gene, ["usage.cai"], None, 0)
    for gene in ("aceB", "thrA"):
        design(f"{gene} CAI ribosomal AGG", gene, ["usage.cai"], None, 0)
    for gene in ("yagF", "waaS", "elp", "yeeJ", "rpoB"):
        design(f"{gene} CAI assembly", gene, ["usage.cai"], assembly, 0)
    design("rpoB CAI Golden Gate", "rpoB", ["usage.cai"], golden, 0)
    design("rpoB CAI BioBrick", "rpoB", ["usage.cai"], biobrick, 0)
    design("waaS CAI+burden assembly", "waaS", ["usage.cai", "burden.growth_cost"], assembly, 0)
    design("rpoB CAI+burden empty", "rpoB", ["usage.cai", "burden.growth_cost"], None, 0)
    design("aceB CAI+burden empty", "aceB", ["usage.cai", "burden.growth_cost"], None, 0)
    structure = ["usage.cai", "structure.local_accessibility"]
    for gene in ("rpsL", "waaS", "elp", "lacZ", "yeeJ"):
        design(f"{gene} CAI+structure assembly", gene, structure, assembly, 2)
    design(
        "waaS CAI+structure+burden assembly",
        "waaS",
        ["usage.cai", "structure.local_accessibility", "burden.growth_cost"],
        assembly,
        2,
    )
    rows.append(
        {
            "kind": "operon",
            "name": "operon rpsL+waaS",
            "genes": ["rpsL", "waaS"],
            "fluxes": [1.0, 2.0],
            "passes": 2,
            "generations": 0,
            "budget": "fast",
        }
    )
    rows.append(
        {
            "kind": "operon",
            "name": "operon lacZ+thrA",
            "genes": ["lacZ", "thrA"],
            "fluxes": [1.0, 2.0],
            "passes": 2,
            "generations": 0,
            "budget": "fast",
        }
    )
    for gene, synthetic in (("yeeJ", False), ("rpoB", False), ("yeeJ", True)):
        label = "synthetic flanks" if synthetic else "CDS flanks"
        rows.append(
            {
                "kind": "primer",
                "name": f"primers {gene} {label}",
                "gene": gene,
                "synthetic_flanks": synthetic,
            }
        )
    return rows


def _fasta(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    return "".join(line.strip() for line in lines if line and not line.startswith(">"))


def _sequences() -> dict[str, str]:
    from codonforge.bench.corpus import load_panel, select_genes

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
        return {
            "schema_version": "1.0.0",
            "constraint_set_id": "cst_empty",
            "constraints": [],
        }
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def _profile() -> dict:
    path = ROOT / "examples" / "mg1655_ribosomal_usage.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _rss_kb() -> int:
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def _run_design(spec: dict) -> dict:
    from codonforge import _core
    from codonforge.modelreg import ArtifactStore
    from codonforge.orchestrator.design import DesignRequest, design

    sequences = _sequences()
    cds = sequences[spec["gene"]]
    constraints = _constraints(spec["constraints"])
    request = DesignRequest(
        host=HOST,
        cds=cds,
        host_profile=_profile(),
        objectives=list(spec["objectives"]),
        constraints=constraints,
        translation_table=11,
        budget=spec["budget"],
        seed=SEED,
        generations=int(spec["generations"]),
        stagnation_window=21,
        segment=2,
        host_digest=DIGEST,
    )
    started = time.perf_counter()
    store = ArtifactStore(Path("/tmp") / f"cf-stress-{spec['name'].replace(' ', '_')}")
    result = design(request, store)
    wall = time.perf_counter() - started
    status = result.status
    tiers = [
        {
            "tier": tier["tier"],
            "state": tier["state"],
            "generation": tier.get("generation"),
        }
        for tier in status["progress"]["tiers"]
    ]
    payload: dict = {
        "state": result.job.state,
        "wall_s": round(wall, 3),
        "peak_rss_kb": _rss_kb(),
        "nt": len(cds),
        "tiers": tiers,
        "warnings": sorted({warning["code"] for warning in status.get("warnings", [])}),
    }
    if result.job.problem:
        payload["problem"] = result.job.problem.get("code")
    if result.front and result.front.get("candidates"):
        chosen = result.front["candidates"][0]
        designed = chosen["sequences"][0]["cds"]
        report = json.loads(
            _core.constraint_violations(
                json.dumps(constraints),
                designed,
                [(0, len(designed))],
                "assembled",
            )
        )
        actionable = sum(1 for item in report["violations"] if item.get("actionable"))
        payload.update(
            {
                "vector": chosen["objective_vector"],
                "rule": result.front["recommended"]["rule"],
                "feasible": (
                    chosen["feasible"]
                    if chosen.get("feasible") is not None
                    else result.job.is_feasible(designed)
                ),
                "actionable_violations": actionable,
                "protein_identical": _core.translate(designed, 11) == _core.translate(cds, 11),
                "in_frame_agg": designed.count("AGG") and sum(
                    1
                    for index in range(0, len(designed), 3)
                    if designed[index : index + 3] == "AGG"
                ),
                "n_candidates": len(result.front["candidates"]),
            }
        )
    session = result.job._fold_session
    if session is not None:
        payload["folds"] = {
            "hits": session.hits,
            "full_folds": session.full_folds,
            "local_edits": session.local_edits,
        }
    return payload


def _run_operon(spec: dict) -> dict:
    from codonforge import _core
    from codonforge.modelreg import ArtifactStore
    from codonforge.orchestrator import Cistron, design_operon

    sequences = _sequences()
    cistrons = [
        Cistron(sequences[gene], float(flux), gene)
        for gene, flux in zip(spec["genes"], spec["fluxes"], strict=True)
    ]
    started = time.perf_counter()
    result = design_operon(
        cistrons,
        ArtifactStore(Path("/tmp") / f"cf-stress-{spec['name'].replace(' ', '_')}"),
        host=HOST,
        host_profile=_profile(),
        objectives=["burden.growth_cost"],
        translation_table=11,
        budget=spec["budget"],
        seed=SEED,
        generations=int(spec["generations"]),
        passes=int(spec["passes"]),
        host_digest=DIGEST,
    )
    wall = time.perf_counter() - started
    return {
        "state": "succeeded",
        "wall_s": round(wall, 3),
        "peak_rss_kb": _rss_kb(),
        "nt": [len(sequence) for sequence in result.sequences],
        "joint_growth_cost": result.joint_growth_cost,
        "independent_joint_growth_cost": result.independent_joint_growth_cost,
        "solo_growth_costs": result.solo_growth_costs,
        "swaps_per_pass": result.swaps_per_pass,
        "protein_identical": [
            _core.translate(designed, 11) == protein
            for designed, protein in zip(result.sequences, result.proteins, strict=True)
        ],
        "cistron_states": [job.job.state for job in result.jobs],
    }


def _run_primer(spec: dict) -> dict:
    from codonforge import _core

    sequence = _sequences()[spec["gene"]]
    flank = 300
    if spec["synthetic_flanks"]:
        template = SYNTHETIC_FLANK + sequence + SYNTHETIC_FLANK
        start = len(SYNTHETIC_FLANK)
        end = start + len(sequence)
        flank_note = "synthetic ACGT repeat, not a biological flank"
    else:
        template = sequence
        start = flank
        end = len(sequence) - flank
        flank_note = "flanks are the same CDS"
    region = end - start
    body = {
        "sequence": template,
        "region": {"start": start, "end": end},
        "config": {
            "tm_target": 60.0,
            "product_size_range": {"min": region, "max": region + 2 * flank},
            "flank_length": flank,
            "primer_length_range": {"min": 18, "max": 25},
        },
    }
    started = time.perf_counter()
    try:
        raw = _core.primer_search(json.dumps(body), lambda _fraction: None)
    except ValueError as exc:
        return {
            "state": "failed",
            "wall_s": round(time.perf_counter() - started, 3),
            "peak_rss_kb": _rss_kb(),
            "nt": len(sequence),
            "region_nt": region,
            "flanks": flank_note,
            "problem": str(exc).split(":", 1)[0],
            "detail": str(exc),
        }
    wall = time.perf_counter() - started
    result = json.loads(raw)
    pairs = result.get("pairs") or []
    recommended = pairs[0] if pairs else {}
    return {
        "state": "succeeded",
        "wall_s": round(wall, 3),
        "peak_rss_kb": _rss_kb(),
        "nt": len(sequence),
        "template_nt": len(template),
        "region_nt": region,
        "flanks": flank_note,
        "n_pairs": len(pairs),
        "survivors": result.get("survivors"),
        "truncated": result.get("truncated"),
        "product_size": recommended.get("product_size"),
        "tm_forward": recommended.get("tm_forward"),
        "tm_reverse": recommended.get("tm_reverse"),
        "hetero_dimer_dG": recommended.get("hetero_dimer_dG"),
        "refused": recommended.get("refused"),
    }


def run_spec(spec: dict) -> dict:
    try:
        if spec["kind"] == "design":
            measured = _run_design(spec)
        elif spec["kind"] == "operon":
            measured = _run_operon(spec)
        elif spec["kind"] == "primer":
            measured = _run_primer(spec)
        else:
            raise ValueError(spec["kind"])
    except Exception as exc:  # a stress row should be recorded, not abort the suite
        measured = {
            "state": "error",
            "problem": type(exc).__name__,
            "detail": str(exc),
            "trace": traceback.format_exc(limit=8),
            "peak_rss_kb": _rss_kb(),
        }
    return {"name": spec["name"], "kind": spec["kind"], "spec": spec, **measured}


def _worker() -> int:
    spec = json.loads(sys.stdin.read())
    json.dump(run_spec(spec), sys.stdout)
    return 0


def _one(spec: dict, timeout: int) -> dict:
    try:
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--worker"],
            input=json.dumps(spec),
            text=True,
            capture_output=True,
            timeout=timeout,
            cwd=ROOT,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return {
            "name": spec["name"],
            "kind": spec["kind"],
            "spec": spec,
            "state": "timeout",
            "wall_s": timeout,
            "timeout_s": timeout,
        }
    if completed.returncode != 0 or not completed.stdout.strip():
        return {
            "name": spec["name"],
            "kind": spec["kind"],
            "spec": spec,
            "state": "error",
            "detail": (completed.stderr or completed.stdout)[-2000:],
            "returncode": completed.returncode,
        }
    return json.loads(completed.stdout)


def _standard_specs(source: Path, limit: int) -> list[dict]:
    rows = json.loads(source.read_text(encoding="utf-8"))["rows"]
    designs = [
        row
        for row in rows
        if row.get("kind") == "design" and isinstance(row.get("wall_s"), (int, float))
    ]
    designs.sort(key=lambda row: row["wall_s"], reverse=True)
    specs = []
    for row in designs[:limit]:
        spec = dict(row["spec"])
        spec["budget"] = "standard"
        spec["generations"] = 4
        spec["name"] = row["name"] + " standard"
        specs.append(spec)
    return specs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--suite", choices=("sweep", "standard"), default="sweep")
    parser.add_argument("--from", dest="source", type=Path)
    parser.add_argument("--out", type=Path, required=False)
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument("--limit", type=int, default=4)
    args = parser.parse_args(argv)
    if args.worker:
        return _worker()
    if args.out is None:
        parser.error("--out is required")
    if args.suite == "standard":
        if args.source is None:
            parser.error("--suite standard needs --from")
        specs = _standard_specs(args.source, args.limit)
    else:
        specs = _cases()
    rows = []
    for spec in specs:
        print(spec["name"], flush=True)
        row = _one(spec, args.timeout)
        rows.append(row)
        print(f"  {row.get('state')} {row.get('wall_s')} s", flush=True)
    document = {
        "seed": SEED,
        "host": HOST,
        "suite": args.suite,
        "timeout_s": args.timeout,
        "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
