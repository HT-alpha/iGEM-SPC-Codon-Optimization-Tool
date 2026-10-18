#!/usr/bin/env python3
"""Harsher wall time, phase split, peak RSS, and output pins.

Heavier than ``tools/challenging_perf.py`` and ``tools/profile_suite.py``:
the next-longest real CDS after yeeJ, a tighter assembly set, CAI with
structure or burden or Potts, a long operon, Primer design on long templates
and wide oligo windows, and a few concurrent API jobs.

Seed is 7. Host is ``escherichia_coli:511145/BL21``. Codon weights are
``examples/mg1655_ribosomal_usage.json``. Translation table 11. Each row is
its own process, so ``ru_maxrss`` is that process. Linux reports kilobytes.
``CF_SEED_PROFILE``, ``CF_PRIMER_PROFILE`` and ``CF_TIER3_PROFILE`` are the
phase counters.

yeeJ (7077 nt) is the longest CDS in ``benchmarks/corpus``. yfhM (4962) and
rpoC (4224) are the next two, taken from that same file.

From the repository root, after ``make build``:

    .venv/bin/python tools/harsh_perf.py --suite all --out /tmp/cf-harsh.json

``--suite`` is one of ``pins``, ``design``, ``seed``, ``objectives``,
``operon``, ``primers``, ``api``, ``gate``, ``all``.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import resource
import subprocess
import sys
import time
import traceback
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmarks" / "validation"))

import cai_independent  # noqa: E402

SEED = 7
HOST = "escherichia_coli:511145/BL21"
DIGEST = "sha256:" + "ab" * 32
ASSEMBLY = "assembly"
GOLDEN = "golden_gate"
BIOBRICK = "biobrick_rfc10"
HARSH = "harsh_assembly"

# v2.9.4 search-escape table, budget fast, generations 0, seed 7.
# A later row may be higher. It may not be lower.
ASSEMBLY_FLOOR = {
    "rpsL": 0.9866866273575934,
    "elastin_vpgvg20": 0.4642192192077335,
    "gfp_BBa_E0040": 0.9725744807921104,
    "mrfp_BBa_E1010": 0.9605695971330996,
    "mcherry_BBa_J06504": 0.9644182772231985,
    "bla": 0.9847528752532632,
    "waaS": 0.9259163095438141,
    "yagF": 0.9453398274862173,
    "yeeJ": 0.9727484991594185,
}

# The 300 nt yeeJ flank region. Byte-identical with test_profile_identity.
PRIMER_PIN = {
    "forward_sequence": "TGACCGACAAATATTATG",
    "reverse_sequence": "GAACTACAGTGGTTTAAA",
    "product_size": 6843,
    "hetero_dimer_dG": -1.283850000000001,
}

FASTA = {
    "rpsL": "benchmarks/validation/sequences/rpsL.fa",
    "elastin_vpgvg20": "benchmarks/validation/sequences/elp_vpgvg.fa",
    "gfp_BBa_E0040": "benchmarks/validation/sequences/BBa_E0040_onestop.fa",
    "mrfp_BBa_E1010": "benchmarks/validation/sequences/BBa_E1010.fa",
    "mcherry_BBa_J06504": "benchmarks/validation/sequences/BBa_J06504.fa",
    "bla": "benchmarks/validation/sequences/bla.fa",
    "waaS": "benchmarks/validation/sequences/waaS.fa",
    "yagF": "benchmarks/validation/sequences/yagF.fa",
    "lacZ": "benchmarks/validation/sequences/lacZ.fa",
    "yeeJ": "benchmarks/validation/sequences/yeeJ.fa",
}


def _fasta(path: Path) -> str:
    sequence = (
        "".join(
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith(">")
        )
        .upper()
        .replace("U", "T")
    )
    if not sequence or len(sequence) % 3 or any(base not in "ACGT" for base in sequence):
        raise ValueError(f"{path} is not an ACGT CDS")
    return sequence


def _corpus_gene(gene: str) -> str:
    found = ""
    current = None
    buf: list[str] = []
    with gzip.open(ROOT / "benchmarks" / "corpus" / "NC_000913.3_cds.fna.gz", "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if current == gene:
                    found = "".join(buf)
                    break
                current = None
                if "[gene=" in line:
                    current = line.split("[gene=", 1)[1].split("]", 1)[0]
                buf = []
            else:
                buf.append(line.strip())
        else:
            if current == gene:
                found = "".join(buf)
    if not found:
        raise ValueError(f"{gene} is not in the corpus CDS file")
    return found.upper()


def _sequences() -> dict[str, str]:
    found = {name: _fasta(ROOT / relative) for name, relative in FASTA.items()}
    for gene in ("rpoB", "rpoC", "yfhM", "thrA", "aceB"):
        found[gene] = _corpus_gene(gene)
    return found


def _profile() -> dict:
    return json.loads(
        (ROOT / "examples" / "mg1655_ribosomal_usage.json").read_text(encoding="utf-8")
    )


def _weights() -> dict[str, float]:
    return cai_independent.load_ribosomal_weights(ROOT / "examples" / "mg1655_ribosomal_usage.json")


def _constraints(name: str | None) -> dict:
    if name is None:
        return {
            "schema_version": "1.0.0",
            "constraint_set_id": "cst_empty",
            "constraints": [],
        }
    files = {
        ASSEMBLY: ROOT / "benchmarks" / "validation" / "constraints_assembly.json",
        GOLDEN: ROOT / "examples" / "golden_gate.json",
        BIOBRICK: ROOT / "examples" / "biobrick_rfc10.json",
    }
    if name == HARSH:
        document = json.loads(files[ASSEMBLY].read_text(encoding="utf-8"))
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
    return json.loads(files[name].read_text(encoding="utf-8"))


def _rss_kb() -> int:
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


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
        "primer_search",
    ):
        if hasattr(_core, name):
            setattr(_core, name, wrap(name, getattr(_core, name)))
    return totals


def _window_count(span: int, primer_min: int, primer_max: int) -> int:
    total = 0
    for length in range(primer_min, primer_max + 1):
        if span >= length:
            total += span - length + 1
    return total


def _design_row(spec: dict, sequences: dict[str, str]) -> dict:
    from dna_design_suite import _core
    from dna_design_suite.codon.orchestrator.design import DesignRequest, design
    from dna_design_suite.modelreg import ArtifactStore

    totals = _wrap_core()
    cds = sequences[spec["gene"]]
    constraints = _constraints(spec["constraints"])
    weights = _weights()
    started = time.perf_counter()
    result = design(
        DesignRequest(
            host=HOST,
            cds=cds,
            host_profile=_profile(),
            objectives=list(spec["objectives"]),
            constraints=constraints,
            translation_table=11,
            budget=spec.get("budget", "fast"),
            seed=SEED,
            generations=int(spec["generations"]),
            stagnation_window=21,
            segment=2,
            host_digest=DIGEST,
        ),
        ArtifactStore(Path("/tmp") / f"cf-harsh-{_sha(spec['name'])[:12]}"),
    )
    wall = time.perf_counter() - started
    payload: dict = {
        "state": result.job.state,
        "wall_s": wall,
        "peak_rss_kb": _rss_kb(),
        "nt": len(cds),
        "stages": {
            name: {"seconds": round(item["seconds"], 6), "calls": item["calls"]}
            for name, item in sorted(totals.items())
        },
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
        actionable = [item for item in report["violations"] if item.get("actionable")]
        independent = cai_independent.score(designed, weights)
        payload.update(
            {
                "vector": chosen["objective_vector"],
                "cds_sha256": _sha(designed),
                "cai": independent["cai"],
                "engine_cai": chosen["objective_vector"][0]
                if spec["objectives"][0] == "usage.cai"
                else None,
                "counted_codons": independent["counted_codons"],
                "non_top_codons": independent["non_top_codons"],
                "protein_identical": _core.translate(designed, 11) == _core.translate(cds, 11),
                "actionable_violations": len(actionable),
                "actionable_ids": sorted({item.get("id") for item in actionable}),
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
    payload["pin"] = _design_pin(spec, payload)
    return payload


def _design_pin(spec: dict, payload: dict) -> dict:
    kind = spec.get("pin")
    if kind == "best_effort":
        return {"ok": True, "detail": str(payload.get("state"))}
    if payload.get("state") != "succeeded":
        return {"ok": False, "detail": payload.get("state")}
    if not payload.get("protein_identical", False):
        return {"ok": False, "detail": "protein changed"}
    cai = payload.get("cai")
    if kind == "empty":
        if spec["gene"] == "rpoC":
            # rpoC starts with GTG. The start codon is locked, and table 11
            # reads GTG as valine, so that one codon stays at its table weight
            # and CAI cannot be 1. Same reason the panel drops tufA.
            counted = payload.get("counted_codons") or 0
            weight = _weights()["GTG"]
            expected = math.exp(math.log(weight) / counted) if counted else None
            ok = (
                payload.get("non_top_codons") == 1
                and expected is not None
                and cai is not None
                and abs(cai - expected) <= 1e-12
                and payload.get("actionable_violations") == 0
            )
            return {"ok": ok, "detail": f"locked GTG start cai={cai} expected={expected}"}
        ok = cai == 1.0 and payload.get("actionable_violations") == 0
        return {"ok": ok, "detail": f"cai={cai}"}
    if kind == "assembly_floor":
        floor = ASSEMBLY_FLOOR[spec["gene"]]
        ok = cai is not None and cai + 1e-12 >= floor and payload.get("actionable_violations") == 0
        return {"ok": ok, "detail": f"cai={cai} floor={floor}"}
    if kind == "dp":
        target = float(spec["dp_cai"])
        gap = None if cai is None else target - cai
        ok = (
            cai is not None
            and abs(cai - target) <= 1e-12
            and payload.get("actionable_violations") == 0
        )
        return {"ok": ok, "detail": f"cai={cai} dp={target} gap={gap}"}
    if kind == "structure":
        # Pinned by test_yeej_structure_vector_and_fold_counts. The CAI axis
        # sits under the assembly-only seed because the structure term moves it.
        expected = [0.9723247970586831, 2.4166863375261896]
        folds = payload.get("folds") or {}
        ok = (
            payload.get("vector") == expected
            and payload.get("actionable_violations") == 0
            and folds.get("hits") == 12
            and folds.get("full_folds") == 1
            and folds.get("local_edits") == 18
        )
        return {"ok": ok, "detail": f"vector={payload.get('vector')} folds={folds}"}
    clean = payload.get("actionable_violations", 0) == 0
    return {"ok": clean, "detail": "recorded"}


def _seed_row(spec: dict, sequences: dict[str, str]) -> dict:
    from dna_design_suite import _core

    cds = sequences[spec["gene"]]
    request = {
        "cds": cds,
        "translation_table": 11,
        "five_prime_flank": "",
        "three_prime_flank": "",
        "constraint_set": _constraints(spec["constraints"]),
        "host_profile": _profile(),
        "host_id": HOST,
        "host_digest": DIGEST,
        "terms": ["usage.cai"],
        "weights": [1.0],
        "root_seed": SEED,
    }
    started = time.perf_counter()
    report = json.loads(_core.resolve_feasible_seed(json.dumps(request)))
    designed = report["cds"]
    independent = cai_independent.score(designed, _weights())
    return {
        "state": "succeeded" if report["feasible"] else "infeasible",
        "wall_s": time.perf_counter() - started,
        "peak_rss_kb": _rss_kb(),
        "nt": len(cds),
        "feasible": report["feasible"],
        "edits": report["edits"],
        "cycles": report["cycles"],
        "cds_sha256": _sha(designed),
        "cai": independent["cai"],
        "pin": _seed_pin(spec, bool(report["feasible"]), independent["cai"]),
    }


def _seed_pin(spec: dict, feasible: bool, cai: float | None) -> dict:
    """Assembly-only rows keep the search-escape floor. Tighter sets may score less."""
    if spec.get("constraints") == HARSH:
        return {"ok": True, "detail": f"feasible={feasible} cai={cai}"}
    if not feasible or cai is None:
        return {"ok": False, "detail": f"feasible={feasible} cai={cai}"}
    if spec["constraints"] == ASSEMBLY and spec["gene"] in ASSEMBLY_FLOOR:
        floor = ASSEMBLY_FLOOR[spec["gene"]]
        return {"ok": cai + 1e-12 >= floor, "detail": f"cai={cai} floor={floor}"}
    return {"ok": True, "detail": f"cai={cai}"}


def _primer_row(spec: dict, sequences: dict[str, str]) -> dict:
    from dna_design_suite import _core

    sequence = sequences[spec["gene"]]
    flank = int(spec["flank"])
    primer_min = int(spec.get("primer_min", 18))
    primer_max = int(spec.get("primer_max", 25))
    if spec.get("region") in {"pin", "tail"}:
        start, end = flank, len(sequence) - flank
    else:
        start, end = int(spec["start"]), int(spec["end"])
    region = end - start
    product_min = int(spec["product_min"]) if "product_min" in spec else region
    body = {
        "sequence": sequence,
        "region": {"start": start, "end": end},
        "config": {
            "tm_target": 60.0,
            "product_size_range": {"min": product_min, "max": region + 2 * flank},
            "flank_length": flank,
            "primer_length_range": {"min": primer_min, "max": primer_max},
        },
    }
    windows = _window_count(min(flank, start), primer_min, primer_max) + _window_count(
        min(flank, len(sequence) - end), primer_min, primer_max
    )
    started = time.perf_counter()
    raw = _core.primer_search(json.dumps(body), lambda _fraction: None)
    wall = time.perf_counter() - started
    result = json.loads(raw)
    recommended = (result.get("pairs") or [{}])[0]
    pin = {"ok": True, "detail": "recorded"}
    if spec.get("pin") == "yeej300":
        ok = (
            recommended.get("forward_sequence") == PRIMER_PIN["forward_sequence"]
            and recommended.get("reverse_sequence") == PRIMER_PIN["reverse_sequence"]
            and recommended.get("product_size") == PRIMER_PIN["product_size"]
            and recommended.get("hetero_dimer_dG") == PRIMER_PIN["hetero_dimer_dG"]
            and recommended.get("refused") is False
        )
        pin = {"ok": ok, "detail": "yeeJ 300 nt flank pair"}
    return {
        "state": "succeeded",
        "wall_s": wall,
        "peak_rss_kb": _rss_kb(),
        "nt": len(sequence),
        "template_nt": len(sequence),
        "region_nt": region,
        "flank": flank,
        "primer_length": [primer_min, primer_max],
        "oligo_windows": windows,
        "survivors": result.get("survivors"),
        "truncated": result.get("truncated"),
        "n_pairs": len(result.get("pairs") or []),
        "response_bytes": len(raw.encode()),
        "forward_sequence": recommended.get("forward_sequence"),
        "reverse_sequence": recommended.get("reverse_sequence"),
        "product_size": recommended.get("product_size"),
        "hetero_dimer_dG": recommended.get("hetero_dimer_dG"),
        "refused": recommended.get("refused"),
        "document_sha256": _sha(raw),
        "pin": pin,
    }


def _operon_row(spec: dict, sequences: dict[str, str]) -> dict:
    from dna_design_suite import _core
    from dna_design_suite.codon.orchestrator import Cistron, design_operon
    from dna_design_suite.modelreg import ArtifactStore

    cistrons = [
        Cistron(sequences[gene], float(flux), gene)
        for gene, flux in zip(spec["genes"], spec["fluxes"], strict=True)
    ]
    started = time.perf_counter()
    result = design_operon(
        cistrons,
        ArtifactStore(Path("/tmp") / f"cf-harsh-operon-{_sha(spec['name'])[:12]}"),
        host=HOST,
        host_profile=_profile(),
        objectives=["burden.growth_cost"],
        translation_table=11,
        budget="fast",
        seed=SEED,
        generations=int(spec.get("generations", 0)),
        passes=int(spec["passes"]),
        host_digest=DIGEST,
    )
    return {
        "state": "succeeded",
        "wall_s": time.perf_counter() - started,
        "peak_rss_kb": _rss_kb(),
        "nt": [len(sequence) for sequence in result.sequences],
        "joint_growth_cost": result.joint_growth_cost,
        "swaps_per_pass": result.swaps_per_pass,
        "protein_identical": [
            _core.translate(designed, 11) == protein
            for designed, protein in zip(result.sequences, result.proteins, strict=True)
        ],
        "cds_sha256": [_sha(sequence) for sequence in result.sequences],
        "pin": {
            "ok": all(
                _core.translate(designed, 11) == protein
                for designed, protein in zip(result.sequences, result.proteins, strict=True)
            ),
            "detail": f"joint={result.joint_growth_cost}",
        },
    }


def _http(method: str, url: str, body: dict | None = None, headers: dict | None = None):
    data = None if body is None else json.dumps(body).encode()
    request = urllib.request.Request(url, data=data, headers=headers or {}, method=method)
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            raw = response.read()
            return response.status, json.loads(raw) if raw else None
    except urllib.error.HTTPError as exc:
        raw = exc.read()
        return exc.code, json.loads(raw) if raw else None


def _api_row(_spec: dict, sequences: dict[str, str]) -> dict:
    """Concurrent accepts, one serialized design worker, refused routes."""
    from dna_design_suite.api import serve
    from dna_design_suite.codon.orchestrator.design import DesignRequest, design
    from dna_design_suite.modelreg import ArtifactStore

    handle = design(
        DesignRequest(
            host=HOST,
            cds="ATGGGGTAA",
            host_profile=_profile(),
            objectives=["usage.cai"],
            constraints={
                "schema_version": "1.0.0",
                "constraint_set_id": "cst_none",
                "constraints": [],
            },
            generations=0,
            budget="fast",
        ),
        ArtifactStore(Path("/tmp/cf-harsh-bg")),
        background=True,
    )
    background_result = handle.result(timeout=60)
    background_ok = bool(handle.job_id) and background_result.job.state == "succeeded"

    server = serve(ArtifactStore(Path("/tmp/cf-harsh-api")))
    base = f"http://127.0.0.1:{server.server_address[1]}"
    headers = {
        "Authorization": "Bearer t_localdev1",
        "Content-Type": "application/json",
    }
    hosts_status, hosts_body = _http("GET", f"{base}/v1/hosts", headers=headers)
    empty_feas_status, _empty_feas = _http("POST", f"{base}/v1/feasibility", {}, headers)
    rpsl = sequences["rpsL"]
    feas_status, feas_body = _http(
        "POST",
        f"{base}/v1/feasibility",
        {
            "cds": rpsl,
            "constraints": _constraints(None),
            "translation_table": 11,
        },
        headers,
    )
    jobs = []
    for index, (gene, objectives) in enumerate(
        (
            ("rpsL", ["usage.cai"]),
            ("waaS", ["usage.cai"]),
            ("rpsL", ["usage.cai", "burden.growth_cost"]),
        )
    ):
        jobs.append(
            {
                "host": HOST,
                "cds": sequences[gene] if gene != "rpsL" else rpsl,
                "host_profile": _profile(),
                "objectives": objectives,
                "constraints": _constraints(ASSEMBLY if gene == "waaS" else None),
                "translation_table": 11,
                "budget": "fast",
                "seed": SEED,
                "generations": 0,
                "segment": 2,
                "host_digest": DIGEST,
                "label": f"{gene}-{index}",
            }
        )
    started = time.perf_counter()
    accepted = []
    import threading

    errors: list[str] = []

    def post(job: dict) -> None:
        label = job.pop("label")
        status, payload = _http("POST", f"{base}/v1/jobs", job, headers)
        accepted.append((label, status, (payload or {}).get("job_id"), time.perf_counter()))

    threads = [threading.Thread(target=post, args=(dict(job),)) for job in jobs]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    accept_wall = time.perf_counter() - started
    deadline = time.perf_counter() + 120
    finished = []
    pending = {
        job_id: label for label, status, job_id, _when in accepted if status == 202 and job_id
    }
    while pending and time.perf_counter() < deadline:
        for job_id, label in list(pending.items()):
            status, payload = _http("GET", f"{base}/v1/jobs/{job_id}", headers=headers)
            state = (payload or {}).get("state")
            if state in {"succeeded", "failed", "cancelled"}:
                vector = None
                if state == "succeeded":
                    _front_status, front = _http(
                        "GET", f"{base}/v1/jobs/{job_id}/front", headers=headers
                    )
                    candidates = (front or {}).get("candidates") or []
                    if candidates:
                        vector = candidates[0].get("objective_vector")
                finished.append(
                    {
                        "label": label,
                        "state": state,
                        "http": status,
                        "vector": vector,
                        "problem": (payload or {}).get("problem", {}).get("code")
                        if isinstance((payload or {}).get("problem"), dict)
                        else None,
                    }
                )
                del pending[job_id]
        if pending:
            time.sleep(0.05)
    dwell_status, dwell_body = _http(
        "POST",
        f"{base}/v1/jobs",
        {
            "host": HOST,
            "cds": rpsl,
            "host_profile": _profile(),
            "objectives": ["usage.cai", "epistasis.dwell_deviation"],
            "constraints": _constraints(None),
            "translation_table": 11,
            "budget": "fast",
            "seed": SEED,
            "generations": 0,
            "host_digest": DIGEST,
        },
        headers,
    )
    dwell_code = None
    if dwell_status == 202 and isinstance(dwell_body, dict) and dwell_body.get("job_id"):
        dwell_deadline = time.perf_counter() + 30
        while time.perf_counter() < dwell_deadline:
            _status, payload = _http(
                "GET", f"{base}/v1/jobs/{dwell_body['job_id']}", headers=headers
            )
            if (payload or {}).get("state") in {"failed", "succeeded", "cancelled"}:
                problem = (payload or {}).get("problem") or {}
                dwell_code = problem.get("code") if isinstance(problem, dict) else None
                break
            time.sleep(0.05)
    else:
        dwell_code = (dwell_body or {}).get("code") if isinstance(dwell_body, dict) else None
    primer_started = time.perf_counter()
    yeej = sequences["yeeJ"]
    flank = 300
    region = len(yeej) - 2 * flank
    primer_status, primer_body = _http(
        "POST",
        f"{base}/v1/primers/candidates",
        {
            "sequence": yeej,
            "region": {"start": flank, "end": len(yeej) - flank},
            "config": {
                "tm_target": 60.0,
                "product_size_range": {"min": region, "max": region + 2 * flank},
                "flank_length": flank,
                "primer_length_range": {"min": 18, "max": 25},
            },
        },
        headers,
    )
    primer_job = (primer_body or {}).get("job_id") if isinstance(primer_body, dict) else None
    primer_pair = None
    if primer_status == 202 and primer_job:
        primer_deadline = time.perf_counter() + 60
        while time.perf_counter() < primer_deadline:
            status, payload = _http("GET", f"{base}/v1/primers/{primer_job}", headers=headers)
            if (payload or {}).get("state") == "succeeded":
                pairs = ((payload or {}).get("primer_result") or {}).get("pairs") or []
                primer_pair = pairs[0] if pairs else None
                break
            if (payload or {}).get("state") == "failed":
                errors.append(str((payload or {}).get("problem")))
                break
            time.sleep(0.05)
    primer_wall = time.perf_counter() - primer_started
    server.shutdown()
    hosts_ok = (
        hosts_status == 200
        and isinstance(hosts_body, dict)
        and any(item.get("installed") for item in hosts_body.get("hosts") or [])
    )
    checks = (feas_body or {}).get("checks_run") if isinstance(feas_body, dict) else None
    feas_ok = (
        empty_feas_status == 400
        and feas_status == 200
        and isinstance(checks, list)
        and len(checks) == 5
    )
    primer_ok = (
        primer_pair is not None
        and primer_pair.get("forward_sequence") == PRIMER_PIN["forward_sequence"]
        and primer_pair.get("reverse_sequence") == PRIMER_PIN["reverse_sequence"]
        and primer_pair.get("product_size") == PRIMER_PIN["product_size"]
        and primer_pair.get("hetero_dimer_dG") == PRIMER_PIN["hetero_dimer_dG"]
    )
    designs_ok = pending == {} and all(item["state"] == "succeeded" for item in finished)
    pin_ok = (
        background_ok
        and hosts_ok
        and feas_ok
        and dwell_code == "CF-LAND-001"
        and designs_ok
        and primer_ok
    )
    return {
        "state": "succeeded" if pin_ok else "failed",
        "wall_s": time.perf_counter() - started,
        "peak_rss_kb": _rss_kb(),
        "accept_wall_s": accept_wall,
        "accepted": [
            {"label": label, "http": status, "job_id": job_id}
            for label, status, job_id, _when in accepted
        ],
        "finished": finished,
        "still_pending": list(pending),
        "primer_wall_s": primer_wall,
        "primer_http": primer_status,
        "primer_forward": None if primer_pair is None else primer_pair.get("forward_sequence"),
        "primer_reverse": None if primer_pair is None else primer_pair.get("reverse_sequence"),
        "primer_product_size": None if primer_pair is None else primer_pair.get("product_size"),
        "primer_hetero_dimer_dG": None
        if primer_pair is None
        else primer_pair.get("hetero_dimer_dG"),
        "refusals": {
            "library_background": "handle" if background_ok else "failed",
            "GET /v1/hosts": {"http": hosts_status, "ok": hosts_ok},
            "POST /v1/feasibility": {
                "empty_http": empty_feas_status,
                "http": feas_status,
                "checks": checks,
            },
            "dwell": {"http": dwell_status, "code": dwell_code},
        },
        "errors": errors,
        "pin": {
            "ok": pin_ok,
            "detail": (
                f"background={background_ok} hosts={hosts_ok} "
                f"feasibility={feas_ok} dwell={dwell_code} primer_ok={primer_ok}"
            ),
        },
    }


def _gate_row(_spec: dict, _sequences: dict[str, str]) -> dict:
    from dna_design_suite.codon.bench.corpus import check_cai_gate, run_cai_only

    stored_gate = json.loads((ROOT / "benchmarks" / "corpus" / "gates.json").read_text())[0]
    stored_rows = json.loads((ROOT / "benchmarks" / "corpus" / "cai_only.json").read_text())
    check_cai_gate(stored_gate, stored_rows)
    started = time.perf_counter()
    fresh = run_cai_only(Path("/tmp/cf-harsh-gate"))
    wall = time.perf_counter() - started
    fresh_rows = json.loads(Path("/tmp/cf-harsh-gate/cai_only.json").read_text())
    check_cai_gate(fresh, fresh_rows)
    stored_by = {row["gene"]: row for row in stored_rows}
    moved = []
    for row in fresh_rows:
        previous = stored_by[row["gene"]]["codonforge_cai"]
        if row["codonforge_cai"] + 1e-12 < previous:
            moved.append(
                {
                    "gene": row["gene"],
                    "was": previous,
                    "now": row["codonforge_cai"],
                }
            )
    return {
        "state": "succeeded" if fresh.get("claimed") and not moved else "failed",
        "wall_s": wall,
        "peak_rss_kb": _rss_kb(),
        "claimed": fresh.get("claimed"),
        "n_claimed": fresh.get("n_claimed"),
        "codonforge_ge_greedy": fresh.get("codonforge_ge_greedy"),
        "behind_greedy": fresh.get("behind_greedy"),
        "definition": fresh.get("definition"),
        "regressed": moved,
        "pin": {
            "ok": fresh.get("claimed") is True
            and fresh.get("n_claimed") == 50
            and fresh.get("definition") == "cai_only"
            and not moved,
            "detail": f"claimed {fresh.get('n_claimed')} of 50, regressed={len(moved)}",
        },
    }


def _run_job(job: dict) -> dict:
    sequences = _sequences()
    kind = job["kind"]
    if kind == "design":
        return _design_row(job, sequences)
    if kind == "seed":
        return _seed_row(job, sequences)
    if kind == "primer":
        return _primer_row(job, sequences)
    if kind == "operon":
        return _operon_row(job, sequences)
    if kind == "api":
        return _api_row(job, sequences)
    if kind == "gate":
        return _gate_row(job, sequences)
    raise ValueError(kind)


def _dp_targets() -> list[dict]:
    document = json.loads(
        (ROOT / "docs" / "benchmarks" / "v2.9.3-quality-rows.json").read_text(encoding="utf-8")
    )
    rows = []
    for item in document["local_dp"]:
        if item.get("nt", 10**9) > 2500:
            continue
        if item.get("constraint_set") not in {GOLDEN, BIOBRICK}:
            continue
        if not item.get("dp_sequence_is_full_optimum"):
            continue
        rows.append(item)
    return rows


def _jobs(suite: str) -> list[dict]:
    rows: list[dict] = []

    def add(item: dict) -> None:
        item.setdefault("suite", suite if suite != "all" else item["suite"])
        rows.append(item)

    if suite in {"pins", "all"}:
        for gene in (
            "rpsL",
            "elastin_vpgvg20",
            "gfp_BBa_E0040",
            "mrfp_BBa_E1010",
            "mcherry_BBa_J06504",
            "bla",
            "waaS",
            "yagF",
            "lacZ",
            "yeeJ",
            "rpoB",
            "rpoC",
            "yfhM",
            "aceB",
            "thrA",
        ):
            add(
                {
                    "suite": "pins",
                    "name": f"pin {gene} empty CAI",
                    "kind": "design",
                    "gene": gene,
                    "objectives": ["usage.cai"],
                    "constraints": None,
                    "generations": 0,
                    "pin": "empty",
                }
            )
        for gene in ASSEMBLY_FLOOR:
            add(
                {
                    "suite": "pins",
                    "name": f"pin {gene} assembly CAI",
                    "kind": "design",
                    "gene": gene,
                    "objectives": ["usage.cai"],
                    "constraints": ASSEMBLY,
                    "generations": 0,
                    "pin": "assembly_floor",
                }
            )
        for item in _dp_targets():
            add(
                {
                    "suite": "pins",
                    "name": f"pin {item['gene']} {item['constraint_set']}",
                    "kind": "design",
                    "gene": item["gene"],
                    "objectives": ["usage.cai"],
                    "constraints": item["constraint_set"],
                    "generations": 0,
                    "pin": "dp",
                    "dp_cai": item["dp_cai"],
                }
            )
        add(
            {
                "suite": "pins",
                "name": "pin yeeJ primer 300",
                "kind": "primer",
                "gene": "yeeJ",
                "flank": 300,
                "region": "pin",
                "pin": "yeej300",
            }
        )
    if suite in {"design", "all"}:
        for gene in ("yeeJ", "yfhM", "rpoC", "rpoB", "lacZ"):
            add(
                {
                    "suite": "design",
                    "name": f"{gene} CAI assembly",
                    "kind": "design",
                    "gene": gene,
                    "objectives": ["usage.cai"],
                    "constraints": ASSEMBLY,
                    "generations": 0,
                    "pin": "assembly_floor" if gene in ASSEMBLY_FLOOR else "recorded",
                }
            )
        for gene in ("yeeJ", "waaS", "elastin_vpgvg20"):
            add(
                {
                    "suite": "design",
                    "name": f"{gene} CAI harsh assembly",
                    "kind": "design",
                    "gene": gene,
                    "objectives": ["usage.cai"],
                    "constraints": HARSH,
                    "generations": 0,
                    "pin": "best_effort",
                }
            )
        add(
            {
                "suite": "design",
                "name": "yeeJ CAI assembly standard x4",
                "kind": "design",
                "gene": "yeeJ",
                "objectives": ["usage.cai"],
                "constraints": ASSEMBLY,
                "generations": 4,
                "budget": "standard",
                "pin": "assembly_floor",
            }
        )
        add(
            {
                "suite": "design",
                "name": "yfhM CAI assembly standard x4",
                "kind": "design",
                "gene": "yfhM",
                "objectives": ["usage.cai"],
                "constraints": ASSEMBLY,
                "generations": 4,
                "budget": "standard",
            }
        )
    if suite in {"seed", "all"}:
        for gene, constraints in (
            ("yeeJ", ASSEMBLY),
            ("yfhM", ASSEMBLY),
            ("rpoC", ASSEMBLY),
            ("elastin_vpgvg20", ASSEMBLY),
            ("waaS", ASSEMBLY),
            ("yeeJ", HARSH),
            ("elastin_vpgvg20", HARSH),
        ):
            add(
                {
                    "suite": "seed",
                    "name": f"seed {gene} {constraints}",
                    "kind": "seed",
                    "gene": gene,
                    "constraints": constraints,
                }
            )
    if suite in {"objectives", "all"}:
        add(
            {
                "suite": "objectives",
                "name": "yeeJ CAI+structure assembly",
                "kind": "design",
                "gene": "yeeJ",
                "objectives": ["usage.cai", "structure.local_accessibility"],
                "constraints": ASSEMBLY,
                "generations": 2,
                "pin": "structure",
            }
        )
        add(
            {
                "suite": "objectives",
                "name": "rpoB CAI+structure assembly",
                "kind": "design",
                "gene": "rpoB",
                "objectives": ["usage.cai", "structure.local_accessibility"],
                "constraints": ASSEMBLY,
                "generations": 2,
            }
        )
        add(
            {
                "suite": "objectives",
                "name": "yeeJ CAI+burden assembly",
                "kind": "design",
                "gene": "yeeJ",
                "objectives": ["usage.cai", "burden.growth_cost"],
                "constraints": ASSEMBLY,
                "generations": 0,
            }
        )
        add(
            {
                "suite": "objectives",
                "name": "waaS CAI+structure+burden assembly",
                "kind": "design",
                "gene": "waaS",
                "objectives": [
                    "usage.cai",
                    "structure.local_accessibility",
                    "burden.growth_cost",
                ],
                "constraints": ASSEMBLY,
                "generations": 2,
            }
        )
        for gene in ("yeeJ", "rpoC", "lacZ", "rpsL"):
            add(
                {
                    "suite": "objectives",
                    "name": f"{gene} CAI+Potts",
                    "kind": "design",
                    "gene": gene,
                    "objectives": ["usage.cai", "epistasis.context_naturalness"],
                    "constraints": None,
                    "generations": 0,
                }
            )
    if suite in {"operon", "all"}:
        add(
            {
                "suite": "operon",
                "name": "operon rpsL+waaS",
                "kind": "operon",
                "genes": ["rpsL", "waaS"],
                "fluxes": [1.0, 2.0],
                "passes": 2,
            }
        )
        add(
            {
                "suite": "operon",
                "name": "operon lacZ+thrA",
                "kind": "operon",
                "genes": ["lacZ", "thrA"],
                "fluxes": [1.0, 2.0],
                "passes": 2,
            }
        )
        add(
            {
                "suite": "operon",
                "name": "operon yeeJ+rpoB",
                "kind": "operon",
                "genes": ["yeeJ", "rpoB"],
                "fluxes": [1.0, 2.0],
                "passes": 2,
            }
        )
    if suite in {"primers", "all"}:
        add(
            {
                "suite": "primers",
                "name": "yeeJ template 300 flank pin",
                "kind": "primer",
                "gene": "yeeJ",
                "flank": 300,
                "region": "pin",
                "pin": "yeej300",
            }
        )
        for gene in ("yfhM", "rpoC"):
            add(
                {
                    "suite": "primers",
                    "name": f"{gene} long template flank 300",
                    "kind": "primer",
                    "gene": gene,
                    "flank": 300,
                    "region": "tail",
                }
            )
        add(
            {
                "suite": "primers",
                "name": "yeeJ oligo pool 12-36",
                "kind": "primer",
                "gene": "yeeJ",
                "flank": 300,
                "region": "pin",
                "primer_min": 12,
                "primer_max": 36,
            }
        )
        add(
            {
                "suite": "primers",
                "name": "yfhM oligo pool 12-36",
                "kind": "primer",
                "gene": "yfhM",
                "flank": 300,
                "start": 300,
                "end": 4300,
                "primer_min": 12,
                "primer_max": 36,
            }
        )
    if suite in {"api", "all"}:
        add({"suite": "api", "name": "api concurrent and refusals", "kind": "api"})
    if suite in {"gate", "all"}:
        add({"suite": "gate", "name": "GATE-BEN-001 cai_only", "kind": "gate"})
    return rows


def _child_env(job: dict) -> dict[str, str]:
    env = os.environ.copy()
    if job["kind"] in {"seed", "design"} and job.get("constraints"):
        env["CF_SEED_PROFILE"] = "1"
    if job["kind"] == "primer":
        env["CF_PRIMER_PROFILE"] = "1"
    objectives = job.get("objectives") or []
    if job["kind"] == "design" and (
        "structure.local_accessibility" in objectives or int(job.get("generations") or 0) > 0
    ):
        env["CF_TIER3_PROFILE"] = "1"
    return env


def _parse_profile(stderr: str) -> dict[str, list]:
    found: dict[str, list] = {}
    for line in stderr.splitlines():
        for prefix in ("CF_SEED_PROFILE ", "CF_PRIMER_PROFILE ", "CF_TIER3_PROFILE "):
            if line.startswith(prefix):
                found[prefix.strip()] = json.loads("[" + line[len(prefix) :] + "]")
    return found


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
            "pin": {"ok": False, "detail": str(exc)},
        }
    json.dump(measured, sys.stdout)
    return 0


def _public_spec(job: dict) -> dict:
    skip = {"suite"}
    return {key: value for key, value in job.items() if key not in skip}


def _one(job: dict, timeout: int) -> dict:
    try:
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--worker"],
            input=json.dumps(job),
            text=True,
            capture_output=True,
            timeout=timeout,
            cwd=ROOT,
            check=False,
            env=_child_env(job),
        )
    except subprocess.TimeoutExpired:
        return {
            "name": job["name"],
            "suite": job["suite"],
            "kind": job["kind"],
            "state": "timeout",
            "wall_s": timeout,
            "pin": {"ok": False, "detail": "timeout"},
        }
    profiles = _parse_profile(completed.stderr)
    if completed.returncode != 0 or not completed.stdout.strip():
        return {
            "name": job["name"],
            "suite": job["suite"],
            "kind": job["kind"],
            "state": "error",
            "detail": (completed.stderr or completed.stdout)[-2000:],
            "returncode": completed.returncode,
            "rust_profile": profiles,
            "pin": {"ok": False, "detail": "worker failed"},
        }
    measured = json.loads(completed.stdout)
    measured["name"] = job["name"]
    measured["suite"] = job["suite"]
    measured["kind"] = job["kind"]
    measured["spec"] = _public_spec(job)
    if profiles:
        measured["rust_profile"] = profiles
    return measured


def _machine() -> dict:
    cpu = "unknown"
    try:
        for line in Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    viennarna = None
    try:
        import RNA

        viennarna = RNA.__version__
    except ImportError:
        viennarna = None
    return {
        "python": sys.version.split()[0],
        "cpu": cpu,
        "viennarna": viennarna,
        "package": "0.1.0",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument(
        "--suite",
        choices=(
            "pins",
            "design",
            "seed",
            "objectives",
            "operon",
            "primers",
            "api",
            "gate",
            "all",
        ),
        default="all",
    )
    parser.add_argument("--out", type=Path)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--only", help="substring of the row name")
    args = parser.parse_args(argv)
    if args.worker:
        return _worker()
    if args.out is None:
        parser.error("--out is required")
    jobs = _jobs(args.suite)
    if args.only:
        jobs = [job for job in jobs if args.only in job["name"]]
    rows = []
    for job in jobs:
        print(job["name"], flush=True)
        row = _one(job, args.timeout)
        rows.append(row)
        pin = (row.get("pin") or {}).get("ok")
        print(
            f"  {row.get('state')} {row.get('wall_s')} s pin={pin} rss={row.get('peak_rss_kb')}",
            flush=True,
        )
    failed = [row["name"] for row in rows if not (row.get("pin") or {}).get("ok", False)]
    document = {
        "seed": SEED,
        "host": HOST,
        "suite": args.suite,
        "timeout_s": args.timeout,
        "machine": _machine(),
        "pin_failures": failed,
        "rows": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.out} failures={len(failed)}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
