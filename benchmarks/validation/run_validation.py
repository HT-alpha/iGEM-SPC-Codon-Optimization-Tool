"""Run the v2.5.0 validation round.

CAI arithmetic is `cai_independent.py`, which does not import CodonForge.
Designs and primer searches shell out to the installed package or its HTTP API.
Numbers in the JSON under results/ are whatever those processes returned.
"""

from __future__ import annotations

import json
import platform
import subprocess
import sys
import time
import urllib.error
import urllib.request
from collections import Counter
from pathlib import Path

from cai_independent import (
    AA,
    constraint_census,
    gc_fraction,
    load_kazusa_weights,
    load_ribosomal_weights,
    most_frequent,
    score,
    score_with_floor,
    translate,
)

ROOT = Path(__file__).resolve().parents[2]
VAL = Path(__file__).resolve().parent
SEQ = VAL / "sequences"
OUT = VAL / "results"
PY = ROOT / ".venv" / "bin" / "python"
HOST = "escherichia_coli:511145/MG1655"
PROFILE = ROOT / "examples" / "mg1655_ribosomal_usage.json"
EMPTY = ROOT / "examples" / "rpsL_constraints.json"
ASSEMBLY = VAL / "constraints_assembly.json"
CRYPTIC = VAL / "constraints_cryptic_probe.json"
KAZUSA = VAL / "kazusa" / "e_coli_316407.csv"
TOKEN = "t_localdemo01"
API = "http://127.0.0.1:8091"

PERF = [
    "rpsL",
    "lacZ",
    "BBa_E0040",
    "BBa_E0040_onestop",
    "BBa_E1010",
    "BBa_J06504",
    "bla",
]
HARD = ["yagF", "waaS", "elp_vpgvg", "yeeJ"]
SITES = {
    "c_bsai": "BsaI",
    "c_bsmbi": "BsmBI",
    "c_ecori": "EcoRI",
    "c_xbai": "XbaI",
    "c_spei": "SpeI",
    "c_psti": "PstI",
    "c_noti": "NotI",
}


def read_fa(name: str) -> tuple[str, str]:
    lines = (SEQ / f"{name}.fa").read_text(encoding="utf-8").splitlines()
    seq = "".join(lines[1:]).replace(" ", "").upper()
    return lines[0][1:], seq


def meta() -> dict[str, object]:
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    cpu = "unknown"
    for line in Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines():
        if line.startswith("model name"):
            cpu = line.split(":", 1)[1].strip()
            break
    return {
        "commit": commit,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "cpu": cpu,
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }


def zero_codons(seq: str, weights: dict[str, float]) -> list[dict[str, object]]:
    found = []
    for index in range(0, len(seq), 3):
        codon = seq[index : index + 3]
        if AA[codon] in {"*", "M", "W"}:
            continue
        if weights[codon] == 0.0:
            found.append({"codon_index": index // 3, "codon": codon, "amino_acid": AA[codon]})
    return found


def run_cli(args: list[str], timeout: int) -> dict[str, object]:
    started = time.perf_counter()
    try:
        proc = subprocess.run(
            [str(PY), "-m", "codonforge.cli.main", *args],
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        return {
            "timeout": True,
            "timeout_s": timeout,
            "wall_s": time.perf_counter() - started,
            "stdout": (exc.stdout or "")[-2000:],
            "stderr": (exc.stderr or "")[-2000:],
        }
    return {
        "timeout": False,
        "returncode": proc.returncode,
        "wall_s": time.perf_counter() - started,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "runtime_warning": "RuntimeWarning" in proc.stderr,
    }


def design(name: str, seq: str, constraints: Path, budget: str, timeout: int) -> dict[str, object]:
    store = Path(f"/tmp/val-store/{name}-{budget}-{constraints.stem}")
    store.mkdir(parents=True, exist_ok=True)
    args = [
        "design",
        "--host",
        HOST,
        "--cds",
        seq,
        "--host-profile",
        str(PROFILE),
        "--constraints",
        str(constraints),
        "--objectives",
        "usage.cai",
        "--budget",
        budget,
        "--seed",
        "7",
        "--output",
        "json",
        "--store",
        str(store),
    ]
    if constraints != EMPTY:
        args.extend(["--generations", "12"])
    raw = run_cli(args, timeout)
    parsed: dict[str, object] = {
        "timeout": raw["timeout"],
        "wall_s": raw["wall_s"],
        "returncode": raw.get("returncode"),
        "runtime_warning": raw.get("runtime_warning"),
        "stderr": raw["stderr"],
    }
    stdout = raw.get("stdout") or ""
    if stdout.strip().startswith("{"):
        document = json.loads(stdout)
        status = document["status"]
        front = document.get("front")
        parsed.update(
            {
                "state": status["state"],
                "problem": status.get("problem"),
                "job_id": status["job_id"],
                "warnings": status.get("warnings") or [],
                "tiers": status["progress"]["tiers"],
            }
        )
        if not front or not front.get("candidates"):
            return parsed
        candidate = front["candidates"][0]
        order = front["objective_order"]
        parsed.update(
            {
                "state": status["state"],
                "problem": status.get("problem"),
                "job_id": status["job_id"],
                "warnings": status.get("warnings") or [],
                "tiers": status["progress"]["tiers"],
                "origin": candidate["provenance"]["origin"],
                "feasible_field": candidate["constraint_report"]["feasible"],
                "n_blocking_violations_field": candidate["constraint_report"][
                    "n_blocking_violations"
                ],
                "objective": dict(zip(order, candidate["objective_vector"], strict=True)),
                "n_candidates": len(front["candidates"]),
                "termination": (front.get("summary") or {}).get("termination_reason"),
                "output_cds": candidate["sequences"][0]["cds"],
            }
        )
    return parsed


def feasibility(seq: str, constraints: Path) -> dict[str, object]:
    raw = run_cli(
        [
            "feasibility",
            "--cds",
            seq,
            "--constraints",
            str(constraints),
            "--translation-table",
            "11",
            "--output",
            "json",
        ],
        120,
    )
    if not str(raw.get("stdout") or "").strip().startswith("{"):
        return {
            "returncode": raw.get("returncode"),
            "timeout": raw["timeout"],
            "stderr": (raw.get("stderr") or "")[-2000:],
        }
    document = json.loads(raw["stdout"])
    counts = Counter(item["constraint_id"] for item in document["violations"])
    return {
        "returncode": raw["returncode"],
        "n_violations": len(document["violations"]),
        "by_id": dict(sorted(counts.items())),
        "undetectable": document.get("undetectable") or [],
        "wall_s": raw["wall_s"],
    }


def score_pair(seq: str, rib: dict[str, float], kaz: dict[str, float]) -> dict[str, object]:
    ribosomal = score(seq, rib)
    row: dict[str, object] = {
        "length_nt": len(seq),
        "gc": gc_fraction(seq),
        "protein": translate(seq) if len(seq) % 3 == 0 else None,
        "ribosomal": ribosomal,
        "kazusa_w": score(seq, kaz),
    }
    if ribosomal["zero_weight_codons"]:
        row["zero_weight_positions"] = zero_codons(seq, rib)
        row["cai_floor_sensitivity"] = score_with_floor(seq, rib)
        row["floor_note"] = (
            "Replaces each 0 weight with the smallest positive weight in the "
            "ribosomal table. Not Sharp & Li."
        )
    return row


def expected_counts(census: dict[str, object]) -> dict[str, int]:
    sites = census["restriction_sites"]
    assert isinstance(sites, dict)
    counts = {
        "c_gc50": int(census["gc_windows_outside_0.4_0.6"]),
        "c_homopolymer": int(census["homopolymer_runs_longer_than_6"]),
        "c_sd": int(census["sd_aggagg_both_strands"]),
        "c_kmer12": int(census["kmer12_over_limit"]),
    }
    for constraint_id, enzyme in SITES.items():
        counts[constraint_id] = int(sites[enzyme])
    return {key: value for key, value in counts.items() if value}


def chisel(seq: str, timeout: int) -> dict[str, object]:
    """DNA Chisel in a child process so a stall can be killed."""
    script = r"""
import json, sys, traceback
seq = sys.stdin.read().strip()
try:
    from dnachisel import (
        AvoidPattern, CodonOptimize, DnaOptimizationProblem, EnforceGCContent,
        EnforceTranslation, EnzymeSitePattern, HomopolymerPattern,
    )
    enzymes = ["BsaI", "BsmBI", "EcoRI", "XbaI", "SpeI", "PstI", "NotI"]
    problem = DnaOptimizationProblem(
        sequence=seq,
        constraints=[
            EnforceTranslation(),
            EnforceGCContent(mini=0.4, maxi=0.6, window=50),
            *[AvoidPattern(EnzymeSitePattern(name)) for name in enzymes],
            AvoidPattern(HomopolymerPattern("A", 7)),
            AvoidPattern(HomopolymerPattern("C", 7)),
            AvoidPattern(HomopolymerPattern("G", 7)),
            AvoidPattern(HomopolymerPattern("T", 7)),
            AvoidPattern("AGGAGG"),
            AvoidPattern("CCTCCT"),
        ],
        objectives=[CodonOptimize(species="e_coli", method="use_best_codon")],
    )
    problem.resolve_constraints()
    problem.optimize()
    print(json.dumps({"ok": True, "sequence": str(problem.sequence)}))
except Exception as exc:
    print(json.dumps({"ok": False, "error_type": type(exc).__name__, "error": str(exc)[:2000]}))
"""
    started = time.perf_counter()
    try:
        proc = subprocess.run(
            [str(PY), "-c", script],
            input=seq,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {"timeout": True, "timeout_s": timeout, "wall_s": time.perf_counter() - started}
    payload = proc.stdout.strip().splitlines()
    parsed = json.loads(payload[-1]) if payload else {"ok": False, "error": proc.stderr[-2000:]}
    parsed["timeout"] = False
    parsed["wall_s"] = time.perf_counter() - started
    parsed["returncode"] = proc.returncode
    if proc.stderr:
        parsed["stderr_tail"] = proc.stderr[-500:]
    return parsed


def summarize_design(
    name: str,
    seq: str,
    raw: dict[str, object],
    rib: dict[str, float],
    kaz: dict[str, float],
) -> dict[str, object]:
    row: dict[str, object] = {
        "wall_s": raw["wall_s"],
        "returncode": raw.get("returncode"),
        "timeout": raw["timeout"],
        "runtime_warning": raw.get("runtime_warning"),
        "stderr": raw.get("stderr") if raw.get("returncode") not in (0, None) or raw["timeout"] else "",
    }
    if not raw.get("output_cds"):
        row.update(
            {
                "state": raw.get("state"),
                "problem": raw.get("problem"),
                "tiers": raw.get("tiers"),
                "warnings": [
                    {"code": item.get("code"), "title": item.get("title")}
                    for item in raw.get("warnings") or []
                ],
                "stderr": raw.get("stderr"),
            }
        )
        return row
    output = str(raw["output_cds"])
    census = constraint_census(output)
    product = feasibility(output, ASSEMBLY)
    row.update(
        {
            "state": raw.get("state"),
            "problem": raw.get("problem"),
            "job_id": raw.get("job_id"),
            "origin": raw.get("origin"),
            "tiers": raw.get("tiers"),
            "warnings": [
                {"code": item.get("code"), "title": item.get("title")}
                for item in raw.get("warnings") or []
            ],
            "feasible_field": raw.get("feasible_field"),
            "n_blocking_violations_field": raw.get("n_blocking_violations_field"),
            "n_blocking_note": (
                "The front stores 0 when feasible and 1 otherwise. It is not a violation count."
            ),
            "objective": raw.get("objective"),
            "n_candidates": raw.get("n_candidates"),
            "termination": raw.get("termination"),
            "same_protein": translate(seq) == translate(output),
            "same_cds": seq == output,
            "independent": score_pair(output, rib, kaz),
            "census": census,
            "product_feasibility": product,
            "census_matches_product": expected_counts(census) == product.get("by_id"),
            "output_cds": output,
        }
    )
    return row


def cmd_cai() -> None:
    rib = load_ribosomal_weights(PROFILE)
    kaz = load_kazusa_weights(KAZUSA)
    inputs = []
    for name in PERF:
        _header, seq = read_fa(name)
        row = {"name": name, **score_pair(seq, rib, kaz)}
        inputs.append(row)
        print(name, "input ribosomal", row["ribosomal"]["cai"], flush=True)

    logged = []
    jsonl = Path("/tmp/bench-out/results.jsonl")
    if jsonl.exists():
        grouped: dict[str, list[dict[str, object]]] = {}
        for line in jsonl.read_text(encoding="utf-8").splitlines():
            item = json.loads(line)
            if item.get("objectives") != "usage.cai" or item.get("path") != "cli":
                continue
            grouped.setdefault(item["name"], []).append(item)
        for name, items in grouped.items():
            cds_values = {item.get("output_cds") for item in items if item.get("output_cds")}
            succeeded = [item for item in items if item.get("state") == "succeeded"]
            sample = succeeded[0] if succeeded else items[0]
            entry: dict[str, object] = {
                "name": name,
                "n_cli_rows": len(items),
                "n_distinct_output_cds": len(cds_values),
                "states": sorted({str(item.get("state")) for item in items}),
                "product_output_cai_values": sorted(
                    {item.get("output_cai") for item in succeeded}
                ),
            }
            if len(cds_values) == 1:
                only = next(iter(cds_values))
                _header, seq = read_fa(name)
                entry["independent_output"] = score_pair(str(only), rib, kaz)
                entry["same_protein"] = translate(seq) == translate(str(only))
                entry["output_cds"] = only
            else:
                entry["output_cds"] = None
            if name == "BBa_E0040":
                entry["error"] = sample.get("error") or sample.get("stderr")
            logged.append(entry)
            print(name, "logged distinct", len(cds_values), flush=True)

    rerun = []
    for name in PERF:
        _header, seq = read_fa(name)
        print("rerun", name, flush=True)
        raw = design(name, seq, EMPTY, "fast", 120)
        entry = {
            "name": name,
            "returncode": raw.get("returncode"),
            "wall_s": raw.get("wall_s"),
            "runtime_warning": raw.get("runtime_warning"),
            "stderr": raw.get("stderr") if "output_cds" not in raw else "",
        }
        if "output_cds" in raw:
            output = str(raw["output_cds"])
            entry.update(
                {
                    "state": raw.get("state"),
                    "origin": raw.get("origin"),
                    "product_cai": (raw.get("objective") or {}).get("usage.cai"),
                    "objective": raw.get("objective"),
                    "warnings": raw.get("warnings"),
                    "tiers": raw.get("tiers"),
                    "independent": score_pair(output, rib, kaz),
                    "same_protein": translate(seq) == translate(output),
                    "output_cds": output,
                }
            )
            match = next((item for item in logged if item["name"] == name), None)
            entry["matches_performance_log_cds"] = bool(
                match and match.get("output_cds") == output
            )
        rerun.append(entry)
        print(" ", entry.get("returncode"), entry.get("product_cai"), flush=True)

    payload = {"meta": meta(), "inputs": inputs, "performance_log": logged, "rerun_seed7_fast": rerun}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "cai.json").write_text(json.dumps(payload), encoding="utf-8")
    print("wrote", OUT / "cai.json", flush=True)


def cmd_hard() -> None:
    rib = load_ribosomal_weights(PROFILE)
    kaz = load_kazusa_weights(KAZUSA)
    cryptic = feasibility(read_fa("yagF")[1][:300], CRYPTIC)
    cases = []
    for name in HARD:
        _header, seq = read_fa(name)
        print("hard", name, len(seq), flush=True)
        native_census = constraint_census(seq)
        native_product = feasibility(seq, ASSEMBLY)
        baseline = most_frequent(seq, rib)
        case: dict[str, object] = {
            "name": name,
            "header": _header,
            "input": {
                "independent": score_pair(seq, rib, kaz),
                "census": native_census,
                "product_feasibility": native_product,
                "census_matches_product": expected_counts(native_census) == native_product.get("by_id"),
            },
            "cryptic_full": feasibility(seq, CRYPTIC),
            "most_frequent_ribosomal": {
                "same_protein": translate(seq) == translate(baseline),
                "independent": score_pair(baseline, rib, kaz),
                "census": constraint_census(baseline),
                "product_feasibility": feasibility(baseline, ASSEMBLY),
            },
        }
        designs = {}
        for budget, timeout in (("fast", 180), ("standard", 240 if name != "yeeJ" else 300)):
            print(" design", budget, flush=True)
            raw = design(name, seq, ASSEMBLY, budget, timeout)
            designs[budget] = summarize_design(name, seq, raw, rib, kaz)
            print(
                "  ",
                designs[budget].get("returncode"),
                designs[budget].get("feasible_field"),
                designs[budget].get("wall_s"),
                flush=True,
            )
        case["codonforge"] = designs
        print(" dnachisel", flush=True)
        melted = chisel(seq, 180 if name != "yeeJ" else 240)
        if melted.get("ok") and melted.get("sequence"):
            out_seq = str(melted["sequence"])
            melted = {
                "ok": True,
                "wall_s": melted["wall_s"],
                "same_protein": translate(seq) == translate(out_seq),
                "independent": score_pair(out_seq, rib, kaz),
                "census": constraint_census(out_seq),
                "product_feasibility": feasibility(out_seq, ASSEMBLY),
                "output_cds": out_seq,
                "note": (
                    "Objective is DNA Chisel use_best_codon for species e_coli "
                    "(python_codon_tables, the Kazusa genome-wide table), not the "
                    "ribosomal table. RepeatedKmerPattern was not added: it matches "
                    "tandem repeats only."
                ),
            }
        case["dnachisel"] = melted
        cases.append(case)
        (OUT / "hard.json").write_text(
            json.dumps({"meta": meta(), "cryptic_probe": cryptic, "cases": cases}),
            encoding="utf-8",
        )
        print(" checkpoint", name, flush=True)
    print("wrote", OUT / "hard.json", flush=True)


def http(method: str, path: str, body: dict[str, object] | None = None) -> tuple[int, dict]:
    data = None if body is None else json.dumps(body).encode()
    request = urllib.request.Request(
        API + path,
        data=data,
        method=method,
        headers={"Authorization": f"Bearer {TOKEN}", "Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return response.status, json.loads(response.read().decode())
    except urllib.error.HTTPError as exc:
        return exc.code, json.loads(exc.read().decode())


def primer_job(body: dict[str, object]) -> dict[str, object]:
    status, posted = http("POST", "/v1/primers/candidates", body)
    if status != 202:
        return {"http": status, "problem": posted}
    job_id = posted["job_id"]
    last = posted
    for _ in range(60):
        code, last = http("GET", f"/v1/primers/{job_id}")
        if code != 200:
            return {"http": code, "problem": last, "job_id": job_id}
        if last.get("state") in {"succeeded", "failed"}:
            break
        time.sleep(0.2)
    return {"http": 200, "job_id": job_id, "status": last}


def check_pair(template: str, pair: dict[str, object]) -> dict[str, object]:
    from Bio.SeqUtils import MeltingTemp as mt

    import primer3

    forward = str(pair["forward_sequence"])
    reverse = str(pair["reverse_sequence"])
    fs, fe = int(pair["forward_start"]), int(pair["forward_end"])
    rs, re = int(pair["reverse_start"]), int(pair["reverse_end"])
    site_f = template[fs:fe]
    site_r = template[rs:re]
    rc = site_r.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    def bio(seq: str, matched: bool) -> float:
        # 1 M Na+, salt term is 0. Oligo Ct = 1 uM, so Ct/4 = 250 nM.
        return float(
            mt.Tm_NN(
                seq,
                nn_table=mt.DNA_NN4,
                Na=1000,
                saltcorr=5,
                dnac1=1000 if matched else 250,
                dnac2=0,
                selfcomp=matched,
            )
        )

    def bio_default(seq: str) -> float:
        return float(mt.Tm_NN(seq, nn_table=mt.DNA_NN4))

    aligned = dict(mv_conc=1000, dv_conc=0, dntp_conc=0, dna_conc=1000)
    hairpin_f = primer3.calc_hairpin(forward, **aligned)
    hairpin_r = primer3.calc_hairpin(reverse, **aligned)
    self_f = primer3.calc_homodimer(forward, **aligned)
    self_r = primer3.calc_homodimer(reverse, **aligned)
    hetero = primer3.calc_heterodimer(forward, reverse, **aligned)
    return {
        "forward_equals_template_slice": forward == site_f,
        "reverse_equals_revcomp": reverse == rc,
        "product_size_ok": int(pair["product_size"]) == re - fs,
        "product_sequence_ok": pair.get("product_sequence") == template[fs:re],
        "length_forward": len(forward),
        "length_reverse": len(reverse),
        "gc_forward_oligo": 100 * gc_fraction(forward),
        "gc_reverse_oligo": 100 * gc_fraction(reverse),
        "reported_gc_forward": pair.get("gc_forward"),
        "reported_gc_reverse": pair.get("gc_reverse"),
        "clamp_forward_is_gc": forward[-1] in "GC",
        "clamp_reverse_is_gc": reverse[-1] in "GC",
        "reported_clamp_forward": pair.get("gc_clamp_forward"),
        "reported_clamp_reverse": pair.get("gc_clamp_reverse"),
        "reported_tm_forward": pair.get("tm_forward"),
        "reported_tm_reverse": pair.get("tm_reverse"),
        "biopython_tm_forward_1M": bio(forward, False),
        "biopython_tm_reverse_1M": bio(reverse, False),
        "biopython_tm_forward_default": bio_default(forward),
        "biopython_tm_reverse_default": bio_default(reverse),
        "primer3_tm_forward_1M": primer3.calc_tm(forward, **aligned),
        "primer3_tm_reverse_1M": primer3.calc_tm(reverse, **aligned),
        "primer3_tm_forward_default": primer3.calc_tm(forward),
        "reported_hairpin_forward": pair.get("hairpin_dG_forward"),
        "primer3_hairpin_forward_kcal": hairpin_f.dg / 1000,
        "reported_hairpin_reverse": pair.get("hairpin_dG_reverse"),
        "primer3_hairpin_reverse_kcal": hairpin_r.dg / 1000,
        "reported_self_forward": pair.get("self_dimer_dG_forward"),
        "primer3_self_forward_kcal": self_f.dg / 1000,
        "reported_self_reverse": pair.get("self_dimer_dG_reverse"),
        "primer3_self_reverse_kcal": self_r.dg / 1000,
        "reported_hetero": pair.get("hetero_dimer_dG"),
        "primer3_hetero_kcal": hetero.dg / 1000,
        "primer3_dg_unit": "primer3 dg is cal/mol; values above are divided by 1000",
        "thermo_note": (
            "DNA Design Suite uses SantaLucia 1998 Table 2 at 1 M Na+ and 1 uM "
            "strand, with no divalent ions. Biopython defaults are 50 mM Na+, "
            "Owczarzy saltcorr=5, 25 nM. primer3 defaults add 1.5 mM divalent "
            "and 0.6 mM dNTP. Hairpin loop initiation is not in the bundled "
            "Table 2, so hairpin ΔG is not expected to match primer3."
        ),
    }


def cmd_primers() -> None:
    templates = {}
    for name in ("BBa_E0040_onestop", "lacZ", "yagF", "waaS"):
        templates[name] = read_fa(name)[1]
    specs = [
        ("BBa_E0040_onestop", 200, 400, 80, 150, 400),
        ("lacZ", 400, 650, 80, 180, 400),
        ("yagF", 400, 650, 80, 180, 400),
        ("waaS", 200, 450, 80, 180, 420),
    ]
    successes = []
    for name, start, end, flank, pmin, pmax in specs:
        body = {
            "sequence": templates[name],
            "region": {"start": start, "end": end},
            "config": {
                "tm_target": 60,
                "product_size_range": {"min": pmin, "max": pmax},
                "flank_length": flank,
                "primer_length_range": {"min": 18, "max": 25},
            },
        }
        print("primer", name, flush=True)
        result = primer_job(body)
        status = result.get("status") or {}
        pair = None
        pairs = (status.get("primer_result") or {}).get("pairs") or []
        if pairs:
            pair = pairs[0]
        entry: dict[str, object] = {
            "name": name,
            "region": [start, end],
            "state": status.get("state"),
            "problem": status.get("problem") or result.get("problem"),
            "n_pairs": len(pairs),
            "warnings": status.get("warnings"),
        }
        if pair:
            entry["pair_id"] = pair.get("pair_id")
            entry["checks"] = check_pair(templates[name], pair)
            entry["pair"] = {
                key: pair.get(key)
                for key in (
                    "forward_sequence",
                    "reverse_sequence",
                    "forward_start",
                    "forward_end",
                    "reverse_start",
                    "reverse_end",
                    "product_size",
                    "tm_forward",
                    "tm_reverse",
                    "gc_forward",
                    "gc_reverse",
                    "gc_clamp_forward",
                    "gc_clamp_reverse",
                    "hairpin_dG_forward",
                    "hairpin_dG_reverse",
                    "self_dimer_dG_forward",
                    "self_dimer_dG_reverse",
                    "hetero_dimer_dG",
                    "parameter_citation",
                    "blast",
                )
            }
        successes.append(entry)
        print(" ", entry["state"], entry["n_pairs"], flush=True)

    gfp = templates["BBa_E0040_onestop"]
    edges = {
        "very_short": primer_job(
            {
                "sequence": "ATGAAACCCGGGTTTACTAA",
                "region": {"start": 6, "end": 12},
                "config": {
                    "tm_target": 60,
                    "product_size_range": {"min": 20, "max": 80},
                    "flank_length": 40,
                    "primer_length_range": {"min": 18, "max": 22},
                },
            }
        ),
        "homopolymer_no_primer": primer_job(
            {
                "sequence": (
                    "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
                    "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
                    "GCTAGCTAGCTAGCTAGCTAGCTAGCTAG"
                ),
                "region": {"start": 30, "end": 55},
                "config": {
                    "tm_target": 60,
                    "product_size_range": {"min": 20, "max": 90},
                    "flank_length": 40,
                    "primer_length_range": {"min": 18, "max": 22},
                },
            }
        ),
        "ambiguous_N": primer_job(
            {
                "sequence": gfp[:180] + "N" + gfp[181:],
                "region": {"start": 200, "end": 400},
                "config": {
                    "tm_target": 60,
                    "product_size_range": {"min": 150, "max": 400},
                    "flank_length": 80,
                    "primer_length_range": {"min": 18, "max": 25},
                },
            }
        ),
        "lowercase": primer_job(
            {
                "sequence": gfp.lower(),
                "region": {"start": 200, "end": 400},
                "config": {
                    "tm_target": 60,
                    "product_size_range": {"min": 150, "max": 400},
                    "flank_length": 80,
                    "primer_length_range": {"min": 18, "max": 25},
                },
            }
        ),
        "fasta_header": primer_job(
            {
                "sequence": ">GFP demo\n" + gfp,
                "region": {"start": 200, "end": 400},
                "config": {
                    "tm_target": 60,
                    "product_size_range": {"min": 150, "max": 400},
                    "flank_length": 80,
                    "primer_length_range": {"min": 18, "max": 25},
                },
            }
        ),
    }
    compact = {}
    for name, result in edges.items():
        status = result.get("status") or {}
        problem = status.get("problem") or result.get("problem")
        compact[name] = {
            "http": result.get("http"),
            "state": status.get("state"),
            "problem": problem,
            "n_pairs": len((status.get("primer_result") or {}).get("pairs") or []),
        }
        print("edge", name, compact[name]["state"], compact[name]["http"], flush=True)

    payload = {"meta": meta(), "designs": successes, "edges": compact}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "primers.json").write_text(json.dumps(payload), encoding="utf-8")
    print("wrote", OUT / "primers.json", flush=True)


def main() -> None:
    if len(sys.argv) != 2 or sys.argv[1] not in {"cai", "hard", "primers"}:
        raise SystemExit("usage: run_validation.py cai|hard|primers")
    OUT.mkdir(parents=True, exist_ok=True)
    {"cai": cmd_cai, "hard": cmd_hard, "primers": cmd_primers}[sys.argv[1]]()


if __name__ == "__main__":
    main()
