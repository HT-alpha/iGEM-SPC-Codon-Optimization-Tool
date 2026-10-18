#!/usr/bin/env python3
"""How close DNA Design Suite's search gets to an optimum, and where it stops moving.

Part B of the v2.9.3 design-quality study. Exact numbers are the output of
this process. The unconstrained CAI optimum is the per-codon maximum of the
ribosomal relative adaptiveness (ties break toward the lexicographically
smaller codon). That assignment is also what ``cai_independent.most_frequent``
returns, stops left as written. A chain DP with the same objective and no
patterns has to land on the same CAI.

With restriction sites, homopolymers and fixed motifs, the same DP is exact:
the state is the last 7 nucleotides, which covers NotI (8) and every shorter
pattern. GC windows and scattered 12-mers are not in that state. Their CAI is
an upper bound. If the DP sequence also passes the full checker, it is an
optimum of the full set, because nothing feasible for the full set can score
above the relaxed optimum.

A 54 nt synthetic CDS is enumerated outright, including GC windows and 12-mers.
``2**16`` codon assignments, stop fixed.

From the repository root, after ``make build``:

    .venv/bin/python tools/optimization_quality.py \\
        --out docs/benchmarks/v2.9.3-quality-rows.json
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
import time
from array import array
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "benchmarks" / "validation"))

import cai_independent  # noqa: E402

import benchmark_common as common  # noqa: E402
from dna_design_suite import _core  # noqa: E402
from dna_design_suite.modelreg import ArtifactStore  # noqa: E402

BASES = "ACGT"
BASE_INDEX = {base: index for index, base in enumerate(BASES)}
SUFFIX = 7
NSTATE = 4**SUFFIX
NEG = -1e300

# Hypervolume is minimisation. CAI is negated. The reference is fixed before
# looking at a front so two generation caps stay comparable. A point that does
# not strictly dominate the reference contributes nothing; the row says so.
HV_REFERENCE = {
    ("usage.cai", "burden.growth_cost"): [0.05, 5.0],
    ("usage.cai", "structure.local_accessibility"): [0.05, 8.0],
}


def options_for(weights: dict[str, float]) -> dict[str, list[tuple[str, float]]]:
    grouped: dict[str, list[tuple[str, float]]] = {}
    for codon, amino in cai_independent.AA.items():
        added = common.log_add(codon, weights)
        if added is None:
            continue
        grouped.setdefault(amino, []).append((codon, added))
    for amino in grouped:
        grouped[amino].sort()
    return grouped


def pattern_free(window: str, patterns: tuple[str, ...]) -> bool:
    return all(pattern not in window for pattern in patterns)


def decode(state: int) -> str:
    chars = []
    for shift in range(SUFFIX - 1, -1, -1):
        chars.append(BASES[(state >> (2 * shift)) & 3])
    return "".join(chars)


def encode(suffix: str) -> int:
    state = 0
    for base in suffix:
        state = (state << 2) | BASE_INDEX[base]
    return state


def append_codon(state: int, codon: str) -> int:
    for base in codon:
        state = ((state << 2) | BASE_INDEX[base]) & (NSTATE - 1)
    return state


_TRANSITION_CACHE: dict[tuple[str, ...], dict[str, list[list[tuple[int, float, str]]]]] = {}


def build_transitions(
    grouped: dict[str, list[tuple[str, float]]], patterns: tuple[str, ...]
) -> dict[str, list[list[tuple[int, float, str]]]]:
    """Per amino acid: (codon_id table unused, list over states of legal moves).

    ``moves[state]`` is a list of ``(new_state, log_add, codon)``. Built once
    per pattern set. The suffix is exactly 7 bases, so this table is used
    after the sequence is at least that long.
    """
    decoded = [decode(state) for state in range(NSTATE)]
    moves: dict[str, list[list[tuple[int, float, str]]]] = {}
    for amino, choices in grouped.items():
        per_state: list[list[tuple[int, float, str]]] = [[] for _ in range(NSTATE)]
        for state, suffix in enumerate(decoded):
            legal = []
            for codon, added in choices:
                if pattern_free(suffix + codon, patterns):
                    legal.append((append_codon(state, codon), added, codon))
            per_state[state] = legal
        moves[amino] = per_state
    return moves


def cached_transitions(
    grouped: dict[str, list[tuple[str, float]]], patterns: tuple[str, ...]
) -> dict[str, list[list[tuple[int, float, str]]]]:
    cached = _TRANSITION_CACHE.get(patterns)
    if cached is None:
        cached = build_transitions(grouped, patterns)
        _TRANSITION_CACHE[patterns] = cached
    return cached


def sequence_cai(seq: str, weights: dict[str, float]) -> float | None:
    return cai_independent.score(seq, weights)["cai"]


def exact_local(cds: str, weights: dict[str, float], patterns: list[str]) -> dict[str, Any]:
    """Maximum counted CAI under a regular pattern set.

    Returns one optimal CDS. Ties take the lexicographically smaller codon at
    the first position where two survivors differ, by keeping the smaller
    codon when the log-sum is equal.
    """
    protein = cai_independent.translate(cds)
    grouped = options_for(weights)
    forbidden = tuple(sorted(set(patterns)))
    if any(amino not in grouped for amino in protein):
        return {"feasible": False, "reason": "a counted amino acid has no positive-weight codon"}

    # Grow from the empty prefix until every survivor has a 7-mer suffix.
    layer: dict[str, tuple[float, str]] = {"": (0.0, "")}
    start_prefix: dict[int, tuple[float, str]] | None = None
    consumed = 0
    for amino in protein:
        if layer and min(len(suffix) for suffix in layer) >= SUFFIX:
            start_prefix = {}
            for suffix, (score, prefix) in layer.items():
                state = encode(suffix)
                current = start_prefix.get(state)
                if (
                    current is None
                    or score > current[0]
                    or (score == current[0] and prefix < current[1])
                ):
                    start_prefix[state] = (score, prefix)
            break
        nxt: dict[str, tuple[float, str]] = {}
        for suffix, (score, prefix) in layer.items():
            for codon, added in grouped[amino]:
                candidate = suffix + codon
                if not pattern_free(candidate, forbidden):
                    continue
                new_suffix = candidate[-SUFFIX:]
                new_score = score + added
                new_prefix = prefix + codon
                current = nxt.get(new_suffix)
                if (
                    current is None
                    or new_score > current[0]
                    or (new_score == current[0] and new_prefix < current[1])
                ):
                    nxt[new_suffix] = (new_score, new_prefix)
        layer = nxt
        consumed += 1
        if not layer:
            return {
                "feasible": False,
                "reason": "no sequence survives the patterns",
                "codons_placed": consumed,
            }
    if start_prefix is None:
        # The whole protein fitted in the growing phase.
        score, prefix = min(layer.values(), key=lambda item: (-item[0], item[1]))
        return _witness(prefix, cds, weights, forbidden, score)

    transitions = cached_transitions(grouped, forbidden)
    scores = array("d", [NEG] * NSTATE)
    for state, (score, _prefix) in start_prefix.items():
        scores[state] = score
    parent_state: list[array] = []
    parent_codon: list[array] = []
    for amino in protein[consumed:]:
        new_scores = array("d", [NEG] * NSTATE)
        back_state = array("i", [-1] * NSTATE)
        back_codon = array("I", [0] * NSTATE)
        moves = transitions[amino]
        for state in range(NSTATE):
            score = scores[state]
            if score <= NEG / 2:
                continue
            for new_state, added, codon in moves[state]:
                value = score + added
                current = new_scores[new_state]
                smaller_codon = back_state[new_state] != -1 and codon < _codon_of(
                    int(back_codon[new_state])
                )
                if value > current or (value == current and smaller_codon):
                    new_scores[new_state] = value
                    back_state[new_state] = state
                    back_codon[new_state] = _codon_id(codon)
        scores = new_scores
        parent_state.append(back_state)
        parent_codon.append(back_codon)
        if max(scores) <= NEG / 2:
            return {"feasible": False, "reason": "no sequence survives the patterns"}

    best_state = max(range(NSTATE), key=lambda state: scores[state])
    if scores[best_state] <= NEG / 2:
        return {"feasible": False, "reason": "no sequence survives the patterns"}
    codons: list[str] = []
    state = best_state
    for layer_index in range(len(parent_state) - 1, -1, -1):
        codon = _codon_of(parent_codon[layer_index][state])
        codons.append(codon)
        state = parent_state[layer_index][state]
        if state < 0:
            return {"feasible": False, "reason": "parent pointer missing"}
    codons.reverse()
    prefix = start_prefix[state][1]
    sequence = prefix + "".join(codons)
    return _witness(sequence, cds, weights, forbidden, scores[best_state])


_CODON_LIST = [a + b + c for a in BASES for b in BASES for c in BASES]
_CODON_ID = {codon: index for index, codon in enumerate(_CODON_LIST)}


def _codon_id(codon: str) -> int:
    return _CODON_ID[codon]


def _codon_of(codon_id: int) -> str:
    return _CODON_LIST[codon_id]


def _witness(
    sequence: str,
    source: str,
    weights: dict[str, float],
    forbidden: tuple[str, ...],
    log_sum: float,
) -> dict[str, Any]:
    cai = sequence_cai(sequence, weights)
    counted = cai_independent.score(sequence, weights)["counted_codons"]
    recomputed = None if cai is None or counted == 0 else abs(math.log(cai) * counted - log_sum)
    return {
        "feasible": True,
        "cds": sequence,
        "sha256": common.sha256_text(sequence),
        "cai": cai,
        "log_sum": log_sum,
        "log_sum_abs_delta": recomputed,
        "same_protein": cai_independent.translate(sequence) == cai_independent.translate(source),
        "patterns_clear": pattern_free(sequence, forbidden),
        "gc": common.gc_fraction(sequence),
    }


def brute_optimum(protein: str, weights: dict[str, float], patterns: list[str]) -> dict[str, Any]:
    grouped = options_for(weights)
    forbidden = tuple(patterns)
    best_score = None
    best_seq = None
    seen = 0
    stacks = [grouped[amino] for amino in protein]
    for choice in itertools.product(*[[codon for codon, _added in item] for item in stacks]):
        seen += 1
        sequence = "".join(choice)
        if not pattern_free(sequence, forbidden):
            continue
        cai = sequence_cai(sequence, weights)
        if (
            best_score is None
            or (cai is not None and cai > best_score)
            or (cai == best_score and (best_seq is None or sequence < best_seq))
        ):
            best_score = cai
            best_seq = sequence
    return {"seen": seen, "cai": best_score, "cds": best_seq}


def self_check(weights: dict[str, float]) -> dict[str, Any]:
    """The DP has to match exhaustive search on a protein small enough to list."""
    protein_cds = "ATG" + "CAT" + "AAT" + "TTT" + "TAT" + "AAA" + "GAT" + "GAA" + "TAA"
    # M H N F Y K D E *
    protein = cai_independent.translate(protein_cds)
    greedy = cai_independent.most_frequent(protein_cds, weights)
    planted = greedy[:6]
    patterns = [planted, "AAAAAAA", "CCCCCCC", "GGGGGGG", "TTTTTTT"]
    brute = brute_optimum(protein, weights, patterns)
    solved = exact_local(protein_cds, weights, patterns)
    relaxed = exact_local(protein_cds, weights, [])
    relaxed_cai = sequence_cai(greedy, weights)
    return {
        "protein": protein,
        "enumerations": brute["seen"],
        "planted_pattern": planted,
        "brute_cai": brute["cai"],
        "dp_cai": solved.get("cai"),
        "dp_matches_brute": (
            solved.get("cai") == brute["cai"] and solved.get("patterns_clear") is True
        ),
        "dp_same_protein": solved.get("same_protein"),
        "unconstrained_dp_cai": relaxed.get("cai"),
        "most_frequent_cai": relaxed_cai,
        "unconstrained_dp_matches_closed_form": relaxed.get("cai") == relaxed_cai,
        "planted_pattern_binds": (
            brute["cai"] is not None and relaxed_cai is not None and brute["cai"] < relaxed_cai
        ),
    }


def enumerate_synthetic(weights: dict[str, float], sets: list[dict[str, Any]]) -> dict[str, Any]:
    """54 nt, 16 binary codon choices, Met and the stop fixed.

    The protein is M FYHQNKDE FYHQNKDE *. The two-codon families are the ones
    whose high-GC synonym can sit near GC 0.5, so the 0.40-0.60 window is not
    an empty set by construction.
    """
    order = "FYHQNKDE"
    pairs = {
        "F": ("TTT", "TTC"),
        "Y": ("TAT", "TAC"),
        "H": ("CAT", "CAC"),
        "Q": ("CAA", "CAG"),
        "N": ("AAT", "AAC"),
        "K": ("AAA", "AAG"),
        "D": ("GAT", "GAC"),
        "E": ("GAA", "GAG"),
    }
    amino = "M" + order + order
    # Bit i selects the second codon of amino acid i+1 (the variable ones).
    variable = list(amino[1:])
    assert len(variable) == 16
    assert all(item in pairs for item in variable)
    low = "ATG" + "".join(pairs[item][0] for item in variable) + "TAA"
    results = []
    assignments = 1 << 16
    sequences = []
    cais = []
    for mask in range(assignments):
        codons = ["ATG"]
        for index, item in enumerate(variable):
            codons.append(pairs[item][(mask >> index) & 1])
        codons.append("TAA")
        sequence = "".join(codons)
        sequences.append(sequence)
        cais.append(sequence_cai(sequence, weights))
    for constraint in sets:
        document = constraint["document"]
        feasible_flags: list[bool] = []
        batch = 4096
        for start in range(0, assignments, batch):
            chunk = sequences[start : start + batch]
            flags = _core.constraint_feasible_batch(
                json.dumps(document),
                chunk,
                [(0, len(chunk[0]))],
                "assembled",
            )
            feasible_flags.extend(bool(flag) for flag in flags)
        best = None
        best_seq = None
        n_feasible = 0
        for sequence, cai, flag in zip(sequences, cais, feasible_flags, strict=True):
            if not flag or cai is None:
                continue
            n_feasible += 1
            if best is None or cai > best or (cai == best and sequence < best_seq):
                best = cai
                best_seq = sequence
        confirmed = None
        if best_seq is not None:
            confirmed = common.violation_summary(best_seq, document)["feasible"]
            if confirmed is not True:
                raise SystemExit(
                    f"enumeration and constraint_violations disagree on {constraint['name']}"
                )
        results.append(
            {
                "constraint_set": constraint["name"],
                "assignments": assignments,
                "n_feasible": n_feasible,
                "exact_cai": best,
                "exact_sha256": None if best_seq is None else common.sha256_text(best_seq),
                "exact_cds": best_seq,
                "exact_gc": None if best_seq is None else common.gc_fraction(best_seq),
                "violation_report_agrees": confirmed,
            }
        )
    return {
        "protein": amino + "*",
        "nt": 54,
        "input_cds": low,
        "input_cai": sequence_cai(low, weights),
        "input_gc": common.gc_fraction(low),
        "note": (
            "Every variable amino acid has two codons. The input uses the "
            "lexicographically smaller one, which is the AT-rich synonym. "
            "Met is ATG and the stop is TAA. 65536 sequences, all of them."
        ),
        "by_constraint": results,
    }


def dominates(left: list[float], right: list[float], senses: list[str]) -> bool:
    better = False
    for x, y, sense in zip(left, right, senses, strict=True):
        if sense == "maximize":
            if x < y:
                return False
            if x > y:
                better = True
        else:
            if x > y:
                return False
            if x < y:
                better = True
    return better


def pareto_report(run: dict[str, Any]) -> dict[str, Any]:
    senses = list(run.get("objective_sense") or [])
    complete = []
    for candidate in run.get("candidates") or []:
        vector = candidate.get("objective_vector")
        if not vector or any(value is None for value in vector):
            continue
        if candidate.get("feasible_flag") is False:
            continue
        complete.append([float(value) for value in vector])
    dominated_pairs = 0
    for left, right in itertools.combinations(complete, 2):
        if dominates(left, right, senses) or dominates(right, left, senses):
            dominated_pairs += 1
    spread = []
    if complete and senses:
        for index, sense in enumerate(senses):
            column = [vector[index] for vector in complete]
            spread.append(
                {
                    "sense": sense,
                    "min": min(column),
                    "max": max(column),
                    "range": max(column) - min(column),
                }
            )
    reference = HV_REFERENCE.get(tuple(run.get("objectives") or []))
    hypervolume = None
    outside = None
    if reference is not None and complete and len(reference) == len(senses):
        minimized = []
        for vector in complete:
            minimized.append(
                [
                    -value if sense == "maximize" else value
                    for value, sense in zip(vector, senses, strict=True)
                ]
            )
        outside = sum(
            1
            for point in minimized
            if any(value >= bound for value, bound in zip(point, reference, strict=True))
        )
        try:
            hypervolume = float(_core.indicator_hypervolume(minimized, reference))
        except ValueError as exc:
            hypervolume = {"error": str(exc)}
    return {
        "n_candidates": len(run.get("candidates") or []),
        "n_feasible_complete": len(complete),
        "dominated_pairs": dominated_pairs,
        "nondominated": dominated_pairs == 0,
        "spread": spread,
        "hypervolume_reference_minimized": reference,
        "points_outside_reference": outside,
        "hypervolume": hypervolume,
        "hypervolume_history": run.get("hypervolume_history"),
        "recommended_rule": run.get("recommended_rule"),
        "objective_sense": senses,
    }


def tier_request(
    cds: str,
    constraints: dict[str, Any],
    profile: dict[str, Any],
    digest: str,
    seed: int,
) -> str:
    return json.dumps(
        {
            "cds": cds,
            "five_prime_flank": "",
            "three_prime_flank": "",
            "translation_table": common.TRANSLATION_TABLE,
            "constraint_set": constraints,
            "terms": ["usage.cai"],
            "weights": [1.0],
            "host_profile": profile,
            "host_id": common.HOST,
            "host_digest": digest,
            "root_seed": seed,
        }
    )


def compact_run(
    run: dict[str, Any], weights: dict[str, float], constraints: dict[str, Any]
) -> dict[str, Any]:
    seq = run.get("cds")
    cai = None
    checker = None
    if seq:
        cai = sequence_cai(seq, weights)
        checker = common.violation_summary(seq, constraints)
    origins = [item["origin"] for item in run.get("archive") or []]
    return {
        "state": run.get("state"),
        "wall_s": run.get("wall_s"),
        "budget": run.get("budget"),
        "seed": run.get("seed"),
        "generations_requested": run.get("generations_requested"),
        "stagnation_window": run.get("stagnation_window"),
        "objectives": run.get("objectives"),
        "tier1_beam": run.get("tier1_beam"),
        "tier1_seeds": run.get("tier1_seeds"),
        "tier2_block": run.get("tier2_block"),
        "tier2_sweeps": run.get("tier2_sweeps"),
        "population": run.get("population"),
        "warning_codes": [item.get("code") for item in run.get("warnings") or []],
        "problem": (
            (run.get("problem") or {}).get("code")
            if isinstance(run.get("problem"), dict)
            else run.get("problem")
        ),
        "cai": cai,
        "sha256": None if not seq else common.sha256_text(seq),
        "gc": None if not seq else common.gc_fraction(seq),
        "feasible": None if checker is None else checker["feasible"],
        "violations": None if checker is None else checker["actionable_by_id"],
        "n_candidates": len(run.get("candidates") or []),
        "origins": origins,
        "objective_vector": (
            (run.get("candidates") or [{}])[0].get("objective_vector")
            if run.get("candidates")
            else None
        ),
        "pareto": pareto_report(run) if run.get("candidates") else None,
        "progress": run.get("progress"),
    }


def gap(exact: float | None, achieved: float | None) -> dict[str, float | None]:
    if exact is None or achieved is None:
        return {"exact": exact, "achieved": achieved, "gap": None}
    return {"exact": exact, "achieved": achieved, "gap": exact - achieved}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--self-check-only", action="store_true")
    args = parser.parse_args()

    weights = common.load_weights()
    profile = common.host_profile(weights)
    digest = common.profile_digest(profile)
    genes = common.load_genes()
    sets = {item["name"]: item for item in common.load_constraint_sets()}
    assembly = sets["assembly"]["document"]
    print("self-check", flush=True)
    checked = self_check(weights)
    print(json.dumps({key: checked[key] for key in checked if key != "protein"}), flush=True)
    if not checked["dp_matches_brute"] or not checked["unconstrained_dp_matches_closed_form"]:
        raise SystemExit("exact DP does not match exhaustive search")
    if args.self_check_only:
        return

    store = ArtifactStore(Path("/tmp/cf-v293-quality-store"))
    by_gene = {gene["gene"]: gene for gene in genes}

    print("unconstrained gap", flush=True)
    unconstrained = []
    for gene in genes:
        closed = cai_independent.most_frequent(gene["cds"], weights)
        closed_cai = sequence_cai(closed, weights)
        # No cross-codon term: the closed form is the optimum. The DP is
        # checked against it on the self-check protein and on short genes.
        dp_cai = None
        dp_match = None
        if gene["nt"] <= 400:
            dp = exact_local(gene["cds"], weights, [])
            dp_cai = dp.get("cai")
            dp_match = dp_cai == closed_cai
        designed = common.run_codonforge(
            gene["cds"],
            common.empty_constraints(),
            profile,
            digest,
            store,
            budget="fast",
            seed=7,
            generations=0,
            stagnation_window=21,
            job_id=f"unc-{gene['gene']}",
        )
        achieved = sequence_cai(designed["cds"], weights) if designed.get("cds") else None
        unconstrained.append(
            {
                "gene": gene["gene"],
                "nt": gene["nt"],
                "closed_form_cai": closed_cai,
                "dp_cai": dp_cai,
                "dp_matches_closed_form": dp_match,
                "dp_note": (
                    "Ran. Empty pattern set."
                    if gene["nt"] <= 400
                    else (
                        "Not run. Unconstrained CAI separates by codon, so the closed form "
                        "is exact; the DP is checked on every gene of at most 400 nt and "
                        "on the self-check protein."
                    )
                ),
                "codonforge_fast_generations_0": compact_run(
                    designed, weights, common.empty_constraints()
                ),
                "gap_fast": gap(closed_cai, achieved),
            }
        )
        print(f"  {gene['gene']} closed {closed_cai} achieved {achieved}", flush=True)

    print("local DP", flush=True)
    local_rows = []
    # lacZ and yeeJ are exact on the unconstrained problem above. The pattern
    # DP's state is 4**7; those lengths are still finite, and a slow gene is
    # recorded rather than skipped without a time.
    for gene in genes:
        if gene["nt"] > 2500:
            local_rows.append(
                {
                    "gene": gene["gene"],
                    "nt": gene["nt"],
                    "skipped": True,
                    "reason": (
                        "pattern DP not run above 2500 nt in this study; "
                        "unconstrained gap is in the other table"
                    ),
                }
            )
            continue
        for name in ("golden_gate", "biobrick_rfc10", "assembly"):
            document = sets[name]["document"]
            described = common.local_patterns(document)
            print(f"  dp {gene['gene']} {name}", flush=True)
            started = time.perf_counter()
            solved = exact_local(gene["cds"], weights, described["patterns"])
            elapsed = time.perf_counter() - started
            full = None
            if solved.get("cds"):
                full = common.violation_summary(solved["cds"], document)
            local_rows.append(
                {
                    "gene": gene["gene"],
                    "nt": gene["nt"],
                    "constraint_set": name,
                    "deferred_constraint_ids": described["deferred_constraint_ids"],
                    "exact_for_full_set": described["exact_for_full_set"],
                    "dp_wall_s": elapsed,
                    "dp_feasible": solved.get("feasible"),
                    "dp_cai": solved.get("cai"),
                    "dp_same_protein": solved.get("same_protein"),
                    "dp_patterns_clear": solved.get("patterns_clear"),
                    "dp_log_sum_abs_delta": solved.get("log_sum_abs_delta"),
                    "full_checker_on_dp_sequence": full,
                    "dp_sequence_is_full_optimum": bool(
                        solved.get("feasible") and full and full["feasible"]
                    ),
                    "reason": solved.get("reason"),
                }
            )

    print("synthetic enumeration", flush=True)
    synthetic = enumerate_synthetic(weights, list(sets.values()))
    synthetic_designs = []
    for name in ("assembly", "golden_gate", "biobrick_rfc10"):
        for generations, budget in ((0, "fast"), (4, "standard"), (12, "standard")):
            print(f"  synthetic {name} {budget} gen {generations}", flush=True)
            run = common.run_codonforge(
                synthetic["input_cds"],
                sets[name]["document"],
                profile,
                digest,
                store,
                budget=budget,
                seed=7,
                generations=generations,
                stagnation_window=21,
                job_id=f"syn-{name}-{budget}-{generations}",
            )
            compact = compact_run(run, weights, sets[name]["document"])
            exact = next(
                item for item in synthetic["by_constraint"] if item["constraint_set"] == name
            )
            compact["gap_to_enumeration"] = gap(exact["exact_cai"], compact["cai"])
            synthetic_designs.append(compact | {"constraint_set": name})

    print("seeds", flush=True)
    seed_rows = []
    for gene_name in ("rpsL", "elastin_vpgvg20", "waaS", "yagF", "mcherry_BBa_J06504"):
        gene = by_gene[gene_name]
        for seed in (1, 2, 7, 11, 42):
            print(f"  seed {gene_name} {seed}", flush=True)
            run = common.run_codonforge(
                gene["cds"],
                assembly,
                profile,
                digest,
                store,
                budget="standard",
                seed=seed,
                generations=4,
                stagnation_window=21,
                job_id=f"seed-{gene_name}-{seed}",
            )
            seed_rows.append(
                compact_run(run, weights, assembly)
                | {"gene": gene_name, "constraint_set": "assembly"}
            )
        for seed in (1, 7, 42):
            run = common.run_codonforge(
                gene["cds"],
                common.empty_constraints(),
                profile,
                digest,
                store,
                budget="fast",
                seed=seed,
                generations=0,
                stagnation_window=21,
                job_id=f"seed-open-{gene_name}-{seed}",
            )
            seed_rows.append(
                compact_run(run, weights, common.empty_constraints())
                | {"gene": gene_name, "constraint_set": "empty"}
            )

    print("generations", flush=True)
    generation_rows = []
    generation_plan = (
        ("fast", 0, 21),
        ("standard", 0, 21),
        ("standard", 4, 21),
        ("standard", 12, 21),
        ("standard", 12, 8),
    )
    for gene_name in ("rpsL", "elastin_vpgvg20", "waaS", "yagF", "yeeJ"):
        gene = by_gene[gene_name]
        for budget, generations, window in generation_plan:
            print(f"  gen {gene_name} {budget} {generations} window {window}", flush=True)
            run = common.run_codonforge(
                gene["cds"],
                assembly,
                profile,
                digest,
                store,
                budget=budget,
                seed=7,
                generations=generations,
                stagnation_window=window,
                job_id=f"gen-{gene_name}-{budget}-{generations}-{window}",
            )
            generation_rows.append(
                compact_run(run, weights, assembly)
                | {"gene": gene_name, "constraint_set": "assembly"}
            )

    print("pareto", flush=True)
    pareto_rows = []
    for gene_name in ("rpsL", "waaS", "elastin_vpgvg20"):
        gene = by_gene[gene_name]
        for objectives in (
            ["usage.cai", "burden.growth_cost"],
            ["usage.cai", "structure.local_accessibility"],
        ):
            for generations in (0, 4, 12):
                print(f"  pareto {gene_name} {objectives} gen {generations}", flush=True)
                run = common.run_codonforge(
                    gene["cds"],
                    assembly,
                    profile,
                    digest,
                    store,
                    objectives=objectives,
                    budget="standard",
                    seed=7,
                    generations=generations,
                    stagnation_window=21,
                    job_id=f"pareto-{gene_name}-{generations}-{objectives[-1].split('.')[-1]}",
                )
                pareto_rows.append(
                    compact_run(run, weights, assembly)
                    | {"gene": gene_name, "constraint_set": "assembly"}
                )

    print("repair bias", flush=True)
    repair_rows = []
    repair_genes = (
        "rpsL",
        "elastin_vpgvg20",
        "gfp_BBa_E0040",
        "waaS",
        "yagF",
        "mcherry_BBa_J06504",
    )
    for gene_name in repair_genes:
        gene = by_gene[gene_name]
        greedy = cai_independent.most_frequent(gene["cds"], weights)
        for label, start in (("native", gene["cds"]), ("cai_greedy", greedy)):
            print(f"  repair {gene_name} from {label}", flush=True)
            started = time.perf_counter()
            report = json.loads(
                _core.resolve_feasible_seed(tier_request(start, assembly, profile, digest, 7))
            )
            repair_wall = time.perf_counter() - started
            repaired = report.get("cds")
            repair_rows.append(
                {
                    "gene": gene_name,
                    "start": label,
                    "repair_wall_s": repair_wall,
                    "repair_feasible": report.get("feasible"),
                    "repair_edits": report.get("edits"),
                    "repair_cycles": report.get("cycles"),
                    "repair_cai": None if not repaired else sequence_cai(repaired, weights),
                    "repair_sha256": None if not repaired else common.sha256_text(repaired),
                    "repair_checker": None
                    if not repaired
                    else common.violation_summary(repaired, assembly),
                    "residual": report.get("residual"),
                }
            )
            run = common.run_codonforge(
                start,
                assembly,
                profile,
                digest,
                store,
                budget="standard",
                seed=7,
                generations=4,
                stagnation_window=21,
                job_id=f"repair-design-{gene_name}-{label}",
            )
            repair_rows[-1]["design"] = compact_run(run, weights, assembly)

    # Attach codonforge gaps to the local DP rows where a matching design exists.
    design_index = {
        (
            row["gene"],
            row["constraint_set"],
            row["budget"],
            row["generations_requested"],
            row["seed"],
        ): row
        for row in generation_rows
    }
    for row in local_rows:
        if row.get("skipped"):
            continue
        match = design_index.get((row["gene"], row["constraint_set"], "standard", 4, 7))
        if match is None:
            continue
        if row.get("dp_sequence_is_full_optimum"):
            row["gap_standard_gen4"] = gap(row.get("dp_cai"), match["cai"])
        else:
            row["relaxed_upper_bound_minus_standard_gen4"] = gap(row.get("dp_cai"), match["cai"])

    document = {
        "study": "v2.9.3 optimization quality",
        "versions": common.tool_versions(),
        "weights": str(common.WEIGHTS_PATH.relative_to(common.ROOT)),
        "host": common.HOST,
        "commands": {
            "this_script": (
                ".venv/bin/python tools/optimization_quality.py "
                "--out docs/benchmarks/v2.9.3-quality-rows.json"
            )
        },
        "self_check": checked,
        "unconstrained": unconstrained,
        "local_dp": local_rows,
        "synthetic_enumeration": synthetic,
        "synthetic_designs": synthetic_designs,
        "seeds": seed_rows,
        "generations": generation_rows,
        "pareto": pareto_rows,
        "repair": repair_rows,
        "search_defaults_from_the_runs": {
            "note": (
                "Taken from the job inputs of these runs, not from the planning "
                "document. Tier 1 beam, tier 2 sweeps and the population are the "
                "JobInputs defaults unless a row says otherwise."
            )
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
