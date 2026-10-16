"""Sharp & Li CAI without importing CodonForge.

Counted codons exclude stops and the single-codon families Met (ATG) and Trp
(TGG). A counted codon whose weight is 0 makes the geometric mean 0. Relative
adaptiveness at 1, including ties, is a top synonym.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

# NCBI translation table 11 elongation codes. For CDS that start with ATG they
# match table 1 at every codon, including the terminal stop.
_TABLE = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)
_CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
AA = dict(zip(_CODONS, _TABLE, strict=True))
STOPS = {codon for codon, aa in AA.items() if aa == "*"}
SINGLE = {codon for codon, aa in AA.items() if aa in {"M", "W"}}


def translate(seq: str) -> str:
    if len(seq) % 3:
        raise ValueError(f"length {len(seq)} is not a multiple of 3")
    return "".join(AA[seq[i : i + 3]] for i in range(0, len(seq), 3))


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def load_ribosomal_weights(path: Path) -> dict[str, float]:
    import json

    document = json.loads(path.read_text(encoding="utf-8"))
    table = document["codon_usage"]["table"]
    weights = {codon: float(table[codon]["relative_adaptiveness"]) for codon in _CODONS}
    return weights


def load_kazusa_weights(path: Path) -> dict[str, float]:
    """Kazusa relative frequencies, converted to w = f / max synonymous f.

    The file stores genome-wide fractions (Nakamura, Gojobori, Ikemura, Nucleic
    Acids Research 2000;28:292), not Sharp & Li highly-expressed weights.
    """
    by_aa: dict[str, list[tuple[str, float]]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            codon = row["codon"].replace("U", "T")
            by_aa.setdefault(row["amino_acid"], []).append(
                (codon, float(row["relative_frequency"]))
            )
    weights: dict[str, float] = {}
    for pairs in by_aa.values():
        peak = max(freq for _, freq in pairs)
        if peak <= 0:
            raise ValueError("Kazusa family has no positive frequency")
        for codon, freq in pairs:
            weights[codon] = freq / peak
    missing = [codon for codon in _CODONS if codon not in weights]
    if missing:
        raise ValueError(f"Kazusa table is missing {missing}")
    return weights


def score(seq: str, weights: dict[str, float]) -> dict[str, float | int | bool | None]:
    """Geometric mean of w over counted codons.

    `cai_is_one_iff_every_counted_codon_is_top` is true when every counted
    codon has w == 1. A zero weight short-circuits the mean to 0.
    """
    logs: list[float] = []
    counted = 0
    non_top = 0
    zero_codons: list[str] = []
    families: dict[str, float] = {}
    for aa in {item for item in AA.values() if item not in {"*", "M", "W"}}:
        families[aa] = max(weights[codon] for codon, amino in AA.items() if amino == aa)
    for index in range(0, len(seq), 3):
        codon = seq[index : index + 3]
        amino = AA[codon]
        if amino in {"*", "M", "W"}:
            continue
        weight = weights[codon]
        counted += 1
        if weight + 1e-12 < families[amino]:
            non_top += 1
        if weight == 0.0:
            zero_codons.append(codon)
        else:
            logs.append(math.log(weight))
    if counted == 0:
        cai = None
    elif zero_codons:
        cai = 0.0
    else:
        cai = math.exp(sum(logs) / len(logs))
    all_top = counted > 0 and non_top == 0 and not zero_codons
    return {
        "cai": cai,
        "counted_codons": counted,
        "non_top_codons": non_top,
        "zero_weight_codons": len(zero_codons),
        "cai_is_one_iff_every_counted_codon_is_top": bool(cai == 1.0) == all_top,
        "all_counted_codons_are_top": all_top,
    }


def score_with_floor(seq: str, weights: dict[str, float]) -> float | None:
    """Sensitivity only: replace 0 with the smallest positive weight.

    This is not Sharp & Li. It shows what a floor would have reported.
    """
    positive = [weight for weight in weights.values() if weight > 0]
    floor = min(positive)
    floored = {codon: (weight if weight > 0 else floor) for codon, weight in weights.items()}
    return score(seq, floored)["cai"]


def most_frequent(seq: str, weights: dict[str, float]) -> str:
    """Highest weight synonym. Ties break toward the lexicographically smaller codon.

    The terminal stop, and any internal stop, is left as written.
    """
    choice: dict[str, str] = {}
    for amino in {item for item in AA.values() if item != "*"}:
        options = [codon for codon, aa in AA.items() if aa == amino]
        choice[amino] = min(options, key=lambda codon: (-weights[codon], codon))
    out = []
    for index in range(0, len(seq), 3):
        codon = seq[index : index + 3]
        amino = AA[codon]
        out.append(codon if amino == "*" else choice[amino])
    return "".join(out)


def gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    return sum(base in "GC" for base in seq) / len(seq)


SITES = {
    "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC",
    "EcoRI": "GAATTC",
    "XbaI": "TCTAGA",
    "SpeI": "ACTAGT",
    "PstI": "CTGCAG",
    "NotI": "GCGGCCGC",
}


def _hits(seq: str, pattern: str) -> int:
    count = 0
    start = 0
    while True:
        found = seq.find(pattern, start)
        if found < 0:
            return count
        count += 1
        start = found + 1


def constraint_census(seq: str) -> dict[str, object]:
    """Same detectable rules as constraints_assembly.json, counted independently.

    GC windows slide by 1. A palindromic site is counted once. Homopolymer
    violations are runs longer than 6. A 12-mer is an offender when it occurs
    more than once, overlaps included. The SD-plus-ATG count is not a product
    rule; it is recorded beside the census.
    """
    windows = []
    if len(seq) >= 50:
        for start in range(0, len(seq) - 50 + 1):
            fraction = gc_fraction(seq[start : start + 50])
            if fraction < 0.4 or fraction > 0.6:
                windows.append(fraction)
    sites = {}
    for enzyme, site in SITES.items():
        rc = revcomp(site)
        sites[enzyme] = _hits(seq, site) + (0 if rc == site else _hits(seq, rc))
    runs = []
    index = 0
    while index < len(seq):
        end = index + 1
        while end < len(seq) and seq[end] == seq[index]:
            end += 1
        if end - index > 6:
            runs.append({"base": seq[index], "start": index, "length": end - index})
        index = end
    kmers: dict[str, int] = {}
    if len(seq) >= 12:
        for start in range(0, len(seq) - 12 + 1):
            kmer = seq[start : start + 12]
            kmers[kmer] = kmers.get(kmer, 0) + 1
    offenders = {kmer: count for kmer, count in kmers.items() if count > 1}
    cryptic = 0
    motif = "AGGAGG"
    start = 0
    while True:
        found = seq.find(motif, start)
        if found < 0:
            break
        after = found + len(motif)
        for spacer in range(4, 13):
            if seq[after + spacer : after + spacer + 3] == "ATG":
                cryptic += 1
                break
        start = found + 1
    return {
        "length_nt": len(seq),
        "gc_global": gc_fraction(seq),
        "gc_windows_outside_0.4_0.6": len(windows),
        "gc_window_min": min(windows) if windows else None,
        "gc_window_max": max(windows) if windows else None,
        "restriction_sites": sites,
        "homopolymer_runs_longer_than_6": len(runs),
        "longest_homopolymer": max((run["length"] for run in runs), default=0),
        "sd_aggagg_both_strands": _hits(seq, "AGGAGG") + _hits(seq, revcomp("AGGAGG")),
        "kmer12_over_limit": len(offenders),
        "kmer12_max_count": max(offenders.values(), default=0),
        "cryptic_sd_atg_sense_not_a_product_rule": cryptic,
    }
