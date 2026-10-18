"""Map a coding sequence onto synthesis-rule classes vendors describe in public.

The numeric cutoffs in the bundled sets are this repository's synthesis and
assembly presets (`examples/synthesis.json`,
`benchmarks/validation/constraints_assembly.json`). They are not GeneArt,
IDT, GenScript, or Twist acceptance thresholds. Those companies' sites were
not opened for this file, and their scoring weights are not published in a
form this screen can replay.

What is cited:

* Raab, Graf, Notka, Schödl and Wagner, Syst Synth Biol 2010, 4:215–225,
  doi:10.1007/s11693-010-9062-3. GeneOptimizer's published parameter classes
  are GC content, repeats, restriction sites and sequence motifs.
* The comparison note records that IDT, GenScript and Twist describe GC,
  repeats and (for structure) hairpins. Hairpins are not scored here.
  Windowed accessibility and global folding are separate, labelled reports.

Wet-lab expression is out of scope. Nothing in the report is an expression
measurement.

`DEP-API-005`.
"""

from __future__ import annotations

from typing import Any

# Preset numbers already used by this repository. Not a manufacturer cutoff.
_SYNTHESIS_GC = (50, 0.40, 0.60)
_ASSEMBLY_SITES = (
    ("EcoRI", "GAATTC"),
    ("XbaI", "TCTAGA"),
    ("SpeI", "ACTAGT"),
    ("PstI", "CTGCAG"),
    ("NotI", "GCGGCCGC"),
    ("BsaI", "GGTCTC"),
    ("BsmBI", "CGTCTC"),
)
_SHINE_DALGARNO = "AGGAGG"

_VENDOR_CLASSES = {
    "gc_window": (
        "GeneArt / GeneOptimizer (Raab et al. 2010): GC content. "
        "IDT, GenScript and Twist public descriptions also name GC. "
        "The numbers below are the DNA Design Suite synthesis preset, not those vendors' cutoffs."
    ),
    "homopolymer": (
        "IDT, GenScript and Twist public descriptions name homopolymer runs. "
        "No vendor cutoff was transcribed. The cap is the DNA Design Suite synthesis preset "
        "(a run longer than 6)."
    ),
    "repeat": (
        "GeneOptimizer (Raab et al. 2010) names repeats. IDT, GenScript and Twist "
        "public descriptions name repeats. The 12-mer cap is the DNA Design Suite "
        "synthesis preset, not a vendor cutoff."
    ),
    "restriction_site": (
        "GeneOptimizer (Raab et al. 2010) names restriction sites. The sites are the "
        "assembly comparison set in this repository (EcoRI, XbaI, SpeI, PstI, NotI, BsaI, BsmBI)."
    ),
    "shine_dalgarno": (
        "The assembly comparison forbids AGGAGG on both strands. That motif is a "
        "ribosome-binding pattern used as a hard constraint in this repository. "
        "It is not a claim about expression."
    ),
}


def screen_sequence(sequence: str, *, sets: list[str] | None = None) -> dict[str, Any]:
    """Flag synthesis-rule classes on one DNA sequence.

    `sets` chooses `dds_synthesis`, `dds_assembly`, or both. The default is both.
    """
    dna = "".join(sequence.upper().split()).replace("U", "T")
    if not dna or any(base not in "ACGT" for base in dna):
        raise ValueError("vendor screen needs an A/C/G/T sequence")
    chosen = sets or ["dds_synthesis", "dds_assembly"]
    unknown = [name for name in chosen if name not in {"dds_synthesis", "dds_assembly"}]
    if unknown:
        raise ValueError(
            f"unknown vendor-screen set {unknown}; known sets are dds_synthesis and dds_assembly"
        )
    rows: list[dict[str, Any]] = []
    if "dds_synthesis" in chosen:
        window, low, high = _SYNTHESIS_GC
        outside = _gc_windows(dna, window, low, high)
        rows.append(
            _row(
                "dds_synthesis",
                "gc_window",
                "c_gc50",
                outside == 0,
                outside,
                (f"{outside} window(s) of {window} nt fall outside GC {low:.2f}–{high:.2f}"),
            )
        )
        runs = _homopolymer_runs(dna, 6)
        rows.append(
            _row(
                "dds_synthesis",
                "homopolymer",
                "c_homopolymer",
                runs == 0,
                runs,
                (f"{runs} homopolymer run(s) longer than 6"),
            )
        )
        repeats = _repeated_kmers(dna, 12)
        rows.append(
            _row(
                "dds_synthesis",
                "repeat",
                "c_kmer12",
                repeats == 0,
                repeats,
                (f"{repeats} distinct 12-mer(s) occur more than once"),
            )
        )
    if "dds_assembly" in chosen:
        window, low, high = _SYNTHESIS_GC
        outside = _gc_windows(dna, window, low, high)
        rows.append(
            _row(
                "dds_assembly",
                "gc_window",
                "c_gc50",
                outside == 0,
                outside,
                (f"{outside} window(s) of {window} nt fall outside GC {low:.2f}–{high:.2f}"),
            )
        )
        runs = _homopolymer_runs(dna, 6)
        rows.append(
            _row(
                "dds_assembly",
                "homopolymer",
                "c_homopolymer",
                runs == 0,
                runs,
                (f"{runs} homopolymer run(s) longer than 6"),
            )
        )
        site_hits = _sites(dna)
        rows.append(
            _row(
                "dds_assembly",
                "restriction_site",
                "c_sites",
                site_hits == 0,
                site_hits,
                f"{site_hits} restriction-site hit(s) on either strand",
            )
        )
        sd = _motif_hits(dna, _SHINE_DALGARNO)
        rows.append(
            _row(
                "dds_assembly",
                "shine_dalgarno",
                "c_sd",
                sd == 0,
                sd,
                f"{sd} AGGAGG hit(s) on either strand",
            )
        )
        repeats = _repeated_kmers(dna, 12)
        rows.append(
            _row(
                "dds_assembly",
                "repeat",
                "c_kmer12",
                repeats == 0,
                repeats,
                (f"{repeats} distinct 12-mer(s) occur more than once"),
            )
        )
    return {
        "schema_version": "1.0.0",
        "kind": "vendor_screen",
        "vendor_score": False,
        "wet_lab_expression": "out_of_scope",
        "length_nt": len(dna),
        "note": (
            "Dry run against DNA Design Suite synthesis and assembly presets, "
            "labelled with the vendor rule classes those presets correspond to. "
            "This is not a GeneArt, IDT, GenScript, or Twist score, and it is "
            "not an expression assay. A gene still has to pass the manufacturer's own screen."
        ),
        "hairpin": {
            "scored": False,
            "reason": (
                "Hairpins are not part of this screen. structure.local_accessibility "
                "is windowed ViennaRNA accessibility. structure global fold is a separate "
                "labelled report and is not a vendor hairpin cutoff."
            ),
        },
        "checks": rows,
        "n_flags": sum(1 for row in rows if not row["passed"]),
    }


def _row(
    rule_set: str,
    kind: str,
    check_id: str,
    passed: bool,
    count: int,
    detail: str,
) -> dict[str, Any]:
    return {
        "set": rule_set,
        "id": check_id,
        "kind": kind,
        "passed": passed,
        "count": count,
        "detail": detail,
        "vendor_cutoff_transcribed": False,
        "numeric_cutoff_source": "dna_design_suite preset",
        "vendor_class": _VENDOR_CLASSES[kind],
    }


def _gc_windows(dna: str, window: int, low: float, high: float) -> int:
    if len(dna) < window:
        fraction = (dna.count("G") + dna.count("C")) / len(dna)
        return int(fraction < low or fraction > high)
    outside = 0
    for start in range(0, len(dna) - window + 1):
        piece = dna[start : start + window]
        fraction = (piece.count("G") + piece.count("C")) / window
        if fraction < low or fraction > high:
            outside += 1
    return outside


def _homopolymer_runs(dna: str, max_run: int) -> int:
    runs = 0
    index = 0
    while index < len(dna):
        base = dna[index]
        end = index + 1
        while end < len(dna) and dna[end] == base:
            end += 1
        if end - index > max_run:
            runs += 1
        index = end
    return runs


def _repeated_kmers(dna: str, k: int) -> int:
    if len(dna) < k:
        return 0
    seen: dict[str, int] = {}
    for start in range(0, len(dna) - k + 1):
        word = dna[start : start + k]
        seen[word] = seen.get(word, 0) + 1
    return sum(1 for count in seen.values() if count > 1)


def _sites(dna: str) -> int:
    reverse = _reverse_complement(dna)
    hits = 0
    for _name, site in _ASSEMBLY_SITES:
        hits += dna.count(site)
        if site != _reverse_complement(site):
            hits += reverse.count(site)
    return hits


def _motif_hits(dna: str, motif: str) -> int:
    reverse = _reverse_complement(dna)
    return dna.count(motif) + reverse.count(motif)


def _reverse_complement(dna: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return dna.translate(table)[::-1]
