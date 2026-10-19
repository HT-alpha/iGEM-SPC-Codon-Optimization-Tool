"""Gibson context for a codon design.

The codon search already counts k-mers on the assembled construct, which is
the 5' flank, the CDS, and the 3' flank. ``examples/gibson.json`` is the
blocking rule that asks for that count. This module does not keep a second
copy of the rule and does not count k-mers itself.

An homology arm is the junction-adjacent slice of a flank the caller
supplied. Bases are not invented. Overlap melting temperature is not
computed.

Nothing here is specific to one gene or one organism.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from dna_design_suite import _core

# Gibson et al., Nature Methods 2009, 6:343–345, is the method. A protocol
# that reprints the ISO-buffer mix describes 20–150 bp overlaps as suitable
# for that mix (Salis lab protocol book, gibson-assembly.md). This module
# did not measure an optimum. 40 sits inside that range. Lengths outside it
# are refused rather than labelled as Gibson arms.
ARM_LENGTH_DEFAULT = 40
ARM_LENGTH_MIN = 20
ARM_LENGTH_MAX = 150

_DETAIL = re.compile(
    r"^(?P<k>\d+)-mer (?P<kmer>[ACGT]+) occurs (?P<count>\d+) times "
    r"\(limit (?P<limit>\d+)\), at \[(?P<positions>[0-9, ]*)\]$"
)

TITLES = {
    "CF-GIB-001": "A homology arm cannot be copied from the supplied flank",
    "CF-GIB-002": "No flanking sequence was supplied on this side",
    "CF-GIB-003": "Gibson assembly input is malformed",
}


class GibsonError(Exception):
    """A refusal with a registered ``CF-GIB-*`` code."""

    def __init__(self, code: str, detail: str) -> None:
        self.code = code
        self.detail = detail
        super().__init__(f"{code}: {detail}")


def constraint_path() -> Path:
    """The one constraint file. The packaged copy is the wheel's fallback."""
    here = Path(__file__).resolve()
    candidates = (
        here.parents[3] / "examples" / "gibson.json",
        here.parents[1] / "api" / "examples" / "gibson.json",
    )
    for path in candidates:
        if path.is_file():
            return path
    raise GibsonError(
        "CF-GIB-003",
        "examples/gibson.json was not found in the checkout or the installed package",
    )


def constraint_set() -> dict[str, Any]:
    """The Gibson constraint document, parsed from that file."""
    try:
        document = json.loads(constraint_path().read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise GibsonError("CF-GIB-003", f"examples/gibson.json is not JSON: {exc}") from exc
    rules = [item for item in document.get("constraints", []) if _is_kmer(item)]
    if len(rules) != 1:
        raise GibsonError(
            "CF-GIB-003",
            "examples/gibson.json must contain exactly one kmer_repeat rule",
        )
    return document


def assess(
    cds: str,
    five_prime_flank: str = "",
    three_prime_flank: str = "",
    arm_length: int = ARM_LENGTH_DEFAULT,
) -> dict[str, Any]:
    """Homology arms and direct repeats for one CDS and its flanks.

    The repeat list is the constraint checker's list for
    ``examples/gibson.json`` on ``five_prime_flank + cds + three_prime_flank``.
    """
    coding = _dna(cds, "cds")
    five = _dna(five_prime_flank, "five_prime_flank")
    three = _dna(three_prime_flank, "three_prime_flank")
    if not coding:
        raise GibsonError("CF-GIB-003", "cds is empty")
    length = _arm_length(arm_length)
    document = constraint_set()
    rule = next(item for item in document["constraints"] if _is_kmer(item))
    spec = rule["rule"]
    k = int(spec["k"])
    five_arm = _arm(five, length, "5'")
    three_arm = _arm(three, length, "3'")
    sites = _sites(document, coding, five, three, five_arm, three_arm)
    problems = _problems(five_arm, three_arm)
    refused = any(item["code"] == "CF-GIB-001" for item in problems)
    if refused:
        fragment = None
        fragment_note = (
            "No insert sequence is reported. A flank you supplied is shorter than "
            "the requested arm, and the missing bases were not invented."
        )
    elif five_arm["designed"] or three_arm["designed"]:
        fragment = (five_arm["sequence"] or "") + coding + (three_arm["sequence"] or "")
        fragment_note = (
            "The insert is each designed arm copied onto the CDS. It is the DNA "
            "for those junctions. The repeat scan used the full flanks, not only "
            "these slices."
        )
    else:
        fragment = None
        fragment_note = "Neither flank was supplied, so no insert was assembled."
    n_actionable = sum(1 for site in sites if site["actionable"])
    n_flank_only = sum(1 for site in sites if site["class"] == "flank_only")
    return {
        "schema_version": "1.0.0",
        "kind": "gibson_assembly",
        "guarantees_assembly": False,
        "guarantee_statement": (
            f"A report with no actionable {k}-mer does not mean the fragments will "
            "assemble, and it does not mean the construct is free of recombination. "
            f"The scan is exact identity of {k} nt on the linear sequence that was supplied."
        ),
        "limitations": _limitations(k),
        "repeat_rule": {
            "source": "examples/gibson.json",
            "constraint_id": rule["id"],
            "scope": rule.get("scope", "construct"),
            "k": k,
            "max_occurrences": spec["max_occurrences"],
            "count_overlapping": spec.get("count_overlapping", True),
            "canonical_strand": spec.get("canonical_strand", False),
        },
        "context": {
            "five_prime_flank_length": len(five),
            "three_prime_flank_length": len(three),
            "cds_length": len(coding),
            "cds_span": [len(five), len(five) + len(coding)],
            "scanned": "five_prime_flank + cds + three_prime_flank",
        },
        "homology_arms": {
            "arm_length_requested": length,
            "arm_length_default": ARM_LENGTH_DEFAULT,
            "arm_length_min": ARM_LENGTH_MIN,
            "arm_length_max": ARM_LENGTH_MAX,
            "method": (
                "Copy the junction-adjacent slice of the supplied flank. "
                "No bases are added. Melting temperature is not computed."
            ),
            "five_prime": five_arm,
            "three_prime": three_arm,
            "fragment": fragment,
            "fragment_note": fragment_note,
            "tm_celsius": None,
        },
        "recombination": {
            "k": k,
            "sites": sites,
            "n_sites": len(sites),
            "n_actionable": n_actionable,
            "n_flank_only": n_flank_only,
            "n_cds_context": sum(1 for site in sites if site["class"] == "cds_context"),
            "n_cds_internal": sum(1 for site in sites if site["class"] == "cds_internal"),
        },
        "problems": problems,
    }


def format_text(report: dict[str, Any]) -> str:
    """The human view of :func:`assess`. It prints fields, it does not recompute them."""
    arms = report["homology_arms"]
    repeats = report["recombination"]
    lines = [
        "gibson",
        f"  {report['guarantee_statement']}",
    ]
    for label, key in (("5'", "five_prime"), ("3'", "three_prime")):
        arm = arms[key]
        if arm["designed"]:
            lines.append(f"  {label} arm  {arm['sequence']}")
        else:
            lines.append(f"  {label} arm  not designed  {arm['code']}  {arm['detail']}")
    if arms["fragment"]:
        lines.append(f"  fragment  {arms['fragment']}")
    else:
        lines.append("  fragment  not reported")
    lines.append(
        f"  repeats  {repeats['n_sites']} {repeats['k']}-mer(s) over the limit "
        f"({repeats['n_actionable']} the codon search can try to edit, "
        f"{repeats['n_flank_only']} only in the flanks)"
    )
    for site in repeats["sites"]:
        lines.append(
            f"    {site['class']}  {site['kmer']}  x{site['count']}  "
            f"actionable={str(site['actionable']).lower()}"
        )
    return "\n".join(lines) + "\n"


def _limitations(k: int) -> list[str]:
    return [
        (
            "The scan is the 5' flank, then the CDS, then the 3' flank, as one "
            "linear string. Sequence that was not supplied is invisible, including "
            "the rest of a vector."
        ),
        (
            "The outer ends are not joined. A repeat that would exist only because "
            "the plasmid is circular is not counted."
        ),
        (
            "Only exact direct repeats are counted. A mismatch, a gap, or the "
            "reverse complement is not a hit. Inverted repeats are not treated as "
            "recombination sites."
        ),
        (
            f"The repeat length is {k} nt, the same length as the synthesis preset, "
            "not a measured Gibson or RecA threshold. Gibson et al., Nature Methods "
            "2009, 6:343–345, is the assembly method. A protocol that reprints the "
            "ISO-buffer mix describes 20–150 bp overlaps as suitable for that mix. "
            f"Clearing exact {k}-mers can reject a design that would have assembled, "
            "and it misses a longer imperfect repeat."
        ),
        (
            "A repeat that lies entirely in a flank cannot be removed by codon "
            "optimization. It is reported and left in place."
        ),
        (
            "A homology arm is a copy of the supplied flank. Overlap melting "
            "temperature, exonuclease chew-back, polymerase fill-in, and fragment "
            "concentration are not evaluated."
        ),
        "The same rule is applied to every gene and every organism.",
    ]


def _problems(five_arm: dict[str, Any], three_arm: dict[str, Any]) -> list[dict[str, Any]]:
    problems = []
    for subject, arm in (("five_prime", five_arm), ("three_prime", three_arm)):
        code = arm.get("code")
        if code:
            status = 200 if code == "CF-GIB-002" else 422
            problems.append(
                {
                    "type": f"https://dna_design_suite.org/problems/{code.lower()}",
                    "title": TITLES[code],
                    "status": status,
                    "code": code,
                    "detail": arm["detail"],
                    "stage": "gibson",
                    "subject": subject,
                }
            )
    return problems


def _sites(
    document: dict[str, Any],
    cds: str,
    five: str,
    three: str,
    five_arm: dict[str, Any],
    three_arm: dict[str, Any],
) -> list[dict[str, Any]]:
    cds_start = len(five)
    cds_end = cds_start + len(cds)
    raw = json.loads(
        _core.constraint_violations(
            json.dumps(document),
            five + cds + three,
            [(cds_start, cds_end)],
            "assembled",
        )
    )
    sites = []
    for item in raw["violations"]:
        parsed = _DETAIL.fullmatch(item["detail"])
        if parsed is None:
            raise GibsonError(
                "CF-GIB-003",
                "the repeat checker returned a detail this report cannot read: " + item["detail"],
            )
        k = int(parsed.group("k"))
        positions = [int(part) for part in parsed.group("positions").split(",") if part.strip()]
        occurrences = [
            {"start_nt": pos, "end_nt": pos + k, "where": _where(pos, k, cds_start, cds_end)}
            for pos in positions
        ]
        places = {occ["where"] for occ in occurrences}
        touches = any(place in {"cds", "junction"} for place in places)
        outside = any(place in {"five_prime_flank", "three_prime_flank"} for place in places)
        if not touches:
            klass = "flank_only"
        elif outside:
            klass = "cds_context"
        else:
            klass = "cds_internal"
        kmer = parsed.group("kmer")
        five_sequence = five_arm["sequence"] or ""
        three_sequence = three_arm["sequence"] or ""
        sites.append(
            {
                "constraint_id": item["constraint_id"],
                "kmer": kmer,
                "k": k,
                "count": int(parsed.group("count")),
                "limit": int(parsed.group("limit")),
                "actionable": bool(item["actionable"]),
                "class": klass,
                "in_five_prime_arm": bool(five_sequence) and kmer in five_sequence,
                "in_three_prime_arm": bool(three_sequence) and kmer in three_sequence,
                "occurrences": occurrences,
                "detail": item["detail"],
            }
        )
    sites.sort(key=lambda site: (site["occurrences"][0]["start_nt"], site["kmer"]))
    return sites


def _where(pos: int, k: int, cds_start: int, cds_end: int) -> str:
    end = pos + k
    if pos < cds_end and end > cds_start:
        if pos >= cds_start and end <= cds_end:
            return "cds"
        return "junction"
    if end <= cds_start:
        return "five_prime_flank"
    return "three_prime_flank"


def _arm(flank: str, arm_length: int, side: str) -> dict[str, Any]:
    if side == "5'":
        source = "suffix of five_prime_flank"
        sequence = flank[-arm_length:] if flank else ""
        subject = "5' flank"
    else:
        source = "prefix of three_prime_flank"
        sequence = flank[:arm_length] if flank else ""
        subject = "3' flank"
    if not flank:
        return {
            "designed": False,
            "sequence": None,
            "length": 0,
            "source": source,
            "code": "CF-GIB-002",
            "detail": (
                f"no {subject} was supplied, so no homology arm was copied and the "
                "repeat scan did not see a neighbor on that side"
            ),
        }
    if len(flank) < arm_length:
        return {
            "designed": False,
            "sequence": None,
            "length": len(flank),
            "source": source,
            "code": "CF-GIB-001",
            "detail": (
                f"the {subject} is {len(flank)} nt and the requested homology arm is "
                f"{arm_length} nt. The arm is copied from the flank. Supply a longer "
                "flank or a shorter arm_length. No bases were added."
            ),
        }
    return {
        "designed": True,
        "sequence": sequence,
        "length": arm_length,
        "source": source,
        "code": None,
        "detail": f"{arm_length} nt copied from the {source}",
    }


def _arm_length(arm_length: int) -> int:
    if isinstance(arm_length, bool) or not isinstance(arm_length, int):
        raise GibsonError("CF-GIB-003", "arm_length must be an integer")
    if not ARM_LENGTH_MIN <= arm_length <= ARM_LENGTH_MAX:
        raise GibsonError(
            "CF-GIB-003",
            f"arm_length must be from {ARM_LENGTH_MIN} to {ARM_LENGTH_MAX} nt. "
            f"{arm_length} is outside that range. The range is the 20–150 nt overlap "
            "a protocol reprinting the ISO-buffer Gibson mix describes as suitable. "
            "This tool did not measure an optimum.",
        )
    return arm_length


def _dna(value: str, label: str) -> str:
    if not isinstance(value, str):
        raise GibsonError("CF-GIB-003", f"{label} must be a string")
    text = "".join(value.upper().split())
    illegal = [base for base in text if base not in "ACGT"]
    if illegal:
        raise GibsonError(
            "CF-GIB-003",
            f"{label} contains {illegal[0]!r}, which is not A, C, G or T. U is not rewritten to T.",
        )
    return text


def _is_kmer(item: dict[str, Any]) -> bool:
    rule = item.get("rule")
    return isinstance(rule, dict) and rule.get("kind") == "kmer_repeat"
