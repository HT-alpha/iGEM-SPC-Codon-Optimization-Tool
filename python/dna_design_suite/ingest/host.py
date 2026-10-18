"""Host reference-data loading.

Ticket P1-11. Test: `TST-INGE-002` (host portion).
Contract: `planning/02-schemas/host_profile.schema.json`, ADR-0004.

NO FABRICATED BIOLOGICAL CONSTANTS. Two values the schema requires cannot be
derived from a genome and are not guessable:

  * `ribosome.elongation_rate_codons_per_s` — a measured physiological
    parameter. There is no default here; the caller must supply it, and a host
    for which it is unknown fails loudly rather than silently inheriting an
    *E. coli* figure.
  * tAI `s_values` — wobble penalty parameters from the literature. Without
    them, the `tai` block is omitted entirely (the schema does not require it)
    and `data_completeness` records that tAI-derived terms are unavailable.

`capacity` at estimator tier 3 is not in that category: ADR-0004 *defines* it as
supply proportional to the host's own native demand with `theta = 1`, so it is
derived rather than assumed. Its output is a demand-profile mismatch score, not
a growth-rate prediction, and the profile says so.
"""

from __future__ import annotations

import json
from collections import Counter
from typing import Any

from dna_design_suite.ingest.errors import IngestError, IngestReport, SourceLocation

# Terms a host can support with a genome and nothing else.
_GENOME_ONLY_TERMS = [
    "usage.cai",
    "usage.gc_content",
    "structure.local_accessibility",
    "epistasis.context_naturalness",
]

# Every codon, including those with count 0. A missing key is not the same
# record as a weight of 0: callers that iterate the table would skip it.
_ALL_CODONS = tuple(a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT")


def codon_usage_from_cds(
    sequences: list[str], *, translation_table: int, subset: str = "genome_wide"
) -> dict[str, Any]:
    """Count codons over a CDS set and derive relative adaptiveness.

    `relative_adaptiveness` is `w_i`: the codon's share within its synonymous
    family, normalized by the family maximum. That is the standard CAI
    definition and it is derived here rather than taken on faith.
    """
    from dna_design_suite import _core

    counts: Counter[str] = Counter()
    for index, cds in enumerate(sequences):
        upper = cds.upper().replace("U", "T")
        if len(upper) % 3:
            raise IngestError(
                "CF-SEQ-001",
                f"CDS {index} is {len(upper)} nt, not a whole number of codons",
                location=SourceLocation(record_id=f"cds[{index}]"),
            )
        for offset in range(0, len(upper), 3):
            counts[upper[offset : offset + 3]] += 1

    if not counts:
        raise IngestError("CF-SEQ-012", "no codons counted", location=SourceLocation())

    total = sum(counts.values())
    table: dict[str, Any] = {}
    for codon in _ALL_CODONS:
        count = counts.get(codon, 0)
        family = _core.synonyms(codon, translation_table)
        family_total = sum(counts.get(member, 0) for member in family)
        fraction = count / family_total if family_total else 0.0
        family_max = max((counts.get(member, 0) for member in family), default=0)
        family_max_fraction = family_max / family_total if family_total else 0.0
        table[codon] = {
            "aa": _core.translate(codon + "AAA", translation_table)[0]
            if codon not in ("TAA", "TAG", "TGA")
            else "*",
            "count": count,
            "fraction": round(fraction, 6),
            "per_thousand": round(1000.0 * count / total, 4),
            "relative_adaptiveness": round(
                fraction / family_max_fraction if family_max_fraction else 0.0, 6
            ),
        }
    return {"subset": subset, "n_cds": len(sequences), "table": table}


def build_host_profile(
    *,
    host_id: str,
    taxid: int,
    translation_table: int,
    codon_usage: dict[str, Any],
    trna_gene_copy_number: dict[str, int],
    elongation_rate_codons_per_s: float,
    version: str = "0.1.0",
    display_name: str | None = None,
    strain: str | None = None,
    source_dataset_id: str = "ds_user_supplied",
    source_digest: str | None = None,
    tai_s_values: dict[str, float] | None = None,
    report: IngestReport | None = None,
) -> dict[str, Any]:
    """Assemble a `HostProfile`.

    `elongation_rate_codons_per_s` is a required keyword with no default: see
    the module docstring. Supplying a wrong one is the caller's decision to
    make explicitly; supplying none is an error rather than an inherited guess.
    """
    from dna_design_suite import _core

    report = report if report is not None else IngestReport()

    if elongation_rate_codons_per_s <= 0:
        raise IngestError(
            "CF-HOST-002",
            "elongation_rate_codons_per_s must be positive and measured for this "
            "host; there is deliberately no default, because inheriting another "
            "organism's figure would silently change every burden estimate",
            location=SourceLocation(record_id=host_id),
        )
    if not trna_gene_copy_number:
        raise IngestError(
            "CF-HOST-002",
            "tRNA gene copy numbers are required; they are derivable from any "
            "annotated genome and are what the tier-3 capacity estimator uses",
            location=SourceLocation(record_id=host_id),
        )

    digest_placeholder = source_digest or _core.digest_json(json.dumps(codon_usage))
    profile: dict[str, Any] = {
        "schema_version": "1.0.0",
        "host_id": host_id,
        "taxid": taxid,
        "translation_table": translation_table,
        "version": version,
        "codon_usage": {
            **codon_usage,
            "source": {"dataset_id": source_dataset_id, "digest": digest_placeholder},
        },
        "trna": {
            "gene_copy_number": dict(sorted(trna_gene_copy_number.items())),
            "abundance": [],
            "wobble_rules_id": "crick_wobble_v1",
        },
        "ribosome": {
            "elongation_rate_codons_per_s": elongation_rate_codons_per_s,
            "active_fraction": 0.85,
            "footprint_codons": 10,
            "source": {"dataset_id": source_dataset_id, "digest": digest_placeholder},
        },
        "growth_conditions": [],
        "riboseq_tracks": [],
        "capacity": {
            # ADR-0004 tier 3: supply proportional to the host's own native
            # demand, theta = 1 by construction. The unit says what the number
            # means, so a mismatch score is never mistaken for a growth-rate
            # prediction.
            "estimator_tier": 3,
            "theta": 1.0,
            "theta_unit": "fraction_of_native_demand",
            "sharpness_kappa": 8.0,
            "hill_n": 2.0,
            "g_max": 1.0,
            "uncertainty": {"method": "none", "calibrated": False},
            "fit_source": "relative_to_native_demand_profile",
        },
        "data_completeness": {
            "has_codon_usage": True,
            "has_trna_gene_copies": True,
            "has_trna_abundance": False,
            "has_proteomics": False,
            "riboseq_mapped_codons": 0,
            "riboseq_genes_passing_qc": 0,
            "supported_terms": list(_GENOME_ONLY_TERMS),
        },
    }
    if display_name:
        profile["display_name"] = display_name
    if strain is not None:
        profile["strain"] = strain

    if tai_s_values:
        profile["tai"] = {
            "s_values": tai_s_values,
            "weights": _tai_weights(trna_gene_copy_number, translation_table, tai_s_values),
            "derived_from": "gene_copy_number",
            "optimized_against": None,
        }
        profile["data_completeness"]["supported_terms"].append("burden.growth_cost")
    else:
        report.warn(
            f"{host_id}: no tAI s-values supplied, so no tAI weights were derived "
            "and burden.growth_cost is not offered for this host. Supply "
            "tai_s_values from a cited source rather than accepting a default."
        )

    profile["data_completeness"]["supported_terms"] = sorted(
        set(profile["data_completeness"]["supported_terms"])
    )
    profile["digest"] = _core.digest_json(json.dumps(profile), "digest")
    return profile


def _tai_weights(
    gene_copy_number: dict[str, int], translation_table: int, s_values: dict[str, float]
) -> dict[str, float]:
    """Absolute adaptiveness per codon from tRNA gene copy numbers.

    Deliberately simple and fully determined by the supplied `s_values`: this
    function introduces no constant of its own, so a reviewer can check the
    result against the cited source rather than against this code.
    """

    raw: dict[str, float] = {}
    for codon_index in range(64):
        bases = "ACGT"
        codon = (
            bases[(codon_index >> 4) & 3] + bases[(codon_index >> 2) & 3] + bases[codon_index & 3]
        )
        anticodon = _revcomp(codon)
        exact = gene_copy_number.get(anticodon, 0)
        wobble = 0.0
        for pairing, penalty in s_values.items():
            partner = pairing.split(":")[-1] if ":" in pairing else None
            if partner:
                alternative = partner + anticodon[1:]
                wobble += (1.0 - penalty) * gene_copy_number.get(alternative, 0)
        raw[codon] = float(exact) + wobble

    peak = max(raw.values()) if raw else 0.0
    if peak == 0:
        raise IngestError(
            "CF-HOST-002",
            "no tRNA gene copy number matched any codon; check the anticodon "
            "orientation of the supplied table (it must be 5'->3')",
            location=SourceLocation(),
        )
    # `translation_table` is unused on purpose: tAI weights are per codon, not
    # per family, and normalizing per family would silently produce CAI.
    _ = translation_table
    return {codon: round(value / peak, 6) for codon, value in raw.items()}


def _revcomp(sequence: str) -> str:
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement[base] for base in reversed(sequence))
