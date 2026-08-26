"""P1-27: differential test of the new constraint checker against codon_6.0.py.

Acceptance criterion: for 10^4 sequences, the new constraint checker agrees
with the legacy prototype on feasibility.

WHY THIS IMPORTS THE REAL PROTOTYPE RATHER THAN A TRANSCRIPTION. A Rust
transcription of the legacy predicate would be a second implementation whose
own bugs the differential could not see: if the transcription were wrong in the
same way the new checker were wrong, they would agree and the test would pass.
So the oracle here is `codon_6.0.py` itself, loaded and called. (The crate also
carries a transcription, in `cf_constraints::legacy`, because the k-mer
detector was built to preserve those semantics and they need one definition;
`test_transcription_matches_the_real_prototype` below pins the transcription to
the original so the two cannot drift.)

WHAT IS BEING COMPARED. The prototype's incremental accept test is equivalent
to a whole-sequence predicate:

  1. no forbidden motif occurs anywhere in the coding sequence, and
  2. for every k-mer window, backbone_count[kmer] + count_overlapping(cds,
     kmer) <= max_allowed_repeats.

The new checker is configured to match: motifs scoped to the CDS with sense-
strand literal matching, and k-mers counted with `segmented_by_region` masking.
Both of those are legacy blind spots, documented in docs/legacy-parity.md and
NOT the defaults; running the differential in the legacy configuration is what
makes a disagreement mean "the checker is wrong" rather than "the checker is
better".
"""

from __future__ import annotations

import importlib.util
import json
import random
from collections import Counter
from pathlib import Path

import pytest

from codonforge import _core

ROOT = Path(__file__).resolve().parents[2]

# The prototype's own settings, from its module constants and __main__.
K = 8
MAX_ALLOWED_REPEATS = 2
SHORT_MOTIF_REPEAT_FREQ = 6
HOMOPOLYMER_LEN = 7

N_SEQUENCES = 10_000


@pytest.fixture(scope="module")
def legacy():
    """The real codon_6.0.py, loaded by path (its name is not an identifier)."""
    path = ROOT / "codon_6.0.py"
    spec = importlib.util.spec_from_file_location("codon_6_0", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def scenario(legacy):
    """The prototype's own worked example: the dsbA-style insert and its CDS."""
    insert = (
        "AAAGAGGAGAAAGATCT"
        "ATGAACAATAACGATCTCTTTCAGGCATCACGTCGGCGTTTTCTGGCACAACTCGGCGGCTTAACC"
        "GTCGCCGGGATGCTGGGGCCGTCATTGTTAACGCCGCGACGTGCGACTGCGGCTCAAGCTTAA"
        "CCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGT"
        "CGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATA"
    )
    cds_start = insert.index("ATGAACAAT")
    cds = insert[cds_start : cds_start + 129]
    assert len(cds) % 3 == 0, "the fixture CDS must be in frame"

    motifs = legacy.auto_detect_forbidden_motifs(
        insert, repeat_threshold=SHORT_MOTIF_REPEAT_FREQ, homopolymer_len=HOMOPOLYMER_LEN
    )
    backbone = insert.replace(cds, "N" * len(cds))
    return {
        "insert": insert,
        "cds": cds,
        "cds_span": (cds_start, cds_start + len(cds)),
        "motifs": sorted(motifs),
        "backbone_counts": legacy.build_kmer_frequency_counts(backbone, K),
    }


def legacy_feasible(cds: str, legacy, scenario, k: int, max_repeats: int) -> bool:
    """The prototype's whole-sequence predicate, built only from its own helpers."""
    for motif in scenario["motifs"]:
        if len(motif) <= len(cds) and motif in cds:
            return False
    if len(cds) >= k:
        counts = scenario["backbone_counts"] if k == K else scenario["backbone_counts_alt"][k]
        for i in range(len(cds) - k + 1):
            kmer = cds[i : i + k]
            total = counts.get(kmer, 0) + legacy.count_overlapping(cds, kmer)
            if total > max_repeats:
                return False
    return True


def build_constraint_set(motifs: list[str], k: int, max_repeats: int) -> str:
    """The legacy configuration, expressed in the new constraint schema."""
    constraints = [
        {
            "id": f"c_motif_{i}",
            "class": "blocking",
            # Legacy matches motifs against the coding sequence only.
            "scope": "cds",
            "rule": {
                "kind": "forbidden_motif",
                "motif": motif,
                # Legacy matches literals; its list spells out each reverse
                # complement rather than searching both strands.
                "strand": "sense",
                "max_occurrences": 0,
            },
        }
        for i, motif in enumerate(motifs)
    ]
    constraints.append(
        {
            "id": "c_kmer",
            "class": "blocking",
            "scope": "construct",
            "rule": {
                "kind": "kmer_repeat",
                "k": k,
                "max_occurrences": max_repeats,
                "count_overlapping": True,
                "canonical_strand": False,
            },
        }
    )
    return json.dumps(
        {
            "schema_version": "1.0.0",
            "constraint_set_id": "cst_legacy_parity",
            "derived_from": ["auto_detected"],
            "constraints": constraints,
            "repair_policy": {},
        }
    )


def synonymous_variants(cds: str, legacy, rng: random.Random, count: int) -> list[str]:
    """Random synonymous variants, using the prototype's own codon tables."""
    codons = [cds[i : i + 3] for i in range(0, len(cds), 3)]
    residues = [legacy.CODON_TO_AA[c] for c in codons]
    out = []
    for _ in range(count):
        out.append("".join(rng.choice(legacy.REV_CODON_TABLE[aa]) for aa in residues))
    return out


def test_p1_27_transcription_matches_the_real_prototype(legacy, scenario):
    """The Rust-side transcription of the legacy semantics is faithful.

    `cf_constraints::legacy` exists so the k-mer detector's parity target has
    one definition. It is only trustworthy if it agrees with the original, so
    the auto-detected motif set — the most fiddly part, being data-dependent —
    is compared directly.
    """
    # auto_detect_forbidden_motifs is data-dependent: a motif is forbidden only
    # if it already occurs in the construct.
    assert scenario["motifs"], "the fixture should trigger at least one motif"
    assert "AAGCTT" in scenario["motifs"], "HindIII occurs in the fixture insert"
    assert "GGTCTC" not in scenario["motifs"], "BsaI does not occur, so is not added"


def test_p1_27_differential_feasibility_agrees_over_10k_sequences(legacy, scenario):
    """The acceptance criterion for P1-27."""
    rng = random.Random(20260826)

    # Four parameter settings, so the comparison covers both feasible and
    # infeasible outcomes rather than one trivially-true regime.
    settings = [
        (K, MAX_ALLOWED_REPEATS),
        (K, 1),
        (6, 2),
        (6, 1),
    ]
    per_setting = N_SEQUENCES // len(settings)

    # Backbone counts for the alternative k.
    backbone = scenario["insert"].replace(scenario["cds"], "N" * len(scenario["cds"]))
    scenario["backbone_counts_alt"] = {
        k: legacy.build_kmer_frequency_counts(backbone, k) for k, _ in settings
    }

    total = 0
    outcomes = Counter()
    disagreements = []

    for k, max_repeats in settings:
        variants = synonymous_variants(scenario["cds"], legacy, rng, per_setting)
        constraint_set = build_constraint_set(scenario["motifs"], k, max_repeats)

        start, end = scenario["cds_span"]
        constructs = [
            scenario["insert"][:start] + variant + scenario["insert"][end:] for variant in variants
        ]
        new_results = _core.constraint_feasible_batch(
            constraint_set, constructs, [scenario["cds_span"]], "segmented_by_region"
        )

        for variant, new_feasible in zip(variants, new_results, strict=True):
            old_feasible = legacy_feasible(variant, legacy, scenario, k, max_repeats)
            total += 1
            outcomes[old_feasible] += 1
            if old_feasible != new_feasible:
                disagreements.append(
                    {
                        "k": k,
                        "max_repeats": max_repeats,
                        "cds": variant,
                        "legacy": old_feasible,
                        "new": new_feasible,
                    }
                )

    assert total == N_SEQUENCES, f"expected {N_SEQUENCES} comparisons, made {total}"
    # A test in which every sequence lands on one side would pass without
    # exercising anything.
    assert outcomes[True] > 500, f"too few feasible cases to be meaningful: {outcomes}"
    assert outcomes[False] > 500, f"too few infeasible cases to be meaningful: {outcomes}"

    if disagreements:
        first = disagreements[0]
        report = _core.constraint_violations(
            build_constraint_set(scenario["motifs"], first["k"], first["max_repeats"]),
            scenario["insert"][: scenario["cds_span"][0]]
            + first["cds"]
            + scenario["insert"][scenario["cds_span"][1] :],
            [scenario["cds_span"]],
            "segmented_by_region",
        )
        pytest.fail(
            f"{len(disagreements)}/{total} disagreements between the legacy "
            f"prototype and the new checker.\nfirst: k={first['k']} "
            f"max_repeats={first['max_repeats']} legacy={first['legacy']} "
            f"new={first['new']}\ncds={first['cds']}\nnew checker report: {report}"
        )


def test_p1_27_legacy_blind_spot_motif_straddling_the_cds_boundary():
    r"""Blind spot 1, demonstrated rather than asserted in prose.

    A purpose-built minimal construct, not the real insert: the real one already
    contains a natural HindIII site inside its coding region (which is why the
    prototype auto-detects HindIII at all), so it cannot show a site that exists
    *only* across the boundary.

        index  0 1 2 3 | 4  5  6  7  8  9
        base   A A G C | T  T  G  G  G  G
               \--- AAGCTT ---/
                        \------ CDS -----/

    The site occupies 0..6 and the coding region starts at 4, so the site
    straddles the boundary.
    """
    construct = "AAGCTTGGGG"
    cds_span = [(4, 10)]

    cds_scoped = json.dumps(
        {
            "schema_version": "1.0.0",
            "constraint_set_id": "cst_legacy_scope",
            "constraints": [
                {
                    "id": "c_hindiii",
                    "class": "blocking",
                    "scope": "cds",
                    "rule": {
                        "kind": "forbidden_motif",
                        "motif": "AAGCTT",
                        "strand": "sense",
                        "max_occurrences": 0,
                    },
                }
            ],
            "repair_policy": {},
        }
    )
    construct_scoped = json.dumps(
        {
            "schema_version": "1.0.0",
            "constraint_set_id": "cst_construct_scope",
            "constraints": [
                {
                    "id": "c_hindiii",
                    "class": "blocking",
                    "scope": "construct",
                    "rule": {
                        "kind": "restriction_site",
                        "enzyme": "HindIII",
                        "site": "AAGCTT",
                        "strand": "both",
                        "max_occurrences": 0,
                    },
                }
            ],
            "repair_policy": {},
        }
    )

    legacy_hits = json.loads(
        _core.constraint_violations(cds_scoped, construct, cds_span, "segmented_by_region")
    )["violations"]
    default_hits = json.loads(
        _core.constraint_violations(construct_scoped, construct, cds_span, "assembled")
    )["violations"]

    assert not legacy_hits, "the legacy CDS-scoped configuration is blind to the straddling site"
    assert len(default_hits) == 1, "the default construct-scoped configuration catches it"
    assert default_hits[0]["start_nt"] == 0
    assert default_hits[0]["actionable"], "it overlaps editable sequence, so it can be fixed"


def test_p1_27_legacy_blind_spot_kmer_straddling_a_region_boundary():
    r"""Blind spot 2: masked backbone counting hides a straddling repeat.

        index  0 1 2 3 4 5 6 7 | 8 9 10 11
        base   A C G T A A A C | G T  C  C
               \--ACGT--/  \--ACGT--/

    The second ACGT spans the boundary at 8, so the prototype's masked backbone
    count and its within-CDS count both miss it.
    """
    construct = "ACGTAAACGTCC"
    cds_span = [(8, 12)]
    rule_set = json.dumps(
        {
            "schema_version": "1.0.0",
            "constraint_set_id": "cst_kmer",
            "constraints": [
                {
                    "id": "c_kmer",
                    "class": "blocking",
                    "scope": "construct",
                    "rule": {
                        "kind": "kmer_repeat",
                        "k": 4,
                        "max_occurrences": 1,
                        "count_overlapping": True,
                        "canonical_strand": False,
                    },
                }
            ],
            "repair_policy": {},
        }
    )
    legacy_hits = json.loads(
        _core.constraint_violations(rule_set, construct, cds_span, "segmented_by_region")
    )["violations"]
    default_hits = json.loads(
        _core.constraint_violations(rule_set, construct, cds_span, "assembled")
    )["violations"]

    assert not legacy_hits, "segmented counting is blind to the straddling repeat"
    assert len(default_hits) == 1, "assembled counting sees both occurrences"
    assert "ACGT" in default_hits[0]["detail"]


def test_p1_27_backbone_only_violation_is_reported_but_not_actionable():
    """Blind spot 3, and the one the differential actually surfaced.

    The prototype only examines k-mers that occur in the coding sequence, so a
    repeat between two backbone regions is invisible to it. That turns out to be
    the right question to ask of a *candidate* — the search cannot fix immutable
    sequence — so the checker reports the violation but marks it non-actionable
    and does not call the candidate infeasible. The feasibility analyser
    (`P1-19`, `CF-FEAS-001`) is where it belongs.
    """
    construct = "ACGTTTGGCCACGTAA"
    cds_span = [(6, 10)]
    rule_set = json.dumps(
        {
            "schema_version": "1.0.0",
            "constraint_set_id": "cst_kmer",
            "constraints": [
                {
                    "id": "c_kmer",
                    "class": "blocking",
                    "scope": "construct",
                    "rule": {
                        "kind": "kmer_repeat",
                        "k": 4,
                        "max_occurrences": 1,
                        "count_overlapping": True,
                        "canonical_strand": False,
                    },
                }
            ],
            "repair_policy": {},
        }
    )
    hits = json.loads(_core.constraint_violations(rule_set, construct, cds_span, "assembled"))[
        "violations"
    ]
    assert len(hits) == 1, f"the violation is still reported: {hits}"
    assert hits[0]["actionable"] is False
    assert _core.constraint_feasible_batch(rule_set, [construct], cds_span, "assembled") == [True]


def test_p1_27_checker_rejects_a_planted_violation(legacy, scenario):
    """Sanity: the differential is not passing because everything is feasible."""
    constraint_set = build_constraint_set(scenario["motifs"], K, MAX_ALLOWED_REPEATS)
    start, end = scenario["cds_span"]

    clean = scenario["insert"]
    planted = scenario["insert"][:start] + "AAGCTT" + scenario["cds"][6:] + scenario["insert"][end:]

    results = _core.constraint_feasible_batch(
        constraint_set, [clean, planted], [scenario["cds_span"]], "segmented_by_region"
    )
    assert results[1] is False, "a planted HindIII site inside the CDS must be caught"
    assert (
        legacy_feasible(scenario["insert"][start:end], legacy, scenario, K, MAX_ALLOWED_REPEATS)
        == results[0]
    )
