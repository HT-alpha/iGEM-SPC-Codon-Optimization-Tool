"""P1-27b — the Rust `baseline_backtracking` port against the real prototype.

Ticket P1-27b.

`P1-27a` established the standard: an oracle must be the actual `codon_6.0.py`,
loaded and called, not a transcription that could share a bug with the thing it
checks. The Rust port in `cf-search` is a *reference solver*, not an oracle, and
it earns that role only by being demonstrably faithful — so this compares it
against the prototype before anything else compares against it.

WHAT IS COMPARED, AND WHY NOT SEQUENCE EQUALITY. The prototype draws its codon
ordering from Python's Mersenne Twister via `random.uniform`. Reproducing that
stream bit-for-bit would mean vendoring MT19937 to imitate a script that never
promised reproducibility, and the resulting equality would be a fact about two
RNGs rather than about two searches. So the comparison is on the properties that
carry meaning:

  1. the acceptance predicate — the prototype's incremental prefix check and the
     port's must agree on every prefix, which is where a transcription error
     would actually live;
  2. solvability — on instances the prototype solves, the port solves too, and
     on instances neither can solve, neither claims to;
  3. validity — every solution either produces satisfies the *other* side's
     predicate.

Item 1 is the strong one. It compares the two implementations position by
position over thousands of prefixes, which is a far tighter net than comparing
two final answers.
"""

from __future__ import annotations

import importlib.util
import random
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]

# The prototype's own settings.
K = 8
MAX_ALLOWED_REPEATS = 2
SHORT_MOTIF_REPEAT_FREQ = 6
HOMOPOLYMER_LEN = 7

N_PREFIXES = 5_000


@pytest.fixture(scope="module")
def legacy():
    """The real codon_6.0.py, loaded by path (its name is not an identifier)."""
    path = ROOT / "codon_6.0.py"
    spec = importlib.util.spec_from_file_location("codon_6_0_b", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def scenario(legacy):
    insert = (
        "AAAGAGGAGAAAGATCT"
        "ATGAACAATAACGATCTCTTTCAGGCATCACGTCGGCGTTTTCTGGCACAACTCGGCGGCTTAACC"
        "GTCGCCGGGATGCTGGGGCCGTCATTGTTAACGCCGCGACGTGCGACTGCGGCTCAAGCTTAA"
        "CCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGT"
        "CGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATA"
    )
    cds_start = insert.index("ATGAACAAT")
    cds = insert[cds_start : cds_start + 129]
    motifs = legacy.auto_detect_forbidden_motifs(
        insert, repeat_threshold=SHORT_MOTIF_REPEAT_FREQ, homopolymer_len=HOMOPOLYMER_LEN
    )
    backbone = insert.replace(cds, "N" * len(cds))
    counts = legacy.build_kmer_frequency_counts(backbone, K)
    return {
        "insert": insert,
        "cds": cds,
        "motifs": sorted(motifs),
        # The prototype indexes this with `[kmer]`, relying on it being a
        # defaultdict. Handing it a plain dict raises KeyError on the first
        # unseen k-mer, so the two forms are kept apart deliberately: the
        # prototype gets its own object, the port gets a plain mapping.
        "backbone_counts": counts,
        "backbone_plain": dict(counts),
    }


def prototype_prefix_accepts(legacy, prefix: str, motifs, k: int, limit: int, counts) -> bool:
    """The prototype's own incremental check, extracted verbatim.

    This is the block inside `optimize_cds_sequence`'s `while stack:` loop that
    decides whether `test_seq` may be extended. It is reproduced here rather
    than called because the prototype interleaves it with stack bookkeeping and
    printing; every expression, including the two `+ 2` window slacks, is copied
    unchanged from lines 167-184 of `codon_6.0.py`.
    """
    test_seq = prefix
    for motif in motifs:
        motif_len = len(motif)
        if len(test_seq) >= motif_len:
            check_tail = test_seq[max(0, len(test_seq) - motif_len - 2) :]
            if motif in check_tail:
                return False

    if len(test_seq) >= k:
        check_start = max(0, len(test_seq) - (k + 2))
        for j in range(check_start, len(test_seq) - k + 1):
            kmer = test_seq[j : j + k]
            existing = counts.get(kmer, 0) + legacy.count_overlapping(test_seq, kmer)
            if existing > limit:
                return False
    return True


# (k, max_allowed_repeats). Two regimes on purpose: at k=8 with a limit of 2
# almost every synonymous prefix is accepted, so a comparison run only there
# would be agreement about the accepting branch alone. A short k with a limit of
# 1 rejects most prefixes and exercises the other side.
REGIMES = ((8, 2), (4, 1))


def test_p1_27b_the_prefix_predicate_matches_the_prototype(legacy, scenario):
    """The port's incremental check agrees with the prototype on every prefix.

    This is where a transcription error would live: the two windowed checks and
    their `+ 2` slack are easy to paraphrase into something that agrees on most
    inputs. Comparing thousands of prefixes across two regimes catches a
    disagreement at the exact position it first appears, which comparing final
    answers would not.
    """
    from codonforge import _core

    rng = random.Random(20260826)
    cds = scenario["cds"]
    protein = "".join(legacy.CODON_TO_AA[cds[i : i + 3]] for i in range(0, len(cds), 3))
    backbone = scenario["insert"].replace(cds, "N" * len(cds))

    checked = 0
    accepted = 0
    for k, limit in REGIMES:
        counts = legacy.build_kmer_frequency_counts(backbone, k)
        plain = dict(counts)
        for _ in range(N_PREFIXES):
            length = rng.randrange(1, len(protein) + 1)
            prefix = "".join(
                rng.choice(legacy.REV_CODON_TABLE[protein[i]]) for i in range(length)
            )
            mine = _core.legacy_prefix_accepts(prefix, scenario["motifs"], k, limit, plain)
            theirs = prototype_prefix_accepts(
                legacy, prefix, scenario["motifs"], k, limit, counts
            )
            assert mine == theirs, (
                f"prefix predicate disagreed at k={k} limit={limit} length={length}: "
                f"port={mine} prototype={theirs} for {prefix}"
            )
            accepted += int(theirs)
            checked += 1

    assert checked == N_PREFIXES * len(REGIMES)
    # Both outcomes must be well represented, or the agreement is agreement
    # about nothing — the same requirement TST-CONS-007 imposes.
    assert 0.1 < accepted / checked < 0.9, (
        f"{accepted}/{checked} prefixes were accepted; the fixture does not exercise "
        "both branches"
    )


def test_p1_27b_the_port_solves_what_the_prototype_solves(legacy, scenario):
    """Solvability agrees on the prototype's own worked instance.

    Not sequence equality: the two draw their codon ordering from different
    RNGs, so the answers differ and should. What must agree is whether an answer
    exists, and whether the one produced is valid under the other side's rules.
    """
    from codonforge import _core

    random.seed(20260826)
    theirs = legacy.optimize_cds_sequence(
        scenario["cds"],
        scenario["backbone_counts"],
        K,
        MAX_ALLOWED_REPEATS,
        scenario["motifs"],
    )
    assert theirs is not None, "the prototype cannot solve its own example"

    solved = 0
    for seed in range(20):
        mine = _core.baseline_backtracking(
            scenario["cds"],
            1,
            scenario["motifs"],
            K,
            MAX_ALLOWED_REPEATS,
            scenario["backbone_plain"],
            seed,
        )
        if mine is None:
            continue
        solved += 1
        # Same protein.
        assert _core.translate(mine, 1) == _core.translate(scenario["cds"], 1)
        # Valid under the prototype's own predicate, evaluated by the prototype.
        assert prototype_prefix_accepts(
            legacy, mine, scenario["motifs"], K, MAX_ALLOWED_REPEATS, scenario["backbone_counts"]
        ), f"seed {seed}: the port returned a sequence the prototype rejects"

    assert solved >= 18, f"the port solved only {solved}/20 of an instance the prototype solves"

    # And the reverse direction: the prototype's answer passes the port's check.
    assert _core.legacy_prefix_accepts(
        theirs, scenario["motifs"], K, MAX_ALLOWED_REPEATS, scenario["backbone_plain"]
    )


def test_p1_27b_neither_claims_to_solve_the_unsolvable(legacy, scenario):
    """An impossible instance is refused by both, and neither reports a proof.

    The prototype prints "Exceeded max restarts" and returns `None`; the port
    returns no solution with `BudgetExhausted`. Both are timeouts, not proofs —
    the same distinction `P1-17`'s `GiveUpReason` draws — so neither may be read
    as evidence that no solution exists.
    """
    from codonforge import _core

    # Every codon of the first amino acid forbidden outright.
    cds = scenario["cds"]
    first_aa = legacy.CODON_TO_AA[cds[0:3]]
    impossible = sorted(legacy.REV_CODON_TABLE[first_aa])

    random.seed(1)
    theirs = legacy.optimize_cds_sequence(
        cds, scenario["backbone_counts"], K, MAX_ALLOWED_REPEATS, impossible
    )
    assert theirs is None

    mine = _core.baseline_backtracking(
        cds, 1, impossible, K, MAX_ALLOWED_REPEATS, scenario["backbone_plain"], 5
    )
    assert mine is None
