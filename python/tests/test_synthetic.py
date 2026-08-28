"""TST-SYN-001 and TST-SYN-004: the synthetic suite's ground truths.

Ticket P1-32.

A synthetic suite is only worth having if its "ground truth" is genuinely known
rather than observed from the thing under test. So every assertion here either
checks the generator's construction against a closed-form prediction, or checks
a search result against exhaustive enumeration.

The tie coverage is the part worth reading. A suite of instances with one
obvious right answer cannot test tie-breaking anywhere, and tie-breaking is
where determinism bugs live: the tier-2 determinism fixture used distinct
weights, so the optimum was unique, the beam's tie-break was never invoked, and
an implementation keying block RNG streams on iteration order passed it. Both
`SYN-01` and `SYN-04` therefore ship degenerate configurations.
"""

from __future__ import annotations

import json

import pytest

from codonforge import _core
from codonforge.testkit import synthetic as syn

# --------------------------------------------------------------------------- #
# TST-SYN-001 — the toy exhaustive space
# --------------------------------------------------------------------------- #


def test_tst_syn_001_the_toy_space_is_the_size_the_plan_states() -> None:
    """8 designable codons, 3 synonyms each, 6 561 sequences.

    The reduced table matters: under table 1 only isoleucine is exactly 3-fold
    degenerate, so a genuine `3⁸` space cannot be built without restricting.
    """
    space = syn.syn01_distinct()
    assert space.size == 6561
    assert len(list(space.enumerate())) == 6561
    assert len(set(space.enumerate())) == 6561, "the enumeration repeats a sequence"


def test_tst_syn_001_every_enumerated_sequence_is_the_same_protein() -> None:
    """A design space that changed the protein would not be a design space."""
    space = syn.syn01_distinct()
    proteins = {_core.translate(cds, 1) for cds in space.enumerate()}
    assert proteins == {space.protein}


def test_tst_syn_001_the_optimum_is_computable_by_brute_force() -> None:
    space = syn.syn01_distinct()
    best = max(space.enumerate(), key=space.score)
    assert space.optima() == [best]
    assert space.score(best) == pytest.approx(space.optimum_value())

    # And it is genuinely the best, not merely a high scorer.
    assert all(space.score(cds) <= space.score(best) for cds in space.enumerate())


# --------------------------------------------------------------------------- #
# the tie coverage, which is the point
# --------------------------------------------------------------------------- #


def test_p1_32_the_toy_space_ships_a_degenerate_configuration() -> None:
    """DEP-SYN-001. Several distinct sequences attain the optimum exactly.

    The count is predicted from the construction — the product of the synonym
    counts of the tied residues — rather than read off the enumeration, so it
    is a ground truth rather than an observation of whatever came out.
    """
    space = syn.syn01_degenerate()
    optima = space.optima()

    assert space.degenerate
    assert len(optima) == 9, f"expected 3 x 3 optima, found {len(optima)}"
    assert len(set(optima)) == len(optima), "the optima are not distinct sequences"

    values = {space.score(cds) for cds in optima}
    assert len(values) == 1, "the 'tied' optima do not actually tie"
    assert all(space.score(cds) <= space.optimum_value() for cds in space.enumerate())


def test_p1_32_the_distinct_and_degenerate_configurations_really_differ() -> None:
    """Otherwise shipping both proves nothing."""
    assert len(syn.syn01_distinct().optima()) == 1
    assert len(syn.syn01_degenerate().optima()) > 1
    assert syn.syn01_distinct().size == syn.syn01_degenerate().size


def test_p1_32_a_tie_is_a_property_of_the_weights_not_of_the_protein() -> None:
    """Both configurations span the same space, so a search sees one variable."""
    distinct, degenerate = syn.syn01_distinct(), syn.syn01_degenerate()
    assert distinct.protein == degenerate.protein
    assert list(distinct.enumerate()) == list(degenerate.enumerate())


# --------------------------------------------------------------------------- #
# TST-SYN-004 — planted violations
# --------------------------------------------------------------------------- #


def _violations(case: syn.PlantedViolation, cds: str) -> list[dict]:
    report = json.loads(
        _core.constraint_violations(
            json.dumps(case.constraint_set), cds, [(0, len(cds))], "assembled"
        )
    )
    return [v for v in report["violations"] if v["actionable"]]


def test_tst_syn_004_planted_violations_are_actually_present() -> None:
    """Each generator produces a sequence with exactly the intended violation."""
    for build in (
        syn.syn04_unique_fix,
        syn.syn04_multiple_equal_fixes,
        syn.syn04_no_single_codon_fix,
    ):
        case = build()
        found = _violations(case, case.cds)
        assert found, f"{build.__name__} planted no detectable violation"
        assert all(v["constraint_id"] in json.dumps(case.constraint_set) for v in found)


def test_tst_syn_004_the_minimal_fixes_are_exactly_the_ones_claimed() -> None:
    """Every claimed fix works, and no unclaimed single substitution does."""
    for build in (syn.syn04_unique_fix, syn.syn04_multiple_equal_fixes):
        case = build()
        codons = [case.cds[i : i + 3] for i in range(0, len(case.cds), 3)]

        for index, codon in case.minimal_fixes:
            patched = codons.copy()
            patched[index] = codon
            assert not _violations(case, "".join(patched)), (
                f"{build.__name__} claims ({index}, {codon}) fixes it; it does not"
            )
            assert _core.translate("".join(patched), 1) == _core.translate(case.cds, 1)

        # Exhaustive: nothing else at distance 1 works.
        working = set()
        for index, codon in enumerate(codons):
            for synonym in _core.synonyms(codon, 1):
                if synonym == codon:
                    continue
                patched = codons.copy()
                patched[index] = synonym
                if not _violations(case, "".join(patched)):
                    working.add((index, synonym))
        assert working == set(case.minimal_fixes)


def test_tst_syn_004_the_unique_case_has_exactly_one_fix() -> None:
    """Constructed on tryptophan, which has one codon and so cannot move."""
    case = syn.syn04_unique_fix()
    assert case.fix_count == 1
    assert case.minimal_distance == 1


def test_p1_32_the_multi_fix_case_makes_tie_breaking_testable() -> None:
    """DEP-SYN-001 for repair.

    With one fix available, a repair that always returns the same answer is
    indistinguishable from one that chooses among equals — so the diversity
    property `TST-CONS-002` asserts would be unfalsifiable.
    """
    case = syn.syn04_multiple_equal_fixes()
    assert case.fix_count >= 3, f"only {case.fix_count} fixes; a tie needs several"
    assert len(set(case.minimal_fixes)) == case.fix_count
    assert all(distance == 1 for distance in [case.minimal_distance])


def test_p1_32_a_case_with_no_single_codon_fix_is_included() -> None:
    """ "Gave up" and "fixed it" are different outcomes, so both are covered."""
    case = syn.syn04_no_single_codon_fix()
    assert case.minimal_fixes == ()
    assert case.minimal_distance == 2

    codons = [case.cds[i : i + 3] for i in range(0, len(case.cds), 3)]
    for index, codon in enumerate(codons):
        for synonym in _core.synonyms(codon, 1):
            patched = codons.copy()
            patched[index] = synonym
            assert _violations(case, "".join(patched)), (
                f"({index}, {synonym}) fixes a case claiming no single-codon fix"
            )


# --------------------------------------------------------------------------- #
# SYN-05 — infeasible sets
# --------------------------------------------------------------------------- #


def test_p1_32_syn05_cases_behave_as_their_expected_code_says() -> None:
    """Including the case that must *not* be reported infeasible.

    A violation confined to immutable flank is non-actionable by
    `DEP-CONS-002`, and a suite that only contained genuine infeasibilities
    could not tell a correct analyser from one that escalates everything.
    """
    for case in syn.syn05_cases():
        construct = case.five_prime_flank + case.cds + case.three_prime_flank
        start = len(case.five_prime_flank)
        report = json.loads(
            _core.constraint_violations(
                json.dumps(case.constraint_set),
                construct,
                [(start, start + len(case.cds))],
                "assembled",
            )
        )
        actionable = [v for v in report["violations"] if v["actionable"]]

        if case.expected_code == "":
            assert not actionable, (
                f"{case.description}: a non-actionable violation was reported as one"
            )
            assert report["violations"], "the case planted no violation at all"
        else:
            assert actionable, f"{case.description}: expected an actionable violation"


# --------------------------------------------------------------------------- #
# SYN-07 — property strategies
# --------------------------------------------------------------------------- #


def test_p1_32_strategies_are_reproducible_from_their_seed() -> None:
    """A property failure that cannot be replayed is a report nobody can act on."""
    a = syn.Strategies(seed=7)
    b = syn.Strategies(seed=7)
    c = syn.Strategies(seed=8)

    first = [a.arbitrary_cds() for _ in range(20)]
    assert first == [b.arbitrary_cds() for _ in range(20)]
    assert first != [c.arbitrary_cds() for _ in range(20)], (
        "two seeds produced identical streams; the seed is not reaching the generator"
    )


def test_p1_32_generated_sequences_are_valid_cds() -> None:
    strategies = syn.Strategies(seed=11)
    for _ in range(50):
        cds = strategies.arbitrary_cds(length=8)
        assert len(cds) % 3 == 0
        protein = _core.translate(cds, 1)
        assert protein.startswith("M")
        assert protein.endswith("*")
        assert "*" not in protein[:-1], "an internal stop was generated"


def test_p1_32_generated_substitutions_are_synonymous() -> None:
    strategies = syn.Strategies(seed=13)
    for _ in range(50):
        cds = strategies.arbitrary_cds()
        substitution = strategies.arbitrary_substitution(cds)
        if substitution is None:
            continue
        index, codon = substitution
        codons = [cds[i : i + 3] for i in range(0, len(cds), 3)]
        assert codon != codons[index]
        codons[index] = codon
        assert _core.translate("".join(codons), 1) == _core.translate(cds, 1)


def test_p1_32_every_generator_is_reachable_from_the_registry() -> None:
    """So a caller can iterate the suite rather than have to name each one."""
    assert set(syn.GENERATORS) == {
        "SYN-01/distinct",
        "SYN-01/degenerate",
        "SYN-04/unique-fix",
        "SYN-04/multiple-equal-fixes",
        "SYN-04/no-single-codon-fix",
        "SYN-05",
    }
    for name, build in syn.GENERATORS.items():
        assert build() is not None, name
