"""The reported candidate never regresses, whatever the stage order.

Ticket P1-30, corrected. Code: `CF-ORCH-007`.

The point of these tests is that they say nothing about which stages exist.
They construct arbitrary stage sequences, including orderings no pipeline
currently produces, and require the same answer from all of them. A test that
asserted "tier 1 beats the generation loop" would pass for an implementation
that special-cased those two names, and would say nothing about the ordering
`P2-02` introduces.
"""

from __future__ import annotations

import hashlib
import itertools
import random

import pytest

from codonforge.orchestrator.archive import ORIGINS, CandidateArchive

# Origins used as arbitrary stage labels. Which ones they are does not matter,
# and that is the point.
STAGES = ("input", "tier1_dp", "tier2_block_ca", "mutation", "crossover", "repair")


def digest(text: str) -> str:
    """A digest-shaped label. Content matters only in being distinct.

    `hashlib` rather than `hash`, which is seeded per process for strings: a
    fixture whose tie-break inputs changed between runs would make a failure
    reproduce only sometimes.
    """
    return "sha256:" + hashlib.sha256(text.encode()).hexdigest()


def archive_of(entries: list[tuple[str, float, bool]]) -> CandidateArchive:
    archive = CandidateArchive()
    for origin, score, feasible in entries:
        label = f"{origin}:{score}:{feasible}"
        archive.record(origin, f"CDS_{label}", score, feasible, digest(label))
    return archive


# --------------------------------------------------------------------------- #
# stage-order independence
# --------------------------------------------------------------------------- #


def test_the_best_candidate_is_independent_of_stage_order() -> None:
    """Every permutation of the same stage outputs reports the same candidate.

    This is the guarantee stated without reference to any particular pipeline:
    the archive is a function of the *set* of candidates, not of the sequence.
    """
    entries = [
        ("input", 0.10, True),
        ("tier1_dp", 0.92, True),
        ("tier2_block_ca", 0.88, True),
        ("mutation", 0.62, True),
        ("crossover", 0.41, True),
    ]

    reported = {archive_of(list(order)).best().cds for order in itertools.permutations(entries)}
    assert len(reported) == 1, (
        f"the reported candidate depends on the order stages ran: {sorted(reported)}"
    )
    assert archive_of(entries).best().origin == "tier1_dp"


def test_a_later_stage_producing_a_worse_candidate_is_never_reported() -> None:
    """The constructed case the decision asks for, at every insertion point.

    A regression is placed after each prefix of an otherwise-improving run, so
    the assertion does not depend on the regression being last, first, or in
    any particular position.
    """
    improving = [
        ("input", 0.10, True),
        ("tier1_dp", 0.50, True),
        ("tier2_block_ca", 0.90, True),
    ]
    regression = ("mutation", 0.20, True)

    for position in range(len(improving) + 1):
        entries = improving[:position] + [regression] + improving[position:]
        best = archive_of(entries).best()
        assert best.score == 0.90, (
            f"a regression inserted at position {position} was reported: "
            f"{best.origin} at {best.score}"
        )
        assert best.origin == "tier2_block_ca"


def test_the_guarantee_holds_for_stage_orders_no_pipeline_produces() -> None:
    """Randomised orders over randomised scores, including future orderings.

    `P2-02` will add selection and may run stages in an order nothing here
    anticipates. The archive is checked against arbitrary orders rather than
    against the one this round happens to build.
    """
    rng = random.Random(20260828)
    for _ in range(200):
        entries = [
            (rng.choice(STAGES), round(rng.uniform(-5, 5), 6), rng.random() > 0.25)
            for _ in range(rng.randint(1, 8))
        ]
        shuffled = list(entries)
        rng.shuffle(shuffled)

        a, b = archive_of(entries).best(), archive_of(shuffled).best()
        assert (a.score, a.feasible, a.digest) == (b.score, b.feasible, b.digest), (
            f"shuffling {entries} changed the reported candidate"
        )


# --------------------------------------------------------------------------- #
# feasibility outranks score
# --------------------------------------------------------------------------- #


def test_a_feasible_candidate_beats_a_higher_scoring_infeasible_one() -> None:
    """A sequence that violates a blocking constraint is not a better answer.

    Selection was on score alone at first, and the *input* is recorded like any
    other candidate — so an infeasible input that happened to score well could
    be reported as the design. Reporting an infeasible sequence is worse than
    reporting a suboptimal one.
    """
    archive = archive_of(
        [
            ("input", 99.0, False),
            ("tier1_dp", 0.01, True),
        ]
    )
    best = archive.best()
    assert best.feasible
    assert best.origin == "tier1_dp"


def test_among_infeasible_candidates_the_best_scoring_is_still_chosen() -> None:
    """When nothing is feasible the job still reports its best attempt."""
    archive = archive_of([("input", 0.2, False), ("tier1_dp", 0.7, False)])
    best = archive.best()
    assert not best.feasible
    assert best.score == 0.7


# --------------------------------------------------------------------------- #
# ties
# --------------------------------------------------------------------------- #


def test_a_tie_is_broken_deterministically_and_not_by_position() -> None:
    """Genuine ties must not resolve by list order.

    A tie broken by position makes the answer depend on the order stages were
    recorded in, which is the order-dependence the archive exists to remove —
    so the tie-break would have quietly reintroduced it in the one case where
    it is hardest to notice.
    """
    tied = [
        ("tier1_dp", 0.75, True),
        ("tier2_block_ca", 0.75, True),
        ("mutation", 0.75, True),
    ]
    reported = {archive_of(list(order)).best().cds for order in itertools.permutations(tied)}
    assert len(reported) == 1, f"a tie resolved by position: {sorted(reported)}"

    # And the rule is the stated one: smallest digest.
    archive = archive_of(tied)
    assert archive.best().digest == min(c.digest for c in archive)


def test_ties_are_broken_the_same_way_across_processes() -> None:
    """The tie-break is a pure function of the digest string.

    Not of a hash, a set iteration, or anything else seeded per process — the
    same reason the tabu snapshot is sorted rather than left to `HashMap`
    order, which was a near-miss source of intermittent failures.
    """
    entries = [
        ("tier1_dp", 1.0, True, "sha256:" + "b" * 64),
        ("mutation", 1.0, True, "sha256:" + "a" * 64),
        ("crossover", 1.0, True, "sha256:" + "c" * 64),
    ]
    archive = CandidateArchive()
    for origin, score, feasible, dig in entries:
        archive.record(origin, "ATGTAA", score, feasible, dig)
    assert archive.best().digest == "sha256:" + "a" * 64


# --------------------------------------------------------------------------- #
# regression reporting, and round-tripping
# --------------------------------------------------------------------------- #


def test_a_regression_is_detected_without_naming_any_stage() -> None:
    archive = archive_of([("tier1_dp", 0.9, True), ("mutation", 0.3, True)])
    assert archive.regressed()
    assert archive.last_recorded().origin == "mutation"
    assert archive.best().origin == "tier1_dp"

    monotone = archive_of([("tier1_dp", 0.3, True), ("mutation", 0.9, True)])
    assert not monotone.regressed()


def test_an_empty_archive_reports_nothing_rather_than_guessing() -> None:
    archive = CandidateArchive()
    assert archive.best() is None
    assert archive.last_recorded() is None
    assert not archive.regressed()


def test_the_archive_round_trips_through_its_document() -> None:
    """It travels in the checkpoint, so it has to survive serialization.

    Three things the reported result depends on have already been found missing
    from the checkpoint, each by a different failing test. This one is checked
    before a fourth test has to find it.
    """
    archive = archive_of([("input", 0.1, True), ("tier1_dp", 0.9, True), ("mutation", 0.4, False)])
    restored = CandidateArchive.from_json(archive.to_json())

    assert [c.to_json() for c in restored] == [c.to_json() for c in archive]
    assert restored.best().to_json() == archive.best().to_json()
    assert restored.regressed() == archive.regressed()


def test_an_unknown_origin_is_refused() -> None:
    """The vocabulary is the front schema's, so a writer cannot be surprised."""
    archive = CandidateArchive()
    with pytest.raises(ValueError, match="provenance origin"):
        archive.record("tier9_magic", "ATGTAA", 1.0, True, digest("x"))
    assert set(ORIGINS) >= set(STAGES)


# --------------------------------------------------------------------------- #
# CF-ORCH-007 over a Pareto front, not a scalar winner
# --------------------------------------------------------------------------- #


def multi(entries: list[tuple[str, list[float | None], bool]]) -> CandidateArchive:
    """An archive of multi-objective candidates. All values minimized."""
    archive = CandidateArchive()
    for origin, objectives, feasible in entries:
        label = f"{origin}:{objectives}:{feasible}"
        archive.record(
            origin,
            f"CDS_{label}",
            score=0.0,
            feasible=feasible,
            digest=digest(label),
            objectives=tuple(objectives),
        )
    return archive


def test_the_front_is_the_non_dominated_set_not_a_scalar_winner() -> None:
    """The guarantee P2 exists to produce.

    Collapsing a trade-off curve to its scalar argmax hands a caller who asked
    for three objectives a single design, and nothing in the document says the
    rest were found and discarded.
    """
    archive = multi(
        [
            ("tier1_dp", [1.0, 5.0], True),
            ("tier2_block_ca", [3.0, 3.0], True),
            ("mutation", [5.0, 1.0], True),
            ("crossover", [9.0, 9.0], True),
        ]
    )
    front = archive.front()

    assert len(front) == 3, f"expected three mutually non-dominated members, got {len(front)}"
    assert {tuple(c.objectives) for c in front} == {(1.0, 5.0), (3.0, 3.0), (5.0, 1.0)}
    assert all(tuple(c.objectives) != (9.0, 9.0) for c in front), "a dominated member was reported"


def test_the_front_is_independent_of_record_order() -> None:
    entries = [
        ("tier1_dp", [1.0, 5.0], True),
        ("tier2_block_ca", [3.0, 3.0], True),
        ("mutation", [5.0, 1.0], True),
        ("crossover", [4.0, 4.0], True),
    ]
    reported = {
        tuple(c.digest for c in multi(list(order)).front())
        for order in itertools.permutations(entries)
    }
    assert len(reported) == 1, f"the reported front depends on record order: {reported}"


def test_a_feasible_front_excludes_infeasible_candidates_that_would_dominate() -> None:
    """Feasibility partitions the front exactly as it outranks score in `best`.

    An infeasible design is not a trade-off against a feasible one; it is not
    an answer. Its objectives dominating is beside the point.
    """
    archive = multi(
        [
            ("input", [0.0, 0.0], False),
            ("tier1_dp", [1.0, 5.0], True),
            ("mutation", [5.0, 1.0], True),
        ]
    )
    front = archive.front()

    assert len(front) == 2
    assert all(c.feasible for c in front)
    assert all(tuple(c.objectives) != (0.0, 0.0) for c in front)


def test_when_nothing_is_feasible_the_front_is_drawn_from_the_attempts() -> None:
    """A job that found no feasible design still reports its best attempts."""
    archive = multi(
        [
            ("input", [1.0, 5.0], False),
            ("mutation", [5.0, 1.0], False),
            ("crossover", [9.0, 9.0], False),
        ]
    )
    front = archive.front()

    assert len(front) == 2
    assert all(not c.feasible for c in front)


def test_an_unknown_term_makes_a_candidate_incomparable_rather_than_imputed() -> None:
    """§4's rule, at the reporting layer.

    A candidate whose provider missed its budget must not be dominated away on
    a value nobody measured, and must not dominate on one either. Imputing a
    mean or a worst case biases the front invisibly.
    """
    # The magnitudes are chosen so imputation *changes the answer*. A first
    # version used [9.0, None] against [1.0, 1.0]: imputing the unknown to zero
    # gives [9.0, 0.0], which still neither dominates nor is dominated, so the
    # front stayed at two and the test was blind to the imputation it was
    # written to catch. At [0.5, None] the imputed candidate dominates.
    archive = multi(
        [
            ("tier1_dp", [1.0, 1.0], True),
            ("mutation", [0.5, None], True),
        ]
    )
    front = archive.front()

    assert len(front) == 2, "an unknown term was imputed and the candidate dominated away"
    assert any(None in tuple(c.objectives) for c in front)

    # And the imputation would genuinely have changed it, so the assertion
    # above is not about an outcome that holds either way.
    assert all(value is not None for value in front[0].objectives) or all(
        value is not None for value in front[1].objectives
    )


def test_a_fully_tied_front_is_reported_whole() -> None:
    """Every member mutually non-dominated: none may be dropped."""
    archive = multi(
        [
            ("tier1_dp", [2.0, 2.0], True),
            ("tier2_block_ca", [2.0, 2.0], True),
            ("mutation", [2.0, 2.0], True),
        ]
    )
    assert len(archive.front()) == 3


def test_a_single_objective_archive_gives_a_front_of_one_agreeing_with_best() -> None:
    """The two answers agree where they should.

    Every Phase 1 job records a scalarization and no objective vector, so the
    front has to reduce to the same candidate `best` reports - otherwise the
    front document and the status document would disagree about the design.
    """
    archive = archive_of(
        [("input", 0.10, True), ("tier1_dp", 0.92, True), ("mutation", 0.62, True)]
    )
    front = archive.front()

    assert len(front) == 1
    assert front[0].digest == archive.best().digest
    assert front[0].origin == "tier1_dp"


def test_the_front_survives_the_checkpoint_round_trip() -> None:
    archive = multi(
        [
            ("tier1_dp", [1.0, 5.0], True),
            ("mutation", [5.0, 1.0], True),
            ("crossover", [9.0, 9.0], True),
        ]
    )
    restored = CandidateArchive.from_json(archive.to_json())

    assert [c.to_json() for c in restored.front()] == [c.to_json() for c in archive.front()]
    assert [tuple(c.objectives) for c in restored] == [tuple(c.objectives) for c in archive]
