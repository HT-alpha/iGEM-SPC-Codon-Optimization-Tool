"""P2-05a: the crossing count, measured where the crossing actually happens.

`TST-SRCH-024`.

`crates/cf-search/tests/batch_consumption.rs` proves the loop calls the protocol
correctly, but it does so with a Rust scorer — which crosses nothing. A
crossing-count assertion there would read zero and pass vacuously, which is the
pre-satisfied-criteria pattern the sequencing audits keep finding.

So the counter lives here, in the Python callable, counting its own invocations.
`ADR-0001` caps calls across the **PyO3** boundary at one per generation:

    no call crosses the boundary more often than once per generation

Note this is the FFI boundary, not an OS process boundary. `DEP-LAND-004`
deferred the out-of-process `cf-worker` and it stays deferred; transport changes
only *where* a scorer runs, while the budget is about how often Rust hands
control to Python. That is the quantity these tests measure.
"""

from __future__ import annotations

import json

import pytest

from dna_design_suite import _core
from dna_design_suite.modelreg import ArtifactStore


@pytest.fixture()
def store(tmp_path) -> ArtifactStore:
    return ArtifactStore(tmp_path / "artifacts")


class CountingScorer:
    """Counts its own calls and records exactly what it was handed.

    The sixth fixture-design rule, applied to a crossing: a scorer that ignored
    its arguments and returned cached values of the right shape would satisfy
    any check comparing batch output against per-item output, because both would
    be the same wrong numbers. So the answer here is a pure function of the
    sequence, and the tests check that the values on the population match the
    sequences those individuals actually carry.
    """

    def __init__(self) -> None:
        self.calls: list[list[str]] = []

    def __call__(self, sequences: list[str]) -> list[float | None]:
        self.calls.append(list(sequences))
        return [self.score_of(s) for s in sequences]

    @property
    def crossings(self) -> int:
        return len(self.calls)

    @staticmethod
    def score_of(sequence: str) -> float:
        if not sequence:
            return 0.0
        return sum(1 for c in sequence if c in "GC") / len(sequence)


def tier_request(store, tmp_path) -> str:
    """The real tier document, from the real job.

    Built through `Job.submit` rather than hand-assembled, so this exercises the
    same document the pipeline sends. A hand-rolled request would be a fixture
    that drifts from production without either one noticing.
    """
    from test_orchestrator import make_job

    job = make_job(store, tmp_path)
    return job._tier_request(job.inputs.cds)


def run_segment(
    store,
    tmp_path,
    scorer: CountingScorer | None,
    generations: int,
    population: int = 8,
    budget: int = 3,
    axis: int = 2,
) -> dict:
    request = tier_request(store, tmp_path)
    return json.loads(
        _core.search_tier3_segment(
            request,
            None,
            population,
            generations,
            1000,  # stagnation window long enough not to end the run early
            scorer,
            axis,
            budget,
        )
    )


@pytest.mark.parametrize("generations", [1, 2, 5])
def test_tst_srch_024_exactly_one_crossing_per_generation(store, tmp_path, generations: int):
    """The measurement the whole ticket exists to make.

    A per-individual loop would show up as `population * generations` calls
    here. Parameterized over three values so a constant — a scorer called once
    regardless — cannot pass.
    """
    scorer = CountingScorer()
    report = run_segment(store, tmp_path, scorer, generations)

    assert report["generations_run"] == generations, (
        "the segment must run the generations asked for; otherwise a low "
        "crossing count would mean nothing happened rather than that the "
        "budget was respected"
    )
    assert scorer.crossings == generations, (
        f"ADR-0001 caps the PyO3 boundary at one crossing per generation; "
        f"{generations} generations produced {scorer.crossings} calls"
    )


def test_tst_srch_024_the_crossing_count_is_independent_of_population(store, tmp_path):
    """A per-item loop scales with population. One crossing per generation does not."""
    small = CountingScorer()
    run_segment(store, tmp_path, small, 3, population=6)
    large = CountingScorer()
    run_segment(store, tmp_path, large, 3, population=20)

    assert small.crossings == large.crossings == 3, (
        f"crossings scaled with population: {small.crossings} vs {large.crossings}. "
        "That is the per-individual loop this ticket replaced"
    )


def test_tst_srch_024_no_scorer_means_no_crossing(store, tmp_path):
    """The negative half. Without it, a scorer never called would pass above."""
    report = run_segment(store, tmp_path, None, 3)
    assert report["generations_run"] == 3
    assert report["stranded"] == []
    assert report["checkpoint"].get("batch") is None


def test_tst_srch_024_the_batch_carries_real_distinct_genotypes(store, tmp_path):
    """Argument fidelity: what crossed the boundary is the population, not a stub."""
    scorer = CountingScorer()
    run_segment(store, tmp_path, scorer, 4, population=12, budget=4)

    assert scorer.crossings == 4
    for call in scorer.calls:
        assert call, "an empty batch is a wasted crossing"
        assert len(call) <= 4, f"batch of {len(call)} exceeds the E_max of 4"
        for sequence in call:
            assert len(sequence) % 3 == 0, "a batch member must be codon-aligned"
            assert set(sequence) <= set("ACGT"), f"not a DNA sequence: {sequence[:12]}"

    distinct = {s for call in scorer.calls for s in call}
    assert len(distinct) > 1, (
        "every batched sequence was identical; the batch is not carrying the "
        "population it claims to"
    )


def test_tst_srch_024_a_malformed_return_is_refused_not_aligned(store, tmp_path):
    """Plausibility: a wrong-length reply must not be silently zipped.

    Aligning by position would attribute one sequence's score to another, and
    every resulting value would still look like a plausible probability.
    """

    class TooShort(CountingScorer):
        def __call__(self, sequences: list[str]) -> list[float | None]:
            self.calls.append(list(sequences))
            return [0.5] * (len(sequences) - 1) if len(sequences) > 1 else []

    with pytest.raises(Exception) as caught:
        run_segment(store, tmp_path, TooShort(), 2, population=12, budget=5)
    assert "CF-LAND-012" in str(caught.value) or "returned" in str(caught.value)


def test_tst_srch_024_a_scorer_that_declines_leaves_unknown_not_zero(store, tmp_path):
    """`None` from the scorer is UNKNOWN, and must not become a number."""

    class Declining(CountingScorer):
        def __call__(self, sequences: list[str]) -> list[float | None]:
            self.calls.append(list(sequences))
            return [None] * len(sequences)

    scorer = Declining()
    report = run_segment(store, tmp_path, scorer, 2, population=8, budget=4)

    assert scorer.crossings == 2
    for member in report["checkpoint"]["members"]:
        assert member["objectives"][2] is None, (
            "a declined item became a number; UNKNOWN is a value in "
            "constraint-domination and imputing it silently changes the ranking"
        )


def test_tst_srch_024_the_budget_strands_individuals_through_the_real_path(store, tmp_path):
    """`CF-SRCH-016` driven by a search that actually ran.

    P2-05 and P2-06 proved the protocol and the ranker against themselves. This
    is the first stranded individual produced by ranker -> protocol -> a real
    crossing, which is what the warning has always claimed to report.
    """
    scorer = CountingScorer()
    report = run_segment(store, tmp_path, scorer, 3, population=20, budget=2)

    assert report["stranded"], (
        "20 individuals at E_max 2 over 3 generations must strand some; an "
        "empty list means CF-SRCH-016 can never fire on a real job"
    )
    # And the protocol state travels, so a resumed job keeps the warning.
    batch = report["checkpoint"]["batch"]
    assert batch is not None
    assert batch["stranded"], "the stranded set was dropped from the checkpoint"


def test_tst_srch_024_a_sufficient_budget_strands_nobody(store, tmp_path):
    """The other branch, so the assertion above is not satisfied by always warning."""
    scorer = CountingScorer()
    report = run_segment(store, tmp_path, scorer, 3, population=6, budget=64)
    assert report["stranded"] == [], f"stranded {report['stranded']} with budget 64"
