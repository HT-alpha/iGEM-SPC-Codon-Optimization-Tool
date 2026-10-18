"""P2-10: the incremental accessibility cache.

`TST-STRU-002` (exactness) and `TST-STRU-006` (the sixth fixture-design rule).

The audit found `TST-STRU-002` — "incremental accessibility equals full
recomputation" — is satisfied perfectly by a cache that calls full recomputation
every time. It verifies correctness and says nothing about the thing being a
cache. So these tests come in three layers, and the second and third are the
ones the sixth rule adds:

1. **Agreement.** The incremental result equals a full recomputation, exactly.
2. **Argument fidelity.** The edit *range* actually determines what is
   recomputed. A cache that ignores the range and reuses everything is stale; one
   that ignores it and refolds everything is correct but is not a cache. Both are
   caught here and neither is caught by (1) alone.
3. **Plausibility.** Values stay inside the range a probability can occupy,
   recomputation is bounded by the folding span rather than the sequence length,
   and using the cache before priming it is a refusal rather than a plausible
   empty answer.
"""

from __future__ import annotations

import os
import random

import pytest

from dna_design_suite.codon.orchestrator.structure_diag import FoldSession
from dna_design_suite.codon.predictors import (
    StructureEngineError,
    available_engines,
    select_engine,
)
from dna_design_suite.codon.predictors.cache import AccessibilityCache
from dna_design_suite.codon.predictors.viennarna import accessibility_penalty

_ENGINE_PRESENT = "viennarna" in available_engines()

if os.environ.get("CF_REQUIRE_STRUCTURE_ENGINE") == "1" and not _ENGINE_PRESENT:
    raise RuntimeError(
        "CF_REQUIRE_STRUCTURE_ENGINE=1 but no structure engine is installed; "
        "the cache tests would silently skip and prove nothing about P2-10"
    )

vienna = pytest.mark.skipif(
    not _ENGINE_PRESENT, reason="ViennaRNA is a user-installed backend (AMD-013)"
)

WINDOW = 6
SPAN = 40


def cache_for(sequence: str) -> AccessibilityCache:
    cache = AccessibilityCache(
        engine=select_engine("viennarna"),
        window_nt=WINDOW,
        span_nt=SPAN,
        temperature_c=37.0,
    )
    cache.prime(sequence)
    cache.stats.reset()
    return cache


def random_sequence(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def full_recompute(sequence: str) -> list[float]:
    engine = select_engine("viennarna")
    return list(engine.local_accessibility(sequence, WINDOW, SPAN, 37.0))


# --------------------------------------------------------------------------- #
# 1. agreement
# --------------------------------------------------------------------------- #


@vienna
def test_tst_stru_002_incremental_equals_full_recomputation_over_many_edits():
    """Exactness under a walk, not under a single lucky edit.

    Edits are applied in sequence without re-priming, so an error that only
    appears after the cache has been updated several times — a drifting offset,
    a region boundary that creeps — has somewhere to accumulate. A test that
    re-primed between edits could not see that class of bug at all.
    """
    sequence = random_sequence(300, seed=11)
    cache = cache_for(sequence)
    rng = random.Random(4)

    for _ in range(40):
        start = rng.randrange(0, len(sequence) - 3)
        start -= start % 3  # codon-aligned, as a synonymous substitution is
        replacement = "".join(rng.choice("ACGT") for _ in range(3))

        incremental = cache.substitute(start, replacement)
        sequence = sequence[:start] + replacement + sequence[start + 3 :]
        expected = full_recompute(sequence)

        assert len(incremental) == len(expected)
        worst = max(abs(a - b) for a, b in zip(incremental, expected, strict=True))
        assert worst < 1e-9, f"incremental drifted from full recomputation by {worst}"


# --------------------------------------------------------------------------- #
# 2. argument fidelity — the edit range must matter
# --------------------------------------------------------------------------- #


@vienna
def test_tst_stru_006_the_edit_position_determines_which_windows_are_recomputed():
    """A cache that ignores where the edit was is either stale or not a cache.

    Two edits far apart must touch disjoint window ranges. If the computed range
    is the same for both, the range argument is being ignored — and an agreement
    test would still pass as long as the ignored range was "everything".
    """
    sequence = random_sequence(400, seed=3)
    cache = cache_for(sequence)

    near_start = cache.affected_window_range(30, 33)
    near_end = cache.affected_window_range(350, 353)

    assert near_start != near_end, "the edit position must change which windows move"
    assert near_start[1] <= near_end[0], (
        f"edits 320 nt apart with a {SPAN} nt span must affect disjoint windows; "
        f"got {near_start} and {near_end}"
    )


@vienna
def test_tst_stru_006_a_distant_window_is_unchanged_by_an_edit():
    """Locality is the claim the whole cache rests on.

    If accessibility were not local, reusing a distant value would be an
    approximation rather than an exact reuse — and the cache would be quietly
    wrong in a way agreement-on-the-happy-path might not surface, because the
    error would be small.

    So this checks the claim directly against the engine, not against the cache:
    a full recomputation before and after an edit must agree far from it.
    """
    sequence = random_sequence(400, seed=5)
    before = full_recompute(sequence)
    edited = sequence[:200] + "GGG" + sequence[203:]
    after = full_recompute(edited)

    far = len(before) - 1  # last window, ~200 nt from the edit
    assert abs(before[far] - after[far]) < 1e-12, (
        "a window 200 nt from a 3 nt edit changed; if this fails the locality "
        "assumption behind the cache is wrong and the cache is unsound, not slow"
    )

    # And the edit did change something, so the test above is not vacuous.
    assert any(abs(a - b) > 1e-9 for a, b in zip(before, after, strict=True)), (
        "the edit changed nothing anywhere, so this fixture cannot detect a locality failure"
    )


@vienna
def test_tst_stru_006_the_cache_actually_reuses_and_actually_recomputes():
    """Both counters must be non-zero. Either alone is a different bug.

    Zero reuse means it is not a cache. Zero recomputation means it is serving
    stale values. `TST-STRU-002` catches the second and is blind to the first.
    """
    cache = cache_for(random_sequence(400, seed=9))
    cache.substitute(200, "GGG")

    assert cache.stats.windows_reused > 0, (
        "nothing was reused, so this refolds everything on every edit and the "
        "cache is a wrapper with counters"
    )
    assert cache.stats.windows_recomputed > 0, (
        "nothing was recomputed, so the edit was never folded and every value "
        "served after it is stale"
    )


# --------------------------------------------------------------------------- #
# 3. plausibility
# --------------------------------------------------------------------------- #


@vienna
def test_tst_stru_006_recomputation_is_bounded_by_the_span_not_the_length():
    """The property that makes it worth having.

    Recomputation should depend on the folding span, which is fixed, and not on
    the sequence length, which is not. A cache whose work grows with the
    construct is a full refold wearing a hat.
    """
    counts = []
    for length in (200, 400, 800):
        cache = cache_for(random_sequence(length, seed=length))
        cache.substitute(length // 2, "GGG")
        counts.append(cache.stats.windows_recomputed)

    assert len(set(counts)) == 1, (
        f"recomputation grew with sequence length: {counts}. It should be a "
        f"function of the {SPAN} nt span alone"
    )
    # And reuse must grow, or the constant above would be consistent with the
    # cache simply doing nothing.
    reuse = []
    for length in (200, 800):
        cache = cache_for(random_sequence(length, seed=length))
        cache.substitute(length // 2, "GGG")
        reuse.append(cache.stats.windows_reused)
    assert reuse[1] > reuse[0] * 2, f"reuse should scale with length, got {reuse}"


@vienna
def test_tst_stru_006_every_value_stays_inside_the_range_a_probability_can_occupy():
    """Confident nonsense is usually out of range.

    An index shift, a misaligned region offset or an off-by-one in the reuse
    loop tends to produce values that are still finite and still plausible one
    at a time. Checking the range catches the subset that escapes, and costs
    nothing.
    """
    sequence = random_sequence(300, seed=13)
    cache = cache_for(sequence)
    rng = random.Random(21)

    for _ in range(20):
        start = rng.randrange(0, len(sequence) - 3)
        values = cache.substitute(start, "".join(rng.choice("ACGT") for _ in range(3)))
        assert all(0.0 <= v <= 1.0 for v in values), (
            f"accessibility outside [0, 1]: min {min(values)}, max {max(values)}"
        )


@vienna
def test_tst_stru_006_an_unprimed_cache_refuses_rather_than_returning_nothing():
    cache = AccessibilityCache(
        engine=select_engine("viennarna"),
        window_nt=WINDOW,
        span_nt=SPAN,
        temperature_c=37.0,
    )
    with pytest.raises(StructureEngineError) as caught:
        cache.current()
    assert caught.value.code == "CF-STRU-005"

    with pytest.raises(StructureEngineError):
        cache.substitute(0, "GGG")


@vienna
def test_tst_stru_006_an_edit_outside_the_sequence_is_refused():
    cache = cache_for(random_sequence(100, seed=2))
    with pytest.raises(StructureEngineError) as caught:
        cache.substitute(99, "GGG")
    assert caught.value.code == "CF-STRU-004"


@vienna
def test_fold_session_reuses_a_duplicate_and_edits_one_codon() -> None:
    """A batch of copies must not fold three times.

    Removing the stored-penalty hit makes `prime` run once per copy, and this
    assertion is what notices. The one-codon neighbour must match a full fold
    of that neighbour, not a stale copy of the parent.
    """
    sequence = random_sequence(240, seed=5)
    start = 60
    replacement = "GGG" if sequence[start : start + 3] != "GGG" else "AAA"
    edited = sequence[:start] + replacement + sequence[start + 3 :]
    calls = {"prime": 0}
    original = AccessibilityCache.prime

    def counted(self: AccessibilityCache, seq: str) -> list[float]:
        calls["prime"] += 1
        return original(self, seq)

    AccessibilityCache.prime = counted  # type: ignore[method-assign]
    try:
        session = FoldSession()
        values = session.evaluate([sequence, sequence, edited])
    finally:
        AccessibilityCache.prime = original  # type: ignore[method-assign]

    assert calls["prime"] == 1
    assert session.hits == 1
    assert session.local_edits == 1
    assert values[0] == values[1]
    assert values[0] == accessibility_penalty(full_recompute(sequence))
    assert values[2] == accessibility_penalty(full_recompute(edited))


@vienna
def test_fold_session_folds_a_distant_sequence_from_scratch() -> None:
    sequence = random_sequence(240, seed=5)
    # Two codon edits is past the length//160 budget of 1, so this is a full fold.
    far = sequence[:3] + "GGG" + sequence[6:30] + "CCC" + sequence[33:]
    if far[3:6] == sequence[3:6]:
        far = sequence[:3] + "AAA" + far[6:]
    session = FoldSession()
    values = session.evaluate([sequence, far])
    assert session.full_folds == 2
    assert session.local_edits == 0
    assert values[1] == accessibility_penalty(full_recompute(far))
