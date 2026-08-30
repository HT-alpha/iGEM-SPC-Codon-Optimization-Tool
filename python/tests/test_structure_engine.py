"""P2-09 and P2-08: the structure seam and the ViennaRNA binding.

`TST-STRU-001`, `TST-STRU-003`, `TST-STRU-004`. `AMD-013`, `AMD-001`,
`DEP-STRU-002`.

The sequencing audit rewrote both stated criteria, because both tested the wrong
thing once `AMD-013` made ViennaRNA the engine rather than the reference:

* `TST-STRU-001` compared "the native kernel" against an installed ViennaRNA.
  With ViennaRNA *as* the engine that compares it to itself.
* `TST-STRU-004` required agreement with "an independent brute-force partition
  function". CodonForge has no parameter set of its own, so no such reference
  can be built, and a test that pretended to would be fabricating exactly the
  thing `DEP-STRU-001` exists to prevent.

What replaces them is shaped around a specific question: **which bugs would
produce plausible numbers?** A binding that returns something is easy. A binding
that returns the right thing has to survive an index check, an
argument-fidelity check (does it actually use the window and temperature it was
given, or silently use defaults?), and physical checks that would fail on
confident nonsense.
"""

from __future__ import annotations

import math
import os

import pytest

from codonforge.predictors import (
    StructureEngineError,
    available_engines,
    select_engine,
)

# `TST-STRU-001` says "skipped, not failed, when ViennaRNA is absent" - correct
# for a contributor's laptop, and a trap for CI. A suite where every reference
# check skips is green and proves nothing about the binding, which is the exact
# shape of the counting bug found two rounds ago. So CI sets
# CF_REQUIRE_STRUCTURE_ENGINE=1 and the skip becomes a failure.
_ENGINE_PRESENT = "viennarna" in available_engines()
_ENGINE_REQUIRED = os.environ.get("CF_REQUIRE_STRUCTURE_ENGINE") == "1"

if _ENGINE_REQUIRED and not _ENGINE_PRESENT:
    raise RuntimeError(
        "CF_REQUIRE_STRUCTURE_ENGINE=1 but no structure engine is installed. "
        "These checks are the only thing standing between the ViennaRNA binding "
        "and a plausible-looking wrong answer; skipping them in CI would make the "
        "suite green without testing P2-08 at all. Install with "
        "`pip install -e '.[structure-tests]'`."
    )

vienna = pytest.mark.skipif(
    not _ENGINE_PRESENT,
    reason="ViennaRNA is not installed; it is a user-installed backend (AMD-013)",
)

# A hairpin with a GC-rich stem, then an unstructured tail. Constructed rather
# than sampled: the fixture rules require a stochastic-looking condition to be
# built deliberately, and "this region is paired" must be true by design, not by
# luck of the folding.
STEM = "GGGGGGAAAACCCCCC"
TAIL = "A" * 20
HAIRPIN = STEM + TAIL


@vienna
def test_tst_stru_001_the_binding_matches_a_directly_computed_reference():
    """The binding's numbers are ViennaRNA's numbers, correctly extracted.

    The reference is recomputed here from `pfl_fold_up` with the 1-based
    indexing written out longhand. That is deliberately not a call to the same
    helper: what this catches is an off-by-one in the extraction, which shifts
    every accessibility by a nucleotide and which nothing downstream would find
    surprising.
    """
    import RNA

    window, span, temperature = 4, len(HAIRPIN), 37.0
    engine = select_engine("viennarna")
    produced = engine.local_accessibility(HAIRPIN, window, span, temperature)

    RNA.cvar.temperature = temperature
    matrix = RNA.pfl_fold_up(HAIRPIN.replace("T", "U"), window, span, span)
    expected = [matrix[end][window] for end in range(window, len(HAIRPIN) + 1)]

    assert len(produced) == len(expected)
    for index, (got, want) in enumerate(zip(produced, expected, strict=True)):
        assert got == pytest.approx(want, rel=0, abs=0), (
            f"position {index} disagrees with the reference: {got} vs {want}"
        )


@vienna
def test_tst_stru_001_the_binding_actually_uses_the_arguments_it_is_given():
    """A passthrough that ignores its arguments would pass a self-comparison.

    This is the check the reference test above cannot make: if the binding
    silently folded at a default temperature or a default window, comparing it
    against a reference computed with *the same defaults* would agree perfectly.
    So each argument is varied and required to change the answer.
    """
    engine = select_engine("viennarna")

    base = engine.local_accessibility(HAIRPIN, 4, len(HAIRPIN), 37.0)
    wider = engine.local_accessibility(HAIRPIN, 8, len(HAIRPIN), 37.0)
    hotter = engine.local_accessibility(HAIRPIN, 4, len(HAIRPIN), 85.0)

    assert len(wider) == len(HAIRPIN) - 8 + 1 != len(base), (
        "a different window must change how many windows fit; if the length is "
        "unchanged the window argument is being ignored"
    )
    assert base[:5] != hotter[:5], (
        "temperature must reach the partition function; if it does not, every "
        "run folds at whatever temperature the library last saw"
    )


@vienna
def test_tst_stru_004_accessibility_discriminates_structure_from_no_structure():
    """A physical property, not an agreement.

    This one holds no matter which engine is bound and would fail against a
    binding that returned confident nonsense — including one that returned a
    constant, which every agreement-style check above would still pass if the
    reference were computed the same wrong way.
    """
    engine = select_engine("viennarna")
    accessible = engine.local_accessibility("A" * len(HAIRPIN), 4, len(HAIRPIN), 37.0)
    folded = engine.local_accessibility(HAIRPIN, 4, len(HAIRPIN), 37.0)

    assert min(accessible) > 0.99, (
        "an unstructured poly-A has nothing to pair with, so every window should "
        "be essentially certain to be unpaired"
    )
    assert min(folded) < 1e-3, (
        "a GC-rich stem should be strongly paired; if the minimum accessibility "
        "over a hairpin is high, the engine is not folding anything"
    )
    assert min(folded) < min(accessible)


@vienna
def test_tst_stru_004_accessibility_is_monotone_in_temperature():
    """Heat opens structure. Three points, strictly increasing.

    Two points could pass on noise. The stem is probed rather than the tail
    because the tail is saturated at 1.0 and cannot move - a fixture that
    measured there would be unable to fail.
    """
    engine = select_engine("viennarna")
    stem_end = len(STEM) - 4  # index of a window lying inside the stem

    readings = [
        min(engine.local_accessibility(HAIRPIN, 4, len(HAIRPIN), t)[:stem_end])
        for t in (10.0, 37.0, 90.0)
    ]
    assert readings[0] < readings[1] < readings[2], (
        f"stem accessibility should rise with temperature, got {readings}"
    )


@vienna
def test_tst_stru_004_an_impossible_window_is_absent_rather_than_reported_as_zero():
    """ViennaRNA returns exactly 0.0 for a run that cannot fit.

    Passing that through would report the 5' edge as maximally *inaccessible*
    when the truth is that no window of that length ends there. The binding
    returns one value per window end instead, so the question is never answered
    for positions where it was not asked.
    """
    engine = select_engine("viennarna")
    window = 6
    values = engine.local_accessibility(HAIRPIN, window, len(HAIRPIN), 37.0)

    assert len(values) == len(HAIRPIN) - window + 1
    # Index k is the run ending at sequence position k + window - 1, so the last
    # entry ends at the final base. Anchoring both ends pins the alignment.
    tail_windows = values[-4:]
    assert all(v > 0.99 for v in tail_windows), (
        "the last windows lie entirely in the poly-A tail and must be accessible; "
        "if they are not, the returned values are misaligned with the sequence"
    )


@vienna
def test_tst_stru_004_the_penalty_orders_designs_the_way_the_objective_needs():
    from codonforge.predictors.viennarna import accessibility_penalty

    engine = select_engine("viennarna")
    folded = accessibility_penalty(engine.local_accessibility(HAIRPIN, 4, len(HAIRPIN), 37.0))
    open_ = accessibility_penalty(
        engine.local_accessibility("A" * len(HAIRPIN), 4, len(HAIRPIN), 37.0)
    )

    assert folded > open_, "the term is minimized, so structure must cost more"
    assert math.isfinite(folded) and math.isfinite(open_)


@vienna
def test_tst_stru_001_the_engine_identity_reaches_the_lockfile_shape():
    engine = select_engine("viennarna")
    identity = engine.identity().as_document()
    assert identity["name"] == "viennarna"
    assert identity["version"], "an unversioned engine cannot be reproduced"
    assert identity["parameter_set"], (
        "DEP-STRU-001's reproducibility half: which parameter table produced "
        "these energies has to be recorded even when licensing is settled"
    )


def test_tst_stru_003_with_no_engine_installed_the_job_is_refused_not_degraded(
    monkeypatch,
):
    """`AMD-013` removed the native fallback, so absence must be visible.

    The original criterion said the native engine is selected when nothing is
    installed. There is no native engine now, and the replacement behaviour is
    a refusal: a proxy would return a number for every candidate, so every
    downstream check would pass while the objective meant nothing.
    """
    import codonforge.predictors.structure as seam

    monkeypatch.setattr(seam, "available_engines", lambda: [])
    with pytest.raises(StructureEngineError) as caught:
        seam.select_engine()

    assert caught.value.code == "CF-STRU-003"
    assert "ViennaRNA" in caught.value.detail, (
        "the refusal must name what to install; a user who hits this is being "
        "asked to go and do something"
    )


def test_tst_stru_003_no_import_of_an_absent_package_is_attempted():
    """Detection uses `find_spec`, not a `try: import`.

    A `try: import` runs the package's top-level code whenever it happens to be
    importable, which is how a detection routine acquires side effects.
    """
    from pathlib import Path

    source = Path(seam_path()).read_text(encoding="utf-8")
    assert "find_spec" in source
    detection = source[source.index("def available_engines") : source.index("def select_engine")]
    assert "import RNA" not in detection, (
        "engine detection must not import the engine; that is what find_spec is for"
    )


def seam_path() -> str:
    import codonforge.predictors.structure as seam

    return seam.__file__


def test_tst_stru_003_a_requested_but_absent_engine_is_named_in_the_refusal(
    monkeypatch,
):
    import codonforge.predictors.structure as seam

    monkeypatch.setattr(seam, "available_engines", lambda: ["viennarna"])
    with pytest.raises(StructureEngineError) as caught:
        seam.select_engine("linearfold")
    assert caught.value.code == "CF-STRU-003"
    assert "linearfold" in caught.value.detail


def test_amd_001_the_engine_offers_no_global_mfe_to_the_search_loop():
    """`AMD-001` keeps global MFE out of the loop, enforced by absence.

    A `global_mfe` method would eventually be called from a per-move path by
    someone with a good local reason. The cheapest enforcement is to give the
    loop nothing to call, and the cheapest way to keep that true is a test that
    fails when the method reappears.
    """
    from codonforge.predictors.viennarna import ViennaRNAEngine

    assert not hasattr(ViennaRNAEngine, "global_mfe"), (
        "AMD-001 excludes global MFE from the search loop; if an exact-tier "
        "re-ranking path needs it, it belongs on a separate interface that the "
        "per-move path cannot reach"
    )
