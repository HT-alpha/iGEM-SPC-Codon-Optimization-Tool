"""TST-CLI-001 through TST-CLI-004: the assembled pipeline, CLI and library.

Ticket P1-30.

`TST-CLI-003` is the load-bearing one. `TST-CLI-001` asserts the three surfaces
agree, and they share one orchestrator, so three surfaces that all skip a tier
agree perfectly and it passes. What has to be asserted separately is that the
pipeline *ran* what the budget planned — and the evidence is the `CF-ORCH-006`
warnings disappearing, because their presence in `P1-29` was the evidence it
had not.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

from test_orchestrator import CDS, FIVE_PRIME, HOST, THREE_PRIME, constraint_set, host_profile

from codonforge import _core
from codonforge.cli import main as cli_main
from codonforge.modelreg import ArtifactStore
from codonforge.orchestrator import DesignError, DesignRequest, design

SEED = 20260828


@pytest.fixture()
def store(tmp_path: Path) -> ArtifactStore:
    return ArtifactStore(tmp_path / "artifacts")


def request_for(
    objectives: list[str] | None = None,
    budget: str = "standard",
    generations: int = 8,
    **overrides: object,
) -> DesignRequest:
    return DesignRequest(
        host=HOST,
        cds=CDS,
        host_profile=host_profile(),
        # Two conflicting objectives by default: with one, tier 1's exact DP is
        # already at the optimum, tier 3's population converges to a single
        # sequence, and every assertion about the population loop is vacuous.
        objectives=objectives or ["usage.cai", "usage.gc_content"],
        constraints=constraint_set(),
        five_prime_flank=FIVE_PRIME,
        three_prime_flank=THREE_PRIME,
        translation_table=1,
        budget=budget,
        seed=SEED,
        generations=generations,
        segment=4,
        host_digest="sha256:" + "1" * 64,
        job_id="job_01JBQ8Z4K7M2N3P5R6T7V8W9XY",
        **overrides,  # type: ignore[arg-type]
    )


def gc(sequence: str) -> float:
    return sum(base in "GC" for base in sequence) / len(sequence)


# --------------------------------------------------------------------------- #
# TST-CLI-003 — the pipeline runs what the budget planned
# --------------------------------------------------------------------------- #


def test_tst_cli_003_the_planned_tiers_run_and_stop_being_reported_as_skipped(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """The CF-ORCH-006 placeholders disappear for a job that runs the tiers.

    `P1-29` reported `tier1_seed` and `tier2_block_ca` as `skipped` with a
    warning naming what was missing. If those warnings survive `P1-30`, the
    wiring is present and not real.
    """
    result = design(request_for(), store)
    status = result.status

    states = {tier["tier"]: tier["state"] for tier in status["progress"]["tiers"]}
    assert states["tier1_seed"] == "done"
    assert states["tier2_block_ca"] == "done"
    assert states["tier3_moea"] == "done"

    placeholders = [w for w in status["warnings"] if w["code"] == "CF-ORCH-006"]
    assert placeholders == [], (
        "a tier the budget planned is still reported as unassembled: "
        f"{[w['subject'] for w in placeholders]}"
    )

    # Present is not the same as effective: the tiers must have done something.
    assert result.job.stage_reports["tier1_seed"]["seeds"], "tier 1 produced no seeds"
    assert result.job.stage_reports["tier2_block_ca"]["sweeps"] >= 1
    assert result.job.seeded_cds != CDS, "the pipeline returned its input unchanged"


def test_tst_cli_003_tier_2_is_not_a_no_op(store: ArtifactStore) -> None:
    """Tier 2 improves a sequence it is given unpolished.

    Separated from the pipeline test on purpose. In the assembled pipeline tier
    2 usually commits nothing, because tier 1's DP is already optimal for an
    additive objective — so a pipeline-level assertion that tier 2 changed
    something would fail for a correct implementation, and one that asserted
    nothing would pass for a tier 2 that did not run at all.
    """
    probe = json.dumps(
        {
            "cds": CDS,
            "translation_table": 1,
            "five_prime_flank": FIVE_PRIME,
            "three_prime_flank": THREE_PRIME,
            "constraint_set": constraint_set(),
            "host_profile": host_profile(),
            "host_id": HOST,
            "host_digest": "sha256:" + "1" * 64,
            "terms": ["usage.cai"],
            "weights": [1.0],
            "root_seed": SEED,
        }
    )
    before = json.loads(_core.search_tier2_polish(probe, 4, 0))
    after = json.loads(_core.search_tier2_polish(probe, 4, 6))

    assert after["blocks_committed"] > 0, "tier 2 committed no block from an unpolished start"
    assert after["objective"] > before["objective"]
    assert after["incumbent"]["cds"] != CDS


def test_tst_cli_003_a_fast_budget_runs_fewer_tiers_and_says_so(
    store: ArtifactStore,
) -> None:
    """The plan is honoured, not ignored in either direction."""
    result = design(request_for(budget="fast"), store)
    states = {t["tier"]: t["state"] for t in result.status["progress"]["tiers"]}
    assert set(states) == {"tier1_seed", "tier2_block_ca"}, (
        "a fast budget must not silently acquire the population loop"
    )
    assert all(state == "done" for state in states.values())


def test_tst_cli_003_an_exact_budget_still_warns_about_the_decision_stage(
    store: ArtifactStore,
) -> None:
    """Assembling two tiers must not silence the warning about a third.

    `decision` is `P2-08`'s. If wiring tiers 1 and 2 had removed the mechanism
    rather than satisfied it for those two, this would be green and wrong.
    """
    result = design(request_for(budget="exact"), store)
    placeholders = [w for w in result.status["warnings"] if w["code"] == "CF-ORCH-006"]
    assert [w["subject"] for w in placeholders] == ["decision"]
    assert placeholders[0]["stage"] == "tier_plan"


# --------------------------------------------------------------------------- #
# TST-CLI-004 — objective direction survives end to end
# --------------------------------------------------------------------------- #


def test_tst_cli_004_a_minimize_term_is_actually_minimized(store: ArtifactStore) -> None:
    """The search must move the way the document says it moved.

    This failed when it was first written. Tier 1 and tier 2 each apply
    `provider.direction()` themselves, and the job layer signed the weights as
    well — so `usage.gc_content`, declared `minimize`, came out maximized while
    `objective_sense` still read `minimize`.
    """
    result = design(request_for(objectives=["usage.gc_content"]), store)

    assert result.front["objective_sense"] == ["minimize"]
    assert result.status["resolved_objectives"][0]["direction"] == "minimize"

    # Asserted on the *tiers'* own output, not only on the front. The front
    # reports the best candidate found, and the generation loop wanders far
    # enough that "better than the input" is satisfied even by a pipeline
    # optimizing in the wrong direction - a signed-weight defect passed that
    # threshold. What has to hold is that tier 1 and tier 2 each moved the way
    # the provider declared.
    tier1 = result.job.stage_reports["tier1_seed"]["seeds"][0]["cds"]
    tier2 = result.job.stage_reports["tier2_block_ca"]["incumbent"]["cds"]
    assert gc(tier1) < gc(CDS), (
        f"tier 1 did not minimize a minimize term: GC went {gc(CDS):.3f} -> {gc(tier1):.3f}"
    )
    assert gc(tier2) <= gc(tier1), (
        f"tier 2 moved a minimize term the wrong way: {gc(tier1):.3f} -> {gc(tier2):.3f}"
    )

    produced = result.front["candidates"][0]["sequences"][0]["cds"]
    assert gc(produced) <= gc(tier1), "the front reported a worse candidate than tier 1 produced"


def test_tst_cli_004_a_maximize_term_is_actually_maximized(store: ArtifactStore) -> None:
    """The other half, so the test above cannot pass by a blanket sign flip."""
    result = design(request_for(objectives=["usage.cai"]), store)
    assert result.front["objective_sense"] == ["maximize"]

    scores = json.loads(
        _core.term_scores(
            json.dumps(
                {
                    "cds": CDS,
                    "translation_table": 1,
                    "five_prime_flank": FIVE_PRIME,
                    "three_prime_flank": THREE_PRIME,
                    "constraint_set": constraint_set(),
                    "host_profile": host_profile(),
                    "host_id": HOST,
                    "host_digest": "sha256:" + "1" * 64,
                    "terms": ["usage.cai"],
                    "weights": [1.0],
                    "root_seed": SEED,
                }
            ),
            CDS,
        )
    )
    produced = result.front["candidates"][0]["objective_vector"][0]
    assert produced > scores["usage.cai"], (
        f"a maximize term was not maximized: {scores['usage.cai']:.4f} -> {produced:.4f}"
    )


def test_tst_cli_004_directions_come_from_the_providers(store: ArtifactStore) -> None:
    """Both senses appear in one front, in the objective order."""
    result = design(request_for(objectives=["usage.cai", "usage.gc_content"]), store)
    assert result.front["objective_order"] == ["usage.cai", "usage.gc_content"]
    assert result.front["objective_sense"] == ["maximize", "minimize"]
    assert [o["direction"] for o in result.status["resolved_objectives"]] == [
        "maximize",
        "minimize",
    ]


def test_tst_cli_004_a_signed_weight_is_refused(store: ArtifactStore) -> None:
    """Applying the direction twice is refused rather than silently obeyed."""
    with pytest.raises(ValueError, match="magnitudes"):
        design(request_for(weights=[-1.0]), store)


# --------------------------------------------------------------------------- #
# TST-CLI-001 — the surfaces agree
# --------------------------------------------------------------------------- #


def test_tst_cli_001_cli_and_library_produce_the_same_front(
    store: ArtifactStore, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Same request, same seed, same front — through both surfaces."""
    request = request_for(objectives=["usage.cai", "usage.gc_content"])
    library = design(request, store)

    path = tmp_path / "request.json"
    path.write_text(request.to_json(), encoding="utf-8")
    code = cli_main(["design", "--json", str(path), "--output", "json", "--store", str(store.root)])
    assert code == 0
    emitted = json.loads(capsys.readouterr().out)

    for key in ("objective_order", "objective_sense", "candidates", "reduction"):
        assert emitted["front"][key] == library.front[key]
    assert (
        emitted["status"]["landscape_lockfile_digest"]
        == (library.status["landscape_lockfile_digest"])
    )
    assert emitted["status"]["progress"] == library.status["progress"]
    assert emitted["status"]["warnings"] == library.status["warnings"]


def test_tst_cli_001_the_front_document_validates(store: ArtifactStore) -> None:
    jsonschema = pytest.importorskip("jsonschema")
    from referencing import Registry, Resource

    schema_dir = Path(__file__).resolve().parents[2] / "schemas" / "v1"
    registry = Registry().with_resources(
        [
            (p.name, Resource.from_contents(json.loads(p.read_text())))
            for p in schema_dir.glob("*.schema.json")
        ]
    )
    schema = json.loads((schema_dir / "pareto_front.schema.json").read_text())

    result = design(request_for(objectives=["usage.cai", "usage.gc_content"]), store)
    jsonschema.Draft202012Validator(schema, registry=registry).validate(result.front)


def test_tst_cli_001_the_same_seed_reproduces_and_a_different_one_explores_differently(
    store: ArtifactStore,
) -> None:
    """Same seed, same everything. Different seed, different search.

    The seed sensitivity is asserted on the *trajectory*, not on the reported
    candidate, and that is not a weaker claim — it is the accurate one. Now
    that the front reports the best candidate rather than the last, tier 1's DP
    dominates it, and that DP is exact for an additive objective: it finds the
    same optimum from any seed, using randomness only to break ties. Requiring
    two seeds to disagree on the answer would be requiring the optimizer to be
    worse than it is.
    """
    a = design(request_for(), store)
    again = design(request_for(), store)
    other = request_for()
    other.seed = SEED + 1
    b = design(other, store)

    assert a.front["candidates"] == again.front["candidates"]
    assert a.job.checkpoint["members"] == again.job.checkpoint["members"]

    assert b.job.checkpoint["members"] != a.job.checkpoint["members"], (
        "two root seeds produced an identical search trajectory; the seed is not "
        "reaching the engine"
    )


# --------------------------------------------------------------------------- #
# TST-CLI-002 — exit codes
# --------------------------------------------------------------------------- #


def test_tst_cli_002_exit_codes_match_the_documented_table(
    store: ArtifactStore, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Infeasible constraints exit 3 and not 1, so a script can tell them apart."""
    constraints = tmp_path / "c.json"
    constraints.write_text(json.dumps(constraint_set()), encoding="utf-8")

    # 0 — a clean run.
    profile = tmp_path / "h.json"
    profile.write_text(json.dumps(host_profile()), encoding="utf-8")
    assert (
        cli_main(
            [
                "design",
                "--host",
                HOST,
                "--cds",
                CDS,
                "--host-profile",
                str(profile),
                "--constraints",
                str(constraints),
                "--generations",
                "4",
                "--store",
                str(store.root),
            ]
        )
        == 0
    )

    # 2 — invalid input.
    assert cli_main(["design", "--cds", CDS, "--store", str(store.root)]) == 2
    assert "CF-CLI-001" in capsys.readouterr().err

    # 3 — infeasible: the construct violates its own constraint set.
    bad = tmp_path / "bad.json"
    bad.write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "constraint_set_id": "cst_bad",
                "constraints": [
                    {
                        "id": "no_atg",
                        "class": "blocking",
                        "scope": "cds",
                        "rule": {
                            "kind": "forbidden_motif",
                            "motif": "ATG",
                            "strand": "sense",
                            "max_occurrences": 0,
                        },
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    assert cli_main(["feasibility", "--cds", CDS, "--constraints", str(bad)]) == 3


def test_tst_cli_002_background_is_refused_rather_than_faked(store: ArtifactStore) -> None:
    """DEP-CLI-004. Returning a finished handle would be worse than refusing."""
    with pytest.raises(DesignError) as excinfo:
        design(request_for(), store, background=True)
    assert excinfo.value.code == "CF-CLI-003"
    assert "P1-29" in excinfo.value.detail


def test_tst_cli_002_an_unknown_objective_is_refused(store: ArtifactStore) -> None:
    """A job that asked for four objectives must not quietly optimize three."""
    with pytest.raises(DesignError) as excinfo:
        design(request_for(objectives=["usage.cai", "kinetics.dwell"]), store)
    assert excinfo.value.code == "CF-LAND-001"
    assert "kinetics.dwell" in excinfo.value.detail


# --------------------------------------------------------------------------- #
# checkpoint/resume and cancellation, through the surface
# --------------------------------------------------------------------------- #


def test_p1_30_resume_through_the_library_matches_an_uninterrupted_run(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """DEP-CLI-001, and P1-29's determinism claim through the public surface."""
    straight = design(request_for(generations=12), store)

    checkpoints = tmp_path / "ckpt"
    stopped = design(
        request_for(generations=12, checkpoint_dir=str(checkpoints)),
        store,
        on_segment=lambda job: setattr(job, "cancel_requested", True),
    )
    assert stopped.status["state"] == "cancelled"
    assert stopped.job.generations_run == 4

    resumed = design(
        request_for(generations=12, checkpoint_dir=str(checkpoints), resume=True), store
    )
    assert resumed.status["state"] == "succeeded"
    # Compared against the uninterrupted run rather than a fixed count: tier 3
    # terminates when it converges, so the number is a property of the search
    # and asserting 12 would be asserting that convergence never happens.
    assert resumed.job.generations_run == straight.job.generations_run
    assert resumed.job.checkpoint == straight.job.checkpoint, (
        "a run resumed through the library diverged from an uninterrupted one"
    )
    assert resumed.front["candidates"] == straight.front["candidates"]


def test_p1_30_cancellation_through_the_surface_stops_cleanly(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """The same clean segment-boundary stop P1-29 validated internally."""
    checkpoints = tmp_path / "ckpt"
    result = design(
        request_for(generations=12, checkpoint_dir=str(checkpoints)),
        store,
        on_segment=lambda job: setattr(job, "cancel_requested", True),
    )

    assert result.status["state"] == "cancelled"
    assert result.status["problem"]["code"] == "CF-ORCH-003"
    assert result.front["summary"]["termination_reason"] == "user_cancel"

    document = json.loads(
        (checkpoints / "job_01JBQ8Z4K7M2N3P5R6T7V8W9XY.checkpoint.json").read_text()
    )
    assert document["generations_run"] == result.job.generations_run
    assert document["budget_consumed"] == result.job.budget.consumed
    assert document["state"]["generation"] == result.job.generations_run


def test_p1_30_resume_without_a_checkpoint_is_refused(store: ArtifactStore, tmp_path: Path) -> None:
    with pytest.raises(DesignError) as excinfo:
        design(request_for(checkpoint_dir=str(tmp_path / "empty"), resume=True), store)
    assert excinfo.value.code == "CF-ORCH-001"


# --------------------------------------------------------------------------- #
# progress structure and warnings provenance through the surface
# --------------------------------------------------------------------------- #


def test_p1_30_progress_stays_per_tier_through_the_cli(
    store: ArtifactStore, capsys: pytest.CaptureFixture[str]
) -> None:
    """Not flattened to a percentage for display convenience.

    A job that ran tier 1 and skipped tier 2 and one that ran both are
    different outcomes, and one number cannot say which happened.
    """
    request = request_for()
    cli_main(["design", "--json", _write(request, store), "--store", str(store.root)])
    printed = capsys.readouterr().out

    assert "tier1_seed" in printed
    assert "tier2_block_ca" in printed
    assert "tier3_moea" in printed
    # And only the generational tier shows a generation count.
    lines = {line.split()[0]: line for line in printed.splitlines() if line.startswith("  tier")}
    assert "/" in lines["tier3_moea"]
    assert "/" not in lines["tier1_seed"], (
        "tier 1 was reported with a generation count it does not have"
    )


def test_p1_30_warning_provenance_survives_the_cli(
    store: ArtifactStore, capsys: pytest.CaptureFixture[str]
) -> None:
    """DEP-ORCH-001's (code, stage, subject) is not collapsed for display."""
    request = request_for(budget="exact")
    code = cli_main(
        ["design", "--json", _write(request, store), "--output", "json", "--store", str(store.root)]
    )
    assert code == 0
    warnings = json.loads(capsys.readouterr().out)["status"]["warnings"]

    assert warnings, "the operator-mix warnings did not reach the CLI"
    for warning in warnings:
        assert warning["stage"]
        assert warning["subject"]
        assert warning["code"].startswith("CF-")
    assert {w["subject"] for w in warnings} >= {"potts_guided", "block_gibbs", "decision"}


def _write(request: DesignRequest, store: ArtifactStore) -> str:
    path = Path(store.root).parent / "request.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(request.to_json(), encoding="utf-8")
    return str(path)


# --------------------------------------------------------------------------- #
# the pipeline is not monotone, and says so
# --------------------------------------------------------------------------- #


def test_p1_30_the_front_reports_the_best_candidate_not_the_last(
    store: ArtifactStore,
) -> None:
    """Found by a negative control that would not fire.

    The generation loop proposes and repairs but does not select — selection is
    NSGA's, which is `P2-02` — so running it after an optimizing tier moves away
    from the optimum rather than past it. Here tier 1 takes CAI from 0.25 to
    0.92 and the loop returns 0.62. Reporting the loop's output published that
    loss as the answer, and made every direction assertion in this file a
    statement about a random walk rather than about the search.
    """
    result = design(request_for(objectives=["usage.cai"], generations=12), store)
    probe = json.dumps(
        {
            "cds": CDS,
            "translation_table": 1,
            "five_prime_flank": FIVE_PRIME,
            "three_prime_flank": THREE_PRIME,
            "constraint_set": constraint_set(),
            "host_profile": host_profile(),
            "host_id": HOST,
            "host_digest": "sha256:" + "1" * 64,
            "terms": ["usage.cai"],
            "weights": [1.0],
            "root_seed": SEED,
        }
    )
    score = lambda cds: json.loads(_core.term_scores(probe, cds))["usage.cai"]  # noqa: E731

    reported = result.front["candidates"][0]["sequences"][0]["cds"]
    tier1_output = result.job.stage_reports["tier1_seed"]["seeds"][0]["cds"]

    # The claim has changed shape, and the change is the point. P1-30 found
    # that the objective-free mutation loop moved *away* from tier 1's optimum
    # and the front reported the loss. P2-04 replaced that loop with tier 3,
    # which is elitist - NSGA-II keeps the best front - so the pipeline no
    # longer degrades and there is nothing for the archive to compensate for
    # here. What must still hold is that the reported candidate is no worse
    # than any stage produced, which is the guarantee rather than the symptom.
    for candidate in result.job.archive:
        if candidate.feasible:
            assert score(reported) >= score(candidate.cds) - 1e-12, (
                f"{candidate.origin} produced a better candidate than the one reported"
            )
    assert score(reported) >= score(tier1_output)


def test_p1_30_a_non_monotone_stage_is_reported_rather_than_hidden(
    store: ArtifactStore,
) -> None:
    """CF-ORCH-007, through the same provenance path as every other warning."""
    result = design(request_for(objectives=["usage.cai"], generations=12), store)
    notes = [w for w in result.status["warnings"] if w["code"] == "CF-ORCH-007"]

    # Tier 3 is elitist, so a later stage no longer regresses and the warning
    # correctly does not fire. The mechanism is not gone - it is exercised
    # directly in `test_archive.py` against constructed regressions, which is
    # where it belongs now that the pipeline does not produce one. Asserting
    # its absence here keeps the claim honest: if a future stage does regress,
    # this test says nothing and the archive tests still do.
    assert notes == [], (
        "a stage regressed. That is permitted - exploration is not monotone - "
        "but it means this pipeline is no longer the monotone one this test "
        "documents, and the change deserves a look rather than a passing test"
    )
    assert len(notes) == 0
    # And the reported candidate now comes from the population loop rather than
    # from an earlier stage the archive had to rescue.
    assert result.front["candidates"][0]["provenance"]["origin"] in (
        "tier1_dp",
        "tier2_block_ca",
        "mutation",
    )


def test_p1_30_the_reported_candidate_is_never_worse_than_any_stage_produced(
    store: ArtifactStore,
) -> None:
    """The archive guarantee, observed end to end through the real pipeline.

    `test_archive.py` proves the rule holds for arbitrary stage orders; this
    checks that the pipeline actually routes every stage's output through it,
    which is the half a unit test of the archive cannot see.
    """
    result = design(request_for(objectives=["usage.cai"], generations=12), store)
    probe = json.dumps(
        {
            "cds": CDS,
            "translation_table": 1,
            "five_prime_flank": FIVE_PRIME,
            "three_prime_flank": THREE_PRIME,
            "constraint_set": constraint_set(),
            "host_profile": host_profile(),
            "host_id": HOST,
            "host_digest": "sha256:" + "1" * 64,
            "terms": ["usage.cai"],
            "weights": [1.0],
            "root_seed": SEED,
        }
    )

    recorded = list(result.job.archive)
    assert len(recorded) >= 4, (
        f"only {len(recorded)} candidates recorded; the stages are not all routed "
        "through the archive, so the guarantee covers less than it claims"
    )
    assert {c.origin for c in recorded} >= {"input", "tier1_dp", "tier2_block_ca", "mutation"}

    reported = result.front["candidates"][0]["sequences"][0]["cds"]
    reported_score = _core.scalarized_score(probe, reported)
    for candidate in recorded:
        if candidate.feasible:
            assert reported_score >= candidate.score - 1e-12, (
                f"{candidate.origin} produced a better candidate "
                f"({candidate.score:.6g}) than the one reported ({reported_score:.6g})"
            )

    # And the reported one is feasible, which outranks score.
    assert result.front["candidates"][0]["constraint_report"]["feasible"]
    assert result.job.is_feasible(reported)


def test_p1_30_an_infeasible_input_is_never_reported_over_a_feasible_design(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """Feasibility outranks score, through the pipeline rather than in isolation.

    The input is recorded like any other candidate, so a constraint set the
    input violates is exactly the case where scoring alone would report an
    infeasible sequence as the design.
    """
    request = request_for(objectives=["usage.cai"])
    request.constraints = {
        "schema_version": "1.0.0",
        "constraint_set_id": "cst_input_violates",
        "constraints": [
            {
                "id": "no_ctgcgt",
                "class": "blocking",
                "scope": "cds",
                "rule": {
                    "kind": "forbidden_motif",
                    "motif": "CTGCGT",
                    "strand": "sense",
                    "max_occurrences": 0,
                },
            }
        ],
    }
    result = design(request, store)

    infeasible = [c for c in result.job.archive if not c.feasible]
    assert infeasible, "the fixture's input must actually violate its constraint set"
    assert {c.origin for c in infeasible} >= {"input"}

    assert result.front["candidates"][0]["constraint_report"]["feasible"]
    assert result.job.best_candidate().origin != "input"


# --------------------------------------------------------------------------- #
# TST-CLI-008 — tier 3 through the actual job and CLI surfaces
# --------------------------------------------------------------------------- #


def test_tst_cli_008_a_real_job_runs_tier_3_and_its_front_reaches_job_status(
    store: ArtifactStore,
) -> None:
    """The pipeline-level test the archive-level ones cannot stand in for.

    `P2-01` and `P2-02` prove the machinery directly, and a pipeline that never
    called it would pass every one of them — the same reason `P1-30` needed
    `TST-CLI-003` beyond `P1-29`'s segment-level tests.
    """
    result = design(request_for(objectives=["usage.cai", "usage.gc_content"]), store)

    # Tier 3 actually ran, and produced a population rather than an incumbent.
    assert result.job.tier3 is not None, "the job never called tier 3"
    assert len(result.job.checkpoint["members"]) == result.job.inputs.population
    assert result.job.tier3["selection"] == "crowding_distance"

    # Its front reached job_status through the schema's own structure.
    summary = result.status["current_front_summary"]
    assert summary["n_nondominated"] == len(result.job.tier3["front"])
    assert summary["n_nondominated"] > 0
    assert len(summary["hypervolume_history"]) == result.job.generations_run
    assert summary["hypervolume"] == summary["hypervolume_history"][-1]

    # `best_per_objective` is a per-axis extremum, not a selected candidate:
    # every value must be attained by some front member on that axis.
    for index, term in enumerate(result.front["objective_order"]):
        attained = {m["objectives"][index] for m in result.job.tier3["front"]}
        assert summary["best_per_objective"][term] == min(attained)


def test_tst_cli_008_the_front_reaches_the_cli_output(
    store: ArtifactStore, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    from codonforge.cli import main as cli_main

    path = tmp_path / "request.json"
    path.write_text(request_for().to_json(), encoding="utf-8")
    assert (
        cli_main(["design", "--json", str(path), "--output", "json", "--store", str(store.root)])
        == 0
    )
    emitted = json.loads(capsys.readouterr().out)

    summary = emitted["status"]["current_front_summary"]
    assert summary["n_nondominated"] > 0
    assert summary["hypervolume_history"]
    assert emitted["front"]["candidates"], "the front document has no candidates"


def test_tst_cli_008_a_job_terminates_on_stagnation_rather_than_the_budget(
    store: ArtifactStore,
) -> None:
    """Termination reaches the job, not just the engine.

    The budget is deliberately far larger than the search needs, so a run that
    used all of it would be one where termination never fired.
    """
    request = request_for(generations=200)
    request.segment = 4
    result = design(request, store)

    assert result.status["state"] == "succeeded"
    assert result.job.generations_run < 200, (
        "the job ran its whole budget; hypervolume-stagnation termination did not "
        "reach the job layer"
    )
    assert result.job.tier3["terminated"] is True
    assert result.job.tier3["termination"] == "hypervolume_stagnation"


def test_tst_cli_008_the_archive_holds_every_front_member_not_one_representative(
    store: ArtifactStore,
) -> None:
    """CF-ORCH-007 receives the whole front.

    Recording a single representative would discard the trade-off before
    anything downstream could see it, which is the failure mode the front/best
    split exists to prevent.
    """
    result = design(request_for(objectives=["usage.cai", "usage.gc_content"]), store)

    from_loop = {c.cds for c in result.job.archive if c.origin == "mutation"}
    reported = {m["cds"] for m in result.job.tier3["front"]}
    assert reported <= from_loop, "front members never reached the archive"
    assert len(from_loop) > 1, (
        "only one candidate came from the population loop, so this cannot tell "
        "'records the front' from 'records a representative'"
    )

    # And every archived candidate carries its objective vector, so the
    # archive's front() has coordinates to work with rather than falling back
    # to the scalarization.
    for candidate in result.job.archive:
        if candidate.origin == "mutation":
            assert len(candidate.objectives) == 2
