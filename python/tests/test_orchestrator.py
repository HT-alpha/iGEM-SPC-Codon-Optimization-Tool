"""TST-ORCH-001, TST-ORCH-002, and the three carried-forward warning obligations.

Ticket P1-29.

The determinism claim is the one worth being careful about. "Checkpoint and
resume work" is easy to satisfy and hard to trust: a resumed job that restarted
its RNG stream counter at zero would be perfectly reproducible *against itself*
and quietly different from an uninterrupted run. So the test compares a resumed
run against an uninterrupted one of the same seed plan, state for state, rather
than checking that resume merely produces something.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from codonforge import _core
from codonforge.modelreg import ArtifactStore, Pin, ResolutionRequest, resolve, seal
from codonforge.orchestrator import ComputeBudget, Job, JobInputs, TierPlanError

HOST = "escherichia_coli:511145/BL21"

# `M L R S G V A P T K L R S *` — sixteen codons of real synonymous choice.
CDS = "ATGCTGCGTAGCGGCGTGGCGCCGACCAAACTGCGTAGCTAA"
FIVE_PRIME = "TTGACAATTAATCATCGGCTCGTATAATGTGTGGA"
THREE_PRIME = "AAAGCCCGGATCCTCTAGAGGG"

ROOT_SEED = 20260828


def constraint_set() -> dict:
    return {
        "schema_version": "1.0.0",
        "constraint_set_id": "cst_orch",
        "constraints": [
            {
                "id": "no_cgcagc",
                "class": "blocking",
                "scope": "construct",
                # Both strands, and a motif that is tripped in the *first*
                # generation rather than one the trajectory eventually
                # wanders into.
                #
                # This fixture has now gone vacuous twice. On the sense strand
                # alone the motif never arose at all, so repair never fired and
                # every assertion about the checkpoint carrying the tabu ledger
                # was empty - found by a negative control that would not fire.
                # Then DEP-DET-003 changed the generation loop's RNG root, the
                # trajectory moved, and a six-base motif that had been tripped
                # by luck around generation 18 stopped being tripped at all.
                # A four-base motif is violated early and often, so the fixture
                # no longer depends on where a particular random walk goes.
                "rule": {
                    "kind": "forbidden_motif",
                    "motif": "TTCG",
                    "strand": "both",
                    "max_occurrences": 0,
                },
            }
        ],
    }


def lock_scaffold(root_seed: int = ROOT_SEED) -> dict:
    return {
        "host_profile": {
            "host_id": HOST,
            "version": "1.0.0",
            "digest": "sha256:" + "1" * 64,
        },
        "inputs": {
            "construct_digest": "sha256:" + "2" * 64,
            "constraint_set_digest": "sha256:" + "3" * 64,
            "input_sequence_digests": ["sha256:" + "4" * 64],
        },
        "rng": {
            # xoshiro256pp, which is what SearchRng actually is since
            # DEP-DET-001. Naming an algorithm the engine does not implement
            # would make this block a description of something that did not
            # happen.
            "algorithm": "xoshiro256pp",
            "root_seed": root_seed,
            # The finest granularity the run uses; see the schema note.
            "stream_assignment": "per_block_index",
            "deterministic_reduction": True,
        },
        "environment": {
            "codonforge_version": "0.1.0",
            "git_sha": "0" * 40,
            "platform": "linux-x86_64",
            "thread_count": 1,
            "budget_tier": "standard",
        },
    }


def manifest_document(version: str = "1.0.0", note: str = "ci") -> dict:
    return {
        "schema_version": "1.0.0",
        "artifact_id": f"cf:epistasis_dwell:{HOST}:{version}",
        "kind": "epistasis_dwell",
        "host_id": HOST,
        "version": version,
        "contract_major": 1,
        "provides_terms": ["epistasis.dwell_deviation"],
        "created_at": "2026-08-28T00:00:00Z",
        "code_version": {
            "package": "codonforge",
            "version": "0.1.0",
            "git_sha": "0" * 40,
            "dirty": False,
        },
        "training": {
            "datasets": [{"dataset_id": "ds_test", "digest": "sha256:" + "9" * 64}],
            "holdout_datasets": [],
            "split_strategy": "by_gene",
            "seed": 1,
        },
        "evaluation": {"metrics": {"spearman": 0.61}, "gates": []},
        "isolation": {"scope": "public", "tenant_id": None, "may_export": True},
        "lifecycle": {"status": "active"},
        "created_by": note,
    }


@pytest.fixture()
def store(tmp_path: Path) -> ArtifactStore:
    return ArtifactStore(tmp_path / "artifacts")


def sealed_lockfile(store: ArtifactStore, root_seed: int = ROOT_SEED):
    manifest = store.publish(manifest_document())
    entries = resolve(
        store,
        ResolutionRequest(
            terms=["epistasis.dwell_deviation", "usage.cai"],
            host_id=HOST,
            contract_major=1,
            pins={"epistasis.dwell_deviation": Pin(manifest.artifact_id, manifest.digest)},
        ),
    )
    return seal(
        entries,
        lockfile_id="lock_01JBQ8Z4K7M2N3P5R6T7V8W9XZ",
        contract_major=1,
        **lock_scaffold(root_seed),
    )


def host_profile() -> dict:
    """A `highly_expressed` codon-usage table over the whole genetic code.

    Deterministic and non-uniform. A uniform table would make CAI indifferent
    between synonyms, so tier 1 and tier 2 would have nothing to optimize and
    every assertion about them changing the sequence would pass by doing
    nothing - the vacuity failure mode this project keeps finding. The
    assertion that the pipeline actually moves the sequence is in the tests
    below rather than here, but the fixture has to make movement possible.
    """
    bases = "TCAG"
    table = {}
    for i, a in enumerate(bases):
        for j, b in enumerate(bases):
            for k, c in enumerate(bases):
                index = i * 16 + j * 4 + k
                table[a + b + c] = {
                    "relative_adaptiveness": round(0.05 + 0.95 * ((index * 17) % 23) / 22.0, 4)
                }
    return {"codon_usage": {"subset": "highly_expressed", "table": table}}


def inputs() -> JobInputs:
    return JobInputs(
        cds=CDS,
        translation_table=1,
        five_prime_flank=FIVE_PRIME,
        three_prime_flank=THREE_PRIME,
        constraint_set=constraint_set(),
        codon_weights=[],
        tabu_leak=0.0,
        host_profile=host_profile(),
        host_id=HOST,
        host_digest="sha256:" + "1" * 64,
        tier1_beam=16,
        tier1_seeds=4,
        tier2_block=4,
        tier2_sweeps=3,
    )


def make_job(
    store: ArtifactStore,
    tmp_path: Path,
    *,
    total: int = 12,
    segment: int = 4,
    job_id: str = "job_01JBQ8Z4K7M2N3P5R6T7V8W9XY",
    terms: list[str] | None = None,
    capabilities: dict[str, str] | None = None,
    radii: dict[str, int | None] | None = None,
    block_size: int | None = None,
    directions: dict[str, str] | None = None,
    root_seed: int = ROOT_SEED,
) -> Job:
    return Job.submit(
        job_id=job_id,
        inputs=inputs(),
        budget=ComputeBudget(tier="standard", total_generations=total, segment=segment),
        lockfile=sealed_lockfile(store, root_seed),
        store=store,
        terms=terms if terms is not None else ["usage.cai"],
        capabilities=capabilities or {"usage.cai": "L1"},
        radii=radii or {"usage.cai": 0},
        block_size=block_size,
        directions=directions
        or {
            "usage.cai": "maximize",
            "usage.gc_content": "minimize",
            "lm.naturalness": "maximize",
            "structure.mfe_global": "minimize",
        },
        checkpoint_dir=tmp_path / "checkpoints",
    )


# --------------------------------------------------------------------------- #
# TST-ORCH-001 — the lockfile is sealed before any search work
# --------------------------------------------------------------------------- #


def test_tst_orch_001_an_unsealed_lockfile_cannot_start_a_job(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """A job that could seal later could also resolve later."""
    from codonforge.modelreg import Lockfile

    unsealed = Lockfile(
        lockfile_id="lock_01JBQ8Z4K7M2N3P5R6T7V8W9XZ",
        created_at="2026-08-28T00:00:00Z",
        contract_major=1,
        resolution=[],
        **lock_scaffold(),
    )
    assert unsealed.digest is None

    with pytest.raises(ValueError, match="sealed"):
        Job.submit(
            job_id="job_01JBQ8Z4K7M2N3P5R6T7V8W9XY",
            inputs=inputs(),
            budget=ComputeBudget(tier="standard", total_generations=4),
            lockfile=unsealed,
            store=store,
            terms=["usage.cai"],
            capabilities={"usage.cai": "L1"},
        )


def test_tst_orch_001_a_job_replays_its_lockfile_before_running(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """A missing pinned artifact stops the job before any generation runs."""
    job = make_job(store, tmp_path)
    digest = job.lockfile.entry("epistasis.dwell_deviation").digest

    # Remove the pinned artifact after submission.
    target = next(
        p for p in Path(store.root).rglob("*.json") if json.loads(p.read_text())["digest"] == digest
    )
    target.unlink()

    from codonforge.modelreg import RegistryError

    with pytest.raises(RegistryError) as excinfo:
        job.run()
    assert excinfo.value.code == "CF-RESOLVE-006"
    assert job.generations_run == 0, "no search work may have happened"
    assert job.checkpoint is None


# --------------------------------------------------------------------------- #
# TST-ORCH-002 — checkpoint and resume
# --------------------------------------------------------------------------- #


def test_tst_orch_002_a_resumed_run_matches_an_uninterrupted_one(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """Resume at generation k and finish; compare with never having stopped.

    The comparison is on the whole checkpoint — incumbent, tabu ledger and
    per-generation history — not on the final sequence alone. A run that
    diverged and reconverged would pass a final-state comparison.
    """
    uninterrupted = make_job(store, tmp_path / "a", total=12, segment=4).run()
    assert uninterrupted.state == "succeeded"
    assert uninterrupted.generations_run == 12

    # The fixture must exercise repair, or the claims about the checkpoint
    # carrying the tabu ledger across a resume are about an empty list.
    assert uninterrupted.checkpoint["repairs_attempted"] > 0, "repair never fired"

    # Stop after the first segment.
    stopped = make_job(store, tmp_path / "b", total=12, segment=3)

    def cancel_after_first(job: Job) -> None:
        if job.generations_run >= 3:
            job.cancel_requested = True

    stopped.run(on_segment=cancel_after_first)
    assert stopped.state == "cancelled"
    assert stopped.generations_run == 3

    # Resume from the persisted checkpoint into a fresh job object.
    document = Job.load_checkpoint(stopped.checkpoint_path())
    resumed = make_job(store, tmp_path / "b", total=12, segment=3).resume(document)
    resumed.run()

    assert resumed.state == "succeeded"
    assert resumed.generations_run == 12
    assert resumed.checkpoint == uninterrupted.checkpoint, (
        "a resumed run diverged from an uninterrupted one of the same seed plan"
    )
    # The ledger has to be live at some boundary, or "resume carries it" and
    # "resume drops it" are the same test. Which boundary depends on where
    # repair happens to fire, so it is found rather than assumed.
    assert any(h["tabu_live"] > 0 for h in uninterrupted.checkpoint["history"]), (
        "the tabu ledger was empty in every generation; the fixture proves nothing "
        "about the checkpoint carrying it"
    )


def test_tst_orch_002_the_history_is_identical_generation_by_generation(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """Not just the endpoint: every generation must match."""
    straight = make_job(store, tmp_path / "a", total=9, segment=9).run()
    assert straight.checkpoint["repairs_attempted"] > 0, "repair never fired"

    piecewise = make_job(store, tmp_path / "b", total=9, segment=3)
    for _ in range(3):
        document = None
        piecewise.run(on_segment=lambda job: setattr(job, "cancel_requested", True))
        document = Job.load_checkpoint(piecewise.checkpoint_path())
        piecewise = make_job(store, tmp_path / "b", total=9, segment=3).resume(document)
    assert piecewise.generations_run == 9

    assert straight.checkpoint["history"] == piecewise.checkpoint["history"]
    assert straight.checkpoint["incumbent"] == piecewise.checkpoint["incumbent"]
    assert straight.checkpoint["tabu"] == piecewise.checkpoint["tabu"]


def test_tst_orch_002_a_checkpoint_from_another_lockfile_is_refused(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """Resuming under a different lockfile would change which models ran."""
    first = make_job(store, tmp_path / "a", total=4, segment=4).run()
    document = Job.load_checkpoint(first.checkpoint_path())

    other = make_job(store, tmp_path / "b", total=4, segment=4, root_seed=999)
    with pytest.raises(ValueError, match="different landscape lockfile"):
        other.resume(document)

    mismatched = dict(document, job_id="job_01JBQ8Z4K7M2N3P5R6T7V8W9ZZ")
    with pytest.raises(ValueError, match="belongs to job"):
        make_job(store, tmp_path / "a", total=4, segment=4).resume(mismatched)


# --------------------------------------------------------------------------- #
# cancellation
# --------------------------------------------------------------------------- #


def test_p1_29_cancellation_stops_cleanly_and_leaves_a_usable_state(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """No partial state: the checkpoint is always a completed generation."""
    job = make_job(store, tmp_path, total=12, segment=3)
    job.run(on_segment=lambda j: setattr(j, "cancel_requested", True))

    assert job.state == "cancelled"
    assert job.generations_run == 3
    assert job.budget.consumed == 3
    # Consistency: the persisted checkpoint agrees with the in-memory job.
    document = Job.load_checkpoint(job.checkpoint_path())
    assert document["generations_run"] == job.generations_run
    assert document["budget_consumed"] == job.budget.consumed
    assert document["state"] == job.checkpoint
    assert document["state"]["generation"] == 3
    assert len(document["state"]["history"]) == 3

    status = job.status()
    assert status["state"] == "cancelled"
    assert status["problem"]["code"] == "CF-ORCH-003"
    loop = next(t for t in status["progress"]["tiers"] if t["tier"] == "tier3_moea")
    assert loop["generation"] == 3

    # And it is genuinely resumable.
    resumed = make_job(store, tmp_path, total=12, segment=3).resume(document)
    resumed.run()
    assert resumed.state == "succeeded"
    assert resumed.generations_run == 12


def test_p1_29_cancelling_before_the_first_segment_runs_nothing(
    store: ArtifactStore, tmp_path: Path
) -> None:
    job = make_job(store, tmp_path, total=8, segment=4)
    job.cancel_requested = True
    job.run()
    assert job.state == "cancelled"
    assert job.generations_run == 0
    assert job.checkpoint is None


# --------------------------------------------------------------------------- #
# the three carried-forward warning obligations
# --------------------------------------------------------------------------- #


def test_p1_29_operator_warnings_reach_job_status(store: ArtifactStore, tmp_path: Path) -> None:
    """P1-23's obligation, via P1-18b's plumbing, discharged here.

    `potts_guided` and `block_gibbs` are absent without a context model, and the
    mix renormalizes from 0.35/0.35/0.20/0.10 to 0.35/0.10 — a real change in
    how the search explores that nothing could previously show a user.
    """
    job = make_job(store, tmp_path, total=4, segment=4).run()
    status = job.status()

    operators = [w for w in status["warnings"] if w["code"] == "CF-SRCH-001"]
    assert len(operators) == 2, f"expected two absent operators, got {status['warnings']}"
    assert {w["subject"] for w in operators} == {"potts_guided", "block_gibbs"}
    for warning in operators:
        assert warning["stage"] == "operators"
        assert "P3-05" in warning["detail"]
        assert "renormalized" in warning["detail"]


def test_p1_29_tier_capability_warnings_reach_job_status(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """P1-25's obligation, and the tier-1/tier-2 inconsistency resolved.

    The planner decides once for the whole job, so the same term gets the same
    answer whichever tier asks — previously tier 1 refused an L0 term outright
    and tier 2 dropped it with a warning.
    """
    job = make_job(
        store,
        tmp_path,
        total=4,
        segment=4,
        terms=["usage.cai", "lm.naturalness"],
        capabilities={"usage.cai": "L1", "lm.naturalness": "L0"},
        radii={"usage.cai": 0, "lm.naturalness": 0},
    )
    assert job.plan.dropped == ["lm.naturalness"]
    assert job.plan.scored == ["usage.cai"]

    job.run()
    status = job.status()

    dropped = [w for w in status["warnings"] if w["code"] == "CF-SRCH-010"]
    assert len(dropped) == 1
    assert dropped[0]["stage"] == "tier_plan"
    assert dropped[0]["subject"] == "lm.naturalness"
    assert [o["term_id"] for o in status["resolved_objectives"]] == ["usage.cai"]


def test_p1_29_every_warning_is_attributable_to_a_stage_and_subject(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """DEP-ORCH-001: all three sources in one document, still distinguishable."""
    job = make_job(
        store,
        tmp_path,
        total=4,
        segment=4,
        terms=["usage.cai", "lm.naturalness"],
        capabilities={"usage.cai": "L1", "lm.naturalness": "L0"},
        radii={"usage.cai": 0, "lm.naturalness": 0},
    ).run()
    status = job.status()

    # Two operators and one dropped term. The two CF-ORCH-006 tier placeholders
    # are gone: P1-30 assembled tier 1 and tier 2, and TST-CLI-003 asserts that
    # disappearance directly.
    assert len(status["warnings"]) == 3
    for warning in status["warnings"]:
        assert warning["stage"], "a warning with no stage cannot be attributed"
        assert warning["subject"], "a warning with no subject names nothing"
        assert warning["code"].startswith("CF-")

    assert {w["stage"] for w in status["warnings"]} == {"operators", "tier_plan"}
    assert len(job.warnings_from("operators")) == 2
    assert len(job.warnings_from("tier_plan")) == 1
    assert {w["subject"] for w in job.warnings_from("tier_plan")} == {"lm.naturalness"}

    # Repeated segments must not repeat the same warning.
    longer = make_job(store, tmp_path / "long", total=12, segment=2).run()
    assert len([w for w in longer.status()["warnings"] if w["code"] == "CF-SRCH-001"]) == 2


def test_p1_29_a_term_wider_than_the_block_fails_rather_than_warns(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """CF-SRCH-008 invalidates the partition, so it is not a warning.

    Stated explicitly because it is sometimes listed alongside the warning
    codes. A term wider than the block does not degrade the search; it makes
    tier 2's independence claim false, so it belongs in `problem`.
    """
    with pytest.raises(TierPlanError) as excinfo:
        make_job(
            store,
            tmp_path,
            terms=["usage.cai", "structure.mfe_global"],
            capabilities={"usage.cai": "L1", "structure.mfe_global": "L1"},
            radii={"usage.cai": 0, "structure.mfe_global": None},
            block_size=30,
        )
    problem = excinfo.value.to_problem_details()
    assert problem["code"] == "CF-SRCH-008"
    assert problem["stage"] == "tier_plan"
    assert problem["subject"] == "structure.mfe_global"
    assert "P1-26" in problem["detail"]


def test_p1_29_the_status_document_validates(store: ArtifactStore, tmp_path: Path) -> None:
    jsonschema = pytest.importorskip("jsonschema")
    from referencing import Registry, Resource

    schema_dir = Path(__file__).resolve().parents[2] / "schemas" / "v1"
    registry = Registry().with_resources(
        [
            (path.name, Resource.from_contents(json.loads(path.read_text())))
            for path in schema_dir.glob("*.schema.json")
        ]
    )
    schema = json.loads((schema_dir / "job_status.schema.json").read_text())

    job = make_job(
        store,
        tmp_path,
        total=4,
        segment=4,
        terms=["usage.cai", "lm.naturalness"],
        capabilities={"usage.cai": "L1", "lm.naturalness": "L0"},
        radii={"usage.cai": 0, "lm.naturalness": 0},
    ).run()
    jsonschema.Draft202012Validator(schema, registry=registry).validate(job.status())


# --------------------------------------------------------------------------- #
# resume against the non-blocking-recalibration guarantee
# --------------------------------------------------------------------------- #


def test_p1_29_a_resumed_job_does_not_observe_a_mid_flight_promotion(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """The lockfile guarantee has to survive a gap, not just a continuous run.

    A checkpointed job is exactly the case where a promotion has time to land.
    Resume replays the sealed digests rather than re-resolving, so the newer
    artifact is invisible to it — while a job that has not sealed sees it, which
    is what makes the first half a guarantee rather than an artefact of nothing
    having changed.
    """
    job = make_job(store, tmp_path, total=8, segment=4)
    sealed_digest = job.lockfile.entry("epistasis.dwell_deviation").digest
    job.run(on_segment=lambda j: setattr(j, "cancel_requested", True))
    assert job.state == "cancelled"

    # A recalibration lands while the job is checkpointed.
    recalibrated = store.publish(
        manifest_document(version="1.1.0", note="recalibrated from new wet-lab data")
    )
    assert recalibrated.digest != sealed_digest

    document = Job.load_checkpoint(job.checkpoint_path())
    resumed = make_job(store, tmp_path, total=8, segment=4).resume(document)
    resumed.run()

    resolved = resumed.resolved_models()
    assert resolved["epistasis.dwell_deviation"].digest == sealed_digest
    assert resolved["epistasis.dwell_deviation"].version == "1.0.0"
    assert resumed.status()["landscape_lockfile_digest"] == job.lockfile.digest


# --------------------------------------------------------------------------- #
# the RNG plan
# --------------------------------------------------------------------------- #


def test_p1_29_the_root_seed_comes_from_the_lockfile_and_nowhere_else(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """A job layer with its own seed source makes the lockfile a fiction."""
    job = make_job(store, tmp_path, total=4, segment=4, root_seed=4242)
    assert job.root_seed() == 4242
    assert job.rng_plan()["algorithm"] == "xoshiro256pp"

    # Two jobs differing only in root seed must diverge; one repeated must not.
    a = make_job(store, tmp_path / "a", total=6, segment=6, root_seed=1).run()
    b = make_job(store, tmp_path / "b", total=6, segment=6, root_seed=2).run()
    again = make_job(store, tmp_path / "c", total=6, segment=6, root_seed=1).run()

    assert a.checkpoint == again.checkpoint
    assert a.checkpoint["history"] != b.checkpoint["history"], (
        "two root seeds produced identical runs; the seed is not reaching the engine"
    )


def test_p1_29_the_budget_is_enforced(store: ArtifactStore, tmp_path: Path) -> None:
    job = make_job(store, tmp_path, total=7, segment=3).run()
    assert job.generations_run == 7, "the last segment must be trimmed to the remainder"
    assert job.budget.remaining == 0
    assert job.budget.exhausted
    assert job.checkpoint["generation"] == 7


def test_dep_det_003_every_stage_derives_from_a_distinct_reserved_root(
    store: ArtifactStore, tmp_path: Path
) -> None:
    """Each stage's stream family must not overlap another's.

    `DEP-DET-002` requires a reserved index per stage. Tiers 1 and 2 were
    namespaced and the generation loop was left deriving from the job root
    directly, keying generations `0, 1, 2, …` — which is the index space a
    population loop's per-individual streams naturally use, so `P2-02` would
    have collided with it.

    Asserted on the roots rather than on a trajectory, because a trajectory
    test passes whichever root is used as long as it is used consistently:
    reverting to the job seed changes the answer but keeps every split-versus-
    whole comparison green. Found by a negative control that would not fire.
    """
    job = make_job(store, tmp_path, total=4, segment=4)
    root = job.root_seed()

    # Asserted on what the job *uses*, not merely on the derivation function:
    # a test that only called `_core.tier3_root_seed` would pass whichever
    # root the job actually passed to the engine.
    tier3 = job.tier3_seed()
    assert tier3 == _core.tier3_root_seed(root)
    assert tier3 != root, (
        "tier 3 derives from the job root itself, so its per-individual streams "
        "share an index space with the generation loop"
    )

    # And it is a function of the seed, not a constant.
    assert _core.tier3_root_seed(root + 1) != tier3

    # The job must pass it to the engine, which is the part the control
    # reverted. Checked by observing the seed the segment call receives.
    seen: list[int] = []
    original = _core.search_run_segment

    def watch(*args):
        seen.append(args[8])
        return original(*args)

    _core.search_run_segment = watch  # type: ignore[assignment]
    try:
        make_job(store, tmp_path / "watch", total=4, segment=4).run()
    finally:
        _core.search_run_segment = original  # type: ignore[assignment]

    assert seen, "the generation loop was never called"
    assert all(value == tier3 for value in seen), (
        f"the job passed {seen[0]} to the generation loop, not the reserved tier-3 root {tier3}"
    )
