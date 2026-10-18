"""TST-MREG-002 and the P1-28 scope boundary.

Ticket P1-28.

The headline test is the non-blocking-recalibration control: a job that sealed
its lockfile must not observe an artifact published after it sealed. That
guarantee belongs to the lockfile and to content addressing, so it is testable
now even though the machinery that makes a *promotion* possible — channels —
arrives in `P3-15`.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dna_design_suite.modelreg import (
    ArtifactStore,
    Lockfile,
    Pin,
    RegistryError,
    ResolutionRequest,
    implemented_rungs,
    resolve,
    seal,
    unimplemented_rungs,
)

SCHEMA_DIR = Path(__file__).resolve().parents[2] / "schemas" / "v1"

HOST = "escherichia_coli:511145/BL21"


def manifest_document(
    *,
    kind: str = "epistasis_dwell",
    host_id: str | None = HOST,
    version: str = "1.0.0",
    terms: tuple[str, ...] = ("epistasis.dwell_deviation",),
    contract_major: int = 1,
    status: str = "active",
    tenant_id: str | None = None,
    note: str = "",
) -> dict:
    """A schema-valid manifest.

    `note` exists so two otherwise identical artifacts hash differently, which
    is how the mid-flight-publication test creates a genuinely new digest for
    the same coordinate.
    """
    host_part = "host_agnostic" if host_id in (None, "host_agnostic") else host_id
    return {
        "schema_version": "1.0.0",
        "artifact_id": f"cf:{kind}:{host_part}:{version}",
        "kind": kind,
        "host_id": None if host_id in (None, "host_agnostic") else host_id,
        "version": version,
        "contract_major": contract_major,
        "provides_terms": list(terms),
        "created_at": "2026-08-26T00:00:00Z",
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
        "evaluation": {"metrics": {"spearman": 0.61 if not note else 0.72}, "gates": []},
        "isolation": {
            "scope": "tenant" if tenant_id else "public",
            "tenant_id": tenant_id,
            "may_export": True,
        },
        "lifecycle": {"status": status},
        "created_by": note or "ci",
    }


def lock_scaffold() -> dict:
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
            "algorithm": "pcg64",
            "root_seed": 20260826,
            "stream_assignment": "per_block_index",
            "deterministic_reduction": True,
        },
        "environment": {
            "codonforge_version": "0.1.0",
            "git_sha": "0" * 40,
            "platform": "linux-x86_64",
            "thread_count": 8,
            "budget_tier": "standard",
        },
    }


@pytest.fixture()
def store(tmp_path: Path) -> ArtifactStore:
    return ArtifactStore(tmp_path / "artifacts")


# --------------------------------------------------------------------------- #
# store and manifest
# --------------------------------------------------------------------------- #


def test_publish_derives_the_digest_and_is_idempotent(store: ArtifactStore) -> None:
    document = manifest_document()
    first = store.publish(document)
    second = store.publish(dict(document))
    assert first.digest == second.digest
    assert first.digest in store

    # A submitted digest is ignored, not trusted: a store that believed the
    # label would let a mislabelled artifact in and then agree with the label
    # at every later verification.
    lying = dict(document)
    lying["digest"] = "sha256:" + "f" * 64
    assert store.publish(lying).digest == first.digest


def test_a_stored_manifest_verifies_against_its_own_content(store: ArtifactStore) -> None:
    manifest = store.publish(manifest_document())
    manifest.verify_digest()

    # Tamper with the file on disk; the store must refuse it on read.
    path = next(Path(store.root).rglob("*.json"))
    document = json.loads(path.read_text())
    document["evaluation"]["metrics"]["spearman"] = 0.99
    path.write_text(json.dumps(document))

    with pytest.raises(RegistryError) as excinfo:
        store.get(manifest.digest)
    assert excinfo.value.code == "CF-RESOLVE-006"


def test_manifests_validate_against_the_shipped_schema(store: ArtifactStore) -> None:
    jsonschema = pytest.importorskip("jsonschema")
    from referencing import Registry, Resource

    registry = Registry().with_resources(
        [
            (path.name, Resource.from_contents(json.loads(path.read_text())))
            for path in SCHEMA_DIR.glob("*.schema.json")
        ]
    )
    schema = json.loads((SCHEMA_DIR / "model_manifest.schema.json").read_text())
    validator = jsonschema.Draft202012Validator(schema, registry=registry)

    manifest = store.publish(manifest_document())
    validator.validate(manifest.document)


# --------------------------------------------------------------------------- #
# resolution ladder
# --------------------------------------------------------------------------- #


def test_builtin_terms_resolve_without_an_artifact(store: ArtifactStore) -> None:
    entries = resolve(
        store,
        ResolutionRequest(
            terms=["usage.cai", "structure.local_accessibility", "structure.local_accessibility"],
            host_id=HOST,
            contract_major=1,
        ),
    )
    assert [e.resolved_by for e in entries] == ["builtin"] * 3
    assert all(e.artifact_id.startswith("builtin:") for e in entries)


def test_a_pin_resolves_by_digest_and_is_checked(store: ArtifactStore) -> None:
    manifest = store.publish(manifest_document())
    entries = resolve(
        store,
        ResolutionRequest(
            terms=["epistasis.dwell_deviation"],
            host_id=HOST,
            contract_major=1,
            pins={"epistasis.dwell_deviation": Pin(manifest.artifact_id, manifest.digest)},
        ),
    )
    assert entries[0].resolved_by == "pin"
    assert entries[0].digest == manifest.digest

    # A pin whose coordinate and digest disagree is refused rather than one of
    # them being preferred.
    with pytest.raises(RegistryError) as excinfo:
        resolve(
            store,
            ResolutionRequest(
                terms=["epistasis.dwell_deviation"],
                host_id=HOST,
                contract_major=1,
                pins={
                    "epistasis.dwell_deviation": Pin(
                        "cf:epistasis_dwell:escherichia_coli:511145/BL21:9.9.9",
                        manifest.digest,
                    )
                },
            ),
        )
    assert excinfo.value.code == "CF-RESOLVE-006"


def test_the_base_rung_resolves_a_host_agnostic_artifact(store: ArtifactStore) -> None:
    store.publish(
        manifest_document(
            kind="lm_naturalness",
            host_id=None,
            version="2.0.1",
            terms=("lm.naturalness",),
        )
    )
    entries = resolve(
        store,
        ResolutionRequest(terms=["lm.naturalness"], host_id=HOST, contract_major=1),
    )
    assert entries[0].resolved_by == "base"
    assert "host_agnostic" in entries[0].artifact_id


def test_base_picks_the_newest_version_numerically(store: ArtifactStore) -> None:
    for version in ("2.9.0", "2.10.0", "2.2.0"):
        store.publish(
            manifest_document(
                kind="lm_naturalness",
                host_id=None,
                version=version,
                terms=("lm.naturalness",),
            )
        )
    entries = resolve(
        store,
        ResolutionRequest(terms=["lm.naturalness"], host_id=HOST, contract_major=1),
    )
    assert entries[0].artifact_id.endswith(":2.10.0"), (
        "SemVer must order numerically; lexicographic order puts 2.10.0 before 2.9.0"
    )


def test_a_deferred_rung_raises_rather_than_falling_through(store: ArtifactStore) -> None:
    """A host-specific artifact exists but needs a rung P1-28 does not have.

    Falling through to `base` would return the host-agnostic model and report
    success, which is the substitution the ladder's ordering exists to prevent.
    """
    store.publish(
        manifest_document(
            kind="lm_naturalness", host_id=HOST, version="2.0.1", terms=("lm.naturalness",)
        )
    )
    store.publish(
        manifest_document(
            kind="lm_naturalness", host_id=None, version="2.0.1", terms=("lm.naturalness",)
        )
    )

    with pytest.raises(RegistryError) as excinfo:
        resolve(
            store,
            ResolutionRequest(terms=["lm.naturalness"], host_id=HOST, contract_major=1),
        )
    assert excinfo.value.code == "CF-RESOLVE-007"
    assert excinfo.value.extra["rung"] == "host_channel"
    assert excinfo.value.extra["ticket"] == "P3-15"

    # And the escape hatch the message names actually works.
    host_specific = store.find(term_id="lm.naturalness", contract_major=1, host_id=HOST)[0]
    entries = resolve(
        store,
        ResolutionRequest(
            terms=["lm.naturalness"],
            host_id=HOST,
            contract_major=1,
            pins={"lm.naturalness": Pin(host_specific.artifact_id, host_specific.digest)},
        ),
    )
    assert entries[0].resolved_by == "pin"


def test_the_scope_boundary_is_explicit() -> None:
    assert set(implemented_rungs()) == {"pin", "base", "builtin"}
    assert unimplemented_rungs() == {
        "tenant_channel": "P3-15",
        "host_channel": "P3-15",
        "nearest_host": "P3-16",
    }


def test_a_missing_term_fails_by_name(store: ArtifactStore) -> None:
    with pytest.raises(RegistryError) as excinfo:
        resolve(
            store,
            ResolutionRequest(terms=["kinetics.not_a_term"], host_id=HOST, contract_major=1),
        )
    assert excinfo.value.code == "CF-RESOLVE-001"
    assert "kinetics.not_a_term" in excinfo.value.detail


def test_contract_major_mismatch_is_refused(store: ArtifactStore) -> None:
    manifest = store.publish(manifest_document(contract_major=2))
    with pytest.raises(RegistryError) as excinfo:
        resolve(
            store,
            ResolutionRequest(
                terms=["epistasis.dwell_deviation"],
                host_id=HOST,
                contract_major=1,
                pins={"epistasis.dwell_deviation": Pin(manifest.artifact_id, manifest.digest)},
            ),
        )
    assert excinfo.value.code == "CF-RESOLVE-002"


def test_a_tenant_artifact_is_invisible_to_another_tenant(store: ArtifactStore) -> None:
    """ADR-0006: reproducibility does not degrade into cross-tenant access."""
    manifest = store.publish(manifest_document(tenant_id="t_alpha"))
    with pytest.raises(RegistryError) as excinfo:
        resolve(
            store,
            ResolutionRequest(
                terms=["epistasis.dwell_deviation"],
                host_id=HOST,
                contract_major=1,
                tenant_id="t_beta",
                pins={"epistasis.dwell_deviation": Pin(manifest.artifact_id, manifest.digest)},
            ),
        )
    assert excinfo.value.code == "CF-RESOLVE-006"


# --------------------------------------------------------------------------- #
# TST-MREG-002 — replay fidelity
# --------------------------------------------------------------------------- #


def _sealed(store: ArtifactStore, pins: dict[str, Pin] | None = None) -> Lockfile:
    entries = resolve(
        store,
        ResolutionRequest(
            terms=["epistasis.dwell_deviation", "usage.cai"],
            host_id=HOST,
            contract_major=1,
            pins=pins,
        ),
    )
    return seal(
        entries,
        lockfile_id="lock_01JBQ8Z4K7M2N3P5R6T7V8W9XZ",
        contract_major=1,
        **lock_scaffold(),
    )


def test_tst_mreg_002_a_missing_digest_fails_and_never_substitutes(
    store: ArtifactStore, tmp_path: Path
) -> None:
    manifest = store.publish(manifest_document(version="1.0.0"))
    lock = _sealed(
        store,
        {"epistasis.dwell_deviation": Pin(manifest.artifact_id, manifest.digest)},
    )
    assert lock.replay(store)["epistasis.dwell_deviation"].digest == manifest.digest

    # A newer artifact for the same coordinate appears, and the pinned one is
    # removed. The replay must stop rather than take the newer one.
    newer = store.publish(manifest_document(version="1.1.0", note="recalibrated"))
    assert newer.digest != manifest.digest
    next(
        p
        for p in Path(store.root).rglob("*.json")
        if json.loads(p.read_text())["digest"] == manifest.digest
    ).unlink()

    with pytest.raises(RegistryError) as excinfo:
        lock.replay(store)
    assert excinfo.value.code == "CF-RESOLVE-006"
    assert manifest.digest in excinfo.value.detail
    assert newer.digest not in excinfo.value.detail


def test_a_lockfile_round_trips_and_seals_deterministically(store: ArtifactStore) -> None:
    manifest = store.publish(manifest_document())
    lock = _sealed(
        store,
        {"epistasis.dwell_deviation": Pin(manifest.artifact_id, manifest.digest)},
    )
    document = lock.to_json()
    reloaded = Lockfile.from_json(document)
    assert reloaded.to_json() == document
    assert reloaded.digest == lock.digest
    assert reloaded.terms() == ["epistasis.dwell_deviation", "usage.cai"]

    jsonschema = pytest.importorskip("jsonschema")
    from referencing import Registry, Resource

    registry = Registry().with_resources(
        [
            (path.name, Resource.from_contents(json.loads(path.read_text())))
            for path in SCHEMA_DIR.glob("*.schema.json")
        ]
    )
    schema = json.loads((SCHEMA_DIR / "landscape_lockfile.schema.json").read_text())
    jsonschema.Draft202012Validator(schema, registry=registry).validate(document)


def test_a_builtin_entry_replays_without_an_artifact(store: ArtifactStore) -> None:
    manifest = store.publish(manifest_document())
    lock = _sealed(
        store,
        {"epistasis.dwell_deviation": Pin(manifest.artifact_id, manifest.digest)},
    )
    resolved = lock.replay(store)
    assert resolved["usage.cai"] is None
    assert resolved["epistasis.dwell_deviation"] is not None


# --------------------------------------------------------------------------- #
# non-blocking recalibration
# --------------------------------------------------------------------------- #


def test_a_sealed_job_does_not_observe_a_mid_flight_publication(
    store: ArtifactStore,
) -> None:
    """The `DEP-ACT-001` guarantee, exercised with the parts `P1-28` owns.

    A job seals a lockfile. A recalibration then publishes a new artifact for
    the same coordinate — the situation `P3-15`'s channel promotion will create
    for real. The job replays and must still see the digest it sealed, with no
    locking anywhere: the store only ever adds digests, and a digest that
    resolved once resolves to the same bytes forever.
    """
    original = store.publish(manifest_document(version="1.0.0"))
    lock = _sealed(store, {"epistasis.dwell_deviation": Pin(original.artifact_id, original.digest)})
    sealed_digest = lock.entry("epistasis.dwell_deviation").digest

    # ---- recalibration lands while the job is running ----
    recalibrated = store.publish(
        manifest_document(version="1.1.0", note="recalibrated from new wet-lab data")
    )
    assert recalibrated.digest != original.digest
    assert recalibrated.digest in store and original.digest in store

    for _ in range(3):
        resolved = lock.replay(store)
        assert resolved["epistasis.dwell_deviation"].digest == sealed_digest
        assert resolved["epistasis.dwell_deviation"].version == "1.0.0"

    # A job that has *not* sealed sees the new one, which is what makes the
    # assertion above a guarantee rather than an artefact of nothing changing.
    fresh = resolve(
        store,
        ResolutionRequest(
            terms=["epistasis.dwell_deviation"],
            host_id=HOST,
            contract_major=1,
            pins={"epistasis.dwell_deviation": Pin(recalibrated.artifact_id, recalibrated.digest)},
        ),
    )
    assert fresh[0].digest == recalibrated.digest
    assert fresh[0].digest != sealed_digest


def test_deprecating_and_yanking_do_not_disturb_a_sealed_job(
    store: ArtifactStore,
) -> None:
    """Lifecycle changes are published as new manifests, so a pin is untouched.

    `yanked` never resolves through a channel and always resolves by digest —
    that asymmetry is what keeps a published result reproducible after its
    model has been withdrawn.
    """
    original = store.publish(manifest_document(version="1.0.0"))
    lock = _sealed(store, {"epistasis.dwell_deviation": Pin(original.artifact_id, original.digest)})

    store.publish(manifest_document(version="1.0.0", status="yanked", note="withdrawn"))
    assert lock.replay(store)["epistasis.dwell_deviation"].digest == original.digest

    # The yanked manifest is excluded from search, and still retrievable by
    # digest.
    yanked = next(m for m in store.manifests() if m.lifecycle_status == "yanked")
    assert yanked not in store.find(
        term_id="epistasis.dwell_deviation", contract_major=1, host_id=HOST
    )
    assert store.require(yanked.digest).lifecycle_status == "yanked"
