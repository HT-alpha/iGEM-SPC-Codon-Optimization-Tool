"""The resolution ladder, Phase 1 subset.

Ticket P1-28. Specification: `planning/04-api/04-model-resolution.md`, ADR-0009.

For each requested term the first rung yielding a compatible artifact wins:

    1  pin              explicit model_pins, or a replay lockfile   IMPLEMENTED
    2  tenant_channel   the tenant's recalibrated adapter           P3-15
    3  host_channel     the normal path for a supported host        P3-15
    4  nearest_host     a related host, within a distance bound     P3-16
    5  base             the host-agnostic artifact                  IMPLEMENTED
    -  builtin          shipped code; not an artifact               IMPLEMENTED

## Why the missing rungs raise instead of falling through

A partially-implemented ladder that let unimplemented rungs fall through would
still return an artifact for most requests, and the artifact would be
plausible. That is the problem. The rungs are *ordered*: a lower rung is by
construction a worse answer, so silently taking one while reporting success
produces exactly the result nobody can explain later — "why did my E. coli job
use the host-agnostic model?" with nothing in the output that says it did.

So `resolve` refuses a request that would need a deferred rung, with
`CF-RESOLVE-007` naming the rung and the ticket that will implement it. A
caller who genuinely wants the host-agnostic model can ask for it: that is what
`allow_base` means.
"""

from __future__ import annotations

import datetime as _dt
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from .errors import CF_NO_ARTIFACT, CF_RUNG_DEFERRED, RegistryError
from .lockfile import (
    DEFERRED_RUNGS,
    IMPLEMENTED_RUNGS,
    Lockfile,
    ResolutionEntry,
    builtin_digest,
)
from .manifest import Manifest
from .store import ArtifactStore

#: Terms computed by shipped code rather than by a model artifact.
#:
#: Exactly the `P1-21` built-ins. A term is on this list because there is
#: nothing to train, not because no artifact happens to exist yet — the
#: difference matters, since the second case must raise `CF-RESOLVE-001`.
BUILTIN_TERMS = (
    "usage.cai",
    "usage.gc_content",
    "constraint.violation_count",
    "structure.local_accessibility",
    "burden.growth_cost",
)


@dataclass(frozen=True)
class Pin:
    """An explicit artifact pin: the top rung."""

    artifact_id: str
    digest: str


@dataclass(frozen=True)
class ResolutionRequest:
    """What a job wants resolved."""

    terms: Sequence[str]
    host_id: str
    contract_major: int
    tenant_id: str | None = None
    pins: dict[str, Pin] | None = None
    #: Permit the host-agnostic rung. Off by default: a caller who has not said
    #: that a host-agnostic model is acceptable for their host has not said it.
    allow_base: bool = True
    builtin_engine_version: str = "0.1.0"


def resolve(store: ArtifactStore, request: ResolutionRequest) -> list[ResolutionEntry]:
    """Resolve every requested term, in the caller's order.

    Fails on the first term that cannot be resolved, naming it. It does not
    return the terms it could resolve: a job that optimizes four objectives
    when the user asked for five, without saying so, is worse than a job that
    fails.
    """
    pins = request.pins or {}
    entries: list[ResolutionEntry] = []

    for term_id in request.terms:
        pin = pins.get(term_id)
        if pin is not None:
            manifest = store.require(pin.digest)
            _check_pin(manifest, pin, term_id, request)
            entries.append(
                ResolutionEntry(
                    term_id=term_id,
                    artifact_id=manifest.artifact_id,
                    digest=manifest.digest,
                    resolved_by="pin",
                )
            )
            continue

        if term_id in BUILTIN_TERMS:
            entries.append(
                ResolutionEntry(
                    term_id=term_id,
                    artifact_id=f"builtin:{term_id}:{request.builtin_engine_version}",
                    digest=builtin_digest(term_id, request.builtin_engine_version),
                    resolved_by="builtin",
                )
            )
            continue

        # Rungs 2-4 would be consulted here. They are not implemented, and a
        # term that needs one must not quietly drop to `base`.
        host_specific = store.find(
            term_id=term_id,
            contract_major=request.contract_major,
            host_id=request.host_id,
            tenant_id=request.tenant_id,
        )
        if host_specific:
            raise RegistryError(
                CF_RUNG_DEFERRED,
                f"{term_id} has a host-specific artifact for {request.host_id} "
                f"({host_specific[-1].artifact_id}), but selecting it needs the "
                f"host_channel rung, which arrives in {DEFERRED_RUNGS['host_channel']}. "
                "Pin the artifact by digest to use it now",
                term_id=term_id,
                rung="host_channel",
                ticket=DEFERRED_RUNGS["host_channel"],
                candidate=host_specific[-1].artifact_id,
            )

        if request.allow_base:
            base = store.find(
                term_id=term_id,
                contract_major=request.contract_major,
                host_agnostic=True,
                tenant_id=request.tenant_id,
            )
            if base:
                manifest = base[-1]
                entries.append(
                    ResolutionEntry(
                        term_id=term_id,
                        artifact_id=manifest.artifact_id,
                        digest=manifest.digest,
                        resolved_by="base",
                    )
                )
                continue

        raise RegistryError(
            CF_NO_ARTIFACT,
            f"no artifact resolves for term {term_id} and host {request.host_id}",
            term_id=term_id,
            host_id=request.host_id,
        )

    return entries


def _check_pin(manifest: Manifest, pin: Pin, term_id: str, request: ResolutionRequest) -> None:
    if manifest.artifact_id != pin.artifact_id:
        raise RegistryError(
            "CF-RESOLVE-006",
            f"pin for {term_id} names {pin.artifact_id} but digest {pin.digest} "
            f"holds {manifest.artifact_id}",
        )
    if not manifest.provides(term_id):
        raise RegistryError(
            CF_NO_ARTIFACT,
            f"pinned artifact {manifest.artifact_id} does not provide {term_id}",
        )
    manifest.check_contract_major(request.contract_major)
    if manifest.isolation_scope == "tenant" and manifest.tenant_id != request.tenant_id:
        # ADR-0006: a lockfile pinning a tenant-A artifact cannot be replayed by
        # tenant B. Reproducibility does not silently degrade into cross-tenant
        # access.
        raise RegistryError(
            "CF-RESOLVE-006",
            f"pinned artifact {manifest.artifact_id} belongs to another tenant",
        )
    # A yanked artifact still resolves by digest, forever, so that a published
    # result stays reproducible. Only channel-based resolution refuses it, and
    # channels are P3-15.


def unimplemented_rungs() -> dict[str, str]:
    """Deferred rungs and the ticket that will provide each.

    Public so a caller — the CLI's `hosts list`, the API's `/v1/models` — can
    say what is not available yet instead of appearing to offer it.
    """
    return dict(DEFERRED_RUNGS)


def implemented_rungs() -> tuple[str, ...]:
    return IMPLEMENTED_RUNGS


def seal(
    entries: Iterable[ResolutionEntry],
    *,
    lockfile_id: str,
    contract_major: int,
    host_profile: dict,
    inputs: dict,
    rng: dict,
    environment: dict,
    job_id: str | None = None,
    tenant_id: str | None = None,
    engines: list[dict] | None = None,
    created_at: str | None = None,
) -> Lockfile:
    """Build and seal a lockfile from a resolution.

    Sealing happens once, before any search work (`P1-29`). After it, the job's
    artifact set is fixed: this is the moment the non-blocking-recalibration
    guarantee starts applying.
    """
    lock = Lockfile(
        lockfile_id=lockfile_id,
        created_at=created_at or _dt.datetime.now(_dt.UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        contract_major=contract_major,
        resolution=list(entries),
        host_profile=host_profile,
        inputs=inputs,
        rng=rng,
        environment=environment,
        job_id=job_id,
        tenant_id=tenant_id,
        engines=list(engines or []),
    )
    lock.seal()
    return lock
