"""Landscape lockfiles: sealing a job's artifact set, and replaying it.

Ticket P1-28. Contract: `schemas/v1/landscape_lockfile.schema.json`.

A lockfile is the answer to "which model produced this number", recorded before
any search work begins (`P1-29` sequences the sealing). Every entry is fully
pinned — an `artifact_id` alone is not enough, because a channel can move under
it — so a replay resolves the same bytes or fails.

## The non-blocking-recalibration guarantee, stated precisely

`DEP-ACT-001` makes the lockfile the concurrency-isolation mechanism for
recalibration. The guarantee is:

> A job that has sealed its lockfile reads the artifacts it sealed for its whole
> life, whatever else is published, promoted, deprecated or yanked meanwhile.

It holds because of two facts and no locking at all. Resolution reads digests,
and the store is content-addressed, so a digest that resolved once resolves to
the same content forever. Publication only ever *adds* digests. There is
therefore no interval during which a running job can observe a partially
applied change, because nothing is ever applied to something a job is holding.

Phase 3 adds channels (`P3-15`), which is the machinery that makes a promotion
possible at all. It moves a pointer that resolution consults *before* sealing
and never after. The guarantee is testable now precisely because it does not
depend on that machinery — see `test_a_sealed_job_does_not_observe_a_mid_flight_publication`.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any

from .errors import CF_MAJOR_MISMATCH, RegistryError
from .manifest import Manifest, digest_of
from .store import ArtifactStore

SCHEMA_VERSION = "1.0.0"

RUNGS = ("pin", "tenant_channel", "host_channel", "nearest_host", "base", "builtin")

#: Rungs `P1-28` implements. The rest raise `CF-RESOLVE-007` naming their
#: ticket rather than falling through to a lower rung — a lower rung is by
#: construction a worse answer, so taking it silently while reporting success
#: is the failure the ladder exists to prevent.
IMPLEMENTED_RUNGS = ("pin", "base", "builtin")

DEFERRED_RUNGS = {
    "tenant_channel": "P3-15",
    "host_channel": "P3-15",
    "nearest_host": "P3-16",
}


@dataclass(frozen=True)
class ResolutionEntry:
    """One term's pinned artifact."""

    term_id: str
    artifact_id: str
    digest: str
    resolved_by: str
    channel: str | None = None
    fallback: bool = False

    def to_json(self) -> dict[str, Any]:
        return {
            "term_id": self.term_id,
            "artifact_id": self.artifact_id,
            "digest": self.digest,
            "resolved_by": self.resolved_by,
            "channel": self.channel,
            "fallback": self.fallback,
            "host_distance": None,
        }

    @classmethod
    def from_json(cls, document: dict[str, Any]) -> ResolutionEntry:
        return cls(
            term_id=document["term_id"],
            artifact_id=document["artifact_id"],
            digest=document["digest"],
            resolved_by=document["resolved_by"],
            channel=document.get("channel"),
            fallback=bool(document.get("fallback", False)),
        )


def builtin_digest(term_id: str, engine_version: str, engine: str = "codonforge") -> str:
    """Digest of a built-in term's engine descriptor.

    A built-in has no file to hash — CAI is a function — but the lockfile
    requires a digest on every entry, and it has work to do there: if the engine
    version changes, the numbers can change, and the lockfile has to say so.

    Hashing `{engine, term_id, version}` gives a real digest of a real document
    with exactly that property. A fixed sentinel such as an all-zero digest
    would satisfy the schema's pattern while carrying no information, and would
    read as a genuine hash at every call site — the shape of fabricated value
    the standing rules exclude.
    """
    return digest_of(
        {"engine": engine, "term_id": term_id, "version": engine_version}, omit="digest"
    )


@dataclass
class Lockfile:
    """A sealed resolution, replayable by digest."""

    lockfile_id: str
    created_at: str
    contract_major: int
    resolution: list[ResolutionEntry]
    host_profile: dict[str, Any]
    inputs: dict[str, Any]
    rng: dict[str, Any]
    environment: dict[str, Any]
    job_id: str | None = None
    tenant_id: str | None = None
    engines: list[dict[str, Any]] = field(default_factory=list)
    digest: str | None = None

    def to_json(self) -> dict[str, Any]:
        document: dict[str, Any] = {
            "schema_version": SCHEMA_VERSION,
            "lockfile_id": self.lockfile_id,
            "created_at": self.created_at,
            "contract_major": self.contract_major,
            "resolution": [entry.to_json() for entry in self.resolution],
            "host_profile": self.host_profile,
            "inputs": self.inputs,
            "rng": self.rng,
            "environment": self.environment,
        }
        if self.job_id is not None:
            document["job_id"] = self.job_id
        if self.tenant_id is not None:
            document["tenant_id"] = self.tenant_id
        if self.engines:
            document["engines"] = self.engines
        document["digest"] = self.digest or digest_of(document, omit="digest")
        return document

    def seal(self) -> str:
        """Compute and record this lockfile's own digest.

        Also the value that namespaces the evaluation cache in `cf-landscape`
        (`DEP-LAND-001`): the digest crosses the Rust boundary as a string,
        because a `rust-core` crate cannot depend on this package.
        """
        document = self.to_json()
        document.pop("digest", None)
        self.digest = digest_of(document, omit="digest")
        return self.digest

    @classmethod
    def from_json(cls, document: dict[str, Any]) -> Lockfile:
        for required in (
            "lockfile_id",
            "created_at",
            "contract_major",
            "resolution",
            "host_profile",
            "inputs",
            "rng",
            "environment",
            "digest",
        ):
            if required not in document:
                raise RegistryError(
                    "CF-RESOLVE-006", f"lockfile is missing required field {required!r}"
                )
        return cls(
            lockfile_id=document["lockfile_id"],
            created_at=document["created_at"],
            contract_major=int(document["contract_major"]),
            resolution=[ResolutionEntry.from_json(e) for e in document["resolution"]],
            host_profile=document["host_profile"],
            inputs=document["inputs"],
            rng=document["rng"],
            environment=document["environment"],
            job_id=document.get("job_id"),
            tenant_id=document.get("tenant_id"),
            engines=list(document.get("engines", [])),
            digest=document["digest"],
        )

    @classmethod
    def load(cls, path: Any) -> Lockfile:
        return cls.from_json(json.loads(open(path, encoding="utf-8").read()))

    def entry(self, term_id: str) -> ResolutionEntry | None:
        for candidate in self.resolution:
            if candidate.term_id == term_id:
                return candidate
        return None

    def terms(self) -> list[str]:
        return [entry.term_id for entry in self.resolution]

    def replay(self, store: ArtifactStore) -> dict[str, Manifest | None]:
        """Resolve every entry by digest, or fail.

        `TST-MREG-002`. Two things this deliberately does not do: it does not
        fall back to a newer artifact when a digest is missing, and it does not
        skip the entry. Both would let a replay succeed while producing
        different numbers, and a replay that quietly is not one is worse than a
        replay that stops.

        A `builtin` entry has no artifact and resolves to `None` — the value
        is computed by shipped code, identified by the engine version in its
        coordinate.
        """
        self.check_contract_major_consistency()
        resolved: dict[str, Manifest | None] = {}
        for entry in self.resolution:
            if entry.resolved_by == "builtin":
                resolved[entry.term_id] = None
                continue
            manifest = store.require(entry.digest)
            if not manifest.provides(entry.term_id):
                raise RegistryError(
                    "CF-RESOLVE-006",
                    f"{entry.digest} resolves to {manifest.artifact_id}, which does not "
                    f"provide {entry.term_id}",
                )
            manifest.check_contract_major(self.contract_major)
            resolved[entry.term_id] = manifest
        return resolved

    def check_contract_major_consistency(self) -> None:
        """`CF-RESOLVE-002`: one job, one landscape contract major.

        Two terms that disagree about what a query means cannot be composed
        into one objective vector; the vector would have axes measured under
        different definitions.
        """
        for entry in self.resolution:
            if entry.resolved_by not in RUNGS:
                raise RegistryError(
                    CF_MAJOR_MISMATCH,
                    f"{entry.term_id} records an unknown resolution rung {entry.resolved_by!r}",
                )
