"""cf-modelreg — model artifact registry.

Ticket P1-28. Contracts: `schemas/v1/model_manifest.schema.json` and
`schemas/v1/landscape_lockfile.schema.json`.

Phase 1 is a deliberate subset of the scheme in ADR-0009 and
`planning/04-api/04-model-resolution.md`, not a parallel mechanism. It owns:

  * manifest identity — coordinate, content digest, `contract_major`,
    `provides_terms`;
  * a content-addressed local artifact store;
  * the `pin`, `base` and `builtin` rungs of the resolution ladder;
  * lockfile sealing and replay, including `CF-RESOLVE-006` on a missing digest;
  * the non-blocking-recalibration guarantee, which is a property of the
    lockfile plus content addressing and therefore holds before channels exist.

It defers channels and the `tenant_channel` / `host_channel` rungs to `P3-15`,
and `nearest_host` with its distance metric, admissibility rules and
uncertainty inflation to `P3-16`. A request that would need a deferred rung
raises `CF-RESOLVE-007` naming the rung and its ticket rather than falling
through to a lower one.

Declared dependencies: cf-types.
"""

from .errors import RegistryError
from .lockfile import Lockfile, ResolutionEntry, builtin_digest
from .manifest import Manifest, digest_of
from .resolve import (
    BUILTIN_TERMS,
    Pin,
    ResolutionRequest,
    implemented_rungs,
    resolve,
    seal,
    unimplemented_rungs,
)
from .store import ArtifactStore

__all__ = [
    "BUILTIN_TERMS",
    "ArtifactStore",
    "Lockfile",
    "Manifest",
    "Pin",
    "RegistryError",
    "ResolutionEntry",
    "ResolutionRequest",
    "builtin_digest",
    "digest_of",
    "implemented_rungs",
    "resolve",
    "seal",
    "unimplemented_rungs",
]
