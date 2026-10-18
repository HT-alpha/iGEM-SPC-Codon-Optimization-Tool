"""A content-addressed artifact store on the local filesystem.

Ticket P1-28.

## Why content-addressed, and why that is the whole non-blocking guarantee

Artifacts are filed under their content digest, and a digest that is present
never changes meaning. That single property is what makes recalibration
non-blocking (`DEP-ACT-001`): a running job holds digests, a recalibration
publishes a *new* digest, and the two cannot collide. Nothing has to be locked,
no job has to be drained, and there is no window during which a job might see
half of an update — because there is no update, only an addition.

Publishing the same digest twice is idempotent. Publishing *different* content
under a digest it does not hash to is refused, which is not a hypothetical:
that is what a corrupted download looks like.

The Phase 3 channel machinery (`P3-15`) will sit on top of this. It moves
pointers; it does not mutate artifacts. That is why the guarantee this module
establishes survives the feature that appears to threaten it.
"""

from __future__ import annotations

import json
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from .errors import RegistryError
from .manifest import Manifest, digest_of


class ArtifactStore:
    """Manifests filed by content digest under a local directory."""

    def __init__(self, root: Path | str) -> None:
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def _path(self, digest: str) -> Path:
        if not digest.startswith("sha256:"):
            raise RegistryError("CF-RESOLVE-006", f"{digest!r} is not a sha256 digest reference")
        body = digest.removeprefix("sha256:")
        # Two-character shard, as git and every container registry do: a flat
        # directory of a few hundred thousand entries is slow to list on every
        # filesystem that matters.
        return self.root / body[:2] / f"{body}.json"

    def publish(self, document: dict[str, Any]) -> Manifest:
        """Store a manifest, computing and stamping its digest.

        The `digest` field is *derived*, never accepted from the caller. A
        store that trusted a submitted digest would let a mislabelled artifact
        in, and every later verification would agree with the label rather than
        the content.
        """
        stamped = dict(document)
        stamped.pop("digest", None)
        digest = digest_of(stamped, omit="digest")
        stamped["digest"] = digest

        manifest = Manifest.from_json(stamped)
        path = self._path(digest)
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists():
            existing = json.loads(path.read_text(encoding="utf-8"))
            if existing != stamped:
                raise RegistryError(
                    "CF-RESOLVE-006",
                    f"digest {digest} already holds different content; "
                    "a content-addressed store cannot have two answers",
                )
            return manifest
        path.write_text(json.dumps(stamped, indent=2, sort_keys=True), encoding="utf-8")
        return manifest

    def get(self, digest: str) -> Manifest | None:
        path = self._path(digest)
        if not path.is_file():
            return None
        manifest = Manifest.from_json(json.loads(path.read_text(encoding="utf-8")))
        manifest.verify_digest()
        return manifest

    def require(self, digest: str) -> Manifest:
        """Fetch by digest, or fail the way a replay must fail.

        `TST-MREG-002`: never substitutes a newer artifact. A silent
        substitution would turn "reproduce my result" into "produce a different
        result without telling me", which is worse than not reproducing it at
        all — the second is visible.
        """
        manifest = self.get(digest)
        if manifest is None:
            raise RegistryError(
                "CF-RESOLVE-006",
                f"pinned artifact digest {digest} is not in the store; replay cannot proceed",
                digest=digest,
            )
        return manifest

    def __contains__(self, digest: object) -> bool:
        return isinstance(digest, str) and self._path(digest).is_file()

    def manifests(self) -> Iterator[Manifest]:
        """Every stored manifest, in digest order.

        Sorted so a listing is stable across runs and filesystems; directory
        order is neither.
        """
        for path in sorted(self.root.rglob("*.json")):
            yield Manifest.from_json(json.loads(path.read_text(encoding="utf-8")))

    def find(
        self,
        *,
        term_id: str,
        contract_major: int,
        host_id: str | None = None,
        host_agnostic: bool = False,
        tenant_id: str | None = None,
    ) -> list[Manifest]:
        """Candidate artifacts for a term, newest version last.

        Filters on lifecycle and isolation the same way the ADR-0009
        compatibility check does: `yanked` never resolves by search (it stays
        retrievable by digest forever, so a published result remains
        reproducible), and a tenant artifact is invisible to everyone else.
        """
        found = []
        for manifest in self.manifests():
            if not manifest.provides(term_id):
                continue
            if manifest.contract_major != contract_major:
                continue
            if manifest.lifecycle_status in ("draft", "yanked"):
                continue
            if manifest.isolation_scope == "tenant" and manifest.tenant_id != tenant_id:
                continue
            if host_agnostic:
                if not manifest.is_host_agnostic():
                    continue
            elif host_id is not None and manifest.host_id != host_id:
                continue
            found.append(manifest)
        found.sort(key=lambda m: (_version_key(m.version), m.artifact_id))
        return found


def _version_key(version: str) -> tuple[int, ...]:
    """Order SemVer numerically, not lexicographically.

    `1.10.0` follows `1.9.0`; string order puts it before. Pre-release and
    build metadata are dropped for ordering purposes, which is enough here
    because resolution by version only ever picks the newest compatible
    artifact and a pre-release is not published to a channel.
    """
    core = version.split("-", 1)[0].split("+", 1)[0]
    parts = []
    for piece in core.split("."):
        try:
            parts.append(int(piece))
        except ValueError:
            parts.append(0)
    return tuple(parts)
