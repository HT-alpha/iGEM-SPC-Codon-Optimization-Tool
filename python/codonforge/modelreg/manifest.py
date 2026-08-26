"""Model artifact manifests.

Ticket P1-28. Contract: `schemas/v1/model_manifest.schema.json`.

ADR-0009: the human-readable coordinate is not enough for reproducibility, so a
content digest accompanies it everywhere the coordinate appears. `DEP-VER-001`
records that as a departure from the architecture proposal, which treated the
coordinate as sufficient.

The digest is computed over the manifest with its own `digest` field omitted —
a document cannot contain its own hash — using the same RFC 8785 canonical form
the Rust core uses. `verify_digest` recomputes rather than trusting the field,
because a manifest whose recorded digest disagrees with its content is precisely
the object a replay must refuse.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any

from codonforge import _core

from .errors import CF_MAJOR_MISMATCH, RegistryError

#: Artifact kinds the coordinate grammar admits.
KINDS = (
    "epistasis_potts",
    "epistasis_dwell",
    "burden_capacity",
    "lm_naturalness",
    "structure_params",
    "kinetics_domainmap",
    "calibrator",
    "host_profile",
)

LIFECYCLE_STATUSES = ("draft", "active", "deprecated", "yanked")


def digest_of(document: dict[str, Any], *, omit: str = "digest") -> str:
    """Content digest of a document, with one field omitted.

    Delegates to the compiled core so the Python and Rust sides cannot disagree
    about what "the digest of this document" means. A second implementation
    here would be a second definition, and the two would drift on the first
    unusual float.
    """
    return _core.digest_json(json.dumps(document), omit)


@dataclass(frozen=True)
class Manifest:
    """One model artifact's identity and provenance."""

    artifact_id: str
    kind: str
    host_id: str | None
    version: str
    contract_major: int
    provides_terms: tuple[str, ...]
    digest: str
    lifecycle_status: str = "active"
    isolation_scope: str = "public"
    tenant_id: str | None = None
    document: dict[str, Any] = field(default_factory=dict, repr=False)

    @classmethod
    def from_json(cls, document: dict[str, Any]) -> Manifest:
        for required in (
            "artifact_id",
            "kind",
            "version",
            "contract_major",
            "provides_terms",
            "digest",
        ):
            if required not in document:
                raise RegistryError(
                    "CF-RESOLVE-001", f"manifest is missing required field {required!r}"
                )
        if document["kind"] not in KINDS:
            raise RegistryError(
                "CF-RESOLVE-001",
                f"unknown artifact kind {document['kind']!r}; expected one of {', '.join(KINDS)}",
            )
        lifecycle = document.get("lifecycle", {})
        status = lifecycle.get("status", "active")
        if status not in LIFECYCLE_STATUSES:
            raise RegistryError("CF-RESOLVE-001", f"unknown lifecycle status {status!r}")
        isolation = document.get("isolation", {})

        return cls(
            artifact_id=document["artifact_id"],
            kind=document["kind"],
            host_id=document.get("host_id"),
            version=document["version"],
            contract_major=int(document["contract_major"]),
            provides_terms=tuple(document["provides_terms"]),
            digest=document["digest"],
            lifecycle_status=status,
            isolation_scope=isolation.get("scope", "public"),
            tenant_id=isolation.get("tenant_id"),
            document=document,
        )

    def verify_digest(self) -> None:
        """Recompute the digest and compare.

        Raised as `CF-RESOLVE-006` rather than a generic validation error: from
        a replay's point of view, a manifest whose content does not match the
        digest it claims is indistinguishable from the pinned artifact being
        unavailable, and both must stop the replay.
        """
        recomputed = digest_of(self.document)
        if recomputed != self.digest:
            raise RegistryError(
                "CF-RESOLVE-006",
                f"{self.artifact_id} records digest {self.digest} but its content "
                f"hashes to {recomputed}",
            )

    def provides(self, term_id: str) -> bool:
        return term_id in self.provides_terms

    def is_host_agnostic(self) -> bool:
        """Whether this artifact is the `base` rung for its kind.

        Host-agnostic means the host is an *input* to the model rather than
        baked into its weights — the LM naturalness model conditions on a host
        embedding, so one artifact serves every host. It is a property the
        manifest states, not one inferred from a missing field.
        """
        return self.host_id in (None, "host_agnostic")

    def check_contract_major(self, expected: int) -> None:
        if self.contract_major != expected:
            raise RegistryError(
                CF_MAJOR_MISMATCH,
                f"{self.artifact_id} speaks landscape contract major "
                f"{self.contract_major}; this job speaks {expected}",
            )
