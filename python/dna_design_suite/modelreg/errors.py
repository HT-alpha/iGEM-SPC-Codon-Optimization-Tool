"""Registry errors, carrying the `CF-RESOLVE-*` codes.

Ticket P1-28.
"""

from __future__ import annotations

from typing import Any

CF_NO_ARTIFACT = "CF-RESOLVE-001"
CF_MAJOR_MISMATCH = "CF-RESOLVE-002"
CF_FALLBACK_USED = "CF-RESOLVE-004"
CF_DEPRECATED = "CF-RESOLVE-005"
CF_DIGEST_MISSING = "CF-RESOLVE-006"
CF_RUNG_DEFERRED = "CF-RESOLVE-007"

_STATUS = {
    CF_NO_ARTIFACT: 422,
    CF_MAJOR_MISMATCH: 422,
    CF_DIGEST_MISSING: 409,
    CF_RUNG_DEFERRED: 501,
}


class RegistryError(Exception):
    """A registry failure that carries a machine-readable `CF-*` code."""

    def __init__(self, code: str, detail: str, **extra: Any) -> None:
        self.code = code
        self.detail = detail
        self.status = _STATUS.get(code, 422)
        self.extra = extra
        super().__init__(f"{code}: {detail}")

    def to_problem_details(self) -> dict[str, Any]:
        return {
            "type": f"https://dna_design_suite.org/problems/{self.code.lower()}",
            "status": self.status,
            "code": self.code,
            "detail": self.detail,
            **self.extra,
        }
