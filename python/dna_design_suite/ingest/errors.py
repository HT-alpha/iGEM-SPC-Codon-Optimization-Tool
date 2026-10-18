"""Typed ingest errors carrying a `CF-*` code.

Tickets P1-09..P1-12.

WHY INGEST DOES NOT VALIDATE SEQUENCES ITSELF. `CodonSequence` already enforces
frame, alphabet, stop placement, translation and lock ranges (`P1-06`), and it
reports each with a specific `CF-SEQ-*` code. Re-implementing those checks here
would create a second, drifting definition of what a valid sequence is — and the
first symptom would be a file that ingest accepts and the core later rejects.

So ingest parses, hands the result to the core, and translates the core's
problem detail into an `IngestError` that keeps the code and adds where in the
file it came from. A malformed record therefore fails *at ingest*, with a code
the caller can branch on, rather than in some later ticket that assumed ingest
had already checked.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any


@dataclass
class SourceLocation:
    """Where in an input a problem was found."""

    uri: str | None = None
    record_id: str | None = None
    line: int | None = None

    def describe(self) -> str:
        parts = []
        if self.uri:
            parts.append(str(self.uri))
        if self.record_id:
            parts.append(f"record {self.record_id!r}")
        if self.line is not None:
            parts.append(f"line {self.line}")
        return ", ".join(parts) if parts else "input"


class IngestError(Exception):
    """An ingest failure that carries a machine-readable `CF-*` code."""

    def __init__(
        self,
        code: str,
        detail: str,
        *,
        location: SourceLocation | None = None,
        status: int = 422,
        extra: dict[str, Any] | None = None,
    ) -> None:
        self.code = code
        self.detail = detail
        self.location = location or SourceLocation()
        self.status = status
        self.extra = extra or {}
        super().__init__(f"{code}: {detail} (at {self.location.describe()})")

    @classmethod
    def from_core(cls, exc: Exception, *, location: SourceLocation | None = None) -> IngestError:
        """Translate a `dna_design_suite._core` error, preserving its code.

        The core raises `ValueError` whose message is the RFC 9457 problem
        detail. Parsing it back rather than re-deriving a code is what keeps the
        two sides from drifting.
        """
        try:
            problem = json.loads(str(exc))
            code = problem["code"]
            detail = problem.get("detail", "")
            status = problem.get("status", 422)
        except (ValueError, KeyError):
            # Not a problem detail: surface it rather than inventing a code.
            return cls("CF-SEQ-001", f"unparseable core error: {exc}", location=location)
        return cls(code, detail, location=location, status=status)

    def to_problem_details(self) -> dict[str, Any]:
        return {
            "type": f"https://dna_design_suite.org/problems/{self.code.lower()}",
            "title": self.detail.split(";")[0][:120],
            "status": self.status,
            "code": self.code,
            "detail": self.detail,
            "location": self.location.describe(),
            **self.extra,
        }


@dataclass
class IngestReport:
    """Non-fatal observations from a successful ingest.

    Anything that changed the meaning of the input without failing it — a
    transcoded RNA sequence, a record skipped, an ambiguity code encountered —
    lands here rather than in a log, because a user whose input was silently
    altered has no way to find out otherwise.
    """

    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def warn(self, message: str) -> None:
        self.warnings.append(message)

    def note(self, message: str) -> None:
        self.notes.append(message)

    def __bool__(self) -> bool:
        return bool(self.warnings or self.notes)


# Codes ingest may raise that are not produced by the core. Everything else is
# passed through from `CodonSequence`.
CF_INGEST_MALFORMED_RECORD = "CF-SEQ-001"
CF_INGEST_INVALID_BASE = "CF-SEQ-002"
CF_INGEST_BAD_IDENTIFIER = "CF-SEQ-013"
CF_HOST_PROFILE_INVALID = "CF-HOST-002"
