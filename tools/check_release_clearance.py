#!/usr/bin/env python3
"""Block a release of structure-kernel functionality without human sign-off.

`AMD-013`. Run by the `legal-clearance` job in `.github/workflows/wheels.yml`,
which every tagged release must pass before a wheel is built.

Why this is a script and not a line in a planning document
----------------------------------------------------------
`AMD-013` cleared ViennaRNA and LinearFold from their primary licence texts, and
that analysis is solid. It also concluded that CodonForge does not practise
LinearDesign's claimed method, and that conclusion is an engineer reading claim
text. Claim construction is a legal exercise, the doctrine of equivalents cannot
be evaluated by reading claims, and grant status for the patent family was never
confirmed. "We avoided the claimed method" only works as a defence if the
reasoning about what was avoided survives legal scrutiny rather than engineering
scrutiny.

A requirement recorded only in `planning/` is a requirement a release process
never consults. This one runs in the release path and fails it.

The gate is expected to fail today. Nobody has signed off. That is the correct
state, not a bug to route around: the fix is a signature, not a code change.
"""

from __future__ import annotations

import datetime as dt
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CLEARANCE = ROOT / "docs" / "legal-clearance.yaml"

# Presence of any of these means the build ships structure-kernel functionality
# and the gate applies. Keyed to what the code *does*, not to a version number
# someone would have to remember to bump.
STRUCTURE_SURFACE = [
    ROOT / "python" / "codonforge" / "predictors" / "structure.py",
    ROOT / "python" / "codonforge" / "predictors" / "viennarna.py",
]

REQUIRED_FIELDS = {
    "signed_off_by": "the name of the person who reviewed this",
    "role": "their role - faculty advisor, PI, or institutional tech-transfer/legal",
    "date": "when they signed off (YYYY-MM-DD)",
    "amendment": "the amendment reviewed; must be AMD-013",
    "scope": "what they actually reviewed, in their words",
}

# A sign-off is about a specific analysis. If the analysis changes materially,
# the old signature does not cover the new one.
MAX_AGE_DAYS = 365


def ships_structure_kernel() -> list[Path]:
    return [path for path in STRUCTURE_SURFACE if path.exists()]


def main() -> int:
    shipped = ships_structure_kernel()
    if not shipped:
        print("[clearance] no structure-kernel surface in this build; gate not applicable")
        return 0

    print("[clearance] this build ships structure-kernel functionality:")
    for path in shipped:
        print(f"             {path.relative_to(ROOT)}")

    if not CLEARANCE.exists():
        print()
        print("BLOCKED: no legal clearance on file.")
        print()
        print("  AMD-013 permits ViennaRNA and LinearFold from their own licence texts,")
        print("  and concludes CodonForge does not practise LinearDesign's claimed")
        print("  method. The first half is a reading of plain licence terms. The second")
        print("  is an engineer's claim-scope analysis, made without confirming grant")
        print("  status, and it is not a legal opinion.")
        print()
        print("  Before a release that includes structure-kernel functionality ships, a")
        print("  faculty advisor, PI, or institutional technology-transfer/legal contact")
        print("  must review AMD-013 and record their sign-off at:")
        print(f"      {CLEARANCE.relative_to(ROOT)}")
        print()
        print("  See docs/legal-clearance.example.yaml for the required fields.")
        print("  This is a hard gate. Do not remove this check to make a release pass.")
        return 1

    try:
        document = yaml.safe_load(CLEARANCE.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as exc:
        print(f"BLOCKED: {CLEARANCE.name} is not valid YAML: {exc}")
        return 1

    problems = []
    for field, description in REQUIRED_FIELDS.items():
        value = document.get(field)
        if value is None or (isinstance(value, str) and not value.strip()):
            problems.append(f"missing {field} ({description})")

    if document.get("amendment") not in (None, "AMD-013"):
        problems.append(
            f"amendment is {document['amendment']!r}; this gate covers AMD-013. A "
            "sign-off on a different analysis does not clear this one"
        )

    signed = document.get("date")
    if signed is not None:
        if isinstance(signed, str):
            try:
                signed = dt.date.fromisoformat(signed)
            except ValueError:
                problems.append(f"date {signed!r} is not YYYY-MM-DD")
                signed = None
        if isinstance(signed, dt.date):
            age = (dt.date.today() - signed).days
            if age > MAX_AGE_DAYS:
                problems.append(
                    f"the sign-off is {age} days old (limit {MAX_AGE_DAYS}); patent "
                    "status and claim scope change over time, and a stale signature "
                    "covers an analysis that may no longer describe this code"
                )
            if age < 0:
                problems.append(f"the sign-off is dated {age * -1} days in the future")

    if problems:
        print("BLOCKED: the clearance record is incomplete:")
        for problem in problems:
            print(f"  - {problem}")
        return 1

    print(
        f"[clearance] signed off by {document['signed_off_by']} "
        f"({document['role']}) on {document['date']} for {document['amendment']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
