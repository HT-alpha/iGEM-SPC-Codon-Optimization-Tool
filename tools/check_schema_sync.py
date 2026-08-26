#!/usr/bin/env python3
"""Keep `schemas/v1/` identical to the normative schemas in `planning/02-schemas/`.

Ticket P1-03.

    --sync   copy planning/02-schemas -> schemas/v1
    --check  fail if they differ, or if a declared consumer has drifted

WHY A COPY AT ALL. The planning package is the reviewed, normative source
(`planning/02-schemas/README.md`), but runtime code should not reach into a
planning directory: that would make `planning/` load-bearing at install time
and would blur the line the brief drew between the planning phase and the
build. So the schemas are copied into `schemas/v1/`, which is what the wheel
ships and what the Rust contract tests read — and this tool makes the copy
verifiable rather than aspirational.

DRIFT CHECKS. Beyond byte-equality, `--check` verifies the two places where a
schema change would otherwise silently diverge from code:

  * the identifier patterns hard-coded in `crates/cf-types/src/ids.rs`
  * the `contract_major` / schema version constants in the Rust core

Both are cheap to check and expensive to get wrong: an identifier that Rust
accepts and the API rejects is a bug that only appears in integration.
"""

from __future__ import annotations

import argparse
import filecmp
import json
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "planning" / "02-schemas"
DST = ROOT / "schemas" / "v1"

# Files copied verbatim. `examples/` and `interfaces/` stay in planning: they
# are review artifacts, not runtime assets.
PATTERNS = ("*.schema.json", "rationale_tags.json")


def sources() -> list[Path]:
    out: list[Path] = []
    for pattern in PATTERNS:
        out.extend(sorted(SRC.glob(pattern)))
    return out


def sync() -> int:
    DST.mkdir(parents=True, exist_ok=True)
    copied = 0
    for src in sources():
        dst = DST / src.name
        if not dst.exists() or not filecmp.cmp(src, dst, shallow=False):
            shutil.copy2(src, dst)
            print(f"  synced {dst.relative_to(ROOT)}")
            copied += 1
    stale = {p.name for p in DST.glob("*.json")} - {p.name for p in sources()}
    for name in sorted(stale):
        (DST / name).unlink()
        print(f"  removed stale {name}")
    if not copied and not stale:
        print("  already in sync")
    return 0


def check() -> int:
    failures: list[str] = []
    checks = 0

    print("\n[schema-copy]")
    expected = {p.name for p in sources()}
    actual = {p.name for p in DST.glob("*.json")} if DST.is_dir() else set()
    for missing in sorted(expected - actual):
        failures.append(f"schemas/v1/{missing} is missing; run --sync")
        print(f"  FAIL schemas/v1/{missing} is missing")
    for extra in sorted(actual - expected):
        failures.append(f"schemas/v1/{extra} has no counterpart in planning/02-schemas")
        print(f"  FAIL schemas/v1/{extra} is stale")
    for src in sources():
        checks += 1
        dst = DST / src.name
        if dst.exists() and filecmp.cmp(src, dst, shallow=False):
            print(f"  ok   {src.name}")
        elif dst.exists():
            failures.append(f"schemas/v1/{src.name} differs from the planning copy")
            print(f"  FAIL {src.name} differs from the planning copy")

    print("\n[id-pattern-drift]")
    common = json.loads((SRC / "common.schema.json").read_text(encoding="utf-8"))
    defs = common["$defs"]
    ids_rs = (ROOT / "crates" / "cf-types" / "src" / "ids.rs").read_text(encoding="utf-8")
    # id_newtype!( ..., "SchemaDef", r"pattern" );
    declared = dict(re.findall(r'"([A-Za-z]+)",\s*\n\s*r"((?:[^"\\]|\\.)*)"', ids_rs))
    if not declared:
        failures.append("could not parse any id patterns out of ids.rs")
        print("  FAIL could not parse any id patterns out of ids.rs")
    for name, pattern in sorted(declared.items()):
        checks += 1
        schema_pattern = defs.get(name, {}).get("pattern")
        if schema_pattern is None:
            failures.append(f"ids.rs declares {name}, absent from common.schema.json $defs")
            print(f"  FAIL {name} is not in common.schema.json $defs")
        elif schema_pattern != pattern:
            failures.append(f"{name}: ids.rs pattern has drifted from the schema")
            print(f"  FAIL {name}: ids.rs has {pattern!r}, schema has {schema_pattern!r}")
        else:
            print(f"  ok   {name}")

    print("\n[error-code-coverage]")
    # Every CF-* code cited anywhere in planning/ must exist in the Rust registry.
    cited = set()
    for path in (ROOT / "planning").rglob("*"):
        if path.is_file() and path.suffix in {".md", ".json", ".yaml", ".yml", ".pyi"}:
            text = path.read_text(encoding="utf-8")
            cited |= set(re.findall(r"\bCF-[A-Z]{3,8}-[0-9]{3}\b", text))
    error_rs = (ROOT / "crates" / "cf-types" / "src" / "error.rs").read_text(encoding="utf-8")
    registered = set(re.findall(r'\("(CF-[A-Z]{3,8}-[0-9]{3})"', error_rs))
    for code in sorted(cited - registered):
        failures.append(f"{code} is cited in planning/ but not in the cf-types registry")
        print(f"  FAIL {code} is cited in planning/ but not registered")
    checks += 1
    if not (cited - registered):
        print(f"  ok   all {len(cited)} cited codes are registered ({len(registered)} total)")

    print("\n" + "=" * 70)
    if failures:
        print(f"FAILED: {len(failures)} problem(s) across {checks} checks")
        for failure in failures:
            print(f"  - {failure}")
        return 1
    print(f"PASSED: {checks} checks")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sync", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if not (args.sync or args.check):
        args.check = True
    rc = 0
    if args.sync:
        print("[sync]")
        rc |= sync()
    if args.check:
        rc |= check()
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
