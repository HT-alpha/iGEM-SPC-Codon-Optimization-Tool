#!/usr/bin/env python3
"""Recompute and re-thread the digests in the worked schema examples.

    python3 planning/tools/refresh_example_digests.py [--check] [-v]

`--check` reports what would change without writing, which is what the
`digests` group of `validate_planning_package.py` uses.

The classification lives in `planning/02-schemas/examples/digest_map.yaml`; this
tool is only the mechanism. It resolves the reference graph topologically, so an
example that is cited by another is finalized before the one citing it, and it
iterates the whole set to a fixed point as a safety net — a self digest changes
when a threaded reference inside it changes, so one pass is not enough for a
document that both cites and is cited.

Canonicalization uses `planning/tools/jcs.py`, the from-spec reference
implementation, so this tool needs no compiled core and runs in a CI job with no
Rust toolchain. `python/tests/test_contract_sequence_schema.py` separately
asserts that the Rust core agrees with it, so a divergence between the two
cannot hide here.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import jcs  # noqa: E402

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit("PyYAML is required: pip install pyyaml")

PLANNING = Path(__file__).resolve().parent.parent
ROOT = PLANNING.parent
EXAMPLES = PLANNING / "02-schemas" / "examples"
MAP_PATH = EXAMPLES / "digest_map.yaml"

MAX_PASSES = 8


# --------------------------------------------------------------------------- #
# RFC 6901 JSON Pointer
# --------------------------------------------------------------------------- #


def unescape(token: str) -> str:
    return token.replace("~1", "/").replace("~0", "~")


def pointer_parts(pointer: str) -> list[str]:
    if pointer in ("", "/"):
        return []
    if not pointer.startswith("/"):
        raise ValueError(f"not a JSON Pointer: {pointer!r}")
    return [unescape(t) for t in pointer[1:].split("/")]


def resolve(document, pointer: str):
    node = document
    for token in pointer_parts(pointer):
        node = node[int(token)] if isinstance(node, list) else node[token]
    return node


def resolve_parent(document, pointer: str):
    parts = pointer_parts(pointer)
    node = document
    for token in parts[:-1]:
        node = node[int(token)] if isinstance(node, list) else node[token]
    return node, parts[-1]


def assign(document, pointer: str, value) -> None:
    parent, last = resolve_parent(document, pointer)
    if isinstance(parent, list):
        parent[int(last)] = value
    else:
        parent[last] = value


# --------------------------------------------------------------------------- #
# the pass
# --------------------------------------------------------------------------- #


def load_examples() -> dict[str, dict]:
    return {
        path.name: json.loads(path.read_text(encoding="utf-8"))
        for path in sorted(EXAMPLES.glob("*.example.json"))
    }


def order_computed(spec: dict) -> list[dict]:
    """Self digests innermost-first.

    A document that embeds another (a job request embedding a construct
    context) must have the inner digest settled before the outer one is taken
    over a payload that contains it. Deeper pointer first does exactly that.
    """
    return sorted(spec["computed"], key=lambda e: -len(pointer_parts(e["at"])))


def apply_once(docs: dict[str, dict], spec: dict, tables: dict) -> list[str]:
    """One pass. Returns the list of paths whose value changed."""
    changed: list[str] = []

    def set_if_different(file: str, at: str, value: str) -> None:
        current = resolve(docs[file], at)
        if current != value:
            assign(docs[file], at, value)
            changed.append(f"{file}{at}")

    for entry in order_computed(spec):
        doc = docs[entry["file"]]
        rule = entry["rule"]
        if rule == "self_document":
            parent, _ = resolve_parent(doc, entry["at"])
            value = jcs.digest_of_document(parent, "digest")
        elif rule == "codon_sequence_document":
            cds = resolve(doc, entry["cds"])
            table = entry["table_const"] if "table_const" in entry else resolve(doc, entry["table"])
            value = jcs.codon_sequence_document(cds, int(table), tables)["digest"]
        elif rule == "raw_bytes":
            value = jcs.digest_of_bytes(resolve(doc, entry["source"]).encode("utf-8"))
        else:
            raise ValueError(f"unknown rule {rule!r}")
        set_if_different(entry["file"], entry["at"], value)

    for entry in spec["reference"]:
        target_file, _, target_pointer = entry["same_as"].partition("#")
        value = resolve(docs[target_file], target_pointer)
        set_if_different(entry["file"], entry["at"], value)

    return changed


def refresh(check_only: bool, verbose: bool) -> int:
    spec = yaml.safe_load(MAP_PATH.read_text(encoding="utf-8"))
    tables = jcs.load_genetic_codes(ROOT)
    docs = load_examples()
    original = {name: json.dumps(doc, sort_keys=True) for name, doc in docs.items()}

    all_changed: list[str] = []
    for iteration in range(MAX_PASSES):
        changed = apply_once(docs, spec, tables)
        if verbose and changed:
            print(f"  pass {iteration + 1}: {len(changed)} value(s) updated")
        all_changed.extend(changed)
        if not changed:
            break
    else:
        print(
            f"  FAIL did not reach a fixed point in {MAX_PASSES} passes; "
            "the reference graph is probably cyclic"
        )
        return 1

    touched = [
        name for name, doc in docs.items() if json.dumps(doc, sort_keys=True) != original[name]
    ]

    if check_only:
        if touched:
            print(f"  FAIL {len(set(all_changed))} digest(s) are stale in {len(touched)} file(s)")
            for path in sorted(set(all_changed)):
                print(f"    {path}")
            print("  run: python3 planning/tools/refresh_example_digests.py")
            return 1
        print(f"  ok   every digest in {len(docs)} example(s) is current")
        return 0

    for name in touched:
        path = EXAMPLES / name
        path.write_text(json.dumps(docs[name], indent=2) + "\n", encoding="utf-8")
        print(f"  rewrote {path.relative_to(ROOT)}")
    if not touched:
        print("  already current")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    return refresh(args.check, args.verbose)


if __name__ == "__main__":
    raise SystemExit(main())
