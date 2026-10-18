#!/usr/bin/env python3
"""Decide whether a CI diff has to run the heavy jobs.

`make check` stays the gate for Rust and Python. A diff that is only
documentation skips clippy, the crate tests, the extension build, pytest,
and the generated-asset check. The repository-invariant job still runs.

An empty diff, or a path that escapes `docs/` with `..`, runs the heavy
gate. A filter that cannot see a file must not go green by skipping tests.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Iterable
from pathlib import Path

# Root files a docs-only change may touch. Anything else outside `docs/`
# is code, planning, schemas, or workflow, and the heavy gate runs.
ROOT_DOCS = frozenset(
    {
        "README.md",
        "CONTRIBUTING.md",
        "CHANGELOG.md",
        "AGENTS.md",
        "LICENSE",
    }
)


def normalize_repo_path(path: str) -> str:
    """Lexical repo-relative path. `docs/../python/x.py` is `python/x.py`."""
    slash = path.replace("\\", "/").strip()
    parts: list[str] = []
    for part in slash.split("/"):
        if part in ("", "."):
            continue
        if part == "..":
            if not parts:
                return ""
            parts.pop()
            continue
        parts.append(part)
    return "/".join(parts)


def is_docs_only(path: str) -> bool:
    """True for `docs/**` (including benchmark write-ups) and `ROOT_DOCS`."""
    normal = normalize_repo_path(path)
    if not normal:
        return False
    if normal in ROOT_DOCS:
        return True
    return normal == "docs" or normal.startswith("docs/")


def requires_heavy(paths: Iterable[str]) -> bool:
    """True unless every path is documentation.

    An empty list runs the gate. That is the fail-safe for a diff the
    workflow could not list.
    """
    materialized = list(paths)
    if not materialized:
        return True
    return any(not is_docs_only(path) for path in materialized)


def paths_from_text(text: str) -> list[str]:
    return [line.strip() for line in text.splitlines() if line.strip()]


def github_output_line(paths: Iterable[str]) -> str:
    heavy = requires_heavy(paths)
    return f"heavy={'true' if heavy else 'false'}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "paths_file",
        nargs="?",
        help="newline-separated paths; stdin when omitted",
    )
    args = parser.parse_args(argv)
    if args.paths_file:
        text = Path(args.paths_file).read_text(encoding="utf-8")
    else:
        text = sys.stdin.read()
    paths = paths_from_text(text)
    line = github_output_line(paths)
    print(line)
    print(f"[ci-path-filter] {line} ({len(paths)} paths)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
