"""Anchor-verifying text edits.

Every edit here fails loudly when its anchor is missing. `str.replace` does
not: it returns the text unchanged and the caller carries on, so an edit whose
anchor drifted — because a formatter rewrapped a signature, say — silently does
nothing and is discovered later from a confusing test failure rather than from
the edit step. That has happened twice in this repository, both times costing a
misdiagnosed negative control, so scripted edits go through here.
"""

from __future__ import annotations

import pathlib
import sys


class AnchorError(RuntimeError):
    """An edit whose anchor was absent, or present more times than expected."""


def sub(path: str | pathlib.Path, old: str, new: str, *, count: int = 1) -> None:
    """Replace `old` with `new`, requiring exactly `count` occurrences."""
    target = pathlib.Path(path)
    text = target.read_text(encoding="utf-8")
    found = text.count(old)
    if found != count:
        raise AnchorError(
            f"{target}: expected {count} occurrence(s) of the anchor, found {found}\n"
            f"--- anchor ---\n{old[:400]}"
        )
    target.write_text(text.replace(old, new, count), encoding="utf-8")


def append(path: str | pathlib.Path, text: str) -> None:
    target = pathlib.Path(path)
    target.write_text(target.read_text(encoding="utf-8") + text, encoding="utf-8")


def main() -> int:
    print(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
