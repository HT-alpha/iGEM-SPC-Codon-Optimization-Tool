"""One cleaning rule for a primer template.

The page applies the same steps before it draws the sequence or posts it.
Doing it again here is idempotent, so an API client that pastes a FASTA
record and a page that already stripped the header describe the same bases.
"""

from __future__ import annotations


def clean_template(text: str) -> tuple[str, str | None]:
    """Uppercase a template, drop FASTA headers, and name an illegal base.

    Returns `(sequence, None)` or `("", problem)`. The position is 1-based in
    the header-stripped sequence. Lowercase bases are accepted by uppercasing
    them. Whitespace is removed after the header lines are dropped.
    """
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    kept: list[str] = []
    for line in lines:
        trimmed = line.strip()
        if not trimmed or trimmed.startswith(">") or trimmed.startswith(";"):
            continue
        kept.append(trimmed)
    sequence = "".join("".join(kept).upper().split())
    for index, base in enumerate(sequence):
        if base not in "ACGT":
            return "", (
                f"CF-PRIM-006 base {base} at position {index + 1} is not A, C, G or T"
            )
    return sequence, None
