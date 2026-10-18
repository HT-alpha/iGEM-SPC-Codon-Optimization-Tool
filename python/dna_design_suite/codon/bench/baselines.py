"""Baselines the harness can run without vendoring a tool. `P2-19`, `AMD-013`.

`native` is the input CDS. `cai_greedy` picks, at each position, the synonym
with the highest `relative_adaptiveness` in the host profile the caller
already supplied. A tie takes the lexicographically smaller codon. A codon
absent from the table has weight 0, which is "no observation", not a guessed
frequency.

An external tool is whatever executable the user points at. This module does
not name LinearDesign's flags: the argv template is the user's, and it must
contain `{input}` so the sequence file is actually passed. Stdout must be a
FASTA record the existing parser accepts, and the translation must match the
input. A tool that prints a different protein is rejected.
"""

from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

from dna_design_suite import _core
from dna_design_suite.codon.bench.reference import BenchError
from dna_design_suite.ingest.fasta import parse_fasta


def native(cds: str) -> str:
    return cds


def cai_greedy(cds: str, translation_table: int, usage: dict[str, float]) -> str:
    """One codon at a time. Deterministic. Uses the caller's usage table."""
    if len(cds) % 3 != 0:
        raise BenchError("CAI-greedy input is not a whole number of codons")
    codons = []
    for index in range(0, len(cds), 3):
        codon = cds[index : index + 3].upper()
        options = _core.synonyms(codon, translation_table)
        chosen = min(options, key=lambda item: (-float(usage.get(item, 0.0)), item))
        codons.append(chosen)
    return "".join(codons)


@dataclass(frozen=True)
class ExternalTool:
    """A user-installed executable. The binary is never shipped."""

    name: str
    executable: str
    argv: tuple[str, ...]
    timeout_seconds: float = 60.0

    def available(self) -> bool:
        path = Path(self.executable)
        return path.is_file() and os.access(path, os.X_OK)

    def command(self, input_path: Path) -> list[str]:
        if not any("{input}" in arg for arg in self.argv):
            raise BenchError(
                f"{self.name}: the argv template does not contain {{input}}, so the "
                "sequence would never be passed to the tool"
            )
        if not self.available():
            raise BenchError(
                f"{self.name}: {self.executable} is not an executable this user installed"
            )
        rendered = [arg.replace("{input}", str(input_path)) for arg in self.argv]
        return [self.executable, *rendered]


def run_external(
    tool: ExternalTool,
    cds: str,
    translation_table: int,
    workdir: Path,
) -> str:
    """Run the tool and accept the result only if it encodes the same protein."""
    workdir.mkdir(parents=True, exist_ok=True)
    source = workdir / "input.fa"
    source.write_text(f">input\n{cds}\n", encoding="utf-8")
    argv = tool.command(source)
    completed = subprocess.run(
        argv,
        check=False,
        capture_output=True,
        text=True,
        timeout=tool.timeout_seconds,
        cwd=workdir,
    )
    if completed.returncode != 0:
        raise BenchError(f"{tool.name} exited {completed.returncode}: {completed.stderr.strip()}")
    parsed = parse_fasta(completed.stdout)
    if len(parsed) != 1:
        raise BenchError(f"{tool.name} stdout is not a single FASTA record ({len(parsed)} records)")
    produced = parsed[0].sequence
    expected = _core.translate(cds, translation_table)
    found = _core.translate(produced, translation_table)
    if found != expected:
        raise BenchError(f"{tool.name} changed the protein; the output is not a synonymous design")
    return produced
