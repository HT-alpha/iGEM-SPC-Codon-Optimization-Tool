"""cf-ingest — context-ingestion layer.

FASTA/GenBank parsing to the canonical representation, construct-context
assembly, host reference-data loading, and the external-tool seed importers.

Declared dependencies: cf-types, cf-modelreg.
Contract: planning/02-schemas/construct_context.schema.json

Tickets: P1-09 (FASTA/GenBank readers), P1-10 (construct context),
P1-11 (host profile loader), P1-12 (seed importers).

Ingest validates by *delegating*, not by re-implementing. `CodonSequence`
already defines what a valid sequence is and reports each failure with a
specific `CF-SEQ-*` code; a second definition here would drift, and the first
symptom would be a file ingest accepts and the core later rejects. So every
sequence goes through `codonforge._core`, and the core's problem detail is
translated into an `IngestError` that keeps the code and adds the location in
the file.
"""

from __future__ import annotations

from codonforge.ingest.errors import IngestError, IngestReport, SourceLocation
from codonforge.ingest.fasta import (
    FastaRecord,
    normalize_nucleotides,
    parse_fasta,
    read_cds,
    read_cds_file,
)
from codonforge.ingest.genbank import extract_cds_sequences, read_construct
from codonforge.ingest.host import build_host_profile, codon_usage_from_cds
from codonforge.ingest.seeds import KNOWN_TOOLS, import_seeds

__all__ = [
    "FastaRecord",
    "IngestError",
    "IngestReport",
    "KNOWN_TOOLS",
    "SourceLocation",
    "build_host_profile",
    "codon_usage_from_cds",
    "extract_cds_sequences",
    "import_seeds",
    "normalize_nucleotides",
    "parse_fasta",
    "read_cds",
    "read_cds_file",
    "read_construct",
]
