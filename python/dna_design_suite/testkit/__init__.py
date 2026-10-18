"""cf-testkit — synthetic instances whose right answer is known in advance.

Ticket P1-32. Specification: `planning/05-testing/03-synthetic-data-suite.md`.

Phase 1 ships `SYN-01`, `SYN-04`, `SYN-05` and `SYN-07`. `SYN-02`, `SYN-03` and
`SYN-06` need the Potts fitter, the shared resource ledger and the burden
estimators, all of which are Phase 2's, and are absent rather than stubbed.

Both `SYN-01` and `SYN-04` ship *degenerate* configurations alongside their
unambiguous ones, because a suite of instances with one obvious right answer
cannot exercise tie-breaking — and tie-breaking is where determinism bugs live.

Declared dependencies: cf-types, cf-metrics.
"""

from .synthetic import (
    GENERATORS,
    REDUCED_TABLE,
    InfeasibleCase,
    PlantedViolation,
    Strategies,
    ToySpace,
    syn01_degenerate,
    syn01_distinct,
    syn04_multiple_equal_fixes,
    syn04_no_single_codon_fix,
    syn04_unique_fix,
    syn05_cases,
)

__all__ = [
    "GENERATORS",
    "REDUCED_TABLE",
    "InfeasibleCase",
    "PlantedViolation",
    "Strategies",
    "ToySpace",
    "syn01_degenerate",
    "syn01_distinct",
    "syn04_multiple_equal_fixes",
    "syn04_no_single_codon_fix",
    "syn04_unique_fix",
    "syn05_cases",
]
