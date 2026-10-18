"""cf-bench — benchmark harness.

`P2-19`. Test-only in the sense that the shipped library does not import it
to design a sequence. It does call the library's scoring path.

A fixed reference point, comparison refusal, in-process native and CAI-greedy
baselines, and an external-tool adapter that runs a binary the user installed.
No tool is vendored. The real-protein panel is 50 CDS from NC_000913.3
(`DEP-BEN-002`). `GATE-BEN-001` is claimed only as the CAI-only record in
`benchmarks/corpus/gates.json`: 50 of 50 proteins. LinearFold's global MFE
is not `structure.local_accessibility`, and that comparison is not the claim.

Declared dependencies: cf-types, cf-metrics, cf-orchestrator, cf-testkit, cf-ingest.
Contract: planning/05-testing/04-benchmark-harness.md. `DEP-BEN-001`.
"""

from dna_design_suite.codon.bench.baselines import ExternalTool, cai_greedy, native, run_external
from dna_design_suite.codon.bench.compare import RunRecord, assert_comparable
from dna_design_suite.codon.bench.exp05 import Decision, Trial, decide
from dna_design_suite.codon.bench.reference import BenchError, fixed_reference_point, hypervolume
from dna_design_suite.codon.bench.results import check_reproduces, render_report, write_result

__all__ = [
    "BenchError",
    "Decision",
    "ExternalTool",
    "RunRecord",
    "Trial",
    "assert_comparable",
    "cai_greedy",
    "check_reproduces",
    "decide",
    "fixed_reference_point",
    "hypervolume",
    "native",
    "render_report",
    "run_external",
    "write_result",
]
