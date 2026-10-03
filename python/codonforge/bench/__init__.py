"""cf-bench — benchmark harness.

`P2-19`. Test-only in the sense that the shipped library does not import it
to design a sequence. It does call the library's scoring path.

A fixed reference point, comparison refusal, in-process native and CAI-greedy
baselines, and an external-tool adapter that runs a binary the user installed.
No tool is vendored. The real-protein panel is 13 CDS from NC_000913.3
(`DEP-BEN-002`). It does not claim `GATE-BEN-001`.

Declared dependencies: cf-types, cf-metrics, cf-orchestrator, cf-testkit, cf-ingest.
Contract: planning/05-testing/04-benchmark-harness.md. `DEP-BEN-001`.
"""

from codonforge.bench.baselines import ExternalTool, cai_greedy, native, run_external
from codonforge.bench.compare import RunRecord, assert_comparable
from codonforge.bench.exp05 import Decision, Trial, decide
from codonforge.bench.reference import BenchError, fixed_reference_point, hypervolume
from codonforge.bench.results import check_reproduces, render_report, write_result

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
