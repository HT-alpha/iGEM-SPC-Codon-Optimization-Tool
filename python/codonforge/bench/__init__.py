"""cf-bench — support layer.

Benchmark harness: reference datasets, baseline adapters, indicator computation,
threshold gates, and the regression-performance suite. Test-only target; never a runtime
dependency of the shipped library.

Declared dependencies: cf-types, cf-metrics, cf-orchestrator, cf-testkit, cf-ingest.
Contract: planning/05-testing/04-benchmark-harness.md
"""
