"""cf-cli — support layer.

`codonforge` command-line entry point. A thin argument-parsing shell over `cf-
orchestrator`; every CLI verb maps 1:1 to a library call.

Declared dependencies: cf-types, cf-orchestrator, cf-active, cf-modelreg.
Contract: planning/04-api/02-cli-and-library.md
"""
