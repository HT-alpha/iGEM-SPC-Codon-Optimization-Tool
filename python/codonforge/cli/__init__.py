"""cf-cli — the `codonforge` command-line entry point.

A thin argument-parsing shell over `cf-orchestrator`: every verb maps 1:1 to a
library call, and this package holds no validation, no defaults and no
objective presets of its own. Each of those would be a second code path, and a
second code path is how three surfaces that are supposed to agree stop
agreeing.

Declared dependencies: cf-types, cf-orchestrator, cf-active, cf-modelreg.
"""

from .main import build_parser, main

__all__ = ["build_parser", "main"]
