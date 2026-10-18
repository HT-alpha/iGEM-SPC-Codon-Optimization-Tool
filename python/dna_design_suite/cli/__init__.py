"""cf-cli — the `dna-design-suite` command-line entry point.

A thin argument-parsing shell over `cf-orchestrator`: every verb maps 1:1 to a
library call, and this package holds no validation, no defaults and no
objective presets of its own. Each of those would be a second code path, and a
second code path is how three surfaces that are supposed to agree stop
agreeing.

Declared dependencies: cf-types, cf-orchestrator, cf-active, cf-modelreg.
"""

import importlib
from typing import Any

__all__ = ["build_parser", "main"]


def __getattr__(name: str) -> Any:
    """Import `main` on first use.

    `python -m dna_design_suite.cli.main` imports the package before it executes the
    module. Importing `main` here would leave `dna_design_suite.cli.main` in
    `sys.modules` before that execution, which is the RuntimeWarning every
    CLI invocation used to print. `importlib.import_module` loads the submodule
    without looking the name up on this package, which is what recurses.
    """
    if name in __all__:
        module = importlib.import_module(f"{__name__}.main")
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
