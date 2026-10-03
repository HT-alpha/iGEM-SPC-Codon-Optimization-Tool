"""cf-api — the HTTP surface.

`P2-18`, `DEP-API-001`. The server is the standard library. Routes the
library can execute are real; the rest answer `CF-CLI-003`.

Declared dependencies: cf-types, cf-orchestrator, cf-active, cf-modelreg.
Contract: planning/04-api/openapi.yaml
"""

from codonforge.api.server import Api, serve

__all__ = ["Api", "serve"]
