"""Open the DNA Design Suite page and wait.

    python -m dna_design_suite.api

Listens on 127.0.0.1 port 8000, opens a browser, and blocks until Ctrl+C.
``--no-browser`` skips the browser. ``--port`` overrides 8000.
"""

from __future__ import annotations

import argparse
import sys

from dna_design_suite.api.server import run_page
from dna_design_suite.modelreg import ArtifactStore


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m dna_design_suite.api",
        description="Open the DNA Design Suite page and wait until Ctrl+C.",
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--store", default=".dna-design-suite/artifacts")
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="do not open a browser; print the address and wait",
    )
    args = parser.parse_args(argv)
    url = f"http://{args.host}:{args.port}/"
    print(f"DNA Design Suite  {url}", flush=True)
    print("Waiting. Press Ctrl+C to stop.", flush=True)
    try:
        run_page(
            ArtifactStore(args.store),
            args.host,
            args.port,
            open_browser=not args.no_browser,
            block=True,
        )
    except OSError as exc:
        print(
            f"CF-CLI-001: could not listen on {args.host}:{args.port}. {exc}",
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
