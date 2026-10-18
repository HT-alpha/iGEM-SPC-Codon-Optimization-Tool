#!/usr/bin/env python3
"""Choose which wheels a run builds.

Ordinary pushes do not call this. The workflow runs on a manual dispatch
and on a version tag. Dispatch defaults to Linux x86_64. `all`, and any
run that is not a dispatch, builds every platform Hugo publishes.
"""

from __future__ import annotations

import argparse
import json
import sys

LINUX_X86_64 = {"os": "ubuntu-latest", "target": "x86_64", "manylinux": "auto"}
LINUX_AARCH64 = {"os": "ubuntu-latest", "target": "aarch64", "manylinux": "auto"}
MACOS_X86_64 = {"os": "macos-latest", "target": "x86_64"}
MACOS_AARCH64 = {"os": "macos-latest", "target": "aarch64"}
WINDOWS_X64 = {"os": "windows-latest", "target": "x64"}

# Order is stable so a serialized matrix does not churn between runs.
ALL_PLATFORMS: list[dict[str, str]] = [
    LINUX_X86_64,
    LINUX_AARCH64,
    MACOS_X86_64,
    MACOS_AARCH64,
    WINDOWS_X64,
]


def select_platforms(event_name: str, targets: str) -> list[dict[str, str]]:
    """Linux x86_64 unless this is a publish.

    A version tag is a publish. A dispatch whose target input is `all` is
    a publish. Any other dispatch value, including an empty or unknown
    one, stays on the cheap default so a bad input cannot start macOS.
    """
    if event_name != "workflow_dispatch" or targets == "all":
        return [dict(platform) for platform in ALL_PLATFORMS]
    return [dict(LINUX_X86_64)]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event-name", required=True)
    parser.add_argument("--targets", default="")
    args = parser.parse_args(argv)
    json.dump(
        select_platforms(args.event_name, args.targets),
        sys.stdout,
        separators=(",", ":"),
    )
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
