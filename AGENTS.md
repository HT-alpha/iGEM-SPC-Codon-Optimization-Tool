# AGENTS.md

Working notes for agents and engineers changing this repository. Everything
here is enforced by a check, learned from a defect that shipped, or both.

## Build, test, run

```sh
make setup     # virtualenv + Python dev dependencies (NOT python3 -m venv; see below)
make build     # maturin develop: compiles the Rust core into the virtualenv
make test      # cargo test --workspace --all-targets, then pytest python/tests
make lint      # cargo clippy -D warnings, ruff check
```

Run a design end to end:

```sh
.venv/bin/python -m codonforge.cli.main design \
    --host escherichia_coli:511145/BL21 --cds <CDS> \
    --host-profile host.json --constraints constraints.json \
    --objectives usage.cai --seed 20260828 --output explain
```

Repository invariants, all of which must pass:

```sh
.venv/bin/python planning/tools/validate_planning_package.py
python3 tools/check_module_layering.py --check
python3 tools/check_schema_sync.py --check
```

## Environment

**Use `virtualenv`, never a bare `python3 -m venv`.** The stdlib module needs
`ensurepip`, which Debian and Ubuntu ship separately in `python3-venv` — a
package that is not always present and sometimes not installable. When it is
missing the failure names a command inside `.venv/bin/` and looks like a broken
checkout. `make setup` handles this; if you add tooling or CI that creates an
environment, keep it handling it.

## Where things live

`planning/` is normative. Schemas, algorithm specifications, ADRs, the test
plan and the ticketed build sequence live there, and code is expected to match
them. When it does not, the disagreement is a bug in one of the two and both
get fixed in the same change, with a `DEP-*` identifier registered in
`planning/00-scope/04-departure-index.md`.

`crates/` is the Rust core, one crate per module in
`planning/01-modules/modules.yaml`. `python/codonforge/` is ingest, the model
registry, orchestration, decision and the CLI. `crates/cf-pyext` is the only
crate allowed to depend on `pyo3` and the only place the boundary is crossed —
enforced by `check_module_layering.py`.

## Rules that have teeth

**A test must fail when the thing it describes is broken.** This has gone wrong
often enough to be worth stating: a determinism fixture with a unique optimum
never invoked tie-breaking, so an implementation keying RNG streams on iteration
order passed it; a repair fixture whose forbidden motif never arose left the
tabu ledger empty in every generation, making every assertion about checkpoint
persistence vacuous. Break the code on purpose and confirm your test notices.

**Anything the reported result depends on must travel in the checkpoint.**
Learned three times — the sealed lockfile, the stage reports, the candidate
archive — each found by a different failing test. A resumed job that recomputes
rather than restores reports a different answer from an uninterrupted one, and
both look internally consistent.

**Derive RNG streams from a structural index, never from iteration order**, and
**order any collection that reaches a serialized document**. `HashMap`
iteration is seeded per process, so an unordered snapshot fails intermittently.

**One producer per value.** `objective_sense`, `resolved_objectives[].direction`
and `progress.tiers` each have exactly one source, and readers render them
rather than recomputing them. Tier 1 and tier 2 once each decided what to do
with an unscoreable term and reached different answers; a value computed twice
can disagree with itself.

**No fabricated constants.** If a value is not measured, derived, or read from
a specification, it does not go in. A capability a caller asserts is not a
capability.

## Reporting a change

Per ticket: what changed, the acceptance criterion and its status,
negative-control evidence (inject the defect, confirm the test fires, revert),
and anything the reader should know that the diff does not show. If a negative
control does not fire on the first try, assume the fixture is inadequate before
concluding the check is redundant — that assumption has been correct every time
so far.

## Cursor Cloud environment

`main` (`ae6822c`, 2026-08-05) is the original stub: `README.md`, `codon_6.0.py`,
`old_codon.py`. It is an ancestor of the product, not a second codebase.
`git rev-list --count origin/main..v2.2.0` is 145, and the other direction is 0.
Every release tag from `v0.1.0-igem` through `v2.3.0` is on the line that runs
through `origin/cursor/v21-open-items-b1dc` (`v2.2.0`, `28720f0`). GitHub's
default branch is still `main`, so a fresh clone does not contain those tags.
Fast-forwarding `main` to the current stable tag would not rewrite history.
Do not do it unless someone has confirmed it. The command, once confirmed, is
`git push origin <stable-tag>:main`. Do not rebase the release line onto `main`.

Two branches from July 2026 are not on this line:
`origin/cursor/build-bioforge-1637` and
`origin/cursor/setup-dev-environment-a1db`.
`origin/cursor/p1-26-reread-recommended-front-design-6c83` has one commit that
is not in the release line. `P1-26` stays held, so that branch is not merged
from here.

ViennaRNA is installed by the user, never vendored. The structure tests skip
without it, and `CF_REQUIRE_STRUCTURE_ENGINE=1` turns that skip into a failure.
This environment has ViennaRNA 2.7.2, from:

```sh
.venv/bin/pip install 'ViennaRNA==2.7.2'
```

The helix tests and the viewer-selection test run JavaScript under Node.
`CF_REQUIRE_NODE=1` turns a missing Node into a failure. `node` here is
v22.14.0 at `/exec-daemon/node`.

The API server has no CLI entry point. From the repo root, after `make build`:

```sh
.venv/bin/python -c '
from pathlib import Path
import time
from codonforge.api import serve
from codonforge.modelreg import ArtifactStore
serve(ArtifactStore(Path("/tmp/codonforge-artifacts")), port=8765)
print("http://127.0.0.1:8765/", flush=True)
while True:
    time.sleep(3600)
'
```

No environment variable is required. Static files are read from
`python/codonforge/api/static` on each request. A new path is served only after
the process restarts, because the route table is built at import.

MCP sign-in does not complete under automation. Authenticating the Adapter
server timed out, and there is no workaround yet, so a run here cannot read or
write that workspace.
