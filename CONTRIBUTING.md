# Contributing to DNA Design Suite

## The two toolchains

DNA Design Suite is a Rust core (`crates/`) exposed to Python (`python/`) through a
single PyO3 extension. You need both toolchains: a stable Rust compiler and
Python 3.11 or later.

```sh
make setup     # virtualenv + Python dev dependencies
make build     # maturin develop: compiles the core into the virtualenv
```

### Why `make setup` avoids `python3 -m venv`

`python3 -m venv` needs `ensurepip`, which Debian and Ubuntu ship in a separate
`python3-venv` package. That package is not always installed, and on some
images it is not installable at all. When it is missing, the error names a
failing command inside `.venv/bin/` and reads like a broken checkout rather
than a missing system package.

`make setup` therefore uses `virtualenv`, which vendors its own bootstrap,
installing it first if it is absent, and only falls back to `python3 -m venv`
if that fails. **This is the supported path, not an implementation detail.** If
you create an environment by hand, use `virtualenv .venv`. If you add tooling,
CI or documentation that creates one, use `virtualenv` there too — a bare
`python3 -m venv` reintroduces a failure that costs the next contributor an
afternoon.

## Layout

| Directory | Contents |
| --- | --- |
| `crates/` | The Rust core, one crate per architectural module |
| `python/dna_design_suite/` | Ingest, model registry, orchestration, decision, CLI |
| `planning/` | The normative planning package: schemas, ADRs, specs, tickets |
| `schemas/v1/` | Runtime copies of the schemas, synced from `planning/02-schemas/` |
| `tools/` | Repository invariant checks |

`crates/cf-pyext` is the only crate permitted to depend on `pyo3`, and it is the
only place the Rust/Python boundary is crossed. That is enforced rather than
assumed — see below.

## Module layering

`planning/01-modules/modules.yaml` is the dependency graph, and it is normative.
A module may depend only on modules it declares, and a Rust core module may not
depend on a Python orchestration one.

```sh
python3 tools/check_module_layering.py --check
```

This will reject an import you did not declare. When it does, the question to
ask is whether the dependency is right — the answer is often to pass a value
rather than acquire a module. A shared constant between two modules that do not
declare each other is a dependency wearing a disguise.

## Before you push

```sh
make lint      # cargo clippy -D warnings, ruff check
make test      # cargo test --workspace, pytest
.venv/bin/python planning/tools/validate_planning_package.py
python3 tools/check_schema_sync.py --check
```

All five must pass. `check_schema_sync` also verifies that every `CF-*` code
cited in `planning/` is registered in `cf-types`, so inventing a code in prose
without defining it is caught.

## Continuous integration

CI is `ubuntu-latest` only. It runs on pull requests and on pushes to
`main`. A feature-branch push is not a second CI run beside the pull
request. `concurrency` puts a pull request and a push of the same branch
in one group and cancels the older run.

`make check` stays the gate when a diff touches Rust or Python. That gate
is clippy, the per-crate tests, the `cf-pyext` build, pytest on 3.11 and
3.12, the generated-asset check, and the repository invariants. A diff
that is only documentation does not run the heavy jobs. Documentation
means `docs/` (including `docs/benchmarks/`) and the root files
`README.md`, `CONTRIBUTING.md`, `CHANGELOG.md`, `AGENTS.md`, and
`LICENSE`. The invariant job still runs. Anything else — `crates/`,
`python/`, `planning/`, `schemas/`, `tools/`, `examples/`, a workflow —
runs the full gate. An empty or unreadable diff runs the full gate too,
so a filter failure cannot go green by skipping the tests.

The classifier is `tools/ci_path_filter.py`. Cargo caches are
`Swatinem/rust-cache`, including clippy. pip caches through
`actions/setup-python`, keyed on `pyproject.toml`.

## Building wheels

Wheels are not built on a push, and there is no weekly schedule. Those
runs were spending private-repository Actions minutes on five platform
builds for commits that were never published. macOS and Windows runners
cost more than Linux, which is why they are opt-in.

To build one: **Actions → wheels → Run workflow**. The default target is
Linux x86_64, plus an sdist, plus a check that the Linux wheel installs
with no Rust toolchain. Choose **all** when you are about to publish.
That adds Linux aarch64, macOS x86_64, macOS aarch64, and Windows x64.

Pushing a version tag `v*.*.*` (for example `v0.1.0`) builds that same
full set. A tag with a hyphen (`v0.1.0-igem`) skips `legal-clearance`,
and every wheel job needs that job, so nothing is built. Creating a
GitHub Release does not start the workflow. Publish from the tag run,
or dispatch the workflow yourself.

The workflow uploads artifacts. It does not publish to PyPI. Download
`wheel-*` and `sdist`, then `twine upload` from a machine that has your
token. The token stays out of the repository. `tools/wheel_targets.py`
is the platform list; the workflow reads it. Wheel compiles use
`sccache` through `PyO3/maturin-action`.

## What a change is expected to carry

**Tests that would fail without it.** The bar here is specifically that a test
must fail when the behaviour it describes is broken. Several times in this
codebase a test has passed against a deliberately broken implementation because
its fixture could not reach the defect — a determinism fixture with a unique
optimum never invoked tie-breaking; a repair fixture whose motif never arose
left the tabu ledger empty in every generation. If you are unsure, break the
code on purpose and check that your test notices.

**Schema and planning updates in the same commit.** If behaviour disagrees with
`planning/`, the planning artifact is wrong too, and both get fixed together
with a registered `DEP-*` identifier in
`planning/00-scope/04-departure-index.md`. A correction that lives only in code
is a correction the next reader will not find.

**A `CF-*` code for anything a user can hit.** Errors and warnings carry a
stable code, a `stage` and a `subject`, so a caller can attribute one without
parsing prose. Register it in `crates/cf-types/src/error.rs`.

## Determinism

Reproducibility is a correctness property here, not a nice-to-have.
`GATE-DET-001` requires bit-identical output across 1, 4 and 16 threads.

Two rules follow, and both have been violated at least once:

**Derive random streams from a structural index, never from iteration or
completion order.** `SearchRng::stream(root_seed, index)` where `index` is a
block, generation or individual number. A stream keyed on the order work
happened to be scheduled in is not reproducible across thread counts.

**Order any collection that reaches a serialized document.** Rust seeds
`HashMap` iteration per process, so an unordered snapshot produces different
bytes on different runs — which surfaces as an intermittent failure, the worst
kind to diagnose. Sort it, and assert the sort in a test rather than relying on
incidental behaviour.

## Checkpoints

**Anything the reported result depends on must travel in the checkpoint.**

This rule exists because it was learned three times: the sealed lockfile, the
stage reports and the candidate archive were each omitted, and each omission
was found by a different failing test. A resumed job that recomputes something
instead of restoring it will silently report a different answer from an
uninterrupted one — and both will look internally consistent.

If you add state that influences the output, add it to
`Job.persist_checkpoint` and `Job.resume` in the same change.
