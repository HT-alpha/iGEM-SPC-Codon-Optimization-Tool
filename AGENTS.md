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
