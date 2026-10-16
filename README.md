# DNA Design Suite

CodonForge is the codon-optimization tool in this suite. PrimerForge is the
primer-design tool. They share one page and one visual language.

Epistasis-aware, systems-level, multi-objective codon optimization.

CodonForge treats a coding sequence as a design problem with several objectives
that genuinely conflict — codon adaptation, local mRNA structure, host resource
burden, and hard sequence constraints — rather than as a table lookup. The
search runs in Rust; the job lifecycle, ingest and CLI run in Python.

**This is Phase 1.** What works is described below, and so is what does not.
The [capability boundaries](#what-does-not-work-yet) section is not a roadmap:
it lists things you might reasonably expect to work, that do not, and what
happens if you try.

---

## Quickstart

### Install a wheel

ViennaRNA is optional and installed by the user. CodonForge does not vendor it.

```sh
pip install ViennaRNA
pip install codonforge
```

`rpsL` below is the *E. coli* K-12 MG1655 CDS `NP_417801.1` from `NC_000913.3`,
125 codons. The usage table is the relative adaptiveness counted from the
ribosomal-protein CDS in that same file (`rpl`/`rps`/`rpm`), labelled
`user_supplied`.

```sh
python -m codonforge.cli.main design \
    --host escherichia_coli:511145/MG1655 \
    --cds "$(tr -d '\n' < examples/rpsL.cds)" \
    --host-profile examples/mg1655_ribosomal_usage.json \
    --constraints examples/rpsL_constraints.json \
    --objectives usage.cai \
    --budget fast \
    --generations 4 \
    --seed 7 \
    --output explain
```

### Open the page

CodonForge and PrimerForge are one page. From a checkout, after `make setup`
and `make build`:

```sh
.venv/bin/python -m codonforge.api
```

The process listens on port 8000, opens a browser at the page, and waits until
Ctrl+C. `--no-browser` skips the window. `--port` chooses a different port;
8000 is the default.

On Windows, in PowerShell, from the repository directory:

```powershell
.\.venv\Scripts\python.exe -m codonforge.api
```

Ctrl+C stops it. The same module command works on macOS and Linux when
`python` is the environment you installed into, including after a wheel
install:

```sh
python -m codonforge.api
```

The form pastes or uploads a coding sequence, takes the host, and can load
the bundled ribosomal usage table in one click. The rpsL demo fills the
sequence, that table, and the example constraints. A result can be copied or
saved as FASTA, GenBank, or CSV. PrimerForge is the other tab on the same page.

### 1. Build from a checkout

CodonForge is a Rust core with a Python extension, so a build needs both
toolchains. `make setup` creates a virtualenv and installs the Python side;
`make build` compiles the Rust core into it.

```sh
git clone <this repository> && cd codonforge
make setup
make build
```

> **On the virtualenv.** `make setup` uses `virtualenv`, installing it first if
> it is missing, and only falls back to the standard library's `python3 -m venv`
> if that fails. This is deliberate and load-bearing: `python3 -m venv` requires
> `ensurepip`, which Debian and Ubuntu ship in a separate `python3-venv` package
> that is not always installed and not always installable. When it is missing
> the failure reads like a broken checkout rather than a missing package.
> `virtualenv` vendors its own bootstrap and needs neither.
>
> If you create the environment yourself, use `virtualenv .venv` rather than
> `python3 -m venv .venv` for the same reason.

### 2. Design a sequence

You need a host codon-usage profile and, optionally, a constraint set. Both are
JSON documents validated against the schemas in [`schemas/v1`](schemas/v1).

```sh
.venv/bin/python -m codonforge.cli.main design \
    --host escherichia_coli:511145/BL21 \
    --cds ATGCTGCGTAGCGGCGTGGCGCCGACCAAACTGCGTAGCTAA \
    --host-profile host.json \
    --constraints constraints.json \
    --objectives usage.cai,usage.gc_content \
    --seed 20260828 \
    --output explain
```

```
candidate  cand_60f6aa72c5038562
job        job_HR9ECVXNGDE8QY9K680EZG3H8M  succeeded
lockfile   sha256:e68367b9b4ee75a84afc2b0c8787c6eead5846eabd7243a9f0a0c88c6baef485
origin     tier2_block_ca

objectives
  usage.cai               0.864322  maximize
  usage.gc_content        0.428571  minimize

tiers
  tier1_seed      done
  tier2_block_ca  done
  tier3_moea      done  8/8

constraints  satisfied
```

Other formats: `--output json` (the schema-validated documents, verbatim),
`fasta`, `genbank`, `csv`, `table`.

### 3. Or from Python

```python
import codonforge as cf
from codonforge.modelreg import ArtifactStore
from codonforge.orchestrator import DesignRequest

result = cf.design(
    DesignRequest(
        host="escherichia_coli:511145/BL21",
        cds="ATGCTGCGTAGCGGCGTGGCGCCGACCAAACTGCGTAGCTAA",
        host_profile=host_profile,
        objectives=["usage.cai"],
        constraints=constraint_set,
        seed=20260828,
    ),
    ArtifactStore("artifacts"),
)

result.front            # a pareto_front.schema.json document
result.status           # a job_status.schema.json document
result.best()           # the reported candidate
```

The CLI builds exactly this `DesignRequest` from its flags and calls exactly
this function. There is no second validation path and no second set of
defaults.

---

## What it does

A design job runs in stages, and the job reports **the best candidate any stage
produced** — not the last stage's output.

| Stage | What it does |
| --- | --- |
| `tier1_seed` | Constrained beam DP over the coding sequence; exact for an additive objective |
| `tier2_block_ca` | Red-black block coordinate ascent, parallel across same-colour blocks |
| `tier3_moea` | Mutation and repair over generations |

Stages are not required to improve on each other, because exploration is not
monotone. When the reported candidate did not come from the last stage to run,
the job says so with a `CF-ORCH-007` warning naming both.

Every run is pinned: the landscape lockfile is sealed before any search work
begins, and the same request with the same seed reproduces the same design. A
job can be checkpointed, cancelled with `SIGINT`, and resumed with `--resume`.

---

## What does not work yet

Each of these produces a specific, named error or warning rather than failing
quietly.

**Background execution.** `cf.design(..., background=True)` raises
`CF-CLI-003`. It needs a worker pool that does not exist. It does not silently
run synchronously, because that would make every timing assumption in your code
wrong.

**A multi-objective Pareto front.** A Phase 1 front is a valid
`pareto_front.schema.json` document with `reduction.method: "none"`, and a
single-objective job produces a front of one. The document shape, the objective
order and the objective senses are real; the population loop that produces a
genuinely multi-objective front is Phase 2.

**Some tiers your budget names.** A budget tier may plan a stage this build
does not execute — an `exact` budget plans a `decision` stage, which is Phase
2's. The job reports it as `skipped` and warns with `CF-ORCH-006` naming the
stage, rather than reporting success over less search than the plan describes.
If a stage the budget planned is missing from `progress.tiers` as `done` and
there is no `CF-ORCH-006` warning for it, that is a bug worth reporting.

**Most objective terms.** This build ships `usage.cai`,
`constraint.violation_count` and `structure.local_accessibility` as rankable
objectives, plus `usage.gc_content` as a **diagnostic**. Requesting anything
else is refused with `CF-LAND-001` rather than silently optimizing a smaller
objective set — a job that asked for four objectives must not quietly optimize
three.

`usage.gc_content` is reported and never ranked on, and that is deliberate. GC
content has no direction: too high and the transcript folds, too low and it
does the same for other reasons. Its provider caps itself at `diagnostic` for
that reason, and a job requesting it as an objective gets `CF-LAND-014` saying
so. Bound GC as a *constraint* instead, with `Rule::GcWindow`.

`structure.local_accessibility` needs a user-installed ViennaRNA
(`pip install ViennaRNA`); CodonForge ships no folding engine (`AMD-013`). It is
scored out of core, once per generation on the whole population, so tiers 1 and
2 do not see it — `CF-SRCH-010` says so on every job that asks for it.

**A REST server.** `codonforge.api.serve` serves the job routes. Hosts, models,
observations, presets, protein input, construct input, and the five-check
feasibility route answer `CF-CLI-003`.

---

## Exit codes

A script needs to tell an infeasible constraint set from a crash.

| Code | Meaning |
| --- | --- |
| `0` | Success |
| `1` | Design failed |
| `2` | Invalid input |
| `3` | Infeasible constraints |
| `4` | Quota or budget |
| `5` | Internal error |

Every error prints the same `CF-*` code the API returns.

---

## Development

```sh
make setup     # virtualenv and Python dev dependencies
make build     # compile the Rust core into the virtualenv
make test      # cargo test + pytest
make lint      # clippy + ruff
```

The planning package in [`planning/`](planning/) is normative: schemas,
algorithm specifications, ADRs, the test plan and the build sequence. It is
mechanically validated —

```sh
.venv/bin/python planning/tools/validate_planning_package.py
python3 tools/check_module_layering.py --check
python3 tools/check_schema_sync.py --check
```

See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the two-toolchain setup, the
module layering rules and what a change is expected to carry.
