# DNA Design Suite

DNA Design Suite is one local page with two tabs.

**CodonForge** rewrites a coding sequence for a named host. The protein stays
the same. Objectives and hard sequence constraints are explicit.

**PrimerForge** designs PCR primer pairs for a region on a DNA template, and
checks a pair you already have.

The search runs in Rust. The job lifecycle, ingest, and command line run in
Python. CodonForge and PrimerForge share the page and the visual language.

Licensed under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0).
The text is [`LICENSE`](LICENSE). `pyproject.toml` and `Cargo.toml` declare
the same identifier. ViennaRNA is optional, installed by you, and is not
part of this licence.

Validation for v2.5.0 through v2.9.1 is in [`docs/benchmarks/`](docs/benchmarks/).
The Software-wiki draft, including one validation table, is
[`docs/wiki-software-draft.md`](docs/wiki-software-draft.md). Current status
is [`docs/project-status.md`](docs/project-status.md). History is
[`CHANGELOG.md`](CHANGELOG.md).

---

## Quick start

The page listens on <http://127.0.0.1:8000/>, opens a browser, and waits
until Ctrl+C. `--no-browser` skips the window. `--port` chooses another port.

A wheel is not on PyPI for this tag. `pip install codonforge` does not
resolve. Install from this checkout. Python 3.11 or newer, and a Rust
toolchain (`rustup`), are required because the core is a compiled extension.

### Windows (PowerShell, Anaconda)

From the repository directory:

```powershell
conda create -n dna-design python=3.11 -y
conda activate dna-design
pip install "maturin>=1.5,<2"
maturin develop
python -m codonforge.api
```

After that environment exists, the command is `python -m codonforge.api`.

### macOS and Linux

```sh
make setup
make build
.venv/bin/python -m codonforge.api
```

`make setup` uses `virtualenv` and falls back to `python3 -m venv` only if
that fails. The reason is in [`CONTRIBUTING.md`](CONTRIBUTING.md): Debian
and Ubuntu ship `ensurepip` in a separate `python3-venv` package. The same
module command works in any environment the package was installed into:

```sh
python -m codonforge.api
```

### First run

CodonForge is the tab that opens. **Load the rpsL demo**, then **Run
design**. The page loads the bundled *E. coli* K-12 MG1655 ribosomal usage
table on its own. A finished result can be copied or saved as FASTA,
GenBank, or CSV.

PrimerForge is the other tab. **Load the rpsL sequence** fills the gene,
selects bases 41–300, and sets the flank to 40 nt. Then **Find primer
pairs**.

`rpsL` is the *E. coli* K-12 MG1655 CDS `NP_417801.1` from `NC_000913.3`,
125 codons. The usage table is the relative adaptiveness counted from the
ribosomal-protein CDS in that same file (`rpl`/`rps`/`rpm`), labelled
`user_supplied`.

The same demo from the command line, after the install above:

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

On Windows PowerShell, pass the CDS with
`(Get-Content -Raw examples\rpsL.cds) -replace '\s',''` in place of the
`tr` substitution.

---

## Honest limits

The phase plan is
[`planning/06-build/00-sequencing-overview.md`](planning/06-build/00-sequencing-overview.md).
This is what the current tree runs. The refusals below are named errors,
recorded in [`CHANGELOG.md`](CHANGELOG.md).

| Phase | What this suite runs |
| --- | --- |
| 1 | Hard constraints and a CAI-optimized synonymous CDS. Tier 1 is dynamic programming for an additive objective. Tier 2 is block coordinate ascent. When that program cannot carry a GC window or a scattered repeat, a deterministic seed repair runs. A job that still cannot is exit 3, `CF-SRCH-006`. |
| 2 | `structure.local_accessibility` (ViennaRNA, scored in tier 3) and `burden.growth_cost`. Two or more complete axes recommend the knee of that front. |
| 3 | `epistasis.context_naturalness`, the selected *E. coli* synonymous Potts model, window 10. |
| 4 | `design_operon`. Each CDS is designed, then a synonymous swap is kept only when the caller-supplied flux-weighted burden falls. |
| 5 | `POST /v1/observations` stores a validated measurement. When a line is identified, it is ordinary least squares of one uncensored numeric readout on the GC fraction of an explicit CDS. |

Stays refused:

- **Dwell.** `epistasis.dwell_deviation` is `CF-LAND-001`. The dwell head is
  `not_fitted`. No 5′-end offset for GEO GSE53767 was fitted, and those files
  were not reachable. Changelog: v2.4.0, v2.7.0, v2.9.1.
- **Host catalogue and recalibration.** `GET /v1/hosts`, `GET /v1/models`,
  `GET /v1/recalibrations`, seed import, protein input, construct input, and
  the five-check feasibility route answer `CF-CLI-003`. An accepted
  observation is stored. It does not start a recalibration. Changelog: v2.9.0,
  v2.9.1.
- **No invented supply, s-value, or distance.** Burden uses
  `trna.gene_copy_number` when that map is non-empty, otherwise the codon-usage
  table. A profile with neither is `CF-BURD-004`. The caller supplies every
  operon flux. A Potts request for a non-*E. coli* host, or a missing MG1655
  artifact, is `CF-LAND-001`. No host-distance penalty is applied. Changelog:
  v2.6.0, v2.7.0, v2.8.0, v2.9.1.
- **`GATE-BEN-001` is CAI-only.** 50 of 50 *E. coli* K-12 proteins, budget
  `fast`, 4 generations, seed 7, baseline `cai_greedy`. CodonForge meets or
  equals that baseline on all 50. The v2.4.2 record was 48 of 50; v2.9.1
  brings `aceB` and `thrA` to 1.0 with no in-frame `AGG`. The claim is not an
  MFE comparison and not a comparison to LinearDesign or CDSfold. Changelog:
  v2.4.2, v2.9.1. The row is in `benchmarks/corpus/gates.json`.

Also held, and named as such in [`docs/project-status.md`](docs/project-status.md):
`P1-26`, LinearDesign's method, bundled thermodynamic tables, EXP-06, and
`P3-14` (language-model scoring). EXP-05 did not separate Lamarckian from
Baldwinian repair; Lamarckian stays the default. PrimerForge does not check
specificity against a genome, does not wrap Primer3, and does not use mFold
for hairpins. Commercial distribution is outside the AMD-013 note in
[`docs/legal-clearance.yaml`](docs/legal-clearance.yaml): a non-commercial
iGEM educational project, and the reviewer states that the email is not a
formal legal opinion.

The suite note is v2.9.1. The package version stays 0.1.0.
`python -m codonforge.cli.main --version` prints 0.1.0.

---

## What it does

A design job runs in stages. A single-objective job reports the best candidate
any stage produced, not the last stage's output. A job with two or more
complete axes reports the knee of that front (`knee_of_complete_front`).
`candidates[0]` is that recommendation either way.

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

`design_operon` designs each coding sequence on its own, then keeps a
synonymous swap only when the shared flux-weighted codon demand scores a
lower `burden.growth_cost`. The caller supplies each flux.

Benchmark write-ups, in order:

| Release | File |
| --- | --- |
| v2.5.0 performance | [`docs/benchmarks/v2.5.0-performance.md`](docs/benchmarks/v2.5.0-performance.md) |
| v2.5.0 validation | [`docs/benchmarks/v2.5.0-validation.md`](docs/benchmarks/v2.5.0-validation.md) |
| v2.6.0 structure and burden | [`docs/benchmarks/v2.6.0-phase2.md`](docs/benchmarks/v2.6.0-phase2.md) |
| v2.9.1 gate, folds, vectors | [`docs/benchmarks/v2.9.1-validation.md`](docs/benchmarks/v2.9.1-validation.md) |
| v2.9.1 stress timings | [`docs/benchmarks/v2.9.1-challenging-perf.md`](docs/benchmarks/v2.9.1-challenging-perf.md) |
| PrimerForge thermodynamics | [`docs/primerforge-validation.md`](docs/primerforge-validation.md) |

---

## What does not work yet

The [honest limits](#honest-limits) section is the map. Each row here is the
error you get if you call the path.

**Background execution.** `cf.design(..., background=True)` raises
`CF-CLI-003`. It needs a worker pool that does not exist. It does not silently
run synchronously, because that would make every timing assumption in your code
wrong.

**Some tiers your budget names.** A budget tier may plan a stage this build
does not execute. An `exact` budget plans a `decision` stage. The job reports
it as `skipped` and warns with `CF-ORCH-006` naming the
stage, rather than reporting success over less search than the plan describes.
If a stage the budget planned is missing from `progress.tiers` as `done` and
there is no `CF-ORCH-006` warning for it, that is a bug worth reporting.

**Most objective terms.** This build ships `usage.cai`,
`constraint.violation_count`, `structure.local_accessibility`,
`burden.growth_cost`, and `epistasis.context_naturalness` as rankable
objectives, plus `usage.gc_content` as a **diagnostic**. The context term is
ranked only when the verified E. coli Potts artifact is installed. A missing
artifact, or a host that is not *E. coli*, is `CF-LAND-001` and is not scored
as zero. Requesting anything else is refused with `CF-LAND-001` rather than
silently optimizing a smaller objective set — a job that asked for four
objectives must not quietly optimize three.

`burden.growth_cost` is a tier-3 demand-profile mismatch in `[0, 1)`, unit
`dimensionless`, minimized. It is not a predicted growth rate. Supply is the
host's `trna.gene_copy_number` when that map is non-empty, otherwise the
synonymous shares of `codon_usage.table`. A host with neither is `CF-BURD-004`.

`usage.gc_content` is reported and never ranked on, and that is deliberate. GC
content has no direction: too high and the transcript folds, too low and it
does the same for other reasons. Its provider caps itself at `diagnostic` for
that reason, and a job requesting it as an objective gets `CF-LAND-014` saying
so. Bound GC as a *constraint* instead, with `Rule::GcWindow`.

`structure.local_accessibility` needs a user-installed ViennaRNA
(`pip install ViennaRNA`); CodonForge ships no folding engine (`AMD-013`). It is
scored out of core, once per generation on the whole population, so tiers 1 and
2 do not see it — `CF-SRCH-010` says so on every job that asks for it.

**A REST server.** `codonforge.api.serve` serves the job routes and
`/v1/observations`. Hosts, models, recalibrations, presets, protein input,
construct input, and the five-check feasibility route answer `CF-CLI-003`.
An accepted observation is stored. It does not start a recalibration.

**Shadow prices and a joint NSGA-III pass.** `design_operon` does not run
them. The shared ledger is the flux-weighted codon demand.

---

## Design a sequence

You need a host codon-usage profile and, optionally, a constraint set. Both are
JSON documents validated against the schemas in [`schemas/v1`](schemas/v1).

```sh
.venv/bin/python -m codonforge.cli.main design \
    --host escherichia_coli:511145/BL21 \
    --cds ATGCTGCGTAGCGGCGTGGCGCCGACCAAACTGCGTAGCTAA \
    --host-profile examples/mg1655_ribosomal_usage.json \
    --constraints examples/rpsL_constraints.json \
    --objectives usage.cai \
    --seed 20260828 \
    --output explain
```

On this tree that command succeeds. The explain header is:

```
objectives
  usage.cai        1.000000  maximize

tiers
  tier1_seed      done
  tier2_block_ca  done
  tier3_moea      done  9/12

constraints  satisfied
```

`usage.cai` of 1.0 here means every counted codon is a top synonym of the
loaded ribosomal table. The same CDS scores 0.733847 on the Kazusa genome-wide
weights, which did not choose the codons. The same run warns `CF-SRCH-001`
(the Potts-guided and block-Gibbs operators are absent, so the mix is
renormalized) and `CF-LAND-013` (`structure.local_accessibility` is in the
report and was not an objective). Adding `usage.gc_content` to `--objectives`
still succeeds and warns `CF-LAND-014`: that term is reported, not ranked.

Other formats: `--output json` (the schema-validated documents, verbatim),
`fasta`, `genbank`, `csv`, `table`.

### Or from Python

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
make check     # the invariants, the linters, and the tests
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
