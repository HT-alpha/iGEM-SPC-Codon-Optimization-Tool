# DNA Design Suite

DNA Design Suite is one local page with two tabs.

The codon optimization tool rewrites a coding sequence for a named host. The
protein stays the same. Objectives and hard sequence constraints are explicit.

The primer design tool picks PCR primer pairs for a region on a DNA template,
and checks a pair you already have.

The search runs in Rust. The job lifecycle, ingest, and command line run in
Python. Both tools share the page and the visual language.

Licensed under the [Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0).
The text is [`LICENSE`](LICENSE). `pyproject.toml` and `Cargo.toml` declare
the same identifier. ViennaRNA is optional, installed by you, and is not
part of this licence.

Validation notes, including v2.5.0 through v2.9.2, are in [`docs/benchmarks/`](docs/benchmarks/).
The Software-wiki draft, including one validation table, is
[`docs/wiki-software-draft.md`](docs/wiki-software-draft.md). Current status
is [`docs/project-status.md`](docs/project-status.md). History is
[`CHANGELOG.md`](CHANGELOG.md).

---

## Quick start

The page listens on <http://127.0.0.1:8000/>, opens a browser, and waits
until Ctrl+C. `--no-browser` skips the window. `--port` chooses another port.

`pip install dna-design-suite` does not resolve: the package is not on PyPI,
there is no release tag, and this repository stays private. Version 0.1.0
is the package version. Python 3.11 or newer is required. One wheel per
operating system covers 3.11 through 3.13 (`abi3-py311`) when a wheel file
exists. CPython 3.9 and 3.10 are not supported.

The install path from a checkout of this repository is source, below.
GitHub Actions is not an install path: there is no release, and this README
does not send you to a workflow artifact. A wheel file is something you
already have on disk. This README does not claim a wheel is attached to a
release, and it does not claim the package is published.

### Install a wheel you already have

Installing one is pip against that file. It is not `pip install dna-design-suite`
from an index, and it is not a download from a GitHub Actions run.

Linux x86_64. `virtualenv` is the bootstrap `make setup` uses. `python3 -m venv` is only the fallback when virtualenv itself cannot be installed, because that module needs `ensurepip`.

```sh
python3 -m virtualenv .venv
.venv/bin/pip install dna-design-suite-0.1.0-cp311-abi3-manylinux_*.whl
.venv/bin/python -m dna_design_suite.api
```

Windows (PowerShell, Anaconda), from the directory that contains the wheel:

```powershell
conda create -n dna-design python=3.11 -y
conda activate dna-design
pip install .\dna-design-suite-0.1.0-cp311-abi3-win_amd64.whl
python -m dna_design_suite.api
```

The filename ends in `win_amd64` on 64-bit Windows. ViennaRNA is still
optional and is not inside the wheel.

### Install from source

Use this when no wheel is available. A Rust toolchain (`rustup`) is
required, because the core is a compiled extension.

Windows (PowerShell, Anaconda), from the repository directory:

```powershell
conda create -n dna-design python=3.11 -y
conda activate dna-design
pip install "maturin>=1.5,<2"
maturin develop
python -m dna_design_suite.api
```

After that environment exists, the command is `python -m dna_design_suite.api`.

macOS and Linux:

```sh
make setup
make build
.venv/bin/python -m dna_design_suite.api
```

`make setup` uses `virtualenv` and falls back to `python3 -m venv` only if
that fails. The reason is in [`CONTRIBUTING.md`](CONTRIBUTING.md): Debian
and Ubuntu ship `ensurepip` in a separate `python3-venv` package. The same
module command works in any environment the package was installed into:

```sh
python -m dna_design_suite.api
```

### Publishing a release

Not done. The package is not on PyPI. There is no release tag. Publishing
is Hugo's step, from a machine that already has a wheel file, a PyPI
account, and a token he supplies: `pip install twine` and `twine upload`
of that file. Do not add the token to the repository. Until that upload,
install from the source checkout.

### First run

DNA Design Suite is the tab that opens. **Load the rpsL demo**, then **Run
design**. The page loads the bundled *E. coli* K-12 MG1655 ribosomal usage
table on its own. A finished result can be copied or saved as FASTA,
GenBank, or CSV.

Primer design is the other tab. **Load the rpsL sequence** fills the gene,
selects bases 41–300, and sets the flank to 40 nt. Then **Find primer
pairs**.

`rpsL` is the *E. coli* K-12 MG1655 CDS `NP_417801.1` from `NC_000913.3`,
125 codons. The usage table is the relative adaptiveness counted from the
ribosomal-protein CDS in that same file (`rpl`/`rps`/`rpm`), labelled
`user_supplied`.

The same demo from the command line, after the install above:

```sh
python -m dna_design_suite.cli.main design \
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
| 2 | `structure.local_accessibility` (ViennaRNA). A construct that fits in one 40 nt span is folded in tiers 1 and 2. A longer gene stays on the tier 3 batch, and that job still warns `CF-SRCH-010`. `burden.growth_cost` is the demand mismatch. Two or more complete axes recommend the knee of that front. |
| 3 | `epistasis.context_naturalness`, the selected *E. coli* synonymous Potts model, window 10, when that file is installed and verifies. |
| 4 | `design_operon`. Each CDS is designed, then a synonymous swap is kept only when the caller-supplied flux-weighted burden falls. Shadow prices are finite differences of that same cost. A joint NSGA-III pass adopts a sequence only when the joint cost falls. |
| 5 | `POST /v1/observations` stores a validated measurement. When the Phase 5 ordinary least squares fit of one uncensored numeric readout on the GC fraction of an explicit CDS is `fitted`, a recalibration record is stored. `not_fitted` is not scheduled. |

Stays refused:

- **Dwell.** `epistasis.dwell_deviation` is `CF-LAND-001`. The dwell head is
  `not_fitted`. No paper states one 5′-end A-site offset for GEO GSE53767
  (Li et al., Cell 2014, center-weight reads of 20–42 nt; the GSE53767 files
  were not reachable). A zero deviation is not substituted. Changelog:
  v2.4.0, v2.7.0, v2.9.1.
- **Hosts.** `GET /v1/hosts` is the example codon-usage files on disk, not an
  organism catalogue. `POST /v1/hosts:import` reads Kazusa text, a delimited
  codon table, or a CDS set and stores `user_supplied` relative adaptiveness
  (count divided by the family maximum; no pseudocount). `genome_wide` is
  refused because `usage.cai` will not score that subset (`CF-HOST-003`).
  Models are the verified Potts file, or an empty list with `installed: false`.
- **No invented supply, s-value, or distance.** Burden uses
  `trna.gene_copy_number` when that map is non-empty, otherwise the codon-usage
  table. A profile with neither is `CF-BURD-004`. The caller supplies every
  operon flux. A Potts request for a non-*E. coli* host, or a missing MG1655
  artifact, is `CF-LAND-001`. No host-distance penalty is applied. Changelog:
  v2.6.0, v2.7.0, v2.8.0, v2.9.1.
- **rpsL assembly CAI.** Budget `fast`, generations 0, seed 7, the ribosomal
  table, and `benchmarks/validation/constraints_assembly.json` stay at
  0.9866866273575934 (SHA-256
  `0624d04f50821da0804a0cf91fa59b9a6bde2ba8ef6156bb79e9680ba09eedfc`).
  Codons 74 (`CAA`) and 88 (`GAT`) are the only counted codons off their
  heaviest synonym. `CAG` creates one PstI site, and the cheapest Leu codon
  that removes it costs more CAI than `CAG` returns. `GAC` creates one
  repeated 12-mer, and every feasible repair of that repeat costs more than
  `GAC` returns. The published difference 0.0038104428470379093 is to that
  infeasible sequence. No search change in this tree chases it.
- **`GATE-BEN-001` is CAI-only.** 50 of 50 *E. coli* K-12 proteins, budget
  `fast`, 4 generations, seed 7, baseline `cai_greedy`. DNA Design Suite meets or
  equals that baseline on all 50. The v2.4.2 record was 48 of 50; v2.9.1
  brings `aceB` and `thrA` to 1.0 with no in-frame `AGG`. The claim is not an
  MFE comparison and not a comparison to LinearDesign or CDSfold. Changelog:
  v2.4.2, v2.9.1. The row is in `benchmarks/corpus/gates.json`.

Also held, and named as such in [`docs/project-status.md`](docs/project-status.md):
`P1-26`, LinearDesign's method, bundled thermodynamic tables, EXP-06, and
`P3-14` (language-model scoring). EXP-05 did not separate Lamarckian from
Baldwinian repair; Lamarckian stays the default. Primer design does not check
specificity against a genome, does not wrap Primer3, and does not use mFold
for hairpins. Commercial distribution is outside the AMD-013 note in
[`docs/legal-clearance.yaml`](docs/legal-clearance.yaml): a non-commercial
iGEM educational project, and the reviewer states that the email is not a
formal legal opinion.

The status note is the current tip. The package version is 0.1.0.
There is no release tag, and the package is not on PyPI.
`python -m dna_design_suite.cli.main --version` prints 0.1.0.

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
| v2.9.2 profiles | [`docs/benchmarks/v2.9.2-profiles.md`](docs/benchmarks/v2.9.2-profiles.md) |
| Primer design thermodynamics | [`docs/primer-thermo-validation.md`](docs/primer-thermo-validation.md) |
| Commercial codon tools | [`docs/benchmarks/commercial-codon-tools.md`](docs/benchmarks/commercial-codon-tools.md) |

---

## What does not work yet

The [honest limits](#honest-limits) section is the map. Each row here is the
error you get if you call the path.

**Background execution.** `design(..., background=True)` returns a handle
with `job_id`, `done()`, and `result(timeout)`. The job runs on a worker
pool. The call does not silently run synchronously: a handle whose `done()`
is already true before the job has started would make every timing
assumption in your code wrong. `CF-CLI-003` remains the answer for an
unknown route. Seed import is `POST /v1/seeds:import` and
`seeds import`.

**The decision stage.** An `exact` budget includes `tier1_seed`,
`tier2_block_ca`, `tier3_moea`, and `decision`. The decision stage records
the same recommendation the front already names and does not replace the
candidate. Its progress row is `done`. A tier this build does not run is
still reported `skipped` with `CF-ORCH-006`. If a stage the budget planned
is missing from `progress.tiers` as `done` and there is no `CF-ORCH-006`
warning for it, that is a bug worth reporting.

**Most objective terms.** This build ships `usage.cai`,
`usage.harmonize_rca`, `usage.match_frequency`,
`constraint.violation_count`, `structure.local_accessibility`,
`burden.growth_cost`, and `epistasis.context_naturalness` as rankable
objectives, plus `usage.gc_content` as a **diagnostic**. The context term is
ranked only when the verified E. coli Potts artifact is installed. A missing
artifact, or a host that is not *E. coli*, is `CF-LAND-001` and is not scored
as zero. Requesting anything else is refused with `CF-LAND-001` rather than
silently optimizing a smaller objective set — a job that asked for four
objectives must not quietly optimize three.

`usage.harmonize_rca` minimizes the mean absolute difference between each
input codon's relative adaptiveness in a source table and the chosen codon's
relative adaptiveness in the expression-host table. The method is Claassens
et al., PLoS ONE 2017 (doi:10.1371/journal.pone.0184355), the same per-codon
replacement Angov et al., PLoS ONE 2008 (doi:10.1371/journal.pone.0002189)
describe. Pass the source table as `host_profile.source_codon_usage`, or on
the command line as `--source-profile`. The input CDS is pinned before the
job id is derived. A missing source table is `CF-HOST-003`. The suite does
not invent one.

`usage.match_frequency` minimizes the squared difference between the gene's
synonymous codon fractions and the host table's fractions. One codon changes
every synonym of that amino acid, so tiers 1 and 2 omit it. A fast budget
warns `CF-SRCH-010` and still reports the score. Standard and exact rank it
in the population loop. It is not a bit-identical port of another tool's
frequency matcher.

`burden.growth_cost` is a tier-3 demand-profile mismatch in `[0, 1)`, unit
`dimensionless`, minimized. It is not a predicted growth rate. Supply is the
host's `trna.gene_copy_number` when that map is non-empty, otherwise the
synonymous shares of `codon_usage.table`. A host with neither is `CF-BURD-004`.

`usage.gc_content` is reported and never ranked on, and that is deliberate. GC
content has no direction: too high and the transcript folds, too low and it
does the same for other reasons. Its provider caps itself at `diagnostic` for
that reason, and a job requesting it as an objective gets `CF-LAND-014` saying
so. Bound GC as a *constraint* instead, with `Rule::GcWindow`.

`structure.local_accessibility` uses a user-installed ViennaRNA
(`pip install ViennaRNA`) unless the request sets `structure_engine` to
`nussinov`. DNA Design Suite ships no thermodynamic engine (`AMD-013`).
Nussinov is an opt-in maximum matching (pair score 1 for A–U, G–C and G–U,
hairpin of at least 3; Nussinov and Jacobson, PNAS 1980). The value is the
mean paired fraction in [0, 1]. It is not a free energy and not an unpaired
probability, and it is not selected when ViennaRNA is missing. A construct
longer than one 40 nt span is scored out of core, once per generation, and
`CF-SRCH-010` says tiers 1 and 2 did not fold it. A construct that fits in
that span is folded during those tiers. There is no claim of parity with a
global minimum free energy or with LinearDesign.

A separate report, `dna_design_suite fold` and `POST /v1/structure:global`,
labels a whole-sequence fold. When ViennaRNA is installed the value is
`RNA.fold`'s global minimum free energy in kcal/mol. When it is not, a
sequence of at most 400 nt gets a Watson-Crick pair count, and the label
says the unit is base pairs. A longer sequence without ViennaRNA is left
unscored. That report is not `structure.local_accessibility` and not the
`nussinov` objective.

A CAI-only job does not fold. The old diagnostic fold is not run.

**The REST server.** `dna_design_suite.api.serve` serves jobs, observations,
the example host files, a verified model file, recalibrations stored from a
fitted observation, objective presets (`cai`, `cai_burden`, `cai_structure`,
`cai_context`, `harmonize`, `match_frequency`), protein input
reverse-translated with the caller's usage table, construct input read from
the GenBank CDS features in the body, and `POST /v1/feasibility` (five
checks, plus a constructive witness when one exists). An empty feasibility
body is `400`. `POST /v1/seeds:import` checks a coding sequence against the
protein. `POST /v1/hosts:import` builds a usage table. `POST /v1/measurements`
stores a lab number and `GET /v1/observations?format=csv` exports it. An
unknown route is `CF-CLI-003`.

**Host tables.** There is no species catalogue. `GET /v1/hosts` lists the
example files installed with the page. `POST /v1/hosts:import`, and
`dna_design_suite hosts import`, accept Kazusa text, a delimited table, a
JSON `codon_usage` block, or a FASTA of coding sequences. Relative
adaptiveness is each value divided by the maximum in its synonymous family.
The stored label stays `user_supplied` unless the caller sets
`highly_expressed`. `genome_wide` is refused (`CF-HOST-003`). Asserting
`highly_expressed` is not a check that the file is a highly expressed
reference.

**Vendor-class screen.** `dna_design_suite vendor-screen` and
`POST /v1/vendor-screen` flag GC windows, homopolymers, repeated 12-mers,
and the assembly restriction sites and `AGGAGG` motif, using the numbers in
this suite's synthesis and assembly presets. The checklist is
[`docs/vendor-screen.md`](docs/vendor-screen.md). `vendor_score` is false.
Wet-lab expression is out of scope. A gene still has to pass the
manufacturer's own screen.

```sh
.venv/bin/python -m dna_design_suite.cli.main hosts import \
    --input benchmarks/validation/kazusa/e_coli_316407.csv
.venv/bin/python -m dna_design_suite.cli.main vendor-screen --cds ATGAAATAA
.venv/bin/python -m dna_design_suite.cli.main fold --cds GGGAAACCC
.venv/bin/python -m dna_design_suite.cli.main design \
    --host escherichia_coli:511145/MG1655 \
    --cds ATGCTGCGTAGCGGCGTGGCGCCGACCAAACTGCGTAGCTAA \
    --host-profile examples/mg1655_ribosomal_usage.json \
    --source-profile examples/mg1655_ribosomal_usage.json \
    --objectives usage.harmonize_rca \
    --budget fast \
    --seed 7
```

**Operon prices.** `design_operon` reports 64 finite differences of
`burden.growth_cost`. A slot the caller's table cannot score is null. The
joint pass does not invent expression rates.

---

## Design a sequence

You need a host codon-usage profile and, optionally, a constraint set. Both are
JSON documents validated against the schemas in [`schemas/v1`](schemas/v1).

```sh
.venv/bin/python -m dna_design_suite.cli.main design \
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
renormalized). A CAI-only job does not fold, so it does not warn
`CF-LAND-013`. Adding `usage.gc_content` to `--objectives`
still succeeds and warns `CF-LAND-014`: that term is reported, not ranked.

Other formats: `--output json` (the schema-validated documents, verbatim),
`fasta`, `genbank`, `csv`, `table`.

### Or from Python

```python
import dna_design_suite as cf
from dna_design_suite.modelreg import ArtifactStore
from dna_design_suite.codon.orchestrator import DesignRequest

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
