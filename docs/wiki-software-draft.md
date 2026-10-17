# DNA Design Suite — Software wiki draft

Paste into the team Software page. Replace the two placeholders before
publishing. Upload the screenshots named below, or follow the checklist if
the files are not in the checkout yet.

- Public repository: `{{PUBLIC_REPOSITORY_URL}}`
- Tag `v2.9.1`: `{{V2_9_1_TAG_URL}}`

The repository this draft was written from is still private. Publishing the
wiki does not publish the git repository.

## Purpose

DNA Design Suite is a local tool for designing the DNA of a coding sequence
and the PCR primers that amplify it.

CodonForge keeps the protein the same and chooses synonymous codons for a
host you name. You can ask it to raise codon adaptation, to respect hard
sequence rules (restriction sites, GC windows, repeats, homopolymers), and,
when you ask, to rank local mRNA accessibility and a host-burden mismatch.
A job with two or more scored axes returns the knee of that front, and the
page shows the range of each axis. Every job seals a lockfile so the same
request and seed can be run again.

PrimerForge is the other tab on the same page. It designs primer pairs for
a region you select, and it checks a pair you paste. Melting temperature
and dimer ΔG use SantaLucia 1998 Table 2 in a two-state calculation written
for this project. It does not search a genome for off-targets.

The search core is Rust. The page, the command line, and the job files are
Python. The licence is Apache-2.0.

## Install

Python 3.11 or newer, and Rust (`rustup`). A PyPI wheel is not published
for v2.9.1. Clone `{{PUBLIC_REPOSITORY_URL}}` and build the extension.

Windows (PowerShell, Anaconda), from the repository directory:

```powershell
conda create -n dna-design python=3.11 -y
conda activate dna-design
pip install "maturin>=1.5,<2"
maturin develop
python -m codonforge.api
```

macOS and Linux:

```sh
make setup
make build
.venv/bin/python -m codonforge.api
```

The page opens at <http://127.0.0.1:8000/> and runs until Ctrl+C.
CodonForge is the first tab. PrimerForge is the second. **Load the rpsL
demo** fills the *E. coli* K-12 MG1655 gene `rpsL` (`NP_417801.1`, 125
codons), the ribosomal codon table, and the example constraint file.

ViennaRNA is optional. CodonForge does not vendor it. Structure scoring
needs `pip install ViennaRNA` (the v2.9.1 timings used 2.7.2). A CAI-only
job runs without it. A structure request with no engine installed is
`CF-STRU-003`.

Command-line design, primer routes, and the Python API are in the
repository `README.md`.

## What the five phases do here

The plan is `planning/06-build/00-sequencing-overview.md`. The tree at
`v2.9.1` runs the following. Details and the changelog citations are in
the README section **Honest limits**.

| Phase | In this suite |
| --- | --- |
| 1 | Constrained CAI design: tier 1 dynamic programming, tier 2 block ascent, and a deterministic seed repair when the program cannot carry a GC window or a scattered repeat. A job that still cannot is exit 3, `CF-SRCH-006`. |
| 2 | Local mRNA accessibility and `burden.growth_cost` on a Pareto front. The recommendation for two or more complete axes is the knee. |
| 3 | `epistasis.context_naturalness` from the selected *E. coli* Potts model (window 10, held-out naturalness −0.8708682912309792 nats/codon). |
| 4 | `design_operon`: a synonymous swap is kept only when the caller-supplied flux-weighted burden falls. |
| 5 | `POST /v1/observations` stores a validated measurement. The secondary model is a least-squares line of one uncensored numeric readout on CDS GC fraction, or `not_fitted`. |

## What this version refuses

- Dwell (`epistasis.dwell_deviation`) is `CF-LAND-001`. The head is
  `not_fitted`. Changelog v2.4.0, v2.7.0, v2.9.1.
- The host catalogue and `GET /v1/recalibrations` are `CF-CLI-003`. Storing
  an observation does not refit a model. Changelog v2.9.0.
- Gene-copy numbers, s-values, and host-distance penalties are not invented.
  Burden supply is the profile's `trna.gene_copy_number` when present,
  otherwise the codon table. Operon fluxes are the caller's. Changelog
  v2.6.0, v2.7.0, v2.8.0.
- `GATE-BEN-001` is the CAI-only comparison against `cai_greedy` on 50 of
  50 proteins. It is not an MFE comparison and not LinearDesign or CDSfold.
  Changelog v2.4.2 (48 of 50) and v2.9.1 (50 of 50).
- PrimerForge does not check genome specificity, does not wrap Primer3, and
  does not score hairpins with mFold.
- Commercial distribution is outside `docs/legal-clearance.yaml`. The
  recorded review is a non-commercial iGEM educational project, and the
  reviewer states that it is not a formal legal opinion.

## Validation

One table. Every figure is copied from the file in **Source**. Budgets
differ by row; they are in **Setting**.

| Check | Setting | Result | Source |
| --- | --- | --- | --- |
| CAI gate `GATE-BEN-001` | 50 *E. coli* K-12 CDS; budget `fast`, 4 generations, seed 7; baseline `cai_greedy`; empty constraints | Claimed. CodonForge meets or equals `cai_greedy` on 50 of 50. `codonforge_ge_greedy: 50`, `behind_greedy: []`. CAI-only. The MFE comparison is not this claim. | `benchmarks/corpus/gates.json`; `CHANGELOG.md` v2.9.1 |
| `aceB` / `thrA` on that gate | Same gate budget. Native `AGG` (weight 0 in the ribosomal table) | Both score CAI 1.0 with no in-frame `AGG`. The previous record scored 0. The other 48 `codonforge_cai` values are unchanged. | `docs/benchmarks/v2.9.1-validation.md` |
| Constrained yagF | Assembly rules; CAI only; budget `fast`, generations 0, seed 7 | 1968 nt, CAI 0.9286425884184984, 0 violations, same protein, 0.602 s | `docs/benchmarks/v2.9.1-validation.md` |
| Constrained waaS | Same | 936 nt, CAI 0.9173459052341191, 0 violations, same protein, 0.194 s | `docs/benchmarks/v2.9.1-validation.md` |
| Constrained elastin-like VPGVG×20 | Same | 306 nt, CAI 0.43802555520567116, 0 violations, same protein, 0.144 s | `docs/benchmarks/v2.9.1-validation.md` |
| Constrained yeeJ | Same | 7077 nt, CAI 0.9580065306093907, 0 violations, same protein, 21.386 s | `docs/benchmarks/v2.9.1-validation.md` |
| PrimerForge, Table 2 duplex | `5'-CGTTGA-3'` and its complement; SantaLucia 1998 Table 2 | ΔH° −41.2 kcal/mol, ΔS° −115.4 cal/K·mol, ΔG°37 −5.409 kcal/mol, Tm 9.805 °C. Engine matches the hand sum. | `docs/primerforge-validation.md` |
| PrimerForge, published self-dimer | `CAGCCACTGATGCCAAGAGA` | ΔG°37 −3.556 kcal/mol. Published OligoAnalyzer cell −3.55 kcal/mol. Difference 0.006 kcal/mol. | `docs/primerforge-validation.md` |
| PrimerForge, yeeJ CDS flanks | 300 nt of the CDS on each side | 20 pairs, product 6843, ΔG −1.28385, not refused, 0.158 s | `docs/benchmarks/v2.9.1-challenging-perf.md` |
| PrimerForge, synthetic flanks | yeeJ; `ACGT` repeated, 300 nt each side | 20 pairs, product 7395, ΔG −20.754845, refused (dimer below −9 kcal/mol) | `docs/benchmarks/v2.9.1-challenging-perf.md` |
| Fold cache, rpsL | CAI + structure, assembly, budget `standard`, 4 generations, seed 7 | Vector unchanged `[0.9731171656422364, 2.581242176257437]`. Wall 1.357 s then 0.561 s. | `docs/benchmarks/v2.9.1-validation.md` |
| Fold cache, waaS | Same | Vector unchanged `[0.9086209272414913, 1.775449070836324]`. Wall 4.211 s then 2.255 s. | `docs/benchmarks/v2.9.1-validation.md` |
| Operon lacZ + thrA | Fluxes 1 and 2; burden; budget `fast`, generations 0, seed 7 | Joint cost 0.4452346719259927, swaps `[26, 0]`. Wall 9.185 s then 1.632 s. | `CHANGELOG.md` v2.9.1; `docs/benchmarks/v2.9.1-challenging-perf.md` |
| Stress, yeeJ standard | CAI, assembly, budget `standard`, 4 generations | Wall 25.446 s, CAI 0.9580065306093907, same vector as generations 0. Tier 3 ran 4/4. The sweep reports no infeasible CDS. | `docs/benchmarks/v2.9.1-challenging-perf.md` |

The generations-12 constrained rerun (different CAI values) is
`docs/benchmarks/v2.5.0-validation.md`. The v2.6.0 fronts, including
three-objective knees, are `docs/benchmarks/v2.6.0-phase2.md`. v2.9.1
re-measured the phase 2 vectors and records them as unchanged.

`burden.growth_cost` is a dimensionless demand mismatch. It is not a
predicted growth rate. CAI 1.0 on an empty constraint set means every
counted codon is a top synonym of the loaded table.

## Screenshots

This checkout includes three captures from 17 October 2026
(`docs/screenshots/`, captions in `docs/screenshots/captions.md`). They are
the rpsL demo on the page. They are not the assembly-constraint numbers in
the validation table. On the wiki, upload each file and keep the caption.

| File | Caption |
| --- | --- |
| `docs/screenshots/codonforge-rpsl.png` | CodonForge tab after **Load the rpsL demo** and **Run design**. The status is the finished job and the front shows the reported candidate. |
| `docs/screenshots/primerforge-rpsl.png` | PrimerForge tab after **Load the rpsL sequence** (bases 41–300, flank 40 nt) and **Find primer pairs**. |
| `docs/screenshots/codonforge-pareto.png` | CodonForge front for a job with more than one scored axis (CAI and local structure, or the explain ranges on the page). |

If a file is missing, capture it on a machine with a display after
`python -m codonforge.api`:

1. CodonForge. Open the page. Confirm the ribosomal table note. Click
   **Load the rpsL demo**, then **Run design**. Wait until the status is
   succeeded and the Pareto front has a candidate. Save
   `codonforge-rpsl.png`.
2. PrimerForge. Open the PrimerForge tab. Click **Load the rpsL sequence**.
   That button selects bases 41–300 and sets the flank to 40 nt. Click
   **Find primer pairs**. Save `primerforge-rpsl.png` with the pair table
   in view. **Select the whole sequence** on this CDS leaves no flank for
   the primers to sit in; use it only when the template is longer than the
   amplicon.
3. Pareto. On CodonForge, load rpsL, tick **Local mRNA structure** as well
   as CAI (ViennaRNA must be installed), set budget to Standard under
   Advanced, and run. Save `codonforge-pareto.png` showing the front table
   or the axis ranges. If ViennaRNA is absent, skip this shot and say so;
   do not substitute a CAI-only panel and caption it as a front.

## Licence and reuse

Apache License, Version 2.0: `{{PUBLIC_REPOSITORY_URL}}` file `LICENSE`,
and <https://www.apache.org/licenses/LICENSE-2.0>.

Build, test, and the module layout are in `README.md` and
`CONTRIBUTING.md`. `make check` is the command CI runs. Schemas and the
algorithm notes are in `planning/`. A refused feature returns a `CF-*` code
on the page and on the command line, so a later project can tell a missing
capability from a failed design.
