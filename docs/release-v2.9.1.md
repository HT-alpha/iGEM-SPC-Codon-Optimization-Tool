# DNA Design Suite v2.9.1

Paste this when tagging `v2.9.1`. This note does not create the tag.

DNA Design Suite is CodonForge and PrimerForge on one local page. The
package version stays 0.1.0. `python -m codonforge.cli.main --version`
prints 0.1.0. The suite note is this tag.

## Install

A wheel is not on PyPI. From a checkout, with Python 3.11+ and Rust:

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

The page is <http://127.0.0.1:8000/>. CodonForge is the first tab.
PrimerForge is the second. ViennaRNA is optional and installed by the user
(`pip install 'ViennaRNA==2.7.2'` matches the structure timings below).

Licence: Apache-2.0 (`LICENSE`).

## What this tag measures

- `GATE-BEN-001` is claimed on 50 of 50 proteins. The claim is CAI-only
  against `cai_greedy` (budget `fast`, 4 generations, seed 7). `aceB` and
  `thrA` score 1.0, equal to `cai_greedy`, with no in-frame `AGG`. The other
  48 scores are unchanged from the previous file. The MFE comparison is not
  this claim.
- A structure job reuses an exact accessibility value for a repeated CDS and
  for a short codon edit. The phase 2 knee vectors are unchanged. On the
  recorded machine, standard 4-generation CAI+structure wall time fell from
  1.357 s to 0.561 s on rpsL and from 4.211 s to 2.255 s on waaS.
- `design_operon` updates the flux-weighted codon ledger one codon at a time.
  On lacZ plus thrA (fluxes 1 and 2, budget `fast`, generations 0, seed 7)
  the joint cost stays 0.4452346719259927 and the swaps stay `[26, 0]`. Wall
  time on that machine fell from 9.185 s to 1.632 s.
- The challenging-sequence sweep returned no infeasible CDS. The slowest job
  is yeeJ under the assembly rules, 21.202 s before the population loop
  (25.446 s at budget `standard`, 4 generations, same CAI).

Sources: `CHANGELOG.md` (v2.9.1), `docs/benchmarks/v2.9.1-validation.md`,
`docs/benchmarks/v2.9.1-challenging-perf.md`, `benchmarks/corpus/gates.json`.

## What stays refused

- `epistasis.dwell_deviation` is `CF-LAND-001`. The dwell head is `not_fitted`.
- Host catalogue, model catalogue, and `GET /v1/recalibrations` are
  `CF-CLI-003`. Stored observations do not schedule a recalibration.
- No gene-copy number, s-value, or host-distance penalty was added.
- `GATE-BEN-001` is not an MFE comparison and not LinearDesign or CDSfold.
- `P1-26` stays held. LinearDesign's method is not implemented.

Full wording: `docs/project-status.md` and the Honest limits section of
`README.md`.
