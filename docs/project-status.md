# Project status

Canonical status for the DNA Design Suite at `v2.9.0`. Older notes in
`docs/igem-submission.md` and `CHANGELOG.md` keep the history of how each
release got here. The CAI gate claim is the one recorded for `v2.4.2`. The
v2.6.0 fronts are in `docs/benchmarks/v2.6.0-phase2.md`. The package version
stays 0.1.0. `python -m codonforge.cli.main --version` prints 0.1.0.

## Shipped and working

CodonForge designs an *E. coli* coding sequence under explicit objectives and
hard constraints. The search is Rust (tier 1 exact dynamic programming, tier
2 block coordinate ascent, tier 3 a multi-objective population). The job
lifecycle, ingest, and CLI are Python. `POST /v1/jobs`, `POST /v1/score`,
`POST /v1/observations`, `GET /v1/observations`, `GET /v1/calibration`,
`GET /v1/healthz`, and `GET /v1/version` are served by the standard-library
HTTP server. Hosts, models, recalibrations, seed import, and the five-check
feasibility analyser answer `CF-CLI-003`.

The page in `python/codonforge/api/static/` is the shared CodonForge and
PrimerForge UI. It is a guided form: a coding sequence, a host, the bundled
MG1655 ribosomal usage table, optional constraints, and objectives with the
terms named in plain language. `python -m codonforge.api` binds port 8000,
opens a browser, and blocks until interrupted. The DNA view is the SVG helix
in `helix.js` (CSS keyframes, B-form geometry measured from PDB 1BNA). Primer
selection uses `document.elementFromPoint` so a drag still names the
nucleotide under the pointer when Chrome retargets the event. A region can
also be typed.

Structure scoring uses ViennaRNA, installed by the user and not vendored.
`structure.local_accessibility` is windowed accessibility (`Codons(14)` from
a 40 nt span). When a job requests it, tier 3 ranks on that penalty and the
recommended candidate is informed by it. A request with no engine installed
is `CF-STRU-003`. A CAI-only job does not need the engine. A `local_structure`
rule with a declared span is measured the same way, and that id is not
reported as `CF-FEAS-008`. A missing span is `CF-STRU-004`.
`min_codon_frequency` is still `CF-FEAS-008`.

`burden.growth_cost` is a ranked, dimensionless demand-profile mismatch.
Supply is gene-copy counts when the host profile has them, otherwise the
codon-usage table. A profile with neither is `CF-BURD-004`. No gene-copy
numbers and no s-values are invented. The number is not a growth rate.

A single-objective recommendation is `feasible_then_score_then_digest`. Two
or more complete axes use `knee_of_complete_front`. The page and
`--output explain` print the axis ranges of those vectors.

Primer melting temperatures and dimer ΔG are a clean-room
two-state calculation from SantaLucia 1998 Table 2
(`docs/primerforge-validation.md`). `POST /v1/primers/candidates` and
`POST /v1/primers/analyze` are the primer routes.

The synonymous Potts surface (`epistasis.context_naturalness`) is fit on
`NC_000913.3` with the 13-protein panel held out. The selected model is
`W = 10`: 2333 L-BFGS iterations, gradient norm 0.09770886731134135,
`converged: true`. Held-out naturalness is −0.8708682912309792 nats/codon
against a codon-usage baseline of −0.991738843659599. `GET /v1/calibration`
returns those figures. A design job ranks the term when the verified artifact
is installed. The artifact host is MG1655. A job whose host id starts with
`escherichia_coli:` pins it, and no distance penalty is applied because none
was measured. Any other host is `CF-LAND-001`. The dwell head is not loaded
by that pin. The job raises `tier2_block` to the Potts window. `P1-26` stays
held.

`design_operon` designs two or more cistrons on one flux-weighted codon
ledger and keeps a synonymous swap only when the shared `burden.growth_cost`
falls. The caller supplies every flux.

`POST /v1/observations` appends a validated observation for the calling
tenant. The outcome model is a least-squares line of one uncensored numeric
readout on the GC fraction of an explicit CDS, or `not_fitted` when that line
is not identified. Recalibration is not scheduled.

The benchmark panel is 50 CDS from `NC_000913.3`. The ribosomal usage table
lists all 64 codons. `AGG` and `TAG` have count 0 in the 54 ribosomal CDS, so
their weights are 0. The CAI-only comparison (budget `fast`, 4 generations,
seed 7) is in `benchmarks/corpus/cai_only.json`.

## Claimed

`GATE-BEN-001` is claimed with this scope, and only this scope:

> CAI-only front, 48 of 50 E. coli K-12 proteins (aceB and thrA excluded, see degenerate_cases), cai_greedy baseline, CodonForge meets or exceeds baseline on all 48.

The 48 scores are the 4-generation run. Equality counts as meeting
`cai_greedy`. `check_cai_gate` checks the claim against those rows.
`benchmarks/corpus/gates.json` has `claimed: true` and `n_claimed: 48`.

`aceB` (`NP_418438.1`) and `thrA` (`NP_414543.1`) are excluded. They contain
`AGG` codons with weight 0 in the ribosomal-protein training corpus.
`usage.cai` is identically 0 for any candidate retaining an `AGG` codon. The
greedy baseline wins by construction on a flat-zero objective regardless of
search budget. A further run at budget `fast`, 20 generations, stagnation
window 21, seed 7, still ran 20 generations, still scored 0 against
`cai_greedy` 1, and still left in-frame `AGG` (aceB designed codons 80, 136,
321, 454, 497, 502, 508; thrA designed codons 243, 406, 421, 538, 680, 731,
755, 792).

## Not claimed

- The dwell head is `not_fitted`. No paper states one 5′-end offset for GEO
  GSE53767, and the GEO files were not reachable. `epistasis.dwell_deviation`
  returns `CF-LAND-001`.
- Recalibration is not scheduled. Stored observations do not train a Potts
  model, a dwell head, or a capacity curve.
- Joint search does not run shadow prices or a multi-cistron NSGA-III. The
  shared ledger is the flux-weighted codon demand the caller supplied.
- EXP-05 does not separate the repair policies. On `gapA`, `pfkA`, `mdh`,
  `eno`, and `pykF`, 20 seeds per arm and 40 generations, the hypervolume
  margin is 0 and the Lamarckian standard deviation is 0. Lamarckian stays
  the default. `decide()` was not edited.
- `P1-26` stays held. No registered term has an interaction radius that a
  block cannot honestly cover.
- Commercial distribution is outside the AMD-013 sign-off. The recorded
  review is a non-commercial iGEM educational project, and the reviewer
  states that it is not a formal legal opinion.
- The original MFE/CAI frontier against LinearDesign or CDSfold was not run.
  LinearDesign stays outside AMD-013. Bundled thermodynamic tables, EXP-06,
  and `P3-14` stay held.

## What would make an unclaimed item claimable

- Dwell head: cite one primary-source A-site offset for GSE53767, or apply
  the Li et al. 2014 center-weight rule to reads that can actually be
  downloaded, then fit and pass `verify_dwell`.
- EXP-05 separation: score an objective on which Lamarckian and Baldwinian
  repair change hypervolume; more seeds cannot separate a margin that is
  exactly 0.
- `P1-26`: register an in-core L1 term whose declared radius is greater than
  a block size that cannot honestly be raised, and a provider that evaluates
  the seam.
- Commercial distribution: a new AMD-013 review whose scope includes it.
  The wheel gate stays `tools/check_release_clearance.py`.
