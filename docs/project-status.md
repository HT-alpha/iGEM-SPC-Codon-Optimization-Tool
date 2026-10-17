# Project status

Canonical status for the DNA Design Suite at `v2.5.0`. Older notes in
`docs/igem-submission.md` and `CHANGELOG.md` keep the history of how each
release got here. The benchmark claim is the one recorded for `v2.4.2`.

## Shipped and working

CodonForge designs an *E. coli* coding sequence under explicit objectives and
hard constraints. The search is Rust (tier 1 exact dynamic programming, tier
2 block coordinate ascent, tier 3 a multi-objective population). The job
lifecycle, ingest, and CLI are Python. `POST /v1/jobs`, `POST /v1/score`,
`GET /v1/calibration`, `GET /v1/healthz`, and `GET /v1/version` are served by
the standard-library HTTP server. Hosts, models, observations,
recalibrations, seed import, and the five-check feasibility analyser answer
`CF-CLI-003`.

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
a 40 nt span). Primer melting temperatures and dimer ΔG are a clean-room
two-state calculation from SantaLucia 1998 Table 2
(`docs/primerforge-validation.md`). `POST /v1/primers/candidates` and
`POST /v1/primers/analyze` are the primer routes.

The synonymous Potts surface (`epistasis.context_naturalness`) is fit on
`NC_000913.3` with the 13-protein panel held out. The selected model is
`W = 10`: 2333 L-BFGS iterations, gradient norm 0.0977, `converged: true`.
Held-out naturalness is −0.870868 nats/codon against a codon-usage baseline
of −0.991739. `GET /v1/calibration` returns those figures.

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
