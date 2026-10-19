# Project status

Canonical status for the DNA Design Suite at the current tip. The package
version is 0.1.0.post1. Tag `v0.1.0` is the 0.1.0 release, and that release
is on PyPI. `pip install dna-design-suite` is the install.
`python -m dna_design_suite.cli.main --version` prints the core version,
0.1.0. Older notes in `docs/igem-submission.md` and `CHANGELOG.md`
keep the history. The unreleased delta on top of the v2.9.2 note is the
top of `CHANGELOG.md`. The CAI gate claim is the one in
`benchmarks/corpus/gates.json`. The v2.6.0 fronts are in
`docs/benchmarks/v2.6.0-phase2.md`. The v2.9.1 re-measurement is in
`docs/benchmarks/v2.9.1-validation.md`. The challenging-sequence timings
are in `docs/benchmarks/v2.9.1-challenging-perf.md`. Where those jobs
spent their time on the v2.9.2 tree is
`docs/benchmarks/v2.9.2-profiles.md`. That file is a measurement of that
tree. It is not a version number for this tip.

## Shipped and working

DNA Design Suite designs an *E. coli* coding sequence under explicit objectives and
hard constraints. The search is Rust (tier 1 exact dynamic programming, tier
2 block coordinate ascent, tier 3 a multi-objective population). The job
lifecycle, ingest, and CLI are Python. `POST /v1/jobs`, `POST /v1/score`,
`POST /v1/observations`, `GET /v1/observations`, `GET /v1/calibration`,
`GET /v1/healthz`, and `GET /v1/version` are served by the standard-library
HTTP server. An unknown route is `CF-CLI-003`.

`POST /v1/seeds:import`, `seeds import`, and `starting_cds` accept a
synonymous coding sequence and refuse a different protein with `CF-SEQ-004`
and the residue positions (`DEP-API-003`). `POST /v1/hosts:import` and
`hosts import` turn Kazusa text, a delimited count table, a JSON
`codon_usage` block, or a CDS set into `user_supplied` relative
adaptiveness (`DEP-API-004`). `genome_wide` is refused (`CF-HOST-003`).
Unobserved families stay at weight 0 and are named. No elongation rate,
tRNA copy number, or pseudocount is invented. `GET /v1/hosts` is still the
example files on disk, not a species catalogue.

`POST /v1/feasibility` runs the five checks and reports `satisfiability`
(`DEP-CONS-007`). A synonymous CDS the checker accepts, with no findings,
is `proven_satisfiable`, and `necessary_conditions_only` is false. An
exact finding with no such witness is `proven_unsatisfiable`, and
`necessary_conditions_only` is false. Anything else stays `unknown`, and
`necessary_conditions_only` stays true. `POST /v1/measurements`
stores the same observation document as `POST /v1/observations`
(`DEP-ACT-004`). The GC-fraction line does not change codon choices.

`usage.harmonize_rca` minimizes the mean absolute difference of relative
codon adaptiveness between the input CDS under a source table and the
candidate under the expression-host table (Claassens et al., PLoS ONE 2017,
doi:10.1371/journal.pone.0184355; Angov et al., PLoS ONE 2008,
doi:10.1371/journal.pone.0002189). A missing source table is `CF-HOST-003`.
`usage.match_frequency` minimizes the summed squared error of synonymous
fractions. Its radius is unbounded, so tiers 1 and 2 omit it. A fast budget
warns `CF-SRCH-010` and does not move the seeds for it. Standard and exact
score it in tier 3 (`DEP-LAND-006`). Neither term is a bit-identical port
of another tool, and neither is a measured expression change.

`structure_engine=nussinov` is an opt-in maximum matching for
`structure.local_accessibility` (`DEP-STRU-006`): pair score 1 for A–U,
G–C and G–U, hairpin of at least 3, mean paired fraction in [0, 1]. It is
not a free energy and not an unpaired probability. `select_engine()` does
not choose it. The default engine is still ViennaRNA, installed by the
user. A missing install on that default is `CF-STRU-003`.

`vendor-screen` and `POST /v1/vendor-screen` flag this suite's synthesis
and assembly presets and label the public rule classes those presets
correspond to (`DEP-API-005`). The numeric cutoffs are the presets, not a
transcribed GeneArt, IDT, GenScript, or Twist threshold. `vendor_score` is
false. The screen is not a manufacturer score and not an expression assay.
Wet-lab expression stays out of scope. The checklist is
`docs/vendor-screen.md`.

`examples/gibson.json`, `gibson`, and `POST /v1/gibson` apply the synthesis
preset’s exact 12-mer rule to the 5′ flank, the CDS, and the 3′ flank, and
they copy a homology arm from each supplied flank (`DEP-API-006`). The
default arm is 40 nt, accepted from 20 to 150. No bases are invented.
Overlap melting temperature is not computed. Omitted vector sequence,
circular wrap, mismatches, and inverted repeats are not scanned. A clear
report does not mean the fragments will assemble. The limits are
`docs/gibson-assembly.md`.

`fold` and `POST /v1/structure:global` report a whole-sequence fold
(`DEP-STRU-007`). ViennaRNA `RNA.fold` supplies kcal/mol when that package
imports. Otherwise a sequence of at most 400 nt gets a Watson-Crick pair
count, labelled as base pairs, and a longer sequence is left null. The
report is not an objective, not windowed accessibility, and not the
`nussinov` structure engine.

The published assembly pins recorded in `CHANGELOG.md` are unchanged.

The page in `python/dna_design_suite/api/static/` is the shared DNA Design Suite and
Primer design UI. It is a guided form: a coding sequence, a host, the bundled
MG1655 ribosomal usage table, optional constraints, and objectives with the
terms named in plain language. `python -m dna_design_suite.api` binds port 8000,
opens a browser, and blocks until interrupted. The DNA view is the SVG helix
in `helix.js` (CSS keyframes, B-form geometry measured from PDB 1BNA). Primer
selection uses `document.elementFromPoint` so a drag still names the
nucleotide under the pointer when Chrome retargets the event. A region can
also be typed.

Structure scoring uses ViennaRNA, installed by the user and not vendored,
unless the request sets `structure_engine` to `nussinov`. With ViennaRNA,
`structure.local_accessibility` is windowed accessibility (`Codons(14)` from
a 40 nt span). When a job requests it, tier 3 ranks on that penalty and the
recommended candidate is informed by it. A ViennaRNA request with no engine
installed is `CF-STRU-003`. A CAI-only job does not need the engine. A `local_structure`
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
(`docs/primer-thermo-validation.md`). `POST /v1/primers/candidates` and
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

> CAI-only front, 50 of 50 E. coli K-12 proteins, cai_greedy baseline, DNA Design Suite meets or exceeds baseline on all 50.

The 50 scores are the 4-generation run (budget `fast`, seed 7). Equality
counts as meeting `cai_greedy`. `check_cai_gate` checks the claim against
those rows. `benchmarks/corpus/gates.json` has `claimed: true` and
`n_claimed: 50`.

`aceB` (`NP_418438.1`) and `thrA` (`NP_414543.1`) contain `AGG`. That codon
has weight 0 in the ribosomal-protein table, so `usage.cai` is 0 while any
counted `AGG` remains. The v2.4.2 search left those codons in and scored 0.
`DEP-SRCH-020` changes the search linearization. On this record both designs
score 1.0, equal to `cai_greedy`, and neither designed CDS contains in-frame
`AGG`. The reported CAI definition is unchanged.

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
- `vendor_score` is false. The vendor-class screen does not transcribe a
  manufacturer cutoff and does not report wet-lab expression.
- The global fold report is not an objective. Without ViennaRNA, a sequence
  longer than 400 nt is left null. The Nussinov paired fraction is not a
  free energy.
- An unknown route stays `CF-CLI-003`.

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
