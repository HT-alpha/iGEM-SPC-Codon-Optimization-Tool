# CodonForge — iGEM submission notes

The stable release is `v2.4.1`. `v2.4.0` records the dwell-head investigation,
the P1-26 inventory, the 40-generation EXP-05 run, and the first CAI-only
`GATE-BEN-001` comparison. `v2.4.1` puts every codon, including `AGG` at
weight 0, into the ribosomal usage table and re-runs that comparison. The
gate stays unclaimed.
`docs/legal-clearance.yaml`
records the AMD-013 sign-off from NG Chi Hong's email of 2026-10-05. The
review accepts this analysis for a non-commercial iGEM educational project.
The reviewer states that it is not a formal legal opinion. Wider commercial
distribution still needs a separate legal review.

## What it does

CodonForge designs a coding sequence for *Escherichia coli* under explicit
objectives and hard sequence constraints. The search runs tier 1 (exact
dynamic programming for an additive objective), tier 2 (block coordinate
ascent), and tier 3 (a multi-objective population). The result is a Pareto
front. Each candidate carries `constraint_report.feasible` as a boolean and,
separately, `unknown_terms` for any axis that was not scored. A recommended
candidate names `informed_by` and `not_informed_by`. Warnings are rows of
`code`, `stage`, and `subject`.

The pairwise synonymous Potts model (`epistasis.context_naturalness`) is fit
on the CDS of `NC_000913.3` (assembly `GCF_000005845.2`). The 13-protein panel
was held out of the fit and then scored. The selected `W = 10` fit ran 2333
L-BFGS iterations in total, reached gradient norm 0.0977, and is
`converged: true`. `GET /v1/calibration` returns that iteration count, that
gradient, and `converged: true`.

The dwell head (A-site occupancy signal from ribosome profiling) is not
fitted. No paper states one 5′-end offset for GEO GSE53767. Li et al., Cell
2014 (doi:10.1016/j.cell.2014.02.033) size-select about 15–45 nt and give the
same score to center residues at least 10 nt from either end of 20–42 nt
reads. Li, Oh and Weissman, Nature 2012 (doi:10.1038/nature10965, GSE35641)
do the same kind of center-weighting, at 12 nt from either end of 25–42 nt
reads, and do not state +15; the 2014 methods say the size selection differs.
Ingolia et al., Science 2009 (doi:10.1126/science.1168978) is yeast: Fig. 2
puts 5′ ends 12–13 nt upstream of the start codon. Mohammad, Green and
Buskirk, eLife 2019 (doi:10.7554/eLife.42591; there is no 2019 Cell Systems
paper of that name) assign their own 3′ ends, about 12 nt from the A site, on
a flash-freeze library, not GSE53767. RiboToolkit (doi:10.1093/nar/gkaa395)
defaults to 26–32 nt and a user-supplied or plastid P-site. RiboCode
(doi:10.1093/nar/gky179) infers a P-site per read length and does not name
this accession. The GEO files were not reachable, so no counts were extracted
and none were replaced with zero. `GET /v1/calibration` returns
`dwell.state: not_fitted` and names these papers. Requesting
`epistasis.dwell_deviation` returns `CF-LAND-001`. It does not report a
deviation of zero.

## What it does not do

- It does not claim `GATE-BEN-001`. The MFE fronts are not the same
  objective. The CAI-only rerun, after the usage table lists `AGG` at weight
  0, is still short on `aceB` (`AGG` at codon 14) and `thrA` (`AGG` at codons
  28 and 243): CodonForge and the native CDS score 0, and `cai_greedy` scores
  1.
- It does not implement LinearDesign's grammar-intersection method.
- It does not ship a GPU path. `P3-14` (LM scoring) is cut.
- `POST /v1/feasibility` returns `CF-CLI-003`. The five-check analyser is not
  this build. Feasibility on a finished candidate is the boolean
  `constraint_report.feasible`.
- The `W = 10` Potts fit is converged: 2333 total iterations, gradient norm
  0.0977, `converged: true`. Held-out naturalness is −0.870868 nats/codon,
  against a codon-usage baseline of −0.991739. `W = 5` was not re-fit and is
  not the selected model.
- `P1-26` (`DEP-SRCH-013`) is held.
- `EXP-06` and bundled thermodynamic tables are not in this snapshot. The
  dwell head is recorded as not fitted, for the reason above.
- AMD-013 sign-off is recorded in `docs/legal-clearance.yaml` from the
  reviewer's email of 2026-10-05. The scope is a non-commercial iGEM
  educational project, and the reviewer states that the email is not a formal
  legal opinion. Wider commercial distribution needs a separate legal review.
  The wheel gate is unchanged: a stable tag runs
  `legal-clearance`; a hyphenated pre-release tag does not build wheels.

## What is novel here

- A pairwise synonymous Potts context model, fit by pseudo-likelihood with a
  zero-sum gauge, scored on proteins that were not in the fit.
- A host-burden ledger with explicit capacity tiers (the systems-level burden
  model). It is in the library. It is not the objective of the EXP-05 runs
  below.
- A Pareto front that keeps feasibility as a boolean beside per-term
  `UNKNOWN`, instead of imputing a missing score or collapsing the front to
  one winner.
- The result page states when the Potts surface is absent, and when it is
  present it states the window, the held-out nats, the baseline, and whether
  the fit converged.

## What is not claimed

- `GATE-BEN-001`.
- LinearDesign's method, or a comparison against it.
- GPU acceleration.
- A Baldwinian default. The EXP-05 re-run kept Lamarckian.

## Install

ViennaRNA is installed by the user. CodonForge does not vendor it.

```sh
pip install ViennaRNA
pip install codonforge
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

`examples/rpsL.cds` is the *E. coli* K-12 MG1655 gene `rpsL`, protein
`NP_417801.1`, from `NC_000913.3`. The usage table is relative adaptiveness
counted from ribosomal-protein CDS in that file (gene names `rpl`, `rps`,
`rpm`), labelled `user_supplied`.

## Reproducibility

A finished job seals a landscape lockfile. The digest covers the resolved
artifacts, the host profile, the inputs, and the RNG record, and it omits
`created_at`. Two hypervolumes are comparable when they share a budget tier,
a lockfile digest, and a reference point (`TST-BEN-002`). The EXP-05 arms
share the budget (`standard`, 20 generations, stagnation window 21 so the
budget is what stops the run) and the reference point `0` on `usage.cai`.
They do not share a lockfile digest: the sealed document includes the seed
and the constraint set, and the repair flag is part of that set.

The Potts manifest digest covers `h` and the couplings. The current digest
is `sha256:ac64cae265510c2cad9c3a85a59b9abf389cae456eba25fc4595e21d20e645b6`.
`GET /v1/calibration` returns the window, the held-out nats, the baseline,
`iteration_count` (2333), `final_gradient_norm` (0.09771), and
`converged: true`.

## EXP-05 re-run

EXP-05 ran on two of three panel proteins (rplL infeasible under constraint
set); hypervolume margin was zero on both at 20 generations, Lamarckian
default retained.

Twenty seeds per arm were started on `rpsJ`, `rplL`, and `rpsL`.
`decide()` was not edited. Wall clock for the 120 jobs was 4.3 seconds.

`rpsJ` and `rpsL` ran 20 generations on every seed. `rplL` runs 0
generations: tier 1 returns `CF-SRCH-006` (no feasible assignment under
the interior motif) on both arms, so that protein does not inform the
repair comparison. Its margin is 0 because both arms fail the same way.

Repair invocations, mean per seed, on the jobs that reached tier 3:

| Gene | Lamarckian | Baldwinian | Hypervolume margin |
| --- | --- | --- | --- |
| `rpsJ` | 13.95 | 24.55 | 0 |
| `rplL` | 0 (tier 1 refused) | 0 (tier 1 refused) | 0 |
| `rpsL` | 10.3 | 15.0 | 0 |

Population diversity stays 0 on the Lamarckian arm and about 0.001 on
the Baldwinian arm for `rpsJ` and `rpsL`. That is not a hypervolume
separation. The finding at this budget is inconclusive, and Lamarckian
stays the default.

## PrimerForge

PrimerForge, added in the DNA Design Suite after the 2026-10-05 sign-off,
designs primer pairs for a region the user selects on a DNA sequence and
also checks one pasted pair. Melting temperatures and dimer ΔG use SantaLucia
1998 Table 2 (`docs/primerforge-validation.md`), as a digested parameter
document, in a clean-room two-state calculation. It does not check specificity
against a genome, does not wrap Primer3, and does not claim mFold parity for
hairpins. The forward CaPIP self-dimer is two helices in one register,
ΔG°37 −3.556 kcal/mol against the published OligoAnalyzer cell −3.55; an
earlier build that kept only the first helix reported −0.810.

## Known limitations

- The design benchmark panel is 50 proteins, all from one *E. coli* chromosome.
  `GATE-BEN-001` stays unclaimed. The CAI-only comparison fails on `aceB` and
  `thrA` because `AGG` has weight 0 and the search leaves it in place.
- The epistatic surface is that same host only.
- `W = 5` was not re-fit. Window selection stays `W = 10`.
- EXP-05 ran on two of three panel proteins (`rplL` infeasible under the
  constraint set). Hypervolume margin was zero on both at 20 generations.
  Lamarckian stays the default.
- Wider commercial distribution is outside the AMD-013 sign-off on file.
