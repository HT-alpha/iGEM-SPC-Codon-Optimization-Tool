# CodonForge — iGEM submission notes

This is the state of the tree at the `v1.0.0` stable release. The
`v0.1.0-igem` pre-release is the parent snapshot. `docs/legal-clearance.yaml`
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
was held out of the fit and then scored. The named Ribo-seq dataset is Li,
Burkhardt, Gross and Weissman, Cell 2014, GEO `GSE53767`. The fit does not
consume those read counts. There is no dwell head.

## What it does not do

- It does not claim `GATE-BEN-001`. The panel is 13 proteins, with native and
  CAI-greedy arms only. There is no LinearDesign or CDSfold arm.
- It does not implement LinearDesign's grammar-intersection method.
- It does not ship a GPU path. `P3-14` (LM scoring) is cut.
- `POST /v1/feasibility` returns `CF-CLI-003`. The five-check analyser is not
  this build. Feasibility on a finished candidate is the boolean
  `constraint_report.feasible`.
- The Potts fit is a 40-iteration early stop. The gradient norm at `W = 10`
  is about 1140, above the convergence tolerance. The held-out figure is the
  value at that stop: −0.8712 nats/codon, against a codon-usage baseline of
  −0.9917. It is not a converged surface.
- `P1-26` (`DEP-SRCH-013`) is held.
- `EXP-06`, bundled thermodynamic tables, and a dwell head are not in this
  snapshot.
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
  present it states the window, the held-out nats, the baseline, and the
  early stop.

## What is not claimed

- `GATE-BEN-001`.
- A converged Potts fit.
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
share the budget (`standard`, 4 generations) and the reference point `0` on
`usage.cai`. They do not share a lockfile digest: the sealed document
includes the seed and the constraint set, and the repair flag is part of
that set.

The Potts manifest digest covers `h` and the couplings. The current digest
is `sha256:87a139a46336379fa4e2a955a30e4b04ac5aa9a53f306023deada44bae528f7e`.
`GET /v1/calibration` returns the window, the held-out nats, the baseline,
`iteration_count`, and `final_gradient_norm`.

## EXP-05 re-run

Twenty seeds per arm on `rpsJ`, `rplL`, and `rpsL`. `decide()` was not
edited. Hypervolume margin is 0 on all three, so the rule keeps Lamarckian.
The default was not changed.

Repair invocations, mean per seed:

| Gene | Lamarckian | Baldwinian |
| --- | --- | --- |
| `rpsJ` | 3.35 | 4.75 |
| `rplL` | 0 | 0 |
| `rpsL` | 2.25 | 2.80 |

`rplL` never fired repair under this motif. `rpsJ` and `rpsL` differ, which
is the evidence Baldwinian is repairing the unrepaired genome again.

## Known limitations

- The design benchmark panel is 13 proteins, all from one *E. coli* chromosome.
- The epistatic surface is that same host only.
- The Potts numbers are an early stop, as above.
- EXP-05 at this budget does not separate the arms on hypervolume.
- Wider commercial distribution is outside the AMD-013 sign-off on file.
