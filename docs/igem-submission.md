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

- It does not claim `GATE-BEN-001`. The panel is 50 `NC_000913.3` proteins,
  with native and CAI-greedy arms. LinearFold returns a global MFE and a
  dot-bracket. `structure.local_accessibility` is windowed ViennaRNA
  accessibility. Those fronts are not the same objective, so the gate stays
  unclaimed.
- It does not implement LinearDesign's grammar-intersection method.
- It does not ship a GPU path. `P3-14` (LM scoring) is cut.
- `POST /v1/feasibility` returns `CF-CLI-003`. The five-check analyser is not
  this build. Feasibility on a finished candidate is the boolean
  `constraint_report.feasible`.
- The `W = 10` Potts fit meets the convergence tolerance. The reduced
  gradient `gradient_norm / sqrt(n_positions)` is 8.72e-5, below 1e-4
  (`n_positions` 1,254,595, raw gradient norm 0.09771). That took 2333
  cumulative L-BFGS iterations from zero initialization (40, then 195 in
  15-iteration chunks, then 400, 800, and 898). Limited-memory history was
  reset between those segments. The last segment stopped because the
  tolerance was met, not because the iteration cap was hit. Held-out
  naturalness at convergence is −0.870868 nats/codon, against −0.871175 at
  the 40-iteration early stop and a codon-usage baseline of −0.991739.
  `W = 5` is the unchanged 40-iteration fit and is not the selected model.
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

Twenty seeds per arm on `rpsJ`, `rplL`, and `rpsL`, at 20 generations.
`decide()` was not edited. Hypervolume margin is 0 on all three, so the
rule keeps Lamarckian. The default was not changed. Wall clock for the
120 jobs was 4.3 seconds.

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
  `GATE-BEN-001` stays unclaimed because LinearFold's MFE is not
  `structure.local_accessibility`.
- The epistatic surface is that same host only.
- `W = 5` in the same artifact is still the 40-iteration early stop. It was
  not re-fit. Window selection stays `W = 10`.
- EXP-05 at 20 generations does not separate the arms on hypervolume.
- Wider commercial distribution is outside the AMD-013 sign-off on file.
