# iGEM batch 2 — sequencing, and the Phase 3 minimum

This note records the sequencing check that decided what batch 2 builds, and
the Phase 3 scope that is planned and not coded. It does not move a gate.

## 1. Sequencing check

### 1.1 API routes the demo actually calls

The result page submits a library `DesignRequest` and then reads the job and
the front. Those routes are already real, and they are the ones this round
keeps:

| Route | Why the page needs it |
| --- | --- |
| `POST /v1/jobs` | Submit the design |
| `GET /v1/jobs/{id}` | Per-tier progress, resolved-objective roles, warning rows |
| `GET /v1/jobs/{id}/front` | The Pareto document the page renders |
| `GET /v1/jobs/{id}/lockfile`, `POST /v1/jobs/{id}/cancel`, `POST /v1/score` | Not drawn by the page. Kept because a reproducibility claim and a benchmark claim both use them, and they already call the library |

The page does not call a host catalog, a model catalog, an observation
store, a preset expander, or the five-check feasibility analyser. Promoting
those would be building surfaces the demo does not exercise.

They stay `CF-CLI-003`, with a reason in the problem document (`DEP-API-001`):

| Surface | Why it is not promoted |
| --- | --- |
| `input.kind` of `protein`, `construct`, or `single_cds` | This build does not reverse-translate a protein and does not expand a construct. The demo pastes a CDS |
| `objectives_preset` | A preset would hide which term ids the caller asked for. The page's `informed_by` / `not_informed_by` lists need those ids |
| `GET /v1/hosts`, `GET /v1/models` | No catalog is installed. A job carries the caller's `host_profile`. An empty list would look like a catalog. Public channels are `P3-15` |
| `POST /v1/observations`, `GET /v1/recalibrations` | Wet-lab ingestion is not this release |
| `POST /v1/seeds:import` | The route is not wired. An external baseline, if one is used, is a user-installed executable invoked by the harness, never a vendored binary |
| `POST /v1/feasibility` | The five checks in `03-algorithms/06-constraint-repair.md` §5 (`CF-FEAS-001` through `CF-FEAS-005`) are not this build. The CLI's `constraint_violations` report is a detector pass, not that analyser. Wrapping the detector and returning `200` would look like the analyser |

`RISK-TIME-004`.

### 1.2 What the benchmark claims

`GATE-BEN-001` needs at least 50 proteins and a competitor arm scored on
`structure.mfe_global` and `usage.cai`, with CodonForge hypervolume at least
0.98× the best competitor on at least 90 % of them. The full-objective suite
additionally wants at least 20 seeds. This build has neither the competitor
nor the 20 seeds.

The claim that is defensible is the one `DEP-BEN-002` writes into
`gates.json` with `claimed: false`:

**50 proteins, 0 seeds, two in-process arms (native and CAI-greedy).
`GATE-BEN-001` is not claimed.** LinearFold emits a global MFE and a
dot-bracket. The scored structure term is windowed accessibility. The
fronts are not compared. The original 13 genes remain the Potts held-out
set.

The sequences are the CDS of `NC_000913.3` (assembly `GCF_000005845.2`,
4140 records in the cited FASTA). The gzip's sha256 is
`72b9a16f2c6c48ac021288ad6922b8b7ccf608e31e4ae807eea65b2cdce64e1a`, and it
matches the file retrieved from the Rice 2018 iGEM PORTAL commit
`4337d757da6b54cb3a0d937a422f26a1bd713953`. `verify_shipped_cds` refuses a
file whose digest is not that value.

| Gene | Why it is on the list |
| --- | --- |
| `lacZ` | Reporter. β-galactosidase, `NP_414878.1` |
| `gapA`, `eno`, `pykF`, `pfkA`, `adk` | Metabolic enzymes |
| `dnaK`, `tufB` | Abundant cellular proteins |
| `rpsA`, `rplA`, `rplL`, `rpsL`, `rpsJ` | Ribosomal proteins |

Expression evidence is the citation, not a number: Li, Burkhardt, Gross and
Weissman, Cell 2014, doi:10.1016/j.cell.2014.02.033, GEO `GSE53767`. The
panel's `expression_evidence.value` is null, and a non-null value is refused.
CAI weights are relative adaptiveness counted from the ribosomal-protein CDS
in the same file (54 records whose gene names start with `rpl`, `rps`, or
`rpm`). That geometric mean is not reported as the engine's `usage.cai`.

`tufA` is excluded because its annotated start is GTG and `translate` reads
internal GTG as valine. GFP is excluded: no nucleotide CDS with an INSDC
accession was retrieved, and a protein sequence was not reverse-translated
(`RISK-TIME-006`).

### 1.3 What the page has to show

All six disclosures are essential. Each one is a field the job or the front
already carries, and dropping it would make the page vaguer than the
document:

- recommended versus structure-optimal, with `informed_by` and `not_informed_by`
- per-term `UNKNOWN` on a hatched rail, not imputed onto the numeric axis
- diagnostic versus axis, from `objective_detail.role` and from which column the term sits in
- warning rows with `code`, `stage`, and `subject`
- `constraint_report.feasible` as a boolean, with `unknown_terms` beside it
- one progress track per tier

A preset picker, a host catalog, a checkpoint URI, and a live hypervolume
chart are not required for that reading. The page does not invent a host
profile. The design language is the Part D brief: Source Serif 4, IBM Plex
Sans, IBM Plex Mono, the named state colours, and the 5/10/20/40 spacing.
This environment has no frontend-design skill and no Visualizer; that
substitution is the one the brief already discloses.

The structure-optimal mark is only for a term whose id starts with
`structure.`. Another excluded axis can have a best scored member, and that
member is named, but it is not painted as structure-optimal.

## 2. Phase 3 plan

The pairwise fit described below was run after this note. The dwell head
was not. Training used 4127 CDS; the 13 panel genes were removed first.
At 40 L-BFGS iterations the gradient was still above the convergence
tolerance. Held-out nats per codon: W=5 −0.9005, W=10 −0.8712, codon-usage
baseline −0.9917. The one-standard-error rule selected W=10. The digest
and the parameters are in `benchmarks/epistasis/manifest.json`.

### 2.1 One host, one Ribo-seq dataset

Host: *Escherichia coli* K-12 MG1655, the same assembly as the panel
(`GCF_000005845.2`, `NC_000913.3`).

Dataset: Li, Burkhardt, Gross and Weissman, "Quantifying absolute protein
synthesis rates reveals principles underlying allocation of cellular
resources", Cell 2014, GEO accession `GSE53767`. GWIPS-viz
(Michel et al., Nucleic Acids Research 2018, the 2018 update) hosts
*Escherichia coli* K12 ribosome-profiling tracks and lists this paper as
reference 121. That is the single dataset. The harness note
`bm_riboseq_ecoli` asks for at least two independent sets reprocessed
identically; the second set is not part of the iGEM minimum.

The pairwise context model below does not consume the read counts. The
dataset is named now so the dwell head, when it is fit, has one source and
not a choice made later.

### 2.2 Pairwise epistatic surface (`P3-05`)

The fit this deliverable means is component A of
`03-algorithms/01-epistatic-surface-fit.md`: a synonymous-restricted Potts
model on the annotated CDS, pseudo-likelihood, L-BFGS, zero-sum gauge inside
synonymous families. The dwell head is component B and needs A-site counts
from `GSE53767`. It is not in the minimum.

Data the fit needs: the 4140 CDS already in the panel file. No expression
rate, and no elongation-rate constant. `ADR-0003`'s zero-Ribo-seq row is the
window grid, `{5, 10}`, chosen by held-out pseudo-likelihood with the
one-standard-error rule. Parameter count is `64 + 4096·W`: 20,544 at `W = 5`,
41,024 at `W = 10`. One process. The objective is concave, so the start is
the zero vector and there is no seed. No GPU.

A wall-clock figure is not stated here, because this batch did not run the
fit. The planning claim that is already on record is that this parameter
count fits on a laptop (`01-epistatic-surface-fit.md` §2.1). That is the
budget. A measured time belongs to the run that produces the artifact.

Output: an artifact of `(h, J_d)` with a manifest; held-out pseudo-likelihood
in nats per codon against the codon-usage baseline on a by-gene split; and
`epistasis.context_naturalness` scored on the 13-protein panel. That score
is a report. It is not a design recommendation until the term is a
registered provider.

### 2.3 `P3-14` is cut

LM scoring (`CodonTransformer` / `CodonBERT`-class, `P3-14`) is not in the
iGEM build. `ADR-0003` already rejected an attention model as the v1 surface
because it gives up the concave optimum and the `O(W·S)` deltas tier 2 is
built on. That rejection is `RISK-EPI-002`. The additional cost here is
model weights, a GPU path, and a licence review that `P2-23` has not signed.

Cutting it means this build has no `lm.naturalness` term. The loss is a
long-context prior the pairwise model does not represent. The term can be
reconsidered after the pairwise surface has a held-out number and the legal
gate is signed. Until then a request for the term is `CF-LAND-001`, which is
the existing refusal for a term this build does not ship (`RISK-TIME-003`).

### 2.4 Minimum viable Phase 3, one sentence

Fit a pairwise synonymous Potts model on the `NC_000913.3` CDS set, score
`epistasis.context_naturalness` on the 13-protein panel against a codon-usage
baseline, and keep `GSE53767` as the single named Ribo-seq dataset for a
dwell head this deliverable does not fit.

## 3. Still held

`P1-26` (`DEP-SRCH-013`), `EXP-06`, bundled thermodynamic parameter tables,
and LinearDesign stay held. The `EXP-05` default stays Lamarckian until a
campaign of at least 20 seeds per arm (`RISK-TIME-005`). The release gate
stays blocking (`RISK-TIME-002`). The toolchain pin stays rustc 1.98.1.
