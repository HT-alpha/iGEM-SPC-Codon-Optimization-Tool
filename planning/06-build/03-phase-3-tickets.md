# Phase 3 — The Epistatic Surface

**What Phase 3 delivers:** the layer the architecture is named for, plus the uncertainty machinery that makes its output usable and the registry machinery that makes it versionable. By the end of Phase 3, `EXP-01` has answered whether the central technical bet holds.

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P3-01` | Ribo-seq ingest | `cf-ingest` | Parquet-backed `RiboSeqTrack` loading; QC gates (periodicity, depth, replicate concordance); flag propagation | `TST-INGE-002` | Phase 2 | M |
| `P3-02` | A-site assignment | `cf-ingest` | Fixed, per-read-length and metagene offset methods; 5′/3′ region exclusion; overlapping-ORF masking | Metagene periodicity improves measurably after assignment on a reference track | `P3-01` | M |
| `P3-03` | Potts PLM fitter | `cf-epistasis` | Vectorized pseudo-likelihood objective and gradient; L-BFGS; zero-sum gauge projection; group-`L1` proximal step | `TST-EPI-001` (parameter recovery on `SYN-02`) | `P3-01` | XL |
| `P3-04` | Window selection | `cf-epistasis` | The grid, by-gene CV, the one-standard-error rule, and the data-volume tiering table of ADR-0003 | `EXP-08` runs end-to-end | `P3-03` | M |
| `P3-05` | Potts scorer and deltas | `cf-epistasis` | Rust scorer; per-position field cache; `O(W·S)` exact single-substitution delta; `O(L·W·S²)` delta table; `rebase` | **`TST-EPI-002`** | `P3-03` | L |
| `P3-06` | Dwell head | `cf-epistasis` | NB GLM with gene×replicate offsets, position spline, `dG_local`, tunnel charge, proline, internal SD, batch terms | Held-out NB log-likelihood beats a codon-only GLM | `P3-02`, `P3-03` | XL |
| `P3-07` | Shrinkage prior + `EXP-02` | `cf-epistasis` | Empirical-Bayes `ψ_d ~ N(γ_d·J_d, τ_d²)`; the three-variant comparison | `EXP-02` decided | `P3-06` | M |
| `P3-08` | Dwell terms | `cf-epistasis` | `epistasis.context_naturalness` and `epistasis.dwell_deviation` (rhythm mismatch + pause suppression + designated-window exemption); reference profiles | `TST-KIN-002` (reconciliation with the kinetics term) | `P3-06` | M |
| `P3-09` | Triplet gating | `cf-epistasis` | Residual score statistic, top-`K` selection, refit, hard held-out gate | `EXP-08` triplet portion decided | `P3-06` | L |
| `P3-10` | Deep ensembles | `cf-epistasis` | `K = 5` bootstrap-resampled dwell heads; ensemble artifact packaging; `fast`-tier single-member path | `TST-UNC-001` | `P3-06` | M |
| `P3-11` | Conformal calibration | `cf-epistasis`, `cf-landscape` | Split conformal per term; Mondrian conformal for per-position confidence; quantile artifacts; `calibrated` flag discipline | `TST-UNC-002`, `GATE-UNC-001` | `P3-10` | L |
| `P3-12` | Shared/idiosyncratic error model | `cf-landscape` | Estimation of the shared vs. per-candidate error decomposition during calibration; the joint sampler | `TST-DEC-001` | `P3-11` | L |
| `P3-13` | Probabilistic dominance | `cf-decision` | Monte-Carlo `P(A ≺ B)` over the joint predictive distribution; `robustness.p_nondominated`; expected-hypervolume-contribution reduction | `TST-DEC-001` | `P3-12` | M |
| `P3-14` | LM naturalness term | `cf-biophys-lm` | CodonTransformer/CodonBERT-class scoring in the worker; L0 with optional L2; MC-dropout or checkpoint-ensemble uncertainty | `TST-LM-001`, `TST-LM-002` | `P2-05` | L |
| `P3-15` | Registry channels | `cf-modelreg` | `stable`/`candidate`/`nightly` channels; atomic promotion; gate enforcement; channel history; `yanked` semantics | `TST-MREG-004`, `TST-MREG-005`, `TST-ACT-004` | Phase 2 | M |
| `P3-16` | Nearest-host resolution | `cf-modelreg` | The distance metric, admissibility rules, uncertainty inflation, `low_confidence_host` tagging, `CF-RESOLVE-004` warning | `TST-MREG-001`, `TST-MREG-003` | `P3-15` | M |
| `P3-17` | `EXP-07` runner | `cf-bench` | Cross-host coverage measurement and inflation-factor fitting | `EXP-07` decided; `κ` values adopted or fallback tightened | `P3-16` | M |
| `P3-18` | `bench_epistasis_heldout` + `EXP-01` | `cf-bench` | Within-dataset, cross-dataset and heterologous-protein evaluation against the codon-usage baseline | `GATE-EPI-001`, `GATE-EPI-002`; `EXP-01` decided | `P3-06` | L |
| `P3-19` | Artifact publication | `cf-modelreg` | Public *E. coli* and *S. cerevisiae* artifacts with full manifests, gate reports and citations | Artifacts resolve from a clean install and reproduce their published metrics | `P3-15`, `P3-18` | M |

## Phase 3 exit checklist

- [ ] Phase 2 checklist still green
- [ ] `TST-EPI-*`, `TST-UNC-*`, `TST-LM-*`, `TST-MREG-*` green — **`TST-EPI-002` (delta correctness) above all**
- [ ] `GATE-EPI-001`, `GATE-EPI-002`, `GATE-UNC-001` pass
- [ ] `EXP-01`, `EXP-02`, `EXP-07`, `EXP-08` decided and their outcomes reflected in the shipped defaults
- [ ] If `EXP-01` (c) failed, `epistasis.dwell_deviation` is demoted to `diagnostic` and the change is documented rather than argued away
- [ ] Public artifacts published with manifests that pass the promotion gates
