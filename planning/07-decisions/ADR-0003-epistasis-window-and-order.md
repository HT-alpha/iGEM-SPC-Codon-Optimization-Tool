# ADR-0003 — Epistasis Window Size and Interaction Order

- **Status:** Accepted, with two parameters deferred to named experiments
- **Resolves:** open question 3.1 (`BRIEF-3.1`)
- **Related:** [`../03-algorithms/01-epistatic-surface-fit.md`](../03-algorithms/01-epistatic-surface-fit.md), `EXP-01`, `EXP-08`, `AMD-006`

## Question as posed

> Is a fixed window (e.g. 10 codons) sufficient, or should window size/interaction order be learned or made host-dependent? What is the fallback when a host has too little ribosome-profiling data to fit reliable higher-order terms?

## Decision

**The window is a host-dependent hyperparameter selected by held-out likelihood over a small grid, with a data-volume-tiered default. The interaction order is pairwise, with third-order terms gated behind a hard held-out improvement test that they must pass per host. Data-poor hosts fall back to a smaller window, a pooled cross-host prior, and the context model alone.**

The brief explicitly permits deferring a decision to a named experiment when it cannot be made without data that does not yet exist. Window size is that case: nobody can know the right `W` for *Pichia pastoris* without fitting it. What planning *can* fix — and what this ADR fixes — is the **selection procedure, the grid, the defaults, and the fallbacks**, so that no further design discussion is needed when the data arrives.

### 1. Window selection procedure (fixed now)

```
SELECT_WINDOW(host) -> W

  grid <- [3, 5, 10, 15, 20]                     # codons
  restrict grid by the data tier in §2
  for W in grid:
      fit the context model on a by-gene training split
      score held-out pseudo-likelihood per codon
  W* <- the SMALLEST W whose held-out score is within 1 standard error
        of the best (one-standard-error rule)
  return W*
```

The one-standard-error rule rather than plain argmax, because the held-out curve is typically flat past the true window and plain argmax then picks the largest grid point — inflating the parameter count by `4096` per extra offset and slowing every delta query, for no measured gain.

### 2. Data-volume tiering (fixed now)

`riboseq_mapped_codons` from `HostProfile.data_completeness` is the tiering variable.

| Ribo-seq mapped codons | Window grid | Order | Prior | Terms offered |
| --- | --- | --- | --- | --- |
| **0** (no Ribo-seq) | `{5, 10}`, selected on native-CDS held-out likelihood | pairwise | pooled cross-host empirical-Bayes prior | `epistasis.context_naturalness` only. **No dwell head.** |
| **< 5 × 10⁵** | `{3, 5}` | pairwise | pooled prior, strong shrinkage | context + dwell, dwell flagged low-confidence |
| **5 × 10⁵ – 5 × 10⁶** | `{5, 10, 15}` | pairwise | weak shrinkage | context + dwell |
| **> 5 × 10⁶** | `{5, 10, 15, 20}` | pairwise, triplets **gated** | none | context + dwell + gated triplets |

Two things to note about this table. First, the zero-Ribo-seq row is not a degraded mode, it is the *common* case: the context model needs only an annotated genome, so every host gets a real epistatic term. Second, the dwell head is simply absent below its data threshold rather than fit badly — `data_completeness.supported_terms` omits it and the API refuses a job requesting it, instead of returning a number nobody should use.

### 3. The default before any fit: `W = 10`

`W = 10` is the default that ships and the value used when a host has not yet been fit. It is not arbitrary: the bacterial ribosome footprint is ~28–30 nt ≈ 10 codons, so a 10-codon window is the span over which the ribosome is physically in contact with the message. It is an anchor, not a measurement, and `EXP-08` replaces it with a fitted value per host.

### 4. Interaction order: pairwise, with a gate

Full third-order coupling needs `64³ × C(W,2) ≈ 1.2 × 10⁷` parameters against `~10⁶` codon observations, which is not a modelling choice. The procedure in [`../03-algorithms/01-epistatic-surface-fit.md`](../03-algorithms/01-epistatic-surface-fit.md) §5 selects at most `K = 2000` triplets by a residual score statistic and then applies a hard gate:

> Accept the triplet extension **iff** held-out log-likelihood per codon improves by `≥ δ` **and** the bootstrap 95 % CI of the improvement excludes zero. Otherwise ship the pairwise model.

`δ = 0.005` nats/codon initially, to be replaced by whatever `EXP-08` finds is distinguishable from noise. The gate is a hard requirement, not a preference: a triplet extension that cannot beat its own held-out baseline is an overfit, and shipping it silently degrades every design produced for that host.

### 5. Data-poor fallback: hierarchical pooling

For hosts below `5 × 10⁵` mapped codons, couplings are shrunk toward a pooled cross-host estimate:

```
J_host,d  ~  Normal( J_pooled,d ,  τ_d² )
```

with `J_pooled` fit across all public hosts and `τ_d` estimated by empirical Bayes. Justification for pooling at all: the couplings being estimated are largely a property of ribosome geometry and wobble decoding, which are far more conserved across organisms than codon *usage* is. Justification for the shrinkage being per-offset (`τ_d`): near-neighbour couplings (E/P/A site geometry) should be more conserved than distant ones, and letting the data say so is cheap.

`EXP-01` tests the pooling assumption directly by holding out one host at a time and measuring whether the pooled prior beats a uniform prior. If it does not, the fallback becomes "reduce `W` to 3 and use the context model only" — which is the honest floor.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Fixed `W = 10` for all hosts | Simple, but the footprint argument is bacterial. Eukaryotic ribosome footprints and pause structure differ, and there is no reason a yeast window equals an *E. coli* one. Fixing it forecloses a question the data can answer cheaply. |
| Learn `W` per *gene* | Overfits, and destroys the shared delta cache (which assumes a fixed window for the whole sequence). |
| Learn the window by group-lasso only, no CV grid | Attractive — the group-`L1` penalty already selects offsets, and the plan *does* use it as evidence. Rejected as the sole mechanism because the selected support depends on `μ`, which itself needs cross-validation, so the grid is unavoidable. Both signals are reported. |
| Attention over a long context (a transformer instead of a Potts model) | Would learn the window implicitly and is the natural v2. Rejected for v1 because it forfeits the two properties this layer is built on: a concave objective with a unique optimum, and `O(W·S)` exact deltas. Both are what make tier 2 affordable and the artifact auditable. Recorded as `RISK-EPI-002`. |

## Consequences

- `ModelManifest.parameters` records `window_codons` and `interaction_order` per artifact, so two hosts can legitimately have different windows and a result stays interpretable.
- Delta cost is `O(W·S)`, so a host fit at `W = 20` costs twice as much per move as one at `W = 10`. The one-standard-error rule is what keeps this from happening for free.
- A host with no Ribo-seq is a fully supported case with a documented, reduced capability set — see `AMD-006`.

## Departures

`DEP-EPI-004` — The architecture proposal treats the epistasis window as a fixed design parameter. This plan makes it a fitted, host-dependent, tiered hyperparameter, and makes higher-order terms conditional on passing a per-host held-out gate rather than being a design feature.
