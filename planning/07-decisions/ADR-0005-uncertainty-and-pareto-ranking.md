# ADR-0005 — Uncertainty Representation and Its Propagation into Pareto Ranking

- **Status:** Accepted
- **Resolves:** open question 3.4 (`BRIEF-3.4`)
- **Related:** [`../02-schemas/pareto_front.schema.json`](../02-schemas/pareto_front.schema.json), [`../03-algorithms/03-nsga-loop-and-operators.md`](../03-algorithms/03-nsga-loop-and-operators.md), `AMD-003`

## Question as posed

> What concrete method (ensembling, Bayesian last layer, conformal prediction) will produce the per-position confidence scores promised in the decision-support layer, and how is that uncertainty propagated into the Pareto ranking rather than shown only as a post hoc annotation?

## Decision

**Method: a per-term-appropriate base estimator, wrapped in split conformal calibration.** Not one method for everything — the terms have genuinely different structure, and forcing a single method on all of them would mean using the wrong one somewhere.

| Term | Base uncertainty | Why |
| --- | --- | --- |
| `epistasis.context_naturalness` | Analytic: the Potts fit is a penalized MLE, so the parameter covariance is available from the Hessian and propagates to the score by the delta method | Free, exact to second order, no ensemble needed |
| `epistasis.dwell_deviation` | Deep ensemble, `K = 5` independent dwell-head fits on bootstrap gene resamples | The dwell head has real model-form uncertainty (covariate choice, dispersion) that a single fit's standard errors understate |
| `burden.growth_cost` | Parametric bootstrap over the capacity fit `(θ, n, g_max, p)` and the tRNA supply vector | The dominant uncertainty is in the *inputs* (ribosome counts, tRNA abundance), not in a learned function |
| `structure.local_accessibility` | Ensemble over folding temperature and window size perturbation | Deterministic given parameters; the honest uncertainty is parameter sensitivity |
| `lm.naturalness` | MC dropout or a small checkpoint ensemble, whichever the chosen model supports | Standard practice for a transformer |
| `usage.cai`, `usage.gc_content` | `method: none`, `calibrated: true` | Deterministic. `none` means *deterministic*, not *unknown* — a distinction the schema makes explicit |

**Calibration: split conformal on a held-out set, Mondrian (class-conditional) for per-position outputs.** The base estimators produce *relative* uncertainty; conformal turns it into an interval with a coverage guarantee. Without it, `sigma` is a number nobody can interpret, and the probabilistic dominance below would be silently miscalibrated. `Uncertainty.calibrated` records whether this step happened; a term that has not passed `GATE-UNC-001` (empirical coverage in `[0.87, 0.93]` for a nominal 90 % interval) is marked `calibrated: false` and treated conservatively downstream.

**Per-position confidence uses Mondrian conformal**, conditioning on position class (5′ ramp / interior / 3′ end) and on local data support (how many training observations the model saw in similar contexts). Plain split conformal gives only *marginal* coverage — 90 % on average over positions — which is exactly the wrong guarantee when the user's question is "is this particular region trustworthy". This distinction is the single most common error in applied conformal work and is called out here so it does not get lost in implementation.

## Propagation into ranking — three mechanisms, not an annotation

The brief's real requirement is that uncertainty affect the ranking. It does, in three places:

### 1. Probabilistic dominance

```
P(A ≺ B)  estimated by Monte Carlo over the JOINT predictive distribution:
    for s in 1..S (default 2000):
        f_A^(s) ~ predictive(A);  f_B^(s) ~ predictive(B)      # correlated, see below
        count += 1 if f_A^(s) dominates f_B^(s)
    return count / S
```

**Correlation is modelled, not assumed away.** Two candidates scored by the *same* model share that model's error: if the epistatic surface is systematically wrong for this protein family, it is wrong for both. Sampling independently would make that shared error look like a difference between candidates. The sampler therefore decomposes each term's error into a **shared component** (model-level, sampled once per Monte Carlo draw and applied to all candidates) and an **idiosyncratic component** (per candidate). The shared/idiosyncratic split is estimated during conformal calibration from the residual structure on the held-out set. Shared error largely cancels in a *comparison*, which is why modelling it makes rankings sharper rather than fuzzier.

The split is carried as `Uncertainty.shared_fraction`, a variance fraction in `[0,1]`. It was added during the `P1-22` sequencing audit: `Uncertainty` previously exposed only a scalar `sigma`, so the composition layer had no field in which to hand the sampler its correlation structure, and the sampler's only options would have been to assume independence or to invent a correlation. Assuming independence is the specific error described in the paragraph above, and it fails in the flattering direction — it makes every ranking look better resolved than the evidence supports. The field is therefore **optional and absence-is-not-zero**: a term whose decomposition was never measured is ineligible for probabilistic dominance and is reported as ineligible rather than sampled under an assumption.

### 2. Robustness as a reported quantity

Each candidate carries `robustness.p_nondominated`: the probability it remains non-dominated under resampling of the whole front. A candidate that is non-dominated only because of a favourable point estimate on a high-variance term will show `p_nondominated ≈ 0.2`, and the user sees it.

### 3. Uncertainty-aware front reduction

The reduced front presented to the user (`reduction.method`) selects reference-direction representatives by **expected hypervolume contribution** under the predictive distribution, not by point-estimate contribution. A candidate wins its niche by being reliably good, not by being optimistically scored.

### Where uncertainty is deliberately *not* used

Inside the tier-3 selection loop, ranking uses **deterministic** dominance on point estimates. Running 2 000 Monte Carlo draws per pairwise comparison per generation would cost more than the evaluations. Probabilistic dominance is applied (a) in the final generations, once the front has stabilized, and (b) in `cf-decision` on the certified front. This is a real approximation and is stated rather than hidden: the search explores using point estimates and *selects* using the predictive distribution.

## Alternatives considered

| Option | Why not as the primary mechanism |
| --- | --- |
| Bayesian last layer everywhere | Cheap and principled for the neural terms, but says nothing about the non-neural ones (burden, structure, CAI), which is most of the objective vector. Available as an alternative base estimator for `lm.naturalness`. |
| Full Bayesian posterior (HMC/VI) over the Potts model | The analytic Hessian already gives a good approximation for a penalized MLE with 41 k parameters, at a fraction of the cost. |
| Conformal alone, no base estimator | Conformal calibrates a *score*; it needs a heuristic uncertainty to calibrate. It is a wrapper, not a replacement. |
| Show uncertainty as an annotation only | This is precisely what the brief rules out, and rightly: an annotation the ranking ignores will be ignored by users too. |

## Consequences

- `Uncertainty` is a required field shape in the landscape contract with mandatory `method` and `calibrated`. An uncalibrated `sigma` entering probabilistic dominance would corrupt rankings invisibly, so the schema forces the distinction to be carried.
- Deep ensembles multiply the dwell head's storage by 5 (~800 KB — irrelevant) and its inference cost by 5 (relevant, and why `fast` tier drops to a single member).
- A fallback-resolved term ([ADR-0009](ADR-0009-artifact-versioning-and-resolution.md)) has its `sigma` inflated and `calibrated` set to `false`, so host fallback automatically depresses those candidates' robustness scores. The two decisions compose without special-casing.
- `AMD-003` records the amendment: the architecture proposal promises "per-position confidence scores" as if it were one number; this plan defines it as a per-term composite with stated semantics, because a single unified confidence over terms with different uncertainty structures would not mean anything.

## Departures

`DEP-UNC-001` — The architecture proposal describes per-position confidence in the decision-support layer. This plan specifies different base estimators per term, a conformal calibration wrapper, an explicitly modelled shared/idiosyncratic error decomposition, and states that the search loop itself uses point estimates while selection uses the predictive distribution.
