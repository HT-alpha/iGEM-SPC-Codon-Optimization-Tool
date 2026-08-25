# The Co-Translational Folding-Kinetics Objective

Owning module: `cf-biophys-kinetics`.
Related: `AMD-002` (the proxy ladder), `RISK-KIN-001`, `EXP-06`.

This layer is specified separately from `BRIEF-2.3`'s six named algorithms because the brief does not require pseudocode for it, but it is one of the four capabilities the architecture proposal introduces without a reference implementation, so it needs a defined interface, a defined metric and a defined falsification test.

---

## 1. State this plainly: this is the weakest bet in the system

Of the four unaddressed capabilities, three have a clear statistical object and a clear data source: the epistatic surface has Ribo-seq and pseudo-likelihood, the burden model has tRNA pools and growth data, the joint search has a well-posed optimization structure. Co-translational folding kinetics has neither an established predictor nor a routinely-measured readout that could train one at scale.

Building it as a physics simulation in v1 would mean shipping an unvalidated, expensive component on the critical path of every design job. The plan therefore uses a **proxy ladder**: a cheap, testable proxy first, and a simulation only if the proxy demonstrably fails a validation it could have passed.

`AMD-002` records this as an amendment to the architecture proposal, which describes a folding-kinetics objective without specifying its computational form.

---

## 2. v1 metric — pause/domain-boundary concordance

### 2.1 Premise

The mechanistic claim in the literature is narrow and testable: translational pausing at particular points in a CDS gives the nascent chain time to fold before the next segment emerges, and disrupting that pause architecture can reduce the soluble/active fraction even when total protein output is unchanged. The proxy therefore scores *where the pauses are*, not *how the protein folds*.

Two geometric facts set the windows:

- The ribosome exit tunnel holds roughly **30–40 residues**, so a domain that finishes at residue `d` is not fully outside the ribosome until the ribosome reaches roughly residue `d + 35`.
- A pause useful for folding therefore sits **downstream** of a domain boundary, not on it.

### 2.2 Inputs

| Input | Source | Fallback when absent |
| --- | --- | --- |
| Predicted per-codon dwell `z(i)` | `cf-epistasis` dwell head | tAI-derived dwell proxy; term marked `fallback` |
| Domain boundary set `Dom` | user annotation → InterPro/Pfam spans → predicted from a structure model → linker/hydrophobicity heuristic | If none available, the term is **not offered** for this protein; `data_completeness.supported_terms` excludes it |
| Secondary-structure element spans `SSE` | same source as `Dom` | as above |
| Reference dwell profile `z_ref` | the native CDS of this protein, if the input was a native gene | omit the third component |

Refusing to offer the term when no domain information exists is deliberate. A folding-kinetics score computed without knowing where the domains are is noise wearing a domain-specific name.

### 2.3 Metric

```
For each boundary d in Dom:
    desired window   Wd = [ d + 25 , d + 45 ]      (codon indices)
    pause mass       M(Wd) = Σ_{i∈Wd} softplus( z(i) − z_pause )

For each secondary-structure element e in SSE:
    core window      Ce = interior of e, trimmed by 3 codons at each end
    disruptive mass  X(Ce) = Σ_{i∈Ce} softplus( z(i) − z_pause )

kinetics.pause_concordance =
      w₁ · Σ_d  M(Wd) / (1 + M(Wd))          # saturating: one good pause suffices
    − w₂ · Σ_e  X(Ce)
    + w₃ · spearman( z(·), z_ref(·) )        # only when a reference profile exists
```

- Direction: maximize. Unit: `dimensionless`.
- `z_pause` is the host's 99th-percentile native dwell, shared with `epistasis.dwell_deviation`.
- The saturation on `M(Wd)` matters: the goal is *a* pause per boundary, not the deepest possible pause. Without it the optimizer piles rare codons into one window.
- The third term is what makes "do not gratuitously rewrite a native gene's pause architecture" expressible.

**Interaction with the epistasis term.** `epistasis.dwell_deviation` penalizes strong pauses everywhere; `kinetics.pause_concordance` rewards them in specific windows. Left alone these two fight. They are reconciled by the third component of `dwell_deviation` ([`01-epistatic-surface-fit.md`](01-epistatic-surface-fit.md) §3.4), which exempts *designated* pause windows from the pause penalty. `cf-biophys-kinetics` publishes its desired windows into the query context so the epistasis term can read them. This is why the one predictor-to-predictor dependency edge in the module graph exists.

### 2.4 Cost

`O(L)` given the dwell profile, which `cf-epistasis` has already computed. The domain map is computed once per protein and cached. This layer is computationally free; its cost is entirely in the domain annotation, which is an ingest-time concern.

---

## 3. The ladder, and what promotes a rung

| Rung | Method | Promotion condition |
| --- | --- | --- |
| **v1** | Pause/boundary concordance (§2) | Ships in Phase 2 as a **diagnostic** term, not a Pareto axis. |
| **v1 → Pareto axis** | Same metric | `EXP-06` shows a significant association between the score and a measured solubility/activity readout across ≥3 protein families, after controlling for total expression. |
| **v2** | Coarse-grained co-translational folding simulation on the ribosome-elongation timeline | Only if v1 is promoted and then found insufficient — i.e. the association exists but leaves large systematic residuals that a kinetic model explains. |
| **abandon** | Term removed from the default objective set | `EXP-06` finds no association and no host or family where one appears. |

The `diagnostic` role in `job_request.objectives[].role` exists precisely so a term can be computed, reported and validated **without** steering the search while its evidence base is still being assembled. Making an unvalidated objective a Pareto axis on day one would let it distort every design, and would make the validation itself circular.

`DEP-KIN-001` — The architecture proposal lists folding kinetics as one of five co-equal objectives. This plan ships it as a diagnostic term first and promotes it to a Pareto axis only on evidence, because it is the one objective with no validated predictor to build on.

---

## 4. Tests

| Test | Assertion |
| --- | --- |
| `TST-KIN-001` | Window geometry: for synthetic proteins with planted domain boundaries, the desired windows land at `d+25 … d+45` and never overlap a boundary itself. |
| `TST-KIN-002` | Reconciliation: a designated pause window is exempt from the `dwell_deviation` pause penalty, verified by asserting that inserting a strong pause inside a designated window does not increase `dwell_deviation` while inserting one 100 codons away does. |
| `EXP-06` | The falsification experiment for `RISK-KIN-001`. |
