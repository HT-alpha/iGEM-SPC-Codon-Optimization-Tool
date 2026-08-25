# Active-Learning Ingestion and Non-Blocking Recalibration

Satisfies: `BRIEF-2.4-b` — "the active-learning ingestion API schema … and how a submitted result triggers model recalibration without blocking other users' jobs".
Owning module: `cf-active`.
Schema: [`../02-schemas/observation.schema.json`](../02-schemas/observation.schema.json).
Related: [ADR-0006](../07-decisions/ADR-0006-active-learning-data-isolation.md), `AMD-004`.

---

## 1. The pipeline

```
POST /v1/observations
        │
        ▼
   ┌────────────┐   reject (422, CF-OBS-*)
   │ VALIDATE   │──────────────────────────▶ nothing is queued
   └────┬───────┘
        │ accept (202)
        ▼
   ┌────────────────────────┐
   │ APPEND to tenant       │   append-only, content-addressed, never mutated
   │ observation dataset    │
   └────┬───────────────────┘
        │
        ▼
   ┌────────────────────────┐   debounced: fires on N new observations,
   │ SCHEDULE recalibration │   or T elapsed, or explicit request
   └────┬───────────────────┘
        │
        ▼   (separate worker, separate hardware, its own queue)
   ┌────────────────────────┐
   │ TRAIN a NEW artifact   │   never modifies an existing artifact
   └────┬───────────────────┘
        │
        ▼
   ┌────────────────────────┐   fails ⇒ artifact stays `draft`,
   │ ACCEPTANCE GATES       │   incumbent keeps serving, tenant is told why
   └────┬───────────────────┘
        │ pass
        ▼
   ┌────────────────────────┐
   │ PROMOTE: move the      │   a pointer move. Existing jobs are unaffected
   │ tenant channel pointer │   because they pinned a lockfile at submission.
   └────────────────────────┘
```

---

## 2. Why this cannot block another user's job

Three structural properties, each of which is independently sufficient to make the interference impossible:

1. **Every job pins a `LandscapeLockfile` at submission.** A job resolves its artifacts once, records digests, and reads only those digests for its whole lifetime. A new artifact version appearing mid-job is not merely non-blocking — it is *invisible*. This is the main argument, and it is why the lockfile exists.
2. **Artifacts are immutable and content-addressed.** Recalibration writes a *new* artifact. Nothing an existing job reads is ever mutated, so there is no lock to contend for and no torn read to guard against.
3. **Recalibration runs in a separate worker pool with its own queue and its own resource limits.** It is scheduled work, not request-path work. A tenant submitting 500 observations cannot starve design jobs, because the two never share a queue.

Promotion is a **pointer move on a channel** — an atomic write of a small JSON document naming an artifact id and digest. Jobs submitted after the move resolve the new artifact; jobs submitted before it do not. There is no window in which a job sees a half-updated model.

`DEP-ACT-001` — The architecture proposal describes an active-learning feedback loop without specifying its concurrency semantics. This plan makes the lockfile the isolation mechanism, which also delivers reproducibility (`BRIEF-2.2-d`) and query-time version selection (`BRIEF-2.4-c`) from the same construct.

---

## 3. Validation, and what each check is defending

Performed synchronously; a failure means nothing is queued.

| Check | Code | Why it matters |
| --- | --- | --- |
| Schema conformance | `CF-OBS-001` | — |
| `as_built_digest` matches the referenced candidate | `CF-OBS-002` | Divergence between the designed and the actually-synthesized sequence is common (vendor edits, cloning scars) and silently corrupts the training set. This check is the single highest-value validation in the whole pipeline. |
| At least one control present | `CF-OBS-003` | A readout with no control cannot be placed on a common scale with any other lab's readout, which is the entire point of pooling. |
| `dynamic_range.censoring` declared | `CF-OBS-004` | A value at the top of an assay's range is censored, not measured, and must enter the likelihood as such. Treating a saturated fluorescence reading as a point observation teaches the model that a design is exactly as good as the plate reader's ceiling. |
| `measures` (latent quantity) declared | `CF-OBS-005` | Two assays measuring the same latent quantity are poolable; two measuring different ones are not. The *assay* does not determine this — a Western blot and a fluorescence reading can both measure protein abundance. |
| `exploration` block present | `CF-OBS-006` | See §5. |
| Host and condition resolvable | `CF-OBS-007` | — |
| Replicate count consistent with `replicate_values` | `CF-OBS-008` | — |

---

## 4. What actually gets recalibrated

Not everything is retrained from scratch. The recalibration kind depends on what the observation can inform:

| Readout `measures` | Recalibration target | Method |
| --- | --- | --- |
| `protein_abundance`, `translation_efficiency` | `epistasis_dwell` | Tenant **adapter**: a low-rank/residual correction on the dwell head's coefficients, plus a per-tenant intercept and batch effects. The base model is untouched. |
| `host_fitness` | `burden_capacity` | Refit `θ`, `n`, `g_max`, `p` (capacity tier 1) on the tenant's own constructs. This is the highest-value recalibration available, because tier-1 capacity is otherwise unattainable. |
| `solubility`, `activity` | `kinetics_domainmap` | Only once `EXP-06` has promoted the kinetics term; before that, observations are stored but not used for fitting. |
| `mrna_abundance` | none in v1 | Transcription is out of scope for a codon-optimization tool; stored for future use. |
| Ribo-seq submitted as an observation | `epistasis_potts` + `epistasis_dwell` | Full tenant refit, since a tenant Ribo-seq dataset is a genuine training set rather than a correction. |

**Adapters rather than full fine-tunes** for the common case, for three reasons: they are cheap enough to run on a handful of observations, they cannot destroy the base model's behaviour outside the tenant's own region, and they make the isolation invariant of [ADR-0006](../07-decisions/ADR-0006-active-learning-data-isolation.md) checkable by inspecting a small artifact rather than a whole set of weights.

---

## 5. The feedback loop is a hazard, and the design treats it as one

Retraining a model on data generated by that model's own top-ranked suggestions is a closed loop. The model becomes *more* confident in the region it already liked and *less* informed everywhere else — its held-out metrics can improve while its usefulness for new proteins degrades. This is not hypothetical; it is the default outcome of naive active learning.

Three mitigations, all mandatory:

1. **Exploration arms.** `observation.exploration` is a required block. The recommended protocol, surfaced by the UI and the CLI, is that a submitted batch contain at least 25 % of constructs drawn by `uniform_random_from_front` or `max_uncertainty` rather than `top_ranked`. The fitter reweights by selection rule, so a batch that is 100 % top-ranked is *usable* but contributes less and is flagged.
2. **Selection-rule-aware fitting.** The recalibration model includes the selection rule as a covariate and uses inverse-propensity weighting where the rule is known. An observation with `selection_rule: unknown` is downweighted.
3. **A frozen, never-recalibrated holdout.** Each tenant has a holdout set that recalibration never trains on and never selects from. Gate `GATE-ACT-001` is evaluated on it. If the closed loop is degrading the model, the frozen holdout is where it becomes visible.

`AMD-004` records this as an amendment: the architecture proposal describes the active-learning loop as straightforwardly beneficial, and it is beneficial only with these guards.

---

## 6. Acceptance gates before promotion

A recalibrated artifact is `draft` until every required gate passes. Gates are defined in [`../05-testing/02-acceptance-thresholds.md`](../05-testing/02-acceptance-thresholds.md).

| Gate | Requirement |
| --- | --- |
| `GATE-ACT-001` | No regression on the tenant's frozen holdout beyond a tolerance band, versus the incumbent. |
| `GATE-ACT-002` | No regression on the **public** holdout beyond tolerance. A tenant adapter that improves tenant fit by destroying general behaviour is rejected. |
| `GATE-ACT-003` | Calibration retained: empirical coverage of the 90 % predictive interval stays within `[0.87, 0.93]`. |
| `GATE-ACT-004` | Isolation: the manifest's `training.datasets` contains no dataset belonging to another tenant. Mechanically checked; see [ADR-0006](../07-decisions/ADR-0006-active-learning-data-isolation.md). |
| `GATE-ACT-005` | Minimum evidence: at least 12 observations from at least 2 distinct batches. Fitting a capacity curve to one plate produces a confident, wrong answer. |

A rejected recalibration leaves the incumbent serving, records the gate report, and surfaces it at `GET /v1/recalibrations`. The tenant is told *which* gate failed and by how much — "your recalibration was rejected" without a reason is useless.

---

## 7. Scheduling

| Trigger | Default |
| --- | --- |
| New-observation count | 12 since the last recalibration for that `(kind, host)` |
| Elapsed time | 7 days, if any new observations exist |
| Manual | `POST /v1/recalibrations` (tenant-initiated) |

Debounced with a 1-hour quiet period: a lab uploading a 96-well plate one row at a time triggers one recalibration, not eight. Concurrent recalibrations for the same `(tenant, kind, host)` are collapsed; the later request supersedes the earlier.

---

## 8. Tests

| Test | Assertion |
| --- | --- |
| `TST-ACT-001` | Validation: each check in §3 fires on a purpose-built bad submission and stays silent on a good one. |
| `TST-ACT-002` | Isolation: a tenant-B artifact's manifest never lists a tenant-A dataset, asserted over a simulated multi-tenant history. |
| `TST-ACT-003` | Non-blocking: a job running while a recalibration promotes a new artifact mid-flight produces a front bit-identical to the same job run with no concurrent promotion. This is the direct, executable test of the `BRIEF-2.4-b` requirement. |
| `TST-ACT-004` | Gate enforcement: an artifact failing any required gate is never reachable through a channel, only by explicit digest. |
