# Experiment Register

Satisfies `BRIEF-4` — "Where a decision genuinely cannot be made without empirical data that doesn't yet exist … specify the experiment or benchmark that would resolve it, rather than picking arbitrarily and hiding the uncertainty" — and `BRIEF-5.9`, which allows an open question to be closed by a *named experiment* provided it runs before the relevant build phase.

Every entry states the question, why planning cannot answer it, the design, the data required, **a decision rule fixed in advance**, which phase it gates, and what happens if it fails. The decision rule is the important part: an experiment without a pre-registered decision rule becomes a justification exercise.

| ID | Question | Gates | Fallback if it fails |
| --- | --- | --- | --- |
| `EXP-01` | Does the epistatic surface generalize across datasets and proteins? | Phase 3 | Restrict to same-dataset use; demote the term |
| `EXP-02` | Should the dwell head shrink toward the context model? | Phase 3 | Fit `ψ` freely |
| `EXP-03` | Does per-timepoint ledger resolution change the design? | Phase 4 (v2 scope) | Do not build the two-phase ledger |
| `EXP-04` | Does the tier-3 burden proxy correlate with measured burden? | Phase 2 | Demote burden to `diagnostic` for tier-3 hosts |
| `EXP-05` | Lamarckian or Baldwinian repair? | Phase 2 (moved from Phase 1; `DEP-EXP-001`) | Keep the Lamarckian default |
| `EXP-06` | Does pause/boundary concordance predict solubility or activity? | Phase 2→3 promotion | Keep as `diagnostic`, or remove |
| `EXP-07` | How much should nearest-host fallback inflate uncertainty? | Phase 3 | Keep conservative defaults, `calibrated: false` |
| `EXP-08` | What window and interaction order per host? | Phase 3 | Ship `W = 10`, pairwise only |
| `EXP-09` | Does relaxed-gradient local search beat block coordinate ascent? | none (optional) | Do not build L2 search |
| `EXP-10` | Do the new layers improve *designs*, not just predictions? | Phase 5 | Publish the negative result; re-scope |
| `EXP-11` | What does surrogate rationing cost in final front quality? | Phase 4 | Raise `E_max`, accept the slowdown |

---

## EXP-01 — Cross-dataset and cross-protein generalization of the epistatic surface

**Question.** Does a surface fit by pseudo-likelihood on public Ribo-seq predict codon-level ribosome density for (a) genes held out from the same dataset, (b) a *different* Ribo-seq dataset for the same organism, and (c) heterologous proteins unlike anything in training?

**Why planning cannot answer it.** This is the central technical bet of the whole architecture (`RISK-EPI-001`). It is an empirical question about biology and about how much of the Ribo-seq signal is protocol artefact.

**Design.**
1. Fit on *E. coli* dataset A, evaluate on held-out genes from A (split **by gene**, never by codon).
2. Evaluate the same fitted model on dataset B — a different lab, different A-site offset rules, different drug treatment. This is the step that separates biology from batch effect.
3. Evaluate on a heterologous-protein subset: genes with low sequence and compositional similarity to the training set.
4. Repeat for *S. cerevisiae*.
5. Baselines, in increasing strength: codon usage only; tAI only; codon usage + position spline; the full model.

**Decision rule (pre-registered).**
- **(a) held-out genes:** Spearman improvement over the codon-usage baseline `≥ 0.05` with a bootstrap 95 % CI excluding zero → `GATE-EPI-001`.
- **(b) cross-dataset:** the improvement retains `≥ 50 %` of its within-dataset magnitude. Below that, the surface is substantially fitting protocol artefacts.
- **(c) heterologous:** the improvement retains `≥ 50 %` of the within-dataset magnitude on the dissimilar-protein subset.

**If it fails.** (b) failing → the surface is usable only when the target host's own Ribo-seq was processed identically; the term is restricted and prominently marked. (c) failing → the surface does not transfer to heterologous design, which is the main use case; demote `epistasis.dwell_deviation` to `diagnostic` and keep `epistasis.context_naturalness`, which does not depend on Ribo-seq at all.

**Gates:** Phase 3 exit. **Data:** ≥2 independent Ribo-seq datasets per organism, both reprocessed through the same pipeline.

---

## EXP-02 — Empirical-Bayes coupling between the dwell head and the context model

**Question.** Does shrinking the dwell couplings `ψ_d` toward a scaled version of the context couplings `J_d` improve held-out dwell prediction, compared with fitting `ψ` freely?

**Why planning cannot answer it.** It depends on how much "frequent in native genes" and "fast to translate" actually coincide in a given host — which is itself an interesting result and is not known a priori.

**Design.** Fit three variants on the same by-gene splits: (i) `ψ` free; (ii) `ψ_d = γ_d · J_d` (hard tie); (iii) `ψ_d ~ N(γ_d·J_d, τ_d²)` (shrinkage). Compare held-out NB log-likelihood per codon. Report the fitted `γ_d` — a per-offset measure of how much frequency and speed agree.

**Decision rule.** Adopt the shrinkage variant iff it beats the free fit by `≥ 0.002` nats/codon on held-out data with a bootstrap CI excluding zero, **and** the advantage grows as training data is subsampled (which is the regime the prior is meant to help).

**If it fails.** Fit `ψ` freely; drop the prior. Cost: data-poor hosts lose a stabilizer, so the `EXP-08` tiering thresholds tighten.

**Gates:** Phase 3.

---

## EXP-03 — Does per-timepoint ledger resolution change the design?

**Question.** Does a two-phase (pre/post-induction) resource ledger produce materially different designs from a single steady-state ledger?

**Why planning cannot answer it.** The question is not whether two-phase is *more accurate* — it obviously is — but whether the added accuracy changes the answer enough to justify requiring kinetic parameters most hosts do not have.

**Design.** For an inducible-system benchmark set (≥20 constructs), run the full pipeline with (i) the steady-state ledger and (ii) a two-phase ledger using literature parameters for *E. coli*. Compare the resulting Pareto fronts by IGD+ and by the fraction of codon positions that differ between matched candidates.

**Decision rule.** Build the two-phase extension iff median matched-candidate codon divergence is `> 5 %` **or** IGD+ between the two fronts exceeds the run-to-run variation from re-seeding alone. Otherwise the difference is within noise and the extension is not built.

**If it "fails"** (i.e. no material difference), that is a *good* outcome: it closes a scope question and removes a v2 feature. Record the negative result.

**Gates:** Phase 4 scoping.

---

## EXP-04 — Does the tier-3 burden proxy correlate with measured burden?

**Question.** For hosts with no proteomics and no tRNA abundance data, does the demand-profile mismatch score (capacity tier 3) predict measured growth cost at all?

**Why planning cannot answer it.** Tier 3 rests on the assumption that a host's tRNA pools are matched to its own native codon demand. That assumption is reasonable (it underlies tAI) but its *predictive* value for burden is untested.

**Design.** Take hosts where a tier-1 fit **is** possible. Compute both the tier-1 prediction and the tier-3 mismatch score, deliberately withholding the abundance/proteomics data from the tier-3 path. Correlate each against measured relative growth rate on a construct panel.

**Decision rule.** Keep `burden.growth_cost` as a Pareto axis for tier-3 hosts iff the tier-3 score achieves Spearman `≥ 0.35` against measured burden and retains `≥ 50 %` of the tier-1 correlation.

**If it fails.** Demote burden to `diagnostic` for tier-3 hosts — an outcome the `role` field already supports, so acting on it requires no redesign. Document that burden optimization requires at least tier-2 host data.

**Gates:** Phase 2 exit. **Data:** a published burden/expression panel (Codon-Harmony-Index-style) for at least one organism.

---

## EXP-05 — Lamarckian vs Baldwinian repair

**Question.** Does replacing the genotype with its repaired form (Lamarckian) beat evaluating the repaired form while inheriting the original (Baldwinian)?

**Design.** Run the standard benchmark suite under both policies, 20 seeds each, with tightly-constrained inputs where repair fires often. Compare final hypervolume, population diversity (mean pairwise Hamming distance), and total repair invocations.

**Decision rule.** Switch the default to Baldwinian iff it wins on hypervolume by more than the seed-to-seed standard deviation **and** its repair-invocation overhead is under 2×.

**Unchanged, deliberately.** The rule turns on hypervolume and on population diversity, both of which need a population — so this experiment could not run in Phase 1 and `EXIT-003` moved to Phase 2 rather than the rule being restated in Phase-1-measurable terms. A proxy over single-objective score under a fixed budget would answer a different question, and a default chosen from it would have to be re-decided once `P2-02` landed. `DEP-EXP-001`.

**If it fails.** Keep Lamarckian (the current default). Cheap experiment, cheap either way.

**Gates:** Phase 1 exit (it needs only the constraint layer and a single-objective search).

---

## EXP-06 — Does pause/boundary concordance predict solubility or activity?

**Question.** Does `kinetics.pause_concordance` associate with measured soluble fraction or specific activity, controlling for total protein output?

**Why planning cannot answer it.** This is the falsification test for the weakest layer in the system (`RISK-KIN-001`, `AMD-002`).

**Design.** Assemble panels for ≥3 protein families with published or newly-measured (soluble fraction, total expression, sequence) triples across synonymous variants. Fit `solubility ~ concordance + total_expression + family` as a mixed model with family as a random effect. Pre-register the model.

**Decision rule.** Promote the term from `diagnostic` to `pareto_axis` iff the concordance coefficient is significant at `α = 0.05` after multiple-comparison correction, has the predicted sign, and generalizes to a held-out family.

**If it fails.** Keep as `diagnostic` if the effect is directionally right but weak; **remove from the default objective set** if there is no effect anywhere. Publishing the negative result is valuable in itself: it is a clean test of a widely-repeated claim.

**Gates:** the Phase 2 → Phase 3 promotion of the kinetics term.

---

## EXP-07 — Calibrating nearest-host uncertainty inflation

**Question.** By how much must a fallback model's predictive interval widen to retain its nominal coverage when applied to a different host?

**Why planning cannot answer it.** `κ_tax` and `κ_cu` in [`../04-api/04-model-resolution.md`](../04-api/04-model-resolution.md) §3.3 are exactly the kind of parameter it would be dishonest to guess and then present as calibrated.

**Design.** For every ordered pair of hosts where both have held-out data, evaluate host `h`'s model on host `h'`'s data. Record the empirical coverage of the nominal 90 % interval, the taxonomic LCA rank, and the RSCU cosine. Fit `κ_tax`, `κ_cu` by minimizing the deviation of empirical coverage from 90 % across pairs.

**Decision rule.** Adopt the fitted factors iff, on held-out pairs, empirical coverage lands within `[0.87, 0.93]`. If no factor achieves this, tighten the fallback admissibility threshold instead — reducing the set of permitted fallbacks rather than pretending the remainder are calibrated.

**If it fails.** Keep the conservative defaults and leave `Uncertainty.calibrated = false` for every fallback, which already depresses those candidates' robustness scores automatically.

**Gates:** Phase 3. **Data:** ≥4 hosts with fitted surfaces and held-out sets.

---

## EXP-08 — Per-host window size and interaction order

**Question.** What window `W` and interaction order does each host's data actually support?

**Design.** Run the selection procedure of [ADR-0003](../07-decisions/ADR-0003-epistasis-window-and-order.md) §1 over the grid `{3, 5, 10, 15, 20}` per host with by-gene splits; report held-out pseudo-likelihood curves and the group-lasso-selected offset support. Then run the triplet selection and gating procedure and report whether the gate ever passes.

**Decision rule.** Adopt the one-standard-error `W*` per host. Adopt triplets only where the held-out gate passes (`≥ δ` improvement with a bootstrap CI excluding zero).

**If it fails to discriminate** (the held-out curve is flat across the whole grid), ship `W = 10` pairwise and record that the data does not support choosing — which is itself a publishable, useful finding.

**Gates:** Phase 3.

---

## EXP-09 — Relaxed-gradient local search vs. block coordinate ascent

**Question.** Does optimizing continuous codon logits with gradients and rounding beat discrete block coordinate ascent, at equal compute?

**Why it is registered despite [ADR-0002](../07-decisions/ADR-0002-oracle-vs-differentiable-surrogate.md) rejecting L2 as a requirement.** The ADR predicts gradient-based search will lose, mainly because rounding discards constraint feasibility. That prediction should be testable rather than assumed.

**Design.** Implement a throwaway relaxed optimizer for the differentiable subset of terms; compare against tier 2 at matched wall-clock on the benchmark suite; measure objective value **after** rounding and repair.

**Decision rule.** Add an L2-based tier only if it beats tier 2 by more than the seed-to-seed standard deviation at equal wall time. Expected outcome: it does not.

**Gates:** none — this is optional work, deliberately unblocking.

---

## EXP-10 — Do the new layers improve designs, not just predictions?

**Question.** The end-to-end question the whole project rests on: do sequences designed with the epistasis, burden and kinetics layers outperform sequences designed by CAI/MFE optimization alone, when actually built and measured?

**Why planning cannot answer it.** Every other experiment here validates a *component*. This one validates the *thesis*. All the component tests could pass while the thesis fails — for instance if the improvements are real but small relative to variation in expression from RBS context and plasmid copy number.

**Design.** A blinded head-to-head panel: for each of ≥8 target proteins, build and measure (i) the CodonForge full-objective design, (ii) a CAI-maximized design, (iii) a LinearDesign output, (iv) the native/wild-type sequence, (v) a CodonForge design with the epistasis and burden layers *ablated*. Measure expression, growth cost and soluble fraction with adequate replication. Randomize plate positions; blind the assay operator to the design arm.

**Decision rule.** Pre-register the primary endpoint (soluble active protein per unit growth cost) and a target effect size. Arm (v) — the ablation — is what isolates the new layers' contribution from the rest of the pipeline, and it is what makes the result interpretable rather than merely positive.

**If it fails.** Publish it. A rigorous negative result on whether context-aware and burden-aware codon design beats CAI optimization is a genuine contribution, and the tool remains useful as a constraint-aware multi-objective optimizer regardless.

**Gates:** Phase 5. This is the project's capstone validation, not a build blocker.

---

## EXP-11 — The cost of surrogate rationing

**Question.** How much final front quality is lost by evaluating only `E_max` individuals per generation on the L0 terms?

**Design.** Run the benchmark suite with `E_max ∈ {16, 32, 64, 128, all}` and compare final hypervolume and wall time. Also measure how often front certification ([`../03-algorithms/03-nsga-loop-and-operators.md`](../03-algorithms/03-nsga-loop-and-operators.md) §6) changes the emitted front, which is the direct symptom of over-aggressive rationing.

**Decision rule.** Set the per-tier `E_max` at the smallest value whose hypervolume is within 2 % of the `all` condition. If no finite value achieves that, the surrogate is inadequate and needs better features before `standard` tier can rely on it.

**Gates:** Phase 4 (`GATE-PERF-001` depends on the chosen `E_max`).
