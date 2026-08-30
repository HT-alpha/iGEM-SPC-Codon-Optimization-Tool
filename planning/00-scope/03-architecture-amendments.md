# Architecture Amendments — Where Detailed Planning Found the Design Impractical or Underspecified

The brief's closing note asks for exactly this document:

> If any part of the CodonForge architecture proposal turns out, under detailed planning, to be impractical or underspecified in a way this brief did not anticipate, Opus 5 Max should say so explicitly and propose a specific amendment rather than silently working around it.

Twelve items. Each states what the architecture says, why detailed planning found it a problem, the specific proposed amendment, and what it costs to ignore. Everything here is already reflected in the plan; this document is the consolidated, reviewable list.

Smaller, local deviations are flagged inline throughout the package as `DEP-*` markers with a one-line rationale, per `BRIEF-4.b`. This file is for the ones a reviewer should decide on deliberately.

| ID | Severity | One-line summary |
| --- | --- | --- |
| `AMD-001` | **Blocking** | Global MFE folding cannot run inside the search loop; it alone costs ~53× the whole job's compute budget. |
| `AMD-002` | **Blocking** | The folding-kinetics layer has no validated predictor and must not ship as a Pareto axis. |
| `AMD-003` | Significant | "Per-position confidence" is not one number and cannot be defined as one. |
| `AMD-004` | **Blocking** | The active-learning loop is a closed feedback loop with a known failure mode. |
| `AMD-005` | Editorial | "Epistasis" here is not protein-level epistasis; the terminology will cost credibility. |
| `AMD-006` | Significant | Ribo-seq exists for a handful of organisms; the layer must be tiered, and the coverage claim narrowed. |
| `AMD-007` | Significant | Rare-codon bans and CAI maximization contradict this architecture's own premise. |
| `AMD-008` | Significant | Five-plus objectives make nearly everything non-dominated; the front is unpresentable without reduction. |
| `AMD-009` | **Blocking** | The obvious structure-prediction dependencies are incompatible with an open-source release. |
| `AMD-010` | Significant | Constraint sets can be jointly unsatisfiable; without analysis the tool returns unexplainable nothing. |
| `AMD-011` | Significant | Determinism under parallelism is not free, and the transparency goal depends on it. |
| `AMD-012` | Editorial | The "loop" in active learning closes on a timescale of months, not minutes. |

---

## AMD-001 — Global MFE folding cannot sit inside the evolutionary loop

**Severity:** Blocking.

**What the architecture implies.** mRNA secondary structure is one of the biophysical objectives, alongside LinearDesign and CDSfold as reference points — both of which optimize a global folding energy.

**Why it is impractical.** Zuker-class MFE folding is `O(n³)` time and `O(n²)` memory. At `n = 2 400` nt that is ~4 s per evaluation; a `standard`-tier job performs ~`6.4 × 10⁴` evaluations, giving **~71 CPU-hours (2.6 × 10⁵ CPU-seconds) for a single 2 kb gene**. The whole job's budget is 600 wall-seconds on 8 cores — 4 800 CPU-seconds — so this one objective is ~53× over the entire budget before anything else is computed. At 15 kb it is ~20 minutes *per evaluation* and ~3.5 GB of memory.

**Proposed amendment.** Make **windowed local accessibility** (RNAplfold-class, `w = 80` nt, `O(n·w²)`, ~13 ms at 2 kb, incrementally updatable in ~0.4 ms per edit) the default structure objective. Reserve global MFE for `exact` tier and for final re-ranking of the ≤20-candidate reduced front, where 20 × 4 s is trivially affordable.

**This is also the better-targeted objective.** The mechanistic claims the structure term serves — accessibility of the ribosome binding region, local unfolding ahead of the elongating ribosome — are *local* claims. A single global MFE for a 15 kb polycistronic mRNA is a poor proxy for either. The trade-off is real, though: a design optimized for local accessibility is not optimized for global thermodynamic stability, which matters for mRNA-therapeutic applications. `exact` tier exists for that case.

**Cost of not amending.** The performance budgets in `BRIEF-2.3-g` are unreachable and the workstation target of [ADR-0008](../07-decisions/ADR-0008-compute-budget-and-scaling.md) is unreachable, so the tool becomes cluster-only.

**Where reflected:** [`../03-algorithms/08-complexity-and-performance-budgets.md`](../03-algorithms/08-complexity-and-performance-budgets.md) §3, [ADR-0008](../07-decisions/ADR-0008-compute-budget-and-scaling.md).

---

## AMD-002 — The folding-kinetics objective has no validated predictor and must not ship as a Pareto axis

**Severity:** Blocking.

**What the architecture says.** Co-translational folding kinetics is one of five co-equal objectives.

**Why it is a problem.** Of the four capabilities the architecture introduces without a reference implementation, three have a clear statistical object and a clear data source. This one has neither an established predictor nor a routinely-measured readout at the scale needed to train one. Shipping it as a Pareto axis on day one means an unvalidated term steers every design — and, worse, makes its own validation circular: designs are selected partly by the score, then evaluated for whether the score predicts anything.

**Proposed amendment.** A proxy ladder with an explicit promotion gate:

1. **v1** ships a pause/domain-boundary concordance proxy as a **`diagnostic`** term — computed, reported, not steering the search.
2. Promotion to a Pareto axis requires `EXP-06`: a significant association between the score and measured solubility/activity across ≥3 protein families, controlling for total expression.
3. A kinetic simulation (v2) is built only if the proxy is promoted *and then* found to leave large systematic residuals.
4. If `EXP-06` finds no association anywhere, the term is removed from the default objective set.

The `role: diagnostic` field in the job request exists precisely to make step 1 expressible.

**Cost of not amending.** An unvalidated objective distorts every design, and the architecture's strongest claim — that it models things existing tools ignore — becomes its weakest published result.

**Where reflected:** [`../03-algorithms/07-folding-kinetics-objective.md`](../03-algorithms/07-folding-kinetics-objective.md), `RISK-KIN-001`, `EXP-06`.

---

## AMD-003 — "Per-position confidence" cannot be a single number

**Severity:** Significant.

**What the architecture says.** The decision-support layer promises per-position confidence scores.

**Why it is underspecified.** The objective vector mixes terms with structurally different uncertainty: a penalized-MLE Potts model (analytic parameter covariance), an ensemble-based dwell head (model-form uncertainty), a deterministic CAI (no uncertainty at all), a burden term whose uncertainty lives in its *inputs*, and a transformer. There is no meaningful scalar that summarizes "confidence at codon 214" across those. A single number would either be dominated by whichever term happens to have the largest variance, or be an unweighted average of incommensurable quantities.

There is a second, subtler issue: split conformal gives **marginal** coverage — 90 % averaged over positions — which is the wrong guarantee when the question is "is *this region* trustworthy".

**Proposed amendment.** Define per-position confidence as a **per-term composite** with stated semantics: each term contributes its own calibrated per-position interval, and the reported `confidence` array is accompanied by its `term_id`. Use **Mondrian (class-conditional) conformal**, conditioning on position class and local training support, so the coverage guarantee is conditional rather than marginal. The user-facing `low_confidence_region` tag fires on the composite, but the underlying per-term values remain inspectable.

**Cost of not amending.** A confidence number nobody can interpret, which users will learn to ignore — defeating the purpose of the decision-support layer.

**Where reflected:** [ADR-0005](../07-decisions/ADR-0005-uncertainty-and-pareto-ranking.md), [`../02-schemas/landscape_response.schema.json`](../02-schemas/landscape_response.schema.json).

---

## AMD-004 — The active-learning loop is a closed feedback loop and needs guards

**Severity:** Blocking.

**What the architecture says.** Wet-lab results feed back to recalibrate the models.

**Why it is a problem.** If the constructs sent to the bench are the model's own top-ranked suggestions, and the resulting measurements retrain that model, the loop is closed. The model becomes more confident in the region it already liked and less informed everywhere else. Its held-out metrics can *improve* while its usefulness for new proteins degrades — so the failure is invisible to the obvious monitoring. This is the default outcome of naive active learning, not an edge case.

**Proposed amendment.** Three mandatory guards:

1. **Exploration arms.** `observation.exploration` is a required field. The recommended protocol is ≥25 % of each submitted batch drawn by `uniform_random_from_front` or `max_uncertainty` rather than `top_ranked`. Batches that are entirely top-ranked are usable but downweighted and flagged.
2. **Selection-rule-aware fitting.** The selection rule enters the recalibration model as a covariate, with inverse-propensity weighting where the rule is known.
3. **A frozen holdout that recalibration never trains on and never selects from**, which is where `GATE-ACT-001` is evaluated. This is the only place the closed-loop degradation becomes visible.

**Cost of not amending.** A model that gets measurably better and practically worse, with no signal that it is happening.

**Where reflected:** [`../04-api/03-active-learning-ingestion.md`](../04-api/03-active-learning-ingestion.md) §5, [`../02-schemas/observation.schema.json`](../02-schemas/observation.schema.json).

---

## AMD-005 — "Epistasis" is the wrong word for what is being modelled

**Severity:** Editorial, but it will cost credibility with exactly the reviewers the project wants.

**What the architecture says.** An "epistatic fitness surface", fit by DCA/Potts methods borrowed from protein contact prediction.

**Why it is a problem.** In molecular genetics, *epistasis* normally means non-additive interaction between mutations at the protein or organism level. What CodonForge models is **context-dependence of synonymous codon effects on translation**: codon-pair effects, ribosome footprint interactions, local decoding context. That is statistical epistasis in the fitness-landscape sense, but a reader from protein science will initially expect something else — and the DCA/Potts borrowing reinforces the wrong reading, because DCA in its home domain infers *structural* contacts.

**Proposed amendment.** Keep `epistasis` as the internal namespace (renaming it now would churn every schema), but in all user-facing and published material use **"context-dependent codon effects"** or **"translational context model"**, with a single sentence disambiguating on first use. Add a glossary entry. The term ids already read well under this reading: `epistasis.context_naturalness`, `epistasis.dwell_deviation`.

**Cost of not amending.** Reviewers spend their first pass resolving a terminology confusion instead of evaluating the method, and the paper invites an objection that has nothing to do with the work.

**Where reflected:** [`02-glossary-and-notation.md`](02-glossary-and-notation.md).

---

## AMD-006 — Ribosome profiling exists for a handful of organisms; the coverage claim needs narrowing

**Severity:** Significant.

**What the architecture implies.** A host-specific epistatic surface fit from ribosome-profiling data.

**Why it is a problem.** Good Ribo-seq exists for *E. coli*, *S. cerevisiae*, human cell lines, and a modest number of others. It does not exist for most hosts a synthetic-biology user cares about — *B. subtilis* strains, *P. pastoris*, *Y. lipolytica*, non-model bacteria, most plants and insects. An architecture whose flagship layer requires Ribo-seq is an architecture for four organisms, and stating otherwise sets up an unmeetable expectation.

**Proposed amendment.** Split the layer, and tier it:

- **A context model** fit by pseudo-likelihood on the host's **native CDS set** — needs only an annotated genome, so *every* host gets a real context-aware term.
- **A dwell head** fit on Ribo-seq, offered only above a data threshold and simply **absent** below it, rather than fit badly.
- Data-poor hosts shrink toward a **pooled cross-host prior** (justified because ribosome geometry and wobble decoding are far more conserved than codon usage) and reduce their window.
- `HostProfile.data_completeness.supported_terms` becomes the machine-readable contract, and the API **refuses** a job requesting an unsupported term rather than silently degrading.

And in the architecture document: state the coverage claim as *"context-aware for any annotated genome; dwell-calibrated for hosts with ribosome profiling"*.

**Cost of not amending.** Either the tool works for four organisms, or it silently fabricates a dwell model for hosts with no data — the second being much worse.

**Where reflected:** [ADR-0003](../07-decisions/ADR-0003-epistasis-window-and-order.md) §2, [`../03-algorithms/01-epistatic-surface-fit.md`](../03-algorithms/01-epistatic-surface-fit.md) §1.

---

## AMD-007 — Rare-codon bans and CAI maximization contradict the architecture's own premise

**Severity:** Significant.

**What the architecture implies.** CAI is one of the benchmark axes (against LinearDesign/CDSfold), and users will expect the familiar "ban rare codons" control.

**Why it is a problem.** The architecture's central argument is that codon-frequency-based optimization is inadequate — that context, pausing and burden matter. But a **blocking** rare-codon ban forecloses exactly the designs the epistasis and kinetics layers exist to find: a rare codon retained to preserve a folding pause, or to break an unfavourable context pair. The tool would be fighting itself, and the constraint would win because blocking constraints always win.

**Proposed amendment.**

- `codon_ban` and `min_codon_frequency` default to **`advisory`** class, contributing to the objective rather than to feasibility.
- `usage.cai` remains a first-class Pareto axis — it is what the benchmarks compare on and it is what users understand — but as an *axis*, so trading it away is visible and deliberate rather than forbidden.
- When a rare codon is retained deliberately, the `rare_codon_retained` rationale tag says so in words, with the term that would degrade if it were replaced.

**Cost of not amending.** A default configuration in which the new layers cannot express their conclusions, and benchmark results that look like a slightly worse CAI optimizer.

**Where reflected:** [`../02-schemas/constraint_set.schema.json`](../02-schemas/constraint_set.schema.json), [`../02-schemas/rationale_tags.json`](../02-schemas/rationale_tags.json).

---

## AMD-008 — Five-plus objectives make nearly everything non-dominated

**Severity:** Significant. Presentational, not computational — which is why it is easy to miss until users see it.

**What the architecture says.** Five objectives (expression/epistasis, structure, folding kinetics, LM naturalness, burden) with a Pareto front as the deliverable.

**Why it is a problem.** In `M` dimensions the probability that a random candidate is dominated falls off sharply with `M`. At `M = 5` a converged population is essentially *all* non-dominated: a 212-member population routinely yields 130–200 non-dominated candidates. The search handles this fine. The **user** does not: "here are 137 equally-good options" is not decision support, it is the absence of it.

**Proposed amendment.** Separate the *search* objective set from the *decision* objective set:

1. `job_request.objectives[].role` allows each term to be a `pareto_axis`, an `epsilon_constraint` (a bound, not an axis), or a `diagnostic` (computed and reported, not optimized). Demoting terms is the standard remedy for many-objective front explosion.
2. Default to **≤4 Pareto axes**; a longer list is accepted with a warning.
3. The delivered front is **reduced** — reference-direction representatives by default, with knee selection, clustering and ε-grid available — and `ParetoFront.reduction` records the raw count, the presented count, the method and a link to the full set.
4. Ranking within a niche uses expected hypervolume contribution under the predictive distribution ([ADR-0005](../07-decisions/ADR-0005-uncertainty-and-pareto-ranking.md)), so what survives reduction is what is *reliably* good.

**Cost of not amending.** The decision-support layer delivers an undifferentiated list, and users fall back to sorting by one column — which is the behaviour multi-objective design was supposed to replace.

**Where reflected:** [`../02-schemas/pareto_front.schema.json`](../02-schemas/pareto_front.schema.json), [`../03-algorithms/03-nsga-loop-and-operators.md`](../03-algorithms/03-nsga-loop-and-operators.md) §2.

---

## AMD-009 — The obvious structure-prediction dependencies are incompatible with an open-source release

**Severity:** Blocking for the open-source goal. This is a legal constraint, not a technical one, and it is the item most likely to be discovered too late.

**What the architecture implies.** Standard tooling for RNA folding, with LinearDesign and CDSfold as comparison points.

**Why it is a problem.** Verified licence terms:

- **ViennaRNA** — permits research, educational and commercial *use and modification*, but requires that "the package and any derived works are not redistributed for any fee, other than media costs". That is not an OSI-approved licence and is incompatible with redistributing a combined work under Apache-2.0 or GPL.
- **LinearFold** — the same licence text, with the same no-fee-redistribution clause. The ecosystem generally treats it as non-commercial out of caution (VaxPress ships it only in a `[nonfree]` extra).
- **LinearDesign** — licence listed as "Other", and **Baidu Research has filed a patent on the LinearDesign algorithm**. Patent encumbrance is a distinct risk from licensing and is not cured by a permissive licence.

Vendoring any of these makes CodonForge's own distribution non-redistributable-for-fee, which forecloses commercial support, paid hosting, and inclusion in most distributions.

**Proposed amendment.**

1. Ship CodonForge under **Apache-2.0** (its explicit patent grant matters given the LinearDesign filing).
2. Define a **`StructureEngine` plugin seam** and ship a **permissively-licensed windowed-accessibility kernel as the default** — which `AMD-001` already makes the preferred engine on performance grounds, so the two amendments reinforce each other.
3. ViennaRNA and LinearFold back ends are **feature-gated and user-installed**, detected at runtime, never vendored. There is deliberately **no `pip install codonforge[nonfree]` extra**, because an extra that pulls them in would make CodonForge's wheel non-redistributable.
4. Benchmarks against LinearDesign/CDSfold run against **user-installed** binaries or against **published output**, via the seed-import path of [ADR-0007](../07-decisions/ADR-0007-external-tool-seed-import.md). CI runs them in a separate, non-distributed job.
5. Add a licence-and-patent review to the Phase 2 exit checklist.
6. **Thermodynamic parameters are input data, not shipped constants.** `DEP-STRU-001`. This clause was missing and is the part of `AMD-009` that was in fact discovered late: the amendment verifies the licence position of the folding *code* and says nothing about the nearest-neighbour parameter set that code needs to compute anything. Those parameters have the same class of problem — the published sets are embedded in tools whose licences are exactly what this amendment rejects — so a kernel with a table baked in is either carrying invented values, which the no-fabricated-constants rule forbids, or redistributing data whose terms nobody here has verified. The kernel therefore ships with **no** parameter table: a `structure_parameters` document supplies one, carries its own provenance and citation, and its digest is recorded in the landscape lockfile like any other artifact. A run with no parameter set is refused with a named error rather than falling back to a plausible default.

**Cost of not amending.** A tool that cannot legally be distributed as claimed, discovered at release rather than at design time.

**Where reflected:** [ADR-0001](../07-decisions/ADR-0001-language-and-runtime-split.md) (licensing note), `RISK-DEP-001`, `RISK-DEP-002`, [`../04-api/02-cli-and-library.md`](../04-api/02-cli-and-library.md) §5.

---

## AMD-010 — Constraint sets can be jointly unsatisfiable, and the tool must say which one to relax

**Severity:** Significant.

**What the architecture says.** A hard-constraint filter and an in-loop repair operator.

**Why it is underspecified.** Nothing in the architecture addresses what happens when the constraint set has no solution — a GC band no synonym assignment can reach, a k-mer limit that is unsatisfiable by pigeonhole, a blocking violation inside immutable backbone, or motif constraints that jointly close off every continuation. This is not rare: users combine an assembly standard, a vendor profile and their own additions without checking they are compatible.

The existing prototype in this repository shows the failure mode exactly. After exhausting its restart budget it prints `Error: Exceeded max restarts (1000) without finding a solution.` The user learns nothing about why, or what to change.

**Proposed amendment.** A **feasibility analyser** run at submission and exposed as `POST /v1/feasibility`, with five checks — locked-span violations, the k-mer counting bound, per-position exhaustion, automaton reachability, GC-window attainability — each of which reports **a constraint to relax and by how much**, not merely a failure. Checks 2 and 5 are exact; 1, 3 and 4 are necessary conditions, and the report says so.

**Cost of not amending.** The most common user-facing failure produces an unexplainable empty result, which is how a codon optimizer becomes shelfware.

**Where reflected:** [`../03-algorithms/06-constraint-repair.md`](../03-algorithms/06-constraint-repair.md) §5, [`../04-api/openapi.yaml`](../04-api/openapi.yaml) (`/v1/feasibility`).

---

## AMD-011 — Determinism under parallelism is not free, and the transparency goal depends on it

**Severity:** Significant.

**What the architecture says.** Transparency, auditability and reproducibility are explicit goals.

**Why it is underspecified.** A parallel evolutionary search is reproducible only if two things hold, and neither is automatic: RNG streams must be assigned by a *thread-count-independent* index (per block, per individual — not per worker), and cross-thread floating-point reductions must use a fixed order, since floating-point addition is not associative. Without both, "same seed, different machine" is a different answer, and every reproducibility claim in the plan is void.

**Proposed amendment.** Make determinism a contract, not an aspiration: `LandscapeLockfile.rng` records the algorithm, root seed, stream-assignment scheme and a `deterministic_reduction` flag. `per_worker` stream assignment is permitted only at `fast` tier. Deterministic reduction costs ~3 % throughput and is required for any published result. `TST-SRCH-004` asserts bit-identical output across 1, 4 and 16 threads.

**Cost of not amending.** Results that cannot be reproduced across machines, discovered when a reviewer tries.

**Where reflected:** [`../02-schemas/landscape_lockfile.schema.json`](../02-schemas/landscape_lockfile.schema.json), [`../03-algorithms/02-block-local-coordinate-ascent.md`](../03-algorithms/02-block-local-coordinate-ascent.md) §6.

---

## AMD-012 — The active-learning loop closes on a timescale of months

**Severity:** Editorial, but it should shape expectations and the roadmap.

**What the architecture implies.** An active-learning feedback loop as a core capability.

**Why it matters.** The loop's latency is set by biology: design, order synthesis, clone, transform, grow, assay. That is weeks to months per iteration, and a lab produces on the order of tens of measurements per round. This is a *scheduled batch recalibration* regime, not an online-learning one — and the design consequences are real: adapters rather than full retrains, hard minimum-evidence gates (`GATE-ACT-005` requires ≥12 observations from ≥2 batches), and acceptance gating rather than continuous updating.

**Proposed amendment.** Describe the layer as **scheduled batch recalibration from wet-lab feedback**, not online active learning. State the expected cadence. Position its value correctly: the highest-value recalibration is not a marginally better epistatic surface, it is **moving a host from capacity tier 3 or 2 to tier 1** ([ADR-0004](../07-decisions/ADR-0004-resource-ledger-granularity.md)), which is otherwise unreachable and which changes the burden term from a proxy into a prediction. That should be the loop's headline claim.

**Cost of not amending.** An expectation of a fast-improving system, and a roadmap that under-weights the one recalibration that actually transforms a host's capability.

**Where reflected:** [`../04-api/03-active-learning-ingestion.md`](../04-api/03-active-learning-ingestion.md) §4 and §7, [ADR-0006](../07-decisions/ADR-0006-active-learning-data-isolation.md).

---

## What was checked and found sound

For balance, the parts of the architecture that detailed planning did **not** find problems with:

- **The layer decomposition itself.** Eight layers map cleanly onto eight modules with an acyclic dependency graph, once two support leaves are added to break the latent cycles.
- **The three-tier search structure.** Cheap seeding → local exploitation → population-based multi-objective exploration is well matched to the cost structure, and the tiers turn out to correspond exactly to the L0/L1 capability split.
- **The swappable-layer goal.** Achievable, and the mechanism (a runtime term registry behind an interface crate) is cheap.
- **Pseudo-likelihood for the context model.** Better than expected: restricting the partition function to synonymous sets makes the objective concave with a unique global optimum, gives `O(W·S)` exact deltas, and keeps the parameter count at ~41 k.
- **Systems-level joint search.** The cross-cistron coupling is entirely through a ~40-number vector, which makes price decomposition natural and reduces the joint problem's wall time to roughly that of a single cistron.
- **Hard-constraint repair inside the loop.** Workable, provided it is bounded, minimal, randomized and tabu-guarded — all four of which are cheap.
