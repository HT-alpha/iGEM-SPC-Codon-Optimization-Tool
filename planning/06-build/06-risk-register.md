# Risk Register

Satisfies `BRIEF-2.6-b` and acceptance criterion `BRIEF-5.11`:

> The risk register identifies, for each of the four previously-unaddressed capabilities (epistasis, systems-level design, burden modeling, folding-kinetics/functional feedback), the specific technical bet being made and what would falsify it.

Every entry names **the bet** (a falsifiable claim), **the falsifier** (what observation would show it wrong), **the mitigation** (what reduces the chance), and **the fallback** (what the system does if the bet loses). A risk without a fallback is a hope.

---

## 1. The four technical bets

### `RISK-EPI-001` — Epistatic surfaces generalize from public Ribo-seq to arbitrary target proteins

| | |
| --- | --- |
| **Capability** | Epistatic fitness surface |
| **The bet** | A synonymous-restricted Potts model plus a dwell head, fit on one organism's public ribosome-profiling data, predicts context-dependent ribosome density for *heterologous* proteins never seen in training, well enough to change design decisions. |
| **Why it might be wrong** | Ribo-seq measures ribosome *density*, which mixes elongation kinetics with initiation rate, drug artefacts, ligation bias and mapping bias. Much of the apparent context signal may be protocol artefact, and artefacts do not transfer across datasets or proteins. |
| **Falsifier** | `EXP-01`(b) or (c): cross-dataset or heterologous-protein Δρ retains less than 50 % of its within-dataset magnitude. |
| **Likelihood / impact** | Medium / **High** — this is the architecture's headline claim. |
| **Mitigation** | Nuisance covariates for the known artefacts (`dG_local`, tunnel charge, proline, internal SD, batch, drug); by-gene splits; two independent datasets reprocessed identically; a context model that needs no Ribo-seq at all. |
| **Fallback** | Demote `epistasis.dwell_deviation` to `diagnostic`; keep `epistasis.context_naturalness`, which is fit on native CDSs and does not depend on Ribo-seq. The layer degrades to a context-aware naturalness term rather than disappearing. |
| **Phase** | 3 |

### `RISK-BURD-001` — Host capacity is estimable without host-specific proteomics

| | |
| --- | --- |
| **Capability** | Host resource/burden model |
| **The bet** | For hosts with no proteomics and no tRNA abundance data, a demand-profile mismatch score computed from the genome alone (capacity tier 3) carries enough signal to rank designs by burden. |
| **Why it might be wrong** | Tier 3 assumes tRNA pools are matched to native codon demand. The assumption underlies tAI and is reasonable *on average*, but burden is a bottleneck phenomenon driven by the tail of the distribution, and an average-case assumption may say nothing about the tail. |
| **Falsifier** | `EXP-04`: the tier-3 score achieves Spearman < 0.35 against measured burden, or retains less than half of the tier-1 correlation. |
| **Likelihood / impact** | Medium / Medium |
| **Mitigation** | Three explicit tiers; tier-3 output labelled as a mismatch score with a different unit and `calibrated: false`; `TST-BURD-005` prevents that labelling from eroding. |
| **Fallback** | Demote burden to `diagnostic` for tier-3 hosts. Burden optimization then requires at least tier-2 host data, and the API says so. |
| **Phase** | 2 |

### `RISK-SYS-001` — Cross-cistron resource contention is material and predictable

| | |
| --- | --- |
| **Capability** | Systems-level joint search |
| **The bet** | In real multi-cistron constructs, cistrons compete measurably for shared isoacceptor pools and ribosomes, and jointly-optimized designs outperform independently-optimized ones by enough to justify the layer. |
| **Why it might be wrong** | A heterologous construct's demand may be small relative to the native proteome's, in which case contention is negligible and joint optimization converges to the independent solution. The layer would then be correct and useless. |
| **Falsifier** | `bench_joint_search` on `bm_constructs_20`: joint and independent designs differ in under 2 % of codon positions and their hypervolumes are within seed-to-seed noise. |
| **Likelihood / impact** | Medium / Medium |
| **Mitigation** | `SYN-03` proves the machinery works where contention *is* material, separating "the layer is broken" from "contention is small here"; plasmid copy number and expression priors are first-class inputs, since high-copy strong-promoter constructs are where contention should appear. |
| **Fallback** | Keep the ledger as a *reporting* capability — burden attribution per cistron is useful even when joint optimization changes nothing — and disable joint search by default, enabling it only above a computed contention threshold. |
| **Phase** | 4 |

### `RISK-KIN-001` — Pause architecture predicts functional protein yield

| | |
| --- | --- |
| **Capability** | Folding-kinetics / functional feedback |
| **The bet** | Concordance between predicted translational pauses and domain boundaries associates with soluble or active protein fraction, controlling for total expression. |
| **Why it might be wrong** | The evidence base is a set of individually-compelling case studies, not a demonstrated general effect. The effect may be real but small relative to expression-level variation, or confined to particular protein classes. |
| **Falsifier** | `EXP-06`: no significant association across ≥3 protein families after correction, or no generalization to a held-out family. |
| **Likelihood / impact** | **High** / Medium |
| **Mitigation** | Ships as a `diagnostic` term, so it does not steer designs while unvalidated (`AMD-002`); saturating reward per boundary so it cannot be gamed; explicit reconciliation with the epistasis pause penalty. |
| **Fallback** | Remain a diagnostic, or be removed from the default objective set. A rigorous negative result is published — it is a clean test of a widely-repeated claim. |
| **Phase** | 2 (build) → 5 (validate) |

---

## 2. Scientific and modelling risks

| Id | Risk | Falsifier / trigger | Mitigation | Fallback |
| --- | --- | --- | --- | --- |
| `RISK-EPI-002` | The Potts model class is too weak; a long-context attention model would do materially better | An attention baseline beats the Potts held-out likelihood by a wide margin in `EXP-08`'s reporting | The contract is capability-based, so a different model can serve the same term | Swap the provider. Costs the concave objective (unique optimum) and `O(W·S)` deltas, so tier 2 needs a new performance plan — hence it is a risk, not a free upgrade |
| `RISK-BURD-002` | Per-cistron expression rates `f_k` are unknown, so the ledger's absolute scale is wrong | `EXP-03`/`bench_joint_search` show results are highly sensitive to `f_k` | `expression_prior_source` is recorded; RBS-strength priors; sensitivity analysis in the benchmark | Report burden as a *relative* ranking across candidates for a fixed `f_k`, which is scale-invariant, rather than as an absolute cost |
| `RISK-BURD-003` | Real burden runs through channels the ledger does not model (chaperones, secretion, membrane insertion, metabolic load) | Residual structure in `EXP-04` correlated with protein class rather than with codon composition | Documented scope limit; protein class recorded in observations | Add channels in v2 if a data source appears; otherwise state the limit plainly |
| `RISK-DATA-001` | Ribo-seq batch effects dominate biological signal | `GATE-EPI-002` fails | Processing metadata mandatory in `RiboSeqTrack`; drug and offset method as explicit covariates; identical reprocessing of all tracks | Restrict artifacts to same-pipeline use and mark them; publish the finding |
| `RISK-SCI-001` | The thesis fails: designs from the full pipeline do not outperform CAI/MFE optimization when built and measured | `EXP-10` primary endpoint not met | Ablation arms isolate each layer's contribution; blinding; pre-registration | Publish the negative result. The tool remains a strong constraint-aware multi-objective optimizer with reproducibility guarantees, which is independently useful |

---

## 3. Engineering risks

| Id | Risk | Trigger | Mitigation | Fallback |
| --- | --- | --- | --- | --- |
| `RISK-PERF-001` | The structure term dominates the 15 kb budget (≈75 % of it) and a 5× estimate error breaks `GATE-PERF-001` | `TST-PERF-001` breaches the budget | Windowed instead of global folding (`AMD-001`); incremental cache; two designated mitigations already specified (strided interior sampling, `E_max` rationing of structure) | Apply both mitigations; if still short, restrict `standard` tier at 15 kb to initiation-window structure and say so in the tier table |
| `RISK-ENG-001` | The two-language boundary raises contributor friction enough to hurt the open-source goal | Contributions concentrate in Python; Rust crates stagnate | The rust/python seam coincides with the infer/train seam, so model contributors never touch Rust; prebuilt wheels; dev container | Publish a "Rust-free contribution guide"; consider Numba fallbacks for the least hot kernels |
| `RISK-ENG-002` | Determinism under parallelism is lost by a well-meaning optimization | `GATE-DET-001` fails | Stream assignment by block/individual index; fixed-order reductions; the gate blocks releases | Revert the optimization. Determinism is not negotiable, because every reproducibility claim rests on it |
| `RISK-ENG-003` | A wrong exact delta silently corrupts tier-2 results | `TST-EPI-002` or its siblings fail | Delta-correctness property tests for every L1 term, run on every PR | Disable the L1 path for that term and fall back to rescoring — correct but ~640× slower |
| `RISK-OPS-001` | Artifact sprawl: hundreds of tenant adapters, unclear retention, growing storage | Registry size or resolution latency growth | Artifacts are small (Potts ≈164 KB); manifests are the queryable surface; channel history is bounded, artifact retention is not | Tier cold artifacts to cheaper storage; never delete anything reachable by a published lockfile |

---

## 4. Legal and licensing risks

### `RISK-DEP-001` — Structure-prediction dependencies are not redistributable

| | |
| --- | --- |
| **Risk** | ViennaRNA and LinearFold both permit research, educational and commercial *use and modification* but require that the package and derived works "are not redistributed for any fee, other than media costs". That is not an OSI-approved licence and it is incompatible with redistributing a combined work under Apache-2.0 or GPL. Parts of the ecosystem treat LinearFold as non-commercial out of caution. |
| **Trigger** | Any attempt to vendor either library, or to ship an extra that installs one. |
| **Impact** | **High** — it would make CodonForge's own distribution non-redistributable-for-fee, foreclosing commercial support, paid hosting and most distribution channels. |
| **Mitigation** | `AMD-009`: a permissively-licensed native windowed-accessibility kernel is the default engine; `StructureEngine` is a runtime-detected plugin seam; **no `[nonfree]` extra exists**; `P2-23` is an explicit licence-review ticket. |
| **Fallback** | Native engine only. `AMD-001` already makes it the preferred engine on performance grounds, so the fallback costs little. |

### `RISK-DEP-002` — LinearDesign is patent-encumbered

| | |
| --- | --- |
| **Risk** | Baidu Research has filed a patent on the LinearDesign algorithm. Patent encumbrance is distinct from licensing and is not cured by a permissive licence. |
| **Trigger** | Implementing the LinearDesign algorithm rather than merely benchmarking against it. |
| **Impact** | **High**, and asymmetric: cheap to avoid now, expensive to unwind later. |
| **Mitigation** | CodonForge does not implement the LinearDesign lattice-DP formulation. Its structure objective is windowed accessibility and its search is an MOEA — a different approach reached for independent reasons. Comparison happens through seed import of *output*. Apache-2.0's explicit patent grant is chosen partly for this context. |
| **Fallback** | Benchmark against published output only. Document the position in `NOTICE`. |
| **Action** | `P2-23` includes a patent-position review; any future proposal to implement a lattice-DP structure-and-codon co-optimizer must go through legal review first. |

---

## 5. Security and privacy risks

| Id | Risk | Trigger | Mitigation | Fallback |
| --- | --- | --- | --- | --- |
| `RISK-SEC-001` | Membership inference from a tenant adapter could reveal training sequences | Tenant model sharing is introduced | Tenant adapters are never served to another tenant and never published; isolation is structural, not statistical | If sharing is ever added, differential privacy becomes mandatory, not optional. The plan states this rather than deferring it silently |
| `RISK-SEC-002` | The community pool becomes a channel for controlled or hazardous sequences | Any donation | Review queue with biosecurity screening against controlled-sequence lists before a donation enters a public dataset (`P5-08`) | Suspend the pool. A design tool that aggregates user sequences must be able to stop doing so |
| `RISK-SEC-003` | A malicious `ConstraintSet` regex causes catastrophic backtracking | User-supplied `custom_regex` | Compile to a DFA (no backtracking) with a size cap; reject patterns exceeding it | Disable `custom_regex` for untrusted callers |

---

## 6. Product risks

| Id | Risk | Trigger | Mitigation | Fallback |
| --- | --- | --- | --- | --- |
| `RISK-UX-001` | The many-objective front is unusable: hundreds of equally-good candidates | `M ≥ 5` on any real job | `AMD-008`: role demotion to ε-constraints and diagnostics, ≤4 default axes, mandatory reduction, uncertainty-aware selection | Default to 3 axes with the rest as ε-constraints, and surface a preference-articulation step |
| `RISK-UX-002` | Users treat a tier-3 burden number or a fallback-host prediction as authoritative | Any tier-3 or fallback job | Different unit, `calibrated: false`, `low_confidence_host` tag, structured warnings, `TST-BURD-005` | Refuse to offer the term rather than qualify it, if the qualification proves ineffective |
| `RISK-UX-003` | The tool returns nothing for over-constrained inputs, with no explanation | Any conflicting constraint set | `AMD-010`: the five-check feasibility analyser, run at submission and available as its own endpoint, always naming a constraint to relax | — |
