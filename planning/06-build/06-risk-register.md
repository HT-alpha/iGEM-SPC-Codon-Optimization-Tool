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
| **Impact** | **High** — it would make DNA Design Suite's own distribution non-redistributable-for-fee, foreclosing commercial support, paid hosting and most distribution channels. |
| **Mitigation** | `AMD-009`: a permissively-licensed native windowed-accessibility kernel is the default engine; `StructureEngine` is a runtime-detected plugin seam; **no `[nonfree]` extra exists**; `P2-23` is an explicit licence-review ticket. |
| **Fallback** | Native engine only. `AMD-001` already makes it the preferred engine on performance grounds, so the fallback costs little. |

### `RISK-DEP-002` — LinearDesign is patent-encumbered

| | |
| --- | --- |
| **Risk** | Baidu Research has filed a patent on the LinearDesign algorithm. Patent encumbrance is distinct from licensing and is not cured by a permissive licence. |
| **Trigger** | Implementing the LinearDesign algorithm rather than merely benchmarking against it. |
| **Impact** | **High**, and asymmetric: cheap to avoid now, expensive to unwind later. |
| **Mitigation** | DNA Design Suite does not implement the LinearDesign lattice-DP formulation. Its structure objective is windowed accessibility and its search is an MOEA — a different approach reached for independent reasons. Comparison happens through seed import of *output*. Apache-2.0's explicit patent grant is chosen partly for this context. |
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

---

## 7. Latent hazards — correct today, wrong under a change that has not happened yet

These are not open bugs. Each is a working implementation whose correctness rests
on a condition the code does not enforce, so the failure arrives with a *future*
change rather than with the next test run. They are recorded here rather than
only in a source comment because the engineer who trips one will be working on
something else entirely and will not know to look at the class that carries the
caveat.

| Id | What is true today | What would break it | Where the caveat lives now |
| --- | --- | --- | --- |
| `RISK-LAT-001` | ViennaRNA's temperature is a **process global** (`RNA.cvar.temperature`). The binding sets it on every call, and the structure term runs outside the parallel block sweep, so nothing races. | **Parallelizing structure evaluation inside one process.** Two threads folding at different temperatures would interleave their writes to that global and each would silently fold at the other's temperature — producing plausible energies, no error, and a `GATE-DET-001` violation that only appears as irreproducibility. A worker *process* per thread is safe; threads are not. | `ViennaRNAEngine` docstring |
| `RISK-LAT-002` | `accessibility_penalty` floors accessibility at `1e-300` before taking a log, so a fully-paired window has a large finite cost rather than infinity. | **Treating the floor as a physical constant.** It is a numerical choice with no derivation behind it: it sets how much worse "completely inaccessible" is than "nearly inaccessible", which is a real modelling decision about how hard the objective pushes on the worst window. Anyone tuning structure weight is implicitly tuning against this number and should know it was picked, not measured. | `accessibility_penalty` docstring |
| `RISK-LAT-003` | Tier 1's `half_maximal_expression` locates `θ` by **linear interpolation** between measured panel points, not by fitting the saturating curve `ADR-0004` §2.1 describes. Defensible for the ≥20-point, ≥1-decade panel the tier requires. | **Any host actually reaching tier 1.** No host does today, so the approximation has never been exercised against real data. The interpolation is sensitive to point spacing near the half-maximum in a way a fitted curve is not, and tier 1 is the tier whose whole value is that it rests on measurement. | `capacity.rs`, `fit_capacity_tier1` |
| `RISK-LAT-004` | `expected_gain` approximates expected hypervolume gain by distance past the current best in an individual's strongest objective. Exact EHVI is exponential in objective count and would cost more than the L0 calls it saves. | **A dense front, or many objectives.** The approximation is good when the front is sparse and degrades as it fills in — which is exactly the late-run regime where rationing matters most. It is component (c) of four, so its errors are bounded by the other three components, but it is not what the name implies. | `rationing.rs`, `expected_gain` |
| `RISK-LAT-005` | `checkout@v5`, `setup-python@v6`, `upload-artifact@v6` and `download-artifact@v7` declare Node 24. `rust-cache@v2` and `maturin-action@v1` do too. `dtolnay/rust-toolchain` is composite. | **Pinning those actions back to a Node 20 major.** The runner warns, and GitHub can turn the warning into a failed step. | `.github/workflows/ci.yml`, `.github/workflows/wheels.yml` |

A hazard leaves this table when the condition becomes enforced rather than
assumed — a thread-safety guard, a fitted curve, a derived constant — not when
someone decides it is unlikely.

---

## 8. Timeline risks for the iGEM build

These are the claims that would make the submission false if they were stated
as finished. Each one is a current fact about this build, not a forecast.

| Id | Risk | Trigger | Mitigation | Fallback |
| --- | --- | --- | --- | --- |
| `RISK-TIME-001` | `GATE-BEN-001` is reported as an MFE comparison that was not run | `gates.json` has `definition` other than `cai_only` while `claimed` is true, or a write-up says LinearFold and windowed accessibility were compared | The MFE fronts are still not the same objective. The record that exists is the CAI-only redefinition (`DEP-BEN-003`): claimed on 50 of 50 proteins after `DEP-SRCH-020`. The v2.4.2 arithmetic delta left `aceB` and `thrA` at CAI 0. `check_cai_gate` refuses a claim that covers a recorded score below `cai_greedy`. `DEP-BEN-002` | The submission names the 50-protein CAI-only result and does not say the MFE gate passed |
| `RISK-TIME-002` | A full release wheel ships before legal sign-off | A stable tag (`v1.2.3`, no hyphen), or a wheel job that skips `legal-clearance` | `tools/check_release_clearance.py` blocks that job. A pre-release tag such as `v0.1.0-igem` skips the job because the tag contains a hyphen; it does not build wheels | The gate stays red until a faculty advisor, PI, or tech-transfer contact signs `docs/legal-clearance.yaml` |
| `RISK-TIME-003` | `P3-14` (LM scoring) is pulled onto the iGEM critical path | `lm.naturalness` becomes a default objective, or weights are vendored to make a demo | Cut. Recorded in [`07-igem-batch2-and-phase3.md`](07-igem-batch2-and-phase3.md). An unregistered term is `CF-LAND-001` | The pairwise Potts surface is the Phase 3 minimum. The LM term waits until that surface is scored and the legal gate is signed |
| `RISK-TIME-004` | `/v1/feasibility` looks implemented | A `200` whose body is a stub, or a tri-state boolean | `POST /v1/feasibility` runs the five existing checks (`DEP-API-002`). An empty body is `400`. The page still reads `constraint_report.feasible` on the front | The page reads `constraint_report.feasible` on the front, which is a boolean, plus `unknown_terms` |
| `RISK-TIME-005` | The `EXP-05` repair default flips on a campaign the rule does not accept | A change to the shipped repair rule, or a write-up that treats a zero hypervolume margin as a win, or a write-up that counts `rplL` or `sucD` as a completed arm | `exp05.decide` was not edited. The 20-generation run on `rpsJ` and `rpsL` had margin 0 (`rplL` is `CF-SRCH-006`). A further run, 20 seeds per arm and 40 generations with stagnation window 41, used `gapA`, `pfkA`, `mdh`, `eno` and `pykF` (332–471 codons; tier 1 accepts the interior motif). `sucD` was left out: the same motif is `CF-SRCH-006` there. Every one of the five reached 40 generations. Hypervolume margin was 0 and the Lamarckian standard deviation was 0 on each. The rule keeps Lamarckian. More seeds cannot separate a margin that is exactly 0 | Lamarckian stays the default |
| `RISK-TIME-006` | A protein or an expression number is invented to fill the panel | A CDS with no accession, a reverse-translated GFP, or a synthesis rate typed out of a figure | The accession guard, the `value: null` refusal, and the GFP exclusion in the panel | The protein stays off the list. The citation (`GSE53767`) is the evidence, not a transcribed rate |
| `RISK-TIME-007` | Primer ΔG is treated as an OligoAnalyzer clone | A hairpin number that matches mFold without saying Table 2 has no loop penalty | `docs/primer-thermo-validation.md` records the Table 2 hand sum and the gap versus the published OligoAnalyzer cells. BLAST is stated as not checked. Primer3 is not vendored | The suite does not claim IDT parity for hairpins. A difference above 0.5 kcal/mol stays written down |
