# Phase 1 — Hard-Constraint Filter + Baseline Optimizer

Satisfies `BRIEF-2.6-a` and acceptance criterion `BRIEF-5.10`:

> The build sequencing is granular enough that Phase 1 could be started immediately by an implementation team with no further design discussion.

Every ticket below names its module, its deliverable, its acceptance tests and its dependencies. Nothing in Phase 1 requires a design decision that is not already made in this package.

**What Phase 1 delivers:** a working, installable tool that takes a protein or CDS plus a constraint set and returns a feasible, CAI-optimized, fully-reproducible design — with the feasibility analyser, the seed importers, the repair operator, the three-tier search skeleton and the determinism contract already in place. It is useful on its own, and it is the correctness scaffold for everything after.

**Phase 1 has no learned models.** Every term is a built-in computation. That is deliberate: it means the search engine, the constraint machinery and the reproducibility machinery are all validated before any model can be blamed for a bug.

---

## Infrastructure

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P1-01` | Workspace scaffold | — | Cargo workspace with the 18 crate/package directories of the repo layout; `pyproject.toml`; `maturin` build of `codonforge._core`; dev container with both toolchains | `maturin develop` then `import codonforge` succeeds on Linux/macOS/Windows | — | M |
| `P1-02` | CI pipelines | — | Per-module test jobs; `clippy` + `ruff` + `mypy`; wheel build matrix (manylinux, macOS arm64/x86, Windows); nightly benchmark job stub | A PR touching one crate runs that crate's isolated slice | `P1-01` | M |
| `P1-03` | Schemas into the tree | — | Copy `planning/02-schemas/*.schema.json` to `schemas/v1/`; wire `validate_planning_package.py` into CI; add a drift check that fails if a Rust struct or `.pyi` stub diverges from its schema | Schema-validation job green; a deliberate drift fails it | `P1-01` | M |
| `P1-04` | Error taxonomy | `cf-types` | The `CF-*` code registry, `ProblemDetails` construction, and a test that every emitted code is documented | `TST-API-004` (code registry portion) | `P1-01` | S |

---

## Core types

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P1-05` | Codon indexing and genetic codes | `cf-types` | `index = 16b₀+4b₁+b₂` codec; the 27 NCBI translation tables defined in the id range 1..33 vendored as JSON (7, 8, 17-20 are retired); forward and reverse (synonymous-family) tables | `TST-CORE-002` | `P1-01` | S |
| `P1-06` | `CodonSequence` | `cf-types` | The type, its five invariants enforced at construction, `substitute()` with a synonymy precondition, flank handling, locked spans | `TST-CORE-001` | `P1-05` | M |
| `P1-07` | Canonical serialization and digests | `cf-types` | RFC 8785 canonical JSON; sha256 digests; stable across key order, float formatting and platform | `TST-CORE-003` | `P1-06` | S |
| `P1-08` | Identifiers and units | `cf-types` | `HostId`, `TenantId`, `JobId` (ULID), `ArtifactId`, `TermId` newtypes with the schema's regex validation; unit vocabulary; objective vector type | `contract_*` for `common.schema.json` | `P1-05` | S |

---

## Ingest

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P1-09` | FASTA / GenBank readers | `cf-ingest` | Parse to `CodonSequence`; frame detection; U→T transcoding; explicit handling or rejection of ambiguity codes and lower case | `TST-INGE-001` | `P1-06` | M |
| `P1-10` | Construct context assembly | `cf-ingest` | Build `ConstructContext` from GenBank annotations; infer cistrons from CDS + upstream RBS; nt→codon span conversion (the *only* place a frame error is possible) | `contract_*` against `construct_context.schema.json` | `P1-09` | M |
| `P1-11` | Host profile loader | `cf-ingest` | Codon usage table ingest (genome-wide and highly-expressed subsets); tRNA gene copy numbers; tAI weight derivation; `data_completeness` population | `TST-INGE-002` (host portion) | `P1-09` | M |
| `P1-12` | Seed importers | `cf-ingest` | Adapters for LinearDesign, CDSfold, CodonTransformer, DNA Chisel, plain FASTA, native gene, previous job; acceptance verdicts; **positional protein diff on mismatch** | `TST-INGE-003` | `P1-09` | M |

---

## Constraints and repair

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P1-13` | Constraint set loading | `cf-constraints` | Parse `ConstraintSet`; expand `assembly.standard` and `vendor_profile` into concrete rules; report every auto-added constraint to the user | `contract_*` against `constraint_set.schema.json` | `P1-08` | S |
| `P1-14` | Motif automaton | `cf-constraints` | Aho–Corasick over IUPAC-expanded patterns, both strands; incremental rescan of an `m_max+3` nt window after a substitution | `TST-CONS-001` (motif classes) | `P1-13` | M |
| `P1-15` | k-mer index | `cf-constraints` | Hash map with **overlapping** counts, matching the prototype's `count_overlapping` semantics; `O(k)` incremental update; optional canonical-strand mode | `TST-CONS-001` (k-mer) | `P1-13` | M |
| `P1-16` | Homopolymer and GC detectors | `cf-constraints` | Run-length index; GC prefix sums with `O(1)` per-window queries | `TST-CONS-001` (remaining classes) | `P1-13` | S |
| `P1-19` | Feasibility analyser | `cf-constraints` | All five checks, each reporting a constraint to relax and by how much | `TST-CONS-006` | `P1-14`–`P1-16` | M |
| `P1-17` | Repair operator | `cf-constraints` | Depth-1 exhaustive search; depth-`d` beam; the `(edits, objective cost, random)` ordering key; escalation; bounded give-up distinguishing "provably no fix" from "budget exhausted" | `TST-CONS-002`, `TST-CONS-003`, `TST-CONS-004` | `P1-14`–`P1-16`, **`P1-19`** | L |
| `P1-18a` | Tabu ledger | `cf-constraints` | Per-run `(position, codon)` ledger, `record`, generation-clock `tick`, and the 10 % leak decision | unit tests | `P1-14`–`P1-16` | S |
| `P1-18b` | Tabu wiring | `cf-search` | `TABU_TICK` at the generation boundary; the leak consulted by the mutation operator; the mutation/repair cycle driver | `TST-CONS-005` | `P1-18a`, `P1-23` | S |

> **`P1-18b` hands two obligations to `P1-29`.** Written here rather than left for whoever builds `P1-29` to rediscover, following the precedent `P1-23` set with the operator warnings.
>
> 1. **Job lifecycle.** `run_generations` takes a generation count, a root seed and a repair budget as inputs. The four things a real job adds — a compute budget with a tier plan instead of a generation count, an RNG plan seeded from the lockfile's `rng.root_seed` rather than from anywhere else, checkpoint/resume, and cooperative cancellation — are all `P1-29` deliverables. The loop is a driver; it does not invent any of them.
> 2. **The operator warnings, still parked.** `P1-23` required that `OperatorMix::warnings()` reach the user, since the mix renormalizes from 0.35/0.35/0.20/0.10 to 0.35/0.10 when no context model is registered. The generation loop is the first execution context with a lifetime, so it now *collects* them into `GenerationReport::warnings` — the plumbing moves earlier. It still cannot discharge the obligation: the destination is `job_status.warnings`, and that document is `P1-29`'s. Confirmed parked, not dropped.

> **Sequencing correction (implementation, 2026-08).** Two corrections here, found by the pre-work audit.
>
> **`P1-19` now precedes `P1-17`.** `P1-17`'s three declared tests do not need the analyser, so the original ordering was defensible on the letter of the ticket. But repair's termination condition is only *principled* in combination with it: a bounded budget makes each call terminate while still allowing repair to spend that budget on every candidate forever when a constraint is unsatisfiable in principle. See [`../03-algorithms/06-constraint-repair.md`](../03-algorithms/06-constraint-repair.md) §5.1.
>
> **`P1-18` split.** Repair's signature takes the tabu ledger, so the ledger type must land with `P1-17` — but `TST-CONS-005` scripts a mutation/repair cycle over generations, which needs `P1-23`'s mutation operator and a generation clock. `P1-18a` (the ledger, unit-testable in isolation) lands now; `P1-18b` (wiring and `TST-CONS-005`) is scheduled with `P1-23`. Same shape as the `P1-27a`/`P1-27b` split.

---

## Landscape skeleton

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P1-20` | Landscape contract types | `cf-landscape-api` | `LandscapeQuery`/`Response` in Rust and Python; the three target shapes; the L0/L1/L2 capability model; `TermRole`; `TermRegistry` | `TST-LAPI-001`, `TST-LAPI-002`, `TST-LAPI-003` | `P1-08` | M |
| `P1-21` | Built-in terms | `cf-landscape` | `usage.cai`, `usage.gc_content`, `constraint.violation_count` — all L1 with exact deltas | `TST-EPI-002`-style delta test applied to built-ins | `P1-20` | M |
| `P1-22` | Composition and cache | `cf-landscape` | Term registry resolution; window-keyed evaluation cache namespaced by lockfile digest; `partial`/`UNKNOWN` propagation; per-term capability degradation; role-gated objective vector | `TST-LAND-001`, `TST-LAND-002`, `TST-LAND-003`, `TST-LAND-004`, `TST-LAND-005` | `P1-21` | M |

### `P1-20` — the three invariants that are not shape properties

`TST-LAPI-001` and `TST-LAPI-002` both check the *shape* of the contract: that documents round-trip byte-identically, and that a capability claim matches a method surface. Composition in `P1-22` additionally relies on three *behavioural* properties of the same contract, and the sequencing audit found none of them exercised:

1. **Objective-vector ordering.** [`landscape_response.schema.json`](../02-schemas/landscape_response.schema.json) states that `order` is explicit and must never be inferred from `terms` key iteration. That is an obligation on whoever builds the vector, which is `P1-22`. A composition that built the vector from registry iteration order would satisfy both declared tests while making two runs of the same job produce differently-ordered vectors — and a Pareto front whose axes are permuted between runs is not comparable to itself.
2. **Seed determinism.** `TST-LAND-001` asserts that a cache hit returns exactly what a cold evaluation returns. That claim is vacuous unless a cold evaluation is reproducible, which rests on the query's `seed` field.
3. **Term independence.** Composition evaluates terms in whatever order the registry yields. If a provider could observe or perturb another provider's evaluation, composition would be order-sensitive and the cache would be unsound.

These are properties of the contract, so they belong to `P1-20` rather than to its consumers. No resequencing is needed — `P1-20` already precedes `P1-21`/`P1-22` — but the acceptance criteria gain `TST-LAPI-003`. This is the `P1-17`/`P1-19` pattern again: declared tests that pass without the guarantee the ticket actually owes its consumers.

### `P1-22` — what "black-box oracle" costs the composition layer

Two audit findings, one an architecture defect and one an acceptance gap.

**The `cf-modelreg` edge could never have compiled.** [`modules.yaml`](../01-modules/modules.yaml) declared `cf-landscape` (runtime `rust-core`) as depending on `cf-modelreg` (runtime `python-orchestration`). A Cargo crate cannot depend on a Python package, so this was a declaration no implementation could satisfy. The layering lint did not catch it because its `edges` check looks only for edges present in the tree but absent from `modules.yaml`, never the reverse; a declared-but-impossible edge is invisible to it. What `cf-landscape` genuinely needs from the registry is the resolved lockfile digest that namespaces the evaluation cache, and a digest is a value. The edge is removed and the digest is supplied by the orchestrator at construction. See `DEP-LAND-001`, and check 7 (`runtime-compat`) in `tools/check_module_layering.py`, which now makes this class of edge a lint failure. A side effect is that the apparent `P1-22` → `P1-28` ordering constraint dissolves: the two tickets are genuinely independent.

**The delta path rescored the base.** Found later, by the `P1-23` audit, and recorded here because it belongs with the rest of `P1-22`'s story. `Composer::compute` computed `provider.score(ctx_base) + delta`, so every delta query paid an `O(L)` pass over the base — the same order as the rescore the delta exists to replace, and for CAI strictly worse, since a full rescore pays that pass once and the delta path paid it twice. Nothing observed it: the numbers were right, `full_rescores` read zero, and only the cost differed. `DEP-LAND-003` adds a `rebase` step to the seam and caches the pinned base per `(base, term)`; `TST-LAND-006` asserts on provider call counts, because a value assertion cannot distinguish the two. Two smaller instances went with it — `CaiTerm::delta_table` recomputed its `O(L)` log-sum per alternative, making a table `O(L·P·S)` against a promised `O(P·S)`, and `GcContentTerm::delta` summed GC across the sequence to obtain a denominator equal to `codon_count * 3`.

**Mixed-capability composition is untested by the declared criteria.** `TST-LAND-001` (cache), `TST-LAND-002` (`partial`/`UNKNOWN`) and `TST-LAND-003` (dependency tree) all pass under a composition that responds to a single L0-only term by rescoring *every* term in full. That is the failure mode the binding decision on exact deltas exists to prevent — the measured advantage is ~640x at 2 kb and ~4800x at 15 kb, and it is forfeited entirely if one L0 term can drag the whole query onto the slow path. `TST-LAND-004` pins degradation to the individual term. `TST-LAND-005` covers the role gate that keeps a diagnostic term out of the objective vector.

---

## Search

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P1-23` | Genotype and operators | `cf-search` | Genotype with an incrementally maintained hash; the four mutation operators, two of them against the `ContextModel` seam; block-uniform crossover; duplicate suppression; seed injection | `TST-SRCH-006`, `TST-SRCH-010` | `P1-20` | M |

### `P1-23` — what ships, what waits, and why the operator mix is smaller than §5.1

**Two of the four mutation operators have no model to sample from.** `potts_guided` and `block_gibbs` both draw from the Potts conditional ([`03-nsga-loop-and-operators.md`](../03-algorithms/03-nsga-loop-and-operators.md) §5.1), which is `cf-epistasis` and arrives in `P3-05`. `cf-search` is also specified as talking only to `cf-landscape-api`, so even in Phase 3 the conditional must arrive through a seam rather than a direct call. Both operators are implemented against `ContextModel` and are registered only when a model is supplied; without one, `OperatorMix::without_context_model` renormalizes the remaining weights and reports the absent operators as `CF-SRCH-001` warnings. See `DEP-SRCH-005`.

Backing them with host codon usage instead would have been the tempting move and the wrong one: `potts_guided` over a context-free distribution *is* `uniform_synonymous` with different weights, so the mix would have looked complete while exploring differently, and the change in behaviour would have been attributed to anything but the operator that caused it.

The renormalization is a real behavioural change — 0.35/0.35/0.20/0.10 over four operators becomes 0.35/0.10 over two — so it is surfaced rather than absorbed. `OperatorMix::warnings()` returns the `CF-SRCH-001` problems; **`P1-29` must put them into `job_status.warnings`**, because nothing in `P1-23` has a job to warn from.

**`P1-18b` is not needed, and the boundary is exact.** `P1-18a` shipped `TabuLedger::is_tabu` and `should_suppress` as pure queries, so the mutation operators consult the ledger now and `TST-CONS-005`'s consultation half is tested against a hand-advanced ledger. What waits for `P1-18b` is the clock: `TABU_TICK` belongs to the generation loop, which does not exist until `P1-23` is in place, and `TST-CONS-005` scripts a mutation/repair cycle *across* generations.

**A defect worth carrying forward to `P1-33`.** `SearchRng::new` initially used `seed | 1` to avoid xorshift's absorbing zero state, which maps every even seed onto its odd successor: seeds `2026` and `2027` produced byte-identical runs. Since the root seed comes from the lockfile, two jobs a user seeded deliberately differently would have agreed and looked like a reproducibility success. Replaced with the splitmix64 finalizer, which is a bijection. `cf-constraints`' `RepairRng` had the same construction and the same fix. `GATE-DET-001` and `P1-33` rest on this layer being right.
| `P1-24` | Tier 1 constrained beam DP | `cf-search` | Beam DP on the product with the motif automaton; flank-primed junctions; diversity-capped beam selection | `TST-SRCH-002` (tier-1 portion), `TST-SRCH-011` | `P1-14`, `P1-23` | L |

### `P1-24` — the cost claim, and two boundaries the pseudocode leaves open

**The inner loop is where a values-only test goes blind.** §2 evaluates `Δ ← Σ_m λ_m · delta_m(x | ctx)` once per `(beam state, synonym)` pair — `O(L·b·S)`, about 120 000 evaluations on a 666-codon protein at `b = 64`. An implementation that issued a landscape query for each would return the same seeds and miss the budget by three orders of magnitude, and nothing about the answer would look wrong. `TST-SRCH-011` therefore asserts on landscape call counts: one `rebase` and one `delta_table` per scored term, and no full-sequence `score` at all. Same lesson as `DEP-LAND-003`, applied before the mistake rather than after.

It also carries the tier-planning refusal that the `P1-25` note describes, one tier early: a term whose provider is L0 has no delta table, so tier 1 refuses it up front with `CF-LAND-002` instead of discovering it mid-run or quietly optimizing a smaller objective than the caller asked for.

**Junctions.** §2's pseudocode starts the automaton at `A.start` with an empty prefix. That is right for a bare coding sequence and wrong for one inside a construct: a motif straddling either flank junction is then invisible while the DP extends. Tier 1 primes on the 5′ flank and continues through the 3′ flank, making both junctions structural (`DEP-SRCH-006`). Junction semantics are exactly what `EXIT-001`/`EXIT-002` concern, so leaving them to a post-check would have been the wrong default.

**And the boundary that limits the pruning.** Neither flank walk prunes on a match lying wholly inside immutable sequence. Such a violation is a property of the input rather than of any candidate — `DEP-CONS-002`, the same reason `ConstraintChecker::is_feasible` ignores backbone-confined violations — so pruning on one would make tier 1 reject every candidate for something no candidate can change. The fixture's own backbone contains a BamHI site, which is how this surfaced. A motif rule with `max_occurrences > 0` is likewise not pruned, since distinguishing "seen fewer than the allowance" from "seen too many" needs a counter in the DP state, which is the state-space multiplication §2 declines to pay for k-mers.
| `P1-25` | Tier 2 red-black block CA | `cf-search` | The partition and its `B ≥ R` assertion; same-colour block updates computed before any is committed; incumbent-seeded beam for monotonicity; 1-opt polish; accept-on-truth; the tier-2 capability precheck | `TST-SRCH-001`, `TST-SRCH-002`, `TST-SRCH-012` | `P1-24` | XL |

### `P1-25` — four things the specification leaves to the implementation

**The partition condition is `B ≥ R`, not `B ≥ W`.** §3.1 states it as `B ≥ W`, which is the condition only when `W` is the largest interaction radius among the *scored terms* — and §4 exists precisely because a term can have a larger one. Nothing declared a radius, so the assertion could only have compared one configured number against another and would have passed while the premise was false. `TermProvider::interaction_radius` is now required, and `Tier2Plan` checks the block size against what providers declare. A term declaring `Unbounded`, or a radius above `B`, is refused with `CF-SRCH-008` naming `P1-26` — the fallback that answers that case is not built, and assuming independence silently would be worse than refusing.

**The partition is over codon positions, not editable positions.** §3.1's prose says "partition the editable positions" while its diagram shows codon ranges, and they differ on a sequence with locks. Codon positions are used, because the independence argument is about codon distance: two positions in different same-colour blocks are then at least `B` codons apart directly, with no reasoning about how many locked positions sit between them. A block containing only locked codons is a no-op.

**`ASF` is substituted, and the substitution is recorded rather than implied.** The augmented Tchebycheff scalarizing function of §3.2 needs the ideal and nadir points, which come from `cf-metrics` and the tier-3 population — neither exists yet. Tier 2 accepts on a weighted sum instead. What matters for correctness is the *two-step* structure, which is identical either way: the DP maximizes an additive linearization, and the accept decision is taken on a freshly evaluated true value, so monotonicity holds in the true objective rather than in the linearization. Substituting the true `ASF` is a `P1-25` follow-up once `cf-metrics` supplies reference directions.

**A block can create a motif across its *right* edge, and nothing was watching.** The block DP prunes matches ending inside its own range; one that starts inside and completes in the next block is invisible to it, and nothing downstream catches it either, because the accept step compares objectives and an infeasible candidate can score well. Found by a negative control that would not fire: the guarantee it was meant to break was not being kept in the first place, and the run it exposed was committing an infeasible incumbent and counting it as an improvement. `update_block` now walks the fixed context to the right of the block before choosing, bounded by the longest compiled pattern, and the accept step refuses an infeasible candidate outright. Structurally the same gap as tier 1's 3′ junction, one boundary further in.

**The incumbent is kept as a `b + 1`th beam state, not in competition with the best `b`.** §3.3 writes `next ← next ∪ {incumbent_path_prefix(i)}` before `TOP_B`, which reads as the incumbent being eligible for truncation — and then "block updates are made non-worsening by seeding the beam with the incumbent" is a tendency rather than a guarantee. Displacing the best state instead is worse: at `b = 1` the beam would always *be* the incumbent and the block could never improve.

**Occurrence counting is affordable at tier 2, and is deliberately not taken.** Tier 1 declined to prune motif rules with `max_occurrences > 0` because a counter multiplies the DP state space. Tier 2's situation genuinely differs — occurrences *outside* a block are fixed and countable once, so the block DP would only need a counter over the remaining allowance, which is cheap. It is scoped rather than implemented: taking it would make the DP and the post-check disagree about which constraints are structural, and an asymmetry between tiers deserves its own ticket rather than arriving as a side effect of this one.

> **Tier-2 capability precheck.** `P1-21` found that `constraint.violation_count` is L1 only when every compiled rule has bounded reach; a k-mer repeat rule counts over the whole construct, so the term is honestly L0 and every tier-2 move rescores it. Composition already degrades correctly and reports `full_rescores`, but by then the tier is running. [`04-model-resolution.md`](../04-api/04-model-resolution.md) §2.1 names the remedy — the planner drops an L0 term from tier 2 with a warning, or fails — and `cf_landscape::planned_capability` exists so the planner can ask before starting rather than discovering it mid-run. Wiring that warning is `P1-25`'s, not `P1-21`'s: there is no tier planner yet to warn from.
| `P1-26` | Seam reconciliation fallback | `cf-search` | Overlapping-block path with the exact 3-state chain DP and seam polish | `TST-SRCH-003` | `P1-25`, **a term with a radius above `B`** | M |

> **`P1-26` held, deliberately.** Its dependency on `P1-25` is satisfied, and it was still not taken in the same round. The trigger cannot occur in Phase 1: `usage.cai` and `usage.gc_content` declare `Positionwise`, and `constraint.violation_count`'s radius is bounded by the motif length or the GC window, all far below any sensible `B`. More decisively, `EXACT_OBJECTIVE_CONTRIBUTION(seam, neighbours)` needs a term that scores a seam region under long-range coupling, and §4 names the real instance — RNA base pairing, which is `P2-11`. Building that interface now means guessing at a term class that does not exist, which is the shape that produced placeholders in earlier rounds, and `TST-SRCH-003` would be satisfied entirely against a synthetic long-range term: it would prove the chain DP's arithmetic and nothing about reconciliation. So the dependency is restated as `P1-25` **plus a term that needs it**, and `P1-25` refuses the condition with `CF-SRCH-008` rather than leaving the guarantee implicit and unverified.
| `P1-27` | Baseline backtracking solver | `cf-search` | Port of the existing `codon_6.0.py` stack-based backtracker as `baseline_backtracking`, kept as a reference implementation and a differential-test oracle for the constraint layer | Differential test: for inputs the prototype solves, the new constraint checker agrees on feasibility for 10⁴ sequences | `P1-14`–`P1-16` | M |

`P1-27` is worth doing rather than discarding the prototype: it is an independently-written implementation of the same constraint semantics, which makes it a genuine oracle for `P1-14`–`P1-16`. It also gives a documented answer to "does the new tool still do what the old script did".

> **Sequencing correction (implementation, 2026-08).** `P1-27` originally declared `P1-17` (the repair operator) as its dependency. That was over-stated: the acceptance criterion exercises only the constraint *checker*, and the legacy `count_overlapping` semantics it exists to pin down are implemented in `P1-15`. The dependency is therefore `P1-14`–`P1-16`, and `P1-27` is scheduled **immediately after `P1-16`, before `P1-17`** — a differential oracle that lands after the code it is meant to check has already been built on top of cannot catch drift at the point the drift would be introduced.
>
> **Split (implementation, 2026-08).** `P1-27` bundled two deliverables with one acceptance criterion: a differential oracle for the constraint checker, and a port of the prototype's backtracking *solver* as `baseline_backtracking`. Only the first is covered by the stated criterion, and the second belongs to `cf-search`, which does not exist until `P1-23`. Split accordingly:
>
> | Id | Deliverable | Module | Deps | Status |
> | --- | --- | --- | --- | --- |
> | `P1-27a` | Legacy feasibility oracle and the 10⁴-sequence differential (`TST-CONS-007`) | `cf-constraints` | `P1-14`–`P1-16` | done |
> | `P1-27b` | `baseline_backtracking` solver port, kept as a reference implementation | `cf-search` | `P1-23`, `P1-24` | done |
>
> **`P1-27b` — what the port is an oracle for, and what it is not.** The ticket says "port", which puts a Rust transcription in `cf-search`; the standard `P1-27a` set says an oracle must be the actual prototype, not a transcription that could share a bug with what it checks. Both are right about different objects. The port is a *reference solver* — depth-first with restarts, against tier 1's breadth-first beam — and it earns that role only by being demonstrably faithful, so `python/tests/test_p1_27b_baseline_fidelity.py` compares it against the loaded `codon_6.0.py` before anything else compares against it. The comparison is on the incremental prefix predicate, prefix by prefix across two `(k, limit)` regimes, because that is where a paraphrase of the two windowed checks and their `+ 2` slack would actually differ; sequence equality is not compared, since the prototype's ordering comes from Python's Mersenne Twister and reproducing that stream would prove something about two RNGs rather than two searches.
>
> **The two solvers are not comparable on score.** `optimize_cds_sequence` has no objective — it returns the first assignment that survives its checks. The comparable surface is feasibility, over the constraints both enforce, which excludes `kmer_repeat`: the prototype counts it and §2 excludes it from the DP because it is not a regular language. A differential spanning that row would report design as drift, so the asymmetry is asserted rather than assumed. Completeness is checked against exhaustive enumeration instead of against the baseline, because the baseline's give-up is a timeout and not a proof — the same distinction `P1-17`'s `GiveUpReason` draws.
>
> Findings from `P1-27a` are recorded in `docs/legacy-parity.md`; two of them ask for a decision before Phase 1 exit.

---

## Registry, orchestration, surfaces

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P1-28` | Minimal registry | `cf-modelreg` | Manifest and lockfile types; local filesystem artifact store; the `pin`, `base` and `builtin` rungs (channels arrive in Phase 3) | `TST-MREG-002` | `P1-08` | M |
| `P1-29` | Job lifecycle | `cf-orchestrator` | In-process job execution in segments; lockfile sealed **before** any search work; compute budget and tier plan; RNG plan from `rng.root_seed`; checkpoint/resume; cooperative cancel; `job_status.warnings` aggregation | `TST-ORCH-001`, `TST-ORCH-002` | `P1-25`, `P1-28` | L |

### `P1-29` — the three obligations, and three things the audit changed

**All three carried-forward obligations are discharged.** `P1-23` required `OperatorMix::warnings()` to reach a user; `P1-18b` moved the plumbing as far as `GenerationReport::warnings` and said plainly it could go no further; `P1-25` added `Tier2Plan::warnings()` with the same problem. All three now land in `job_status.warnings` through one path, each tagged with `stage` and `subject` (`DEP-ORCH-001`) so a caller can attribute a warning without reading its prose. A fourth source joined them, described below.

**One decision moved up a layer.** `04-model-resolution.md` §2.1 says "the planner" either drops an unaffordable term with a warning or fails. Until `P1-29` there was no planner, so each tier decided for itself and disagreed: tier 1 refused an L0 term outright with `CF-LAND-002`, tier 2 dropped it with `CF-SRCH-010`. The same term in the same job got two answers depending on which tier looked first. `TierPlan` now decides once and hands each tier a list it can score; the tiers' own refusals remain as backstops that should never fire.

**`CF-SRCH-008` is not a warning.** It is sometimes listed with the warning codes, and it does not belong there: a term wider than the block does not degrade the search, it makes tier 2's independence claim false. It refuses the job and appears in `job_status.problem`.

**The run is a loop over segments, and that is what makes cancellation safe.** ADR-0001 caps boundary crossings at one per generation, so a Python job loop checking a cancellation flag mid-generation would be over the limit. Rust advances whole segments and returns a state; the job layer decides between segments whether to continue, checkpoint or stop. Cancellation and checkpointing therefore happen at the same points, so there is no such thing as a checkpoint taken mid-cancellation — the persisted state is always a completed generation. The state returned *is* the checkpoint, so there is no second representation to keep in sync.

**Resumed runs are identical, not merely reproducible.** Generation streams derive from the *absolute* generation index, so generation `k` sees the same randomness whether it is the first of a resumed segment or the fiftieth of an uninterrupted run. A counter restarting at zero on resume would have made a resumed job reproducible against itself and quietly different from an uninterrupted one — both looking correct. The tabu ledger travels in the checkpoint for the same reason: it carries state across generations by design, and dropping it would re-propose exactly the moves repair had just reverted.

**A gap this ticket surfaces rather than hides.** The job runner drives the mutation/repair generation loop; tier 1 and tier 2 exist in `cf-search` and are not yet wired into a job — `P1-30`'s `design` command is where the pipeline is assembled. A job whose budget plans those tiers reports them as `skipped` with a `CF-ORCH-006` warning naming what is missing, rather than reporting success over less search than the plan describes.
| `P1-30` | CLI and library | `cf-cli` | `design`, `score`, `feasibility`, `seeds import`, `front`, `replay`, `explain`, `hosts list`; `cf.design()` with an identical argument object | `TST-CLI-001`, `TST-CLI-002`, `TST-ORCH-003` | `P1-29` | L |
| `P1-31` | Output writers | `cf-cli` | FASTA, GenBank, CSV and JSON emitters; the aligned `explain` rendering | Round-trip test: emitted GenBank re-ingests to the same digest | `P1-30` | S |

### `P1-28` — what "minimal" defers, precisely

`P1-28` is not a parallel mechanism alongside the versioning scheme in [ADR-0009](../07-decisions/ADR-0009-artifact-versioning-and-resolution.md) and [`04-model-resolution.md`](../04-api/04-model-resolution.md); it is that scheme's Phase 1 subset. The split, stated so that no later ticket assumes a deferred rung already works:

| Capability | Phase 1 (`P1-28`) | Deferred to |
| --- | --- | --- |
| Manifest identity: coordinate, content digest, `contract_major`, `provides_terms` | owned | — |
| Content-addressed local artifact store | owned | — |
| `pin` rung — resolve by exact digest | owned | — |
| `base` rung — the host-agnostic fallback artifact | owned | — |
| `builtin` rung — terms with no artifact at all | owned | — |
| Lockfile sealing, replay, and `CF-RESOLVE-006` on a missing digest | owned | — |
| Non-blocking recalibration *as observed by a sealed job* | owned | — |
| `tenant_channel` and `host_channel` rungs | **deferred** | `P3-15` |
| Channel promotion, `yanked`, staged rollout | **deferred** | `P3-15` |
| `nearest_host` rung and the host-distance metric | **deferred** | `P3-16` |
| Fallback admissibility, `sigma` inflation, `calibrated: false` | **deferred** | `P3-16` |

The non-blocking-recalibration guarantee is testable now despite channels being deferred, because the guarantee is a property of the *lockfile*, not of the promotion machinery: a job that sealed a digest reads that digest for its whole life regardless of what else enters the store. `P3-15` adds the ability for a promotion to occur; `P1-28` establishes that a sealed job would not see one.

The hazard specific to a partial ladder is silent fall-through — a request for a deferred rung quietly resolving one rung lower and returning a plausible artifact. Every deferred rung therefore raises `CF-RESOLVE-007` naming the ticket that will implement it, rather than falling through. The error also names the artifact it *would* have selected and points at the pin that reaches it today, so the refusal is informative rather than merely correct.

Implementation surfaced one schema defect. `resolution[].resolved_by` has always admitted `builtin`, but `resolution[].artifact_id` required an `ArtifactId`, whose grammar begins `cf:<kind>:` and has no spelling for a term with no artifact — so a lockfile recording a `usage.cai` resolution was unwriteable. Fixed as `DEP-VER-002`. The built-in's `digest` is the hash of its engine descriptor rather than a sentinel: it has real work to do there, since a change to a built-in's engine version can change the numbers and the lockfile has to say so.

---

## Validation

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P1-32` | Synthetic suite | `cf-testkit` | `SYN-01` (toy exhaustive), `SYN-04` (planted violations), `SYN-05` (infeasible sets), `SYN-07` (property strategies) | `TST-SYN-001`, `TST-SYN-004` | `P1-06` | M |
| `P1-33` | Determinism harness | `cf-search` | Per-block/per-individual RNG stream assignment; fixed-order float reduction; the 1/4/16-thread comparison harness | `TST-SRCH-004`, `GATE-DET-001` | `P1-25` | M |
| `P1-34` | `EXP-05` runner | `cf-bench` | Lamarckian vs Baldwinian repair comparison on tightly-constrained inputs | `EXP-05` decided and the default confirmed or changed | `P1-17`, `P1-32` | S |
| `P1-35` | Docs | — | README with a working quickstart; contributor guide covering the two-toolchain setup; `AGENTS.md` updated from its current stub state to describe the real build, test and run commands | A new contributor can build, test and run a design from the README alone | `P1-30` | M |

---

## Phase 1 exit decisions

Behaviour questions that implementation surfaced and that need an explicit
answer before Phase 1 closes. Each currently ships on the **legacy-parity**
setting, so deciding is a configuration change rather than a code change, and
not deciding is a decision to keep the legacy behaviour.

| Id | Decision | Currently | Evidence |
| --- | --- | --- | --- |
| `EXIT-001` | Should motif and restriction-site constraints default to `scope: construct`, catching a site that straddles the coding-region boundary? The prototype scopes them to the CDS and misses it. | legacy parity (`scope: cds` where the differential requires it; `spec::expand` already emits `construct` for assembly standards) | `docs/legacy-parity.md` §1 |
| `EXIT-002` | Should k-mer counting default to `KmerMasking::Assembled`, catching a repeat that straddles a region boundary? The prototype's masked backbone counting misses it. | legacy parity available; `Assembled` is the code default and the differential runs in `SegmentedByRegion` | `docs/legacy-parity.md` §2 |
| `EXIT-003` | `EXP-05`: Lamarckian vs Baldwinian repair. Gated on `P1-34`, which needs `P1-17` (done) and `P1-32`. | Lamarckian | `planning/08-experiments/experiment-register.md` |

`EXIT-001` and `EXIT-002` change which designs are called feasible relative to
the tool this repository already ships, so a construct the prototype accepted
may be rejected. That is the intended improvement, but it is a user-visible
behaviour change and should be decided rather than absorbed.

## Phase 1 exit checklist

- [ ] `EXIT-001`, `EXIT-002`, `EXIT-003` decided
- [ ] `TST-CORE-*`, `TST-CONS-*`, `TST-SYN-*`, `TST-INGE-*`, `TST-LAPI-*`, `TST-LAND-*`, `TST-CLI-*`, `TST-ORCH-*` green
- [ ] `TST-SRCH-001`, `-002`, `-003`, `-004`, `-006` green
- [ ] `GATE-DET-001` (bit-identical across thread counts) green
- [ ] `EXP-05` decided
- [ ] A design job on the repository's own dsbA example runs end-to-end, is feasible, and replays bit-identically from its lockfile
- [ ] Wheels build on all three platforms with no Rust toolchain required at install time
- [ ] Every `DEP-*` marker attributed to Phase 1 re-checked and still justified
