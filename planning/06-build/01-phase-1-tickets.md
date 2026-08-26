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
| `P1-17` | Repair operator | `cf-constraints` | Depth-1 exhaustive search; depth-`d` beam; the `(edits, objective cost, random)` ordering key; escalation; bounded give-up | `TST-CONS-002`, `TST-CONS-003`, `TST-CONS-004` | `P1-14`–`P1-16` | L |
| `P1-18` | Tabu ledger | `cf-constraints` | Per-run `(position, codon)` ledger with generation-clock decay and the deliberate 10 % leak | `TST-CONS-005` | `P1-17` | S |
| `P1-19` | Feasibility analyser | `cf-constraints` | All five checks, each reporting a constraint to relax and by how much | `TST-CONS-006` | `P1-14`–`P1-16` | M |

---

## Landscape skeleton

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P1-20` | Landscape contract types | `cf-landscape-api` | `LandscapeQuery`/`Response` in Rust and Python; the three target shapes; the L0/L1/L2 capability model; `TermRegistry` | `TST-LAPI-001`, `TST-LAPI-002` | `P1-08` | M |
| `P1-21` | Built-in terms | `cf-landscape` | `usage.cai`, `usage.gc_content`, `constraint.violation_count` — all L1 with exact deltas | `TST-EPI-002`-style delta test applied to built-ins | `P1-20` | M |
| `P1-22` | Composition and cache | `cf-landscape` | Term registry resolution; window-keyed evaluation cache namespaced by lockfile digest; `partial`/`UNKNOWN` propagation | `TST-LAND-001`, `TST-LAND-002`, `TST-LAND-003` | `P1-21` | M |

---

## Search

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P1-23` | Genotype and operators | `cf-search` | `Vec<u8>` genotype with cached hash; the four mutation operators; block-uniform crossover; duplicate suppression; seed injection | `TST-SRCH-006` | `P1-20` | M |
| `P1-24` | Tier 1 constrained beam DP | `cf-search` | Beam DP on the product with the motif automaton; diversity-capped beam selection | `TST-SRCH-002` (tier-1 portion) | `P1-14`, `P1-23` | L |
| `P1-25` | Tier 2 red-black block CA | `cf-search` | The partition (`B ≥ W` assertion); parallel same-colour block updates; incumbent-seeded beam for monotonicity; 1-opt polish; right-boundary correction | `TST-SRCH-001`, `TST-SRCH-002` | `P1-24` | XL |
| `P1-26` | Seam reconciliation fallback | `cf-search` | Overlapping-block path with the exact 3-state chain DP and seam polish | `TST-SRCH-003` | `P1-25` | M |
| `P1-27` | Baseline backtracking solver | `cf-search` | Port of the existing `codon_6.0.py` stack-based backtracker as `baseline_backtracking`, kept as a reference implementation and a differential-test oracle for the constraint layer | Differential test: for inputs the prototype solves, the new constraint checker agrees on feasibility for 10⁴ sequences | `P1-14`–`P1-16` | M |

`P1-27` is worth doing rather than discarding the prototype: it is an independently-written implementation of the same constraint semantics, which makes it a genuine oracle for `P1-14`–`P1-16`. It also gives a documented answer to "does the new tool still do what the old script did".

> **Sequencing correction (implementation, 2026-08).** `P1-27` originally declared `P1-17` (the repair operator) as its dependency. That was over-stated: the acceptance criterion exercises only the constraint *checker*, and the legacy `count_overlapping` semantics it exists to pin down are implemented in `P1-15`. The dependency is therefore `P1-14`–`P1-16`, and `P1-27` is scheduled **immediately after `P1-16`, before `P1-17`** — a differential oracle that lands after the code it is meant to check has already been built on top of cannot catch drift at the point the drift would be introduced.

---

## Registry, orchestration, surfaces

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P1-28` | Minimal registry | `cf-modelreg` | Manifest and lockfile types; local filesystem artifact store; the `pin`, `base` and `builtin` rungs (channels arrive in Phase 3) | `TST-MREG-002` | `P1-08` | M |
| `P1-29` | Job lifecycle | `cf-orchestrator` | In-process job execution; lockfile sealed **before** any search work; RNG plan; checkpoint/resume; cooperative cancel | `TST-ORCH-001`, `TST-ORCH-002` | `P1-25`, `P1-28` | L |
| `P1-30` | CLI and library | `cf-cli` | `design`, `score`, `feasibility`, `seeds import`, `front`, `replay`, `explain`, `hosts list`; `cf.design()` with an identical argument object | `TST-CLI-001`, `TST-CLI-002`, `TST-ORCH-003` | `P1-29` | L |
| `P1-31` | Output writers | `cf-cli` | FASTA, GenBank, CSV and JSON emitters; the aligned `explain` rendering | Round-trip test: emitted GenBank re-ingests to the same digest | `P1-30` | S |

---

## Validation

| Id | Title | Module | Deliverable | Acceptance | Deps | Size |
| --- | --- | --- | --- | --- | --- | --- |
| `P1-32` | Synthetic suite | `cf-testkit` | `SYN-01` (toy exhaustive), `SYN-04` (planted violations), `SYN-05` (infeasible sets), `SYN-07` (property strategies) | `TST-SYN-001`, `TST-SYN-004` | `P1-06` | M |
| `P1-33` | Determinism harness | `cf-search` | Per-block/per-individual RNG stream assignment; fixed-order float reduction; the 1/4/16-thread comparison harness | `TST-SRCH-004`, `GATE-DET-001` | `P1-25` | M |
| `P1-34` | `EXP-05` runner | `cf-bench` | Lamarckian vs Baldwinian repair comparison on tightly-constrained inputs | `EXP-05` decided and the default confirmed or changed | `P1-17`, `P1-32` | S |
| `P1-35` | Docs | — | README with a working quickstart; contributor guide covering the two-toolchain setup; `AGENTS.md` updated from its current stub state to describe the real build, test and run commands | A new contributor can build, test and run a design from the README alone | `P1-30` | M |

---

## Phase 1 exit checklist

- [ ] `TST-CORE-*`, `TST-CONS-*`, `TST-SYN-*`, `TST-INGE-*`, `TST-LAPI-*`, `TST-LAND-*`, `TST-CLI-*`, `TST-ORCH-*` green
- [ ] `TST-SRCH-001`, `-002`, `-003`, `-004`, `-006` green
- [ ] `GATE-DET-001` (bit-identical across thread counts) green
- [ ] `EXP-05` decided
- [ ] A design job on the repository's own dsbA example runs end-to-end, is feasible, and replays bit-identically from its lockfile
- [ ] Wheels build on all three platforms with no Rust toolchain required at install time
- [ ] Every `DEP-*` marker attributed to Phase 1 re-checked and still justified
