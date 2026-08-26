# Module Dependency Graph

Satisfies: `BRIEF-2.1-b`.
Normative source: [`modules.yaml`](modules.yaml). Rendered here for review.
Checked by: `planning/tools/validate_planning_package.py` (acyclicity + rank monotonicity + dangling-edge detection).

## 1. Graph

```
rank 0   cf-types ────────────────────────────────────────────┐  cf-metrics
            │                                                 │      │
rank 1   cf-landscape-api   cf-constraints   cf-modelreg   cf-testkit │
            │        │            │              │  │              │ │
rank 2   ┌──┴────────┴────┬───────┴──────┬───────┘  └── cf-ingest  │ │
         │                │              │                  │      │ │
     cf-epistasis     cf-burden   cf-biophys-structure       │      │ │
         │                │        cf-biophys-lm             │      │ │
         │                │              │                   │      │ │
rank 3   ├── cf-biophys-kinetics         │                   │      │ │
         │                │              │                   │      │ │
         └────────────────┴──────────────┘                   │      │ │
                (registered at runtime, NOT imported)        │      │ │
                                                             │      │ │
rank 3                    cf-landscape ──────────────────────┘      │ │
                               │                                    │ │
rank 4                    cf-search ────────────────────────────────┘ │
                          cf-decision ────────────────────────────────┘
                               │
rank 5           cf-orchestrator      cf-active
                               │        │
rank 6              cf-api    cf-cli   cf-bench
```

Solid edges are compile-time/import-time dependencies. The line marked *registered at runtime* is deliberately **not** an edge: predictors are discovered through `TermRegistry`, so `cf-landscape` never names them.

## 2. Adjacency (normative)

| Module | rank | depends_on |
| --- | --- | --- |
| `cf-types` | 0 | — |
| `cf-metrics` | 0 | — |
| `cf-landscape-api` | 1 | `cf-types` |
| `cf-constraints` | 1 | `cf-types` |
| `cf-modelreg` | 1 | `cf-types` |
| `cf-testkit` | 1 | `cf-types` |
| `cf-ingest` | 2 | `cf-types`, `cf-modelreg` |
| `cf-epistasis` | 2 | `cf-types`, `cf-landscape-api` |
| `cf-burden` | 2 | `cf-types`, `cf-landscape-api` |
| `cf-biophys-structure` | 2 | `cf-types`, `cf-landscape-api` |
| `cf-biophys-kinetics` | 3 | `cf-types`, `cf-landscape-api`, `cf-epistasis` |
| `cf-biophys-lm` | 2 | `cf-types`, `cf-landscape-api`, `cf-modelreg` |
| `cf-landscape` | 3 | `cf-types`, `cf-landscape-api`, `cf-constraints` |
| `cf-search` | 4 | `cf-types`, `cf-landscape-api`, `cf-constraints`, `cf-metrics` |
| `cf-decision` | 4 | `cf-types`, `cf-metrics`, `cf-landscape-api` |
| `cf-orchestrator` | 5 | `cf-types`, `cf-ingest`, `cf-landscape`, `cf-search`, `cf-constraints`, `cf-decision`, `cf-modelreg` |
| `cf-active` | 5 | `cf-types`, `cf-ingest`, `cf-modelreg` |
| `cf-api` | 6 | `cf-types`, `cf-orchestrator`, `cf-active`, `cf-modelreg` |
| `cf-cli` | 6 | `cf-types`, `cf-orchestrator`, `cf-active`, `cf-modelreg` |
| `cf-bench` | 6 | `cf-types`, `cf-metrics`, `cf-orchestrator`, `cf-testkit`, `cf-ingest` |

`cf-biophys-kinetics -> cf-epistasis` is the only predictor-to-predictor edge. It exists because the v1 kinetics proxy (AMD-002) consumes the epistatic model's per-codon dwell prediction rather than re-deriving one. It is acyclic and deliberate; if the kinetics objective is later replaced by a simulation-based one, the edge disappears.

## 3. The four cycles this design breaks

The architecture proposal describes layers that mutually reference each other. Each such reference is a latent cycle. All four are broken by inverting the dependency onto a leaf or onto an artifact.

### 3.1 `cf-landscape` ↔ predictors

*Naive form:* the composition layer calls each predictor, and each predictor needs the composition layer's caching/uncertainty utilities.

*Break:* `cf-landscape-api` (rank 1) owns the `Predictor` trait, the cache key type, and the uncertainty struct. `cf-landscape` holds `Vec<Box<dyn Predictor>>` populated by a registry at process start. Predictor crates depend on the API, never on the composer.

*Consequence enforced in CI:* `cargo tree -p cf-landscape` must not contain any `cf-epistasis`, `cf-burden`, or `cf-biophys-*`. Test `TST-LAND-003`.

### 3.2 `cf-constraints` ↔ `cf-search`

*Naive form:* "repair must not fight the loop" reads as "repair needs the objective".

*Break:* repair takes an optional `&DeltaTable` argument (type owned by `cf-landscape-api`). The search engine builds it and passes it in. Repair is a pure function; it is fully testable with a hand-written delta table and no search engine present.

### 3.3 `cf-active` ↔ predictors

*Naive form:* active learning retrains the epistatic surface and the burden model, so it imports their training code.

*Break:* `cf-active` writes a `RecalibrationJob` descriptor into `cf-modelreg`. Per-predictor training entry points are separate build targets that poll the registry. The only shared surface is the artifact + manifest. This also means a tenant can run recalibration on separate hardware, which is required for ADR-0006 tenant isolation.

### 3.4 `cf-search` ↔ `cf-decision`

*Naive form:* the search loop terminates on hypervolume stagnation (a decision-support concept); decision support ranks the search engine's population.

*Break:* `cf-metrics` (rank 0, no dependencies) owns non-dominated sorting, crowding, reference directions and hypervolume. Both sides depend on it; neither depends on the other. The search engine emits a serialized `ParetoFront` document; `cf-decision` reads it.

## 4. Runtime-boundary edges

Three edges cross a process boundary and therefore carry a serialization cost that the design must respect (see [ADR-0001](../07-decisions/ADR-0001-language-and-runtime-split.md) and [ADR-0002](../07-decisions/ADR-0002-oracle-vs-differentiable-surrogate.md)):

| Edge | Boundary | Protocol | Frequency |
| --- | --- | --- | --- |
| `cf-search` → `cf-landscape` (in-core terms) | none, same crate graph | direct call | 10⁶–10⁸ per job |
| `cf-landscape` → `cf-biophys-lm` | rust-core → python-worker | batch-yield, Arrow IPC over shared memory | once per generation (10²–10³ per job) |
| `cf-orchestrator` → `cf-search` | python → rust-core | PyO3, GIL released for the duration of a tier | 10¹ per job |

**Rule:** no edge may cross the rust↔python boundary more than once per generation. This rule is what forces the batch-yield protocol instead of a per-move callback, and it is the single most consequential structural decision in the plan. A per-move Python callback at 10⁶ moves/job × ~5 µs GIL round-trip = ~5 s of pure overhead per *thousand* moves, i.e. the search engine would spend >99 % of its time in FFI.

## 5. Verification

`planning/tools/validate_planning_package.py --check deps` asserts:

1. every name in `depends_on` is a declared module (no dangling edges);
2. the relation is a DAG (Kahn's algorithm; reports the cycle if not);
3. `rank(dependency) < rank(module)` for every edge (catches a cycle that a future edit would introduce before it becomes one);
4. every module declares a `public_contract` path that exists;
5. every module declares at least one test id, and every test id appears in `planning/05-testing/01-test-plan.md` (satisfies acceptance criterion §5.8).
