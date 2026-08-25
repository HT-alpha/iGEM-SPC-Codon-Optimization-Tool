# Repository Layout

Satisfies: `BRIEF-2.1-a`.
Normative machine-readable companion: [`modules.yaml`](modules.yaml).
Language/runtime rationale: [ADR-0001](../07-decisions/ADR-0001-language-and-runtime-split.md).

## 1. Design rules the layout enforces

| Rule | Mechanism |
| --- | --- |
| One module per architectural layer | Each layer in the architecture proposal maps to exactly one directory with its own manifest (`Cargo.toml` or `pyproject.toml` entry) and its own test target. |
| Independently testable | Every module has a sibling `tests/` directory and can be tested with the rest of the workspace absent, using `cf-testkit` fakes for its dependencies. CI runs a per-module test job. |
| Independently swappable | Predictor layers are reachable only through `cf-landscape-api` trait objects resolved at runtime from a registry. Nothing above the predictors names a predictor crate. |
| No circular dependencies | `modules.yaml` declares the edges; `planning/tools/validate_planning_package.py` fails the build if the relation is not a DAG or if a module depends on an equal-or-higher rank. |
| Reproducible artifacts | Trained weights never live in the source tree. They are content-addressed blobs referenced by manifests under `registry/`. |

## 2. Directory tree

```
codonforge/
├── Cargo.toml                     # Rust workspace root
├── pyproject.toml                 # Python distribution `codonforge`
├── AGENTS.md
├── README.md
├── LICENSE                        # Apache-2.0 (see ADR-0001 §6)
├── NOTICE                         # third-party attributions
│
├── crates/                        # ---- rust-core runtime -----------------
│   ├── cf-types/
│   │   ├── src/{lib,codon,seqrec,tables,ids,units,error}.rs
│   │   ├── data/genetic_codes/*.json       # NCBI translation tables 1..33
│   │   └── tests/
│   ├── cf-metrics/
│   │   ├── src/{lib,nds,crowding,refdirs,hypervolume,igd}.rs
│   │   └── tests/
│   ├── cf-landscape-api/
│   │   ├── src/{lib,query,response,capability,batch,registry}.rs
│   │   └── tests/
│   ├── cf-constraints/
│   │   ├── src/{lib,spec,detect,repair,tabu,feasibility}.rs
│   │   └── tests/
│   ├── cf-epistasis/
│   │   ├── src/{lib,potts,dwell,delta,artifact}.rs
│   │   └── tests/
│   ├── cf-burden/
│   │   ├── src/{lib,ledger,wobble,fixedpoint,penalty}.rs
│   │   └── tests/
│   ├── cf-biophys-structure/
│   │   ├── src/{lib,window,engine,accessibility}.rs
│   │   ├── src/engines/{native,viennarna_ffi,linearfold_ffi}.rs   # latter two feature-gated
│   │   └── tests/
│   ├── cf-biophys-kinetics/
│   │   ├── src/{lib,pause_profile,domain_map,concordance}.rs
│   │   └── tests/
│   ├── cf-landscape/
│   │   ├── src/{lib,compose,cache,uncertainty,batch_driver}.rs
│   │   └── tests/
│   ├── cf-search/
│   │   ├── src/{lib,tier1_seed,tier2_blockca,tier3_moea}.rs
│   │   ├── src/operators/{mutation,crossover,dedup,seed_inject}.rs
│   │   ├── src/joint/{ledger_share,price_update,convergence}.rs
│   │   └── tests/
│   └── cf-pyext/                  # PyO3 glue -> `codonforge._core`
│       └── src/lib.rs
│
├── python/codonforge/             # ---- python runtimes -------------------
│   ├── __init__.py
│   ├── _core.pyi                  # typed stubs for the Rust extension
│   ├── ingest/{fasta,genbank,construct,host,riboseq,seeds}/
│   ├── modelreg/{manifest,registry,resolve,lockfile,store}.py
│   ├── biophys_lm/{model,score,train,adapters}.py
│   ├── epistasis/train/{plm,dwell,features,cv}.py
│   ├── burden/calibrate/{capacity,proteome_prior}.py
│   ├── decision/{dominance,robustness,reduce,rationale,export}.py
│   ├── orchestrator/{jobs,workers,pin,checkpoint,budget}.py
│   ├── active/{ingest,validate,schedule,gate,exploration}.py
│   ├── api/{app,routes,auth,tenancy,errors}.py
│   ├── cli/{main,design,front,models,observe,bench}.py
│   ├── bench/{harness,datasets,baselines,indicators,gates}.py
│   └── testkit/{synth_host,synth_potts,toy_space,planted,strategies}.py
│
├── worker/                        # ---- python-worker runtime -------------
│   └── cf_worker/{main,batch,gpu,terms}.py
│
├── registry/                      # ---- artifact registry (no weights) ----
│   ├── manifests/<kind>/<host_id>/<semver>.json
│   ├── channels/<tenant>/<kind>/<host_id>.json
│   └── datasets/<dataset_id>.json
│
├── schemas/                       # normative JSON Schemas, copied from
│   └── v1/*.schema.json           # planning/02-schemas at Phase-1 kickoff
│
├── planning/                      # this package (frozen after Phase 1 start)
├── benchmarks/
│   ├── datasets/<dataset_id>/manifest.json
│   └── results/<run_id>/
├── docs/
└── .github/workflows/{ci,bench-nightly,bench-weekly,release}.yml
```

## 3. Layer-to-module map

| Architecture layer | Module | Runtime | Swap seam |
| --- | --- | --- | --- |
| Context ingestion | `cf-ingest` | python-orchestration | `SourceAdapter` per file format; `SeedImporter` per external tool |
| Epistatic fitness surface | `cf-epistasis` | rust-core (fit in Python) | `TermId = epistasis.*`; any L0/L1 provider may replace it |
| Biophysical: structure | `cf-biophys-structure` | rust-core | `StructureEngine` trait; native / ViennaRNA / LinearFold |
| Biophysical: folding kinetics | `cf-biophys-kinetics` | rust-core | `TermId = kinetics.*` |
| Biophysical: LM naturalness | `cf-biophys-lm` | python-worker | `TermId = lm.*`; HF checkpoint id in manifest |
| Host resource / burden model | `cf-burden` | rust-core | `TermId = burden.*`; `CapacityEstimator` trait |
| Search engine (three tiers) | `cf-search` | rust-core | `Tier` trait; tiers can be disabled per budget |
| Hard-constraint repair | `cf-constraints` | rust-core | `RepairStrategy` trait |
| Pareto decision support | `cf-decision` | python-orchestration | `Reducer` / `RationaleTagger` plugins |
| Active-learning feedback | `cf-active` | python-orchestration | `ReadoutAdapter` per assay type |

## 4. Independent testability contract

Every module ships three test tiers, named so CI can select them:

- `unit_*` — no I/O, no fixtures larger than 4 KB, `< 100 ms` each.
- `contract_*` — validates the module's public contract against the JSON Schema in `planning/02-schemas/`, using golden examples. Guarantees a module can be swapped without breaking callers.
- `property_*` — invariants under randomized input (`proptest` in Rust, `hypothesis` in Python).

A module's test target may depend only on `cf-testkit` plus its own declared `depends_on`. CI enforces this by building each module in an isolated workspace slice.

## 5. Departures from the architecture proposal

| ID | Departure | Rationale |
| --- | --- | --- |
| `DEP-LAY-001` | The architecture names "biophysical predictors" as one layer; the layout splits it into three sibling modules. | Structure, kinetics and LM naturalness have different runtimes (two Rust, one GPU-Python) and completely different training/validation regimes. One module would force the GPU dependency on every install. |
| `DEP-LAY-002` | Two support modules with no architectural counterpart (`cf-landscape-api`, `cf-metrics`) are introduced as leaves. | They exist purely to break the two dependency cycles the architecture implies (landscape↔predictors, search↔decision). Without them the "swappable layer" goal is unachievable. |
| `DEP-LAY-003` | `registry/` holds manifests only; weights are external content-addressed blobs. | Keeps the repository clonable without LFS and makes tenant isolation (ADR-0006) a storage-level property rather than a repository-level one. |
