# Schema Conventions and Evolution Rules

Satisfies `BRIEF-4.a` — "every schema should be given in a concrete serializable form (e.g. JSON Schema, protobuf, or typed dataclass definitions) so it can be reviewed line by line".

---

## 1. Inventory

| Schema | Contract for | Brief |
| --- | --- | --- |
| [`common.schema.json`](common.schema.json) | Shared primitives: identifiers, digests, units, uncertainty, provenance, problem details | — |
| [`sequence.schema.json`](sequence.schema.json) | The canonical codon-indexed, frame-safe sequence | 2.2-a |
| [`landscape_query.schema.json`](landscape_query.schema.json) | Search engine → landscape | 2.2-b |
| [`landscape_response.schema.json`](landscape_response.schema.json) | Landscape → search engine | 2.2-b |
| [`host_profile.schema.json`](host_profile.schema.json) | tRNA pools, codon usage, capacity, data completeness | 2.2-c |
| [`riboseq_track.schema.json`](riboseq_track.schema.json) | Ribosome-profiling training data | 2.2-c |
| [`construct_context.schema.json`](construct_context.schema.json) | Multi-cistron layout, RBS/promoter metadata | 2.2-c |
| [`constraint_set.schema.json`](constraint_set.schema.json) | Hard and advisory constraints, repair policy | 2.2-c |
| [`pareto_front.schema.json`](pareto_front.schema.json) | The design job's deliverable | 2.2-c |
| [`rationale_tags.json`](rationale_tags.json) | The closed per-position rationale vocabulary | 2.2-c |
| [`model_manifest.schema.json`](model_manifest.schema.json) | Trained artifact provenance and gates | 2.2-d |
| [`landscape_lockfile.schema.json`](landscape_lockfile.schema.json) | Reproducibility and version pinning | 2.2-d, 2.4-b, 2.4-c |
| [`observation.schema.json`](observation.schema.json) | Wet-lab readout ingestion | 2.4-b |
| [`seed_set.schema.json`](seed_set.schema.json) | External tool output as seeds | 3.6 |
| [`job_request.schema.json`](job_request.schema.json) | `POST /v1/jobs` and `cf.design()` | 2.4-a |
| [`job_status.schema.json`](job_status.schema.json) | `GET /v1/jobs/{id}` | 2.4-a |
| [`interfaces/predictors.pyi`](interfaces/predictors.pyi) | Typed interface stubs for the predictor plugin seam | 2.1, 2.2-b |

Every schema has a worked example in [`examples/`](examples/sequence.example.json), and the validator fails if one does not: **an unexercised contract is an unreviewed contract.**

**The examples' digests are genuine, not placeholders.** Every digest occurrence is classified in [`examples/digest_map.yaml`](examples/digest_map.yaml) as computed (derived from a payload in the repository), a reference (must equal a digest computed elsewhere), or external (payload not in the repository, with a stated reason). `planning/tools/refresh_example_digests.py` recomputes and threads them in dependency order; the `digests` check group fails on any digest the map does not account for, which is what stops a new placeholder from quietly reappearing.

---

## 2. Conventions

| Convention | Rule | Why |
| --- | --- | --- |
| Dialect | JSON Schema draft 2020-12 | Widest tooling support and the dialect OpenAPI 3.1 uses, so the API document can `$ref` these files directly |
| `$id` | Absolute `https://codonforge.org/schemas/v1/<name>.schema.json`; `$ref`s are relative | Relative refs resolve against the `$id` base, so the files work both on disk and when served |
| `additionalProperties` | `false` on every closed object | Catches typos at the boundary. Extensibility is provided by explicit `metadata` / `parameters` objects where it is genuinely needed |
| Descriptions | Every non-obvious field carries a description that says **why the field exists**, not what its type is | The type is already in the schema. The reason is what a reviewer needs |
| Spans | Half-open `[start, end)`. Codon-indexed unless the field name ends in `_nt` | Removes an entire class of off-by-one bugs by convention |
| Digests | `sha256:` prefix, over RFC 8785 canonical JSON for documents and raw bytes for blobs | Canonicalization is what makes a digest stable across languages and key orderings |
| Enums over booleans | `enum` where a third state is plausible | `censoring: none/left/right/interval` rather than `is_censored`; the third state always arrives eventually |
| Units | Mandatory on every numeric term in a landscape response | A unitless objective vector is not comparable across model versions |
| Nullable | `type: [X, "null"]` for "known to be absent"; field omission for "not computed" | These are different, and `UNKNOWN` is a first-class value in the search |

---

## 3. Evolution rules

Three independent version axes, deliberately not coupled:

| Axis | Field | Bump when |
| --- | --- | --- |
| Document schema | `schema_version` (SemVer) | Additive change → minor; removing or retyping a field → major |
| Landscape contract | `contract_major` | The *meaning* of a term changes: new/renamed term, changed unit, changed sign convention |
| API surface | URL path (`/v1`) | An incompatible change to the HTTP surface |

**Additive changes are minor and clients must ignore unknown fields.** Adding a new `TermId` bumps nothing — term ids are namespaced strings, and adding one is not a breaking change. This is precisely why `requested_terms` has no wildcard: a wildcard would make adding a term a breaking change for every existing client and every existing lockfile.

**Any change to a schema goes through this package first**: edit here, run `validate_planning_package.py`, then propagate to `schemas/v1/` in the implementation tree. A drift check in CI (`P1-03`) fails if a Rust struct or `.pyi` stub diverges from its schema. The schemas stay normative rather than degrading into documentation of whatever the code happens to do.

---

## 4. Validation

```bash
python3 planning/tools/validate_planning_package.py -v
```

Checks, in order: every schema is valid draft 2020-12; every example validates against its schema; every schema has an example; the sequence invariants hold on the worked examples; the module dependency graph is a DAG with no rank inversion and no dangling edge; every `ADR-`/`AMD-`/`EXP-`/`RISK-`/`GATE-`/`TST-`/`DEP-` identifier cited anywhere is defined in its owning document; every relative Markdown link resolves; the OpenAPI document parses and its external `$ref`s resolve; every brief requirement maps to an artifact that exists; every rationale tag used in an example is in the closed vocabulary with the right scope.

Requires `jsonschema>=4.18` and `pyyaml`.
