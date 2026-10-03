# REST API Design

Satisfies: `BRIEF-2.4-a`.
Normative document: [`openapi.yaml`](openapi.yaml).
Owning module: `cf-api`.

---

## 1. Surface at a glance

| Method | Path | Purpose |
| --- | --- | --- |
| `POST` | `/v1/jobs` | Submit a design job |
| `GET` | `/v1/jobs` | List jobs |
| `GET` | `/v1/jobs/{id}` | Poll status, progress and live hypervolume |
| `GET` | `/v1/jobs/{id}/front` | Retrieve the Pareto front |
| `GET` | `/v1/jobs/{id}/lockfile` | Retrieve the reproducibility lockfile |
| `GET` | `/v1/jobs/{id}/events` | SSE progress stream |
| `POST` | `/v1/jobs/{id}/cancel` | Cooperative cancel, front still emitted |
| `POST` | `/v1/score` | Score existing sequences, no design |
| `POST` | `/v1/feasibility` | Constraint feasibility analysis |
| `POST` | `/v1/seeds:import` | Import external tool output as seeds |
| `POST` | `/v1/observations` | Submit a wet-lab readout |
| `GET` | `/v1/observations/{id}` | Ingestion state of a submission |
| `GET` | `/v1/recalibrations` | Recalibration runs and their gate reports |
| `GET` | `/v1/hosts`, `/v1/hosts/{id}`, `/v1/hosts/{id}/models` | Host data and resolvable models |
| `GET` | `/v1/models`, `/v1/models/{id}` | Artifact manifests |
| `GET` | `/v1/version`, `/v1/healthz` | Meta |

---

## 2. Design decisions worth defending

### 2.1 Async by default, with the lockfile available before the work starts

`POST /v1/jobs` returns `202` immediately. What matters is *when* `landscape_lockfile_digest` becomes readable: as soon as the job leaves `queued`, and therefore long before the expensive search. A client can inspect exactly which model versions will be used, and abort, before paying for a run. Resolution is cheap (a registry lookup per term); the search is not.

### 2.2 `GET /front` is separate from `GET /jobs/{id}`

A front with per-position rationale for 20 candidates of a 15 kb construct is on the order of a megabyte. Polling status every two seconds must not transfer it. Status carries a *summary* (`current_front_summary`, including the hypervolume history for a live plot); the front is fetched once.

### 2.3 `/v1/score` exists so that three consumers cannot diverge

"What is my current sequence worth", "what is LinearDesign's output worth on our axes", and "what is the before column of the comparison table" are the same computation. Giving them one endpoint means the benchmark numbers in [`../05-testing/04-benchmark-harness.md`](../05-testing/04-benchmark-harness.md) are produced by the same code path a user hits, which is the only way a benchmark claim stays honest.

### 2.4 `/v1/feasibility` is a first-class endpoint, not a job side effect

Finding out that a constraint set is unsatisfiable should take milliseconds and should name the constraint to relax. Discovering it after a 15-minute job returns nothing is the failure mode this endpoint exists to remove — and it is exactly what the current prototype does when it exhausts its restart budget.

### 2.5 Cancellation still yields a result

`POST /cancel` stops at the next generation boundary and emits the best certified front found so far. A user who cancels because "this is taking too long" almost always still wants what the search has already found.

### 2.6 Errors carry a stable code

Every error is `application/problem+json` with a `code` matching `^CF-[A-Z]{3,8}-[0-9]{3}$`. Clients branch on `code`; prose is for humans and may be reworded at any time.

| Code family | Meaning |
| --- | --- |
| `CF-SEQ-*` | Input sequence problems: not in frame, unknown base, internal stop, protein mismatch on a seed |
| `CF-HOST-*` | Unknown host, or a requested term the host has no data to support |
| `CF-FEAS-*` | The five feasibility checks of [`../03-algorithms/06-constraint-repair.md`](../03-algorithms/06-constraint-repair.md) §5 |
| `CF-RESOLVE-*` | Model resolution: no artifact, contract-major mismatch, fallback used but not permitted |
| `CF-BUDGET-*` | Requested budget not satisfiable; downgraded (warning) or refused (error) |
| `CF-LAND-*` | Landscape provider errors, including a cache-evicted `base_digest` |
| `CF-JOINT-*` | Joint-search decomposition failed to converge; fell back |
| `CF-OBS-*` | Observation validation: as-built digest mismatch, missing control, undeclared censoring |
| `CF-AUTH-*` | Authentication and tenancy |

### 2.7 Warnings are structured, not prose in a log

`job_status.warnings` is an array of `ProblemDetails` with `status: 200`. Anything that changed the meaning of the result without failing the job appears there: a host fallback was used, a requested term was dropped, the budget tier was downgraded to fit the hardware, a seed was rejected, the joint decomposition fell back. These are the conditions a user must see; burying them in server logs is how a user ends up trusting a number produced by another organism's model.

---

## 3. Tenancy and isolation

The bearer token carries the tenant. **No endpoint accepts a tenant id in a request body or path.** Cross-tenant addressing is therefore not expressible, which is the API-level half of the isolation invariant in [ADR-0006](../07-decisions/ADR-0006-active-learning-data-isolation.md); the storage-level half is artifact-store partitioning.

Consequences that fall out of this and are worth stating so nobody re-litigates them later:

- `GET /v1/models` returns public artifacts plus the caller's tenant artifacts, never another tenant's.
- `model_pins` may name a public artifact or one of the caller's own; naming a foreign artifact returns `CF-AUTH-003`, not `404`, because the distinction between "does not exist" and "not yours" is not information a caller is entitled to.
- A lockfile from tenant A cannot be replayed by tenant B if it pins a tenant-A artifact. The replay fails with `CF-RESOLVE-006` naming the unavailable digest.

---

## 4. Versioning

Three version axes, deliberately independent:

| Axis | Where | Changes when |
| --- | --- | --- |
| **URL major** (`/v1`) | path | An incompatible change to the HTTP surface. Old majors are served in parallel for a documented window. |
| **`schema_version`** | every document body | A schema evolves. Additive changes bump minor; clients ignore unknown fields. |
| **`contract_major`** | landscape query/response | The meaning of a landscape term changes. Artifacts declare which they satisfy; mixing majors is refused at resolution. |

An additive term (say `structure.mfe_global` gaining a sibling) bumps none of them: term ids are namespaced strings and adding one is not a breaking change. This is why `requested_terms` has no wildcard — a wildcard would make adding a term a breaking change for every existing client.

---

## 5. Rate limiting and quotas

Design jobs are the expensive resource, so limits are expressed in job-seconds rather than requests: a per-tenant concurrent-job cap and a rolling compute-seconds quota, both returned in `X-CodonForge-Quota-*` headers. Read endpoints get a conventional request-rate limit. Exceeding either returns `429` with `Retry-After` and code `CF-QUOTA-001`.

---

## 6. This build

`DEP-API-001`. The server is the Python standard library. `POST /v1/jobs` takes the library `DesignRequest`. A body whose `input.kind` is `protein`, `construct`, or `single_cds`, or whose `objectives_preset` is not `none`, returns `CF-CLI-003` rather than a translated job. `GET /v1/hosts`, `/v1/models`, `/v1/observations`, `/v1/recalibrations`, `/v1/seeds:import`, and `POST /v1/feasibility` return `CF-CLI-003` with a reason. Feasibility of a finished candidate is `constraint_report.feasible` on the front. `GET /` serves the result page, which renders those documents and does not embed scores. A job id the caller cannot see is `CF-AUTH-003` for both a missing id and another tenant's id.

## 7. Tests

| Test | Assertion |
| --- | --- |
| `TST-API-001` | Every response body validates against its declared schema, driven from the OpenAPI document — so a schema change that the API does not follow fails CI. |
| `TST-API-002` | Idempotency: the same `Idempotency-Key` with the same body returns the original job; with a *different* body it returns `409 CF-IDEM-001` rather than silently ignoring the difference. |
| `TST-API-003` | Tenancy: a token for tenant A cannot read, pin or replay any tenant-B artifact, job, observation or lockfile. Enumerates every endpoint. |
| `TST-API-004` | Every error path returns `application/problem+json` with a `code` matching the pattern, and every code used is documented in §2.6. |
