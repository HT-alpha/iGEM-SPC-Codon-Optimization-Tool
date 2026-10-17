# Query-Time Model Version Selection

Satisfies: `BRIEF-2.4-c` — "how host-specific and dataset-specific model versions are selected at query time".
Decision record: [ADR-0009](../07-decisions/ADR-0009-artifact-versioning-and-resolution.md).
Owning module: `cf-modelreg`.

---

## 1. Resolution happens once per job, not once per query

"At query time" is the requirement; resolving *per query* would be a mistake. A tier-2 sweep issues on the order of 10⁶ landscape queries, and a registry lookup on each one would be both slow and — much worse — non-deterministic, because a channel could move mid-job and different queries in the same run would use different models.

So: the ladder below runs **once, when the job leaves `queued`**, and its full result is written into a [`LandscapeLockfile`](../02-schemas/landscape_lockfile.schema.json). Every subsequent query carries `model_pins` derived from that lockfile. The lockfile is what makes results reproducible, what makes recalibration non-blocking, and what makes "which model produced this number" answerable months later.

---

## 2. The resolution ladder

For each requested `term_id`, the first rung that yields a compatible artifact wins.

| Rung | Source | `resolved_by` | Notes |
| --- | --- | --- | --- |
| 1 | Explicit `model_pins` in the request, or a `replay_lockfile_digest` | `pin` | Highest priority. A replay that cannot find a pinned digest **fails**; it never substitutes. |
| 2 | Tenant channel for `(kind, host_id)` | `tenant_channel` | The tenant's own recalibrated adapter, if one has been promoted. |
| 3 | Public channel for `(kind, host_id)` | `host_channel` | The normal path for a supported host. |
| 4 | Nearest-host artifact within the distance threshold | `nearest_host` | **Only if `allow_host_fallback = true`.** Sets `fallback: true` and inflates uncertainty. |
| 5 | Host-agnostic base artifact | `base` | e.g. the LM naturalness model, which is genuinely host-agnostic (host conditioning is an input, not a weight set). |
| — | Built-in, non-learned computation | `builtin` | CAI, GC, k-mer counts, the native structure kernel. Not an artifact; still recorded, with the engine version. |

If no rung yields an artifact, the job is **rejected at submission** with `CF-RESOLVE-001` naming the term and the host. It is not silently dropped from the objective set: a job that optimizes four objectives when the user asked for five, without saying so, is worse than a job that fails.

### 2.0 Rungs not yet implemented

Phase 1 (`P1-28`) implements rungs 1, 5 and the built-in row. Rungs 2, 3 and 4 need the channel machinery (`P3-15`) and the host-distance metric (`P3-16`). Until then, a request that would resolve through an unimplemented rung raises `CF-RESOLVE-007` naming the rung and the ticket that will provide it.

`DEP-LAND-005` — `epistasis.context_naturalness` does not take the nearest-host rung. A job whose host id starts with `escherichia_coli:` pins the public MG1655 artifact. No distance penalty is applied, because none was measured. Any other host is `CF-LAND-001`.

This is deliberate rather than defensive. A partially-implemented ladder whose missing rungs fall through to the next one still returns an artifact, and the artifact is plausible — a `base` model where the caller asked for a host-specific one is exactly the kind of substitution that produces a result nobody can explain later. The whole point of rung ordering is that a lower rung is a *worse* answer, so silently taking it while reporting success is the failure mode the ladder exists to prevent.

### 2.1 Compatibility filter, applied at every rung

```
COMPATIBLE(artifact, request) :=
      artifact.contract_major == request.contract_major
  AND artifact.lifecycle.status in {active, deprecated}      # yanked never resolves
  AND term_id in artifact.provides_terms
  AND (artifact.isolation.scope == public
       OR artifact.isolation.tenant_id == request.tenant_id)
  AND artifact.capability[term_id] >= minimum required by the requested tiers
```

The last clause matters: if the job requests tier-2 block coordinate ascent, a term whose only available provider is L0 forces the planner either to drop that term from tier 2 (default, with a warning) or to fail. Discovering that mid-run would waste the whole tier.

`deprecated` still resolves — deprecation is advice, not removal — but adds a `CF-RESOLVE-005` warning naming the successor. `yanked` never resolves through a channel but stays retrievable by digest forever, so a published result remains reproducible.

---

## 3. Nearest-host fallback

### 3.1 Distance

```
d(h, h') = w_tax · rank(LCA(h, h'))
         + w_cu  · (1 − cos(RSCU(h), RSCU(h')))
         + w_trna· (1 − cos(gcn_profile(h), gcn_profile(h')))

rank:  strain 0 | species 1 | genus 2 | family 3 | order 4 | class 5 | phylum 6 | domain 7 | root 8
defaults: w_tax = 0.5, w_cu = 3.0, w_trna = 3.0
```

Codon usage and tRNA gene-copy profiles are weighted far above taxonomy on purpose. Taxonomic proximity is a poor predictor of translational similarity — two enterobacteria can differ markedly in GC content and therefore in codon usage — whereas the RSCU vector is a direct measurement of the thing being transferred. Taxonomy is retained as a tie-break and as an interpretable label for the user ("used the species-level model").

### 3.2 Admissibility

A fallback is admissible only if **all** of:

- `allow_host_fallback = true` in the request (default is `false`);
- `rank(LCA) ≤ 3` (family or closer);
- `cos(RSCU) ≥ 0.90`.

Otherwise the term is unavailable and the job is rejected. The default being `false` is deliberate: silently scoring an unusual host with another organism's model is exactly the kind of failure a user must opt into.

### 3.3 Uncertainty inflation

```
uncertainty_inflation_factor = 1 + κ_tax · rank(LCA) + κ_cu · (1 − cos(RSCU))
```

`κ_tax` and `κ_cu` are **not guessed**. `EXP-07` measures them: for every pair of hosts where both have held-out data, evaluate host `h`'s model on host `h'`'s data, measure the actual degradation in predictive interval coverage, and fit the inflation factors so that a fallback's stated 90 % interval really covers 90 %. Until `EXP-07` runs, conservative literature-free defaults (`κ_tax = 0.15`, `κ_cu = 2.0`) apply and `Uncertainty.calibrated` is set to `false` for any fallback term — which propagates into probabilistic dominance ([ADR-0005](../07-decisions/ADR-0005-uncertainty-and-pareto-ranking.md)) and downgrades those candidates' robustness scores automatically.

Every fallback also attaches the `low_confidence_host` rationale tag and a `CF-RESOLVE-004` warning in `job_status.warnings`.

---

## 4. Algorithm

```
RESOLVE(request) -> LandscapeLockfile

  lock <- new lockfile(contract_major = request.contract_major)

  if request.replay_lockfile_digest is set:
      src <- FETCH_LOCKFILE(request.replay_lockfile_digest)
      for entry in src.resolution:
          if not AVAILABLE(entry.digest):
              fail CF-RESOLVE-006 naming entry.artifact_id and entry.digest
          lock.resolution.append(entry)                # verbatim, no re-resolution
      return lock

  for term in request.objectives ∪ request.diagnostics:
      cand <- none
      for rung in [pin, tenant_channel, host_channel, nearest_host, base, builtin]:
          if rung == nearest_host and not request.allow_host_fallback: continue
          c <- LOOKUP(rung, term, request.host_id, request.tenant_id)
          if c is none or not COMPATIBLE(c, request): continue
          if rung == nearest_host and not ADMISSIBLE(c): continue
          cand <- c; break
      if cand is none:
          fail CF-RESOLVE-001 naming term and host
      lock.resolution.append(entry(term, cand, rung, fallback, host_distance))

  # a single job may not mix contract majors: term values would not be
  # comparable, and the objective vector would be meaningless
  assert all entries share contract_major else fail CF-RESOLVE-002

  lock.engines      <- RESOLVE_ENGINES(request.budget)
  lock.host_profile <- PIN_HOST_PROFILE(request.host_id)
  lock.inputs       <- DIGESTS(construct, constraints, seeds, sequences)
  lock.rng          <- RNG_PLAN(request.seed, request.budget)
  lock.environment  <- ENVIRONMENT_SNAPSHOT()
  return SEAL(lock)                                    # content-addressed
```

---

## 5. Channels

A channel is a small JSON document naming one artifact id and digest for a `(scope, kind, host_id)` triple.

```
registry/channels/public/epistasis_dwell/escherichia_coli:562%2FK-12_MG1655.json
registry/channels/t_igem2026lab01/burden_capacity/escherichia_coli:562%2FK-12_MG1655.json
```

- Channel names: `stable` (default), `candidate`, `nightly`.
- Promotion is an atomic write of that document. There is no half-promoted state.
- Promotion requires every required gate in the artifact's `evaluation.gates` to be `passed: true`, and requires `code_version.dirty == false` for public channels.
- Channel history is retained, so "which artifact was `stable` on 2026-07-01" is answerable.

---

## 6. Caching

Resolution results are cached per `(tenant, host, term, channel)` with a short TTL (60 s) purely to keep submission latency low. The cache is **never** consulted during a job — jobs read the lockfile. A stale cache can therefore only affect which artifact a *new* job pins, never the consistency of a running one.

---

## 7. Tests

| Test | Assertion |
| --- | --- |
| `TST-MREG-001` | Ladder order: a fixture registry with artifacts at every rung resolves in the documented order, and each `resolved_by` value is recorded correctly. |
| `TST-MREG-002` | Replay fidelity: replaying a lockfile with a missing digest fails with `CF-RESOLVE-006` and never substitutes a newer artifact. |
| `TST-MREG-003` | Fallback admissibility and inflation: inadmissible fallbacks are refused; admissible ones set `fallback: true`, inflate `sigma` by the computed factor, and set `calibrated: false`. |
| `TST-MREG-004` | Contract-major mixing is refused (`CF-RESOLVE-002`). |
| `TST-MREG-005` | A `yanked` artifact never resolves through a channel but always resolves by digest. |
