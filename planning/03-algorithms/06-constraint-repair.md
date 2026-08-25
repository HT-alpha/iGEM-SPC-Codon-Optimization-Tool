# Algorithm (f) — The In-Loop Constraint Repair Operator

Satisfies: `BRIEF-2.3-f`, including "how it interacts with (not fights) the evolutionary loop".
Owning module: `cf-constraints`.
Contract: [`../02-schemas/constraint_set.schema.json`](../02-schemas/constraint_set.schema.json) (`repair_policy`).

---

## 1. What "fighting the loop" actually means

A naive repair operator degrades a multi-objective search in four specific, separately-fixable ways. Naming them is the whole design:

| Failure | Mechanism | Fix in this design |
| --- | --- | --- |
| **Undoing progress** | Repair rewrites a large region, discarding the improvement that variation just made. | Minimal edit distance: repair changes the fewest codons that clear the violation. §3 |
| **Diversity collapse** | Deterministic repair maps many distinct individuals onto the same canonical feasible point, so the population converges to one lineage in a few generations. | Randomized tie-breaking among equally-good repairs. §3.3 |
| **Thrashing** | Mutation proposes edit `E`, repair reverts it, next generation mutation proposes `E` again — forever. | A tabu ledger of `(position, codon)` pairs reverted by repair, consulted by the mutation operator for `tabu_tenure_generations`. §4 |
| **Starvation** | Unbounded repair search consumes the generation's compute budget. | Hard bounds on window, depth, beam and escalations, with a defined give-up that leaves the individual infeasible rather than searching harder. §3.4 |

The give-up path is not a failure mode, it is the design: an individual that cannot be cheaply repaired stays in the population as **infeasible**, where NSGA-III's constraint-domination rule ([`03-nsga-loop-and-operators.md`](03-nsga-loop-and-operators.md) §4) already handles it correctly. Infeasible individuals are useful — they are stepping stones across constraint barriers. Forcing feasibility on every individual destroys them.

---

## 2. Incremental violation detection

Repair is only affordable if detection is incremental. Each constraint class gets an incremental detector, and the per-substitution cost is what makes the whole operator viable:

| Constraint class | Structure | Cost after one codon substitution |
| --- | --- | --- |
| `forbidden_motif`, `restriction_site`, `custom_regex` | Aho–Corasick automaton over all patterns (both strands, IUPAC-expanded) | `O(m_max)` — rescan a window of `m_max + 3` nt |
| `homopolymer` | run-length index | `O(max_run)` |
| `kmer_repeat` | hash map `k-mer → count`, maintained incrementally | `O(k)` — a 3 nt change removes `k+2` and adds `k+2` k-mers |
| `gc_window` | prefix sums over GC | `O(window/3)` windows touched, `O(1)` each |
| `local_structure` | incremental window folding cache | `O(windows touched)` — the expensive one; see [`08`](08-complexity-and-performance-budgets.md) |
| `codon_ban`, `min_codon_frequency` | direct lookup | `O(1)` |

Full detection over a fresh sequence is `O(n)` for everything except `local_structure`.

**Overlapping counting is the default** for `kmer_repeat`, matching `count_overlapping()` in the existing `codon_6.0.py` prototype. Python's `str.count` is non-overlapping and undercounts tandem repeats, which are exactly the structures a synthesis vendor rejects. The schema exposes `count_overlapping` so the semantics are explicit rather than inherited by accident.

---

## 3. The repair operator

### 3.1 Signature and purity

```
repair(sequence, constraint_set, delta_table: Option<&DeltaTable>, rng, tabu)
    -> RepairOutcome { sequence, edits, residual_violations, objective_cost }
```

`cf-constraints` has no dependency on `cf-search` or on any predictor. The exact delta table is *passed in* by the caller. That is what makes repair independently testable (hand-write a delta table, assert the chosen repair) and what keeps the dependency graph acyclic ([`../01-modules/02-module-dependency-graph.md`](../01-modules/02-module-dependency-graph.md) §3.2).

### 3.2 Algorithm

```
REPAIR(x, C, Δ, rng, tabu, policy) -> outcome

  V <- DETECT_BLOCKING(x, C)
  if V is empty: return x unchanged
  edits <- {}

  for v in ORDER_BY_SPAN_THEN_SEVERITY(V):
      if not STILL_VIOLATED(v, x): continue          # an earlier repair fixed it

      r <- policy.max_window_codons
      for escalation in 0 .. policy.max_escalations:
          P <- editable, unlocked codon positions within  span(v) expanded by r
          fix <- SEARCH_MINIMAL_FIX(x, v, P, C, Δ, rng, policy)
          if fix found:
              x <- APPLY(x, fix)
              edits <- edits ∪ fix
              break
          r <- 2·r                                    # escalate the window

      # no fix within budget: leave it. NSGA-III handles infeasibility.

  residual <- DETECT_BLOCKING(x, C)
  TABU_RECORD(tabu, edits)                            # see §4
  return outcome(x, edits, residual, objective_cost = Σ Δ over edits)


SEARCH_MINIMAL_FIX(x, v, P, C, Δ, rng, policy) -> Option<edit set>

  # Depth 1 first: overwhelmingly the common case for restriction sites and
  # homopolymers, and it is exhaustive at ~70 candidates.
  cands <- []
  for p in P:
      for y in Syn(a_p) \ {x_p}:
          if (p, y) in tabu: continue
          x' <- x with x_p := y
          if STILL_VIOLATED(v, x'): continue
          if CREATES_NEW_BLOCKING(x', C, near p): continue
          cands.append( (p, y, cost = Δ_scalarized(p, y)) )
  if cands non-empty:
      return ARGMIN_WITH_RANDOM_TIEBREAK(cands, rng, policy)

  # Depth >= 2: bounded beam over edit sets.
  beam <- [ (edits = {}, x, score = 0) ]
  for depth in 2 .. policy.max_depth:
      next <- []
      for st in beam:
          for p in P, y in Syn(a_p) \ {st.x_p}:
              if (p, y) in tabu: continue
              st' <- extend(st, p, y)
              nviol <- COUNT_VIOLATIONS(st'.x, C, near span(v))
              next.append( st' with key = (nviol, st'.cost, rng.tiebreak()) )
      beam <- TOP_B(next, policy.beam_width, key)
      for st in beam:
          if not STILL_VIOLATED(v, st.x) and not CREATES_NEW_BLOCKING(st.x, C):
              return st.edits
  return None
```

### 3.3 Ordering key, and why every component of it is there

Repairs are ranked lexicographically by `(number of edits, scalarized objective cost, random)`:

1. **Number of edits first.** This is the "do not undo the search's work" property. A one-codon fix that costs a little objective is preferred over a three-codon fix that costs nothing, because the three-codon fix discards more of what variation just produced.
2. **Objective cost second**, from the caller-supplied exact delta table. This is what makes repair *objective-aware* rather than blind. When `policy.objective_aware = false` (or `Δ` is `None`, e.g. in unit tests), this component is dropped and repair falls back to edits-then-random.
3. **Random third.** Not a tie-break of convenience — this is the diversity mechanism. With deterministic tie-breaking, every individual carrying the same violation receives the same repair, and after a few generations the population shares a large identical region.

### 3.4 Bounds

| Bound | Default | Consequence of exceeding it |
| --- | --- | --- |
| `max_window_codons` | 12 | Widen and retry (`max_escalations` times) |
| `max_escalations` | 2 | Give up on this violation |
| `max_depth` | 3 | Give up on this violation |
| `beam_width` | 32 | Truncate the beam |

Worst-case work per violation is therefore bounded by `max_escalations × max_depth × beam_width × |P| × S × detect_cost`, a constant known before the run — which is what keeps repair from starving the generation.

### 3.5 Lamarckian by default

The repaired genotype replaces the original (`repair_policy.lamarckian = true`). Rationale: every solution the user can actually build lies on the feasible manifold, so carrying infeasible genotypes forward spends population slots on regions that cannot yield a deliverable. The Baldwinian alternative — evaluate the repaired form but inherit the unrepaired genotype — preserves diversity better and is available as a config flag; it is not the default because it doubles the effective repair cost (repair runs every generation on the same individual) without a demonstrated benefit. `EXP-05` measures the difference on the benchmark set and the default is revisited if Baldwinian wins.

---

## 4. The tabu ledger — the concrete anti-thrashing mechanism

```
struct TabuLedger {
    entries: HashMap<(CodonIndex, Codon), GenerationsRemaining>,
}

TABU_RECORD(tabu, edits):
    for (p, y_new) in edits:
        y_old = the codon that was there before repair
        tabu[(p, y_old)] = policy.tabu_tenure_generations     # discourage re-proposing it

TABU_TICK(tabu):                       # called once per generation by cf-search
    for key in tabu: tabu[key] -= 1
    remove entries at 0

# consulted by the mutation operator:
MUTATE(...):
    ...
    if (p, y) in tabu and rng.uniform() < 0.9: resample      # 10% leak, deliberate
```

The **10 % leak is deliberate**. A hard tabu can make a genuinely good region unreachable if the constraint that triggered the tabu was only violated in a particular context that no longer applies. The leak lets the population re-test a tabooed edit occasionally, at a rate low enough not to re-establish the thrashing cycle.

The ledger is per-run, shared across the population, and decays on a generation clock — so it is a soft, self-clearing memory rather than a permanent restriction.

---

## 5. Feasibility analysis before the search starts

The existing prototype's failure mode is instructive: after exhausting 1 000 restarts it prints `Error: Exceeded max restarts (1000) without finding a solution.` and stops. The user learns nothing about *why*, or which constraint to relax. CodonForge runs a feasibility analysis at submission and returns a specific, actionable diagnosis instead.

```
ANALYSE_FEASIBILITY(protein, flanks, C) -> FeasibilityReport

  1. LOCKED-SPAN CHECK
     If any blocking violation lies entirely inside locked/immutable sequence
     (backbone, flanks, locked codons), it is unfixable by construction.
     -> CF-FEAS-001, naming the constraint and the span.

  2. COUNTING BOUND (k-mer repeats)
     A construct of length n contains n-k+1 overlapping k-mers drawn from an
     alphabet of 4^k. If  max_occurrences · 4^k  <  n-k+1,  the constraint is
     unsatisfiable by pigeonhole for ANY sequence.
     -> CF-FEAS-002, with the minimum feasible k or max_occurrences.

  3. PER-POSITION EXHAUSTION
     For each editable position, does at least one synonym avoid creating an
     unavoidable violation given the immutable context?
     -> CF-FEAS-003, naming the position.

  4. AUTOMATON REACHABILITY
     Run tier-1 beam DP with regular constraints only and no objective. If the
     beam empties at position i, the motif constraints alone are jointly
     unsatisfiable from that prefix.
     -> CF-FEAS-004, naming position i and the patterns that closed off every
        continuation.

  5. GC-WINDOW ATTAINABILITY
     For each window, compute min and max achievable GC over synonym choices
     (a small per-window DP). If the required band is outside [min, max],
     the window is unsatisfiable.
     -> CF-FEAS-005, reporting the attainable band.
```

Every check reports a **constraint to relax and by how much**, not just a failure. Steps 2 and 5 are exact; steps 1, 3 and 4 are necessary-condition checks, so passing all five does not guarantee feasibility — the report says so.

`DEP-CONS-001` — The architecture proposal does not mention feasibility analysis. It is added here because a hard-constraint layer without one produces unexplainable empty results, which is the most common way a codon optimizer becomes unusable in practice.

---

## 6. Complexity

`|P|` = candidate positions ≈ `2·max_window_codons + span` ≈ 26; `S` ≈ 2.9; `detect` ≈ 50 flops (incremental, all classes except `local_structure`).

| Operation | Cost | Value |
| --- | --- | --- |
| Depth-1 exhaustive search | `O(|P| · S · detect)` | `26 × 2.9 × 50 ≈ 3.8 × 10³` flops |
| Depth-`d` beam search | `O(d · beam · |P| · S · detect)` | `3 × 32 × 26 × 2.9 × 50 ≈ 3.6 × 10⁵` flops (worst case) |
| Per individual per generation | `× V` violations, typically `V ∈ {0,1,2}` | `≈ 0–7 × 10⁵` flops |
| Per generation | `× N = 212` | `≈ 1.5 × 10⁸` flops worst case, `~0` typical |
| Feasibility analysis | `O(n · m_max) + O(L · b · S)` | milliseconds, once per job |

Worst case (every individual carrying two hard violations needing depth-3 repair) is ~1.5×10⁸ flops per generation, comparable to one exact L0 evaluation batch. Typical case, after the first few generations, is near zero because the population has moved onto the feasible manifold. The bounds in §3.4 are what guarantee the worst case is a constant rather than a search.

---

## 7. Tests

| Test | Assertion |
| --- | --- |
| `TST-CONS-001` | Detection parity: incremental detectors agree exactly with a naive `O(n)` rescan over 10⁵ random edits, for every constraint class. |
| `TST-CONS-002` | Minimality: on planted single-site violations with a known unique minimal fix, repair returns exactly that fix. On planted violations with several minimal fixes, repair returns one of them, and over 1 000 runs with different seeds it returns more than one distinct fix (the diversity property of §3.3). |
| `TST-CONS-003` | Termination and bounds: repair never exceeds its declared work bound, and always terminates, over an adversarial suite including overlapping and mutually-conflicting violations. |
| `TST-CONS-004` | Protein and lock invariance: repair never changes the translated protein and never touches a locked codon, over 10⁶ randomized cases. |
| `TST-CONS-005` | Anti-thrashing: in a scripted mutation/repair cycle designed to thrash, the tabu ledger breaks the cycle within `tabu_tenure` generations. |
| `TST-CONS-006` | Feasibility analyser: each of the five checks fires on a purpose-built infeasible input and stays silent on a feasible one. |
