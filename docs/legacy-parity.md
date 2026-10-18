# Legacy parity: `codon_6.0.py` and the DNA Design Suite constraint checker

Ticket `P1-27`. Test: `python/tests/test_p1_27_legacy_differential.py` (`TST-CONS-007`).

The prototype in `codon_6.0.py` is an independently-written implementation of
the same constraint semantics the new checker implements, which makes it a
genuine oracle rather than a formality. `P1-27` runs a differential over 10 000
synonymous variants and requires exact agreement on feasibility.

**Result: 10 000 / 10 000 agreement.** The comparison is not one-sided — 5 468
of the variants are infeasible and 4 532 feasible under the four parameter
settings used, so agreement is a real signal rather than "both said yes".

Reaching that took resolving three places where the two disagreed. All three are
legacy behaviours, all three are preserved so the differential is meaningful,
and none of them is the default in the new checker. They are recorded here for a
decision rather than diverged from unilaterally.

---

## 1. Motifs are matched against the coding sequence only

`optimize_cds_sequence` checks each forbidden motif against `test_seq`, the
partial *coding sequence*. It never sees the backbone. A restriction site that
straddles the boundary between the coding region and the vector is therefore
invisible to the prototype, even though the synthesised molecule contains it.

```
  index  0 1 2 3 | 4  5  6  7  8  9
  base   A A G C | T  T  G  G  G  G
         \--- AAGCTT ---/
                  \------ CDS -----/
```

**Preserved as:** `scope: cds` on the constraint, which is an ordinary schema
value and needs no special machinery.

**Default in the new checker:** `scope: construct`, which catches it.

**Assessment:** almost certainly a prototype bug. A BsaI site spanning the
junction breaks a Golden Gate assembly exactly as much as one three bases
further in. **Recommend defaulting assembly-standard constraints to construct
scope**, which is what `spec::expand` already does.

---

## 2. Backbone k-mer counts are masked, so straddling repeats vanish

`optimize_whole_construct` builds backbone counts from the insert with every
coding region replaced by `N`, and `build_kmer_frequency_counts` skips any k-mer
containing an `N`. The within-CDS count is then added separately. A k-mer
spanning the boundary is counted by neither side.

```
  index  0 1 2 3 4 5 6 7 | 8 9 10 11
  base   A C G T A A A C | G T  C  C
         \--ACGT--/  \--ACGT--/
```

**Preserved as:** `KmerMasking::SegmentedByRegion` on the *view*, not on the
constraint — it describes what the checker is looking at, not what the user
asked for.

**Default in the new checker:** `KmerMasking::Assembled`.

**Assessment:** a prototype bug of the same kind as (1), and for the same
reason. The masking exists to stop the CDS being double-counted, and dropping
junction k-mers is collateral rather than intent.

---

## 3. Only k-mers occurring in the coding sequence are examined

This is the one the differential actually surfaced, and unlike the first two it
is **not** a bug.

The prototype's k-mer loop iterates over windows of the coding sequence:

```python
for i in range(len(test_seq) - k + 1):
    kmer = test_seq[i:i+k]
    total = backbone_kmer_counts[kmer] + count_overlapping(test_seq, kmer)
```

A k-mer that repeats between two *backbone* regions and never appears in the
coding sequence is never looked at. The first version of the new checker counted
it, and disagreed on 1 452 of 10 000 variants — every one of them a repeat
confined to the immutable vector, for example `TGGGCCTT` at positions 184 and
255 with the coding region at 17–146.

The prototype is asking the right question. The search engine cannot edit the
backbone, so a violation confined to it is a property of the *input*, not of the
candidate. Counting it against every candidate would report every design as
infeasible and give the search nothing to optimize toward.

**Resolution:** `Violation::actionable` records whether any editable position is
involved. `ConstraintChecker::is_feasible` considers only actionable blocking
violations; `ConstraintChecker::unfixable` returns the rest. Nothing is hidden —
the violation is still reported, and it is what the feasibility analyser
(`P1-19`, check 1, `CF-FEAS-001`) reports once at submission with the constraint
to relax.

This distinction was implicit in the planning package — `P1-19`'s locked-span
check assumes it — but was not stated in the constraint contract. It now is.

---

## Decisions requested

| # | Legacy behaviour | New default | Recommendation |
| --- | --- | --- | --- |
| 1 | Motifs scoped to the CDS | Construct scope | **Diverge.** Confirm. |
| 2 | Junction k-mers uncounted | Assembled counting | **Diverge.** Confirm. |
| 3 | Backbone-only repeats unexamined | Reported, non-actionable | **Keep legacy's question**, but report the violation. No confirmation needed — this preserves the legacy outcome. |

Items 1 and 2 change which designs are called feasible, so a construct the
prototype accepted may be rejected by DNA Design Suite. That is the intended
improvement, but it is a behaviour change against the tool this repository
already ships, and it should be an explicit decision rather than a side effect.

Until confirmed, both legacy behaviours remain reachable through configuration
and are exercised by the differential test, so reverting either is a
configuration change rather than a code change.
