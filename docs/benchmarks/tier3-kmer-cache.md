# Tier 3: the time after the fold cache

The v2.9.2 profile (`v2.9.2-profiles.md`) measured yeeJ CAI + structure,
assembly constraints, budget `fast`, generations 2, seed 7, with
`maturin develop` (`opt-level = 1`). Tier 3 was 3.202 s of a 5.855 s job.
The fold cache was already at 1 full fold, 21 local edits, 12 hits. This
note is where the rest of that tier 3 time went, measured on a release
build (`maturin develop --release`) so the before and after share one
profile. The documented vector is unchanged.

Host `escherichia_coli:511145/BL21`. Weights
`examples/mg1655_ribosomal_usage.json`. Constraints
`benchmarks/validation/constraints_assembly.json`. ViennaRNA 2.7.2.
`CF_TIER3_PROFILE=1` prints one JSON row per phase of `run_generations`.

## Where the 3 seconds went

Population 12, two generations. Dominance sort and hypervolume are noise
at N = 24. The L0 batch is two Python crossings, one per generation
(`ADR-0001`). Repair of the infeasible offspring is the multi-second cost:
each synonym trial calls `ConstraintChecker::check` on 7077 nt, and
`check_kmers` was building a fresh 12-mer index plus a `Vec` of positions
for every distinct 12-mer.

| phase | calls | before | after |
| --- | ---: | ---: | ---: |
| `breed_check` | 24 | 0.050 s | 0.020 s |
| `feasibility_gate` | 24 | 0.050 s | 0.028 s |
| `repair` | 6 | 2.013 s | 0.746 s |
| `evaluate` | 24 | 0.050 s | 0.015 s |
| `l0_batch` | 2 | 0.656 s | 0.658 s |
| `select` | 2 | 0.000 s | 0.000 s |
| `hypervolume` | 2 | 0.000 s | 0.000 s |
| wall | | 5.125 s | 3.060 s |

The vector on both sides is `[0.9575892496987667, 2.408954059933615]`.
Fold counts on both sides are hits 12, full folds 1, local edits 21.
`python/tests/test_profile_identity.py` pins those.

Skipping the pre-mutation check when the parent was feasible would change
which children enter repair: two feasible parents can breed an infeasible
child, and the violation list is what `violation_directed` reads. That
edit was not made.

## What shipped

1. `check_kmers` records positions only for k-mers already over the limit.
   The offender set is the one the count map already named. The violation
   list, including `at [positions]`, stays the same.
2. A thread-local `KmerIndex` is updated with `substitute` for a single
   span of at most 16 bases, then rewound before the function returns. A
   longer edit, or two spans, rebuilds. `KmerIndex::substitute` matches a
   rebuild (`property_tst_cons_001_incremental_matches_rebuild`). The
   rewind is what `unit_a_codon_edit_reuses_the_kmer_index_and_restores_it`
   checks: leaving the index on the probe makes the next synonym look like
   two spans.

Repair is still about as long as the fold. Motif, GC, and homopolymer
scans inside the same check were left alone.

## Primer response

The core search document for yeeJ with 300 nt flanks is still 166429 bytes.
The HTTP job that embeds it was 169009 bytes, of which the 20 amplicons
are 136549 characters. The default `GET /v1/primers/{id}` is 31980 bytes.
`?product_sequence=1` returns the 169009 byte document. One pair goes from
8480 bytes to 1613. The recommended pair is unchanged.
