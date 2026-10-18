# DNA Design Suite — full performance report (v2.9.4)

Standalone export of the 17 October 2026 measurement. The repository copy is `docs/benchmarks/v2.9.4-full-performance-report.md`.

Measured on 17 October 2026 against commit `b8065f71d99f0cddc5968e3904b97d5f07b937e4`
(PR #33, “Harsh performance suite, and a fix for rpoC assembly stalling”).
Package version is 0.1.0. No tag was cut. Optimizer behaviour was left as
that commit shipped it.

The machine-readable rows sit beside this file:
[`DNA-Design-Suite-performance-rows-2026-10-17.json`](DNA-Design-Suite-performance-rows-2026-10-17.json).
An HTML rendering is `DNA-Design-Suite-performance-report-2026-10-17.html`.

## Executive summary

On a release build (Cargo `opt-level` 3, thin LTO), DNA Design Suite and Primer design
match every published output pin from the v2.9.4 harsh suite: 76 rows, no pin
failures, and the same CAI values, CDS hashes, primer pairs, ΔG values, operon
joint costs, and the CAI-only `GATE-BEN-001` claim of 50 of 50. The opt-level 1
build that earlier notes used returns those same outputs. Headline wall times
on this 4-vCPU Xeon put the long path in seed repair. yeeJ (7077 nt) is 0.526 s
with empty constraints (CAI 1), 2.727 s under the assembly rules (CAI
0.9727484991594185, of which seed repair is 1.764 s), 4.203 s under the tighter harsh set, and
4.012 s at budget `standard` with 4 generations, which does not change that
assembly CDS. CAI plus structure on yeeJ is 3.458 s at the pinned vector
`[0.9723247970586831, 2.4166863375261896]`. Primer design on the 7077 nt yeeJ
template with 300 nt flanks is 0.052 s and returns the pinned pair
`TGACCGACAAATATTATG` / `GAACTACAGTGGTTTAAA`, product 6843, ΔG −1.283850000000001.
A single design stays near 53 MiB RSS. Three API designs serialized by
`Api.run_lock` peak at 63112 KiB. waaS under the unpublished harsh constraint
set still stops at `CF-SRCH-006` in 0.059 s.

## Environment

| item | value |
| --- | --- |
| commit | `b8065f71d99f0cddc5968e3904b97d5f07b937e4` |
| package | codonforge 0.1.0 (`codonforge._core.version()`) |
| OS | Ubuntu 24.04.4 LTS, kernel 6.12.94+ |
| CPU | `Intel(R) Xeon(R) Processor`, family 6 model 207, 4 vCPU, KVM |
| RAM | 15 GiB |
| rustc | 1.98.1 (`48a229cea`, 2026-09-01) |
| Python | 3.12.3 |
| maturin | 1.15.0 |
| ViennaRNA | 2.7.2, imported as `RNA` |
| host id | `escherichia_coli:511145/BL21` |
| weights | `examples/mg1655_ribosomal_usage.json` |
| translation table | 11 |
| seed | 7 |
| stagnation window | 21 where the harness sets it |

Headline numbers below are Cargo profile `release`: `opt-level` 3 (the Cargo
default; `Cargo.toml` does not override it), `lto = "thin"`,
`codegen-units = 1`, `panic = "abort"`. The installed extension
`python/codonforge/_core.abi3.so` was byte-identical to
`target/release/libcf_pyext.so` while these rows were taken.

The comparison column marked opt-level 1 is Cargo profile `dev` with
`opt-level = 1`, which is what `make build` (`maturin develop`) produces.
`docs/benchmarks/v2.9.2-profiles.md`, `v2.9.3-profiles.md`, and
`v2.9.4-harsh-perf.md` were measured on that profile. Both profiles were
built and run on this machine on the same day.

Each harness row is its own process. Wall time is `time.perf_counter` around
the call the harness names (the design, the seed, or `primer_search`). It
does not include interpreter startup. Peak RSS is `ru_maxrss` for that
process. Linux reports kilobytes. CAI in the harsh tables is
`benchmarks/validation/cai_independent.py` unless a cell says “engine vector”.

## Method and commands to reproduce

From the repository root, after `make setup` and
`.venv/bin/pip install 'ViennaRNA==2.7.2'`.

Opt-level 1, the profile the older notes use:

```sh
.venv/bin/maturin develop
.venv/bin/python tools/harsh_perf.py --suite all --out /tmp/cf-perf/opt1-harsh.json
.venv/bin/python tools/profile_suite.py --suite all --out /tmp/cf-perf/opt1-profiles.json
.venv/bin/python tools/challenging_perf.py --suite sweep --out /tmp/cf-perf/opt1-sweep.json
.venv/bin/python tools/challenging_perf.py --suite standard \
  --from /tmp/cf-perf/opt1-sweep.json --out /tmp/cf-perf/opt1-standard.json
```

Release, the headline profile:

```sh
.venv/bin/maturin develop --release
.venv/bin/python tools/harsh_perf.py --suite all --out /tmp/cf-perf/rel-harsh.json
.venv/bin/python tools/profile_suite.py --suite all --out /tmp/cf-perf/rel-profiles.json
.venv/bin/python tools/challenging_perf.py --suite sweep --out /tmp/cf-perf/rel-sweep.json
.venv/bin/python tools/challenging_perf.py --suite standard \
  --from /tmp/cf-perf/rel-sweep.json --out /tmp/cf-perf/rel-standard.json
```

`harsh_perf.py --suite all` is pins, design, seed, objectives, operon,
primers, the API job path, and `GATE-BEN-001`. It exits 1 when any pin's
`ok` is false. The timeout is 300 s per row. `profile_suite.py` is stage
timings, seed counters, primer regions, structure, operon, cProfile, and
Potts. `challenging_perf.py` is the v2.9.1 sweep (budget `fast`; structure
jobs use 2 generations) and a `standard` rerun of the four slowest design
jobs at 4 generations.

On this machine the release pass took 91 s (harsh), 37 s (profile), 26 s
(sweep), and 14 s (standard). The opt-level 1 pass took 103 s, 44 s, 29 s,
and 16 s. Every suite exited 0. No row timed out.

Sequences are the shipped validation FASTA files plus `rpoB`, `rpoC`,
`yfhM`, `thrA`, and `aceB` read from
`benchmarks/corpus/NC_000913.3_cds.fna.gz`. No sequence was invented.
Assembly rules are `benchmarks/validation/constraints_assembly.json`.
Golden Gate and BioBrick are `examples/golden_gate.json` and
`examples/biobrick_rfc10.json`. The harsh set is that assembly file plus a
30 nt GC window in [0.35, 0.65], a homopolymer cap of 5, and a 10-mer cap
of 2 (`constraint_set_id` `cst_harsh_assembly`). It is a stress set. It is
not a published feasibility claim.

## DNA Design Suite

Budget `fast` with generations 0 runs tier 1 and tier 2. A structure
objective is L0, so those jobs use 2 generations. Budget `standard` with
generations 4 adds tier 3. Every succeeded design in these tables translated
back to the input protein and had zero actionable violations under the
constraint set it was given.

The same yeeJ assembly job (fast, generations 0) was launched by four
release harnesses. Walls were 2.727 s (harsh design), 2.460 s (harsh pin),
2.532 s (profile stages), and 2.371 s (challenging sweep). The CDS hash was
`72bd74389f99f3aba3dfcb5422c05c91db606adaacb65f67a5a80f881dc6e540` every
time. Treat a tenth of a second on a multi-second row as run noise. The
tables use one named row each.

### Empty constraints, CAI, budget fast, generations 0

Source: release `harsh_perf.py` pins. CAI 1 means every counted codon is on
its heaviest synonym. rpoC starts with GTG. The start stays locked, and
table 11 reads GTG as valine, so empty-constraint CAI is the closed form
`exp(log(0.245665) / 1363) = 0.998970606154696`. That row is pinned to the
closed form.

| gene | nt | wall s | opt-1 s | RSS KiB | CAI |
| --- | ---: | ---: | ---: | ---: | ---: |
| rpsL | 375 | 0.046 | 0.051 | 47488 | 1 |
| elastin_vpgvg20 | 306 | 0.069 | 0.045 | 47564 | 1 |
| gfp_BBa_E0040 | 717 | 0.064 | 0.065 | 47476 | 1 |
| mrfp_BBa_E1010 | 675 | 0.062 | 0.065 | 47488 | 1 |
| mcherry_BBa_J06504 | 711 | 0.065 | 0.066 | 47488 | 1 |
| bla | 861 | 0.078 | 0.081 | 47740 | 1 |
| waaS | 936 | 0.071 | 0.074 | 47620 | 1 |
| aceB | 1602 | 0.132 | 0.141 | 48068 | 1 |
| yagF | 1968 | 0.164 | 0.176 | 48128 | 1 |
| thrA | 2463 | 0.193 | 0.241 | 48240 | 1 |
| lacZ | 3075 | 0.239 | 0.262 | 48400 | 1 |
| rpoB | 4029 | 0.307 | 0.333 | 48704 | 1 |
| rpoC | 4224 | 0.322 | 0.379 | 48740 | 0.998970606154696 |
| yfhM | 4962 | 0.380 | 0.443 | 49012 | 1 |
| yeeJ | 7077 | 0.526 | 0.563 | 50740 | 1 |

aceB and thrA contain native `AGG`. With empty constraints the CAI-only run
removes it and scores 1. The v2.9.1 note recorded the same outcome.

### Assembly constraints, CAI

Assembly is GC windows, homopolymers, repeats, and the Golden Gate / BioBrick
sites in `constraints_assembly.json`. The floor column is the v2.9.4
search-escape CAI (`docs/benchmarks/v2.9.4-search-escape.md`). Every floor
gene meets it. rpoB, rpoC, yfhM, and lacZ are recorded results from this
tree. They were not in that floor table.

Fast rows are budget `fast`, generations 0. Standard rows are budget
`standard`, generations 4. Where both were run, the CDS hash matches, so
the extra generations did not move the recommended sequence.

| gene | nt | budget | wall s | opt-1 s | seed s | RSS KiB | CAI | vs floor |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | --- |
| rpsL | 375 | fast | 0.070 | 0.083 | 0.022 | 47644 | 0.9866866273575934 | equal |
| elastin_vpgvg20 | 306 | fast | 0.421 | 0.578 | 0.367 | 47500 | 0.4642192192077335 | equal |
| gfp_BBa_E0040 | 717 | fast | 0.128 | 0.167 | 0.056 | 47888 | 0.9725744807921104 | equal |
| mrfp_BBa_E1010 | 675 | fast | 0.124 | 0.164 | 0.054 | 47740 | 0.9605695971330996 | equal |
| mcherry_BBa_J06504 | 711 | fast | 0.145 | 0.183 | 0.071 | 47816 | 0.9644182772231985 | equal |
| bla | 861 | fast | 0.174 | 0.217 | 0.088 | 47820 | 0.9847528752532632 | equal |
| waaS | 936 | fast | 0.243 | 0.313 | 0.152 | 47948 | 0.9259163095438141 | equal |
| yagF | 1968 | fast | 0.597 | 0.863 | 0.397 | 48688 | 0.9453398274862173 | equal |
| lacZ | 3075 | fast | 0.636 | 0.860 | 0.345 | 49772 | 0.9665531384697095 | recorded |
| rpoB | 4029 | fast | 0.665 | 0.826 | 0.293 | 51020 | 0.9830194072110617 | recorded |
| rpoC | 4224 | fast | 1.197 | 1.469 | 0.768 | 51368 | 0.9593415607455033 | recorded |
| yfhM | 4962 | fast | 1.519 | 1.570 | 0.924 | 51664 | 0.9605558574962288 | recorded |
| yeeJ | 7077 | fast | 2.727 | 2.825 | 1.764 | 53332 | 0.9727484991594185 | equal |
| yfhM | 4962 | standard ×4 | 3.081 | 3.716 | 0.831 | 52184 | 0.9605558574962288 | same CDS |
| yeeJ | 7077 | standard ×4 | 4.012 | 4.749 | 1.679 | 53704 | 0.9727484991594185 | same CDS |

Seed seconds are `resolve_feasible_seed` wrapped inside that same design
row, so they sit under the wall time beside them. A direct seed call, with
no tier 1 or tier 2 around it, is a bit lower: yeeJ assembly 1.655 s, yfhM
0.794 s, rpoC 0.785 s, waaS 0.151 s, elastin 0.377 s. yeeJ fast hash
`72bd74389f99f3aba3dfcb5422c05c91db606adaacb65f67a5a80f881dc6e540` (733
edits, 3 cycles). waaS hash
`d227b303a3c477cea7f5ed35c3a0fbcb783cb9169fa38a3737e0fe7699419f16`. rpoC hash
`7e80dcaa96a48b33e0078568ff2d6e5638a3cdedc6d6b001d5401ebe66397ca7`. yfhM hash
`69f1b9d902f6ffd54bb3391eac1f6fbfa6f27901a49ff2089a7f40befdd621be`. lacZ hash
`7d83594101155870cad004b71cd93765bad9fa27fd5d9d7dd96dbd2d551a8fad`.

Elastin’s engine vector is `0.4642192192077335`. The independent checker
prints `0.46421921920773357`, one unit in the last place higher. The pin
accepts that. See the correctness section.

Inside the yeeJ fast design (release, 2.727 s): `resolve_feasible_seed`
1.764 s, tier 2 polish 0.253 s, tier 1 seed 0.112 s. The standard run adds
tier 3 at 1.485 s, of which `repair` is 5 calls and 0.556 s. The recommended
CDS does not change.

### Golden Gate and BioBrick

Genes at or under 2500 nt where the v2.9.3 local dynamic program is a full
optimum. Release harsh pins. Gap is DP CAI minus achieved CAI. Every gap is 0.
Wall times are all under 0.2 s. Opt-level 1 is within a few hundredths.

| gene | nt | set | wall s | CAI | gap |
| --- | ---: | --- | ---: | ---: | ---: |
| rpsL | 375 | golden_gate | 0.044 | 1 | 0 |
| rpsL | 375 | biobrick_rfc10 | 0.045 | 0.9904970702046313 | 0 |
| elastin_vpgvg20 | 306 | golden_gate | 0.044 | 1 | 0 |
| elastin_vpgvg20 | 306 | biobrick_rfc10 | 0.044 | 1 | 0 |
| gfp_BBa_E0040 | 717 | golden_gate | 0.066 | 1 | 0 |
| gfp_BBa_E0040 | 717 | biobrick_rfc10 | 0.065 | 0.9952923787903273 | 0 |
| mrfp_BBa_E1010 | 675 | golden_gate | 0.071 | 1 | 0 |
| mrfp_BBa_E1010 | 675 | biobrick_rfc10 | 0.067 | 0.9793498378226791 | 0 |
| mcherry_BBa_J06504 | 711 | golden_gate | 0.068 | 1 | 0 |
| mcherry_BBa_J06504 | 711 | biobrick_rfc10 | 0.067 | 0.980266655240815 | 0 |
| bla | 861 | golden_gate | 0.084 | 0.9986126172738556 | 0 |
| bla | 861 | biobrick_rfc10 | 0.080 | 1 | 0 |
| waaS | 936 | golden_gate | 0.072 | 1 | 0 |
| waaS | 936 | biobrick_rfc10 | 0.073 | 0.9960797791188929 | 0 |
| yagF | 1968 | golden_gate | 0.167 | 1 | 0 |
| yagF | 1968 | biobrick_rfc10 | 0.166 | 0.9962844172253809 | 0 |

rpoB is 4029 nt, past the length where that DP table claims a full optimum.
The challenging sweep still designs it. Release: Golden Gate 0.302 s, CAI
0.999709780691428; BioBrick 0.303 s, CAI 0.9930619119141477. Both vectors
match the v2.9.1 challenging note. Zero actionable violations.

### Harsh assembly

The stress set described above. Budget `fast`, generations 0. waaS returns
`CF-SRCH-006`. The suite records that as best-effort. The pin stays ok
because this set was never a feasibility claim.

| gene | nt | wall s | opt-1 s | seed s | RSS KiB | result |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| waaS | 936 | 0.059 | 0.084 | | 36620 | `CF-SRCH-006` |
| elastin_vpgvg20 | 306 | 0.601 | 0.836 | 0.543 | 47572 | CAI 0.39695841137981935, feasible |
| yeeJ | 7077 | 4.203 | 5.164 | 3.161 | 54316 | CAI 0.9667370523709706, feasible |

yeeJ harsh hash
`e8afe2becefe5a7e49beca4891ffd398d5a5a95e586bdca172c1522e735bae5e`.
2108 edits. The same hash as `docs/benchmarks/v2.9.4-harsh-rows.json`.

### Multi-objective

Vectors are engine objective vectors, in the order of the objectives column.
Release. The harness column says which suite produced the row. Harsh and
profile agree on the shared vectors (yeeJ structure, rpsL / lacZ Potts).

| case | harness | gen | wall s | opt-1 s | RSS KiB | vector |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| yeeJ CAI + structure, assembly | harsh | 2 | 3.458 | 4.236 | 53976 | `[0.9723247970586831, 2.4166863375261896]` |
| rpoB CAI + structure, assembly | harsh | 2 | 1.204 | 1.490 | 51336 | `[0.9826612532623645, 2.5240523442669134]` |
| lacZ CAI + structure, assembly | profile | 2 | 1.384 | 1.664 | 45976 | `[0.9646814771830261, 2.7308022282793587]` |
| waaS CAI + structure, assembly | profile | 2 | 0.795 | 0.913 | 43656 | `[0.918127802287275, 1.8122103971546997]` |
| elastin CAI + structure, assembly | profile | 2 | 1.348 | 1.563 | 44200 | `[0.46212767144714184, 3.0970951335253227]` |
| rpsL CAI + structure, assembly | profile | 2 | 0.239 | 0.261 | 43356 | `[0.985816320550702, 2.6050213909252715]` |
| yeeJ CAI + burden, assembly | harsh | 0 | 4.061 | 4.809 | 54492 | `[0.9727484991594185, 0.5759101658578069]` |
| waaS CAI + burden, assembly | profile | 0 | 0.299 | 0.400 | 43492 | `[0.8995843357003364, 0.4991567438774198]` |
| rpoB CAI + burden, empty | profile | 0 | 0.636 | 0.727 | 45452 | `[1.0, 0.6171608892818251]` |
| aceB CAI + burden, empty | sweep | 0 | 0.257 | 0.271 | 44216 | `[0.9530633221854043, 0.510773980693139]` |
| waaS CAI + structure + burden | harsh | 2 | 0.859 | 1.019 | 48320 | `[0.8995843357003364, 1.8212043757005512, 0.4991567438774198]` |
| yeeJ CAI + Potts, empty | harsh | 0 | 1.638 | 2.103 | 53292 | `[0.9788608339049295, -0.6223316995686843]` |
| rpoC CAI + Potts, empty | harsh | 0 | 0.815 | 1.006 | 52900 | `[0.9414664807691052, -0.59740611181569]` |
| lacZ CAI + Potts, empty | harsh | 0 | 0.673 | 0.885 | 52856 | `[0.8859755189212043, -0.5677755453484393]` |
| waaS CAI + Potts, empty | profile | 0 | 0.405 | 0.472 | 47200 | `[0.9593146794501609, -0.5143818504117195]` |
| rpsL CAI + Potts, empty | harsh | 0 | 0.342 | 0.392 | 52876 | `[0.9303092945130915, -0.5887977805472148]` |

Potts is `epistasis.context_naturalness`. It is available and it ran. The
second component is the engine’s naturalness term (negative). Empty
constraints, so these rows do not pay assembly seed repair. yeeJ’s 1.638 s
still includes the front. Wrapped Rust stages on that row are tier 2
0.348 s and tier 1 0.078 s.

Fold counts on the structure rows, release profile suite: rpsL 12 hits,
3 full folds, 10 local edits; waaS 13 / 2 / 18; elastin 16 / 4 / 4; lacZ
14 / 1 / 21; yeeJ 12 / 1 / 18. The yeeJ counts are the published pin.

yeeJ CAI plus burden keeps the assembly-only CDS (same hash) and adds the
burden axis. The wall is higher than assembly-only because tier 2 has a
second term to score: tier 2 is 1.295 s against 0.253 s on the CAI-only
design row.

Standard, 4 generations, challenging suite, release. These are the four
slowest design jobs in that sweep.

| case | wall s | RSS KiB | vector |
| --- | ---: | ---: | --- |
| yeeJ CAI + structure, assembly | 4.951 | 50540 | `[0.9719613556538297, 2.414045979917117]` |
| yeeJ CAI, assembly | 3.833 | 50460 | `[0.9727484991594185]` |
| elastin CAI + structure, assembly | 1.876 | 44780 | `[0.46212767144714184, 3.0970951335253227]` |
| lacZ CAI + structure, assembly | 1.857 | 45824 | `[0.9649191059215455, 2.7217215882027954]` |

The yeeJ CAI row matches the fast CDS. The two structure fronts move by a
few thousandths relative to the 2-generation rows. The v2.9.1 standard rerun
picked rpoB assembly instead of elastin structure, because rpoB assembly was
slow then (7.156 s) and is 0.665 s here.

### Operon

Joint growth cost, two passes, fluxes 1 and 2, budget `fast`. Release harsh
suite. Proteins match on every cistron.

| operon | nt | wall s | opt-1 s | RSS KiB | joint | swaps per pass |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| rpsL + waaS | 375 + 936 | 0.382 | 0.350 | 47796 | 0.4475618327684962 | `[22, 0]` |
| lacZ + thrA | 3075 + 2463 | 1.611 | 1.755 | 48548 | 0.4452346719259927 | `[26, 0]` |
| yeeJ + rpoB | 7077 + 4029 | 3.517 | 3.937 | 51080 | 0.44521127470895777 | `[42, 1]` |

The first two joints and swap counts match `v2.9.3-profiles.md` and the
v2.9.1 challenging note. yeeJ + rpoB is the longer case added in the harsh
suite. Its joint matches the published harsh row.

### Phase split on the slow seeds

`CF_SEED_PROFILE` on the release seed rows. Inclusive: `climb` contains
`kmer_replace`, `ascend`, `gc_delta`, and the sweeps. `is_feasible` is the
integer offender check. `checker_check` (the full checker fallback) is 0
calls on both yeeJ seeds.

| function | yeeJ assembly calls | seconds | yeeJ harsh calls | seconds |
| --- | ---: | ---: | ---: | ---: |
| `climb` | 4 | 1.035 | 4 | 2.072 |
| `kmer_replace` | 147250 | 0.335 | 195161 | 0.831 |
| `ascend` | 4 | 0.358 | 4 | 0.543 |
| `fix_kmers` | 11 | 0.098 | 20 | 0.232 |
| `gc_delta` | 270430 | 0.086 | 342977 | 0.174 |
| `polish` | 12 | 0.098 | 12 | 0.148 |
| `gc_sweeps` | 9 | 0.067 | 18 | 0.115 |
| `is_feasible` | 22823 | 0.0011 | 31228 | 0.0014 |
| seed wall | | 1.655 | | 3.161 |

Call counts on the harsh yeeJ seed match the “after” column of
`v2.9.4-harsh-perf.md` (195161 `kmer_replace`, 342977 `gc_delta`, 31228
`is_feasible`). The seconds there were opt-level 1 (`kmer_replace` 0.938 s,
`ascend` 0.724 s, `gc_delta` 0.361 s, climb 2.268 s). Release does the same
calls faster.

cProfile of the release yeeJ assembly design (`profile_suite.py`, cumulative
seconds inside that process): `design` 2.367 s, of which
`resolve_feasible_seed` is 1.599 s, tier 2 polish 0.190 s, and
`structure_diag.annotate` → ViennaRNA `pfl_fold_up` 0.436 s. That
annotation runs because structure is not an objective on a CAI-only job
(`design.py` skips it when the term is already scored). It is about a
sixth of the design call on yeeJ, after the seed. The profiled wrapper also
spends 0.211 s in `_sequences` loading the corpus. That cost is the harness,
measured inside cProfile, and is outside `design()`.

## Primer design

Release harsh suite unless the row says profile. Primer length 18–25 except
the oligo pools, which request 12–36. Tm target 60 °C. The 40-oligo cap
still truncates the ranked set at 1600 survivors. ΔG is
`hetero_dimer_dG` in kcal/mol.

| case | template | windows | wall s | opt-1 s | RSS KiB | pair | product | ΔG | survivors |
| --- | --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: |
| yeeJ flank 300, pin | 7077 | 4472 | 0.052 | 0.067 | 37000 | `TGACCGACAAATATTATG` / `GAACTACAGTGGTTTAAA` | 6843 | −1.283850000000001 | 1600 |
| yfhM flank 300 | 4962 | 4472 | 0.056 | 0.065 | 40320 | `GCCTTCTGAAGTTACTAA` / `ATCTGGTTCAGTAAGTTC` | 4653 | −4.680499999999995 | 1600 |
| rpoC flank 300 | 4224 | 4472 | 0.041 | 0.048 | 38032 | `GGTAAAAGATTACGAGTG` / `GATAACGTTCTCTTTCAG` | 3865 | −6.250615000000003 | 1160 |
| yeeJ oligo 12–36 | 7077 | 13850 | 0.066 | 0.076 | 37508 | `CCAATACGGTGCC` / `CGTGTAATGCCCC` | 6734 | −2.123515000000003 | 1600, truncated |
| yfhM oligo 12–36 | 4962 | 13850 | 0.080 | 0.076 | 38844 | `GTTACGCGTAGCC` / `GCCGCTACGTAAC` | 4573 | −5.792715000000001 | 1600, truncated |

The yeeJ flank-300 document SHA-256 is
`5f1fa115bcefde58d0c3d92a35e7adbc5fcbf04bf78b6efe43bf75f946aa149f`.
That is the pin shared with `test_profile_identity` and with the v2.9.2 /
v2.9.3 profile notes. The API primer job returns the same pair, the same
product, and the same ΔG (HTTP 202, search 0.106 s inside the API row).

The oligo pairs and the yfhM / rpoC flank pairs match the recorded (not
pinned) pairs in `v2.9.4-harsh-perf.md`. The wide pool is the slowest
primer row at 0.080 s.

Phase split on the release yeeJ oligo pool (0.066 s): `tm_prefilter` 12513
calls, 0.024 s; `rank` 0.018 s; `dimer` 1680 calls, 0.011 s. On the pinned
18–25 search, `rank` is 0.019 s and `dimer` is 0.017 s. GC and homopolymer
screens are under 0.001 s.

Profile-suite regions, release. All nine document SHA-256 values match
`v2.9.3-profile-rows.json`. The synthetic `ACGT` flanks are refused because
the dimer is below −9 kcal/mol. They are a repeat, not a biological flank.
Two rows share the name `yeeJ region 1000 flank 300` and differ in geometry;
the product size tells them apart.

| case | wall s | v2.9.3 wall s | product | ΔG | refused | SHA-256 prefix |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| yeeJ region 1000, flank 50 | 0.002 | 0.003 | 1072 | −3.507365 | no | `e57132f535a1` |
| yeeJ region 1000, flank 100 | 0.043 | 0.055 | 1180 | −2.1353599999999986 | no | `8c8d3c988daf` |
| yeeJ region 1000, flank 200 | 0.044 | 0.057 | 1249 | −2.753985 | no | `cc92049d4ea1` |
| yeeJ region 1000, flank 300 (mid window) | 0.050 | 0.077 | 1249 | −2.753985 | no | `e40da98e67dc` |
| yeeJ region 300, flank 300 | 0.057 | 0.071 | 822 | −1.8601799999999997 | no | `04f69cb6ad99` |
| yeeJ region 1000, flank 300 (from nt 300) | 0.051 | 0.076 | 1394 | −1.7191900000000047 | no | `250193a8e2fa` |
| yeeJ region 6477, flank 300 | 0.050 | 0.064 | 6843 | −1.283850000000001 | no | `5f1fa115` |
| rpsL whole CDS, flank 40 | 0.001 | 0.002 | 371 | −3.9297299999999993 | no | `ccc17c9e987e` |
| yeeJ synthetic `ACGT` flanks | 0.039 | 0.054 | 7395 | −20.754845000000046 | yes | `0754bf08ea42` |

Challenging sweep, release, 300 nt CDS flanks, lengths 18–25: yeeJ 0.059 s,
product 6843, ΔG −1.283850000000001, not refused; rpoB 0.060 s, product
3908, ΔG −3.387435, not refused; yeeJ with synthetic flanks 0.051 s,
product 7395, ΔG −20.754845000000046, refused. The rpoB pair’s ΔG matches
the v2.9.1 challenging note.

## API and GATE-BEN-001

One release harsh row posts three designs at once (rpsL empty CAI, waaS
assembly CAI, rpsL CAI + burden), then one yeeJ primer search, then the
refusals the suite requires.

| check | result |
| --- | --- |
| three design accepts | HTTP 202, accept wall 0.006 s |
| three designs finished | all `succeeded`; vectors `[1.0]`, `[0.9259163095438141]`, `[1.0, 0.6171608892818251]` |
| overlap | accepts overlap; `Api.run_lock` runs the designs one at a time |
| row wall | 0.958 s |
| peak RSS | 63112 KiB, the highest RSS in the release harsh suite |
| primer | the pinned yeeJ pair, product 6843, ΔG −1.283850000000001 |
| `design(..., background=True)` | `CF-CLI-003` |
| `GET /v1/hosts`, `POST /v1/feasibility` | HTTP 501, `CF-CLI-003` |
| `epistasis.dwell_deviation` without a fitted head | `CF-LAND-001` |

`GATE-BEN-001` is the stored CAI-only claim, then a fresh `run_cai_only`.
Release: claimed 50 of 50, `codonforge_ge_greedy` 50, `behind_greedy` empty,
zero genes below the stored DNA Design Suite CAI. Wall 5.985 s, RSS 50332 KiB.
Opt-level 1: 6.314 s, same claim. This is the CAI comparison against
`cai_greedy` on 50 *E. coli* K-12 proteins. It is not an MFE comparison and
not LinearDesign or CDSfold.

## Comparison to prior baselines

Outputs were compared field by field. Wall times from older notes were
measured on machines that report the same CPU string (`Intel(R) Xeon(R)
Processor`) and the same ViennaRNA 2.7.2. They are still different runs.
A few tenths of a second on a long row is what this machine showed between
two launches of the same job.

### Harsh v2.9.4, same commit, same outputs

`docs/benchmarks/v2.9.4-harsh-rows.json` was produced by `harsh_perf.py` on
this commit at opt-level 1. This release run and this opt-level 1 run match
that file on state, CAI, objective vector, CDS SHA-256, primer sequences,
product size, ΔG, document SHA-256, operon joint cost, swap counts,
feasibility, problem code, and the gate claim. The diff count on those
fields is 0 across 76 rows. What moved is wall time.

| case | published harsh s (opt-1) | this opt-1 s | this release s |
| --- | ---: | ---: | ---: |
| yeeJ empty CAI | 0.533 | 0.563 | 0.526 |
| yeeJ assembly | 2.785 | 2.825 | 2.727 |
| yeeJ assembly standard ×4 | 4.452 | 4.749 | 4.012 |
| yeeJ harsh assembly | 4.845 | 5.164 | 4.203 |
| yeeJ CAI + structure | 3.956 | 4.236 | 3.458 |
| yeeJ CAI + burden | 4.484 | 4.809 | 4.061 |
| yeeJ CAI + Potts | 1.796 | 2.103 | 1.638 |
| yfhM assembly | 1.536 | 1.570 | 1.519 |
| yfhM assembly standard ×4 | 3.543 | 3.716 | 3.081 |
| rpoC assembly | 1.401 | 1.469 | 1.197 |
| rpoB CAI + structure | 1.385 | 1.490 | 1.204 |
| operon yeeJ + rpoB | 3.706 | 3.937 | 3.517 |
| GATE-BEN-001 | 5.899 | 6.314 | 5.985 |
| yeeJ primer flank 300 | 0.061 | 0.067 | 0.052 |
| yeeJ oligo 12–36 | 0.072 | 0.076 | 0.066 |

Release is faster than both opt-level 1 runs on the long rows, by roughly
10–20% against the published note and more against this machine’s opt-level
1 pass (yeeJ harsh 5.164 s → 4.203 s). The call counts did not change. This
is the effect of `opt-level` 3 and thin LTO on the same seed repair, not a
new search.

### v2.9.4 search-escape CAI

`docs/benchmarks/v2.9.4-search-escape.md` is the quality floor. This run
meets every assembly CAI in that table at equality (the harsh pin
`assembly_floor`). waaS is 0.9259163095438141. yeeJ is 0.9727484991594185.
rpsL is 0.9866866273575934. Elastin’s engine vector is
0.4642192192077335. Fast and standard return the same CDS on yeeJ, yfhM,
and the floor genes that were checked at both budgets, which is what that
note reported: the gain is in the seed, and tier 3 does not replace it.

### v2.9.3 and v2.9.2 profiles

Those notes used opt-level 1, before the v2.9.4 seed escape
(`DEP-SRCH-021`, `DEP-SRCH-022`) and the rpoC pair escape (`DEP-SRCH-023`).
Assembly CAI is higher now, and the seed does more work, so a higher wall
on an assembly row is a different search. Empty-constraint CAI, Potts
vectors on rpsL / waaS / lacZ, the rpsL structure vector, the two shorter
operons, and all nine primer document hashes are the same numbers as
`v2.9.3-profiles.md` (and the primer hashes match v2.9.2, which v2.9.3
already rechecked).

Profile-suite walls, assembly CAI, generations 0. v2.9.3 is the published
opt-level 1 number. “Now” is this release profile run.

| case | v2.9.2 s | v2.9.3 s | v2.9.3 vector | now s | now vector |
| --- | ---: | ---: | --- | ---: | --- |
| yeeJ assembly | 3.048 | 1.175 | `[0.9580065306093907]` | 2.532 | `[0.9727484991594185]` |
| rpoB assembly | 2.057 | 0.589 | `[0.9798209521813181]` | 0.644 | `[0.9830194072110617]` |
| lacZ assembly | 0.646 | 0.418 | `[0.9516318473376221]` | 0.617 | `[0.9665531384697095]` |
| yagF assembly | 0.433 | 0.317 | `[0.9286425884184984]` | 0.566 | `[0.9453398274862173]` |
| waaS assembly | 0.155 | 0.127 | `[0.9173459052341191]` | 0.228 | `[0.9259163095438141]` |
| elastin assembly | 0.099 | 0.113 | `[0.43802555520567116]` | 0.391 | `[0.4642192192077335]` |
| rpsL assembly | 0.060 | 0.057 | `[0.9866866273575934]` | 0.067 | `[0.9866866273575934]` |
| rpoB CAI + burden, empty | 0.701 | 0.688 | `[1.0, 0.6171608892818251]` | 0.636 | same vector |
| waaS CAI + burden, assembly | 0.257 | 0.225 | `[0.9173459052341191, 0.5636056357592945]` | 0.299 | `[0.8995843357003364, 0.4991567438774198]` |

rpsL did not move. The other assembly CAIs rose, which is the search-escape
result, and the seed got slower where the extra repairs do real work
(elastin 0.113 s → 0.391 s on this profile run). waaS burden moved with the
new feasible CDS: both axes changed. rpoB CAI + burden with empty
constraints is the same vector and a similar wall.

Structure vectors that changed with the seed, from the v2.9.3 profile note
to this release profile run:

| gene | v2.9.3 | now |
| --- | --- | --- |
| rpsL | `[0.985816320550702, 2.6050213909252715]` | same |
| waaS | `[0.9096294894349499, 1.7845726655476362]` | `[0.918127802287275, 1.8122103971546997]` |
| elastin | `[0.4236181040759199, 3.3760051578497814]` | `[0.46212767144714184, 3.0970951335253227]` |
| lacZ | `[0.9497890800680886, 2.7181800037278987]` | `[0.9646814771830261, 2.7308022282793587]` |
| yeeJ | `[0.9575892496987667, 2.408954059933615]` | `[0.9723247970586831, 2.4166863375261896]` |

Potts, unchanged: rpsL `[0.9303092945130915, -0.5887977805472148]`, waaS
`[0.9593146794501609, -0.5143818504117195]`, lacZ
`[0.8859755189212043, -0.5677755453484393]`. Tier 2 on lacZ was 3.401 s in
v2.9.2, 0.218 s in v2.9.3, and 0.153 s in this release profile row.

### v2.9.1 challenging sweep

That note is opt-level 1, before the v2.9.2 seed rewrite. yeeJ assembly CAI
was 21.202 s at CAI 0.9580065306093907. The same sweep row here, release,
is 2.371 s at CAI 0.9727484991594185. rpoB assembly was 5.275 s at
0.9798209521813181 and is 0.641 s at 0.9830194072110617. yeeJ CAI + structure
was 24.124 s at `[0.9575892496987667, 2.408954059933615]` and is 3.392 s in
this release sweep at the v2.9.4 structure pin. Empty-constraint CAI was
already 1 and already sub-second; yeeJ empty was 0.549 s and is 0.576 s.
Golden Gate and BioBrick on rpoB are the same CAI as that note.

Those drops are the cumulative search work from v2.9.2 through v2.9.4,
measured again here. They are not a property of the release switch alone.
The opt-level 1 column of the harsh table is the right like-for-like
comparison for this commit.

## Bottlenecks, ranked

Ranked from the release counters on this machine. Inclusive times overlap,
so the rank is by the outermost cost that a change would have to move.

1. **Seed repair on long, constrained CDS.** yeeJ harsh seed 3.161 s of a
   4.203 s design. `climb` is 2.072 s. Inside it, `kmer_replace` is 0.831 s
   (195161 calls) and `ascend` is 0.543 s. `gc_delta` is 0.174 s.
   `is_feasible` is 0.0014 s over 31228 calls, which is the v2.9.4
   predicate change holding. The same order appears on ordinary yeeJ
   assembly (climb 1.035 s, ascend 0.358 s, `kmer_replace` 0.335 s) and on
   yfhM and rpoC, where `kmer_replace` and `ascend` are the two largest
   leaves. This is the same ranking the harsh note left in place.

2. **Tier 3 repair when generations are turned on.** yeeJ assembly standard
   adds 1.485 s of tier 3 (repair 0.556 s) and keeps the generations-0 CDS.
   yfhM standard spends 1.722 s in tier 3 (repair 0.758 s) for the same
   outcome. The feasibility gate and the breed check are a few hundredths.

3. **Structure L0 folding, once structure is an objective.** On yeeJ CAI +
   structure, tier 3 is 1.464 s beside a 1.646 s seed. `l0_batch` is 0.637 s
   and repair is 0.753 s. Fold traffic is 1 full fold, 18 local edits, 12
   cache hits. rpoB structure is dominated by `l0_batch` (0.405 s) plus
   repair (0.344 s).

4. **ViennaRNA accessibility annotation on CAI-only fronts.** 0.436 s of
   `pfl_fold_up` on the yeeJ assembly cProfile, after the search has
   finished. It is how the front reports a structure detail for a job that
   did not optimize structure. It does not run when structure is already
   an objective.

5. **`GATE-BEN-001` as a batch.** 5.985 s for 50 sequential CAI designs.
   That is the longest single row. It is fifty short jobs, not one slow
   gene. yeeJ alone at empty constraints is 0.526 s.

6. **Potts tier 2 is no longer the lacZ spike from v2.9.2.** lacZ tier 2 is
   0.153 s. yeeJ CAI + Potts still costs 1.638 s end to end, with only
   0.43 s in the wrapped tier 1 and tier 2 calls. The rest sits in the
   front, including the annotation in item 4.

7. **Primer design is small next to seed repair.** The slowest primer row is
   0.080 s (yfhM, lengths 12–36, 13850 windows). Time goes to `tm_prefilter`,
   `rank`, and `dimer`, in that order on the wide pool. The 40-oligo cap
   still stops the survivor list at 1600. Raising it would change which
   pairs are scored. It was not raised.

8. **API serialization.** Three accepts return in 0.006 s. The designs then
   run one at a time under `Api.run_lock`, and the row takes 0.958 s. Peak
   RSS is 62 MiB because one process holds the server and the jobs.

## Correctness and pin status

| check | status |
| --- | --- |
| `harsh_perf.py --suite all`, release | pass, `pin_failures` empty, exit 0 |
| `harsh_perf.py --suite all`, opt-level 1 | pass, `pin_failures` empty, exit 0 |
| release vs opt-level 1 on CAI, vectors, CDS hashes, primer pairs, ΔG, operon joints, gate | identical |
| release vs published `v2.9.4-harsh-rows.json` on those fields | identical, 76 rows |
| assembly CAI floors from `v2.9.4-search-escape.md` | pass, at equality |
| Golden Gate / BioBrick DP gaps, genes ≤ 2500 nt with a full optimum | pass, gap 0, 16 rows |
| yeeJ structure vector and fold counts | pass |
| yeeJ primer 300 nt flanks | pass, pair, product 6843, ΔG, document SHA-256 |
| nine profile-suite primer documents vs v2.9.3 | pass, SHA-256 multiset equal |
| `GATE-BEN-001` | pass, CAI-only, 50 of 50, 0 regressions vs the stored scores |
| protein identity and actionable violations on succeeded designs | pass |
| refusals (`CF-CLI-003`, HTTP 501, `CF-LAND-001`) | pass |
| waaS harsh assembly | `CF-SRCH-006` in 0.059 s, as recorded, pin ok |
| rpoC assembly | feasible, CAI 0.9593415607455033, the post-#33 CDS |

Four rows differ by one unit in the last place between the independent CAI
checker and `usage.cai`: elastin assembly (`0.46421921920773357` vs
`0.4642192192077335`), elastin harsh, rpoC Potts, and lacZ Potts. The same
pairs are in the published harsh JSON. The pins use a `1e-12` tolerance on
floors and exact equality on the structure vector, and they pass. No search
result was changed for this report.

No optimizer change, no new constraint, and no new sequence went into this
measurement. A correctness bug that would have justified an optimizer edit
did not show up.

## Honest limits

This run measured wall time, the harness phase counters, peak RSS, feasibility,
CAI and objective vectors, and primer pair / ΔG pins, on one VM, one shot
per row.

- Rows are single launches. The yeeJ assembly fast wall on release spanned
  2.371 s to 2.727 s across harnesses. There is no median and no instruction
  or cache-miss profile.
- The published harsh, v2.9.3, and v2.9.2 numbers come from earlier runs
  that report the same CPU string. Cross-note wall ratios of a few percent
  are not a controlled A/B on one binary.
- Release versus opt-level 1 is one rebuild of each profile, run back to
  back on this VM, with the machine otherwise quiet enough to finish, and
  still a shared KVM guest.
- `GATE-BEN-001` stays the CAI-only gate. LinearDesign, CDSfold, and DNA
  Chisel were not rerun. The v2.9.3 external comparison already records
  that boundary.
- The harsh constraint set is a stress set. waaS failing it is a measured
  `CF-SRCH-006`, not a new feasibility theorem and not a defect against the
  assembly file.
- Checkpoint resume, wheel cold-start, and the browser UI were not timed.
  The API row is in-process `serve` traffic for three jobs and one primer
  search.
- Seed 7 is the published pin seed. The repair does not read the job seed;
  that was established in the search-escape note and was not re-swept here.
- Genes are the shipped panel and validation FASTA set. Nothing was added
  to make a longer CDS.
- The 40-oligo primer cap stayed at 40. A larger pool would change the
  ranked pairs and the ΔG that gets reported.
- ViennaRNA was present (2.7.2). Structure rows need it. Without it those
  rows fail or skip, and this report would not have the structure vectors.
- Package version stays 0.1.0. This document does not tag a release.
