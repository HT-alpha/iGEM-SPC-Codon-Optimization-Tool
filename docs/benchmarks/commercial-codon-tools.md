# Codon optimization: DNA Design Suite and other tools

Package version 0.1.0. Optimizer commit
`db867d4aa31a76229098393dbd0448b01a341032`. Nothing here is tagged.

This note compares DNA Design Suite with commercial codon-optimization
services and with the open tools a third party can actually run. The only
head-to-head numbers are DNA Design Suite against DNA Chisel 3.2.16, on one
machine, with one codon table and one constraint language. Commercial
services were not executed. Their sections are a capability matrix from
published methods and from public descriptions, plus a record of which
sites this environment could open on 18 October 2026.

## Where the earlier comparison lives

Pull request 31 measured the same problem and wrote it down. The files on
`main` are:

- `docs/benchmarks/v2.9.3-external-comparison.md`
- `docs/benchmarks/v2.9.3-optimization-quality.md`
- `docs/benchmarks/v2.9.3-external-rows.json`

That measurement used commit `2a9981a69057486d71b66f8830f64700dcf5b7ec`.
DNA Design Suite's CAI was higher than DNA Chisel's on every assembly gene
where DNA Chisel returned a sequence, and DNA Chisel was never higher on
Golden Gate or BioBrick. `cai_greedy` at CAI 1 was infeasible under the
assembly set. LinearDesign and CDSfold were recorded there as a different
objective: folding free energy, not `usage.cai`, and not given the
constraint files. This note does not replace those tables. It re-runs the
DNA Chisel arm on the current tree for four of those genes, adds four
genes from the 255-gene corpus, and states what cannot be said about the
commercial tools.

## What was run

Host `escherichia_coli:511145/MG1655`, translation table 11, weights
`examples/mg1655_ribosomal_usage.json`. Every family maximum in that table
is 1. The smallest positive weight is 0.002392. The reported CAI is
`benchmarks/validation/cai_independent.py`: Sharp and Li, geometric mean of
relative adaptiveness, Met, Trp and stops left out. The design engine's
`usage.cai` was scored on the same sequences. Across the 99 scored rows the
largest absolute difference is 5.551115123125783e-16.

DNA Design Suite jobs are budget `standard`, 4 generations, stagnation
window 21, seed 7, objective `usage.cai`. That is the job in the pull
request 31 file. It was not retuned per gene.

DNA Chisel 3.2.16 receives the same CDS, the same ribosomal weights through
`CodonOptimize(method="use_best_codon")`, and the same three constraint
sets, mapped by `tools/external_comparison.py`. `numpy` and `random` are
seeded. One seed is seed 7, except the rpsL assembly probe at seeds 1, 7
and 11. `MaximizeCAI` scores `N * log(CAI)` and replaces a zero frequency
with 0.001 (`dnachisel/.../MaximizeCAI.py`). The smallest positive
ribosomal weight is above 0.001, so that floor does not reorder the
positive codons. DNA Chisel still adds Met, Trp and the stop into its own
objective. The CAI column does not.

Constraint sets:

- assembly: `benchmarks/validation/constraints_assembly.json`. GC in every
  50 nt window between 0.40 and 0.60, BsaI, BsmBI, EcoRI, XbaI, SpeI, PstI
  and NotI on both strands, no homopolymer longer than 6, `AGGAGG` forbidden
  on both strands, no repeated 12-mer on the forward strand.
- Golden Gate: `examples/golden_gate.json`. BsaI and BsmBI, both strands.
- BioBrick RFC[10]: `examples/biobrick_rfc10.json`. EcoRI, XbaI, SpeI, PstI,
  NotI, both strands.

The constraint mapping is the one documented in the pull request 31 note:
`EnforceGCContent` steps by 1, `AvoidPattern` rejects the reverse
complement, a homopolymer of 7 or more is `HomopolymerPattern`, and
`UniquifyAllKmers(k=12, include_reverse_complement=False)` matches the
forward-strand 12-mer rule. Three 54 nt inputs that must fail were
evaluated and not repaired. All three fail both checkers
(`negative_cases_fail_both` is true in the JSON).

Named sequences, the ones the task asked to reuse:

| gene | file | nt |
| --- | --- | ---: |
| rpsL | `benchmarks/validation/sequences/rpsL.fa` | 375 |
| waaS | `benchmarks/validation/sequences/waaS.fa` | 936 |
| lacZ | `benchmarks/validation/sequences/lacZ.fa` | 3075 |
| yeeJ | `benchmarks/validation/sequences/yeeJ.fa` | 7077 |

Corpus sample, four genes, not chosen by score. From
`benchmarks/corpus/generality_sample.json` (digest
`f1cf5a7a172bdf4974d942fe7a2e00aaf069eb5f60f627d92b1028b211b214dd`), each
length bin contributes the first record with at least 150 nt whose name is
not one of the four above, ordered by `(nt, gene)`.

| gene | bin | protein_id | nt |
| --- | --- | --- | ---: |
| hokC | short | YP_025292.1 | 153 |
| gfcA | medium | NP_415507.2 | 306 |
| dppC | long | NP_417999.1 | 903 |
| rpoD | very long | NP_417539.1 | 1842 |

The other 251 corpus genes were not sent to DNA Chisel. The 255-gene DNA
Design Suite distributions are already in `docs/benchmarks/v2.10-generality.md`.

Build and machine. Release extension, the same command the earlier note
used:

```sh
.venv/bin/maturin develop --release
.venv/bin/pip install 'dnachisel==3.2.16' 'ViennaRNA==2.7.2'
.venv/bin/python benchmarks/commercial_comparison/run_chisel_vs_dds.py \
    --out docs/benchmarks/commercial-comparison-rows.json
```

The process that wrote the JSON exited 0. Machine: 4 CPUs, Intel Xeon,
16387280 KiB RAM, Linux 6.12.94+, Python 3.12.3, Biopython 1.88, numpy
2.5.4, ViennaRNA 2.7.2. DNA Design Suite package and core both report
0.1.0.

DNA Design Suite wall time is `time.perf_counter` around `design()` in the
already loaded process. DNA Chisel wall time is `time.perf_counter` around
a fresh interpreter, because each solve runs in a child process so a stall
can be killed. The earlier note measured about half a second of that
startup on a Golden Gate solve both tools finish at CAI 1. This run does
not subtract a startup constant. The short DNA Chisel walls sit near 0.55 s,
which is the same pattern. On yeeJ assembly the recorded walls are
1.949 s and 5.614 s, and that gap is larger than half a second.

`accessibility_penalty` is
`dna_design_suite.codon.orchestrator.structure_diag.evaluate` after the
design: mean −log unpaired probability, window 6 nt, span 40 nt, 37.0 °C.
Neither optimizer was asked to optimize it. It is not global minimum free
energy.

## Fairness

A commercial codon optimizer and DNA Design Suite are the same product
category only at the level of "synonymous CDS in, CDS out". They are not
the same optimization problem unless the codon table, the objective, and
the constraints are the same. Those three were aligned for DNA Chisel.
They were not aligned for any commercial tool, because this environment
could not submit a sequence and the vendors do not publish the table and
the weights.

Consequences, each of which blocks a numeric winner:

- A CAI computed on `examples/mg1655_ribosomal_usage.json` is not the CAI
  a vendor prints from a Kazusa genome table or from a proprietary
  expression-weighted table. Rescoring would require the vendor's sequence.
- Expression fold-changes in cell papers are a different measurement from
  CAI. This note does not convert one into the other.
- An unpublished weighted sum of CAI, GC, hairpins and synthesis complexity
  can lose on our CAI and still be the sequence the vendor will synthesize.
  That is not a loss on our objective, and it is not a win on theirs.
- LinearDesign and CDSfold remain outside this comparison. The pull request
  31 file says why: their objective is folding free energy plus a codon
  term, they were not given the constraint files, and on that 15 GiB
  machine lacZ was killed around 15 GB RSS. They were not rerun.
- DNA Chisel can also match a codon distribution (`match_codon_usage`) or
  harmonize relative adaptiveness (`harmonize_rca`). This run asked it for
  `use_best_codon`, which is the CAI objective, because that is the
  objective DNA Design Suite was given. A harmonization contest was not run.
- The corpus arm is four genes. It does not extend a claim to all 255.

`cai_greedy` picks the highest weight at each codon and breaks ties toward
the lexicographically smaller codon. It never sees a constraint. Under
assembly its CAI of 1 is a feasible design only for hokC. On the other
seven genes it is the unconstrained ceiling, and it fails the checker.

## Sites this environment could open

Egress from the measurement VM is an allowlist. On 18 October 2026:

| URL | Result |
| --- | --- |
| `https://github.com/Edinburgh-Genome-Foundry/DnaChisel` README | fetched, HTTP 200 |
| `https://pypi.org/pypi/dnachisel/json` | version 3.2.16 |
| `https://www.idtdna.com/pages/tools/codon-optimization-tool` | connection failed |
| `https://www.genscript.com/tools/gensmart-codon-optimization` | connection failed |
| Thermo Fisher GeneArt gene-synthesis host | connection failed |
| `https://www.twistbioscience.com/` | connection failed |
| `http://www.jcat.de/` | HTTP 403, body `Forbidden` (9 bytes) |
| `http://genomes.urv.es/OPTIMIZER/` | HTTP 403, body `Forbidden` |
| `http://cool.syncti.org/` | DNS did not resolve |
| `http://bioinfo.bti.a-star.edu.sg/COOL/` | HTTP 403, body `Forbidden` |

The 403 responses are the body this VM received. They are not a measurement
that the public server is down. No sequence was pasted into a commercial or
academic web form. No number below is a scraped vendor score.

## Capability matrix

"Run here" means this study executed the tool on the shared CDS. "Offline"
means a person can run it without an account after installing it.

| Tool | Class | How a user runs it | Host tables | Objectives the publisher describes | Constraints the publisher describes | Batch or API | Offline | Published method | Run here |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DNA Design Suite 0.1.0 | local library, Apache-2.0 | CLI, Python, local page | caller-supplied JSON; example files on disk; no downloaded catalogue (`CF-CLI-003` on seed import) | `usage.cai`, violation count, `structure.local_accessibility`, `burden.growth_cost`, `epistasis.context_naturalness` when the E. coli Potts file is installed. `usage.gc_content` is diagnostic only (`CF-LAND-014`) | motifs, restriction sites, homopolymers, k-mer repeats, GC windows, local structure, codon bans, minimum frequency, regex. Schema `schemas/v1/constraint_set.schema.json` | CLI, Python, local HTTP | yes; ViennaRNA is installed by the user | the specifications in `planning/` | yes |
| DNA Chisel 3.2.16 | local library, MIT, Edinburgh Genome Foundry | Python, CLI, and a web app named in the README | `python_codon_tables` names (`e_coli`, `h_sapiens`, and the others listed in `CodonOptimize`), a caller dict, or a TaxID download | `use_best_codon` (CAI), `match_codon_usage`, `harmonize_rca` | the built-in specifications in the README, including GC windows, patterns, homopolymers, k-mers, translation, hairpins | Python and CLI | yes | Zulkower and Rosser, Bioinformatics 2020, 36(16):4508–4509. README: `https://github.com/Edinburgh-Genome-Foundry/DnaChisel` | yes |
| Thermo Fisher GeneArt / GeneOptimizer | commercial, attached to gene synthesis | vendor account and order flow | vendor hosts | multiparameter codon adaptation | the 2010 paper describes GC content, repeats, restriction sites and other motifs inside one sliding-window search | the order flow; no local API was available here | no | Raab, Graf, Notka, Schödl and Wagner, Syst Synth Biol 2010, 4(3):215–225, doi:10.1007/s11693-010-9062-3 | no |
| IDT Codon Optimization Tool | commercial web tool | `https://www.idtdna.com/pages/tools/codon-optimization-tool` | an organism list on the form; the live form was not opened | not published as an equation | the public tool is known as an organism choice plus synthesis-oriented filters such as restriction sites; the current checkbox list was not re-read | web form; no API was called | no | no methods paper that can be replayed on our table | no |
| GenScript GenSmart / OptimumGene | commercial; GenSmart is the free web optimizer, OptimumGene is the name used on synthesis orders | GenScript web tool and order flow | vendor hosts | vendor descriptions list codon adaptation together with GC, repeats and mRNA structure; the weights are not published | vendor descriptions list restriction sites and sequence motifs; not a file this study can diff | web form and the order flow; no API was called | no | no replayable methods paper in this study | no |
| Twist Bioscience | commercial, order flow | Twist gene-order portal | not verified from this VM | not verified from this VM | not verified from this VM | the order flow | no | no paper cited here | no |
| JCat | academic web server | `http://www.jcat.de/` | the server's host list | codon adaptation to the selected host; the paper reports CAI | rho-independent terminators and restriction sites, as described in the paper | one sequence per form submission; this VM received HTTP 403 | no | Grote et al., Nucleic Acids Res 2005, 33:W526–W531, doi:10.1093/nar/gki376 | no |
| OPTIMIZER | academic web server | `http://genomes.urv.es/OPTIMIZER/` | precomputed genome and highly expressed tables | most frequent codon, random draw from the frequencies, or a Monte Carlo search that raises CAI | as described in the paper; not re-tested | web form; this VM received HTTP 403 | no | Puigbò, Guzmán, Romeu and Garcia-Vallvé, Nucleic Acids Res 2007, 35:W126–W131, doi:10.1093/nar/gkm219 | no |
| COOL | academic web server | the 2014 paper's server | the server's tables | a genetic algorithm over codon usage, codon context, hidden stops, GC and repeats | those same terms, as objectives rather than as our JSON | `cool.syncti.org` did not resolve; the A*STAR URL returned HTTP 403 | no | Chin, Chung and Lee, Bioinformatics 2014, 30(15):2210–2212, doi:10.1093/bioinformatics/btu192 | no |

DNA Chisel's `AvoidHairpins` docstring says the hairpin test follows IDT
guidelines, with a default stem of 20 and a window of 200. That is DNA
Chisel's implementation of a synthesis guideline. It is not a measurement
of the IDT codon optimizer, and this run did not enable `AvoidHairpins`.

GeneOptimizer's 2010 paper is a real algorithm description: a sliding window
over synonymous codons, used because the full codon combination space is too
large, trading codon adaptation against sequence parameters such as GC and
unwanted motifs. The paper does not ship a binary, a codon table file, or a
constraint file that this repository can replay. Later GeneArt expression
papers measure protein output in cells. Those numbers are not placed next to
our CAI.

IDT, GenScript and Twist cells that would require the live help page are
either left as "not verified from this VM" or limited to properties that
have been stable in the public descriptions of those products: they are
account or order-flow services, they are not local Apache-licensed
libraries, and they do not accept this repository's constraint JSON. The
GenScript row's objective list is that kind of public description. It was
not re-fetched on 18 October 2026, and it is not a measured feature test.

JCat, OPTIMIZER and COOL are academic servers, not commercial products.
They are in the table because a reader will ask. They were not runnable
from this VM.

## Performance

Full grid, including `cai_greedy`, the native CDS, GC, the feasibility
string, wall time and the post-hoc accessibility score:
`docs/benchmarks/commercial-comparison-tables.md`. Raw rows, including the
sequences:
`docs/benchmarks/commercial-comparison-rows.json`.

Twenty-four optimizer pairs (8 genes × 3 constraint sets). Both tools
returned a sequence on every pair. Every one of those 48 sequences passes
our checker and DNA Chisel's checker. Every one translates to the input
protein, including the stop. Actionable GC-window counts on those sequences
are 0. Actionable motif, site, homopolymer and 12-mer counts are 0.

DNA Design Suite's CAI is higher on 14 pairs and equal on 10. DNA Chisel's
CAI is higher on none.

### Assembly

This is the set where the tools can disagree, because the GC window and the
scattered 12-mer are not a small regular pattern.

| gene | nt | DNA Design Suite CAI | DNA Chisel CAI | both feasible | DNA Design Suite wall s | DNA Chisel wall s |
| --- | ---: | ---: | ---: | --- | ---: | ---: |
| rpsL | 375 | 0.9866866273575934 | 0.9641706286088719 | yes | 0.09200841199998422 | 0.5724054179999598 |
| waaS | 936 | 0.9259163095438141 | 0.874865013241299 | yes | 0.6451901809999754 | 0.7489718560000256 |
| lacZ | 3075 | 0.9665531384697095 | 0.8852227588514382 | yes | 0.7530243009999822 | 1.5047314999999912 |
| yeeJ | 7077 | 0.9727484991594185 | 0.9011770220723586 | yes | 1.9491075590000264 | 5.613802587000009 |
| hokC | 153 | 1.0 | 1.0 | yes | 0.01675851700002795 | 0.5765092009999648 |
| gfcA | 306 | 0.9910026611127586 | 0.9853950011557764 | yes | 0.04517637499998273 | 0.584888478000039 |
| dppC | 903 | 0.9420137707572441 | 0.8746543572307005 | yes | 0.24883572699997103 | 0.6919853640000042 |
| rpoD | 1842 | 0.979989656136133 | 0.9433657287975374 | yes | 0.5275929240000323 | 0.8509378439999864 |

`cai_greedy` at CAI 1 passes assembly only for hokC. On rpsL it fails one
12-mer and one PstI site. On waaS it fails 279 GC windows. On yeeJ it fails
287 twelve-mers. The native yeeJ CAI on this ribosomal table is 0, because
a counted codon has weight 0. The optimized yeeJ sequences are not 0.

### Golden Gate and BioBrick

On the regular site sets the two tools often return the same CAI. Where
the CAI differs, DNA Design Suite is higher, and both sequences still pass
both checkers. The pairs that tie include every Golden Gate gene except
lacZ, yeeJ and rpoD, and the BioBrick genes rpsL, waaS, hokC and gfcA.
Golden Gate lacZ is 0.9984294131415632 versus 0.9976272559213156. BioBrick
lacZ is 0.9928685731868527 versus 0.9874846513175503. The other differing
pairs are in the tables file. Several of the ties are CAI 1, which is the
unconstrained ceiling, so the restriction sites did not bind.

### What changed since pull request 31

For rpsL, waaS, lacZ and yeeJ, each tool, and each of the three constraint
sets, the new sequence was compared with `v2.9.3-external-rows.json` by
SHA-256. Twenty-one of those twenty-four sequences match. DNA Chisel matches
on all twelve. The rpsL assembly seed probe matches the three published
hashes (`0c84699b…`, `5cd715fa…`, `5e64b006…`) at CAI 0.9602135976257714,
0.9641706286088719 and 0.95775539455032. DNA Design Suite's rpsL assembly
sequence is the published one, CAI 0.9866866273575934, above all three.

The three new sequences are DNA Design Suite under assembly:

| gene | pull request 31 CAI | this run CAI | still feasible |
| --- | ---: | ---: | --- |
| waaS | 0.9173459052341192 | 0.9259163095438141 | yes |
| lacZ | 0.9516318473376221 | 0.9665531384697095 | yes |
| yeeJ | 0.9580065306093907 | 0.9727484991594185 | yes |

DNA Chisel's CAI on those three genes is unchanged, so the gap widened
because DNA Design Suite's search result moved, not because DNA Chisel was
given a weaker specification. The commits between `2a9981a` and `db867d4`
include the seed-repair work. This run did not bisect which commit moved
the CAI. Walls are also lower than the numbers in the earlier file. Those
walls are a different binary on the same machine class (4 CPU Xeon, about
15 GiB). They are not a paired timing study, and DNA Chisel's wall still
includes process startup.

The post-hoc accessibility numbers do not follow the CAI order. On assembly
lacZ, DNA Chisel's penalty is 2.650197418888983 and DNA Design Suite's is
2.7314410684084045. On assembly waaS the order is the other way. Nobody
asked either tool to optimize that score.

## What the search still does not settle

The quality note's pattern DP ignores the GC window and the 12-mer. On
rpsL that bound is 0.9904970702046313 and the feasible design is still
0.9866866273575934, because that sequence did not change. The bound
sequence itself fails a 12-mer, so the true assembly optimum sits somewhere
in between. On waaS the published bound is 0.9903106518103836 against a new
achieved CAI of 0.9259163095438141, and the bound sequence was infeasible
(273 GC windows and 4 twelve-mers in that note). The difference is not the
optimality gap. lacZ and yeeJ are longer than the 2500 nt limit that note
set for the pattern DP. No new bound was computed here.

Seven of the eight DNA Design Suite assembly rows carry the warning
`CF-SRCH-017` on `seed_repair`: tier 1 does not carry the GC window or the
scattered 12-mer, and the feasible sequence is what seed repair and four
generations produce from there. hokC does not carry that warning. Its
unconstrained codon choice is already assembly-feasible, and both tools
return CAI 1. The 17 October 2026 performance report
records waaS on the harsher assembly set as `CF-SRCH-006` in 0.059 s. That
set was not repeated against DNA Chisel.

## Advantages of DNA Design Suite

- On the 24 pairs where the codon table, the protein and the constraints
  match, the CAI is higher than DNA Chisel's 14 times and equal 10 times.
  It is never lower. All 48 sequences pass both checkers, with zero
  actionable GC-window and motif violations, and the protein matches the
  input including the stop.
- Under assembly, which is the set that binds, DNA Design Suite is higher
  on seven genes and tied at CAI 1 on hokC, where the unconstrained codon
  choice is already feasible. `cai_greedy` at CAI 1 is infeasible on the
  other seven.
- The regular-constraint sequences from pull request 31 are reproduced by
  SHA-256, including the cases the quality note identified as the pattern-DP
  optimum for genes of at most 2500 nt. lacZ and yeeJ sit past that limit
  and DNA Design Suite remains higher there, on the same sequences as before.
- DNA Chisel's rpsL assembly result depends on the seed. Seeds 1, 7 and 11
  reproduce three different sequences between CAI 0.95775539455032 and
  0.9641706286088719. DNA Design Suite's published rpsL assembly sequence
  stays at 0.9866866273575934.
- The same specification, run on the current tree, raised assembly CAI for
  waaS, lacZ and yeeJ without losing feasibility and without a per-gene
  change of settings.
- The tool runs locally, offline, under Apache-2.0, from a CLI, a Python
  API and a local page. The caller supplies the codon table. Constraints
  are a JSON document the checker and the search share. A job names its
  objective, and a request for an objective this build does not rank is
  refused rather than silently dropped (`CF-LAND-014` for GC as an
  objective, `CF-LAND-001` for an unknown term).
- The three negative fixtures still fail both checkers, so the constraint
  mapping used for the comparison is the mapping that rejects a bad
  sequence.

## Disadvantages of DNA Design Suite

- No commercial tool was run. Nothing in this file shows a higher expression
  level, a better vendor synthesis score, or a better CAI on GeneOptimizer's,
  IDT's, GenScript's or Twist's own table. Those sequences were not obtained.
- Assembly CAI is not shown to be optimal. rpsL remains 0.0038104428470379093
  under a relaxed bound that is itself infeasible. waaS remains under its
  published relaxed bound of 0.9903106518103836. The elastin case in the
  pull request 31 file, feasible CAI 0.4380255552056712 where DNA Chisel
  raised `NoSolutionError`, was not rerun; that CAI is far from 1, and the
  quality note says the true gap is unknown.
- Tier 1 cannot carry a GC window or a scattered k-mer (`CF-SRCH-017`). A
  tighter set can end as `CF-SRCH-006` with no sequence, as waaS did on the
  harsh assembly preset in the 17 October 2026 performance report. A user
  who wanted a best-effort infeasible sequence does not get one.
- There is no codon-harmonization objective and no objective that matches a
  codon-frequency distribution instead of maximizing CAI. DNA Chisel ships
  both (`harmonize_rca`, `match_codon_usage`). OPTIMIZER's paper and COOL's
  paper treat frequency matching and codon context as the point. Maximizing
  CAI piles onto the top codon. That is the objective this study measured,
  and it is the wrong objective for a caller who wanted harmonization.
- Cryptic splice sites, polyadenylation signals and rho-independent
  terminators are not built-in models. A caller who already knows the motif
  can forbid it. JCat's paper and the GeneOptimizer paper treat classes of
  those signals as part of the product. This study did not show DNA Design
  Suite removing a signal it was not given.
- Host tables are files on disk. There is no curated catalogue of hundreds
  of species. DNA Chisel will load `e_coli` or `h_sapiens` by name, and the
  academic servers ship their own lists. A vendor "E. coli" table is
  usually not this ribosomal relative-adaptiveness file, so the two CAI
  numbers are not interchangeable.
- GC content is a constraint (`Rule::GcWindow`), not a ranked objective.
  Global minimum free energy is not an objective. `structure.local_accessibility`
  needs a user-installed ViennaRNA, and a gene longer than one 40 nt span
  is not folded in tiers 1 and 2 (`CF-SRCH-010`). The accessibility column
  in the tables is a post-hoc score. LinearDesign remains a different
  problem, as pull request 31 already concluded.
- The context term is an E. coli Potts model when that file is installed,
  and `epistasis.dwell_deviation` stays `CF-LAND-001`.
- The recorded DNA Chisel walls include interpreter startup. On the short
  genes that startup is most of the half-second. The yeeJ assembly gap is
  larger than that, on one shot, on one release build. It is not a median
  and not a cold-start comparison of two installed packages.
- Diversity is the other side of the seed result. DNA Chisel returned three
  feasible rpsL assembly sequences from three seeds. Vendor tools often
  return more than one candidate. This protocol reports DNA Design Suite's
  single recommended sequence.
- The package version is 0.1.0, the package is not on PyPI, and a
  synthesized gene still has to pass the manufacturer's own screen. DNA
  Design Suite does not place an order.

## Reproduce

```sh
.venv/bin/maturin develop --release
.venv/bin/pip install 'dnachisel==3.2.16' 'ViennaRNA==2.7.2'
.venv/bin/python benchmarks/commercial_comparison/run_chisel_vs_dds.py \
    --out docs/benchmarks/commercial-comparison-rows.json
```

The runner is `benchmarks/commercial_comparison/run_chisel_vs_dds.py`.
It refuses to invent a score for a tool that does not run. Commercial
tools are outside it on purpose.
