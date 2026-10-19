# DNA Design Suite and DNA Chisel on the GATE-BEN-001 proteins

One head-to-head, measured at commit `499f0f3c9c10721d2c0aa225c7c641bd27fa1ebd`. The sequences and the scores are in `docs/benchmarks/gate-ben-001-dnachisel-rows.json`, written by this run. Rewriting this note from that JSON does not re-run either optimizer. CAI in this note is the Sharp and Li geometric mean of a fixed weight table. It is a sequence statistic. It is not a predicted expression level, and this note makes no wet-lab claim.

The published `GATE-BEN-001` claim is a different comparison: the same 50 proteins and the same ribosomal weights, an empty constraint list, and DNA Design Suite against `cai_greedy`. This note does not replace that claim and does not extend it. The hard constraints here are `benchmarks/validation/constraints_assembly.json`, given to both tools.

## Versions, seeds, and the machine

| piece | pin |
| --- | --- |
| git commit measured | `499f0f3c9c10721d2c0aa225c7c641bd27fa1ebd` |
| DNA Design Suite package | 0.1.0.post1 |
| DNA Design Suite core | 0.1.0 |
| DNA Chisel | 3.2.16 |
| scoring Biopython | 1.81 |
| project-environment Biopython | 1.88 |
| numpy | 2.5.4 |
| Python, design | 3.12.3 |
| Python, scoring | 3.12.3 |
| DNA Design Suite seed | 7 |
| DNA Chisel numpy seed | 7 |
| DNA Chisel `PYTHONHASHSEED` | 7 |

Machine: Intel(R) Xeon(R) Processor, 4 CPUs, 16387280 KiB RAM, Linux-6.12.94+-x86_64-with-glibc2.39. Genes ran one at a time. DNA Design Suite time is `time.perf_counter` around `design` in an already imported process, after one untimed warmup call. DNA Chisel solver time is `time.perf_counter` inside a fresh interpreter, around problem construction, `resolve_constraints`, and `optimize`. DNA Chisel wall time is the parent clock around that process, so it includes interpreter startup. Neither clock is a packaging claim about another machine.

## What both tools were given

Fifty CDS from `bm_ecoli_characterized_50`, accession `NC_000913.3`, selected-CDS sha256 `f03c0e882509ec34d8ad0fc327759dabc9b67cdf803a8892bdaeabf475a60bb3`. Host `escherichia_coli:511145/MG1655`, translation table 11. Weights are `examples/mg1655_ribosomal_usage.json`, sha256 `25b7b0578bf0c5f60f85cefbdc75fb9722d6abf84eac0d8167d06ff7b44980dd`. They match `ribosomal_usage()` on the shipped CDS gzip. Every amino-acid family in that table has maximum 1. The smallest positive weight is 0.002392. Weight 0: AGG, TAG.

DNA Design Suite settings are the `GATE-BEN-001` runner's settings, with the assembly file in place of the empty list: budget `fast`, generations 4, seed 7, segment 2, stagnation window 8, objective `usage.cai` only. Budget `fast` runs tier 1 and tier 2. It does not run the generational tier, so the generation count is not extra search.

DNA Chisel settings are library defaults apart from the seed and a silent logger: `max_random_iters` 1000, `mutations_per_iteration` 2, `optimization_stagnation_tolerance` 100, `randomization_threshold` 10000. The objective is `CodonOptimize(method="use_best_codon")` on this weight table. DNA Chisel's own score is `N * log` over every codon, including Met, Trp, and stops, and it replaces a zero with 0.001 before the logarithm. The smallest positive weight is above 0.001, so that floor does not reorder the positive codons. The CAI column below does not use that score.

## Constraint translation

| id | rule | DNA Chisel specification |
| --- | --- | --- |
| `c_gc50` | gc_window 50 nt, 0.4 to 0.6, every start | EnforceGCContent(mini=0.4, maxi=0.6, window=50). The implementation returns len(sequence) - window + 1 windows at step 1. Bounds are inclusive. |
| `c_bsai` | restriction_site BsaI GGTCTC, both strands, max 0 | AvoidPattern(GGTCTC). With no location the strand is 0, and initialization copies that strand onto the full sequence, so the reverse complement is rejected. A palindrome is searched once. |
| `c_bsmbi` | restriction_site BsmBI CGTCTC, both strands, max 0 | AvoidPattern(CGTCTC). With no location the strand is 0, and initialization copies that strand onto the full sequence, so the reverse complement is rejected. A palindrome is searched once. |
| `c_ecori` | restriction_site EcoRI GAATTC, both strands, max 0 | AvoidPattern(GAATTC). With no location the strand is 0, and initialization copies that strand onto the full sequence, so the reverse complement is rejected. A palindrome is searched once. |
| `c_xbai` | restriction_site XbaI TCTAGA, both strands, max 0 | AvoidPattern(TCTAGA). With no location the strand is 0, and initialization copies that strand onto the full sequence, so the reverse complement is rejected. A palindrome is searched once. |
| `c_spei` | restriction_site SpeI ACTAGT, both strands, max 0 | AvoidPattern(ACTAGT). With no location the strand is 0, and initialization copies that strand onto the full sequence, so the reverse complement is rejected. A palindrome is searched once. |
| `c_psti` | restriction_site PstI CTGCAG, both strands, max 0 | AvoidPattern(CTGCAG). With no location the strand is 0, and initialization copies that strand onto the full sequence, so the reverse complement is rejected. A palindrome is searched once. |
| `c_noti` | restriction_site NotI GCGGCCGC, both strands, max 0 | AvoidPattern(GCGGCCGC). With no location the strand is 0, and initialization copies that strand onto the full sequence, so the reverse complement is rejected. A palindrome is searched once. |
| `c_homopolymer` | homopolymer base any, max run 6 | AvoidPattern(HomopolymerPattern(base, 7)) for each of ACGT. A match of that length is a run longer than the limit. Strand 0 also sees the reverse complement, which is the complementary base. |
| `c_sd` | forbidden_motif AGGAGG, both strands, max 0 | AvoidPattern(AGGAGG) on both strands, same path as a restriction site. |
| `c_kmer12` | kmer_repeat k=12, max 1, overlapping True, canonical_strand False | UniquifyAllKmers(k=12, include_reverse_complement=False). Occurrences are every start in the reference span. canonical_strand false keeps a k-mer distinct from its reverse complement. DNA Chisel's docstring says this constraint may fail to resolve a k-mer that already occurs more than twice; that outcome is recorded as no sequence. |

Every constraint in the shared file was sent. None of those rules was dropped.

The file's `repair_policy` is empty. It is search configuration, not a sequence predicate, and it was not sent.

Schema kinds that are not in the shared file were not given to either tool:

- `local_structure`: No DNA Chisel specification is this predicate. AvoidHairpins is a hairpin rule, not a windowed unpaired probability or a pinned minimum free energy, and it does not name the folding engine. The shared file does not contain this kind, so neither tool was given it.
- `codon_ban`: DNA Chisel has no hard specification that removes a listed codon and keeps the amino acid. A zero in the codon table changes the CAI objective and still leaves the codon legal under translation. The shared file does not contain this kind, so neither tool was given it.
- `min_codon_frequency`: AvoidRareCodons compares a frequency in DNA Chisel's codon table. The shared numbers are relative adaptiveness, not frequencies, so the same cutoff would be a different rule. The shared file does not contain this kind, so neither tool was given it.
- `custom_regex`: DNA Chisel SequencePattern can carry a regular expression, so a translation of this kind exists. The shared file does not contain one, so neither tool was given one.

## How the sequences were scored

CAI is `Bio.SeqUtils.CodonAdaptationIndex.calculate` from Biopython 1.81. The reference weights are the relative-adaptiveness values in `examples/mg1655_ribosomal_usage.json` (sha256 25b7b0578bf0c5f60f85cefbdc75fb9722d6abf84eac0d8167d06ff7b44980dd). `calculate` skips ATG and TGG. TAA, TAG, and TGA are omitted from the index, so they are skipped as well. A counted codon of weight 0 is reported as 0, because the geometric mean of a set that contains 0 is 0 and `calculate` raises on `log(0)`. The deprecated `Bio.SeqUtils.CodonUsage` class, which divides by `n - 1`, is not used. `SharpEcoliIndex` is not used.

Global GC is `Bio.SeqUtils.gc_fraction` from the same Biopython. Each hard GC window is the same function on every 50 nt slice, step 1, and the window fails when that value is below 0.40 or above 0.60. The endpoints are included.

Each hard constraint is counted by a string scan in `independent_score.py`: overlapping literal matches for sites and motifs, the reverse complement as well when the rule says both strands (a palindrome is counted once), one violation per homopolymer run longer than 6, and one violation per distinct 12-mer whose overlapping starts exceed 1. A kind with no scan raises. The product's constraint counter is stored beside these counts and is not the headline.

Each CDS is translated with `Bio.Seq.Seq.translate(table=11)`, including the terminal stop as `*`. The design matches when that string equals the translation of the input CDS. A different stop codon still matches. An amino-acid change does not.

On every returned sequence the product's pass/fail agreed with the plain checker. The tables still use the plain checker. Across 98 sequences the largest absolute gap between `usage.cai` and this Biopython CAI is 0.0. That gap is recorded. It is not the headline.

## Losses and ties

A loss is a metric on which DNA Chisel is strictly better on at least one gene, or a total on which DNA Chisel is strictly better. A tie is equality, with CAI equality at an absolute difference of at most 1e-12. Global GC is reported as higher, lower, or equal, and is not called a win or a loss. Genes are named in panel order.

CAI where both tools returned a sequence, whether or not the sequence passes the hard constraints: DNA Design Suite higher on 49, DNA Chisel higher on 0, tie on 0 (of 49).

DNA Chisel CAI is higher on no gene in that set, and there is no CAI tie.

Closest CAI pair among genes both tools returned: sucD, absolute difference 0.00855795017482086, DNA Design Suite 1.0, DNA Chisel 0.9914420498251791.

No CAI is reported where a tool returned no sequence: aceF.

CAI where both sequences are feasible (translation matches and every hard constraint passes): DNA Design Suite higher on 49, DNA Chisel higher on 0, tie on 0 (of 49). Paired means on that set: DNA Design Suite 0.9802699218235877, DNA Chisel 0.9144569845305929.

Among feasible pairs, DNA Chisel CAI is higher on no gene and there is no CAI tie.

A reported CAI of 0 is the geometric mean when a counted codon has weight 0. DNA Chisel's own objective replaces that weight with 0.001 before the logarithm, so the codon stays available to its search. The shared constraints do not ban it. In-frame weight-0 codons in returned sequences:

- rpoB, dnachisel: AGG × 1; CAI 0.0; feasible True

Feasible sequences: both 49, DNA Design Suite only 0, DNA Chisel only 0, neither 1.

Feasibility is a tie: there is no gene where one tool is feasible and the other is not. Neither tool is feasible on aceF.

Translation back to the input protein: DNA Design Suite 49 of 49 returned sequences, DNA Chisel 49 of 49.

Translation is a tie on every returned sequence: each one translates back to its input protein.

Returned a sequence: both 49, DNA Design Suite only 0, DNA Chisel only 0, neither 1.
Neither tool returned a sequence: aceF.

Stored failure text, from the measurement JSON:

- aceF, dna_design_suite, NoSequence: state failed
- aceF, dnachisel, NoSolutionError: While solving EnforceGCContent[0-1893](mini:0.40, maxi:0.60, window:50) in 824-930: Random search did not find a solution in the given number of attempts. Try to increase the number of attempts with: problem.max_random_iters = 5000 # or even 10000, 20000, etc. If the problem persists, you may be in presence of a complex or unsolvable problem.

Per hard constraint, on genes where both tools returned a sequence:

| constraint | both pass | suite fails, Chisel passes | Chisel fails, suite passes | both fail | suite violation sum | Chisel violation sum |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `c_gc50` | 49 | 0 | 0 | 0 | 0 | 0 |
| `c_bsai` | 49 | 0 | 0 | 0 | 0 | 0 |
| `c_bsmbi` | 49 | 0 | 0 | 0 | 0 | 0 |
| `c_ecori` | 49 | 0 | 0 | 0 | 0 | 0 |
| `c_xbai` | 49 | 0 | 0 | 0 | 0 | 0 |
| `c_spei` | 49 | 0 | 0 | 0 | 0 | 0 |
| `c_psti` | 49 | 0 | 0 | 0 | 0 | 0 |
| `c_noti` | 49 | 0 | 0 | 0 | 0 | 0 |
| `c_homopolymer` | 49 | 0 | 0 | 0 | 0 | 0 |
| `c_sd` | 49 | 0 | 0 | 0 | 0 | 0 |
| `c_kmer12` | 49 | 0 | 0 | 0 | 0 | 0 |

`c_gc50` violation sums tie at 0.

`c_bsai` violation sums tie at 0.

`c_bsmbi` violation sums tie at 0.

`c_ecori` violation sums tie at 0.

`c_xbai` violation sums tie at 0.

`c_spei` violation sums tie at 0.

`c_psti` violation sums tie at 0.

`c_noti` violation sums tie at 0.

`c_homopolymer` violation sums tie at 0.

`c_sd` violation sums tie at 0.

`c_kmer12` violation sums tie at 0.

Global GC, both returned a sequence, not an objective: DNA Design Suite higher on 35, DNA Chisel higher on 10, equal on 4.

DNA Chisel global GC higher, not an objective: pykF, rplA, rplL, icd, sucD, pgi, gpmA, glpK, leuA, crr.

Global GC equal: eno, rpsL, mdh, metK.

Solver time, sum of the 50 genes: DNA Design Suite 18.616265151001244 s, DNA Chisel 13.283595282001215 s. DNA Design Suite faster on 8 genes, DNA Chisel faster on 42, equal on 0.

DNA Design Suite solver time faster: lacZ, aceF, sucD, sdhA, ppc, glpK, guaB, rpoB.

DNA Chisel solver time faster: gapA, eno, pykF, pfkA, adk, dnaK, tufB, rpsA, rplA, rplL, rpsL, rpsJ, aceA, aceB, aceE, lpd, gltA, icd, mdh, sucA, sucC, pgi, tpiA, pgk, gpmA, zwf, tktA, ptsG, aspC, ilvC, leuA, thrA, lysC, metK, purA, ndk, atpA, atpD, nuoC, crr, fadB, fabB.

On total solver time, DNA Design Suite is slower (18.616265151001244 s against 13.283595282001215 s).

Wall time, which includes DNA Chisel's fresh interpreter and does not include a fresh interpreter for DNA Design Suite: sums 18.616265151001244 s and 26.193951076000985 s. DNA Chisel's wall is lower on 9 genes. That comparison is not the same boundary on both sides; it is reported anyway.

DNA Chisel wall time lower: tufB, rpsA, lpd, mdh, sucA, sucC, atpA, atpD, fadB.

## Per gene

CAI and GC are the independent scores. Feasible means the translation matches and every hard constraint passes. Solver seconds are the clocks defined above.

| gene | nt | suite CAI | Chisel CAI | CAI | suite GC | Chisel GC | suite feasible | Chisel feasible | suite violations | Chisel violations | suite solver s | Chisel solver s |
| --- | ---: | ---: | ---: | --- | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: |
| lacZ | 3075 | 0.9665531384697095 | 0.8852227588514381 | suite higher | 0.5291056910569105 | 0.5160975609756098 | True | True | 0 | 0 | 0.8405141189998631 | 1.0948165570000583 |
| gapA | 996 | 0.9910560502192561 | 0.9525130693799817 | suite higher | 0.46887550200803213 | 0.46686746987951805 | True | True | 0 | 0 | 0.19404966700017212 | 0.08153227500019966 |
| eno | 1299 | 0.9844322624589411 | 0.9729019581768249 | suite higher | 0.4765204003079292 | 0.4765204003079292 | True | True | 0 | 0 | 0.2310203710001133 | 0.11051584699998784 |
| pykF | 1413 | 0.9832102450598392 | 0.9529529851285281 | suite higher | 0.4720452937013447 | 0.47841472045293704 | True | True | 0 | 0 | 0.2853080840000075 | 0.1675964659998499 |
| pfkA | 963 | 0.9756481568201967 | 0.9291547242760946 | suite higher | 0.4994807892004154 | 0.49117341640706125 | True | True | 0 | 0 | 0.3536124239999481 | 0.10690456799989079 |
| adk | 645 | 0.990191249861141 | 0.9273664773354529 | suite higher | 0.4883720930232558 | 0.48527131782945737 | True | True | 0 | 0 | 0.11819352700013042 | 0.05683120400021835 |
| dnaK | 1917 | 0.9738771227240254 | 0.9442130207707905 | suite higher | 0.48930620761606675 | 0.4846113719353156 | True | True | 0 | 0 | 0.42511344400008966 | 0.26956821500016304 |
| tufB | 1185 | 0.9849780727246936 | 0.9106853775816897 | suite higher | 0.4936708860759494 | 0.4869198312236287 | True | True | 0 | 0 | 0.3996249980000357 | 0.12442403000000013 |
| rpsA | 1674 | 0.9798524648357165 | 0.9366401387600073 | suite higher | 0.47371565113500597 | 0.4731182795698925 | True | True | 0 | 0 | 0.5371542640000371 | 0.23005067999997664 |
| rplA | 705 | 0.9820873720549521 | 0.9432325455985238 | suite higher | 0.475177304964539 | 0.4765957446808511 | True | True | 0 | 0 | 0.26239498399991135 | 0.059847192000006544 |
| rplL | 366 | 0.9728619943579154 | 0.9527743357878348 | suite higher | 0.46174863387978143 | 0.4644808743169399 | True | True | 0 | 0 | 0.11518215600017356 | 0.02904625699989083 |
| rpsL | 375 | 0.9866866273575934 | 0.9641706286088719 | suite higher | 0.48 | 0.48 | True | True | 0 | 0 | 0.09482435099994291 | 0.037162046000048576 |
| rpsJ | 312 | 0.9894463530334181 | 0.9486387847582942 | suite higher | 0.5096153846153846 | 0.5 | True | True | 0 | 0 | 0.048884944000064934 | 0.029503391999924133 |
| aceA | 1305 | 0.9839208596069 | 0.9155768236133295 | suite higher | 0.5065134099616858 | 0.496551724137931 | True | True | 0 | 0 | 0.30933311600006164 | 0.23923738600001343 |
| aceB | 1602 | 0.9813437623830314 | 0.9222569640521062 | suite higher | 0.50187265917603 | 0.5006242197253433 | True | True | 0 | 0 | 0.4277524629999334 | 0.35283433900008276 |
| aceE | 2664 | 0.9767140881225882 | 0.9135698911332043 | suite higher | 0.4954954954954955 | 0.4906156156156156 | True | True | 0 | 0 | 0.6907008979999318 | 0.5012702629999239 |
| aceF | 1893 | — | — | missing | — | — | False | False | — | — | 0.10675761100014824 | 0.11736487200005286 |
| lpd | 1425 | 0.9685717372981407 | 0.932679112549741 | suite higher | 0.4870175438596491 | 0.4849122807017544 | True | True | 0 | 0 | 0.5050370459998703 | 0.18543428700013465 |
| gltA | 1284 | 0.9840070152766143 | 0.9374333180237151 | suite higher | 0.4898753894080997 | 0.48909657320872274 | True | True | 0 | 0 | 0.28507988799992745 | 0.1895321290000993 |
| icd | 1251 | 0.9610658647053114 | 0.9175388736447924 | suite higher | 0.49160671462829736 | 0.4924060751398881 | True | True | 0 | 0 | 0.2982237270000496 | 0.16762552599993796 |
| mdh | 939 | 0.9845385983651036 | 0.9639740179768829 | suite higher | 0.49307774227902024 | 0.49307774227902024 | True | True | 0 | 0 | 0.3470889990001069 | 0.10684145000004719 |
| sucA | 2802 | 0.9706632300423321 | 0.9163196766550785 | suite higher | 0.5060670949321913 | 0.49964311206281226 | True | True | 0 | 0 | 1.0173672379999061 | 0.753341219000049 |
| sucC | 1167 | 0.9752324104469563 | 0.8981918435711392 | suite higher | 0.5021422450728363 | 0.49271636675235647 | True | True | 0 | 0 | 0.4224046179999732 | 0.17984420900006626 |
| sucD | 870 | 1.0 | 0.9914420498251791 | suite higher | 0.47701149425287354 | 0.4781609195402299 | True | True | 0 | 0 | 0.08618949700007761 | 0.09903872799986857 |
| sdhA | 1767 | 0.9744750059980389 | 0.9120597788258858 | suite higher | 0.5110356536502547 | 0.5019807583474816 | True | True | 0 | 0 | 0.28477458800011846 | 0.3691217450000295 |
| pgi | 1650 | 0.9838637131017953 | 0.9403078947772378 | suite higher | 0.4824242424242424 | 0.48484848484848486 | True | True | 0 | 0 | 0.32994001399993067 | 0.2834073730000455 |
| tpiA | 768 | 0.9923411136139978 | 0.9605203031506071 | suite higher | 0.4947916666666667 | 0.4921875 | True | True | 0 | 0 | 0.09427984400008427 | 0.060048966000067594 |
| pgk | 1164 | 0.9768981374212649 | 0.9560972003997259 | suite higher | 0.4879725085910653 | 0.4802405498281787 | True | True | 0 | 0 | 0.19757591700022203 | 0.12969954100003633 |
| gpmA | 753 | 0.9758995326648556 | 0.923798643217927 | suite higher | 0.47941567065073043 | 0.4807436918990704 | True | True | 0 | 0 | 0.1501087860001462 | 0.06973873600009028 |
| ppc | 2652 | 0.9673482753180728 | 0.8801563151065669 | suite higher | 0.5120663650075414 | 0.503393665158371 | True | True | 0 | 0 | 0.720721838999907 | 0.7649371840000185 |
| zwf | 1476 | 0.9822282365547873 | 0.9486253483428443 | suite higher | 0.49796747967479676 | 0.497289972899729 | True | True | 0 | 0 | 0.33776252799998474 | 0.28664989000003516 |
| tktA | 1992 | 0.9785193781139279 | 0.9301436597050051 | suite higher | 0.5040160642570282 | 0.5 | True | True | 0 | 0 | 0.4408724520001215 | 0.3647771729999931 |
| ptsG | 1434 | 0.9884347045135841 | 0.951222525904854 | suite higher | 0.49860529986053 | 0.49511854951185497 | True | True | 0 | 0 | 0.36654683199981264 | 0.25084732699997403 |
| glpK | 1509 | 0.9786954484471359 | 0.9165452299530452 | suite higher | 0.49039098740888004 | 0.49502982107355864 | True | True | 0 | 0 | 0.20673882500000218 | 0.30271121600003426 |
| aspC | 1191 | 0.9702180679759324 | 0.9233005442280394 | suite higher | 0.5029387069689337 | 0.4987405541561713 | True | True | 0 | 0 | 0.45492537900008756 | 0.17045301699999982 |
| ilvC | 1476 | 0.985154970507772 | 0.9245019249186401 | suite higher | 0.48848238482384826 | 0.4878048780487805 | True | True | 0 | 0 | 0.30015651600001547 | 0.23999591000006149 |
| leuA | 1572 | 0.9803733709144183 | 0.9273289965462365 | suite higher | 0.48282442748091603 | 0.4847328244274809 | True | True | 0 | 0 | 0.4999863609998556 | 0.315225952999981 |
| thrA | 2463 | 0.980488168033101 | 0.9084873561814071 | suite higher | 0.5038570848558668 | 0.5006090133982948 | True | True | 0 | 0 | 0.7606097750001481 | 0.6944358600001124 |
| lysC | 1350 | 0.9785786059383567 | 0.9061556076782342 | suite higher | 0.5051851851851852 | 0.5044444444444445 | True | True | 0 | 0 | 0.3743857530000696 | 0.27942480000001524 |
| metK | 1155 | 0.9795084549179818 | 0.9287464628868662 | suite higher | 0.5012987012987012 | 0.5012987012987012 | True | True | 0 | 0 | 0.29297606000000087 | 0.15248201700001118 |
| purA | 1299 | 0.9850517505359427 | 0.9364170221953664 | suite higher | 0.5042340261739799 | 0.5019245573518091 | True | True | 0 | 0 | 0.25473749299999326 | 0.19006500899990897 |
| guaB | 1467 | 0.9842000261097892 | 0.9562116390549741 | suite higher | 0.4948875255623722 | 0.49011588275391954 | True | True | 0 | 0 | 0.25495761299998776 | 0.2565019650000977 |
| ndk | 432 | 0.9755776177351668 | 0.9079783510535328 | suite higher | 0.48842592592592593 | 0.4861111111111111 | True | True | 0 | 0 | 0.07800800400013941 | 0.04341558999999506 |
| atpA | 1542 | 0.985090689899935 | 0.9606581875957805 | suite higher | 0.5006485084306096 | 0.5 | True | True | 0 | 0 | 0.6265249780001341 | 0.2252051410000604 |
| atpD | 1383 | 0.9764850192169936 | 0.9388596825855706 | suite higher | 0.4953000723065799 | 0.4895155459146782 | True | True | 0 | 0 | 0.48388831800002663 | 0.21009536300016407 |
| nuoC | 1791 | 0.9669916042982467 | 0.8996839319014238 | suite higher | 0.5041876046901173 | 0.5019542155220548 | True | True | 0 | 0 | 0.42935521499998686 | 0.34972749500002465 |
| crr | 510 | 0.986362231036357 | 0.9618046554499088 | suite higher | 0.44509803921568625 | 0.4470588235294118 | True | True | 0 | 0 | 0.202416926000069 | 0.04455886499999906 |
| fadB | 2190 | 0.9800139603180573 | 0.9259991377397305 | suite higher | 0.5041095890410959 | 0.5018264840182648 | True | True | 0 | 0 | 0.8343258089998926 | 0.5763454589998673 |
| fabB | 1221 | 0.9904700027348534 | 0.9533324665401405 | suite higher | 0.49713349713349714 | 0.4914004914004914 | True | True | 0 | 0 | 0.29781225700003233 | 0.2162137499999517 |
| rpoB | 4029 | 0.9830194072110617 | 0.0 | suite higher | 0.49094068006949615 | 0.4869694713328369 | True | True | 0 | 0 | 0.9400646349999988 | 1.1473468000001503 |

## Reproduce

```sh
make setup && make build
.venv/bin/pip install 'dnachisel==3.2.16'
python3 -m virtualenv benchmarks/gate_ben_001_chisel/.score-venv
benchmarks/gate_ben_001_chisel/.score-venv/bin/pip install 'biopython==1.81'
benchmarks/gate_ben_001_chisel/.score-venv/bin/python \
    benchmarks/gate_ben_001_chisel/independent_score.py --self-test
.venv/bin/python benchmarks/gate_ben_001_chisel/run_benchmark.py \
    --out docs/benchmarks/gate-ben-001-dnachisel-rows.json \
    --doc docs/benchmarks/gate-ben-001-dnachisel.md
```

Negative control, before this measurement. `independent_score.py --self-test` passed. The reverse-complement count in `_motif_hits` was then replaced with `hits += 0` and the self-test exited 1 with `BsaI reverse complement GAGACC was not counted`. That edit was reverted and the self-test passed again. The fixture is the antisense BsaI site, so a checker that only reads the forward strand does not pass.
