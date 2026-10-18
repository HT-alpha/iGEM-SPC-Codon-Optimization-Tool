# Primer design thermodynamic validation

Parameters are SantaLucia J Jr, Proc Natl Acad Sci USA 1998;95(4):1460–1465,
Table 2 (PMC19045). Tm uses Eq. 3 of that paper. Strand concentration and
sodium concentration are the fields in `crates/cf-primer-thermo/data/santalucia1998.json`
(1 µM, 1 M Na+). Hairpin loop initiation and dangling-end parameters are not
in Table 2 and are not added.

## Tm and ΔG against Table 2

Duplex in Fig. 1: `5'-CGTTGA-3'` with its complement. Hand sum of Table 2:

| Quantity | Hand sum of Table 2 | Engine |
| --- | --- | --- |
| ΔH° | −41.2 kcal/mol | −41.2 kcal/mol |
| ΔS° | −115.4 cal/K·mol | −115.4 cal/K·mol |
| ΔG°37 from Eq. 2 | −5.409 kcal/mol | −5.409 kcal/mol |
| Tm at 1 µM, Eq. 3, Ct/4 | 9.805 °C | 9.805 °C |

Table 1's rounded ΔG°37 sum for the same stacks is −5.35 kcal/mol. The 0.06
kcal/mol gap is Table 2 being published to 0.1 kcal and 0.1 cal. It is under
the 0.5 kcal/mol reporting threshold.

## Self-dimer of `CAGCCACTGATGCCAAGAGA`

Nurhafitri and Suwirmen, Bioeksperimen Jurnal Penelitian Biologi (2026),
Table 2, list an OligoAnalyzer self-dimer of −3.55 kcal/mol on the row that
contains this forward primer.

The minimum-ΔG register pairs the oligo with its reverse as:

```
5'- CAGCCACTG
    |||...|||
3'- GTCACCGAC
```

Two Watson–Crick helices, three mismatches between them. Table 2 stacks, in
order:

| Helix | Top dinucleotide | Table 2 key | ΔH° | ΔS° |
| --- | --- | --- | --- | --- |
| CAG | CA | CA/GT | −8.5 | −22.7 |
| CAG | AG | CT/GA | −7.8 | −21.0 |
| CTG | CT | CT/GA | −7.8 | −21.0 |
| CTG | TG | CA/GT | −8.5 | −22.7 |

Initiation is charged once, for the outer terminal pairs. Both are G·C
(Table 2 "Init. w/term. G·C": ΔH° 0.1, ΔS° −2.8), so two ends contribute
ΔH° 0.2 and ΔS° −5.6. The paired span is not self-complementary, so the
symmetry row (ΔS° −1.4) is not added. Sodium is 1 M, so Eq. 7 contributes 0.

ΔH° = −32.4 kcal/mol. ΔS° = −93.0 cal/K·mol.
ΔG°37 = −32.4 − 310.15 × (−93.0/1000) = −3.556 kcal/mol.

The published cell is −3.55 kcal/mol. The difference is 0.006 kcal/mol.

The previous engine stopped at the first helix (`CAG`, stacks CA/GT and
CT/GA, initiation at both of its ends) and reported −0.810 kcal/mol. That
dropped the second helix in the same register and charged initiation for a
fragment. The gap was that truncation, not a different parameter table.

## Other published cells

`TTGTGGAATGGCATGGCTCT` with itself. The minimum register is two helices,
`ATG` and `CAT`, stacks AT/TA, CA/GT, CA/GT, AT/TA, initiation A·T at both
outer ends (ΔH° 2.3 and ΔS° 4.1, twice). ΔH° −26.8, ΔS° −78.0,
ΔG°37 = −2.608 kcal/mol. A separate single helix, the self-complementary
tetramer `CATG` (stacks CA/GT, AT/TA, CA/GT, two G·C initiations, plus the
symmetry ΔS° −1.4), is −1.421 kcal/mol. The merged HTML also contains
−1.47, which is 0.05 kcal/mol from that tetramer. The row does not label
−1.47 as the reverse primer's minimum, and −1.47 is not the minimum under
Table 2.

The hetero-dimer is one helix of six pairs, stacks AT/TA, CA/GT, GC/CG,
GG/CC, CA/GT, initiation A·T at both ends. ΔG°37 = −5.796 kcal/mol. One
merged HTML cell is −5.38. The difference is 0.42 kcal/mol, under the
0.5 kcal/mol line. Hairpins remain outside this comparison: Owczarzy et al.,
Nucleic Acids Research 2008;36:W163 (PMC2447751), states that OligoAnalyzer
hairpins are mFold, and Table 2 has no hairpin-loop initiation.

A `GCGCGCGCGCGCGCGCGCGC` self-dimer is in the refuse zone (ΔG < −9 kcal/mol).
`CGCAAAGCG` with itself is the threshold check for the register fix: the
single helix `CGC` is −2.459 kcal/mol (not flagged); both helices with one
initiation are below −4 kcal/mol and are flagged.

## Potts surface at convergence

The selected `W = 10` Potts fit is 2333 total iterations, gradient norm
0.0977, `converged: true`. The held-out nats improved by 0.0003 nats/codon
over the early-stop model (−0.8712 to −0.8709), confirming the early-stop
model was already a good approximation.
