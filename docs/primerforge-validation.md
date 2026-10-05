# PrimerForge thermodynamic validation

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

## Dimer ΔG against published OligoAnalyzer examples

Nurhafitri and Suwirmen, Bioeksperimen Jurnal Penelitian Biologi (2026),
Table 2, report OligoAnalyzer ΔG for:

| Pair | Sequence | What they report |
| --- | --- | --- |
| 2_CaPIPF self-dimer | `CAGCCACTGATGCCAAGAGA` | self-dimer among −3.55, −1.47, −5.38, −1.57 kcal/mol (the HTML table merges the two rows) |
| 2_CaPIPR self-dimer | `TTGTGGAATGGCATGGCTCT` | same merged row |
| 2_CaPIP hetero-dimer | the two sequences above | heterodimer among −5.38 and −1.57 kcal/mol |

The article's Table 2 is not machine-readable as separate cells: the HTML
joins the forward and reverse rows. Those numbers are therefore not used as
a pass/fail gate. Owczarzy et al., Nucleic Acids Research 2008;36:W163
(PMC2447751) states that OligoAnalyzer hairpins are mFold, which adds loop
penalties Table 2 does not contain. Dimer figures from OligoAnalyzer also
include dangling-end terms (not in Table 2). A difference above 0.5 kcal/mol
against OligoAnalyzer is that missing physics, not a hidden residual.

Engine two-state ΔG from Table 2 at 1 M Na+, 37°C, for those sequences:

| Sequences | Engine ΔG°37 |
| --- | --- |
| `CAGCCACTGATGCCAAGAGA` with itself | −0.810 kcal/mol |
| `TTGTGGAATGGCATGGCTCT` with itself | −1.421 kcal/mol |
| the two sequences as a hetero-dimer | −5.796 kcal/mol |

The hetero-dimer is 0.42 kcal/mol from one of the merged published figures
(−5.38). The self-dimer of the reverse primer is 0.05 kcal/mol from −1.47.
The forward self-dimer (−0.810) is not within 0.5 kcal/mol of −3.55. That
published cell is not uniquely assigned to this sequence in the merged HTML
row, and Table 2 has no dangling-end bonus. The gap is reported here; it is
not averaged into the Tm check.

A `GCGCGCGCGCGCGCGCGCGC` self-dimer is in the refuse zone (ΔG < −9 kcal/mol).
