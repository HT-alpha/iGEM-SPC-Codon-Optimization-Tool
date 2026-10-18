# Screenshot captions

Captured 17 October 2026 from `python -m dna_design_suite.api` on this tree
(ViennaRNA 2.7.2). These are page captures of the rpsL demo. They are not
the assembly-constraint benchmark vectors in
`docs/benchmarks/`.

| File | Caption |
| --- | --- |
| `dds-rpsl.png` | DNA Design Suite after **Load the rpsL demo** and **Run design** (budget Fast, objective CAI). Status succeeded. `usage.cai` reads 1.0 and the constraints chip says satisfied. `examples/rpsL_constraints.json` checks nothing, so 1.0 means every counted codon is a top synonym of the loaded ribosomal table. |
| `dds-primers-rpsl.png` | Primer design after **Load the rpsL sequence** (bases 41–300, flank 40 nt) and **Find primer pairs**. The table’s first row in this capture is `pair_00_00000029`, product 331, worst dimer ΔG −5.35. |
| `dds-pareto.png` | The same rpsL demo with **Local mRNA structure** ticked and budget **Standard**. The front table lists `usage.cai` and `structure.local_accessibility` as Pareto axes, with `constraint_report.feasible` true. The recommended row in this capture is `cand_e5729347b7ec219c` (CAI 0.994611, structure 2.664833). |
