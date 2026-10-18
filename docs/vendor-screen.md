# Vendor screen checklist

This is a dry run. It maps DNA Design Suite constraints onto synthesis-rule
classes that public descriptions of GeneArt, IDT, GenScript, and Twist name.
It does not transcribe those companies' numeric cutoffs, and it does not
report a manufacturer score. Wet-lab expression is out of scope.

The citation for the GeneOptimizer classes is Raab, Graf, Notka, Schödl and
Wagner, Syst Synth Biol 2010, 4:215–225, doi:10.1007/s11693-010-9062-3.
That paper's parameter classes are GC content, repeats, restriction sites,
and sequence motifs. IDT, GenScript, and Twist public descriptions also name
GC, repeats, and hairpins. Those companies' sites were not fetched for the
numbers below. Hairpins are not scored by this screen. Windowed accessibility
and the global fold are separate, labelled reports.

| Class | DNA Design Suite rule | Numeric cutoff in this screen | Where the number comes from |
| --- | --- | --- | --- |
| GC content | `c_gc50`, every 50 nt | GC fraction 0.40–0.60 | `examples/synthesis.json` |
| Homopolymer | `c_homopolymer` | a run longer than 6 | `examples/synthesis.json` |
| Repeats | `c_kmer12` | a 12-mer that occurs more than once | `examples/synthesis.json` |
| Restriction sites | `c_sites` | EcoRI, XbaI, SpeI, PstI, NotI, BsaI, BsmBI, both strands | assembly comparison set |
| Sequence motif | `c_sd` | `AGGAGG` on either strand | assembly comparison set |

`dds_synthesis` checks the first three rows. `dds_assembly` checks GC,
homopolymers, the restriction sites, the motif, and the 12-mer. Every row
sets `vendor_cutoff_transcribed` to false and `numeric_cutoff_source` to
`dna_design_suite preset`. The report sets `vendor_score` to false and
`wet_lab_expression` to `out_of_scope`.

A sequence this screen calls clear can still fail the manufacturer's own
screen. A sequence it flags is a candidate to edit before that screen, not
a prediction of a vendor's accept or reject.

```sh
.venv/bin/python -m dna_design_suite.cli.main vendor-screen --cds ATGAAATAA
.venv/bin/python -m dna_design_suite.cli.main fold --cds GGGAAACCC
.venv/bin/python -m dna_design_suite.cli.main hosts import \
    --input benchmarks/validation/kazusa/e_coli_316407.csv
```

The page buttons are **Vendor-class screen** and **Global fold**. Upload a
usage table accepts the same files as `hosts import`.
