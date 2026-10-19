# Gibson assembly in the codon tool

The codon tool can do three things for Gibson assembly. All three use the
sequence you supply. None of them is tuned to a particular gene or organism.
None of them is evidence that a fragment will assemble, or that a construct
will be free of recombination in a cell.

## What you paste

Put the coding sequence in the CDS field. Put the upstream neighbor, such as
a promoter, in the 5′ flank. Put the downstream neighbor, such as a
terminator, in the 3′ flank. Those flanks are not rewritten.

The page button **Gibson** loads `examples/gibson.json`. The command is
`design --assembly gibson` with `--five-prime-flank` and
`--three-prime-flank`. `POST /v1/gibson` and the `gibson` command report
arms and repeats for a sequence you already have. They do not optimize it.

## 1. Repeats between the CDS and its neighbors

`examples/gibson.json` is one blocking rule: an exact 12-mer may occur at
most once on the 5′ flank, the CDS, and the 3′ flank read as one string.
Overlapping occurrences count. A 12-mer and its reverse complement are not
the same repeat.

This is the same 12-mer length as `examples/synthesis.json`. The checker
that already runs during a design counts it. A copy in the promoter and a
copy in the CDS is one repeat, and seed repair can change the CDS because
that copy is editable. A copy that never touches the CDS is left in place.
The job does not fail a candidate for a repeat it cannot edit.

12 nt is not a measured Gibson overlap and not a measured RecA length.
Gibson et al., Nature Methods 2009, 6:343–345, is the assembly method. A
protocol that reprints the ISO-buffer mix describes 20–150 bp overlaps as
suitable for that mix. Removing exact 12-mers can reject a design that
would have assembled, and it will not see a longer repeat that differs by
one base.

## 2. Homology arms

Each arm is a copy of the flank at the junction. The 5′ arm is the end of
the 5′ flank. The 3′ arm is the start of the 3′ flank. The default length
is 40 nt. Lengths from 20 to 150 are accepted, which is that protocol’s
range. This repository did not measure an optimal length.

If the flank is shorter than the length you asked for, the arm is not
designed. The code is `CF-GIB-001`. The missing bases are not filled in.
If you leave a flank empty, that side is `CF-GIB-002`: no arm, and the
repeat scan did not see a neighbor there.

The insert string is the arm or arms copied onto the ends of the CDS. It is
what you would order for those junctions. The repeat scan used the full
flanks, not only the arms. Melting temperature is not computed.

`gibson` prints the report and exits 3 when a supplied flank was too short
to copy. `POST /v1/gibson` returns 200 with the same problem inside the
report when the input itself was readable, and 422 with `CF-GIB-003` when
the sequence or the arm length is malformed.

## 3. Direct repeats anywhere in the supplied construct

The same 12-mer rule is the global scan. The report classes each repeat:

- `cds_context` — at least one copy touches the CDS and at least one copy
  does not. This is a repeat between the gene and a neighbor.
- `cds_internal` — every copy touches the CDS.
- `flank_only` — no copy touches the CDS. Codon optimization cannot remove
  it.

`actionable` is the checker’s flag. A flank-only repeat is not actionable.

## What is not supported

- Vector or part sequence you did not paste.
- A repeat that exists only because a circular plasmid joins the outer ends.
- Mismatches, gaps, and reverse complements. Inverted repeats are not
  treated as recombination sites.
- Overlap melting temperature, chew-back, fill-in, or how much DNA to add.
- A promise that the assembly will work.

`guarantees_assembly` on the report is false.
