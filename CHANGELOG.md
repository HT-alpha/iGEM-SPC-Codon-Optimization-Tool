# Changelog

## v2.4.0

No new fitted model and no change to the repair default.

- The dwell head stays `not_fitted`. Li et al. Cell 2014
  (doi:10.1016/j.cell.2014.02.033, GSE53767) score center residues at least
  10 nt from either end of 20–42 nt reads. They do not state a 5′-end offset.
  Li, Oh and Weissman, Nature 2012 (doi:10.1038/nature10965, GSE35641) use
  the same kind of center-weighting, 12 nt from either end, on 25–42 nt
  reads, and the 2014 methods say the size selection differs. Ingolia et al.
  Science 2009 is yeast. Mohammad, Green and Buskirk, eLife 2019
  (doi:10.7554/eLife.42591) assign their own 3′ ends. RiboToolkit and
  RiboCode do not state an offset for GSE53767. The GEO files were not
  reachable. `GET /v1/calibration` names these papers. Nothing was fitted.
- `P1-26` stays held (`DEP-SRCH-013`). The registry is still `usage.cai` and
  `usage.gc_content` (positionwise), `constraint.violation_count` (motif
  window, or unbounded for a k-mer rule), and `structure.local_accessibility`
  (`Codons(14)`, L0). Burden declares no radius. The Potts window `W = 10`
  is not a registered radius; raising the block from 8 to 10 would be the
  ordinary remedy if it were. The dwell head registers no term.
- EXP-05, 20 seeds per arm, 40 generations, stagnation window 41, on
  `gapA`, `pfkA`, `mdh`, `eno` and `pykF`. Each finished 40 generations.
  Hypervolume margin 0, Lamarckian standard deviation 0, on every protein.
  `decide()` was not edited and does not switch the default. `sucD` was not
  included: its interior motif is `CF-SRCH-006`, the same refusal as `rplL`.
- `GATE-BEN-001` (`DEP-BEN-003`) is a CAI-only comparison of 50 proteins.
  `claimed` is false. `aceB` and `thrA` contain `AGG`, missing from the
  ribosomal-protein weight table, so CodonForge and the native CDS score 0
  and `cai_greedy` scores 1. The other 48 are at least the greedy value.
  CDSfold, RNAplfold and LinearDesign were not added.

`P1-26`, EXP-06, LinearDesign, bundled thermodynamic tables and `P3-14`
remain held. Clearance is unchanged.

## v2.3.1

The PrimerForge sequence viewer selects the dragged region again.

In Chrome, `setPointerCapture` on the viewer retargets `pointermove` to the
viewer, so `event.target` is no longer the nucleotide and `indexAt` ignored
the move. A drag stayed on the one base where it started. The nucleotide is
now taken from `document.elementFromPoint` (`viewer.js`). Capture stays, so a
drag that leaves a nucleotide still reports moves. The same code is in
v2.2.0 and v2.3.0.

Confirmed in Chrome 148, on the first 700 nt of the shipped `lacZ` CDS: a
drag from index 300 to 345 selected `Region 300–346 (46 nt).` A drag back
from 80 to 40 selected `Region 40–81 (41 nt).` The viewer's Forward button
still pages a sequence longer than 10,000 nt. Before the change the same
300-to-345 drag selected `Region 300–301 (1 nt).`

`tools/viewer_controls.sh` puts that ignore-the-point behaviour back.
`test_a_retargeted_move_names_the_nucleotide_under_the_pointer` fails, and
the script reverts.

Clearance, `NOTICE`, and `test_amd_013` are unchanged. `P1-26`, EXP-06,
LinearDesign, bundled thermodynamic tables, and `P3-14` stay held.

## v2.3.0

A DNA double helix on the landing page and in the loading and result states of
both tabs. Frontend only: no backend change, no new dependency, no library. The
page gains one script file, `helix.js`, and one entry in the server's static
file list so it can be fetched.

- The helix is built from B-form geometry at run time (SVG, no image) and
  rotates by CSS keyframes. Rise 3.4 A and twist 36 deg per base pair follow
  Watson and Crick 1953. The phosphate radius (9.3 A) and strand offset
  (145 deg, so the grooves alternate) are measured on PDB 1BNA by
  `tools/measure_bform_1bna.py`, recorded in `docs/helix-1bna-measurement.json`.
  The rungs carry the first 20 nt of the `lacZ` CDS from the shipped corpus.
- Backbone in the `recommended` token (`DEP-UI-001`). Rungs in IGV's default
  nucleotide palette muted by one stated rule, in two halves per rung (A green
  and T red, G amber and C blue). IGV's G is orange; the brief's "yellow" is
  not what IGV ships.
- Four states, set only by the existing polling code. Idle: centred, one turn
  per 9 s. Running: in the corner, one turn per 5 s, rungs pulsing bottom to
  top. Result: finishes its turn, docks at 80 x 120 px, result fades in. Error:
  rotation stops, rungs turn `refusal-ink`, the message fades in 1.5 s later.
- `prefers-reduced-motion: reduce` is a stylesheet rule. The helix stays a
  still image at its idle pose.
- Cost: rotation is `transform` and `opacity` keyframes, about 1 ms of main
  thread per second in Chrome. A first version that animated a registered custom
  property cost about 115 ms per second and was replaced. A hidden tab does no
  work. There is no animation loop.
- `helix.js` is 11.4 kB and the SVG it builds is 17.1 kB, 28.5 kB together. The
  helix's keyframes and rules add 10 kB to the shared stylesheet.
- `python/tests/test_helix.py` has 18 tests, and `tools/helix_controls.sh` injects
  24 defects, one per property they guard; all 24 are caught.
- `app.js` is now loaded as a module. PrimerForge's poll loop gains the same
  `catch` the CodonForge one has, so a network failure there ends in the error
  state and a message instead of an unhandled rejection.
- Unchanged: every disclosure on both tabs, the tab colours, the typography, the
  spacing, the Pareto table. The helix takes its place in the page flow above
  the placeholder notes.

Found while testing and not changed: dragging across the PrimerForge sequence
viewer selects one nucleotide in Chrome, because `setPointerCapture` retargets
`pointermove` to the viewer and `indexAt()` then finds no `[data-index]`. It
does the same on v2.2.0.

Still not claimed: everything v2.2.0 does not claim. `P1-26`, EXP-06,
LinearDesign, bundled thermodynamic tables and `P3-14` remain held, and the
legal clearance gate is unchanged.

## v2.2.0

Close-out notes on top of v2.1.0 (`0f22940`). The converged Potts fit, the
50-protein corpus, and the Node 24 action pins are in that tag. This release
records them in the submission notes and states what is still not claimed.

Confirmed from v2.1.0:

- Selected `W = 10` Potts fit: 2333 iterations, gradient norm
  0.09770886731134135, `converged: true`. Digest
  `sha256:ac64cae265510c2cad9c3a85a59b9abf389cae456eba25fc4595e21d20e645b6`.
  Held-out nats improved by 0.0003 nats/codon (−0.8712 to −0.8709).
- Fifty `NC_000913.3` proteins, each with a protein_id from that accession.
  Selected-CDS digest
  `sha256:f03c0e882509ec34d8ad0fc327759dabc9b67cdf803a8892bdaeabf475a60bb3`,
  also recorded on `gates.json`. The Potts held-out set is still the
  original 13.
- `actions/checkout@v5`, `actions/setup-python@v6`, `upload-artifact@v6`,
  and `download-artifact@v7` declare Node 24.

New in this tag:

- The submission and `docs/primerforge-validation.md` state the converged
  fit and the 0.0003 nats observation.
- EXP-05 is two proteins. `rplL` is infeasible under the constraint set
  (`CF-SRCH-006`). That is a property of the gene under the motif, not a
  search bug. Margin was zero on `rpsJ` and `rpsL` at 20 generations.
  Lamarckian stays the default.
- The dwell head is not fitted. GSE53767 states no 5′-end A-site offset.
  `GET /v1/calibration` returns `dwell.state: not_fitted`.
- `GATE-BEN-001` stays unclaimed: global MFE is not windowed accessibility.
- PrimerForge added no landscape term. The tier-2 block is unchanged, and
  `P1-26` stays held.

Still not claimed:

- `GATE-BEN-001`.
- An EXP-05 result that separates Baldwinian from Lamarckian on hypervolume.
- A fitted dwell head.
- `P1-26`. Commercial distribution remains outside the AMD-013 sign-off.

## v2.1.0

- The selected `W = 10` Potts model is converged. Cumulative L-BFGS
  iterations from zero initialization: 2333. Final gradient norm
  0.09770886731134135 (reduced gradient 8.72e-5, tolerance 1e-4). Held-out
  nats/codon −0.8708682912309792, compared with −0.8711749285771091 at the
  40-iteration early stop. `GET /v1/calibration` reports `converged: true`.
  Digest `sha256:ac64cae265510c2cad9c3a85a59b9abf389cae456eba25fc4595e21d20e645b6`.
- The design panel is 50 `NC_000913.3` CDS (the original 13 plus 37, each
  cited by its protein_id). `GATE-BEN-001` stays unclaimed: LinearFold's
  global MFE is not `structure.local_accessibility`. The Potts held-out set
  is still those original 13. Selected-CDS digest
  `sha256:f03c0e882509ec34d8ad0fc327759dabc9b67cdf803a8892bdaeabf475a60bb3`.
- `actions/checkout@v5` and `actions/setup-python@v6` replace the Node 20
  majors. Wheel upload and download use Node 24 majors as well.
- PrimerForge polls primer tiers, drags a live selection, links a hovered
  pair to its binding sites, slides the detail panel in, recolours
  threshold chips from the values already on the pair, and names the
  filter that removed the most candidates on `CF-PRIM-001`.
- `epistasis.dwell_deviation` is refused with `CF-LAND-001` when the fitted
  head is absent. Li et al. Cell 2014 does not state a +15 nt A-site offset,
  so no GSE53767 counts were extracted and none were replaced with zero.
- `P1-26` stays held. The landscape registry still has no in-core L1 term
  whose radius exceeds the default tier-2 block and cannot be covered by
  raising that block. See `DEP-SRCH-013`.
- EXP-05 was re-run at 20 generations, 20 seeds per arm, on `rpsJ`, `rplL`,
  and `rpsL`. `decide()` is unchanged. Hypervolume margin is 0, so
  Lamarckian stays the default. `rplL` does not reach tier 3 (`CF-SRCH-006`
  on both arms).

## v2.0.0

The page is the DNA Design Suite. CodonForge remains the codon-optimization
tool. PrimerForge is the new primer-design tab.

New since v1.1.0:

- A shared shell with two tabs, Codon Optimization and Primer Design. A
  finished codon design can be handed to the primer tab as the input
  sequence.
- Primer search over the flanks of a selected region, and a synchronous
  check of one pasted pair (`POST /v1/primers/candidates`,
  `GET /v1/primers/{id}`, `GET /v1/primers/{id}/pair/{pair_id}`,
  `POST /v1/primers/analyze`).
- A sequence viewer with drag selection, a pair table, and a per-pair
  detail panel. Thermodynamics are SantaLucia 1998 Table 2.

Still not claimed:

- `GATE-BEN-001`.
- A converged Potts fit. The fit remains a 40-iteration early stop.
- An EXP-05 result that separates Baldwinian from Lamarckian on
  hypervolume. The re-run kept Lamarckian.
- Genome specificity (BLAST), Primer3, and mFold parity for hairpins.
