# Changelog

## Unreleased

Package version stays 0.1.0. This note is the DNA Design Suite delta on top
of v2.9.2.

Wheels build on push. `.github/workflows/wheels.yml` uploads a wheel for
Linux x86_64, Linux aarch64, macOS x86_64, macOS aarch64, and Windows x64,
plus an sdist. One `cp311-abi3` wheel covers CPython 3.11 through 3.13.
The workflow does not publish to PyPI. The README says how to install a
downloaded wheel, and how Hugo uploads it with twine from his own token.

Tier 3 on yeeJ with structure was spending its time in repair, not in the
fold. On a release build of this machine, before the k-mer change:
`repair` 6 calls, 2.013 s; `l0_batch` 2 calls, 0.656 s; dominance sort and
hypervolume under 0.001 s; wall 5.125 s. After recording positions only for
k-mers already over the limit, and reusing one k-mer index across a
codon-sized probe: `repair` 0.746 s; `l0_batch` 0.658 s; wall 3.060 s. The
vector stays `[0.9575892496987667, 2.408954059933615]`. Fold counts stay
1 full fold, 21 local edits, 12 hits. The v2.9.2 note's 3.202 s tier 3 was
`opt-level = 1`, which is what `make build` produces. These two numbers are
the release profile, measured the same way on both sides. Merged with
v2.9.3, the same release job is wall 1.977 s. Tier 3 is still 1.480 s
(`repair` 0.752 s, `l0_batch` 0.656 s). The wall drop is the v2.9.3 seed
(0.157 s) and `explain_baseline` (0.003 s).

Primer search still returns `product_sequence` on every pair. That document
is still 166429 bytes and the same SHA-256. `GET /v1/primers/{id}` and
`GET /v1/primers/{id}/pair/{pair_id}` omit the field unless
`?product_sequence=1` (also `true` or `yes`). On the yeeJ 300 nt flank job
the HTTP job document goes from 169009 bytes to 31980, and one pair goes
from 8480 bytes to 1613. The page slices the template it already submitted.
The recommended pair is still `TGACCGACAAATATTATG` / `GAACTACAGTGGTTTAAA`,
product 6843, ΔG −1.283850000000001.

## v2.9.2

DNA Design Suite note for this tree. The package version stays 0.1.0.
`python -m codonforge.cli.main --version` prints 0.1.0. The profile and the
ranked backlog are `docs/benchmarks/v2.9.2-profiles.md`.

Seed repair under the assembly rules was scanning the whole CDS, including a
fresh 12-mer index, on every pattern trial. yeeJ (7077 nt, seed 7, ribosomal
weights, generations 0) went from 20.814 s to 3.048 s for the design and from
19.983 s to 1.772 s for `resolve_feasible_seed`. The CDS is unchanged
(SHA-256 `686adf6a7ed03ac848270c94afbb9f7bbfd7873dbf926bc03390001d80eb9b4c`)
and the CAI vector is still `[0.9580065306093907]`. rpoB seed went from
4.790 s to 1.585 s, same vector `[0.9798209521813181]`. The homopolymer
tie-break is a local delta; 34064 full scans on yeeJ went from 0.815 s to
0.002 s.

Primer search computes hairpin, self-dimer, and Tm once per oligo, and the
Pareto comparison no longer allocates a copy of each objective vector. A
capped yeeJ search (flank 300, 1600 pairs) went from 0.143 s to 0.064 s.
The 300 nt CDS-flank document is the same bytes.

## v2.9.1

DNA Design Suite note for this tree: CodonForge and PrimerForge share one
page. The three measured changes are the CAI-only gate, the structure fold
cache, and the operon ledger. The numbers are in
`docs/benchmarks/v2.9.1-validation.md` and
`docs/benchmarks/v2.9.1-challenging-perf.md`. A pasteable tag note is
`docs/release-v2.9.1.md`. The package version stays 0.1.0.
`python -m codonforge.cli.main --version` prints 0.1.0.

Search on `usage.cai` no longer treats a zero-weight codon as a flat
landscape. When the incumbent has any counted weight of 0, tiers 1 and 2
linearize that term as the change in the sum of log-weights, and a zero
weight sits below every positive log (`DEP-SRCH-020`). The reported value
stays the geometric mean, which is still 0 while any counted weight is 0.
With no zero weight the linearization is the arithmetic delta, bit for bit.

On the CAI-only gate (budget `fast`, 4 generations, seed 7) `aceB` and
`thrA` now score 1.0, equal to `cai_greedy`, with no in-frame `AGG`. The
other 48 scores are unchanged. `GATE-BEN-001` is claimed on 50 of 50
proteins. The MFE comparison is still not this claim. A recorded score
below `cai_greedy` is still refused by `check_cai_gate`.

A structure job reuses an exact accessibility value for a repeated CDS and
for a sequence a few codons away (`AccessibilityCache.substitute`). The
phase 2 knee vectors are unchanged. On this machine, wall time for the
documented standard 4-generation CAI+structure runs fell from 1.357 s to
0.561 s on rpsL and from 4.211 s to 2.255 s on waaS.

`design_operon` updates the flux-weighted codon ledger one codon at a time.
The edited vector matches a full recount. On lacZ plus thrA (fluxes 1 and 2,
burden objective, budget `fast`, generations 0, two passes, seed 7) the joint
cost stays 0.4452346719259927 and the swaps stay `[26, 0]`. Wall time on this
machine fell from 9.185 s to 1.632 s.

The challenging-sequence timings are in
`docs/benchmarks/v2.9.1-challenging-perf.md`. No design in that sweep
returned an infeasible CDS. The slowest job is yeeJ under the assembly
rules, about 21 s before the population loop.

## v2.9.0

`POST /v1/observations` validates a wet-lab submission and appends it.
`GET /v1/observations` returns that tenant's log. The eight checks are
`CF-OBS-001` through `CF-OBS-008`. A failed check stores nothing.

The stored document's `ingest_status.state` is `validated`. Its digest is
recomputed after that field is written. A caller-supplied digest that
disagrees with the submitted document is `CF-OBS-001`.

The secondary model is ordinary least squares of one uncensored numeric
readout on the GC fraction of an explicit CDS. Two rows of
`protein_abundance`, CDS `ATGAAATAA` (value 1) and `ATGGGGGGG` (value 4),
fit with residual sum of squares 0. One row, a constant GC fraction, a
censored readout, or two different latent quantities stay `state: not_fitted`
and carry no coefficient. No codon weight is invented.

Recalibration is not scheduled. `GET /v1/recalibrations` stays `CF-CLI-003`.
`DEP-ACT-003`.

The package version stays 0.1.0.
`python -m codonforge.cli.main --version` prints 0.1.0.

## v2.8.0

`design_operon` designs each cistron, then accepts a synonymous swap only
when `burden.growth_cost` of the caller-supplied flux-weighted codon counts
falls and the sequence stays feasible. The protein is unchanged. Fluxes are
the caller's. No expression rate, gene-copy number, or s-value is invented.

On two copies of `ATGCTGTAA` (translation table 11, fluxes 1 and 1, leucine
relative adaptiveness CTG 1, TTA 0.8, TTG 0.6, CTT 0.4, CTC 0.2, CTA 0.1),
budget `fast`, generations 0, seed 7, the independent designs stay
`ATGCTGTAA` / `ATGCTGTAA` with joint growth cost 0.7968828039139021. One
joint pass moves the first cistron to `ATGTTATAA`. The joint cost is
0.6453291265499066. The solo costs of those two sequences are
0.8597495762920573 and 0.7968828039139021, and their sum is
1.6566323802059593. A second pass accepts no swap. This host table is the
fixture, not a measured expression rate. `DEP-SRCH-019`.

The package version stays 0.1.0.

## v2.7.0

`epistasis.context_naturalness` is an L1 objective: the mean log probability
of the selected synonymous Potts model. Higher is better. The unit is
dimensionless. The couplings stay in a process cache keyed by digest. A tier
request carries `context_digest`. The lockfile resolution entry is the
record, pinned, and the digest is the artifact digest.

The artifact is `cf:epistasis_potts:escherichia_coli:511145/MG1655:0.1.0`,
digest `sha256:ac64cae265510c2cad9c3a85a59b9abf389cae456eba25fc4595e21d20e645b6`,
lifecycle `draft`. The selected window is 10: held-out naturalness
−0.8708682912309792 nats/codon, codon-usage baseline −0.991738843659599,
2333 iterations, gradient norm 0.09770886731134135, converged. The window-5
fit is not selected (held-out −0.900526979129995, not converged). A job whose
host id starts with `escherichia_coli:` may pin the MG1655 artifact. No
distance penalty is applied, because none was measured. Any other host is
`CF-LAND-001`. A missing artifact is `CF-LAND-001` and is not scored as zero.
`epistasis.dwell_deviation` stays `CF-LAND-001`.

The interaction radius is `Codons(10)`. The job sets `tier2_block` to
`max(8, window)`. `P1-26` stays held. `DEP-LAND-005`.

The package version stays 0.1.0.
`python -m codonforge.cli.main --version` prints 0.1.0.

## v2.6.0

Structure and host burden are ranked objectives, and a job that asks for more
than one complete axis recommends the knee of that front.

`structure.local_accessibility` stays the windowed accessibility penalty
(span 40 nt, window 6 nt, 37°C, ViennaRNA, lower is better). Tier 3 scores it
once per generation. Tiers 1 and 2 do not, and say so with `CF-SRCH-010`. A
job that asks for the term, or for a `local_structure` rule, and has no
ViennaRNA installed is refused with `CF-STRU-003`. A CAI-only job still runs.
The measured value is the structure axis of the recommended vector. When
burden is also requested, that axis is the fold, not the burden number.

`local_structure` rules with a declared span are measured with ViennaRNA.
Those ids are no longer `CF-FEAS-008`. A missing span is `CF-STRU-004`.
`min_codon_frequency` is still `CF-FEAS-008`. The Rust checker still lists
`local_structure` as undetectable; the Python job is what drops the warning
after it has measured the rule.

`burden.growth_cost` is a demand-profile mismatch on the synonymous codon
shares. Supply is `trna.gene_copy_number` when that map is non-empty,
otherwise `codon_usage.table`. Neither is invented. The aggregation is the
power mean of order 4, then `g = u^2 / (1 + u^2)` with the tier-3 ceiling 1.
The unit is dimensionless. It is not a growth rate. A codon pile-up scores
worse than the host mix. `DEP-BURD-005`.

The recommendation rule for one objective is still
`feasible_then_score_then_digest`. For two or more complete axes it is
`knee_of_complete_front`. `candidates[0]` is that candidate. `--output explain`
and the page print the range of each axis over the complete feasible vectors.
`CF-ORCH-007` still compares the scalar best across stages.

On the elastin-like repeat (`elp_vpgvg.fa`, 306 nt) under
`constraints_assembly.json` and `examples/mg1655_ribosomal_usage.json`, seed
7, budget `fast`, generations 0, the feasible CDS has independent ribosomal
CAI 0.4380255552056712. The v2.5.0 feasible value on the same table was
0.3699530182864034. The same CDS is what budget `standard` and 4 generations
returns for CAI only. The protein is unchanged. Measured fronts for rpsL,
BBa_E0040_onestop, waaS, and the elastin-like repeat are in
`docs/benchmarks/v2.6.0-phase2.md`.

The package version stays 0.1.0.
`python -m codonforge.cli.main --version` prints 0.1.0.

## v2.5.1

A design with hard constraints no longer stops at tier 1. When the dynamic
program cannot carry a GC window or a scattered repeat, a deterministic seed
repair runs and later stages start from a sequence that satisfies the rules
it can check. The four validation genes (yagF, waaS, the elastin-like repeat,
and yeeJ) return feasible CDS under the assembly rule set. A job that still
cannot is exit 3, `CF-SRCH-006`, and the unchanged input is not printed as
the design. Tiers that ran are `done`; tiers that did not run are `skipped`.

`custom_regex` is enforced when the pattern compiles. `local_structure` and
`min_codon_frequency` are reported as `CF-FEAS-008` and are not treated as
satisfied.

The recommended candidate says what CAI 1.0 means on the loaded table, scores
the same CDS on the Kazusa E. coli genome-wide weights, and names a counted
codon whose weight is 0. Presets for BioBrick RFC[10], Golden Gate, and a
general synthesis rule set are on the page. `examples/rpsL_constraints.json`
is named as a set that checks nothing.

PrimerForge strips a FASTA header the same way on the page and the API,
accepts lowercase by uppercasing it, and names an illegal base with its
position. The 3′ GC clamp is labelled informational. Self-dimer and
heterodimer ΔG stay on SantaLucia 1998 Table 2 without dangling ends, which
is not primer3's number.

The page can trim a registry CDS at the first standard-code stop. Untrimmed
BBa_E0040 still exits 2 with `CF-SEQ-003` when that box is off.

`python -m codonforge.cli.main` no longer prints a `RuntimeWarning`.

## v2.5.0

The DNA Design Suite page is a guided form for people who are not
programmers. The API routes are the same ones as v2.4.2.

CodonForge asks for a coding sequence (paste or FASTA upload), a host, a
codon-usage table, optional constraints, and objectives, then runs. The
bundled `examples/mg1655_ribosomal_usage.json` is a button, and it loads on
its own when the page opens. The rpsL demo fills the CDS, that table, and
`examples/rpsL_constraints.json`. Budget, generations, seed, and the access
token start on the values of that demo (fast, 4, 7, and a local token) and
sit under Advanced. CAI is the objective already ticked.

While a job runs, each search stage keeps its own row. A finished job shows
the sequence, the scores taken from the front document, whether
`constraint_report.feasible` is true or false, and each warning with its
code, stage, and subject plus a sentence on what to do. UNKNOWN stays the
word UNKNOWN. FASTA, GenBank, and CSV are downloaded from that document.
`CF-*` refusals, including Phase 1 routes that answer `CF-CLI-003`, name a
next step instead of leaving the reader with the JSON alone.

PrimerForge is the other tab. A region can be dragged, typed as a start and
end, or set to the whole sequence. The same helix, tokens, and typefaces are
used on both tabs.

The page is started with `python -m codonforge.api`. It binds port 8000,
opens a browser, and blocks until Ctrl+C. README has the PowerShell command
as well as the macOS and Linux one.

The CAI gate claim is unchanged from v2.4.2.

## v2.4.2

`GATE-BEN-001` is claimed on 48 of 50 proteins. The claim is the CAI-only
front against `cai_greedy`, not the MFE comparison.

`aceB` (`NP_418438.1`) and `thrA` (`NP_414543.1`) were run again at budget
`fast`, 20 generations, stagnation window 21, seed 7. Both ran 20
generations. Both still score 0 where `cai_greedy` scores 1, and both
designed sequences still contain in-frame `AGG` (`aceB` at codons 80, 136,
321, 454, 497, 502 and 508; `thrA` at 243, 406, 421, 538, 680, 731, 755 and
792). `usage.cai` is identically 0 while an `AGG` remains, because that
codon has weight 0 in the ribosomal-protein corpus. The greedy baseline
wins by construction on that flat objective. The two proteins are
`degenerate_cases` in `benchmarks/corpus/gates.json` and are excluded from
the claim.

`check_cai_gate` refuses a claim that covers a protein whose recorded
CodonForge score is below `cai_greedy`. `tools/gate_controls.sh` writes an
all-50 claim into `gates.json` and requires that refusal.

The canonical status note is `docs/project-status.md`.

## v2.4.1

The ribosomal usage table lists all 64 codons. `AGG` and `TAG` were the two
absent keys. In the 54 ribosomal-protein CDS of `NC_000913.3`, `AGG` occurs
0 times and `TAG` occurs 0 times, so both weights are 0. The other 62 weights
are unchanged. The CDS gzip digest is unchanged. A test requires every codon
to be present; omitting the zero-count codons makes that test fail.

The CAI-only `GATE-BEN-001` comparison was run again with that table. The 50
scores are the same as v2.4.0. `claimed` stays false.

- `aceB` (`NP_418438.1`): `AGG` at codon 14. Native CAI 0, CodonForge 0,
  `cai_greedy` 1.
- `thrA` (`NP_414543.1`): `AGG` at codons 28 and 243. Native CAI 0,
  CodonForge 0, `cai_greedy` 1.

`usage.cai` is 0 while any counted codon has weight 0. The 4-generation
search did not remove those codons. `cai_greedy` replaces each codon with
the highest-weight synonym, so it never keeps `AGG`. This is a search
shortfall on a flat zero, not a missing table entry.

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
