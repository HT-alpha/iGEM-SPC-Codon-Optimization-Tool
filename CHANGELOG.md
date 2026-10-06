# Changelog

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
