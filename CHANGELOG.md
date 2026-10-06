# Changelog

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
