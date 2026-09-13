# CodonForge frontend — design-phase brief (Part D)

**Substitution, disclosed.** This environment has no frontend-design skill and
no Visualizer tool. That absence was confirmed twice independently in prior
rounds; it is treated as settled. This package is therefore:

- a written design-language brief and information architecture,
- two inspectable HTML/CSS mockups,
- values taken only from a real `design()` emission (the documents in
  [`mockups/fixtures/`](mockups/fixtures/)), not from invented scores.

It is not a production UI, not a component library, and not a substitute for
the tool that was requested and is not here.

---

## 1. Point of view

A codon-optimization result is a lab result. The user is comparing designs the
way they already compare traces: side by side, with the measurement's identity
visible, and with missing data *looking like missing data*. The visual
language is therefore taken from genome browsers, structure viewers, and
instrument software — tools whose job is not to look friendly but to keep a
measurement from being mistaken for a different measurement.

That is the same provenance-transparency standard the backend already holds:
`CF-ORCH-006` for a skipped tier, `unknown_terms` for an unscored axis,
`(code, stage, subject)` for a warning. The UI does not get to be vaguer than
the document it renders.

---

## 2. Design language

### Typography

| Role | Face | Why |
| --- | --- | --- |
| Page title | **Source Serif 4** | NCBI Bookshelf / journal HTML. A result document, not a marketing page. |
| UI copy | **IBM Plex Sans** | IBM Research / Carbon-adjacent scientific products. Distinct from Inter, system-ui, and the default Tailwind sans stack. |
| Sequences, term ids, warning codes, digests, candidate ids | **IBM Plex Mono** | Genome browsers (IGV, NCBI Graphical Sequence Viewer) put nucleotides and accessions in a mono face so a 1-character shift is visible. A warning `code` is an identifier, not a sentence. |

Weights are 400 and 600 only. No display italics, no light/thin.

Mockup HTML does not fetch those faces from Google Fonts. This environment's
headless Chrome hangs on `fonts.googleapis.com` TLS; the named stack still
resolves locally when the faces are installed, and otherwise to Noto / system
as already listed in `codonforge-ui.css`. The design language is the named
faces. The mockup files stay inspectable offline.

### Named colour system

Named after the *state*, not after a palette step. None of these are Tailwind
`green-600` / `red-500` / `slate-400`.

| Token | Hex / treatment | Means | Drawn from |
| --- | --- | --- | --- |
| `feasible` | `#2c5f4a` on `#d5e4db` | Constraint report is true. Two-valued. | Chlorophyll / plus-feature annotation in genome browsers, used here as a *state* colour (ChimeraX: categorical colours for categories, never a sequential ramp for unordered states). |
| `infeasible` | `#8f2d2d` on `#f0d4d0` | Constraint report is false. | Restriction-enzyme map convention (NEB / SnapGene site markup): a blocking site is red because it is a stop, not a low score. |
| `UNKNOWN` | 45° hatch `#c9c3b6` / `#ddd8cc`, ink `#5c574e` | This term was not scored. Not a penalty. | **IGV** empty track; **Jalview** gap (Waterhouse et al. 2009) — a gap is not a residue with score 0; **AlphaFold DB** missing-residue treatment vs pLDDT (Jumper et al. 2021 / Varadi et al. 2022): low confidence and *no coordinates* are different visual channels. |
| `diagnostic-only` | `#8d6b24` on `#efe3c4` | Computed and reported; not an axis. | UCSC Genome Browser custom-track gold (Kent et al. 2002): an overlay that is not the coordinate genome. |
| `hard-refusal` | `#5c1515` on `#f4d6d0` | The job stopped (`job_status.problem`). | Instrument dashboards (Agilent / Thermo run-abort banners): a stopped run is a full-width state, not a toast. |
| `warning-with-provenance` | left rule `#c4922a`, wash `#f5ead0` | A warning that still carries `(code, stage, subject)`. | Same instrument-software warning row; NCBI's problem pages put the code in mono next to the prose. |

`recommended` `#1a4a73` is NCBI header blue — identity of *the named answer*,
not a gradient and not purple. `structure-opt` `#3a5f4d` is a second categorical
mark, not a step on the same ramp.

`--ink` `#1c1917` and `--paper` `#f3efe4` follow Jalview / paper-canvas, not
white-on-slate SaaS.

### Spacing and layout rhythm

The scale is **5 / 10 / 20 / 40 px**, taken from IGV's compact track padding
(~4–6 px between tracks, doubled for section breaks), not Tailwind's 4/8/16.
Layout is a **left label rail** (13.5 rem) plus a full-bleed value column —
the UCSC / IGV track-label pattern — not a centered card. Section rules are
1 px `--rule` or 2 px `--ink`. **Border radius is 0.** Genome-browser tracks
and instrument tables are rectangular.

---

## 3. Information architecture

Surfaces, in the order a user actually moves through a job.

### 3.1 Submission and status

| User question | Source in the documents | Presentation |
| --- | --- | --- |
| What was asked? | `job_status.resolved_objectives[].term_id, role, direction` | Table. `role` is the effective registration (`pareto_axis` / `epsilon_constraint` / `diagnostic`), not `max_role`. |
| Which budget / which tiers? | `pareto_front.summary.budget_tier`; `job_status.progress.tiers[]` | One track per planned tier, each with `state` (`done` / `skipped` / `pending` / …). Not a single percentage. `CF-ORCH-006` is how a skipped planned tier is said out loud. |
| How far through the population loop? | `progress.tiers[tier3_moea].generation / max_generations` | Shown on that track only. |
| Can I resume? | Checkpoint path is a CLI/job input (`--checkpoint-dir`), not a field invented on `job_status` | The status mockup of this run has **no checkpoint URI** because the run did not pass one. The IA records the capability; the mockup does not fabricate a document the job did not emit. |
| Did it finish, and how? | `state`, `termination_reason`, `started_at` / `finished_at` | State chip (feasible-green for `succeeded`; refusal colours for `failed` + `problem`). |

### 3.2 Pareto front

| User question | Source | Presentation |
| --- | --- | --- |
| What is the set? | `candidates[]`, `reduction` | Table + scatter. `reduction.method` is shown (`none` in this run). |
| Which one does CodonForge name? | `recommended.candidate_id`, `rule`, `differs_from_last_stage` | Callout, not a silent sort. |
| Which axes informed that name? | `recommended.informed_by`, `recommended.not_informed_by` | Named lists, same discipline as `CF-ORCH-006`. In this run: informed_by `usage.cai`; not_informed_by `structure.local_accessibility`. |
| Is someone else better on an excluded axis? | Members with a real score on that axis | A second, explicitly labelled mark: *structure-optimal among scored*. Not a second `recommended`. |
| Is it feasible? | `constraint_report.feasible` (boolean) | Two-valued chip. Never a third state. |
| Was every axis scored? | `unknown_terms` ↔ `objective_vector` nulls | Hatch + the word `UNKNOWN`. See §4. |

`objective_detail[].role` is shown per cell so a diagnostic term cannot be
read as an axis.

### 3.3 Warnings

Every warning is a row: **code** (mono), **stage**, **subject**, **title**,
**detail**. Collapsing that triple into a sentence is how `DEP-ORCH-001`
would be undone at the last step. This run's real warnings are `CF-SRCH-010`
(structure is L0), two `CF-SRCH-001` (operators awaiting P3-05), and
`CF-ORCH-007` (recommended came from `tier1_dp`, last stage was `mutation`).

### 3.4 Feasibility, as built

Internal model: `violation` is two-valued; `UNKNOWN` is per-term
`Option<f64>`. The schema matches that (`unknown_terms`, not tri-state
`feasible`). The UI matches that: a candidate can be feasible *and* UNKNOWN
on structure. That is the recommended candidate in this run.

---

## 4. Mixed scored / UNKNOWN on the structure axis

**The fact, from this job.** Sixteen front members. Twelve carry
`unknown_terms: ["structure.local_accessibility"]` and `objective_vector[1]
= null` — they were seeded by tiers 1/2, which cannot score an L0 term.
Four were reached by the tier-3 batch path and carry a real penalty
(2.60 … 5.93). The recommended member is one of the twelve UNKNOWN. The
lowest (best) scored accessibility is `cand_f65bc87b4ba42abe` at
`2.59685`, which is **not** the recommended candidate, because
`recommended.not_informed_by` names `structure.local_accessibility`.

**Treatment.** UNKNOWN members are drawn in a **hatched rail that is not on
the numeric y-axis**. They are not placed at the worst observed score, not
sorted as worst, not coloured with the infeasible or low-end sequential
ramp. The cell in the table is the word `UNKNOWN` on the same hatch. A
caption states the cause already present in the documents: origin
`tier1_dp` / `mutation`, plus `CF-SRCH-010`.

**Why this and not the alternatives.**

| Alternative | Why not |
| --- | --- |
| Impute 0, or the axis maximum | Invents a fold. The backend refuses that (`CF-STRU-006`); the UI does not get to do it in pixels. |
| Sort UNKNOWN last / first | Ranks on a value that does not exist. Dominance already treats `None` as incomparable. |
| Grey sequential (looks like “low”) | AlphaFold's lesson: pLDDT 50 and *no coordinates* must not share a channel. |
| Hide UNKNOWN members | Silently drops the recommended candidate, which is UNKNOWN on structure. |

**Provenance standard.** The backend already says the recommendation was not
informed by structure (`not_informed_by`) and that those members were not
scored (`unknown_terms`). The visual is that statement, not a reinterpretation
of it as a penalty.

---

## 5. Non-goals

- No production React/Vue/Svelte components.
- No REST chrome beyond what `job_status.links` already contains (`self`).
- No fabricated candidates, scores, warnings, checkpoints, or refusals.
- No tri-state feasibility.
- No decision that mixed UNKNOWN is “worse”; that question is closed in §4.

---

## 6. Anti-pattern checklist

| Forbidden | What was used instead | Why |
| --- | --- | --- |
| Unmodified default Tailwind spacing / colour scales | 5/10/20/40 spacing; named state tokens above | Tailwind's `green-600` would paint UNKNOWN and infeasible with the same sequential language the backend just separated. |
| Centered card with shadow | Full-bleed label rail, 0 radius, no shadow | UCSC/IGV tracks; a card would frame a result as a product empty-state. |
| Generic sans with flat hierarchy | Source Serif 4 title, Plex Sans UI, Plex Mono identifiers | Sequence and `CF-*` codes have to be readable as codes. |
| Purple-to-blue gradients | Flat categorical fills; NCBI blue only as `recommended` | Gradients imply a continuous quantity. Role and UNKNOWN are not. |
| Emoji as icons | CSS chips, SVG squares/diamonds, hatch pattern | A warning code is the icon. |
| Rounded-everything | `border-radius: 0` | Instrument tables and genome tracks are rectangular. |

---

## 7. Mockups

| File | What it must show |
| --- | --- |
| [`mockups/pareto_front.html`](mockups/pareto_front.html) | Recommended vs structure-optimal; mixed UNKNOWN hatch; `informed_by` / `not_informed_by`; two-valued feasible; `role`. |
| [`mockups/warnings_status.html`](mockups/warnings_status.html) | Tier tracks; warning `(code, stage, subject)`; resolved objectives' `role`; no fabricated checkpoint or refusal row. |

Regenerate from the fixtures with `python3 planning/04-api/mockups/render.py`.
Do not hand-edit numbers in the HTML.
