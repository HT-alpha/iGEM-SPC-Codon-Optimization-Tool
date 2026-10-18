#!/usr/bin/env python3
"""Render the two Part D mockups from the real job documents.

Values are read from fixtures/ — a front and status emitted by a real
`design()` run of request_for(objectives=['usage.cai',
'structure.local_accessibility']). Nothing here invents a score, a
warning, or a candidate.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FRONT = json.loads((ROOT / "fixtures" / "pareto_front.json").read_text())
STATUS = json.loads((ROOT / "fixtures" / "job_status.json").read_text())

CAI, STRU = FRONT["objective_order"]
REC = FRONT["recommended"]
REC_ID = REC["candidate_id"]

scored = [c for c in FRONT["candidates"] if STRU not in c["unknown_terms"]]
# Minimize: the structure-optimal scored member is the lowest penalty.
SOPT = min(scored, key=lambda c: c["objective_vector"][1])
SOPT_ID = SOPT["candidate_id"]


def n(value: float | None, digits: int = 6) -> str:
    if value is None:
        return "UNKNOWN"
    return f"{value:.{digits}g}"


def cds_of(candidate: dict) -> str:
    return candidate["sequences"][0]["cds"]


HEAD = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<link rel="stylesheet" href="codonforge-ui.css">
</head>
<body>
<nav class="nav">
  <a href="pareto_front.html">Pareto front</a>
  <a href="warnings_status.html">Warnings / status</a>
  <span class="note">Design-phase mockup. Not a production component. Values from fixtures/ (a real design() document). Faces: IBM Plex Sans / Mono and Source Serif 4 when installed; Noto / system otherwise.</span>
</nav>
"""

FOOT = """
<p class="fixture">
Source: <code>planning/04-api/mockups/fixtures/pareto_front.json</code>
and <code>job_status.json</code>, emitted by
<code>design(request_for(objectives=['usage.cai','structure.local_accessibility']))</code>
on 2026-09-13. No values were edited.
</p>
</div>
</body>
</html>
"""


def scatter_svg() -> str:
    # Plot box. UNKNOWN lives in a hatched rail ABOVE the numeric y-axis,
    # not at the worst score.
    W, H = 720, 340
    left, right, top, bottom = 58, 700, 54, 300
    rail_top, rail_bot = 10, 46

    cais = [c["objective_vector"][0] for c in FRONT["candidates"]]
    strs = [c["objective_vector"][1] for c in scored]
    xmin, xmax = min(cais) - 0.04, max(cais) + 0.04
    ymin, ymax = min(strs) - 0.4, max(strs) + 0.4  # higher penalty is visually up

    def x(cai: float) -> float:
        return left + (cai - xmin) / (xmax - xmin) * (right - left)

    def y(acc: float) -> float:
        # Minimize: better (lower) is toward the bottom of the numeric pane.
        return top + (acc - ymin) / (ymax - ymin) * (bottom - top)

    marks = []
    for c in FRONT["candidates"]:
        cx = x(c["objective_vector"][0])
        unknown = STRU in c["unknown_terms"]
        cy = (rail_top + rail_bot) / 2 if unknown else y(c["objective_vector"][1])
        rec = c["candidate_id"] == REC_ID
        sopt = c["candidate_id"] == SOPT_ID
        if rec:
            marks.append(
                f'<rect x="{cx-6:.1f}" y="{cy-6:.1f}" width="12" height="12" '
                f'fill="#d5e2ee" stroke="#1a4a73" stroke-width="2"/>'
            )
        elif sopt:
            marks.append(
                f'<path d="M{cx:.1f},{cy-7:.1f} L{cx+7:.1f},{cy:.1f} '
                f'L{cx:.1f},{cy+7:.1f} L{cx-7:.1f},{cy:.1f} Z" '
                f'fill="#d7e4db" stroke="#3a5f4d" stroke-width="2"/>'
            )
        elif unknown:
            marks.append(
                f'<rect x="{cx-4:.1f}" y="{cy-4:.1f}" width="8" height="8" '
                f'fill="url(#hatch)" stroke="#5c574e" stroke-width="1"/>'
            )
        else:
            marks.append(
                f'<rect x="{cx-3.5:.1f}" y="{cy-3.5:.1f}" width="7" height="7" '
                f'fill="#1c1917" stroke="#1c1917"/>'
            )

    xticks = ""
    for tick in (0.40, 0.55, 0.70, 0.85):
        px = x(tick)
        xticks += (
            f'<line x1="{px:.1f}" y1="{bottom}" x2="{px:.1f}" y2="{bottom+4}" stroke="#1c1917"/>'
            f'<text x="{px:.1f}" y="{bottom+16}" text-anchor="middle" font-size="11" '
            f'font-family="IBM Plex Mono, monospace">{tick:.2f}</text>'
        )
    yticks = ""
    for tick in (2.6, 4.0, 5.2, 6.0):
        py = y(tick)
        yticks += (
            f'<line x1="{left-4}" y1="{py:.1f}" x2="{left}" y2="{py:.1f}" stroke="#1c1917"/>'
            f'<text x="{left-8}" y="{py+4:.1f}" text-anchor="end" font-size="11" '
            f'font-family="IBM Plex Mono, monospace">{tick:.1f}</text>'
        )

    return f'''
<svg class="plot" viewBox="0 0 {W} {H}" width="100%" role="img"
     aria-label="CAI versus local accessibility. UNKNOWN members sit in a hatched rail, not at a worst-case score.">
  <defs>
    <pattern id="hatch" patternUnits="userSpaceOnUse" width="8" height="8">
      <rect width="8" height="8" fill="#ddd8cc"/>
      <path d="M0 8 L8 0" stroke="#c9c3b6" stroke-width="2"/>
    </pattern>
  </defs>
  <rect x="0" y="{rail_top}" width="{W}" height="{rail_bot-rail_top}" fill="url(#hatch)" stroke="#5c574e"/>
  <text x="{left}" y="{rail_top+12}" font-size="11" font-family="IBM Plex Sans, sans-serif" fill="#5c574e">
    UNKNOWN on structure.local_accessibility — unscored, not penalized (tier 1/2 seeds)
  </text>
  <line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#1c1917"/>
  <line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#1c1917"/>
  <text x="{W/2}" y="{H-6}" text-anchor="middle" font-size="12" font-family="IBM Plex Sans, sans-serif">
    usage.cai (maximize)
  </text>
  <text x="14" y="{(top+bottom)/2}" transform="rotate(-90 14 {(top+bottom)/2})" text-anchor="middle"
        font-size="12" font-family="IBM Plex Sans, sans-serif">
    structure.local_accessibility (minimize)
  </text>
  {xticks}{yticks}
  {''.join(marks)}
</svg>
'''


def rows() -> str:
    out = []
    for c in FRONT["candidates"]:
        cai, acc = c["objective_vector"]
        unknown = STRU in c["unknown_terms"]
        cls = ""
        tag = ""
        if c["candidate_id"] == REC_ID:
            cls = "rec"
            tag = '<span class="chip recommended">recommended</span> '
        elif c["candidate_id"] == SOPT_ID:
            cls = "sopt"
            tag = '<span class="chip structure">structure-optimal among scored</span> '
        acc_cell = (
            f'<td class="unk">UNKNOWN</td>'
            if unknown
            else f'<td class="mono">{n(acc)}</td>'
        )
        feas = (
            '<span class="chip feasible">feasible</span>'
            if c["constraint_report"]["feasible"]
            else '<span class="chip infeasible">infeasible</span>'
        )
        role = c["objective_detail"].get(CAI, {}).get("role", "")
        out.append(
            f'<tr class="{cls}">'
            f'<td class="mono">{tag}{c["candidate_id"]}</td>'
            f'<td class="mono">{c["provenance"]["origin"]}</td>'
            f'<td>{feas}</td>'
            f'<td class="mono">{n(cai)}</td>'
            f'{acc_cell}'
            f'<td class="mono">{role}</td>'
            f'<td class="mono">{", ".join(c["unknown_terms"]) or "—"}</td>'
            f'</tr>'
        )
    return "\n".join(out)


def front_html() -> str:
    rec = next(c for c in FRONT["candidates"] if c["candidate_id"] == REC_ID)
    informed = ", ".join(f'<code>{t}</code>' for t in REC["informed_by"])
    excluded = ", ".join(f'<code>{t}</code>' for t in REC["not_informed_by"])
    return (
        HEAD.format(title="DNA Design Suite — Pareto front")
        + f"""
<header class="mast">
  <div class="k">Document</div>
  <div class="v">
    <h1>Pareto front</h1>
    <p class="sub">
      Job <code>{FRONT["job_id"]}</code> ·
      front <code>{FRONT["front_id"]}</code> ·
      lockfile <code>{FRONT["landscape_lockfile_digest"]}</code>
    </p>
  </div>
</header>
<div class="page">
  <p class="note">
    The front is a set. <code>recommended</code> names one member under
    <code>feasible_then_score_then_digest</code>, and says which axes that
    scalarization used. A different member can be better on an axis the
    recommendation did not consider. That distinction is the point of this view.
  </p>

  <div class="legend">
    <span class="chip feasible">feasible</span>
    <span class="chip infeasible">infeasible (none in this job)</span>
    <span class="chip unknown">UNKNOWN (unscored, not worse)</span>
    <span class="chip recommended">recommended</span>
    <span class="chip structure">structure-optimal among scored</span>
    <span class="chip diagnostic">diagnostic-only (none in this job)</span>
  </div>

  <div class="callout">
    <strong>Recommended</strong>
    <code>{REC_ID}</code>
    · rule <code>{REC["rule"]}</code>
    · differs_from_last_stage = <code>{str(REC["differs_from_last_stage"]).lower()}</code>
    <br>
    informed_by {informed}
    <br>
    not_informed_by {excluded}
    <br>
    Origin <code>{rec["provenance"]["origin"]}</code> ·
    {CAI} {n(rec["objective_vector"][0])} (maximize) ·
    {STRU} <span class="chip unknown">UNKNOWN</span>
  </div>

  <div class="callout structure">
    <strong>Structure-optimal among members scored on that axis</strong>
    <code>{SOPT_ID}</code>
    · origin <code>{SOPT["provenance"]["origin"]}</code>
    · {STRU} {n(SOPT["objective_vector"][1])} (minimize) ·
    {CAI} {n(SOPT["objective_vector"][0])}
    <br>
    This is not the recommended candidate. The recommendation was not informed
    by <code>{STRU}</code>. Treating UNKNOWN as a bad score would have ranked
    the recommended member as if it had been folded and found poor; it was
    never folded.
  </div>

  <h2>Axes</h2>
  {scatter_svg()}
  <p class="note">
    Hatched marks sit in a rail that is not on the numeric axis. They are not
    placed at the worst observed accessibility. That is the Jalview-gap /
    IGV-empty-track treatment: absence of a measurement is not a measurement
    of absence. Backend counterpart: <code>Option&lt;f64&gt;</code> is
    incomparable, not imputed.
  </p>

  <h2>Members ({FRONT["reduction"]["presented_count"]} of
      {FRONT["reduction"]["raw_nondominated_count"]},
      reduction={FRONT["reduction"]["method"]})</h2>
  <table class="data">
    <thead>
      <tr>
        <th>candidate_id</th>
        <th>origin</th>
        <th>feasible</th>
        <th>{CAI} ↑</th>
        <th>{STRU} ↓</th>
        <th>role</th>
        <th>unknown_terms</th>
      </tr>
    </thead>
    <tbody>
      {rows()}
    </tbody>
  </table>

  <h2>Recommended CDS</h2>
  <div class="seq">{cds_of(rec)}</div>
  <h2>Structure-optimal CDS (among scored)</h2>
  <div class="seq">{cds_of(SOPT)}</div>

  <p class="note">
    Feasibility is two-valued: every member here is <code>feasible: true</code>
    with <code>n_blocking_violations: 0</code>. UNKNOWN is a per-term list,
    not a third value of <code>constraint_report.feasible</code>.
  </p>
"""
        + FOOT
    )


def status_html() -> str:
    tiers = []
    for t in STATUS["progress"]["tiers"]:
        frac = t.get("fraction")
        gen = t.get("generation")
        mx = t.get("max_generations")
        gen_s = f"{gen}/{mx}" if gen is not None else "—"
        pct = 0 if frac is None else int(round(frac * 100))
        tiers.append(
            f'<div class="tier-track">'
            f'<span>{t["tier"]}</span>'
            f'<span>{t["state"]}</span>'
            f'<div class="bar"><span style="width:{pct}%"></span></div>'
            f'<span>{gen_s}</span>'
            f"</div>"
        )

    warns = []
    for w in STATUS["warnings"]:
        warns.append(
            f'<div class="warn-row">'
            f'<div class="code">{w["code"]}</div>'
            f'<div><span class="mono">{w.get("stage", "?")}</span></div>'
            f'<div><span class="mono">{w.get("subject", "?")}</span></div>'
            f'<div><strong>{w["title"]}</strong><br><span class="note">{w.get("detail","")}</span></div>'
            f"</div>"
        )

    objs = []
    for o in STATUS["resolved_objectives"]:
        objs.append(
            f'<tr>'
            f'<td class="mono">{o["term_id"]}</td>'
            f'<td class="mono">{o["role"]}</td>'
            f'<td class="mono">{o["direction"]}</td>'
            f'</tr>'
        )

    hv = STATUS["current_front_summary"]["hypervolume_history"]
    hv_s = ", ".join(f"{x:.4g}" for x in hv)

    return (
        HEAD.format(title="DNA Design Suite — warnings and status")
        + f"""
<header class="mast">
  <div class="k">Document</div>
  <div class="v">
    <h1>Job status</h1>
    <p class="sub">
      <code>{STATUS["job_id"]}</code>
      · <span class="chip feasible">{STATUS["state"]}</span>
      · started {STATUS["started_at"]}
      · finished {STATUS["finished_at"]}
    </p>
  </div>
</header>
<div class="page">
  <p class="note">
    Provenance is a row, not a tooltip. A warning without
    <code>(code, stage, subject)</code> is the same failure as a skipped
    tier reported as success. This job has no checkpoint URI because the
    run did not pass <code>--checkpoint-dir</code>; the capability exists
    in the job lifecycle and is not invented here.
  </p>

  <h2>Identity</h2>
  <div class="rail"><div class="k">lockfile</div><div class="v mono">{STATUS["landscape_lockfile_digest"]}</div></div>
  <div class="rail"><div class="k">budget_tier</div><div class="v mono">{FRONT["summary"]["budget_tier"]}</div></div>
  <div class="rail"><div class="k">termination</div><div class="v mono">{FRONT["summary"]["termination_reason"]}</div></div>
  <div class="rail"><div class="k">evaluations</div><div class="v mono">{FRONT["summary"]["n_evaluations"]}</div></div>
  <div class="rail"><div class="k">feasibility</div><div class="v">
    <span class="chip feasible">{FRONT["summary"]["n_feasible"]} feasible</span>
    of {FRONT["reduction"]["presented_count"]} presented.
    Two-valued: constraint checking is in-core. UNKNOWN terms are listed
    per candidate, not as a third feasibility state.
  </div></div>

  <h2>Tier plan</h2>
  <p class="note">Each tier is a track. A percentage cannot distinguish a job that ran two tiers from one that skipped one.</p>
  <div class="rail"><div class="k">tier / state / fill / gen</div><div class="v">{''.join(tiers)}</div></div>

  <h2>Resolved objectives</h2>
  <table class="data">
    <thead><tr><th>term_id</th><th>role</th><th>direction</th></tr></thead>
    <tbody>{''.join(objs)}</tbody>
  </table>
  <p class="note">
    <code>role</code> is the effective registration, not <code>max_role</code>.
    Both terms here are <code>pareto_axis</code>. Diagnostic-only is a real
    token in the design language; this job does not have one.
  </p>

  <h2>Warnings</h2>
  {''.join(warns)}

  <div class="callout warning">
    <strong>CF-ORCH-007</strong> is why the recommended candidate can come from
    tier 1 while the last stage was mutation. The Pareto-front view names that
    candidate and the axes that did (and did not) inform it.
  </div>

  <h2>Hypervolume history (from current_front_summary)</h2>
  <p class="mono">{hv_s}</p>
  <p class="note">
    <code>best_per_objective.usage.cai</code> in the status document is
    <code>{STATUS["current_front_summary"]["best_per_objective"]["usage.cai"]}</code>
    — negated, because the population loop minimizes. This view does not plot
    that number as CAI. The front document's <code>objective_vector</code> plus
    <code>objective_sense</code> is the presentation source.
  </p>

  <h2>Colour tokens not observed in this job</h2>
  <div class="legend">
    <span class="chip infeasible">infeasible</span>
    <span class="chip diagnostic">diagnostic-only</span>
    <span class="chip refusal">hard-refusal</span>
  </div>
  <p class="note">
    Shown as tokens, not as fabricated rows. A succeeded job does not get a
    fake problem block to exercise the refusal colour.
  </p>
"""
        + FOOT
    )


def main() -> None:
    (ROOT / "pareto_front.html").write_text(front_html())
    (ROOT / "warnings_status.html").write_text(status_html())
    print("wrote", ROOT / "pareto_front.html")
    print("wrote", ROOT / "warnings_status.html")


if __name__ == "__main__":
    main()
