/* Renders job_status and pareto_front documents. Numbers come from the API.
   Polls GET /v1/jobs/{id} while a job is running. Codon tints, when shown,
   use relative_adaptiveness from the submitted host profile. The front
   document's objective_vector, unknown_terms, and objective_sense are per
   term, not per codon, so this page does not invent a position mask.
   The DNA helix (helix.js) is told the job state from the polling code below
   and from nowhere else. */

import { toCsv, toFasta, toGenbank } from "/formats.js";
import { createHelix } from "/helix.js";
import { nucleotideIndex } from "/viewer.js";

const tokenInput = document.getElementById("token");
const requestInput = document.getElementById("request");
const submitButton = document.getElementById("submit");
const submitError = document.getElementById("submit-error");
const statusRoot = document.getElementById("status-root");
const frontRoot = document.getElementById("front-root");
const codonHelix = createHelix(document.getElementById("codon-helix"));
const primerHelix = createHelix(document.getElementById("primer-helix"));

const tierFraction = new Map();
const openStages = new Set();
let tierJobId = "";
let lastSubmitted = null;
const POLL_MS = 2000;
const READOUT_IDLE = "Coordinates are read here, under the plot, so a label never covers a point or a row.";

function escapeText(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;");
}

function escapeAttr(value) {
  return escapeText(value).replaceAll('"', "&quot;");
}

async function call(method, path, body) {
  const headers = { "Content-Type": "application/json" };
  const token = tokenInput.value.trim();
  if (token) headers.Authorization = `Bearer ${token}`;
  const response = await fetch(path, {
    method,
    headers,
    body: body === undefined ? undefined : JSON.stringify(body),
  });
  const raw = await response.text();
  let payload = null;
  if (raw) payload = JSON.parse(raw);
  return { status: response.status, payload };
}

function chip(kind, label) {
  return `<span class="chip ${kind}">${escapeText(label)}</span>`;
}

function roleChip(role) {
  if (role === "diagnostic") return chip("diagnostic", role);
  if (role === "pareto_axis") return chip("recommended", role);
  return `<span class="mono">${escapeText(role || "")}</span>`;
}

function codeLink(problem) {
  const code = (problem && problem.code) || "";
  if (!code) return "";
  const href = (problem && problem.type) || `https://dna_design_suite.org/problems/${String(code).toLowerCase()}`;
  return `<a href="${escapeAttr(href)}">${escapeText(code)}</a>`;
}

/* What to try, keyed by the registered code. The code, stage, and subject
   stay on the row; this sentence is the action, not a replacement for them. */
const NEXT_STEP = {
  "CF-AUTH-001": "The access token was missing or not in the expected shape. Open Advanced and use a token that starts with t_ followed by 8 to 32 letters or digits. The value already filled in works on this page.",
  "CF-AUTH-003": "That job was started with a different token. Use the same token, or run the design again.",
  "CF-CLI-001": "The request was incomplete or not in the shape this server accepts. Check that the coding sequence uses only A, C, G and T, that its length divides by 3, and that a codon-usage table is loaded.",
  "CF-SEQ-013": "The host name is not in the shape genus_species:taxid/strain. The rpsL demo uses escherichia_coli:511145/MG1655. The codon table is chosen separately from that name.",
  "CF-CLI-002": "The job finished without a candidate. A constraint that no synonymous sequence can satisfy is a common cause. Remove a rule and run again.",
  "CF-CLI-003": "This version does not provide that feature. Codon design accepts DNA, a protein, or one GenBank CDS. Primer design, the other tab, picks PCR primers. Dwell is not scored.",
  "CF-IDEM-001": "The same submit was sent twice with different contents. Reload the page and run the design once.",
  "CF-LAND-001": "That objective is not available. Dwell is not scored: no A-site offset for GSE53767 is stated, so this version does not invent one. Choose CAI, whole-gene GC fraction, constraint violations, or local mRNA structure.",
  "CF-LAND-004": "That term is reported, not ranked. Leave it unchecked, or keep CAI as the objective.",
  "CF-FEAS-001": "A rule is broken in a part of the sequence the search is not allowed to change. Remove that rule, or clear the fixed flanks under Advanced.",
  "CF-FEAS-002": "A repeat rule cannot be satisfied for a sequence of this length. Raise the allowed count, or drop the rule.",
  "CF-FEAS-003": "At one codon, every synonym is forbidden. Relax a motif rule and run again.",
  "CF-FEAS-004": "Taken together, the motif rules leave no way through the gene. Remove one of them.",
  "CF-FEAS-005": "No synonymous coding sequence can land inside that GC window. Widen the minimum and maximum, or lengthen the window.",
  "CF-FEAS-006": "A constraint was not understood. Use the tick boxes, or paste a file that matches the example.",
  "CF-SRCH-001": "One edit operator is not available in this version, so the search used the others. The design is still real. This is a notice, not a failed job.",
  "CF-SRCH-006": "No synonymous sequence satisfied the hard constraints. The detail names the rules that remain. Widen a GC window, raise a repeat allowance, or drop a site, motif, or pattern the codons cannot avoid.",
  "CF-SRCH-017": "The first pass cannot carry a GC window or a scattered repeat. A repair step built a sequence that satisfies the rules it can check, and the later stages started from that sequence.",
  "CF-FEAS-008": "That rule is not searched in this version, so the design does not treat it as satisfied. local_structure needs the Phase 2 structure engine. min_codon_frequency is not enforced.",
  "CF-SRCH-010": "Local mRNA structure cannot be updated inside one of the search stages, so that stage did not use it. Where the score was not computed, the result says UNKNOWN. UNKNOWN is not zero.",
  "CF-STRU-003": "No RNA folding program is installed on the computer running this server. DNA Design Suite does not include one. Install ViennaRNA for the thermodynamic score. The Nussinov option on the form is an approximate paired fraction, and it is used only when you select it. CAI does not need either.",
  "CF-ORCH-002": "The generation budget ran out before every stage finished. The sequence is from the stages that did finish. Raise generations under Advanced if you want the population stage to run longer.",
  "CF-ORCH-006": "The budget you chose names a stage this version does not run. The stages that did run are listed, and the sequence is from those. Fast and Standard do not ask for the extra stage. Choosing Exact asks for a final ranking that is reported as skipped.",
  "CF-ORCH-007": "A later stage did not beat an earlier one. The sequence shown is the better of the two, and the notice names both stages. The search is allowed to try a change that does not win.",
  "CF-PRIM-001": "No primer pair passed the filters. Select a longer region, widen the product-size range, or move the melting-temperature target.",
  "CF-PRIM-002": "There is not enough sequence beside the region for a primer. Select a region further from the end, or paste a longer template.",
  "CF-PRIM-003": "The selected region does not sit inside the sequence. Select again, or use Select the whole sequence.",
  "CF-PRIM-006": "A letter in the template is not A, C, G or T, or a minimum is larger than its maximum. The detail names the base and its position in the header-stripped sequence. Lowercase bases are accepted.",
  "CF-PRIM-007": "Many primers were close to the melting-temperature target, so only the closest ones were checked for dimers. The pairs shown are those. This is a notice, not a failed search.",
};

const TIER_PLAIN = {
  tier1_seed: "First pass. Picks a codon at each position.",
  tier2_block_ca: "Polish. Revisits the gene in blocks.",
  tier3_moea: "A population of alternatives, for the number of generations you set.",
  decision: "A final ranking this version does not run.",
  primer_generate: "Builds candidate primers on each side of the region.",
  primer_prune: "Drops primers whose melting temperature is far from the target.",
  primer_evaluate: "Checks hairpins and dimers.",
  primer_rank: "Sorts the pairs that remain.",
};

const STATE_PLAIN = {
  queued: "Waiting to start. This page checks again every 2 seconds.",
  running: "Working. Each row below is one stage. There is no single percentage, because a skipped stage is not part of a bar.",
  succeeded: "Finished. The sequence, the scores, and the constraint check are under Pareto front.",
  failed: "Stopped before a design was ready. The code names the reason. The sentence under it says what to try.",
  cancelled: "Cancelled. No finished design was requested.",
};

const TERM_PLAIN = {
  "usage.cai": "CAI",
  "usage.harmonize_rca": "Harmonized RCA",
  "usage.match_frequency": "Codon-frequency match",
  "usage.gc_content": "GC fraction",
  "constraint.violation_count": "Constraint violations",
  "structure.local_accessibility": "Local mRNA structure",
  "burden.growth_cost": "Codon-usage burden",
};

function tradeoffNote(front) {
  const order = front.objective_order || [];
  const sense = front.objective_sense || [];
  if (order.length < 2 || sense.length !== order.length) return "";
  const complete = (front.candidates || []).filter((candidate) => {
    const vector = candidate.objective_vector || [];
    return vector.length === order.length && vector.every((value) => value !== null && value !== undefined);
  });
  if (complete.length < 2) return "";
  const parts = order.map((term, index) => {
    const values = complete.map((candidate) => Number(candidate.objective_vector[index]));
    const low = Math.min(...values);
    const high = Math.max(...values);
    return { term, low, high, sense: sense[index] };
  });
  if (parts.every((part) => part.low === part.high)) return "";
  const text = parts
    .map((part) => `${part.term} ${part.low.toFixed(6)} .. ${part.high.toFixed(6)} (${part.sense})`)
    .join("; ");
  return `<p class="note"><strong>Trade-off on this front.</strong> ${complete.length} candidates with complete vectors. ${escapeText(text)}.</p>`;
}

function nextStep(problem) {
  const text = problem && NEXT_STEP[problem.code];
  if (!text) return "";
  return `<p class="next-step"><strong>What to do.</strong> ${escapeText(text)}</p>`;
}

function problemCallout(problem) {
  if (!problem || !problem.code) {
    return `<div class="callout refusal"><strong>No CF-* code.</strong> <span class="note">The response did not name a registered code.</span></div>`;
  }
  const human = problem.detail || problem.title || "";
  return `<div class="callout refusal"><p class="human">${escapeText(human)}</p><p class="note"><strong>${codeLink(problem)}</strong> ${escapeText(problem.title || "")}</p>${nextStep(problem)}</div>`;
}

function lateReveal(html) {
  return `<div class="reveal reveal-late">${html}</div>`;
}

function isRefusal(item) {
  return typeof item.status === "number" && item.status >= 400;
}

function captureOpenStages() {
  openStages.clear();
  statusRoot.querySelectorAll(".warn-group.open").forEach((group) => {
    if (group.dataset.stage) openStages.add(group.dataset.stage);
  });
}

function renderWarnings(document) {
  const items = [...(document.warnings || [])];
  if (document.problem && document.problem.code) items.push(document.problem);
  if (!items.length) return '<p class="note">No warnings on this status.</p>';
  const groups = new Map();
  for (const item of items) {
    const stage = item.stage || "unspecified";
    if (!groups.has(stage)) groups.set(stage, []);
    groups.get(stage).push(item);
  }
  return [...groups.entries()]
    .map(([stage, group]) => {
      const refusal = group.some(isRefusal);
      const open = openStages.has(stage);
      const codes = [...new Set(group.map((item) => item.code).filter(Boolean))].join(", ");
      const firstDetail = (group.find((item) => item.detail) || {}).detail || "";
      const rows = group
        .map(
          (warning) =>
            `<div class="warn-row"><div class="code">${codeLink(warning)}</div><div class="mono">${escapeText(warning.stage || "")}</div><div class="mono">${escapeText(warning.subject || "")}</div><div><strong>${escapeText(warning.title || "")}</strong><div class="note">${escapeText(warning.detail || "")}</div>${nextStep(warning)}</div></div>`,
        )
        .join("");
      const kind = refusal ? "refusal" : "notice";
      const label = refusal ? "hard refusal" : "notice";
      return `<section class="warn-group ${kind}${open ? " open" : ""}" data-stage="${escapeAttr(stage)}">
        <button type="button" class="warn-toggle" aria-expanded="${open ? "true" : "false"}">${escapeText(stage)} · ${group.length} · ${escapeText(label)}${codes ? ` · ${escapeText(codes)}` : ""}${firstDetail ? ` — ${escapeText(firstDetail)}` : ""}</button>
        <div class="warn-body">
          <div class="warn-row head"><div>code</div><div>stage</div><div>subject</div><div>title</div></div>
          ${rows}
        </div>
      </section>`;
    })
    .join("");
}

function bindWarnings() {
  statusRoot.querySelectorAll(".warn-group").forEach((group) => {
    const body = group.querySelector(".warn-body");
    const button = group.querySelector("button");
    if (!body || !button) return;
    const apply = (open, animate) => {
      if (!animate) body.style.transition = "none";
      body.style.maxHeight = open ? `${body.scrollHeight}px` : "0px";
      if (!animate) {
        body.offsetHeight;
        body.style.transition = "";
      }
    };
    apply(group.classList.contains("open"), false);
    button.addEventListener("click", () => {
      const open = !group.classList.contains("open");
      group.classList.toggle("open", open);
      button.setAttribute("aria-expanded", open ? "true" : "false");
      if (open) openStages.add(group.dataset.stage);
      else openStages.delete(group.dataset.stage);
      apply(open, true);
    });
  });
}

function renderStatus(document) {
  lastJobStatus = document;
  captureOpenStages();
  if (tierJobId !== (document.job_id || "")) {
    tierFraction.clear();
    tierJobId = document.job_id || "";
  }
  const progress = document.progress || { tiers: [] };
  const tiers = (progress.tiers || [])
    .map((tier) => {
      const fraction = typeof tier.fraction === "number" ? tier.fraction : null;
      const width = fraction === null ? 0 : Math.max(0, Math.min(1, fraction)) * 100;
      const start = tierFraction.has(tier.tier) ? tierFraction.get(tier.tier) : width;
      tierFraction.set(tier.tier, width);
      const generation =
        tier.generation === undefined || tier.generation === null
          ? ""
          : `${tier.generation}/${tier.max_generations ?? "—"}`;
      const gloss = TIER_PLAIN[tier.tier] ? `<p class="tier-plain">${escapeText(TIER_PLAIN[tier.tier])}</p>` : "";
      return `<div class="tier-block"><div class="tier-track"><span>${escapeText(tier.tier)}</span><span>${escapeText(tier.state)}</span><div class="bar"><span style="width:${start}%" data-target="${width}"></span></div><span>${escapeText(generation)}</span></div>${gloss}</div>`;
    })
    .join("");
  const objectives = (document.resolved_objectives || [])
    .map(
      (term) =>
        `<tr><td class="mono">${escapeText(term.term_id)}</td><td>${roleChip(term.role)}</td><td class="mono">${escapeText(term.direction)}</td></tr>`,
    )
    .join("");
  const failed = document.state === "failed" ? lateReveal(problemCallout(document.problem)) : document.problem && document.problem.code ? problemCallout(document.problem) : "";
  const stateKind = document.state === "succeeded" ? "feasible" : document.state === "failed" ? "refusal" : "warning";
  const plainState = STATE_PLAIN[document.state] || "";
  const running = document.state === "queued" || document.state === "running";
  const stop = running
    ? `<button type="button" id="cancel-job" class="secondary">Stop this job</button>`
    : "";
  statusRoot.innerHTML = `
    <div class="rail"><div class="k">Job</div><div class="v mono">${escapeText(document.job_id || "")} ${stop}</div></div>
    <div class="rail"><div class="k">State</div><div class="v">${chip(stateKind, document.state || "")}${plainState ? `<p class="field-help">${escapeText(plainState)}</p>` : ""}</div></div>
    ${failed}
    <h2>Tiers</h2>
    ${tiers || '<p class="note">No tier has been reported.</p>'}
    <h2>Resolved objectives</h2>
    <div class="hscroll"><table class="data"><thead><tr><th>term</th><th>role</th><th>direction</th></tr></thead><tbody>${objectives}</tbody></table></div>
    <h2>Warnings</h2>
    ${renderWarnings(document)}
  `;
  bindWarnings();
  const cancelButton = statusRoot.querySelector("#cancel-job");
  if (cancelButton && document.job_id) {
    cancelButton.addEventListener("click", () => {
      cancelJob(document.job_id);
    });
  }
  requestAnimationFrame(() => {
    statusRoot.querySelectorAll(".bar > span").forEach((span) => {
      span.style.width = `${span.dataset.target}%`;
    });
  });
}

function termValue(candidate, term, order) {
  const detail = (candidate.objective_detail || {})[term];
  if (detail && detail.value !== undefined && detail.value !== null) return detail.value;
  const index = order.indexOf(term);
  if (index < 0) return null;
  const vector = candidate.objective_vector || [];
  return vector[index] === undefined ? null : vector[index];
}

function isStructureTerm(term) {
  return String(term).startsWith("structure.");
}

function structureOptimal(front) {
  const recommended = front.recommended || {};
  const excluded = recommended.not_informed_by || [];
  const order = front.objective_order || [];
  const senses = front.objective_sense || [];
  const found = [];
  for (const term of excluded) {
    const index = order.indexOf(term);
    const sense = index >= 0 ? senses[index] : "minimize";
    let best = null;
    for (const candidate of front.candidates || []) {
      const value = termValue(candidate, term, order);
      if (value === null || Number.isNaN(Number(value))) continue;
      const numeric = Number(value);
      if (
        best === null ||
        (sense === "minimize" && numeric < best.value) ||
        (sense === "maximize" && numeric > best.value)
      ) {
        best = { candidate, value: numeric, term, sense };
      }
    }
    if (best) found.push(best);
  }
  return { excluded, found };
}

function optimalCallout(item) {
  const structure = isStructureTerm(item.term);
  const title = structure
    ? "Structure-optimal among members scored on an axis the recommendation did not use"
    : "Best scored member on an axis the recommendation did not use";
  const klass = structure ? "callout structure pulse-once" : "callout";
  return `<div class="${klass}"><strong>${escapeText(title)}</strong> <code>${escapeText(item.candidate.candidate_id)}</code> · <code>${escapeText(item.term)}</code> ${escapeText(item.value)} (${escapeText(item.sense)}). This is not a second recommended candidate.</div>`;
}

function roleFor(front, term) {
  for (const candidate of front.candidates || []) {
    const detail = (candidate.objective_detail || {})[term];
    if (detail && detail.role) return detail.role;
  }
  return "";
}

function feasibleCell(candidate) {
  const report = candidate.constraint_report;
  if (!report || typeof report.feasible !== "boolean") {
    return `<td class="note">constraint_report.feasible absent</td>`;
  }
  return `<td>${report.feasible ? chip("feasible", "true") : chip("infeasible", "false")}</td>`;
}

function codonWeights(request) {
  const table = request && request.host_profile && request.host_profile.codon_usage && request.host_profile.codon_usage.table;
  if (!table || typeof table !== "object") return null;
  const weights = {};
  for (const [codon, entry] of Object.entries(table)) {
    if (entry && typeof entry.relative_adaptiveness === "number") weights[codon] = entry.relative_adaptiveness;
  }
  return Object.keys(weights).length ? weights : null;
}

function activeWeights() {
  if (lastSubmitted) return codonWeights(lastSubmitted);
  try {
    return codonWeights(JSON.parse(requestInput.value));
  } catch (error) {
    return null;
  }
}

function topQuartileThreshold(weights) {
  const values = Object.values(weights).filter((value) => typeof value === "number").sort((a, b) => a - b);
  if (!values.length) return null;
  const index = Math.min(values.length - 1, Math.ceil(0.75 * values.length) - 1);
  return values[index];
}

function caiUnknown(candidate, front) {
  if ((candidate.unknown_terms || []).includes("usage.cai")) return true;
  const index = (front.objective_order || []).indexOf("usage.cai");
  if (index < 0) return false;
  const value = (candidate.objective_vector || [])[index];
  return value === null || value === undefined;
}

function colourCodons(cds, weights, threshold, withhold) {
  if (cds.length % 3 !== 0) return escapeText(cds);
  let html = "";
  for (let index = 0; index + 3 <= cds.length; index += 3) {
    const codon = cds.slice(index, index + 3);
    const weight = weights ? weights[codon] : undefined;
    const high = !withhold && weights && threshold !== null && typeof weight === "number" && weight >= threshold;
    html += `<span class="codon${high ? " top" : ""}">${escapeText(codon)}</span>`;
  }
  return html;
}

function renderSequences(front, recommendedId, optimalIds) {
  const weights = activeWeights();
  const threshold = weights ? topQuartileThreshold(weights) : null;
  const wanted = [];
  for (const candidate of front.candidates || []) {
    const id = candidate.candidate_id;
    if (id === recommendedId) wanted.push({ candidate, role: "recommended" });
    else if (optimalIds.has(id)) wanted.push({ candidate, role: "structure-optimal" });
  }
  if (!wanted.length) return "";
  return wanted
    .map(({ candidate, role }) => {
      const withhold = caiUnknown(candidate, front);
      const unknown = (candidate.unknown_terms || []).join(", ");
      const blocks = (candidate.sequences || [])
        .map((sequence) => {
          const cds = sequence.cds || "";
          const note = withhold
            ? "usage.cai is UNKNOWN for this candidate, so codon tints are withheld."
            : weights
              ? `Top-quartile relative adaptiveness from the submitted codon table (threshold ${threshold}). The front document has no per-codon score.`
              : "No codon tint. objective_vector, unknown_terms, and objective_sense do not name codons, and this page has no submitted codon table.";
          return `<div class="seq-actions"><span class="mono">${escapeText(sequence.cistron_id || "cds")}</span><button type="button" data-copy="${escapeAttr(cds)}">Copy raw sequence</button></div>
            <div class="seq-block">${colourCodons(cds, weights, threshold, withhold)}</div>
            <p class="note">${escapeText(note)}</p>`;
        })
        .join("");
      if (role === "recommended" && (candidate.sequences || [])[0]) {
        window.optimizedCds = candidate.sequences[0].cds;
        refreshOptimizedHandoff();
      }
      const unknownLine = unknown
        ? `<p class="note"><span class="unk">UNKNOWN</span> ${escapeText(unknown)}. These are unscored terms, not codon positions. This document has no position mask, so codons are not hatched.</p>`
        : "";
      return `<h2>Sequence · ${escapeText(role)}</h2><p class="mono">${escapeText(candidate.candidate_id)}</p>${unknownLine}${blocks}`;
    })
    .join("");
}

function readoutFor(candidate, order) {
  const parts = order.map((term, index) => {
    const value = (candidate.objective_vector || [])[index];
    const shown = value === null || value === undefined ? "UNKNOWN" : value;
    return `${term} ${shown}`;
  });
  return `${candidate.candidate_id} · ${parts.join(" · ")}`;
}

function renderPlot(front, recommendedId, optimalIds) {
  const order = front.objective_order || [];
  if (order.length < 2) {
    return `<p class="note">This front has fewer than two axes, so there is no scatter. UNKNOWN terms stay in the table.</p>`;
  }
  const senseY = (front.objective_sense || [])[1];
  const points = [];
  const unknown = [];
  for (const candidate of front.candidates || []) {
    const x = (candidate.objective_vector || [])[0];
    const y = (candidate.objective_vector || [])[1];
    const xKnown = !(x === null || x === undefined);
    const yKnown = !(y === null || y === undefined);
    if (xKnown && yKnown) points.push({ candidate, x: Number(x), y: Number(y) });
    else unknown.push({ candidate, x: xKnown ? Number(x) : null, y: yKnown ? Number(y) : null });
  }
  const xs = points.map((point) => point.x).concat(unknown.filter((point) => point.x !== null).map((point) => point.x));
  if (!xs.length && !unknown.length) return "";
  const minX = xs.length ? Math.min(...xs) : 0;
  const maxX = xs.length ? Math.max(...xs) : 1;
  const ys = points.map((point) => point.y);
  const minY = ys.length ? Math.min(...ys) : 0;
  const maxY = ys.length ? Math.max(...ys) : 1;
  const scaleX = (value) => 40 + ((value - minX) / (maxX - minX || 1)) * 640;
  const scaleY = (value) => 20 + (1 - (value - minY) / (maxY - minY || 1)) * 220;
  const mark = (id, x, y) => {
    const fill = id === recommendedId ? "#1a4a73" : optimalIds.has(id) ? "#3a5f4d" : "#1c1917";
    return `<rect data-candidate="${escapeAttr(id)}" x="${x - 8}" y="${y - 8}" width="16" height="16" fill="transparent"/><rect data-candidate="${escapeAttr(id)}" x="${x - 4}" y="${y - 4}" width="8" height="8" fill="${fill}"/>`;
  };
  const marks = points.map((point) => mark(point.candidate.candidate_id, scaleX(point.x), scaleY(point.y))).join("");
  const rail = unknown.length
    ? `<rect x="20" y="250" width="680" height="40" fill="url(#hatch)"/>${unknown
        .map((point, index) => {
          const xPos = point.x === null ? 40 + index * 16 : scaleX(point.x);
          return mark(point.candidate.candidate_id, xPos, 270);
        })
        .join("")}`
    : "";
  const missingX = unknown.filter((point) => point.x === null).length;
  const caption = unknown.length
    ? `<p class="note">${unknown.length} member(s) are UNKNOWN on at least one of <code>${escapeText(order[0])}</code> and <code>${escapeText(order[1])}</code> and sit in the hatched rail. They are not given a number on the ${escapeText(senseY || "numeric")} axis.${missingX ? " A member whose first axis is UNKNOWN is spaced along the rail in member order, which is not a score." : ""}</p>`
    : `<p class="note">Every plotted member has a number on both axes. An UNKNOWN value would sit in the hatched rail, not on this scale.</p>`;
  return `<svg class="plot" viewBox="0 0 720 300" width="100%" role="img" aria-label="First axis versus second. UNKNOWN members sit in a hatched rail, not at a worst-case score.">
    <defs><pattern id="hatch" patternUnits="userSpaceOnUse" width="8" height="8"><rect width="8" height="8" fill="#ddd8cc"/><path d="M0 8 L8 0" stroke="#c9c3b6" stroke-width="2"/></pattern></defs>
    ${marks}${rail}
  </svg>${caption}<p id="plot-readout" class="note plot-readout">${escapeText(READOUT_IDLE)}</p>`;
}

function renderFront(front) {
  const candidates = front.candidates || [];
  if (!candidates.length) {
    const reason = (front.summary && front.summary.termination_reason) || "";
    const match = String(reason).match(/CF-[A-Z]{3,8}-[0-9]{3}/);
    const problem = match
      ? { code: match[0], type: `https://dna_design_suite.org/problems/${match[0].toLowerCase()}`, title: "The front has no candidates", detail: reason }
      : null;
    frontRoot.innerHTML = lateReveal(problem ? problemCallout(problem) : `<div class="callout refusal"><strong>Empty front.</strong> <span class="note">No candidates were returned${reason ? `: ${escapeText(reason)}` : ""}. The document did not name a CF-* code for this outcome.</span></div>`);
    return false;
  }
  const recommended = front.recommended || {};
  const order = front.objective_order || [];
  const optimal = structureOptimal(front);
  const optimalIds = new Set(
    optimal.found
      .filter((item) => isStructureTerm(item.term))
      .map((item) => item.candidate.candidate_id),
  );
  const informed = (recommended.informed_by || []).map((term) => `<code>${escapeText(term)}</code>`).join(" ");
  const excluded = (recommended.not_informed_by || []).map((term) => `<code>${escapeText(term)}</code>`).join(" ");
  const optimalText = optimal.found.length
    ? optimal.found.map((item) => optimalCallout(item)).join("")
    : `<div class="callout"><strong>No scored member on an excluded axis.</strong> ${excluded || "not_informed_by is empty."} UNKNOWN is not placed on the numeric axis.</div>`;
  const rows = candidates
    .map((candidate) => {
      const cells = order
        .map((term, index) => {
          const value = (candidate.objective_vector || [])[index];
          if (value === null || value === undefined) return `<td class="unk">UNKNOWN</td>`;
          const detail = (candidate.objective_detail || {})[term];
          const role = detail && detail.role ? ` ${roleChip(detail.role)}` : "";
          return `<td class="mono">${escapeText(value)}${role}</td>`;
        })
        .join("");
      const extra = Object.entries(candidate.objective_detail || {})
        .filter(([term]) => !order.includes(term))
        .map(([term, detail]) => {
          const role = detail.role || "";
          const shown = detail.value === null || detail.value === undefined ? "UNKNOWN" : detail.value;
          const shownClass = shown === "UNKNOWN" ? "unk" : "mono";
          return `<div>${escapeText(term)} ${roleChip(role)} <span class="${shownClass}">${escapeText(shown)}</span></div>`;
        })
        .join("");
      const klass = candidate.candidate_id === recommended.candidate_id ? "rec" : optimalIds.has(candidate.candidate_id) ? "sopt" : "";
      const unknown = (candidate.unknown_terms || []).join(", ") || "—";
      const readout = escapeAttr(readoutFor(candidate, order));
      return `<tr class="${klass}" data-candidate="${escapeAttr(candidate.candidate_id)}" data-readout="${readout}"><td class="mono">${escapeText(candidate.candidate_id)}</td>${feasibleCell(candidate)}<td class="mono">${escapeText(unknown)}</td>${cells}<td>${extra}</td></tr>`;
    })
    .join("");
  const plot = renderPlot(front, recommended.candidate_id, optimalIds);
  const sequences = renderSequences(front, recommended.candidate_id, optimalIds);
  window.lastFront = front;
  frontRoot.innerHTML = `
    <div class="reveal">
    ${resultSummary(front)}
    <p class="note">reduction.method <code>${escapeText((front.reduction || {}).method || "")}</code> · presented ${(front.reduction || {}).presented_count ?? ""} of ${(front.reduction || {}).raw_nondominated_count ?? ""}</p>
    <div class="legend">${chip("feasible", "feasible")}${chip("infeasible", "infeasible")}${chip("unknown", "UNKNOWN")}${chip("recommended", "recommended")}${chip("structure", "structure-optimal among scored")}${chip("diagnostic", "diagnostic-only")}</div>
    <div class="callout pulse-once"><strong>Recommended</strong> <code>${escapeText(recommended.candidate_id || "")}</code> · rule <code>${escapeText(recommended.rule || "")}</code><br>informed_by ${informed || "—"}<br>not_informed_by ${excluded || "—"}</div>
    ${tradeoffNote(front)}
    ${optimalText}
    ${sequences}
    ${plot}
    <p class="note">A column in this table is a Pareto axis. A diagnostic term is in the last column and is not an axis. UNKNOWN is the word in the hatched cell, not a number.</p>
    <div class="hscroll"><table class="data"><thead><tr><th>candidate</th><th>constraint_report.feasible</th><th>unknown_terms</th>${order
      .map((term) => {
        const role = roleFor(front, term);
        return `<th>${escapeText(term)}${role ? ` ${roleChip(role)}` : ""}</th>`;
      })
      .join("")}<th>other terms</th></tr></thead><tbody>${rows}</tbody></table></div>
    </div>
  `;
  return true;
}

function resultSummary(front) {
  const recommended = front.recommended || {};
  const candidates = front.candidates || [];
  const candidate =
    candidates.find((item) => item.candidate_id === recommended.candidate_id) || candidates[0];
  if (!candidate) return "";
  const order = front.objective_order || [];
  const sense = front.objective_sense || [];
  const scores = order
    .map((term, index) => {
      const value = (candidate.objective_vector || [])[index];
      const label = TERM_PLAIN[term] || term;
      if (value === null || value === undefined) {
        return `<div class="rail"><div class="k">${escapeText(label)} <span class="mono">${escapeText(term)}</span></div><div class="v"><span class="unk">UNKNOWN</span> <p class="field-help">Not scored. This is not zero, and it is not a bad result.</p></div></div>`;
      }
      const direction = sense[index] || "";
      return `<div class="rail"><div class="k">${escapeText(label)} <span class="mono">${escapeText(term)}</span></div><div class="v"><span class="mono">${escapeText(Number(value).toFixed(6))}</span> <span class="note">${escapeText(direction)}</span></div></div>`;
    })
    .join("");
  const report = candidate.constraint_report;
  let feasible;
  if (!report || typeof report.feasible !== "boolean") {
    feasible = `<p class="note">constraint_report.feasible absent</p>`;
  } else if (report.feasible) {
    feasible = `${chip("feasible", "constraints satisfied")}<p class="field-help">constraint_report.feasible is true. Every hard constraint that was checked passed.</p>`;
  } else {
    feasible = `${chip("infeasible", "constraints not satisfied")}<p class="field-help">constraint_report.feasible is false. At least one hard rule still fails. The sequence is shown so you can read it.</p>`;
  }
  const cds = (candidate.sequences || [])[0] && candidate.sequences[0].cds;
  const lengthNote = cds ? `<p class="note">${cds.length} bases, ${Math.floor(cds.length / 3)} codons.</p>` : "";
  const violations = report && Array.isArray(report.violations) ? report.violations : [];
  const violationList = violations.length
    ? `<ul class="finding-list">${violations
        .map(
          (item) =>
            `<li><span class="mono">${escapeText(item.constraint_id || "")}</span> ${escapeText(item.detail || "")}</li>`,
        )
        .join("")}</ul>`
    : "";
  const copy = cds
    ? `<button type="button" data-copy="${escapeAttr(cds)}">Copy coding sequence</button>`
    : "";
  return `<section class="result-summary">
    <h2>Your design</h2>
    <p class="field-help">This is the recommended candidate. The table further down lists every design the search kept.</p>
    <p class="mono">${escapeText(candidate.candidate_id)}</p>
    ${lengthNote}
    ${feasible}
    ${violationList}
    ${submittedGcWindows()}
    ${scores}
    ${rationaleNotes(candidate)}
    <div id="compare-root"><p class="note">The comparison with what you submitted appears here when this page still has that sequence.</p></div>
    <div class="seq-actions">
      ${copy}
      <button type="button" data-export="fasta">Download FASTA</button>
      <button type="button" data-export="genbank">Download GenBank</button>
      <button type="button" data-export="csv">Download CSV</button>
    </div>
    <p class="field-help">FASTA and CSV include every design on this front. GenBank is the one named above. The numbers are the ones in the result document.</p>
  </section>`;
}

function submittedGcWindows() {
  const constraints = lastSubmitted && lastSubmitted.constraints;
  const rules = constraints && constraints.constraints;
  if (!Array.isArray(rules)) return "";
  const windows = rules.filter((item) => item && item.rule && item.rule.kind === "gc_window");
  if (!windows.length) return "";
  const lines = windows.map((item) => {
    const rule = item.rule;
    return `every ${rule.window_nt} nt, GC fraction between ${rule.min} and ${rule.max} (${item.id || "gc_window"})`;
  });
  return `<p class="field-help"><strong>GC window on this job.</strong> ${escapeText(lines.join("; "))}. A GC window is a sliding stretch and a hard band. It is separate from the single whole-gene GC fraction.</p>`;
}

function rationaleNotes(candidate) {
  const annotations = (candidate.rationale && candidate.rationale.annotations) || [];
  if (!annotations.length) return "";
  const items = annotations
    .map((annotation) => `<li><span class="mono">${escapeText(annotation.tag || "")}</span> ${escapeText(annotation.text || "")}</li>`)
    .join("");
  return `<div class="field-help"><p>Why this sequence.</p><ul>${items}</ul></div>`;
}

function recommendedIndex(front) {
  const id = front.recommended && front.recommended.candidate_id;
  const index = (front.candidates || []).findIndex((item) => item.candidate_id === id);
  return index >= 0 ? index : 0;
}

function sequenceDocumentFor(front, index) {
  const candidate = (front.candidates || [])[index] || (front.candidates || [])[0];
  const cds = candidate && candidate.sequences && candidate.sequences[0] ? candidate.sequences[0].cds : "";
  const five = document.getElementById("flank-5");
  const three = document.getElementById("flank-3");
  const table = document.getElementById("translation-table");
  return {
    cds,
    five_prime_flank: five ? five.value.toUpperCase().replace(/[^ACGT]/g, "") : "",
    three_prime_flank: three ? three.value.toUpperCase().replace(/[^ACGT]/g, "") : "",
    translation_table: table ? Number(table.value) || 11 : 11,
  };
}

function exportFront(kind) {
  const front = window.lastFront;
  if (!front) return;
  const stamp = front.job_id || "design";
  if (kind === "fasta") downloadText(`${stamp}.fasta`, toFasta(front), "text/plain");
  else if (kind === "csv") downloadText(`${stamp}.csv`, toCsv(front), "text/csv");
  else if (kind === "genbank") {
    const index = recommendedIndex(front);
    downloadText(`${stamp}.gb`, toGenbank(front, sequenceDocumentFor(front, index), index), "text/plain");
  }
}

function downloadText(filename, text, type) {
  const blob = new Blob([text], { type });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.appendChild(link);
  link.click();
  link.remove();
  URL.revokeObjectURL(url);
}

function clearHot() {
  frontRoot.querySelectorAll(".is-hot").forEach((node) => node.classList.remove("is-hot"));
  const readout = document.getElementById("plot-readout");
  if (readout) readout.textContent = READOUT_IDLE;
}

function setHot(id) {
  frontRoot.querySelectorAll("[data-candidate]").forEach((node) => {
    node.classList.toggle("is-hot", node.dataset.candidate === id);
  });
  const readout = document.getElementById("plot-readout");
  const row = frontRoot.querySelector(`tr[data-candidate="${id}"]`);
  if (readout && row && row.dataset.readout) readout.textContent = row.dataset.readout;
}

async function copyRaw(text, button) {
  try {
    if (navigator.clipboard && navigator.clipboard.writeText) await navigator.clipboard.writeText(text);
    else throw new Error("clipboard unavailable");
  } catch (error) {
    const area = document.createElement("textarea");
    area.value = text;
    area.setAttribute("readonly", "");
    document.body.appendChild(area);
    area.select();
    document.execCommand("copy");
    area.remove();
  }
  if (button) button.textContent = "Copied raw sequence";
}

function showFrontFailure(response) {
  const payload = (response && response.payload) || {};
  const body = payload.code
    ? problemCallout(payload)
    : `<div class="callout refusal"><strong>Front request failed.</strong> <span class="note">HTTP ${escapeText(response.status)}. The body did not name a CF-* code.</span></div>`;
  frontRoot.innerHTML = lateReveal(body);
}

async function loadFront(jobId) {
  const front = await call("GET", `/v1/jobs/${jobId}/front?reduction=default&max_candidates=20`);
  if (front.status === 200 && front.payload && !front.payload.code) {
    codonHelix.setState(renderFront(front.payload) ? "result" : "error");
    fillComparison(front.payload);
  } else {
    showFrontFailure(front);
    codonHelix.setState("error");
  }
}

let pollToken = 0;
let lastJobStatus = null;
let compareToken = 0;

function gcCount(sequence) {
  const letters = String(sequence || "").toUpperCase();
  if (!letters.length) return null;
  const gc = (letters.match(/[GC]/g) || []).length;
  return { gc, length: letters.length };
}

function codonChanges(before, after) {
  if (!before || !after || before.length !== after.length || before.length % 3 !== 0) return null;
  let changed = 0;
  for (let index = 0; index < before.length; index += 3) {
    if (before.slice(index, index + 3) !== after.slice(index, index + 3)) changed += 1;
  }
  return { changed, total: before.length / 3 };
}

function inputDnaForComparison() {
  if (!lastSubmitted) return "";
  if (typeof lastSubmitted.cds === "string" && lastSubmitted.cds) return lastSubmitted.cds;
  const incoming = lastSubmitted.input;
  if (incoming && incoming.kind === "single_cds" && typeof incoming.sequence === "string") return incoming.sequence;
  return "";
}

function scoreTermsForComparison() {
  const resolved = (lastJobStatus && lastJobStatus.resolved_objectives) || [];
  const terms = resolved.map((term) => term.term_id).filter(Boolean);
  if (terms.length) return terms;
  return Array.isArray(lastSubmitted && lastSubmitted.objectives) ? lastSubmitted.objectives : [];
}

async function fillComparison(front) {
  const slot = document.getElementById("compare-root");
  if (!slot) return;
  const token = ++compareToken;
  const recommended = front.recommended || {};
  const candidates = front.candidates || [];
  const candidate = candidates.find((item) => item.candidate_id === recommended.candidate_id) || candidates[0];
  const designed = candidate && candidate.sequences && candidate.sequences[0] ? candidate.sequences[0].cds || "" : "";
  const incoming = lastSubmitted && lastSubmitted.input;
  const kind = incoming && incoming.kind;
  const input = inputDnaForComparison();
  let body = "";
  if (kind === "protein") {
    body = `<p class="field-help">You pasted a protein. The server reverse-translates it with the codon table on the job: the synonym with the highest relative adaptiveness, and the alphabetically earlier codon when those weights tie. Then it searches. This page does not keep a separate DNA “before”.</p>`;
  } else if (kind === "construct") {
    body = `<p class="field-help">You pasted a GenBank record. The server used one CDS feature. The coding sequence below is the design. This page does not rebuild the rest of the record.</p>`;
  } else if (input && designed && input === designed) {
    body = `<p class="field-help">The recommended coding sequence is the same letters as the DNA you submitted.</p>`;
  } else if (input && designed) {
    const diff = codonChanges(input, designed);
    const before = gcCount(input);
    const after = gcCount(designed);
    const diffLine = diff
      ? `<p>${diff.changed} of ${diff.total} codon positions are not the same letters.</p>`
      : `<p>The two sequences are different lengths (${input.length} and ${designed.length} bases), so this page does not line them up codon by codon.</p>`;
    const countLine = (count) => count ? `${count.gc} G or C of ${count.length} bases` : "";
    body = `<div class="compare">
      <div class="pane"><h3>Submitted DNA</h3><p>${escapeText(countLine(before))}</p><p class="note">Counted on this page. This is not the usage.gc_content score.</p></div>
      <div class="pane"><h3>Recommended coding sequence</h3><p>${escapeText(countLine(after))}</p><p class="note">Counted on this page. This is not the usage.gc_content score.</p></div>
    </div>${diffLine}`;
  } else {
    body = `<p class="note">No submitted DNA is available to put beside this design.</p>`;
  }
  slot.innerHTML = `<h2>Before and after</h2>${body}<div id="score-compare"></div>`;
  const terms = scoreTermsForComparison();
  if (!input || !designed || !lastSubmitted || !lastSubmitted.host_profile || !terms.length) return;
  const scoreSlot = document.getElementById("score-compare");
  if (scoreSlot) scoreSlot.innerHTML = `<p class="note">Scoring the submitted DNA and the recommended sequence on the same terms.</p>`;
  const scored = await call("POST", "/v1/score", {
    host: lastSubmitted.host,
    cds: input,
    host_profile: lastSubmitted.host_profile,
    objectives: terms,
    constraints: lastSubmitted.constraints || EMPTY_CONSTRAINTS,
    five_prime_flank: lastSubmitted.five_prime_flank || "",
    three_prime_flank: lastSubmitted.three_prime_flank || "",
    translation_table: lastSubmitted.translation_table || 11,
    budget: lastSubmitted.budget || "fast",
    seed: Number.isFinite(lastSubmitted.seed) ? lastSubmitted.seed : 0,
    generations: lastSubmitted.generations || 4,
    sequences: [input, designed],
  });
  if (token !== compareToken) return;
  const live = document.getElementById("score-compare");
  if (!live) return;
  if (!scored.payload || scored.payload.code || !Array.isArray(scored.payload.results)) {
    const problem = scored.payload && scored.payload.code ? scored.payload : null;
    live.innerHTML = problem
      ? problemCallout(problem)
      : `<p class="note">The score request did not return numbers for both sequences.</p>`;
    return;
  }
  const [beforeScore, afterScore] = scored.payload.results;
  const rows = terms
    .map((term) => {
      const label = TERM_PLAIN[term] || term;
      const show = (record) => {
        if (!record || record[term] === undefined || record[term] === null) return "UNKNOWN";
        const numeric = Number(record[term]);
        return Number.isFinite(numeric) ? numeric.toFixed(6) : "UNKNOWN";
      };
      return `<tr><td>${escapeText(label)} <span class="mono">${escapeText(term)}</span></td><td class="mono">${escapeText(show(beforeScore))}</td><td class="mono">${escapeText(show(afterScore))}</td></tr>`;
    })
    .join("");
  live.innerHTML = `<p class="field-help">These numbers are from POST /v1/score, on the terms the job resolved. UNKNOWN means that term was not measured.</p>
    <div class="hscroll"><table class="data"><thead><tr><th>Score</th><th>Submitted DNA</th><th>Recommended</th></tr></thead><tbody>${rows}</tbody></table></div>`;
}

async function poll(jobId) {
  const token = ++pollToken;
  codonHelix.setState("running");
  try {
    let status = await call("GET", `/v1/jobs/${jobId}`);
    while (
      token === pollToken &&
      status.payload &&
      status.payload.job_id &&
      ["queued", "running"].includes(status.payload.state)
    ) {
      renderStatus(status.payload);
      refreshJobs();
      await new Promise((resolve) => setTimeout(resolve, POLL_MS));
      if (token !== pollToken) return;
      status = await call("GET", `/v1/jobs/${jobId}`);
    }
    if (token !== pollToken || !status.payload) return;
    if (status.payload.code && !status.payload.job_id) {
      submitError.innerHTML = lateReveal(problemCallout(status.payload));
      frontRoot.innerHTML = lateReveal(problemCallout(status.payload));
      codonHelix.setState("error");
      return;
    }
    renderStatus(status.payload);
    refreshJobs();
    if (status.payload.state !== "succeeded") {
      const problem = status.payload.problem || { code: "", title: "", detail: "" };
      frontRoot.innerHTML = lateReveal(status.payload.problem ? problemCallout(problem) : `<div class="callout refusal"><strong>Job ${escapeText(status.payload.state || "unfinished")}.</strong> <span class="note">No front was requested, because the job did not succeed. The status document did not include a problem code.</span></div>`);
      codonHelix.setState("error");
      return;
    }
    await loadFront(jobId);
  } catch (error) {
    if (token !== pollToken) return;
    submitError.innerHTML = lateReveal(`<div class="callout refusal"><strong>The status request did not complete.</strong> <span class="note">${escapeText(error)}</span></div>`);
    codonHelix.setState("error");
  }
}

async function showJob(jobId) {
  const status = await call("GET", `/v1/jobs/${jobId}`);
  if (!status.payload) return;
  if (status.payload.code && !status.payload.job_id) {
    submitError.innerHTML = lateReveal(problemCallout(status.payload));
    frontRoot.innerHTML = lateReveal(problemCallout(status.payload));
    codonHelix.setState("error");
    return;
  }
  renderStatus(status.payload);
  if (status.payload.state === "queued" || status.payload.state === "running") {
    poll(jobId);
    return;
  }
  if (status.payload.state !== "succeeded") {
    frontRoot.innerHTML = lateReveal(status.payload.problem ? problemCallout(status.payload.problem) : `<div class="callout refusal"><strong>Job ${escapeText(status.payload.state || "unfinished")}.</strong> <span class="note">No front is available.</span></div>`);
    codonHelix.setState("error");
    return;
  }
  await loadFront(jobId);
}

function jobFromHash() {
  const hash = location.hash.replace(/^#/, "");
  if (!hash.startsWith("job_") || !tokenInput.value.trim()) return;
  showJob(hash);
}

async function loadCalibration() {
  const note = document.getElementById("calibration");
  if (!note) return;
  const response = await fetch("/v1/calibration");
  const payload = await response.json();
  if (!response.ok || !payload || (payload.state !== "fitted" && payload.state !== "absent")) {
    const code = payload && payload.code ? payload.code : "calibration";
    const detail = payload && payload.detail ? payload.detail : "The calibration response was not fitted or absent.";
    note.textContent = code + ": " + detail;
    return;
  }
  if (payload.state !== "fitted") {
    note.textContent = "epistasis.context_naturalness absent";
    return;
  }
  const held = Number(payload.held_out_nats_per_codon).toFixed(4).replace("-", "");
  const baseline = Number(payload.baseline_nats_per_codon).toFixed(4).replace("-", "");
  const stop = payload.converged
    ? "converged in " + payload.iteration_count + " iterations"
    : payload.iteration_count + "-iteration early stop";
  note.textContent =
    "epistasis.context_naturalness W=" +
    payload.window +
    ", −" +
    held +
    " nats/codon vs −" +
    baseline +
    " baseline, " +
    stop +
    ".";
}

loadCalibration();

window.addEventListener("hashchange", jobFromHash);
jobFromHash();

frontRoot.addEventListener("pointerover", (event) => {
  const node = event.target.closest("[data-candidate]");
  if (!node) return;
  setHot(node.dataset.candidate);
});

frontRoot.addEventListener("pointerout", (event) => {
  const node = event.target.closest("[data-candidate]");
  if (!node) return;
  const related = event.relatedTarget;
  const element = related && related.nodeType === 1 ? related : related && related.parentElement;
  const next = element && element.closest ? element.closest("[data-candidate]") : null;
  if (next) return;
  clearHot();
});

document.body.addEventListener("click", (event) => {
  const button = event.target.closest("[data-copy]");
  if (!button) return;
  copyRaw(button.dataset.copy || "", button);
});

frontRoot.addEventListener("click", (event) => {
  const exporter = event.target.closest("[data-export]");
  if (exporter) exportFront(exporter.dataset.export);
});

submitButton.addEventListener("click", async () => {
  submitError.innerHTML = "";
  let body;
  try {
    prepareRequest();
    body = JSON.parse(requestInput.value);
  } catch (error) {
    const message = error && error.message ? error.message : "The request is not JSON.";
    submitError.innerHTML = `<div class="callout refusal"><strong>${escapeText(message)}</strong> <span class="note">Nothing was posted.</span></div>`;
    const statusHeading = document.getElementById("status");
    if (statusHeading && message) document.getElementById("submit-error").scrollIntoView({ block: "nearest" });
    return;
  }
  lastSubmitted = body;
  try {
    const posted = await call("POST", "/v1/jobs", body);
    if (!posted.payload) return;
    if (posted.payload.code && !posted.payload.job_id) {
      submitError.innerHTML = lateReveal(problemCallout(posted.payload));
      frontRoot.innerHTML = lateReveal(problemCallout(posted.payload));
      codonHelix.setState("error");
      return;
    }
    history.replaceState(null, "", `#${posted.payload.job_id}`);
    renderStatus(posted.payload);
    const statusHeading = document.getElementById("status");
    if (statusHeading) statusHeading.scrollIntoView({ block: "start" });
    poll(posted.payload.job_id);
  } catch (error) {
    submitError.innerHTML = lateReveal(`<div class="callout refusal"><strong>The request did not complete.</strong> <span class="note">${escapeText(error)}</span></div>`);
    codonHelix.setState("error");
  }
});

window.optimizedCds = "";

function refreshOptimizedHandoff() {
  const slot = document.getElementById("use-optimized");
  if (!slot) return;
  slot.innerHTML = window.optimizedCds
    ? `<p class="note">A coding sequence from Codon design is ready.</p><button type="button" id="use-optimized-sequence">Use that coding sequence</button>`
    : `<p class="note empty-state">No codon design yet. A finished gene on the Codon design tab can be copied into this box.</p>`;
}

document.querySelectorAll("button.tab").forEach((button) => {
  button.addEventListener("click", () => {
    document.querySelectorAll("button.tab").forEach((item) => {
      const on = item === button;
      item.classList.toggle("active", on);
      item.setAttribute("aria-selected", on ? "true" : "false");
    });
    document.getElementById("tab-codon").hidden = button.dataset.tab !== "codon";
    document.getElementById("tab-primer").hidden = button.dataset.tab !== "primer";
    if (button.dataset.tab === "primer") drawViewer();
  });
});

const viewer = document.getElementById("primer-viewer");
const primerSequence = document.getElementById("primer-sequence");
let selection = null;
let dragAnchor = null;
let viewStart = 0;
let hoverSites = null;
let lastPrimerResult = null;

function rulerStep(length) {
  if (length > 5000) return 100;
  if (length > 500) return 50;
  return 10;
}

function visibleWindow(sequence) {
  if (sequence.length <= 10000) return { start: 0, text: sequence };
  const start = Math.max(0, Math.min(viewStart, sequence.length - 2000));
  return { start, text: sequence.slice(start, start + 2000) };
}

function drawViewer() {
  if (!viewer || !primerSequence) return;
  const cleaned = cleanedTemplate(primerSequence.value);
  const sequence = cleaned.sequence;
  const templateProblem = primerTemplateProblem(cleaned);
  const windowed = visibleWindow(sequence);
  viewStart = windowed.start;
  const step = rulerStep(sequence.length);
  let html = "";
  for (let index = 0; index < windowed.text.length; index += 1) {
    const absolute = windowed.start + index;
    if (absolute % step === 0) html += `<span class="note">${absolute}</span> `;
    html += `<span class="nt" data-index="${absolute}">${windowed.text[index]}</span>`;
  }
  if (sequence.length > 10000) {
    html = `<p class="note">Showing ${windowed.start}–${windowed.start + windowed.text.length} of ${sequence.length}. <button type="button" id="viewer-back" data-jump="-2000">Back</button> <button type="button" id="viewer-forward" data-jump="2000">Forward</button></p>` + html;
  }
  if (templateProblem) {
    html = `<p class="callout refusal">${escapeText(templateProblem)}</p>` + html;
  }
  viewer.innerHTML = html || `<span class="note">Paste a sequence.</span>`;
  applyHighlight();
}

function applyHighlight() {
  if (!viewer) return;
  const flank = Number(document.getElementById("primer-flank").value) || 200;
  viewer.querySelectorAll("[data-index]").forEach((node) => {
    const absolute = Number(node.dataset.index);
    let klass = "nt";
    if (hoverSites) {
      if (absolute >= hoverSites.forwardStart && absolute < hoverSites.forwardEnd) klass += " bind-f";
      else if (absolute >= hoverSites.reverseStart && absolute < hoverSites.reverseEnd) klass += " bind-r";
      else if (absolute >= hoverSites.ampStart && absolute < hoverSites.ampEnd) klass += " amplicon";
    } else if (selection && absolute >= selection.start && absolute < selection.end) {
      klass += " selected";
    } else if (selection && absolute >= selection.start - flank && absolute < selection.start) {
      klass += " flank";
    } else if (selection && absolute >= selection.end && absolute < selection.end + flank) {
      klass += " flank";
    }
    if (node.className !== klass) node.className = klass;
  });
  const note = document.getElementById("primer-selection");
  if (note) {
    note.textContent = selection
      ? `Bases ${selection.start + 1}–${selection.end} (${selection.end - selection.start} nt). Both ends are included.`
      : "No region selected.";
  }
  syncRegionFields();
}

function showSelectionTip(event) {
  const tip = document.getElementById("selection-tip");
  if (!tip || !selection) return;
  const bases = selection.end - selection.start;
  tip.hidden = false;
  tip.textContent = bases + " bp selected";
  tip.style.left = event.clientX + 12 + "px";
  tip.style.top = event.clientY + 16 + "px";
}

function indexAt(event) {
  return nucleotideIndex(event.target, document.elementFromPoint(event.clientX, event.clientY));
}

if (viewer) {
  viewer.addEventListener("pointerdown", (event) => {
    const index = indexAt(event);
    if (index === null) {
      const jump = event.target.closest("[data-jump]");
      if (jump) {
        viewStart = Math.max(0, viewStart + Number(jump.dataset.jump));
        drawViewer();
      }
      return;
    }
    dragAnchor = index;
    hoverSites = null;
    selection = { start: index, end: index + 1 };
    viewer.setPointerCapture(event.pointerId);
    applyHighlight();
    showSelectionTip(event);
  });
  viewer.addEventListener("pointermove", (event) => {
    if (dragAnchor === null) return;
    const index = indexAt(event);
    if (index === null) return;
    selection = { start: Math.min(dragAnchor, index), end: Math.max(dragAnchor, index) + 1 };
    applyHighlight();
    showSelectionTip(event);
  });
  window.addEventListener("pointerup", () => {
    dragAnchor = null;
    const tip = document.getElementById("selection-tip");
    if (tip) tip.hidden = true;
  });
}

if (primerSequence) primerSequence.addEventListener("input", drawViewer);

document.getElementById("tab-primer").addEventListener("click", (event) => {
  if (event.target.id !== "use-optimized-sequence") return;
  primerSequence.value = window.optimizedCds;
  selection = null;
  drawViewer();
});

function primerConfig() {
  return {
    tm_target: Number(document.getElementById("primer-tm").value),
    product_size_range: {
      min: Number(document.getElementById("primer-pmin").value),
      max: Number(document.getElementById("primer-pmax").value),
    },
    flank_length: Number(document.getElementById("primer-flank").value),
    primer_length_range: {
      min: Number(document.getElementById("primer-lmin").value),
      max: Number(document.getElementById("primer-lmax").value),
    },
  };
}

function zoneClass(flag) {
  if (!flag) return "zone-green";
  return flag.severity === "red" ? "zone-red" : "zone-amber";
}

function clientFlags(pair) {
  const target = Number(document.getElementById("primer-tm").value);
  const productMin = Number(document.getElementById("primer-pmin").value);
  const productMax = Number(document.getElementById("primer-pmax").value);
  const flags = [];
  [
    ["hairpin_dG_forward", pair.hairpin_dG_forward],
    ["hairpin_dG_reverse", pair.hairpin_dG_reverse],
    ["self_dimer_dG_forward", pair.self_dimer_dG_forward],
    ["self_dimer_dG_reverse", pair.self_dimer_dG_reverse],
    ["hetero_dimer_dG", pair.hetero_dimer_dG],
  ].forEach(([property, observed]) => {
    const value = Number(observed);
    if (Number.isNaN(value)) return;
    if (value < -9) flags.push({ property, observed: value, threshold: "ΔG < -9 kcal/mol", severity: "red" });
    else if (value < -4) flags.push({ property, observed: value, threshold: "ΔG < -4 kcal/mol", severity: "amber" });
  });
  [
    ["tm_forward_vs_target", pair.tm_forward],
    ["tm_reverse_vs_target", pair.tm_reverse],
  ].forEach(([property, observed]) => {
    const delta = Math.abs(Number(observed) - target);
    if (Number.isNaN(delta)) return;
    if (delta > 5) flags.push({ property, observed: delta, threshold: "|Tm − target| > 5 C", severity: "red" });
    else if (delta > 2) flags.push({ property, observed: delta, threshold: "|Tm − target| > 2 C", severity: "amber" });
  });
  if (pair.product_size != null && (pair.product_size < productMin || pair.product_size > productMax)) {
    flags.push({
      property: "product_size",
      observed: pair.product_size,
      threshold: "outside " + productMin + "–" + productMax,
      severity: "red",
    });
  }
  return flags;
}

function worstSeverity(flags) {
  if (flags.some((flag) => flag.severity === "red")) return "red";
  if (flags.some((flag) => flag.severity === "amber")) return "amber";
  return "green";
}

let submittedPrimerTemplate = "";

function productSequenceOf(pair, template) {
  if (pair && pair.product_sequence) return pair.product_sequence;
  if (!pair || !template || pair.forward_start == null || pair.reverse_end == null) return "";
  const start = pair.forward_start;
  const end = pair.reverse_end;
  if (!Number.isInteger(start) || !Number.isInteger(end) || start < 0 || end < start || end > template.length) {
    return "";
  }
  return template.slice(start, end);
}

function detailPanel(pair) {
  const flags = clientFlags(pair)
    .map((flag) => `<div class="zone-chip ${zoneClass(flag)}">${escapeText(flag.property)} ${escapeText(flag.observed)} · ${escapeText(flag.threshold)}</div>`)
    .join("");
  return `<div class="pair-detail">
    <p><strong>Forward</strong> <span class="mono">${escapeText(pair.forward_sequence)}</span> <button type="button" data-copy="${escapeAttr(pair.forward_sequence || "")}">Copy forward primer</button></p>
    <p>GC ${escapeText(pair.gc_forward)}% · Tm ${escapeText(Number(pair.tm_forward).toFixed(2))} °C
      <button type="button" class="param-help" title="${escapeAttr(pair.parameter_citation || "")}">?</button></p>
    <p class="note">${escapeText(pair.parameter_set_id || "")} ${escapeText(pair.parameter_digest || "")}</p>
    <pre class="mono mono-reveal">${escapeText(pair.hairpin_diagram_forward || "")}</pre>
    <p>Hairpin ΔG ${escapeText(pair.hairpin_dG_forward)} · self-dimer ΔG ${escapeText(pair.self_dimer_dG_forward)} · 3′ clamp, informational, not a filter: ${pair.gc_clamp_forward ? "G or C" : "not G or C"}</p>
    <p><strong>Reverse</strong> <span class="mono">${escapeText(pair.reverse_sequence)}</span> <button type="button" data-copy="${escapeAttr(pair.reverse_sequence || "")}">Copy reverse primer</button></p>
    <pre class="mono mono-reveal">${escapeText(pair.hairpin_diagram_reverse || "")}</pre>
    <p>Hairpin ΔG ${escapeText(pair.hairpin_dG_reverse)} · self-dimer ΔG ${escapeText(pair.self_dimer_dG_reverse)} · 3′ clamp, informational, not a filter: ${pair.gc_clamp_reverse ? "G or C" : "not G or C"}</p>
    <p><strong>Pair</strong> hetero-dimer ΔG ${escapeText(pair.hetero_dimer_dG)} · 3′ ΔG ${escapeText(pair.hetero_dimer_3prime_dG)} · Tm mismatch ${escapeText(pair.tm_mismatch)} · product ${escapeText(pair.product_size)}</p>
    <p class="note">Self-dimer and heterodimer ΔG use SantaLucia 1998 Table 2 Watson–Crick stacks only: no dangling ends and no mismatch stacks, so they are not primer3 calc_homodimer. Melting temperature is the same SantaLucia model and matches primer3-py at 1 M sodium.</p>
    <pre class="mono mono-reveal">${escapeText(pair.hetero_dimer_alignment || "")}</pre>
    <div class="flag-host">${flags || `<p class="zone-chip zone-green">No property is in the amber or red zone.</p>`}</div>
    <p class="note">Specificity against genome not checked.</p>
    ${(() => {
      const product = productSequenceOf(pair, submittedPrimerTemplate);
      return product ? `<button type="button" data-copy="${escapeAttr(product)}">Copy product sequence</button><div class="seq-block">${escapeText(product)}</div>` : "";
    })()}
  </div>`;
}

function renderPrimerResult(result) {
  const root = document.getElementById("primer-results");
  lastPrimerResult = result;
  const pairs = (result && result.pairs) || [];
  if (!pairs.length) {
    const detail = (result && result.detail) || "No filter counts were returned with this empty result.";
    root.innerHTML = lateReveal(problemCallout({ code: "CF-PRIM-001", title: "No primer pairs survived pruning", detail }));
    return false;
  }
  const recommended = pairs[0];
  const rows = pairs
    .map((pair, index) => {
      const worst = Math.min(pair.hairpin_dG_forward, pair.hairpin_dG_reverse, pair.self_dimer_dG_forward, pair.self_dimer_dG_reverse, pair.hetero_dimer_dG);
      const severity = worstSeverity(clientFlags(pair));
      const zone = severity === "red" ? chip("refusal", "red") : severity === "amber" ? chip("warning", "amber") : chip("feasible", "kept");
      return `<tr data-pair="${index}"><td class="mono">${escapeText(pair.pair_id)}</td><td>${escapeText(Number(pair.tm_mismatch).toFixed(2))}</td><td>${escapeText(Number(worst).toFixed(2))}</td><td>${escapeText(Math.abs(pair.gc_forward - pair.gc_reverse).toFixed(1))}</td><td>${escapeText(pair.product_size)}</td><td class="zone-chip">${zone}</td></tr>`;
    })
    .join("");
  root.innerHTML = `<div class="reveal primer-results-layout">
    <div class="callout pulse-once"><strong>Recommended</strong> <code>${escapeText(recommended.pair_id)}</code><br>
    informed_by <code>primer.tm_mismatch</code> <code>primer.worst_dimer_penalty</code> <code>primer.gc_balance</code> <code>primer.product_deviation</code></div>
    <p class="note">Specificity against genome not checked. ${escapeText(result.parameter_citation || "")}</p>
    <div class="hscroll"><table class="data"><thead><tr><th>pair</th><th>Tm mismatch</th><th>worst dimer ΔG</th><th>GC balance</th><th>product</th><th>zone</th></tr></thead><tbody>${rows}</tbody></table></div>
    <div id="primer-detail">${detailPanel(recommended)}</div>
  </div>`;
  root.querySelectorAll("tr[data-pair]").forEach((row) => {
    const pair = pairs[Number(row.dataset.pair)];
    row.addEventListener("click", () => {
      root.querySelectorAll("tr.open").forEach((other) => other.classList.remove("open"));
      row.classList.add("open");
      document.getElementById("primer-detail").innerHTML = detailPanel(pair);
    });
    row.addEventListener("pointerenter", () => {
      const ampStart = Math.min(pair.forward_start, pair.reverse_start);
      const ampEnd = Math.max(pair.forward_end, pair.reverse_end);
      hoverSites = {
        forwardStart: pair.forward_start,
        forwardEnd: pair.forward_end,
        reverseStart: pair.reverse_start,
        reverseEnd: pair.reverse_end,
        ampStart,
        ampEnd,
      };
      applyHighlight();
    });
    row.addEventListener("pointerleave", () => {
      hoverSites = null;
      applyHighlight();
    });
  });
  const first = root.querySelector("tr[data-pair]");
  if (first) first.classList.add("open");
  return true;
}

function primerProgress(payload) {
  const tiers = (payload.progress && payload.progress.tiers) || [];
  const moving = tiers.some((tier) => Number(tier.fraction) > 0);
  if (!tiers.length || !moving) {
    return `<div class="tier-track"><span>primers</span><span>${escapeText(payload.state || "")}</span><div class="bar indeterminate"><span></span></div><span></span></div>`;
  }
  return tiers
    .map((tier) => {
      const width = Math.round((Number(tier.fraction) || 0) * 100);
      return `<div class="tier-track"><span>${escapeText(tier.tier)}</span><span>${escapeText(tier.state)}</span><div class="bar"><span style="width:${width}%"></span></div><span></span></div>`;
    })
    .join("");
}

function showPrimerFailure(error) {
  document.getElementById("primer-error").innerHTML = lateReveal(`<div class="callout refusal"><strong>The request did not complete.</strong> <span class="note">${escapeText(error)}</span></div>`);
  primerHelix.setState("error");
}

async function pollPrimers(jobId) {
  const bar = document.getElementById("primer-progress");
  primerHelix.setState("running");
  try {
    let status = await call("GET", `/v1/primers/${jobId}`);
    while (status.payload && ["queued", "running"].includes(status.payload.state)) {
      bar.innerHTML = primerProgress(status.payload);
      await new Promise((resolve) => setTimeout(resolve, 2000));
      status = await call("GET", `/v1/primers/${jobId}`);
    }
    if (!status.payload) return;
    bar.innerHTML = primerProgress(status.payload);
    if (status.payload.problem) {
      document.getElementById("primer-error").innerHTML = lateReveal(problemCallout(status.payload.problem));
      document.getElementById("primer-results").innerHTML = lateReveal(problemCallout(status.payload.problem));
      primerHelix.setState("error");
      return;
    }
    primerHelix.setState(renderPrimerResult(status.payload.primer_result) ? "result" : "error");
  } catch (error) {
    showPrimerFailure(error);
  }
}

document.getElementById("primer-search").addEventListener("click", async () => {
  document.getElementById("primer-error").innerHTML = "";
  document.getElementById("primer-results").innerHTML = "";
  if (!selection) {
    document.getElementById("primer-error").innerHTML = `<div class="callout refusal"><strong>No region selected.</strong> <span class="note">Drag across the sequence, type a start and end, or use Select the whole sequence.</span></div>`;
    return;
  }
  const cleaned = cleanedTemplate(primerSequence.value);
  const templateProblem = primerTemplateProblem(cleaned);
  if (templateProblem) {
    document.getElementById("primer-error").innerHTML = `<div class="callout refusal"><strong>${escapeText(templateProblem)}</strong></div>`;
    return;
  }
  const sequence = cleaned.sequence;
  submittedPrimerTemplate = sequence;
  let posted;
  try {
    posted = await call("POST", "/v1/primers/candidates", {
      sequence,
      region: selection,
      config: primerConfig(),
    });
  } catch (error) {
    showPrimerFailure(error);
    return;
  }
  if (!posted.payload) return;
  if (posted.payload.code && !posted.payload.job_id) {
    document.getElementById("primer-error").innerHTML = lateReveal(problemCallout(posted.payload));
    document.getElementById("primer-results").innerHTML = lateReveal(problemCallout(posted.payload));
    primerHelix.setState("error");
    return;
  }
  pollPrimers(posted.payload.job_id);
});

document.getElementById("manual-check").addEventListener("click", async () => {
  const posted = await call("POST", "/v1/primers/analyze", {
    forward: document.getElementById("manual-forward").value.toUpperCase().replace(/[^ACGT]/g, ""),
    reverse: document.getElementById("manual-reverse").value.toUpperCase().replace(/[^ACGT]/g, ""),
  });
  const root = document.getElementById("manual-detail");
  if (!posted.payload) return;
  if (posted.payload.code && !posted.payload.forward_sequence) {
    root.innerHTML = problemCallout(posted.payload);
    return;
  }
  root.innerHTML = `<div class="reveal">${detailPanel(posted.payload)}</div>`;
});

function refreshPrimerZones() {
  if (!lastPrimerResult || !lastPrimerResult.pairs) return;
  const rows = document.querySelectorAll("#primer-results tr[data-pair]");
  if (!rows.length) return;
  rows.forEach((row) => {
    const pair = lastPrimerResult.pairs[Number(row.dataset.pair)];
    const severity = worstSeverity(clientFlags(pair));
    const node = row.querySelector(".chip");
    if (!node) return;
    node.className = "chip " + (severity === "red" ? "refusal" : severity === "amber" ? "warning" : "feasible");
    node.textContent = severity === "green" ? "kept" : severity;
  });
  const open = document.querySelector("#primer-results tr.open");
  const host = document.querySelector("#primer-detail .flag-host");
  if (open && host) {
    const pair = lastPrimerResult.pairs[Number(open.dataset.pair)];
    const flags = clientFlags(pair);
    host.innerHTML = flags.length
      ? flags
          .map((flag) => `<div class="zone-chip ${zoneClass(flag)}">${escapeText(flag.property)} ${escapeText(flag.observed)} · ${escapeText(flag.threshold)}</div>`)
          .join("")
      : `<p class="zone-chip zone-green">No property is in the amber or red zone.</p>`;
  }
}

["primer-tm", "primer-pmin", "primer-pmax"].forEach((id) => {
  const input = document.getElementById(id);
  if (!input) return;
  input.addEventListener("input", () => {
    refreshPrimerZones();
    applyHighlight();
  });
});

drawViewer();
refreshOptimizedHandoff();

const EMPTY_CONSTRAINTS = {
  schema_version: "1.0.0",
  constraint_set_id: "cst_none",
  constraints: [],
};

let loadedProfile = null;
let loadedSource = null;
let loadedConstraints = null;

function syncRegionFields() {
  const start = document.getElementById("region-start");
  const end = document.getElementById("region-end");
  if (!start || !end || !selection) return;
  start.value = String(selection.start + 1);
  end.value = String(selection.end);
}

function selectionFromFields() {
  if (!primerSequence) return;
  const sequence = cleanedTemplate(primerSequence.value).sequence;
  const startBox = document.getElementById("region-start");
  const endBox = document.getElementById("region-end");
  const startBase = Number(startBox.value);
  const endBase = Number(endBox.value);
  if (!Number.isInteger(startBase) || !Number.isInteger(endBase)) return;
  const start = Math.max(0, startBase - 1);
  const end = Math.min(sequence.length, endBase);
  if (end <= start) return;
  selection = { start, end };
  applyHighlight();
}

function cleanedTemplate(raw) {
  const lines = String(raw).replaceAll("\r\n", "\n").replaceAll("\r", "\n").split("\n");
  const kept = [];
  for (const line of lines) {
    const trimmed = line.trim();
    if (!trimmed || trimmed.startsWith(">") || trimmed.startsWith(";")) continue;
    kept.push(trimmed);
  }
  const sequence = kept.join("").toUpperCase().replaceAll(/\s+/g, "");
  const illegal = [];
  for (let index = 0; index < sequence.length; index += 1) {
    const base = sequence[index];
    if (!"ACGT".includes(base)) illegal.push({ base, position: index + 1 });
  }
  return { sequence, illegal };
}

function primerTemplateProblem(cleaned) {
  if (!cleaned.illegal.length) return "";
  const listed = cleaned.illegal
    .slice(0, 4)
    .map((item) => `base ${item.base} at position ${item.position}`)
    .join(", ");
  const more = cleaned.illegal.length > 4 ? `, and ${cleaned.illegal.length - 4} more` : "";
  return `CF-PRIM-006 ${listed}${more} is not A, C, G or T.`;
}

function firstStopTrim(sequence) {
  const stops = ["TAA", "TAG", "TGA"];
  for (let index = 0; index + 3 <= sequence.length; index += 3) {
    if (stops.includes(sequence.slice(index, index + 3))) {
      return sequence.slice(0, index + 3);
    }
  }
  return sequence;
}

function cleanedCds(raw) {
  const lines = String(raw).replaceAll("\r\n", "\n").replaceAll("\r", "\n").split("\n");
  const kept = [];
  let header = false;
  for (const line of lines) {
    const trimmed = line.trim();
    if (!trimmed) continue;
    if (trimmed.startsWith(">") || trimmed.startsWith(";")) {
      header = true;
      continue;
    }
    kept.push(trimmed);
  }
  let sequence = kept.join("").toUpperCase().replaceAll(/\s+/g, "");
  const notes = [];
  if (header) notes.push("Removed the FASTA header.");
  const rna = (sequence.match(/U/g) || []).length;
  if (rna) {
    sequence = sequence.replaceAll("U", "T");
    notes.push(`Read ${rna} U as T, so an RNA paste is treated as DNA.`);
  }
  const digits = (sequence.match(/[0-9]/g) || []).length;
  if (digits) {
    sequence = sequence.replaceAll(/[0-9]/g, "");
    notes.push("Removed position numbers.");
  }
  const trimBox = document.getElementById("trim-stop");
  if (trimBox && trimBox.checked && sequence && sequence.length % 3 === 0 && !/[^ACGT]/.test(sequence)) {
    const trimmed = firstStopTrim(sequence);
    if (trimmed.length !== sequence.length) {
      const codon = trimmed.length / 3 - 1;
      notes.push(`Trimmed at the first in-frame stop (standard code): ${sequence.length} bases became ${trimmed.length}, through codon ${codon} inclusive.`);
      sequence = trimmed;
    }
  }
  const illegal = sequence.match(/[^ACGT]/g);
  return { sequence, notes, illegal: illegal ? [...new Set(illegal)] : [] };
}

function inputKind() {
  const node = document.getElementById("input-kind");
  return node ? node.value : "dna";
}

function cleanedProtein(raw) {
  const lines = String(raw || "").split(/\r?\n/);
  const kept = [];
  const notes = [];
  let header = false;
  for (const line of lines) {
    const trimmed = line.trim();
    if (!trimmed) continue;
    if (trimmed.startsWith(">") || trimmed.startsWith(";")) {
      header = true;
      continue;
    }
    kept.push(trimmed);
  }
  if (header) notes.push("Removed the FASTA header.");
  const sequence = kept.join("").replaceAll(/\s+/g, "").toUpperCase();
  const illegal = sequence.match(/[^ACDEFGHIKLMNPQRSTVWY*]/g);
  return { sequence, notes, illegal: illegal ? [...new Set(illegal)] : [] };
}

function applyInputKind() {
  const kind = inputKind();
  const box = document.getElementById("cds-text");
  const locus = document.getElementById("locus-wrap");
  const locusHelp = document.getElementById("help-locus");
  const trim = document.getElementById("trim-stop");
  const trimLabel = trim && trim.closest("label");
  if (locus) locus.hidden = kind !== "construct";
  if (locusHelp) locusHelp.hidden = kind !== "construct";
  if (trimLabel) trimLabel.hidden = kind !== "dna";
  if (box) {
    if (kind === "protein") box.placeholder = "MKT... or a FASTA protein record";
    else if (kind === "construct") box.placeholder = "LOCUS ... ORIGIN ... //";
    else box.placeholder = "ATG... or a FASTA record";
  }
  refreshSequenceNote();
  syncRequest();
}

function refreshSequenceNote() {
  const note = document.getElementById("cds-note");
  const box = document.getElementById("cds-text");
  if (!note || !box) return;
  if (!box.value.trim()) {
    note.textContent = "No sequence yet. The rpsL demo is the E. coli K-12 MG1655 gene rpsL, 125 codons.";
    return;
  }
  const kind = inputKind();
  if (kind === "protein") {
    const parsed = cleanedProtein(box.value);
    const parts = [...parsed.notes];
    if (!parsed.sequence) parts.push("No amino-acid letters were found.");
    else if (parsed.illegal.length) parts.push(`These letters are not amino acids in this form: ${parsed.illegal.join(", ")}.`);
    else {
      const residues = parsed.sequence.replaceAll("*", "");
      parts.push(`${residues.length} amino acids${parsed.sequence.endsWith("*") ? ", with a stop written at the end" : ""}.`);
      if (/^[ACGT]+$/.test(residues) && residues.length >= 30 && residues.length % 3 === 0) {
        parts.push("These letters are also DNA bases. If this is a gene, switch to DNA coding sequence.");
      }
    }
    note.textContent = parts.join(" ");
    return;
  }
  if (kind === "construct") {
    const text = box.value;
    const parts = [];
    if (!/LOCUS/i.test(text)) parts.push("No LOCUS line yet.");
    if (!/ORIGIN/i.test(text)) parts.push("No ORIGIN line yet.");
    const locus = document.getElementById("locus");
    parts.push(locus && locus.value.trim() ? `CDS name ${locus.value.trim()}.` : "No CDS name. That is fine when the file has one CDS.");
    note.textContent = parts.join(" ");
    return;
  }
  const parsed = cleanedCds(box.value);
  const parts = [...parsed.notes];
  if (!parsed.sequence) {
    parts.push("No A, C, G or T bases were found.");
  } else if (parsed.illegal.length) {
    parts.push(`These letters are not DNA bases: ${parsed.illegal.join(", ")}. For a protein, change “What you are pasting” to Protein sequence.`);
  } else if (parsed.sequence.length % 3 !== 0) {
    parts.push(`${parsed.sequence.length} bases. That is not a multiple of 3, so it is not a whole number of codons.`);
  } else {
    parts.push(`${parsed.sequence.length} bases, ${parsed.sequence.length / 3} codons.`);
    const trim = document.getElementById("trim-stop");
    if (parsed.sequence.endsWith("TAATAA") && trim && !trim.checked) {
      parts.push("The sequence ends with TAATAA. Trim at the first in-frame stop is off unless you tick it.");
    }
  }
  note.textContent = parts.join(" ");
}

function setProfileNote(text) {
  const note = document.getElementById("profile-note");
  if (note) note.textContent = text;
}

async function loadBundledProfile() {
  const response = await fetch("/examples/mg1655_ribosomal_usage.json");
  if (!response.ok) {
    setProfileNote("The bundled E. coli table could not be loaded. Upload a codon-usage JSON instead.");
    return;
  }
  loadedProfile = await response.json();
  setProfileNote("Loaded examples/mg1655_ribosomal_usage.json. E. coli K-12 MG1655 ribosomal proteins, labelled user-supplied.");
  syncRequest();
}

function constraintFromBuilder() {
  const rules = [];
  if (document.getElementById("constraint-bsai").checked) {
    rules.push({
      id: "c_bsai",
      class: "blocking",
      scope: "cds",
      rule: { kind: "restriction_site", enzyme: "BsaI", site: "GGTCTC", strand: "both", max_occurrences: 0 },
    });
  }
  if (document.getElementById("constraint-bsmbi").checked) {
    rules.push({
      id: "c_bsmbi",
      class: "blocking",
      scope: "cds",
      rule: { kind: "restriction_site", enzyme: "BsmBI", site: "CGTCTC", strand: "both", max_occurrences: 0 },
    });
  }
  if (document.getElementById("constraint-homo").checked) {
    rules.push({
      id: "c_homopolymer",
      class: "blocking",
      scope: "cds",
      rule: { kind: "homopolymer", base: "any", max_run: Number(document.getElementById("constraint-homo-run").value) },
    });
  }
  if (document.getElementById("constraint-gc").checked) {
    rules.push({
      id: "c_gc_window",
      class: "blocking",
      scope: "cds",
      rule: {
        kind: "gc_window",
        window_nt: Number(document.getElementById("constraint-gc-window").value),
        min: Number(document.getElementById("constraint-gc-min").value),
        max: Number(document.getElementById("constraint-gc-max").value),
      },
    });
  }
  return {
    schema_version: "1.0.0",
    constraint_set_id: "cst_form",
    constraints: rules,
  };
}

function constraintsForRequest() {
  const mode = document.getElementById("constraint-mode").value;
  if (mode === "none") return EMPTY_CONSTRAINTS;
  if (mode === "example") {
    if (!loadedConstraints) throw new Error("The rpsL constraint file has not loaded yet. Choose it again, or use No extra rules.");
    return loadedConstraints;
  }
  if (mode === "builder") return constraintFromBuilder();
  const text = document.getElementById("constraint-json").value.trim();
  if (!text) throw new Error("Paste a constraint JSON, or choose No extra rules.");
  try {
    return JSON.parse(text);
  } catch (error) {
    throw new Error("The constraint box is not JSON. Nothing was posted.");
  }
}

const HOST_ID = /^[a-z0-9]+(?:_[a-z0-9]+)*:[0-9]+(?:\/[A-Za-z0-9._-]+)?$/;

function selectedHost() {
  const named = document.getElementById("host-name");
  const select = document.getElementById("host-select");
  const custom = select && select.value === "custom" ? document.getElementById("host-custom").value.trim() : "";
  const host = (named && named.value.trim()) || custom;
  if (!host) {
    throw new Error("Type a host name such as escherichia_coli:511145/MG1655. The menu chooses the codon table, which is separate from this name.");
  }
  if (!HOST_ID.test(host)) {
    throw new Error(`“${host}” is not a host name this server accepts. Use genus_species:taxid/strain, for example escherichia_coli:511145/MG1655.`);
  }
  return host;
}

function objectiveFields() {
  const preset = document.getElementById("objective-preset");
  if (preset && preset.value && preset.value !== "custom") return { objectives_preset: preset.value };
  return { objectives: selectedObjectives() };
}

function selectedObjectives() {
  const chosen = [];
  if (document.getElementById("obj-cai").checked) chosen.push("usage.cai");
  if (document.getElementById("obj-harmonize") && document.getElementById("obj-harmonize").checked) {
    chosen.push("usage.harmonize_rca");
  }
  if (document.getElementById("obj-match") && document.getElementById("obj-match").checked) {
    chosen.push("usage.match_frequency");
  }
  if (document.getElementById("obj-gc").checked) chosen.push("usage.gc_content");
  if (document.getElementById("obj-violation").checked) chosen.push("constraint.violation_count");
  if (document.getElementById("obj-structure").checked) chosen.push("structure.local_accessibility");
  if (!chosen.length) throw new Error("Choose at least one objective. CAI is the usual first choice.");
  return chosen;
}

function flankValue(id) {
  return document.getElementById(id).value.toUpperCase().replace(/[^ACGT]/g, "");
}

function buildRequest() {
  const kind = inputKind();
  const box = document.getElementById("cds-text");
  const raw = box.value;
  if (!loadedProfile || !loadedProfile.codon_usage) {
    throw new Error("Load a codon-usage table. The list is the files installed on this computer. The first button loads examples/mg1655_ribosomal_usage.json.");
  }
  const profile = JSON.parse(JSON.stringify(loadedProfile));
  const preset = document.getElementById("objective-preset");
  const harmonize = (preset && preset.value === "harmonize")
    || (document.getElementById("obj-harmonize") && document.getElementById("obj-harmonize").checked && preset && preset.value === "custom");
  if (harmonize) {
    if (!loadedSource || !loadedSource.codon_usage) {
      throw new Error("Harmonization needs a source organism table. Upload one next to that score. The expression host table is the one already loaded, and this page will not invent the source.");
    }
    profile.source_codon_usage = loadedSource.codon_usage;
  }
  const body = {
    host: selectedHost(),
    host_profile: profile,
    constraints: constraintsForRequest(),
    five_prime_flank: flankValue("flank-5"),
    three_prime_flank: flankValue("flank-3"),
    translation_table: Number(document.getElementById("translation-table").value) || 11,
    budget: document.getElementById("budget").value,
    seed: Number(document.getElementById("seed").value) || 0,
    generations: Number(document.getElementById("generations").value) || 4,
    ...objectiveFields(),
  };
  const starting = document.getElementById("starting-cds");
  if (starting && starting.value.trim()) {
    body.starting_cds = starting.value.toUpperCase().replace(/[^ACGT]/g, "");
  }
  const engine = document.getElementById("structure-engine");
  if (engine && engine.value === "nussinov") {
    body.structure_engine = "nussinov";
  }
  if (kind === "protein") {
    const parsed = cleanedProtein(raw);
    if (!parsed.sequence) throw new Error("Paste a protein sequence, or switch to DNA and load the rpsL demo.");
    if (parsed.illegal.length) {
      throw new Error(`These letters are not amino acids in this form: ${parsed.illegal.join(", ")}. The server also refuses a residue the translation table does not name.`);
    }
    body.input = { kind: "protein", protein: parsed.sequence };
    return body;
  }
  if (kind === "construct") {
    const text = raw.trim();
    if (!text) throw new Error("Paste a GenBank record, or switch to DNA and load the rpsL demo.");
    if (!/LOCUS/i.test(text) && !/ORIGIN/i.test(text)) {
      throw new Error("That text does not look like a GenBank record. It needs a LOCUS line and an ORIGIN. A coding sequence belongs under DNA coding sequence.");
    }
    const locus = document.getElementById("locus").value.trim();
    const input = { kind: "construct", genbank: text };
    if (locus) input.locus = locus;
    body.input = input;
    return body;
  }
  const parsed = cleanedCds(raw);
  if (!parsed.sequence) throw new Error("Paste a coding sequence, or load the rpsL demo.");
  if (parsed.illegal.length) {
    throw new Error(`These letters are not DNA bases: ${parsed.illegal.join(", ")}. For a protein, change “What you are pasting” to Protein sequence.`);
  }
  if (parsed.sequence.length % 3 !== 0) {
    throw new Error(`The sequence is ${parsed.sequence.length} bases long, which is not a whole number of codons. Check that the copy is complete.`);
  }
  body.input = { kind: "single_cds", sequence: parsed.sequence };
  return body;
}

function syncRequest() {
  const edit = document.getElementById("edit-json");
  if (!requestInput || (edit && edit.checked)) return;
  try {
    requestInput.value = JSON.stringify(buildRequest(), null, 2);
  } catch (error) {
    requestInput.value = "";
  }
}

function prepareRequest() {
  const edit = document.getElementById("edit-json");
  if (edit && edit.checked) {
    if (!requestInput.value.trim()) throw new Error("The request box is empty. Untick Edit this JSON myself, or paste a request.");
    return;
  }
  requestInput.value = JSON.stringify(buildRequest(), null, 2);
}

function showConstraintMode() {
  const mode = document.getElementById("constraint-mode").value;
  document.getElementById("constraint-builder").hidden = mode !== "builder";
  document.getElementById("constraint-json-wrap").hidden = mode !== "json" && mode !== "example";
  const jsonBox = document.getElementById("constraint-json");
  jsonBox.readOnly = mode === "example";
  syncRequest();
}

async function loadConstraintPreset(file) {
  const response = await fetch(`/examples/${file}`);
  if (!response.ok) throw new Error(`${file} could not be loaded.`);
  const preset = await response.json();
  document.getElementById("constraint-mode").value = "json";
  const box = document.getElementById("constraint-json");
  box.value = JSON.stringify(preset, null, 2);
  box.readOnly = false;
  showConstraintMode();
}

async function loadKazusaProfile() {
  const response = await fetch("/examples/kazusa_e_coli_316407.json");
  if (!response.ok) {
    setProfileNote("The Kazusa table could not be loaded. Upload a codon-usage JSON instead.");
    return;
  }
  loadedProfile = await response.json();
  setProfileNote("Loaded examples/kazusa_e_coli_316407.json. Genome-wide E. coli frequencies (taxid 316407) converted to relative adaptiveness. This replaces the ribosomal table for this search. The result still names the table that was loaded.");
  syncRequest();
}

async function loadExampleConstraints() {
  const response = await fetch("/examples/rpsL_constraints.json");
  if (!response.ok) throw new Error("The rpsL constraint file could not be loaded.");
  loadedConstraints = await response.json();
  document.getElementById("constraint-json").value = JSON.stringify(loadedConstraints, null, 2);
  syncRequest();
}

const hostsById = new Map();

function hostOptionLabel(host) {
  const usage = host.profile && host.profile.codon_usage;
  const source = usage && typeof usage.source === "string" ? usage.source.trim() : "";
  const filename = host.filename || host.id;
  if (!source) return filename;
  const sentence = source.split(". ")[0];
  return `${filename} — ${sentence}`;
}

function applySelectedHost() {
  const select = document.getElementById("host-select");
  const wrap = document.getElementById("host-custom-wrap");
  if (!select) return;
  if (wrap) wrap.hidden = select.value !== "custom";
  if (select.value === "custom" || !select.value) {
    syncRequest();
    return;
  }
  const host = hostsById.get(select.value);
  if (host && host.profile && host.profile.codon_usage) {
    loadedProfile = host.profile;
    const source = host.profile.codon_usage.source;
    setProfileNote(
      source
        ? `Loaded ${host.filename} from GET /v1/hosts. ${source}`
        : `Loaded ${host.filename} from GET /v1/hosts.`,
    );
  }
  syncRequest();
}

async function selectInstalledHost(filename, fallback) {
  const host = [...hostsById.values()].find((item) => item.filename === filename);
  if (host) {
    document.getElementById("host-select").value = host.id;
    applySelectedHost();
    return;
  }
  await fallback();
}

let hostsPromise = null;

async function loadHosts() {
  if (!hostsPromise) {
    hostsPromise = loadHostsOnce().catch((error) => {
      hostsPromise = null;
      throw error;
    });
  }
  return hostsPromise;
}

async function loadHostsOnce() {
  const select = document.getElementById("host-select");
  const response = await call("GET", "/v1/hosts");
  if (!select) return;
  if (!response.payload || response.payload.code || !Array.isArray(response.payload.hosts)) {
    const detail = response.payload && (response.payload.detail || response.payload.code);
    setProfileNote(detail ? `The codon-table list did not load. ${detail}` : "The codon-table list did not load.");
    if (!loadedProfile) await loadBundledProfile();
    if (loadedProfile && select && ![...select.options].some((option) => option.value === "mg1655_ribosomal_usage")) {
      const option = document.createElement("option");
      option.value = "mg1655_ribosomal_usage";
      option.textContent = "mg1655_ribosomal_usage.json";
      select.replaceChildren(option);
      const custom = document.createElement("option");
      custom.value = "custom";
      custom.textContent = "Other host (type the name and upload a table)";
      select.appendChild(custom);
      select.value = "mg1655_ribosomal_usage";
    }
    return;
  }
  const previous = select.value;
  hostsById.clear();
  select.replaceChildren();
  for (const host of response.payload.hosts) {
    if (!host || !host.id || host.installed === false) continue;
    hostsById.set(host.id, host);
    const option = document.createElement("option");
    option.value = host.id;
    option.textContent = hostOptionLabel(host);
    select.appendChild(option);
  }
  const custom = document.createElement("option");
  custom.value = "custom";
  custom.textContent = "Other host (type the name and upload a table)";
  select.appendChild(custom);
  if (!hostsById.size) {
    setProfileNote("No codon-usage file is installed. Upload a usage JSON.");
    select.value = "custom";
    applySelectedHost();
    return;
  }
  const preferred = hostsById.has(previous)
    ? previous
    : hostsById.has("mg1655_ribosomal_usage")
      ? "mg1655_ribosomal_usage"
      : [...hostsById.keys()][0];
  select.value = preferred;
  applySelectedHost();
}

async function loadModels() {
  const note = document.getElementById("model-status");
  if (!note) return;
  const response = await call("GET", "/v1/models");
  if (!response.payload || response.payload.code) {
    const code = response.payload && response.payload.code ? response.payload.code : "";
    const detail = response.payload && response.payload.detail ? response.payload.detail : "The model list did not load.";
    note.textContent = code ? `${code}: ${detail}` : detail;
    return;
  }
  const models = Array.isArray(response.payload.models) ? response.payload.models : [];
  if (response.payload.installed && models.length) {
    note.textContent = `Installed model: ${models.map((item) => item.id).filter(Boolean).join(", ")}.`;
    return;
  }
  note.textContent = "No context model is installed. epistasis.context_naturalness will not be scored. CAI does not need it.";
}

async function loadEngineStatus() {
  const note = document.getElementById("engine-status");
  if (!note) return;
  try {
    const response = await fetch("/v1/healthz");
    const payload = await response.json();
    const engines = Array.isArray(payload.structure_engines) ? payload.structure_engines : [];
    if (engines.includes("viennarna")) {
      note.textContent = "ViennaRNA is installed on this computer. Local mRNA structure can be measured. A score that is still missing is UNKNOWN, not zero.";
      return;
    }
    if (engines.length) {
      note.textContent = `Structure engines reported by this server: ${engines.join(", ")}. ViennaRNA is not installed. DNA Design Suite does not include it, and a structure score that is not measured is UNKNOWN, not zero.`;
      return;
    }
    note.textContent = "ViennaRNA is not installed. DNA Design Suite does not include a folding program. CAI does not need one. If you ask for local mRNA structure, a score that is not measured is UNKNOWN, not zero.";
  } catch (error) {
    note.textContent = "The page could not read whether ViennaRNA is installed.";
  }
}

async function loadVersion() {
  const note = document.getElementById("suite-version");
  if (!note) return;
  try {
    const response = await fetch("/v1/version");
    const payload = await response.json();
    if (payload && payload.version) note.textContent = `DNA Design Suite ${payload.version}. Codon design rewrites a gene for a host. Primer design picks a PCR pair.`;
  } catch (error) {
    note.textContent = "DNA Design Suite. Codon design rewrites a gene for a host. Primer design picks a PCR pair.";
  }
}

function stateChipKind(state) {
  if (state === "succeeded") return "feasible";
  if (state === "failed") return "refusal";
  return "warning";
}

async function refreshJobs() {
  const root = document.getElementById("jobs-root");
  if (!root) return;
  const listed = await call("GET", "/v1/jobs");
  if (!listed.payload || listed.payload.code) {
    root.innerHTML = listed.payload && listed.payload.code
      ? problemCallout(listed.payload)
      : `<p class="note empty-state">The job list did not load.</p>`;
    return;
  }
  const items = [...(listed.payload.items || [])];
  items.sort((left, right) => String(right.created_at || "").localeCompare(String(left.created_at || "")) || String(left.job_id || "").localeCompare(String(right.job_id || "")));
  if (!items.length) {
    root.innerHTML = `<p class="note empty-state">No jobs yet. Run design starts one in the background. This list is filled from GET /v1/jobs.</p>`;
    return;
  }
  root.innerHTML = `<ul class="job-list">${items
    .map((item) => {
      const state = item.state || "";
      const running = state === "queued" || state === "running";
      const stop = running ? `<button type="button" class="secondary" data-cancel="${escapeAttr(item.job_id)}">Stop</button>` : "";
      return `<li><button type="button" class="secondary" data-job="${escapeAttr(item.job_id)}">${escapeText(item.job_id || "")}</button> ${chip(stateChipKind(state), state || "unknown")} ${stop}</li>`;
    })
    .join("")}</ul>`;
}

async function cancelJob(jobId) {
  const posted = await call("POST", `/v1/jobs/${jobId}/cancel`);
  if (posted.payload && posted.payload.code && !posted.payload.job_id) {
    submitError.innerHTML = lateReveal(problemCallout(posted.payload));
  } else if (posted.payload && posted.payload.job_id) {
    renderStatus(posted.payload);
  }
  await refreshJobs();
}

async function checkFeasibility() {
  const root = document.getElementById("feasibility-root");
  if (!root) return;
  if (inputKind() !== "dna") {
    root.innerHTML = `<p class="note empty-state">This check needs DNA in the box. A protein or a GenBank file is turned into DNA on the server when you run, and the finished design says whether the rules were met.</p>`;
    return;
  }
  let body;
  try {
    const parsed = cleanedCds(document.getElementById("cds-text").value);
    if (!parsed.sequence || parsed.illegal.length || parsed.sequence.length % 3 !== 0) {
      throw new Error("Paste a DNA coding sequence first. The note under the box says what is wrong.");
    }
    body = {
      cds: parsed.sequence,
      constraints: constraintsForRequest(),
      translation_table: Number(document.getElementById("translation-table").value) || 11,
      five_prime_flank: flankValue("flank-5"),
      three_prime_flank: flankValue("flank-3"),
    };
  } catch (error) {
    root.innerHTML = `<div class="callout refusal"><p class="human">${escapeText(error.message)}</p><p class="note">Nothing was posted.</p></div>`;
    return;
  }
  root.innerHTML = `<p class="note">Checking the five feasibility tests.</p>`;
  const response = await call("POST", "/v1/feasibility", body);
  if (!response.payload || response.payload.code) {
    root.innerHTML = response.payload && response.payload.code
      ? problemCallout(response.payload)
      : `<p class="note empty-state">The feasibility check did not return a report.</p>`;
    return;
  }
  const report = response.payload;
  const findings = Array.isArray(report.findings) ? report.findings : [];
  const findingHtml = findings.length
    ? `<ul class="finding-list">${findings
        .map((item) => {
          const relax = item.suggested_relaxation ? ` ${item.suggested_relaxation}` : "";
          return `<li><strong>${escapeText(item.detail || "")}</strong> <span class="mono">${escapeText(item.code || "")}</span>${escapeText(relax)}</li>`;
        })
        .join("")}</ul>`
    : `<p class="field-help">No finding was returned.</p>`;
  const caveat = report.necessary_conditions_only
    ? "The report marks these checks as necessary conditions only, so a clear result is not a promise that a design exists."
    : "The finished candidate still carries constraint_report.feasible.";
  const verdict = report.feasible
    ? `${chip("feasible", "no finding")}<p class="field-help">These checks found nothing to report. ${escapeText(caveat)}</p>`
    : `${chip("infeasible", "a rule blocks every synonym")}<p class="field-help">At least one check found a rule the codons cannot satisfy. The sentence is the check’s own detail.</p>`;
  const checkPlain = {
    locked_span: "a rule in sequence the search cannot change",
    counting_bound: "a repeat limit this length cannot meet",
    per_position: "a codon with no allowed synonym",
    automaton_reachability: "motif rules that leave no way through",
    gc_attainability: "a GC window no synonymous sequence can meet",
  };
  const checks = Array.isArray(report.checks_run)
    ? report.checks_run
        .map((name) => (checkPlain[name] ? `${name} (${checkPlain[name]})` : name))
        .join("; ")
    : "";
  root.innerHTML = `${verdict}${findingHtml}${checks ? `<p class="note">Checks that ran: ${escapeText(checks)}.</p>` : ""}`;
}

function currentDna() {
  if (inputKind() !== "dna") {
    throw new Error("This check needs DNA in the box. A protein or a GenBank file is turned into DNA when you run a design, and it is not screened here.");
  }
  const parsed = cleanedCds(document.getElementById("cds-text").value);
  if (!parsed.sequence || parsed.illegal.length) {
    throw new Error("Paste a DNA sequence of A, C, G and T first. The note under the box says what is wrong.");
  }
  return parsed.sequence;
}

async function readUsageFile(file) {
  const text = await file.text();
  const table = Number(document.getElementById("translation-table").value) || 11;
  const response = await call("POST", "/v1/hosts:import", {
    text,
    subset: "user_supplied",
    translation_table: table,
    source: file.name,
  });
  if (!response.payload || response.payload.code || !response.payload.codon_usage) {
    const detail = response.payload && (response.payload.detail || response.payload.title);
    throw new Error(detail || "The usage file was not accepted.");
  }
  return response.payload;
}

async function runVendorScreen() {
  const root = document.getElementById("vendor-root");
  const cds = currentDna();
  root.innerHTML = `<p class="note">Screening synthesis-rule classes against this suite’s presets. This is not a manufacturer score.</p>`;
  const response = await call("POST", "/v1/vendor-screen", { cds });
  if (!response.payload || response.payload.code) {
    root.innerHTML = response.payload && response.payload.code
      ? problemCallout(response.payload)
      : `<p class="note empty-state">The vendor screen did not return a report.</p>`;
    return;
  }
  const report = response.payload;
  const checks = Array.isArray(report.checks) ? report.checks : [];
  const rows = checks.map((row) => {
    const mark = row.passed ? chip("feasible", "clear") : chip("infeasible", "flagged");
    return `<tr><td>${escapeText(row.set || "")}</td><td>${escapeText(row.kind || "")}</td><td>${mark}</td><td class="mono">${escapeText(String(row.count ?? ""))}</td><td>${escapeText(row.detail || "")}</td></tr>`;
  }).join("");
  root.innerHTML = `<div class="callout"><p class="human">${escapeText(report.note || "")}</p><p class="note">vendor_score is ${escapeText(String(report.vendor_score))}. Wet-lab expression is ${escapeText(String(report.wet_lab_expression))}. Flags: ${escapeText(String(report.n_flags ?? ""))}. Hairpins are not scored.</p></div><table><thead><tr><th>Set</th><th>Class</th><th>Result</th><th>Count</th><th>Detail</th></tr></thead><tbody>${rows}</tbody></table>`;
}

async function runGlobalFold() {
  const root = document.getElementById("fold-root");
  const cds = currentDna();
  root.innerHTML = `<p class="note">Folding the whole sequence. This is not the windowed accessibility objective.</p>`;
  const response = await call("POST", "/v1/structure:global", { cds, temperature_c: 37 });
  if (!response.payload || response.payload.code) {
    root.innerHTML = response.payload && response.payload.code
      ? problemCallout(response.payload)
      : `<p class="note empty-state">The fold did not return a report.</p>`;
    return;
  }
  const report = response.payload;
  const missing = report.value === null || report.value === undefined;
  const value = missing ? "not computed" : String(report.value);
  const unit = report.unit ? ` ${report.unit}` : "";
  root.innerHTML = `<div class="callout"><p class="human">${escapeText(report.label || "")}</p><p class="note">Engine ${escapeText(report.engine || "")}. ${escapeText(report.quantity || "quantity")} ${escapeText(value)}${escapeText(unit)}. Length ${escapeText(String(report.length_nt ?? ""))} nt. Wet-lab expression is ${escapeText(String(report.wet_lab_expression))}.</p></div>`;
}

async function loadSequenceFile(path, kind) {
  const response = await fetch(path);
  if (!response.ok) throw new Error(`${path} could not be loaded.`);
  document.getElementById("input-kind").value = kind;
  applyInputKind();
  document.getElementById("cds-text").value = (await response.text()).trim();
  refreshSequenceNote();
  syncRequest();
}

async function loadDemo() {
  await loadSequenceFile("/examples/rpsL.cds", "dna");
  if (!hostsById.size) await loadHosts();
  await selectInstalledHost("mg1655_ribosomal_usage.json", loadBundledProfile);
  document.getElementById("host-name").value = "escherichia_coli:511145/MG1655";
  document.getElementById("obj-cai").checked = true;
  document.getElementById("obj-gc").checked = false;
  document.getElementById("obj-violation").checked = false;
  document.getElementById("obj-structure").checked = false;
  const preset = document.getElementById("objective-preset");
  if (preset) preset.value = "cai";
  document.getElementById("budget").value = "fast";
  document.getElementById("generations").value = "4";
  document.getElementById("seed").value = "7";
  document.getElementById("constraint-mode").value = "example";
  await loadExampleConstraints();
  showConstraintMode();
  refreshSequenceNote();
  syncRequest();
  document.getElementById("submit").focus();
}

function bindGuide() {
  const kind = document.getElementById("input-kind");
  if (kind) kind.addEventListener("change", applyInputKind);
  const locus = document.getElementById("locus");
  if (locus) locus.addEventListener("input", () => {
    refreshSequenceNote();
    syncRequest();
  });
  const cds = document.getElementById("cds-text");
  if (cds) cds.addEventListener("input", () => {
    refreshSequenceNote();
    syncRequest();
  });
  const file = document.getElementById("cds-file");
  if (file) {
    file.addEventListener("change", () => {
      const chosen = file.files && file.files[0];
      if (!chosen) return;
      const reader = new FileReader();
      reader.onload = () => {
        const text = String(reader.result || "");
        cds.value = text;
        const name = (chosen.name || "").toLowerCase();
        if (name.endsWith(".gb") || name.endsWith(".gbk") || name.endsWith(".genbank") || (/LOCUS/i.test(text) && /ORIGIN/i.test(text))) {
          document.getElementById("input-kind").value = "construct";
        }
        applyInputKind();
      };
      reader.readAsText(chosen);
    });
  }
  document.getElementById("load-demo").addEventListener("click", () => {
    loadDemo().catch((error) => {
      submitError.innerHTML = `<div class="callout refusal"><p class="human">${escapeText(error.message)}</p></div>`;
    });
  });
  const gfp = document.getElementById("load-gfp");
  if (gfp) {
    gfp.addEventListener("click", () => {
      loadSequenceFile("/examples/BBa_E0040.fa", "dna")
        .then(() => {
          document.getElementById("trim-stop").checked = false;
          refreshSequenceNote();
          syncRequest();
        })
        .catch((error) => {
          submitError.innerHTML = `<div class="callout refusal"><p class="human">${escapeText(error.message)}</p></div>`;
        });
    });
  }
  const feasibility = document.getElementById("check-feasibility");
  if (feasibility) {
    feasibility.addEventListener("click", () => {
      checkFeasibility().catch((error) => {
        document.getElementById("feasibility-root").innerHTML = `<div class="callout refusal"><p class="human">${escapeText(error.message)}</p></div>`;
      });
    });
  }
  document.getElementById("load-profile").addEventListener("click", () => {
    selectInstalledHost("mg1655_ribosomal_usage.json", loadBundledProfile).catch((error) => {
      setProfileNote(error.message);
    });
  });
  document.getElementById("profile-file").addEventListener("change", (event) => {
    const chosen = event.target.files && event.target.files[0];
    if (!chosen) return;
    readUsageFile(chosen).then((profile) => {
      loadedProfile = profile;
      const subset = profile.codon_usage && profile.codon_usage.subset;
      const source = profile.codon_usage && profile.codon_usage.source;
      setProfileNote(`Loaded ${chosen.name}` + (subset ? ` (${subset})` : "") + (source ? `. ${source}` : "."));
      syncRequest();
    }).catch((error) => {
      setProfileNote(error.message);
    });
  });
  const sourceFile = document.getElementById("source-file");
  if (sourceFile) {
    sourceFile.addEventListener("change", (event) => {
      const chosen = event.target.files && event.target.files[0];
      if (!chosen) return;
      readUsageFile(chosen).then((profile) => {
        loadedSource = profile;
        const note = document.getElementById("source-note");
        if (note) note.textContent = `Source table loaded from ${chosen.name}. Harmonization uses these relative adaptiveness values for the input codons.`;
        syncRequest();
      }).catch((error) => {
        const note = document.getElementById("source-note");
        if (note) note.textContent = error.message;
      });
    });
  }
  const vendorButton = document.getElementById("vendor-screen");
  if (vendorButton) {
    vendorButton.addEventListener("click", () => {
      runVendorScreen().catch((error) => {
        document.getElementById("vendor-root").innerHTML = `<div class="callout refusal"><p class="human">${escapeText(error.message)}</p></div>`;
      });
    });
  }
  const foldButton = document.getElementById("global-fold");
  if (foldButton) {
    foldButton.addEventListener("click", () => {
      runGlobalFold().catch((error) => {
        document.getElementById("fold-root").innerHTML = `<div class="callout refusal"><p class="human">${escapeText(error.message)}</p></div>`;
      });
    });
  }
  const host = document.getElementById("host-select");
  host.addEventListener("change", applySelectedHost);
  document.getElementById("host-custom").addEventListener("input", syncRequest);
  const hostName = document.getElementById("host-name");
  if (hostName) hostName.addEventListener("input", syncRequest);
  const preset = document.getElementById("objective-preset");
  if (preset) preset.addEventListener("change", syncRequest);
  applySelectedHost();
  document.getElementById("constraint-mode").addEventListener("change", () => {
    if (document.getElementById("constraint-mode").value === "example") {
      loadExampleConstraints().catch((error) => {
        document.getElementById("cds-note").textContent = error.message;
      });
    }
    showConstraintMode();
  });
  document.getElementById("load-constraints").addEventListener("click", () => {
    document.getElementById("constraint-mode").value = "json";
    loadExampleConstraints()
      .then(() => {
        document.getElementById("constraint-json").readOnly = false;
      })
      .catch((error) => {
        submitError.innerHTML = `<div class="callout refusal"><strong>${escapeText(error.message)}</strong></div>`;
      });
  });
  document.querySelectorAll("[data-preset]").forEach((button) => {
    button.addEventListener("click", () => {
      loadConstraintPreset(button.dataset.preset).catch((error) => {
        submitError.innerHTML = `<div class="callout refusal"><strong>${escapeText(error.message)}</strong></div>`;
      });
    });
  });
  const kazusa = document.getElementById("load-kazusa");
  if (kazusa) {
    kazusa.addEventListener("click", () => {
      selectInstalledHost("kazusa_e_coli_316407.json", loadKazusaProfile).catch((error) => setProfileNote(error.message));
    });
  }
  const trim = document.getElementById("trim-stop");
  if (trim) {
    trim.addEventListener("change", () => {
      refreshSequenceNote();
      syncRequest();
    });
  }
  ["obj-cai", "obj-gc", "obj-violation", "obj-structure", "obj-harmonize", "obj-match"].forEach((id) => {
    const node = document.getElementById(id);
    if (!node) return;
    node.addEventListener("change", () => {
      const presetBox = document.getElementById("objective-preset");
      if (presetBox) presetBox.value = "custom";
    });
  });
  ["constraint-bsai", "constraint-bsmbi", "constraint-homo", "constraint-gc", "constraint-homo-run", "constraint-gc-window", "constraint-gc-min", "constraint-gc-max", "constraint-json", "obj-cai", "obj-gc", "obj-violation", "obj-structure", "obj-harmonize", "obj-match", "budget", "generations", "seed", "translation-table", "flank-5", "flank-3"].forEach((id) => {
    const node = document.getElementById(id);
    if (!node) return;
    node.addEventListener("input", syncRequest);
    node.addEventListener("change", syncRequest);
  });
  document.getElementById("edit-json").addEventListener("change", (event) => {
    requestInput.readOnly = !event.target.checked;
    if (!event.target.checked) syncRequest();
  });
  document.getElementById("load-primer-demo").addEventListener("click", () => {
    fetch("/examples/rpsL.cds")
      .then((response) => {
        if (!response.ok) throw new Error("The rpsL sequence could not be loaded.");
        return response.text();
      })
      .then((text) => {
        primerSequence.value = text.trim();
        document.getElementById("primer-flank").value = "40";
        selection = { start: 40, end: 300 };
        drawViewer();
      })
      .catch((error) => {
        document.getElementById("primer-error").innerHTML = `<div class="callout refusal"><strong>${escapeText(error.message)}</strong></div>`;
      });
  });
  document.getElementById("select-all").addEventListener("click", () => {
    const sequence = cleanedTemplate(primerSequence.value).sequence;
    if (!sequence) {
      document.getElementById("primer-error").innerHTML = `<div class="callout refusal"><strong>Paste a template first.</strong> <span class="note">Or load the rpsL sequence.</span></div>`;
      return;
    }
    selection = { start: 0, end: sequence.length };
    drawViewer();
  });
  ["region-start", "region-end"].forEach((id) => {
    document.getElementById(id).addEventListener("change", selectionFromFields);
  });
  const jobsRoot = document.getElementById("jobs-root");
  if (jobsRoot) {
    jobsRoot.addEventListener("click", (event) => {
      const cancel = event.target.closest("[data-cancel]");
      if (cancel) {
        cancelJob(cancel.dataset.cancel);
        return;
      }
      const open = event.target.closest("[data-job]");
      if (!open) return;
      history.replaceState(null, "", `#${open.dataset.job}`);
      showJob(open.dataset.job);
    });
  }
  const checkSeed = document.getElementById("check-seed");
  if (checkSeed) checkSeed.addEventListener("click", () => {
    importSeed().catch((error) => {
      document.getElementById("seed-note").textContent = error.message;
    });
  });
  const importUsage = document.getElementById("import-usage");
  if (importUsage) importUsage.addEventListener("click", () => {
    importUsageTable().catch((error) => {
      document.getElementById("usage-note").textContent = error.message;
    });
  });
  const usageFile = document.getElementById("usage-file");
  if (usageFile) usageFile.addEventListener("change", () => {
    const chosen = usageFile.files && usageFile.files[0];
    if (!chosen) return;
    const reader = new FileReader();
    reader.onload = () => {
      document.getElementById("usage-text").value = String(reader.result || "");
    };
    reader.readAsText(chosen);
  });
  const saveMeasurement = document.getElementById("save-measurement");
  if (saveMeasurement) saveMeasurement.addEventListener("click", () => {
    saveLabMeasurement().catch((error) => {
      document.getElementById("measure-note").textContent = error.message;
    });
  });
  const exportMeasurements = document.getElementById("export-measurements");
  if (exportMeasurements) exportMeasurements.addEventListener("click", () => {
    downloadMeasurements().catch((error) => {
      document.getElementById("measure-note").textContent = error.message;
    });
  });
  const measureDate = document.getElementById("measure-date");
  if (measureDate && !measureDate.value) {
    measureDate.value = new Date().toISOString().slice(0, 10);
  }
  showConstraintMode();
  refreshSequenceNote();
  refreshOptimizedHandoff();
  loadHosts().catch((error) => setProfileNote(error.message));
  loadModels().catch((error) => {
    const note = document.getElementById("model-status");
    if (note) note.textContent = error.message;
  });
  loadEngineStatus();
  loadVersion();
  refreshJobs().catch(() => {});
}

async function importSeed() {
  const cds = cleanedCds(document.getElementById("cds-text").value).sequence;
  const starting = document.getElementById("starting-cds").value;
  const result = await call("POST", "/v1/seeds:import", {
    tool: document.getElementById("seed-tool").value,
    payload: starting,
    target_cds: cds,
    translation_table: Number(document.getElementById("translation-table").value) || 11,
    host_id: selectedHost(),
  });
  const note = document.getElementById("seed-note");
  if (result.status >= 400) {
    note.textContent = `${result.payload.code}: ${result.payload.detail}`;
    return;
  }
  const seed = (result.payload.seeds || [])[0];
  const acceptance = seed && seed.acceptance;
  if (!acceptance || !acceptance.accepted) {
    const diffs = (acceptance && acceptance.protein_diff) || [];
    const where = diffs.map((item) => `${item.position}:${item.expected}->${item.found}`).join(", ");
    note.textContent = `Rejected (${acceptance ? acceptance.reason : "empty"}). ${where}`;
    return;
  }
  note.textContent = "Accepted. The protein matches. Run design to start from this CDS.";
}

async function importUsageTable() {
  const text = document.getElementById("usage-text").value;
  const result = await call("POST", "/v1/hosts:import", {
    text,
    format: "auto",
    label: document.getElementById("usage-label").value,
    translation_table: Number(document.getElementById("translation-table").value) || 11,
    organism: selectedHost(),
  });
  const note = document.getElementById("usage-note");
  if (result.status >= 400) {
    note.textContent = `${result.payload.code}: ${result.payload.detail}`;
    return;
  }
  loadedProfile = result.payload;
  const missing = result.payload.unobserved_families || [];
  const extra = missing.length ? ` Unobserved families: ${missing.join(", ")}. CAI will refuse this table until those families have a count.` : "";
  setProfileNote(`Imported a user-supplied table.${extra}`);
  note.textContent = result.payload.codon_usage.source;
  syncRequest();
}

async function saveLabMeasurement() {
  const cds = document.getElementById("measure-cds").value || document.getElementById("cds-text").value;
  const date = document.getElementById("measure-date").value;
  const result = await call("POST", "/v1/measurements", {
    cds,
    host_id: document.getElementById("measure-host").value,
    measures: document.getElementById("measure-kind").value,
    assay: document.getElementById("measure-assay").value,
    value: Number(document.getElementById("measure-value").value),
    unit: document.getElementById("measure-unit").value,
    condition_id: document.getElementById("measure-condition").value,
    control_role: document.getElementById("measure-control").value,
    censoring: document.getElementById("measure-censoring").value,
    measured_on: date ? `${date}T00:00:00Z` : undefined,
  });
  const note = document.getElementById("measure-note");
  if (result.status >= 400) {
    note.textContent = `${result.payload.code}: ${result.payload.detail}`;
    return;
  }
  note.textContent = result.payload.statement || "Stored. No sequence model was trained.";
}

async function downloadMeasurements() {
  const headers = {};
  const token = tokenInput.value.trim();
  if (token) headers.Authorization = `Bearer ${token}`;
  const response = await fetch("/v1/observations?format=csv", { headers });
  const text = await response.text();
  if (!response.ok) {
    throw new Error(text);
  }
  const blob = new Blob([text], { type: "text/csv" });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = "measurements.csv";
  link.click();
  URL.revokeObjectURL(url);
  document.getElementById("measure-note").textContent = "Downloaded measurements.csv.";
}

bindGuide();
