/* Renders job_status and pareto_front documents. Numbers come from the API.
   Polls GET /v1/jobs/{id} while a job is running. Codon tints, when shown,
   use relative_adaptiveness from the submitted host profile. The front
   document's objective_vector, unknown_terms, and objective_sense are per
   term, not per codon, so this page does not invent a position mask. */

const tokenInput = document.getElementById("token");
const requestInput = document.getElementById("request");
const submitButton = document.getElementById("submit");
const submitError = document.getElementById("submit-error");
const statusRoot = document.getElementById("status-root");
const frontRoot = document.getElementById("front-root");

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
  const href = (problem && problem.type) || `https://codonforge.org/problems/${String(code).toLowerCase()}`;
  return `<a href="${escapeAttr(href)}">${escapeText(code)}</a>`;
}

function problemCallout(problem) {
  if (!problem || !problem.code) {
    return `<div class="callout refusal"><strong>No CF-* code.</strong> <span class="note">The response did not name a registered code.</span></div>`;
  }
  return `<div class="callout refusal"><strong>${codeLink(problem)}</strong> ${escapeText(problem.title || "")}<div class="note">${escapeText(problem.detail || "")}</div></div>`;
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
      const rows = group
        .map(
          (warning) =>
            `<div class="warn-row"><div class="code">${codeLink(warning)}</div><div class="mono">${escapeText(warning.stage || "")}</div><div class="mono">${escapeText(warning.subject || "")}</div><div><strong>${escapeText(warning.title || "")}</strong><div class="note">${escapeText(warning.detail || "")}</div></div></div>`,
        )
        .join("");
      const kind = refusal ? "refusal" : "notice";
      const label = refusal ? "hard refusal" : "notice";
      return `<section class="warn-group ${kind}${open ? " open" : ""}" data-stage="${escapeAttr(stage)}">
        <button type="button" class="warn-toggle" aria-expanded="${open ? "true" : "false"}">${escapeText(stage)} · ${group.length} · ${escapeText(label)}${codes ? ` · ${escapeText(codes)}` : ""}</button>
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
      return `<div class="tier-track"><span>${escapeText(tier.tier)}</span><span>${escapeText(tier.state)}</span><div class="bar"><span style="width:${start}%" data-target="${width}"></span></div><span>${escapeText(generation)}</span></div>`;
    })
    .join("");
  const objectives = (document.resolved_objectives || [])
    .map(
      (term) =>
        `<tr><td class="mono">${escapeText(term.term_id)}</td><td>${roleChip(term.role)}</td><td class="mono">${escapeText(term.direction)}</td></tr>`,
    )
    .join("");
  const failed = document.state === "failed" ? problemCallout(document.problem) : document.problem && document.problem.code ? problemCallout(document.problem) : "";
  const stateKind = document.state === "succeeded" ? "feasible" : document.state === "failed" ? "refusal" : "warning";
  statusRoot.innerHTML = `
    <div class="rail"><div class="k">Job</div><div class="v mono">${escapeText(document.job_id || "")}</div></div>
    <div class="rail"><div class="k">State</div><div class="v">${chip(stateKind, document.state || "")}</div></div>
    ${failed}
    <h2>Tiers</h2>
    ${tiers || '<p class="note">No tier has been reported.</p>'}
    <h2>Resolved objectives</h2>
    <div class="hscroll"><table class="data"><thead><tr><th>term</th><th>role</th><th>direction</th></tr></thead><tbody>${objectives}</tbody></table></div>
    <h2>Warnings</h2>
    ${renderWarnings(document)}
  `;
  bindWarnings();
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
      ? { code: match[0], type: `https://codonforge.org/problems/${match[0].toLowerCase()}`, title: "The front has no candidates", detail: reason }
      : null;
    frontRoot.innerHTML = `<div class="reveal">${problem ? problemCallout(problem) : `<div class="callout refusal"><strong>Empty front.</strong> <span class="note">No candidates were returned${reason ? `: ${escapeText(reason)}` : ""}. The document did not name a CF-* code for this outcome.</span></div>`}</div>`;
    return;
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
  frontRoot.innerHTML = `
    <div class="reveal">
    <p class="note">reduction.method <code>${escapeText((front.reduction || {}).method || "")}</code> · presented ${(front.reduction || {}).presented_count ?? ""} of ${(front.reduction || {}).raw_nondominated_count ?? ""}</p>
    <div class="legend">${chip("feasible", "feasible")}${chip("infeasible", "infeasible")}${chip("unknown", "UNKNOWN")}${chip("recommended", "recommended")}${chip("structure", "structure-optimal among scored")}${chip("diagnostic", "diagnostic-only")}</div>
    <div class="callout pulse-once"><strong>Recommended</strong> <code>${escapeText(recommended.candidate_id || "")}</code> · rule <code>${escapeText(recommended.rule || "")}</code><br>informed_by ${informed || "—"}<br>not_informed_by ${excluded || "—"}</div>
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
  frontRoot.innerHTML = `<div class="reveal">${body}</div>`;
}

async function loadFront(jobId) {
  const front = await call("GET", `/v1/jobs/${jobId}/front?reduction=default&max_candidates=20`);
  if (front.status === 200 && front.payload && !front.payload.code) renderFront(front.payload);
  else showFrontFailure(front);
}

let pollToken = 0;

async function poll(jobId) {
  const token = ++pollToken;
  try {
    let status = await call("GET", `/v1/jobs/${jobId}`);
    while (
      token === pollToken &&
      status.payload &&
      status.payload.job_id &&
      ["queued", "running"].includes(status.payload.state)
    ) {
      renderStatus(status.payload);
      await new Promise((resolve) => setTimeout(resolve, POLL_MS));
      if (token !== pollToken) return;
      status = await call("GET", `/v1/jobs/${jobId}`);
    }
    if (token !== pollToken || !status.payload) return;
    if (status.payload.code && !status.payload.job_id) {
      submitError.innerHTML = problemCallout(status.payload);
      frontRoot.innerHTML = `<div class="reveal">${problemCallout(status.payload)}</div>`;
      return;
    }
    renderStatus(status.payload);
    if (status.payload.state !== "succeeded") {
      const problem = status.payload.problem || { code: "", title: "", detail: "" };
      frontRoot.innerHTML = `<div class="reveal">${status.payload.problem ? problemCallout(problem) : `<div class="callout refusal"><strong>Job ${escapeText(status.payload.state || "unfinished")}.</strong> <span class="note">No front was requested, because the job did not succeed. The status document did not include a problem code.</span></div>`}</div>`;
      return;
    }
    await loadFront(jobId);
  } catch (error) {
    if (token !== pollToken) return;
    submitError.innerHTML = `<div class="callout refusal"><strong>The status request did not complete.</strong> <span class="note">${escapeText(error)}</span></div>`;
  }
}

async function showJob(jobId) {
  const status = await call("GET", `/v1/jobs/${jobId}`);
  if (!status.payload) return;
  if (status.payload.code && !status.payload.job_id) {
    submitError.innerHTML = problemCallout(status.payload);
    frontRoot.innerHTML = `<div class="reveal">${problemCallout(status.payload)}</div>`;
    return;
  }
  renderStatus(status.payload);
  if (status.payload.state === "queued" || status.payload.state === "running") {
    poll(jobId);
    return;
  }
  if (status.payload.state !== "succeeded") {
    frontRoot.innerHTML = `<div class="reveal">${status.payload.problem ? problemCallout(status.payload.problem) : `<div class="callout refusal"><strong>Job ${escapeText(status.payload.state || "unfinished")}.</strong> <span class="note">No front is available.</span></div>`}</div>`;
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

frontRoot.addEventListener("click", (event) => {
  const button = event.target.closest("[data-copy]");
  if (!button) return;
  copyRaw(button.dataset.copy || "", button);
});

submitButton.addEventListener("click", async () => {
  submitError.innerHTML = "";
  let body;
  try {
    body = JSON.parse(requestInput.value);
  } catch (error) {
    submitError.innerHTML = `<div class="callout refusal"><strong>The request is not JSON.</strong> <span class="note">Nothing was posted.</span></div>`;
    return;
  }
  lastSubmitted = body;
  try {
    const posted = await call("POST", "/v1/jobs", body);
    if (!posted.payload) return;
    if (posted.payload.code && !posted.payload.job_id) {
      submitError.innerHTML = problemCallout(posted.payload);
      frontRoot.innerHTML = `<div class="reveal">${problemCallout(posted.payload)}</div>`;
      return;
    }
    history.replaceState(null, "", `#${posted.payload.job_id}`);
    renderStatus(posted.payload);
    poll(posted.payload.job_id);
  } catch (error) {
    submitError.innerHTML = `<div class="callout refusal"><strong>The request did not complete.</strong> <span class="note">${escapeText(error)}</span></div>`;
  }
});

window.optimizedCds = "";

function refreshOptimizedHandoff() {
  const slot = document.getElementById("use-optimized");
  if (!slot) return;
  slot.innerHTML = window.optimizedCds
    ? `<button type="button" id="use-optimized-sequence">Use optimized sequence</button>`
    : "";
}

document.querySelectorAll("button.tab").forEach((button) => {
  button.addEventListener("click", () => {
    document.querySelectorAll("button.tab").forEach((item) => item.classList.toggle("active", item === button));
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
  const sequence = primerSequence.value.toUpperCase().replace(/[^ACGT]/g, "");
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
      ? `Region ${selection.start}–${selection.end} (${selection.end - selection.start} nt).`
      : "No region selected.";
  }
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
  const node = event.target.closest("[data-index]");
  return node ? Number(node.dataset.index) : null;
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

function detailPanel(pair) {
  const flags = clientFlags(pair)
    .map((flag) => `<div class="zone-chip ${zoneClass(flag)}">${escapeText(flag.property)} ${escapeText(flag.observed)} · ${escapeText(flag.threshold)}</div>`)
    .join("");
  return `<div class="pair-detail">
    <p><strong>Forward</strong> <span class="mono">${escapeText(pair.forward_sequence)}</span></p>
    <p>GC ${escapeText(pair.gc_forward)}% · Tm ${escapeText(Number(pair.tm_forward).toFixed(2))} °C
      <button type="button" class="param-help" title="${escapeAttr(pair.parameter_citation || "")}">?</button></p>
    <p class="note">${escapeText(pair.parameter_set_id || "")} ${escapeText(pair.parameter_digest || "")}</p>
    <pre class="mono mono-reveal">${escapeText(pair.hairpin_diagram_forward || "")}</pre>
    <p>Hairpin ΔG ${escapeText(pair.hairpin_dG_forward)} · self-dimer ΔG ${escapeText(pair.self_dimer_dG_forward)} · 3′ clamp ${pair.gc_clamp_forward ? "G or C" : "not G or C"}</p>
    <p><strong>Reverse</strong> <span class="mono">${escapeText(pair.reverse_sequence)}</span></p>
    <pre class="mono mono-reveal">${escapeText(pair.hairpin_diagram_reverse || "")}</pre>
    <p>Hairpin ΔG ${escapeText(pair.hairpin_dG_reverse)} · self-dimer ΔG ${escapeText(pair.self_dimer_dG_reverse)} · 3′ clamp ${pair.gc_clamp_reverse ? "G or C" : "not G or C"}</p>
    <p><strong>Pair</strong> hetero-dimer ΔG ${escapeText(pair.hetero_dimer_dG)} · 3′ ΔG ${escapeText(pair.hetero_dimer_3prime_dG)} · Tm mismatch ${escapeText(pair.tm_mismatch)} · product ${escapeText(pair.product_size)}</p>
    <pre class="mono mono-reveal">${escapeText(pair.hetero_dimer_alignment || "")}</pre>
    <div class="flag-host">${flags || `<p class="zone-chip zone-green">No property is in the amber or red zone.</p>`}</div>
    <p class="note">Specificity against genome not checked.</p>
    ${pair.product_sequence ? `<button type="button" data-copy="${escapeAttr(pair.product_sequence)}">Copy product sequence</button><div class="seq-block">${escapeText(pair.product_sequence)}</div>` : ""}
  </div>`;
}

function renderPrimerResult(result) {
  const root = document.getElementById("primer-results");
  lastPrimerResult = result;
  const pairs = (result && result.pairs) || [];
  if (!pairs.length) {
    const detail = (result && result.detail) || "No filter counts were returned with this empty result.";
    root.innerHTML = problemCallout({ code: "CF-PRIM-001", title: "No primer pairs survived pruning", detail });
    return;
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

async function pollPrimers(jobId) {
  const bar = document.getElementById("primer-progress");
  let status = await call("GET", `/v1/primers/${jobId}`);
  while (status.payload && ["queued", "running"].includes(status.payload.state)) {
    bar.innerHTML = primerProgress(status.payload);
    await new Promise((resolve) => setTimeout(resolve, 2000));
    status = await call("GET", `/v1/primers/${jobId}`);
  }
  if (!status.payload) return;
  bar.innerHTML = primerProgress(status.payload);
  if (status.payload.problem) {
    document.getElementById("primer-error").innerHTML = problemCallout(status.payload.problem);
    document.getElementById("primer-results").innerHTML = problemCallout(status.payload.problem);
    return;
  }
  renderPrimerResult(status.payload.primer_result);
}

document.getElementById("primer-search").addEventListener("click", async () => {
  document.getElementById("primer-error").innerHTML = "";
  document.getElementById("primer-results").innerHTML = "";
  if (!selection) {
    document.getElementById("primer-error").innerHTML = `<div class="callout refusal"><strong>No region selected.</strong> <span class="note">Drag across the sequence viewer.</span></div>`;
    return;
  }
  const sequence = primerSequence.value.toUpperCase().replace(/[^ACGT]/g, "");
  const posted = await call("POST", "/v1/primers/candidates", {
    sequence,
    region: selection,
    config: primerConfig(),
  });
  if (!posted.payload) return;
  if (posted.payload.code && !posted.payload.job_id) {
    document.getElementById("primer-error").innerHTML = problemCallout(posted.payload);
    document.getElementById("primer-results").innerHTML = problemCallout(posted.payload);
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
