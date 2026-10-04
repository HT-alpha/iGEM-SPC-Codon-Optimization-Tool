/* Renders job_status and pareto_front documents. Numbers come from the API. */

const tokenInput = document.getElementById("token");
const requestInput = document.getElementById("request");
const submitButton = document.getElementById("submit");
const submitError = document.getElementById("submit-error");
const statusRoot = document.getElementById("status-root");
const frontRoot = document.getElementById("front-root");

function escapeText(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;");
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

function renderStatus(document) {
  const progress = document.progress || { tiers: [] };
  const tiers = (progress.tiers || [])
    .map((tier) => {
      const fraction = typeof tier.fraction === "number" ? tier.fraction : null;
      const width = fraction === null ? 0 : Math.max(0, Math.min(1, fraction)) * 100;
      const generation =
        tier.generation === undefined || tier.generation === null
          ? ""
          : `${tier.generation}/${tier.max_generations ?? "—"}`;
      return `<div class="tier-track"><span>${escapeText(tier.tier)}</span><span>${escapeText(tier.state)}</span><div class="bar"><span style="width:${width}%"></span></div><span>${escapeText(generation)}</span></div>`;
    })
    .join("");
  const objectives = (document.resolved_objectives || [])
    .map(
      (term) =>
        `<tr><td class="mono">${escapeText(term.term_id)}</td><td>${roleChip(term.role)}</td><td class="mono">${escapeText(term.direction)}</td></tr>`,
    )
    .join("");
  const warnings = (document.warnings || [])
    .map(
      (warning) =>
        `<div class="warn-row"><div class="code">${escapeText(warning.code || "")}</div><div class="mono">${escapeText(warning.stage || "")}</div><div class="mono">${escapeText(warning.subject || "")}</div><div><strong>${escapeText(warning.title || "")}</strong><div class="note">${escapeText(warning.detail || "")}</div></div></div>`,
    )
    .join("");
  const warningBlock = warnings
    ? `<div class="warn-row head"><div>code</div><div>stage</div><div>subject</div><div>title</div></div>${warnings}`
    : '<p class="note">No warnings on this status.</p>';
  const stateKind = document.state === "succeeded" ? "feasible" : document.state === "failed" ? "refusal" : "warning";
  statusRoot.innerHTML = `
    <div class="rail"><div class="k">Job</div><div class="v mono">${escapeText(document.job_id || "")}</div></div>
    <div class="rail"><div class="k">State</div><div class="v">${chip(stateKind, document.state || "")}</div></div>
    <h2>Tiers</h2>
    ${tiers || '<p class="note">No tier has been reported.</p>'}
    <h2>Resolved objectives</h2>
    <table class="data"><thead><tr><th>term</th><th>role</th><th>direction</th></tr></thead><tbody>${objectives}</tbody></table>
    <h2>Warnings</h2>
    ${warningBlock}
  `;
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
  const klass = structure ? "callout structure" : "callout";
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

function renderFront(front) {
  const recommended = front.recommended || {};
  const order = front.objective_order || [];
  const senses = front.objective_sense || [];
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
  const rows = (front.candidates || [])
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
      return `<tr class="${klass}"><td class="mono">${escapeText(candidate.candidate_id)}</td>${feasibleCell(candidate)}<td class="mono">${escapeText(unknown)}</td>${cells}<td>${extra}</td></tr>`;
    })
    .join("");
  const plot = renderPlot(front, recommended.candidate_id, optimalIds);
  frontRoot.innerHTML = `
    <p class="note">reduction.method <code>${escapeText((front.reduction || {}).method || "")}</code> · presented ${(front.reduction || {}).presented_count ?? ""} of ${(front.reduction || {}).raw_nondominated_count ?? ""}</p>
    <div class="legend">${chip("feasible", "feasible")}${chip("infeasible", "infeasible")}${chip("unknown", "UNKNOWN")}${chip("recommended", "recommended")}${chip("structure", "structure-optimal among scored")}${chip("diagnostic", "diagnostic-only")}</div>
    <div class="callout"><strong>Recommended</strong> <code>${escapeText(recommended.candidate_id || "")}</code> · rule <code>${escapeText(recommended.rule || "")}</code><br>informed_by ${informed || "—"}<br>not_informed_by ${excluded || "—"}</div>
    ${optimalText}
    ${plot}
    <p class="note">A column in this table is a Pareto axis. A diagnostic term is in the last column and is not an axis. UNKNOWN is the word in the hatched cell, not a number.</p>
    <table class="data"><thead><tr><th>candidate</th><th>constraint_report.feasible</th><th>unknown_terms</th>${order
      .map((term) => {
        const role = roleFor(front, term);
        return `<th>${escapeText(term)}${role ? ` ${roleChip(role)}` : ""}</th>`;
      })
      .join("")}<th>other terms</th></tr></thead><tbody>${rows}</tbody></table>
  `;
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
  const marks = points
    .map((point) => {
      const id = point.candidate.candidate_id;
      const fill = id === recommendedId ? "#1a4a73" : optimalIds.has(id) ? "#3a5f4d" : "#1c1917";
      return `<rect x="${scaleX(point.x) - 4}" y="${scaleY(point.y) - 4}" width="8" height="8" fill="${fill}"/>`;
    })
    .join("");
  const rail = unknown.length
    ? `<rect x="20" y="250" width="680" height="40" fill="url(#hatch)"/>${unknown
        .map((point, index) => {
          const xPos = point.x === null ? 40 + index * 16 : scaleX(point.x);
          return `<rect x="${xPos - 4}" y="266" width="8" height="8" fill="#1c1917"/>`;
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
  </svg>${caption}`;
}

async function poll(jobId) {
  const deadline = Date.now() + 60000;
  let status = await call("GET", `/v1/jobs/${jobId}`);
  while (status.payload && ["queued", "running"].includes(status.payload.state) && Date.now() < deadline) {
    renderStatus(status.payload);
    await new Promise((resolve) => setTimeout(resolve, 400));
    status = await call("GET", `/v1/jobs/${jobId}`);
  }
  if (status.payload) renderStatus(status.payload);
  if (!status.payload || status.payload.state !== "succeeded") return;
  const front = await call(
    "GET",
    `/v1/jobs/${jobId}/front?reduction=default&max_candidates=20`,
  );
  if (front.status === 200) renderFront(front.payload);
  else submitError.textContent = `${front.status} ${front.payload && front.payload.code}`;
}

async function showJob(jobId) {
  const status = await call("GET", `/v1/jobs/${jobId}`);
  if (!status.payload) return;
  if (status.payload.code) {
    submitError.textContent = `${status.status} ${status.payload.code}: ${status.payload.detail || ""}`;
    return;
  }
  renderStatus(status.payload);
  if (status.payload.state !== "succeeded") return;
  const front = await call(
    "GET",
    `/v1/jobs/${jobId}/front?reduction=default&max_candidates=20`,
  );
  if (front.status === 200) renderFront(front.payload);
  else submitError.textContent = `${front.status} ${front.payload && front.payload.code}`;
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
  if (payload.state === "fitted") {
    note.textContent =
      "epistasis.context_naturalness fitted. W=" +
      payload.window +
      ". Held-out " +
      payload.held_out_nats_per_codon +
      " nats/codon. Codon-usage baseline " +
      payload.baseline_nats_per_codon +
      " nats/codon.";
    return;
  }
  note.textContent = "epistasis.context_naturalness absent";
}

loadCalibration();

window.addEventListener("hashchange", jobFromHash);
jobFromHash();

submitButton.addEventListener("click", async () => {
  submitError.textContent = "";
  let body;
  try {
    body = JSON.parse(requestInput.value);
  } catch (error) {
    submitError.textContent = "The request is not JSON.";
    return;
  }
  try {
    const posted = await call("POST", "/v1/jobs", body);
    if (!posted.payload) return;
    if (posted.payload.code) {
      submitError.textContent = `${posted.status} ${posted.payload.code}: ${posted.payload.detail || ""}`;
      return;
    }
    history.replaceState(null, "", `#${posted.payload.job_id}`);
    renderStatus(posted.payload);
    poll(posted.payload.job_id);
  } catch (error) {
    submitError.textContent = String(error);
  }
});
