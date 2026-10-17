import { mkdirSync, readFileSync, writeFileSync } from "node:fs";

const DEBUG = "http://127.0.0.1:9222";
const OUT = "/workspace/benchmarks/validation/results/primer_page.json";

function sleep(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

function fastaSequence(name) {
  const lines = readFileSync(
    `/workspace/benchmarks/validation/sequences/${name}.fa`,
    "utf8",
  ).split("\n");
  return lines.slice(1).join("").trim();
}

const gfp = fastaSequence("BBa_E0040_onestop");
const version = await fetch(`${DEBUG}/json/version`).then((response) => response.json());
const ws = new WebSocket(version.webSocketDebuggerUrl);
await new Promise((resolve) => ws.addEventListener("open", resolve, { once: true }));
let seq = 0;
const pending = new Map();
ws.addEventListener("message", (event) => {
  const message = JSON.parse(event.data);
  if (message.id && pending.has(message.id)) pending.get(message.id)(message);
});
function send(method, params = {}, sessionId) {
  const id = ++seq;
  return new Promise((resolve, reject) => {
    const timer = setTimeout(() => reject(new Error(`timeout ${method}`)), 30000);
    pending.set(id, (message) => {
      clearTimeout(timer);
      if (message.error) reject(new Error(JSON.stringify(message.error)));
      else resolve(message.result);
    });
    const payload = { id, method, params };
    if (sessionId) payload.sessionId = sessionId;
    ws.send(JSON.stringify(payload));
  });
}

const { targetId } = await send("Target.createTarget", { url: "about:blank" });
const { sessionId } = await send("Target.attachToTarget", { targetId, flatten: true });
await send("Page.enable", {}, sessionId);
await send("Runtime.enable", {}, sessionId);

mkdirSync("/opt/cursor/artifacts", { recursive: true });

async function shot(name) {
  const { data } = await send("Page.captureScreenshot", { format: "png" }, sessionId);
  writeFileSync(`/opt/cursor/artifacts/${name}.png`, Buffer.from(data, "base64"));
}

async function evaluate(expression) {
  const result = await send(
    "Runtime.evaluate",
    { expression, awaitPromise: true, returnByValue: true },
    sessionId,
  );
  if (result.exceptionDetails) {
    throw new Error(JSON.stringify(result.exceptionDetails));
  }
  return result.result.value;
}

await send(
  "Emulation.setDeviceMetricsOverride",
  { width: 1280, height: 900, deviceScaleFactor: 1, mobile: false },
  sessionId,
);
await send("Page.navigate", { url: "http://127.0.0.1:8000/" }, sessionId);
for (let i = 0; i < 30; i += 1) {
  const ready = await evaluate("document.getElementById('profile-note') ? document.getElementById('profile-note').innerText : ''");
  if (ready && !ready.includes("No codon table")) break;
  await sleep(200);
}

const report = {};
report.no_region = await evaluate(`
  document.getElementById("tab-btn-primer").click();
  document.getElementById("primer-search").click();
  document.getElementById("primer-error").innerText.trim()
`);

const pasted = `>GFP demo\n${gfp.toLowerCase()}`;
report.page_search = await evaluate(`
  const sequence = ${JSON.stringify(pasted)};
  const box = document.getElementById("primer-sequence");
  box.value = sequence;
  box.dispatchEvent(new Event("input"));
  document.getElementById("primer-flank").value = "80";
  document.getElementById("primer-pmin").value = "150";
  document.getElementById("primer-pmax").value = "400";
  document.getElementById("region-start").value = "201";
  document.getElementById("region-end").value = "400";
  document.getElementById("region-start").dispatchEvent(new Event("change"));
  document.getElementById("primer-search").click();
  ({
    bases: document.querySelectorAll("#primer-viewer .nt").length,
    first: document.querySelector("#primer-viewer .nt") ? document.querySelector("#primer-viewer .nt").textContent : "",
    selection: document.getElementById("primer-selection").innerText,
  })
`);

let text = "";
for (let i = 0; i < 40; i += 1) {
  text = await evaluate(`(() => {
    const message = document.getElementById("primer-error").innerText.trim();
    const results = document.getElementById("primer-results").innerText.trim();
    return JSON.stringify({ error: message, results: results.slice(0, 900), done: results.includes("Forward") || message.includes("did not complete") || message.includes("CF-") });
  })()`);
  const parsed = JSON.parse(text);
  if (parsed.done) {
    report.result = parsed;
    break;
  }
  await sleep(250);
}
if (!report.result) report.result = JSON.parse(text);

report.illegal_base = await evaluate(`(() => {
  const oligo = document.getElementById("primer-sequence");
  oligo.value = ">GFP demo\\nATGNCCC";
  oligo.dispatchEvent(new Event("input"));
  return document.getElementById("primer-viewer").innerText.trim();
})()`);
await shot("primer-illegal-base");

report.codon = await evaluate(`(() => {
  document.getElementById("tab-btn-codon").click();
  const help = document.getElementById("help-constraints").innerText;
  const presets = [...document.querySelectorAll("[data-preset]")].map((button) => button.dataset.preset);
  const kazusa = Boolean(document.getElementById("load-kazusa"));
  const trim = document.getElementById("trim-stop");
  const box = document.getElementById("cds-text");
  const tail = "ATG" + "GCT".repeat(237) + "TAATAA";
  box.value = tail;
  trim.checked = true;
  trim.dispatchEvent(new Event("change"));
  const note = document.getElementById("cds-note").innerText;
  const request = document.getElementById("request").value;
  let sent = "";
  try { sent = JSON.parse(request).cds || ""; } catch (error) { sent = ""; }
  return {
    help,
    presets,
    kazusa,
    note,
    untrimmed: tail.length,
    sent: sent.length,
  };
})()`);
await evaluate("document.getElementById('step-seq').scrollIntoView()");
await shot("trim-at-first-stop");

report.failed_design = await evaluate(`(async () => {
  document.getElementById("tab-btn-codon").click();
  document.getElementById("cds-text").value = "ATG" + "GGG".repeat(12) + "TAA";
  document.getElementById("trim-stop").checked = false;
  document.getElementById("trim-stop").dispatchEvent(new Event("change"));
  document.getElementById("constraint-mode").value = "json";
  document.getElementById("constraint-mode").dispatchEvent(new Event("change"));
  document.getElementById("constraint-json").value = JSON.stringify({
    schema_version: "1.0.0",
    constraint_set_id: "cst_impossible_gc",
    constraints: [{
      id: "c_gc",
      class: "blocking",
      scope: "cds",
      rule: {kind: "gc_window", window_nt: 12, min: 0, max: 0.1}
    }]
  });
  document.getElementById("constraint-json").dispatchEvent(new Event("input"));
  document.getElementById("budget").value = "fast";
  document.getElementById("edit-json").checked = false;
  document.getElementById("submit").click();
  for (let i = 0; i < 20; i += 1) {
    await new Promise((resolve) => setTimeout(resolve, 500));
    const status = document.getElementById("status-root").innerText;
    const front = document.getElementById("front-root").innerText;
    if (status.includes("CF-SRCH-006") || front.includes("CF-SRCH-006")) {
      return {status: status.slice(0, 900), front: front.slice(0, 900)};
    }
  }
  return {
    status: document.getElementById("status-root").innerText.slice(0, 900),
    front: document.getElementById("front-root").innerText.slice(0, 900),
  };
})()`);
await evaluate("document.getElementById('status').scrollIntoView()");
await shot("failed-design-cf-srch-006");

ws.close();
writeFileSync(OUT, JSON.stringify(report, null, 2));
console.log(JSON.stringify(report, null, 2));
