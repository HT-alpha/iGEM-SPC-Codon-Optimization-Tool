import { readFileSync, writeFileSync } from "node:fs";

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
for (let i = 0; i < 20; i += 1) {
  const ready = await evaluate("Boolean(document.getElementById('tab-btn-primer'))");
  if (ready) break;
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
  document.getElementById("primer-flank").value = "80";
  document.getElementById("primer-pmin").value = "150";
  document.getElementById("primer-pmax").value = "400";
  document.getElementById("region-start").value = "201";
  document.getElementById("region-end").value = "400";
  document.getElementById("region-start").dispatchEvent(new Event("change"));
  const cleaned = sequence.toUpperCase().replace(/[^ACGT]/g, "");
  document.getElementById("primer-search").click();
  cleaned
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

ws.close();
writeFileSync(OUT, JSON.stringify(report, null, 2));
console.log(JSON.stringify(report, null, 2));
