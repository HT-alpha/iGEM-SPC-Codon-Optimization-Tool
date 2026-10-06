/* B-form DNA double helix. Builds the SVG once from geometry, then lets CSS
   (codonforge-ui.css, "DNA helix") move it. There is no requestAnimationFrame
   loop and no JS media-query check: the animations are CSS keyframes, the
   browser does not render them in a hidden tab, and prefers-reduced-motion is
   a CSS rule that removes them.

   The module touches only the container it is given, so the Codon and Primer
   tabs each create their own instance. The page's polling code is the only
   caller of setState().

   Geometry
   - Rise 3.4 A and twist 36 deg per base pair, so 10 bp and 34 A per turn.
     Watson and Crick, Nature 171:737 (1953). PDB 1BNA measures 3.35 A and
     35.8 deg per step (docs/helix-1bna-measurement.json).
   - Phosphate radius 9.3 A and a 145 deg offset between the two phosphate
     tracks (minor groove arc 145, major groove arc 215), both measured on 1BNA
     (Drew et al., PNAS 78:2179, 1981) by tools/measure_bform_1bna.py. Watson
     and Crick give 10 A for the radius; 1BNA's spread is +/-0.5 A.
   - Right-handed. Strand A runs 5' to 3' up the page, strand B 5' to 3' down
     it, so strand B at each rung is the complement of strand A (antiparallel).
   - Scale is one number (SCALE units per angstrom). Nothing else is in pixels.

   The 3D look is a depth cue, not a renderer: opacity and stroke width follow
   z, and z is the sine of the strand's angle about the axis.

   How the motion stays in CSS. Rotating a helix by an angle theta is the same
   curve slid along its axis by theta/360 of a pitch, so the whole backbone
   group is one translateY keyframe. Each rung stays at its height and its two
   halves are unit segments placed by cos/sin of the rung's mid angle b:
   translateX(cx + m cos b) scaleX(h sin b), with m and h set on the SVG below.
   Every half runs the same keyframes and differs only by a negative delay, its
   own phase. All of these animations have one length and start together, and
   setState() changes their playback rates in one task, so the backbone and the
   rungs turn as one helix. Each half also carries its turn-0 pose inline, so
   with the animations removed the picture is still the right geometry.
   Keyframes touch transform and opacity only, which keeps the cost off the
   main thread's style and layout work. */

export const B_FORM = Object.freeze({
  risePerBp: 3.4,
  twistPerBp: 36,
  phosphateRadius: 9.3,
  strandOffset: 145,
});

/* First 20 nt of the lacZ CDS (NP_414878.1) in benchmarks/corpus/
   NC_000913.3_cds.fna.gz, the panel's reporter gene. Two full turns. A test
   compares this string with the shipped file. */
export const DEFAULT_SEQUENCE = "ATGACCATGATTACGGATTC";

export const COMPLEMENT = Object.freeze({ A: "T", T: "A", G: "C", C: "G" });
export const STATES = Object.freeze(["idle", "running", "result", "error"]);

const VIEW = Object.freeze({ width: 100, height: 150 });
const SCALE = 2;
const STEP_DEG = 20;
const START_ANGLE = 208;

/* One turn per 9 s idle and per 5 s running. The stylesheet declares the 9 s
   animation; running multiplies its playback rate, which changes speed without
   moving the phase. */
const IDLE_PERIOD_MS = 9000;
const RUNNING_PERIOD_MS = 5000;

const rad = (deg) => (deg * Math.PI) / 180;
const round2 = (n) => Math.round(n * 100) / 100;
const round1 = (n) => Math.round(n * 10) / 10;
const round3 = (n) => Math.round(n * 1e3) / 1e3;
const round5 = (n) => Math.round(n * 1e5) / 1e5;

export function helixGeometry(options = {}) {
  const sequence = String(options.sequence || DEFAULT_SEQUENCE).toUpperCase();
  if (!/^[ACGT]+$/.test(sequence)) {
    throw new RangeError("sequence must be A, C, G and T only");
  }
  const { risePerBp, twistPerBp, phosphateRadius, strandOffset } = B_FORM;
  const count = sequence.length;
  const rise = risePerBp * SCALE;
  const radius = phosphateRadius * SCALE;
  const bpPerTurn = 360 / twistPerBp;
  const pitch = rise * bpPerTurn;
  const centre = VIEW.width / 2;
  const yBottom = (VIEW.height + (count - 1) * rise) / 2;

  const rungs = [...sequence].map((base, index) => ({
    index,
    y: round1(yBottom - index * rise),
    baseA: base,
    baseB: COMPLEMENT[base],
    angleA: START_ANGLE + index * twistPerBp,
    midAngle: START_ANGLE + index * twistPerBp + strandOffset / 2,
  }));
  const rungMid = radius * Math.cos(rad(strandOffset / 2));
  const rungHalf = radius * Math.sin(rad(strandOffset / 2));

  /* Where a rung half sits at turn 0. b is the rung's mid angle. Strand A ends
     at angle b - offset/2 and strand B at b + offset/2; depth is 1 at the front
     and 0 at the back. A half is a unit line, [0, 1] for A and [0, -1] for B,
     so scaleX(h sin b) puts A's end at m + h sin b and B's at m - h sin b. */
  const halfPose = (midAngle, side) => {
    const b = rad(midAngle);
    const endAngle = b + (side === "a" ? -1 : 1) * rad(strandOffset / 2);
    return {
      x: round3(centre + rungMid * Math.cos(b)),
      scale: round3(rungHalf * Math.sin(b)),
      opacity: round3(0.3 + (0.7 * (1 - Math.sin(endAngle))) / 2),
    };
  };
  for (const rung of rungs) {
    rung.a = halfPose(rung.midAngle, "a");
    rung.b = halfPose(rung.midAngle, "b");
    rung.phase = round5(((((rung.midAngle % 360) + 360) % 360) / 360));
  }

  /* The backbone group slides down by up to one pitch, so it carries one extra
     pitch above the window. At a slide of exactly one pitch the picture is the
     starting picture again, which is what makes the loop seamless. */
  const stepBp = STEP_DEG / twistPerBp;
  const from = -2;
  const to = count - 1 + bpPerTurn + 2;
  const pieces = Math.ceil((to - from) / stepBp);
  const segments = [];
  for (const [strand, offset] of [
    ["a", 0],
    ["b", strandOffset],
  ]) {
    const angle = (j) => rad(START_ANGLE + offset + twistPerBp * j);
    const point = (j) => [centre + radius * Math.cos(angle(j)), yBottom - j * rise];
    const slope = (j) => [-radius * Math.sin(angle(j)) * rad(twistPerBp), -rise];
    for (let k = 0; k < pieces; k += 1) {
      const j1 = from + k * stepBp;
      const j2 = j1 + stepBp;
      const [x1, y1] = point(j1);
      const [x2, y2] = point(j2);
      const [dx1, dy1] = slope(j1);
      const [dx2, dy2] = slope(j2);
      const h = stepBp / 3;
      const front = (1 - Math.sin(angle((j1 + j2) / 2))) / 2;
      const d =
        `M${round1(x1)} ${round1(y1)}C${round1(x1 + dx1 * h)} ${round1(y1 + dy1 * h)} ` +
        `${round1(x2 - dx2 * h)} ${round1(y2 - dy2 * h)} ${round1(x2)} ${round1(y2)}`;
      segments.push({ strand, front: round1(front), d });
    }
  }
  /* Painter's order, back to front. A strand's depth at a given point of the
     pattern never changes as the pattern slides, so this order holds at every
     angle of rotation. */
  segments.sort((p, q) => p.front - q.front);

  return {
    viewBox: [0, 0, VIEW.width, VIEW.height],
    sequence,
    centre,
    yBottom,
    startAngle: START_ANGLE,
    rise,
    radius,
    pitch,
    bpPerTurn,
    strandOffset,
    rungs,
    segments,
    /* Mid-point and half-length of a rung, as amplitudes of cos and sin of the
       rung's mid angle b: m = R cos(offset/2) cos b, h = R sin(offset/2) sin b. */
    rungMid,
    rungHalf,
  };
}

const SVG_NS = "http://www.w3.org/2000/svg";

function svgEl(doc, name, attrs) {
  const node = doc.createElementNS(SVG_NS, name);
  for (const [key, value] of Object.entries(attrs)) node.setAttribute(key, String(value));
  return node;
}

function buildSvg(doc, geometry) {
  const svg = svgEl(doc, "svg", {
    class: "hx-svg",
    viewBox: geometry.viewBox.join(" "),
    preserveAspectRatio: "xMidYMid meet",
    "aria-hidden": "true",
    focusable: "false",
  });
  const vars = {
    "--hx-cx": geometry.centre,
    "--hx-m": round5(geometry.rungMid),
    "--hx-h": round5(geometry.rungHalf),
    "--hx-pitch": geometry.pitch,
    "--hx-delta": `${geometry.strandOffset}deg`,
    "--hx-dphase": round5(geometry.strandOffset / 360),
  };
  for (const [name, value] of Object.entries(vars)) svg.style.setProperty(name, String(value));

  const rungs = svgEl(doc, "g", { class: "hx-rungs" });
  for (const rung of geometry.rungs) {
    const group = svgEl(doc, "g", { class: "hx-rung" });
    group.style.setProperty("--hx-i", String(rung.index));
    group.style.setProperty("--hx-f", String(rung.phase));
    for (const [strand, base, pose, reach] of [
      ["a", rung.baseA, rung.a, 1],
      ["b", rung.baseB, rung.b, -1],
    ]) {
      const half = svgEl(doc, "path", {
        class: `hx-${strand} hx-n${base}`,
        d: `M0 ${rung.y}h${reach}`,
      });
      half.style.transform = `translateX(${pose.x}px) scaleX(${pose.scale})`;
      half.style.opacity = String(pose.opacity);
      group.append(half);
    }
    rungs.append(group);
  }

  const backbone = svgEl(doc, "g", { class: "hx-bb" });
  for (const segment of geometry.segments) {
    backbone.append(
      svgEl(doc, "path", {
        d: segment.d,
        "stroke-width": round1(1.8 + 1.8 * segment.front),
        opacity: round1(0.35 + 0.65 * segment.front),
      }),
    );
  }
  svg.append(rungs, backbone);
  return { svg, rungs, backbone };
}

export function createHelix(container, options = {}) {
  const { svg, backbone } = buildSvg(container.ownerDocument, helixGeometry(options));
  container.replaceChildren(svg);

  let state = "idle";
  let pendingTurn = null;
  container.dataset.helixState = state;
  container.dataset.helixPose = "centre";

  const TURNING = ["hx-move", "hx-fade", "hx-slide"];
  const turning = () =>
    svg.getAnimations({ subtree: true }).filter((a) => TURNING.includes(a.animationName));

  function setRate(rate) {
    for (const animation of turning()) animation.updatePlaybackRate(rate);
  }

  function cancelPending() {
    if (!pendingTurn) return;
    backbone.removeEventListener("animationiteration", pendingTurn.onTurn);
    clearTimeout(pendingTurn.timer);
    pendingTurn = null;
  }

  /* Run `work` when the rotation next completes a whole turn. With no running
     rotation (reduced motion, or a browser without the animation) there is no
     turn to wait for, so it runs now. The timer is a backstop for a browser
     that never delivers the iteration event; it is longer than one running
     turn, so it does not fire first in ordinary use. */
  function atTurnEnd(work) {
    cancelPending();
    if (!turning().length) {
      work();
      return;
    }
    const finish = () => {
      cancelPending();
      work();
    };
    const onTurn = (event) => {
      if (event.target === backbone && event.animationName === "hx-slide") finish();
    };
    backbone.addEventListener("animationiteration", onTurn);
    pendingTurn = { onTurn, timer: setTimeout(finish, RUNNING_PERIOD_MS * 1.5) };
  }

  function setState(next) {
    if (!STATES.includes(next)) throw new RangeError(`unknown helix state: ${next}`);
    if (next === state) return;
    cancelPending();
    state = next;
    container.dataset.helixState = next;
    if (next === "idle") {
      setRate(1);
      container.dataset.helixPose = "centre";
    } else if (next === "running") {
      setRate(IDLE_PERIOD_MS / RUNNING_PERIOD_MS);
      container.dataset.helixPose = "corner";
    } else if (next === "result") {
      atTurnEnd(() => {
        setRate(1);
        container.dataset.helixPose = "badge";
      });
    }
  }

  return { setState };
}
