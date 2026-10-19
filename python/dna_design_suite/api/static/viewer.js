/* Which nucleotide a pointer event is on.

   pointerdown calls setPointerCapture on the viewer, so a drag keeps
   delivering pointermove after the pointer leaves the nucleotide it started
   on. Chrome then retargets those events to the viewer: event.target is no
   longer the nucleotide, and closest("[data-index]") on it is null. The
   nucleotide is the element at the pointer's coordinates, which the caller
   gets from document.elementFromPoint and passes as `under`.

   Measured on Chrome 148 before this split: a drag from index 300 to 345
   delivered one move whose target was the nucleotide and then moves whose
   target was the viewer, while the element at the coordinates was nt 304
   through nt 345. The selection stayed on the one nucleotide where it began.

   A wrapped line has gaps: the gutter between groups of ten, the margin
   between lines, and the position numbers beside a line. Those are not
   nucleotides, so elementFromPoint misses. nearestLineIndex and nearestBase
   map that pointer onto the line and base it is nearest, which is how a drag
   continues from one line to the next. */

function closestIndex(node) {
  if (!node || typeof node.closest !== "function") return null;
  const hit = node.closest("[data-index]");
  if (!hit || hit.dataset.index === undefined || hit.dataset.index === "") return null;
  const index = Number(hit.dataset.index);
  return Number.isInteger(index) ? index : null;
}

export function nucleotideIndex(target, under) {
  const direct = closestIndex(target);
  if (direct !== null) return direct;
  return closestIndex(under);
}

export const BASES_PER_LINE = 60;
export const BASES_PER_GROUP = 10;

export function sequenceLines(length, basesPerLine = BASES_PER_LINE) {
  const total = Number(length);
  if (!Number.isInteger(total) || total <= 0) return [];
  const size = Number.isInteger(basesPerLine) && basesPerLine > 0 ? basesPerLine : BASES_PER_LINE;
  const lines = [];
  for (let start = 0; start < total; start += size) {
    lines.push({ start, end: Math.min(total, start + size) });
  }
  return lines;
}

/* selection is half-open, the same {start, end} the primer request posts.
   The label is 1-based and inclusive, which is what the readout and the
   start and end boxes show. */
export function selectionLabel(selection) {
  if (!selection) return null;
  const start = selection.start;
  const end = selection.end;
  if (!Number.isInteger(start) || !Number.isInteger(end) || start < 0 || end <= start) return null;
  return { startBase: start + 1, endBase: end, length: end - start };
}

export function nearestLineIndex(pointerY, lines, slack = 8) {
  if (!Array.isArray(lines) || lines.length === 0 || !Number.isFinite(pointerY)) return null;
  let best = null;
  let bestDy = Infinity;
  for (let index = 0; index < lines.length; index += 1) {
    const line = lines[index];
    if (!line) continue;
    const top = Number(line.top);
    const bottom = Number(line.bottom);
    if (!Number.isFinite(top) || !Number.isFinite(bottom)) continue;
    const dy = pointerY < top ? top - pointerY : pointerY > bottom ? pointerY - bottom : 0;
    if (dy < bestDy) {
      bestDy = dy;
      best = index;
    }
  }
  if (best === null || bestDy > slack) return null;
  return best;
}

export function nearestBase(pointerX, bases) {
  if (!Array.isArray(bases) || !Number.isFinite(pointerX)) return null;
  let chosen = null;
  let bestScore = Infinity;
  for (const base of bases) {
    if (!base || !Number.isInteger(base.index)) continue;
    const left = Number(base.left);
    const right = Number(base.right);
    if (!Number.isFinite(left) || !Number.isFinite(right)) continue;
    const dx = pointerX < left ? left - pointerX : pointerX > right ? pointerX - right : 0;
    const center = (left + right) / 2;
    const score = dx === 0 ? Math.abs(pointerX - center) * 0.001 : dx;
    if (score < bestScore) {
      bestScore = score;
      chosen = base.index;
    }
  }
  return chosen;
}

function escapeBase(base) {
  if (base === "&") return "&amp;";
  if (base === "<") return "&lt;";
  if (base === ">") return "&gt;";
  return base;
}

function baseClass(base) {
  if (base === "A" || base === "C" || base === "G" || base === "T") return `nt nt-${base.toLowerCase()}`;
  return "nt";
}

function renderScale() {
  const ticks = [];
  for (let tick = 1; tick <= BASES_PER_LINE; tick += BASES_PER_GROUP) {
    ticks.push(`<span class="scale-tick">${tick}</span>`);
    if (tick + BASES_PER_GROUP <= BASES_PER_LINE) ticks.push('<span class="nt-gap"></span>');
  }
  return `<div class="seq-ruler-head" aria-hidden="true"><span class="seq-ruler"></span><span class="seq-bases">${ticks.join("")}</span><span class="seq-ruler"></span></div>`;
}

/* origin is the absolute index of sequence[0], so a paged window keeps the
   same base numbers as the whole template. */
export function renderSequenceLines(sequence, origin = 0) {
  const text = String(sequence ?? "");
  const startAt = Number.isInteger(origin) && origin >= 0 ? origin : 0;
  const lines = sequenceLines(text.length);
  if (!lines.length) return "";
  const rows = lines.map((line) => {
    const parts = [];
    for (let index = line.start; index < line.end; index += BASES_PER_GROUP) {
      const groupEnd = Math.min(line.end, index + BASES_PER_GROUP);
      let chars = "";
      for (let cursor = index; cursor < groupEnd; cursor += 1) {
        const base = text[cursor];
        chars += `<span class="${baseClass(base)}" data-index="${startAt + cursor}">${escapeBase(base)}</span>`;
      }
      parts.push(`<span class="nt-group">${chars}</span>`);
      if (groupEnd < line.end) parts.push('<span class="nt-gap"></span>');
    }
    const absoluteStart = startAt + line.start;
    const absoluteEnd = startAt + line.end;
    return `<div class="seq-line" data-line-start="${absoluteStart}" data-line-end="${absoluteEnd}"><span class="seq-ruler seq-ruler-start"><span class="line-pos">${absoluteStart + 1}</span><span class="roi-mark roi-mark-start"></span></span><span class="seq-bases">${parts.join("")}</span><span class="seq-ruler seq-ruler-end"><span class="line-pos">${absoluteEnd}</span><span class="roi-mark roi-mark-end"></span></span></div>`;
  });
  return renderScale() + rows.join("");
}
