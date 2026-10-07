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
   through nt 345. The selection stayed on the one nucleotide where it began. */

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
