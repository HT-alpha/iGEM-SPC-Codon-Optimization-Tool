#!/usr/bin/env python3
"""Measure B-form helix parameters from a PDB coordinate file (1BNA).

The DNA helix in the web UI (python/dna_design_suite/api/static/helix.js) takes two
numbers the brief does not give: the radius of the phosphate backbone and the
angular offset between the two strands. This script derives both, plus the
twist and rise the brief does give, from the coordinates of PDB entry 1BNA
(Drew et al., PNAS 78:2179, 1981), so none of them is chosen by eye.

    python3 tools/measure_bform_1bna.py 1bna.pdb            # print JSON
    python3 tools/measure_bform_1bna.py 1bna.pdb --check    # compare to docs/

The coordinate file is not vendored. docs/helix-1bna-measurement.json records
the sha256 of the file the numbers came from. Standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

Vec = tuple[float, float, float]

RECORD = Path(__file__).resolve().parent.parent / "docs" / "helix-1bna-measurement.json"


def read_atoms(text: str) -> dict[tuple[str, int, str], Vec]:
    atoms: dict[tuple[str, int, str], Vec] = {}
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        key = (line[21], int(line[22:26]), line[12:16].strip())
        atoms[key] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return atoms


def sub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a: Vec, b: Vec) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: Vec, b: Vec) -> Vec:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def unit(a: Vec) -> Vec:
    n = math.sqrt(dot(a, a))
    return (a[0] / n, a[1] / n, a[2] / n)


def mean_sd(values: list[float]) -> dict[str, float]:
    n = len(values)
    mean = sum(values) / n
    sd = math.sqrt(sum((v - mean) ** 2 for v in values) / (n - 1))
    return {"mean": round(mean, 3), "sd": round(sd, 3), "n": n}


def principal_axis(points: list[Vec]) -> Vec:
    """Dominant eigenvector of the covariance of `points` by power iteration."""
    n = len(points)
    centre = tuple(sum(p[i] for p in points) / n for i in range(3))
    cov = [[0.0] * 3 for _ in range(3)]
    for p in points:
        d = sub(p, centre)  # type: ignore[arg-type]
        for i in range(3):
            for j in range(3):
                cov[i][j] += d[i] * d[j]
    v: Vec = (0.3, 0.5, 0.8)
    for _ in range(500):
        w = tuple(sum(cov[i][j] * v[j] for j in range(3)) for i in range(3))
        v = unit(w)  # type: ignore[arg-type]
    return v


def nelder_mead(f, x0: list[float], step: list[float], iters: int = 600) -> list[float]:
    """Minimise f by the Nelder-Mead simplex method."""
    n = len(x0)
    simplex = [x0[:]] + [[x0[j] + (step[j] if j == i else 0.0) for j in range(n)] for i in range(n)]
    values = [f(x) for x in simplex]
    for _ in range(iters):
        order = sorted(range(n + 1), key=lambda i: values[i])
        simplex = [simplex[i] for i in order]
        values = [values[i] for i in order]
        centre = [sum(simplex[i][j] for i in range(n)) / n for j in range(n)]
        worst = simplex[-1]
        reflected = [centre[j] + (centre[j] - worst[j]) for j in range(n)]
        f_reflected = f(reflected)
        if f_reflected < values[0]:
            expanded = [centre[j] + 2 * (centre[j] - worst[j]) for j in range(n)]
            f_expanded = f(expanded)
            simplex[-1], values[-1] = (
                (expanded, f_expanded) if f_expanded < f_reflected else (reflected, f_reflected)
            )
        elif f_reflected < values[-2]:
            simplex[-1], values[-1] = reflected, f_reflected
        else:
            contracted = [centre[j] + 0.5 * (worst[j] - centre[j]) for j in range(n)]
            f_contracted = f(contracted)
            if f_contracted < values[-1]:
                simplex[-1], values[-1] = contracted, f_contracted
            else:
                best = simplex[0]
                simplex = [best] + [
                    [best[j] + 0.5 * (point[j] - best[j]) for j in range(n)]
                    for point in simplex[1:]
                ]
                values = [values[0]] + [f(point) for point in simplex[1:]]
    return simplex[0]


def radial_spread(points: list[Vec], origin: Vec, axis: Vec) -> float:
    """Sum of squared deviations of the points' distances from the axis."""
    radii = []
    for p in points:
        v = sub(p, origin)
        along = dot(v, axis)
        w = (v[0] - along * axis[0], v[1] - along * axis[1], v[2] - along * axis[2])
        radii.append(math.sqrt(dot(w, w)))
    mean = sum(radii) / len(radii)
    return sum((r - mean) ** 2 for r in radii)


def solve3(a: list[list[float]], b: list[float]) -> list[float]:
    """Solve a 3x3 linear system by Gauss-Jordan elimination with pivoting."""
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for col in range(3):
        pivot = max(range(col, 3), key=lambda r: abs(m[r][col]))
        m[col], m[pivot] = m[pivot], m[col]
        for r in range(3):
            if r != col:
                f = m[r][col] / m[col][col]
                m[r] = [x - f * y for x, y in zip(m[r], m[col], strict=True)]
    return [m[i][3] / m[i][i] for i in range(3)]


def measure(pdb_text: str, chain_a: str = "A", chain_b: str = "B") -> dict[str, object]:
    atoms = read_atoms(pdb_text)
    res_a = sorted({k[1] for k in atoms if k[0] == chain_a})
    res_b = sorted({k[1] for k in atoms if k[0] == chain_b})
    if len(res_a) != len(res_b):
        raise ValueError("the two chains are not the same length")
    n = len(res_a)
    pairs = [(res_a[k], res_b[n - 1 - k]) for k in range(n)]

    c1a = [atoms[(chain_a, i, "C1'")] for i, _ in pairs]
    c1b = [atoms[(chain_b, j, "C1'")] for _, j in pairs]
    centres = [tuple((a[i] + b[i]) / 2 for i in range(3)) for a, b in zip(c1a, c1b, strict=True)]
    axis = principal_axis(centres)  # type: ignore[arg-type]
    if dot(axis, sub(centres[-1], centres[0])) < 0:  # type: ignore[arg-type]
        axis = (-axis[0], -axis[1], -axis[2])
    origin = tuple(sum(c[i] for c in centres) / n for i in range(3))
    ref: Vec = (1.0, 0.0, 0.0) if abs(axis[0]) < 0.9 else (0.0, 1.0, 0.0)
    ex = unit(cross(axis, ref))
    ey = cross(axis, ex)

    # The principal axis of the base-pair centres leans on the ends of a short
    # duplex (1.2 turns for 1BNA). An axis is the line the backbone winds
    # around at a fixed distance, so refine it to make the phosphate radii and
    # the C1' radii as constant as possible. Both strands share one axis.
    phosphates = [
        atoms[(c, r, "P")]
        for c, rs in ((chain_a, res_a), (chain_b, res_b))
        for r in rs
        if (c, r, "P") in atoms
    ]
    sugars = c1a + c1b

    def frame(x: list[float]) -> tuple[Vec, Vec]:
        o = tuple(origin[i] + x[0] * ex[i] + x[1] * ey[i] for i in range(3))
        d = unit(tuple(axis[i] + x[2] * ex[i] + x[3] * ey[i] for i in range(3)))
        return o, d  # type: ignore[return-value]

    def objective(x: list[float]) -> float:
        o, d = frame(x)
        return radial_spread(phosphates, o, d) + radial_spread(sugars, o, d)

    best = nelder_mead(objective, [0.0, 0.0, 0.0, 0.0], [0.5, 0.5, 0.05, 0.05])
    origin, axis = frame(best)
    ref = (1.0, 0.0, 0.0) if abs(axis[0]) < 0.9 else (0.0, 1.0, 0.0)
    ex = unit(cross(axis, ref))
    ey = cross(axis, ex)

    def cylindrical(p: Vec) -> tuple[float, float, float]:
        d = sub(p, origin)  # type: ignore[arg-type]
        x, y, z = dot(d, ex), dot(d, ey), dot(d, axis)
        return z, math.hypot(x, y), math.degrees(math.atan2(y, x))

    z_centres = [cylindrical(c)[0] for c in centres]  # type: ignore[arg-type]
    rises = [z_centres[k + 1] - z_centres[k] for k in range(n - 1)]

    spans = []
    for a, b in zip(c1a, c1b, strict=True):
        v = sub(b, a)
        v = sub(v, tuple(dot(v, axis) * axis[i] for i in range(3)))  # type: ignore[arg-type]
        spans.append(unit(v))
    twists = []
    for k in range(n - 1):
        twists.append(
            math.degrees(
                math.atan2(dot(cross(spans[k], spans[k + 1]), axis), dot(spans[k], spans[k + 1]))
            )
        )
    mean_twist = sum(twists) / len(twists)
    handedness = "right" if mean_twist > 0 else "left"
    twists = [abs(t) for t in twists]

    radii = []
    tracks: dict[str, list[tuple[float, float]]] = {chain_a: [], chain_b: []}
    for chain, residues in ((chain_a, res_a), (chain_b, res_b)):
        for r in residues:
            p = atoms.get((chain, r, "P"))
            if p is None:
                continue  # the 5' residue of each strand has no phosphate
            z, rad, phi = cylindrical(p)
            radii.append(rad)
            tracks[chain].append((z, phi))

    # Both strands wind at one shared rate. Fit phi = phi0[strand] + k * z with
    # a common k, unwrapping each strand against the C1' winding rate first.
    sign = 1.0 if handedness == "right" else -1.0
    k_guess = sign * (sum(twists) / len(twists)) / (sum(rises) / len(rises))
    rows: list[tuple[int, float, float]] = []
    for idx, chain in enumerate((chain_a, chain_b)):
        pts = sorted(tracks[chain])
        z0, p0 = pts[0]
        for z, phi in pts:
            predicted = p0 + k_guess * (z - z0)
            phi += 360.0 * round((predicted - phi) / 360.0)
            rows.append((idx, z, phi))
    ata = [[0.0] * 3 for _ in range(3)]
    atb = [0.0] * 3
    for idx, z, phi in rows:
        x = [1.0 if idx == 0 else 0.0, 1.0 if idx == 1 else 0.0, z]
        for i in range(3):
            atb[i] += x[i] * phi
            for j in range(3):
                ata[i][j] += x[i] * x[j]
    phi0_a, phi0_b, slope = solve3(ata, atb)
    offset = (phi0_b - phi0_a) % 360.0
    if sign < 0:
        offset = (-offset) % 360.0

    return {
        "pdb_id": "1BNA",
        "handedness": handedness,
        "base_pairs": n,
        "phosphate_radius_angstrom": mean_sd(radii),
        "twist_deg_per_bp": mean_sd(twists),
        "rise_angstrom_per_bp": mean_sd(rises),
        "bp_per_turn": round(360.0 / (sum(twists) / len(twists)), 2),
        "phosphate_track_offset_deg": round(offset, 1),
        "minor_groove_arc_deg": round(min(offset, 360.0 - offset), 1),
        "major_groove_arc_deg": round(max(offset, 360.0 - offset), 1),
        "phosphate_winding_deg_per_angstrom": round(abs(slope), 3),
    }


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("pdb", type=Path, help="path to 1bna.pdb")
    parser.add_argument("--check", action="store_true", help="compare to the committed record")
    args = parser.parse_args(argv)

    raw = args.pdb.read_bytes()
    result = measure(raw.decode("ascii"))
    result["source_sha256"] = hashlib.sha256(raw).hexdigest()
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if not args.check:
        sys.stdout.write(text)
        return 0
    if json.loads(RECORD.read_text(encoding="utf-8")) != result:
        sys.stderr.write(f"measurement differs from {RECORD}\n{text}")
        return 1
    print(f"matches {RECORD.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
