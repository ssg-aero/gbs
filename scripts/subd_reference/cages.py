"""Cages of reference of the Catmull-Clark subdivision (gbs/subd), written to tests/in/subd/:
each configuration of the rules at least once, with coordinates perturbed so that no symmetry
hides a swapped parameter.

- box: a closed hexahedron, eight extraordinary vertices of valence 3;
- torus: closed, regular everywhere;
- ngons: a pentagonal prism capped by five triangles, extraordinary vertices of valence 5 and
  faces of 3 and 5 sides;
- grid: an open grid of quadrangles, regular boundaries and corners;
- fans: two fans of quadrangles around a boundary vertex, of 3 and 4 faces (irregular creases);
- creases: a hexahedron split once, an equator of creases (regular crease vertices), a crease
  from a corner (a dart at its end, a corner where it meets the equator), two creases at a
  right angle on the top (a crease vertex of 1 and 3 faces on its sides), a tagged corner;
- hexprism: a hexagonal prism, whose hexagons become extraordinary vertices of valence 6 when
  refined.

    python cages.py [output directory]
"""
from __future__ import annotations

import json
import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]

rng = np.random.default_rng(20261103)


def perturbed(vertices, amount):
    v = np.asarray(vertices, float)
    return (v + amount * rng.uniform(-1.0, 1.0, v.shape)).tolist()


def cage(vertices, faces, creases=(), corners=()):
    return {"vertices": [list(map(float, p)) for p in vertices], "faces": [list(map(int, f)) for f in faces],
            "creases": [list(map(int, c)) for c in creases], "corners": list(map(int, corners))}


def box():
    v = [[-1, -1, 1], [1, -1, 1], [-1, 1, 1], [1, 1, 1], [-1, 1, -1], [1, 1, -1], [-1, -1, -1], [1, -1, -1]]
    f = [[0, 1, 3, 2], [2, 3, 5, 4], [4, 5, 7, 6], [6, 7, 1, 0], [1, 7, 5, 3], [6, 0, 2, 4]]
    v = np.array(v, float) * [1.3, 0.9, 0.7]
    return cage(perturbed(v, 0.15), f)


def torus(n=8, m=6, R=2.0, r=0.7):
    v = [[(R + r * math.cos(2 * math.pi * j / m)) * math.cos(2 * math.pi * i / n),
          (R + r * math.cos(2 * math.pi * j / m)) * math.sin(2 * math.pi * i / n),
          r * math.sin(2 * math.pi * j / m)] for i in range(n) for j in range(m)]
    idx = lambda i, j: (i % n) * m + j % m
    f = [[idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)] for i in range(n) for j in range(m)]
    return cage(perturbed(v, 0.1), f)


def ngons():
    b = [[math.cos(2 * math.pi * i / 5), math.sin(2 * math.pi * i / 5), 0.0] for i in range(5)]
    t = [[p[0], p[1], 1.2] for p in b]
    v = b + t + [[0.0, 0.0, 1.9]]
    f = [[4, 3, 2, 1, 0]] + [[i, (i + 1) % 5, 5 + (i + 1) % 5, 5 + i] for i in range(5)]
    f += [[5 + i, 5 + (i + 1) % 5, 10] for i in range(5)]
    return cage(perturbed(v, 0.1), f)


def grid(nx=4, ny=3):
    v = [[i, j, 0.3 * math.sin(1.3 * i) * math.cos(0.9 * j)] for j in range(ny + 1) for i in range(nx + 1)]
    idx = lambda i, j: j * (nx + 1) + i
    f = [[idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)] for j in range(ny) for i in range(nx)]
    return cage(perturbed(v, 0.1), f)


def fans():
    vertices, faces = [], []
    for k, x0 in ((3, 0.0), (4, 4.0)):
        c = len(vertices)
        vertices.append([x0, 0.0, 0.0])
        for i in range(k + 1):
            a = math.pi * i / k
            vertices.append([x0 + math.cos(a), math.sin(a), 0.4 * math.sin(2 * a)])
        for i in range(k):
            a = math.pi * (i + 0.5) / k
            vertices.append([x0 + 1.4 * math.cos(a), 1.4 * math.sin(a), 0.3 * math.cos(a)])
        e = lambda i: c + 1 + i
        d = lambda i: c + 2 + k + i
        faces += [[c, e(i), d(i), e(i + 1)] for i in range(k)]
    return cage(perturbed(vertices, 0.05), faces)


def split(vertices, faces):
    """Each face of n sides into n quadrangles by its centroid and the middles of its sides (the
    topology of one subdivision step, linear positions)."""
    v = [np.array(p, float) for p in vertices]
    edges, out = {}, []
    centre = []
    for f in faces:
        centre.append(len(v))
        v.append(np.mean([v[i] for i in f], axis=0))

    def mid(a, b):
        key = (min(a, b), max(a, b))
        if key not in edges:
            edges[key] = len(v)
            v.append((v[a] + v[b]) / 2)
        return edges[key]

    for fi, f in enumerate(faces):
        n = len(f)
        for i in range(n):
            out.append([f[i], mid(f[i], f[(i + 1) % n]), centre[fi], mid(f[i - 1], f[i])])
    return [p.tolist() for p in v], out


def creases():
    b = box()
    v0 = [[-1, -1, 1], [1, -1, 1], [-1, 1, 1], [1, 1, 1], [-1, 1, -1], [1, 1, -1], [-1, -1, -1], [1, -1, -1]]
    v, f = split(v0, b["faces"])
    v = np.array(v)
    edges = {tuple(sorted((q[i], q[(i + 1) % 4]))) for q in f for i in range(4)}
    find = lambda p: int(np.argmin(np.linalg.norm(v - np.array(p, float), axis=1)))
    equator = [e for e in edges if abs(v[e[0]][2]) < 1e-9 and abs(v[e[1]][2]) < 1e-9]
    down = [tuple(sorted((find([1, 1, 1]), find([1, 1, 0]))))]                       # dart at (1, 1, 1)
    top = [tuple(sorted((find([0, 0, 1]), find([1, 0, 1])))), tuple(sorted((find([0, 0, 1]), find([0, 1, 1]))))]
    corners = [find([-1, -1, -1])]
    v = v * [1.3, 0.9, 0.7]
    return cage(perturbed(v, 0.08), f, equator + down + top, corners)


def hexprism():
    b = [[1.4 * math.cos(math.pi * i / 3), math.sin(math.pi * i / 3), 0.0] for i in range(6)]
    v = b + [[p[0], p[1], 1.5] for p in b]
    f = [[5, 4, 3, 2, 1, 0], [6, 7, 8, 9, 10, 11]] + [[i, (i + 1) % 6, 6 + (i + 1) % 6, 6 + i] for i in range(6)]
    return cage(perturbed(v, 0.1), f)


CAGES = {"box": box, "torus": torus, "ngons": ngons, "grid": grid, "fans": fans, "creases": creases, "hexprism": hexprism}

if __name__ == "__main__":
    out = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "tests" / "in" / "subd"
    out.mkdir(parents=True, exist_ok=True)
    for name, make in CAGES.items():
        with open(out / f"{name}.cage.json", "w") as f:
            json.dump(make(), f)
        print("written", out / f"{name}.cage.json")
