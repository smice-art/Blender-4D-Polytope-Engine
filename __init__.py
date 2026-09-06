# SPDX-License-Identifier: GPL-3.0-or-later
#
# 4D Polytopes -- Blender add-on. Inpired and derivative of the "Math Art"
# project (https://github.com/dwkrider/Math-Art) by David Krider,
# licensed GPL-3.0-or-later; redistribution must keep this license,
# make source available, and credit the upstream project. This note
# is informational, not legal advice.

import math
import time
import itertools
import colorsys
from math import cos, sin, pi, sqrt

import bpy
from bpy.props import (FloatProperty, EnumProperty, IntProperty,
                       BoolProperty, PointerProperty)
from bpy.app.handlers import persistent
import bmesh

try:
    import numpy as np
except ImportError:
    np = None


PHI = (1.0 + 5.0 ** 0.5) / 2.0


def _seed_unit(v):
    l = math.sqrt(sum(x * x for x in v)) or 1.0
    return tuple(x / l for x in v)


def icosa_faces(V):
    """Recover the 20 triangular faces of an icosahedron from its 12
    vertices alone, by walking common neighbours at the minimal edge
    length."""
    n = len(V)
    emin = None
    d2 = {}
    for i in range(n):
        for j in range(i + 1, n):
            d = sum((V[i][k] - V[j][k]) ** 2 for k in range(3))
            d2[(i, j)] = d
            emin = d if emin is None else min(emin, d)
    adj = {i: set() for i in range(n)}
    for (i, j), d in d2.items():
        if d < emin * 1.2:
            adj[i].add(j)
            adj[j].add(i)
    fs = set()
    for i in range(n):
        for j in adj[i]:
            for k in adj[i] & adj[j]:
                fs.add(tuple(sorted((i, j, k))))
    faces = []
    for f in fs:
        a, b, c = (V[i] for i in f)
        nx = ((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]),
              (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2]),
              (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
        cen = [(a[k] + b[k] + c[k]) / 3 for k in range(3)]
        if sum(nx[k] * cen[k] for k in range(3)) < 0:
            faces.append([f[0], f[2], f[1]])
        else:
            faces.append(list(f))
    return faces


def seed_poly(kind, unit=False):
    """(vertices, faces) of one of the five Platonic solids."""
    if kind == 'TETRA':
        V = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
        F = [(0, 1, 2), (0, 2, 3), (0, 3, 1), (1, 3, 2)]
    elif kind == 'OCTA':
        V = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0),
             (0, 0, 1), (0, 0, -1)]
        F = [(0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4),
             (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5)]
    elif kind == 'CUBE':
        V = [(x, y, z) for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
        F = [(0, 1, 3, 2), (4, 6, 7, 5), (0, 4, 5, 1),
             (2, 3, 7, 6), (0, 2, 6, 4), (1, 5, 7, 3)]
    elif kind == 'ICOSA':
        V = []
        for a in (-1, 1):
            for b in (-PHI, PHI):
                V += [(0, a, b), (a, b, 0), (b, 0, a)]
        F = icosa_faces(V)
    elif kind == 'DODECA':
        IV = []
        for a in (-1, 1):
            for b in (-PHI, PHI):
                IV += [(0, a, b), (a, b, 0), (b, 0, a)]
        IF = icosa_faces(IV)
        V = [tuple(sum(IV[i][k] for i in f) / 3 for k in range(3))
             for f in IF]
        F = []
        for vi in range(len(IV)):
            adj = [fi for fi, f in enumerate(IF) if vi in f]
            nrm = _seed_unit(IV[vi])
            u0 = V[adj[0]]
            u = [u0[k] - sum(u0[j] * nrm[j] for j in range(3)) * nrm[k]
                 for k in range(3)]
            w = (nrm[1] * u[2] - nrm[2] * u[1],
                 nrm[2] * u[0] - nrm[0] * u[2],
                 nrm[0] * u[1] - nrm[1] * u[0])

            def ang(fi):
                d = V[fi]
                return math.atan2(sum(d[k] * w[k] for k in range(3)),
                                  sum(d[k] * u[k] for k in range(3)))
            F.append(sorted(adj, key=ang))
    else:
        raise ValueError(kind)
    V = [list(v) for v in V]
    if unit:
        r = math.sqrt(sum(c * c for c in V[0]))
        if r > 1e-12:
            V = [[c / r for c in v] for v in V]
    return V, [list(f) for f in F]


def _perms_even(coords):
    out = set()
    for p in itertools.permutations(range(4)):
        parity = 0
        q = list(p)
        for i in range(4):
            while q[i] != i:
                j = q[i]
                q[i], q[j] = q[j], q[i]
                parity += 1
        if parity % 2 == 0:
            out.add(tuple(coords[i] for i in p))
    return out


def _perms_all(coords):
    return set(itertools.permutations(coords))


def _signs(base):
    out = set()
    for v in base:
        idx = [i for i in range(4) if abs(v[i]) > 1e-12]
        for signs in itertools.product((1, -1), repeat=len(idx)):
            w = list(v)
            for k, i in enumerate(idx):
                w[i] = abs(w[i]) * signs[k]
            out.add(tuple(w))
    return out


def polytope_vertices(kind):
    """Standard coordinate set for one of the six regular convex
    4-polytopes, normalised onto the unit 3-sphere."""
    if kind == 'CELL5':
        s5 = sqrt(5)
        V = [(1, 1, 1, -1 / s5), (1, -1, -1, -1 / s5),
             (-1, 1, -1, -1 / s5), (-1, -1, 1, -1 / s5),
             (0, 0, 0, s5 - 1 / s5)]
    elif kind == 'CELL8':
        V = list(itertools.product((-1, 1), repeat=4))
    elif kind == 'CELL16':
        V = list(_signs(_perms_all((1, 0, 0, 0))))
    elif kind == 'CELL24':
        V = list(_signs(_perms_all((1, 1, 0, 0))))
    elif kind == 'CELL600':
        V = set()
        V |= _signs({(0.5, 0.5, 0.5, 0.5)})
        V |= _signs(_perms_all((1, 0, 0, 0)))
        V |= _signs(_perms_even((PHI / 2, 0.5, 1 / (2 * PHI), 0)))
        V = list(V)
    elif kind == 'CELL120':
        s5 = sqrt(5)
        V = set()
        V |= _signs(_perms_all((2, 2, 0, 0)))
        V |= _signs(_perms_all((s5, 1, 1, 1)))
        V |= _signs(_perms_all((PHI, PHI, PHI, PHI ** -2)))
        V |= _signs(_perms_all((PHI ** 2, PHI ** -1, PHI ** -1, PHI ** -1)))
        V |= _signs(_perms_even((PHI ** 2, PHI ** -2, 1, 0)))
        V |= _signs(_perms_even((s5, PHI ** -1, PHI, 0)))
        V |= _signs(_perms_even((2, 1, PHI, PHI ** -1)))
        V = list(V)
    else:
        raise ValueError(kind)
    out = []
    for v in V:
        n = sqrt(sum(x * x for x in v))
        out.append(tuple(x / n for x in v))
    return out


def _in_flat(u, v, w, x, tol=1e-6):
    """Is x in the 2-flat of R^4 through u, v, w?"""
    a = [v[k] - u[k] for k in range(4)]
    b = [w[k] - u[k] for k in range(4)]
    c = [x[k] - u[k] for k in range(4)]
    la = math.sqrt(sum(t * t for t in a)) or 1.0
    a = [t / la for t in a]
    d = sum(b[k] * a[k] for k in range(4))
    b = [b[k] - d * a[k] for k in range(4)]
    lb = math.sqrt(sum(t * t for t in b)) or 1.0
    b = [t / lb for t in b]
    for e in (a, b):
        d = sum(c[k] * e[k] for k in range(4))
        c = [c[k] - d * e[k] for k in range(4)]
    return math.sqrt(sum(t * t for t in c)) < tol


_FACE_CACHE = {}

_FACE_SIZE = {'CELL5': 3, 'CELL8': 4, 'CELL16': 3, 'CELL24': 3,
              'CELL120': 5, 'CELL600': 3, 'DUOPRISM': 4}


def polytope_faces(kind, V, E, cache_key=None):
    """The 2D faces (polygon vertex cycles): walks along edges staying
    inside a common 2-flat."""
    ckey = cache_key or kind
    if ckey in _FACE_CACHE:
        return _FACE_CACHE[ckey]
    want = _FACE_SIZE[kind]
    from collections import defaultdict
    adj = defaultdict(set)
    for i, j in E:
        adj[i].add(j)
        adj[j].add(i)
    seen = set()
    faces = []
    for u in sorted(adj):
        for v in adj[u]:
            for w in adj[v]:
                if w == u:
                    continue
                cyc = [u, v, w]
                closed = False
                while len(cyc) <= want:
                    prev, cur = cyc[-2], cyc[-1]
                    nxt = None
                    for x in adj[cur]:
                        if x != prev and _in_flat(V[u], V[v], V[w], V[x]):
                            nxt = x
                            break
                    if nxt is None:
                        break
                    if nxt == u:
                        closed = True
                        break
                    cyc.append(nxt)
                if closed and len(cyc) == want:
                    key = frozenset(cyc)
                    if key not in seen:
                        seen.add(key)
                        faces.append(cyc)
    _FACE_CACHE[ckey] = faces
    return faces


def polytope_edges(V):
    """Edges = vertex pairs at the minimal nonzero distance."""
    n = len(V)
    d2min = None
    for i in range(n):
        for j in range(i + 1, n):
            d2 = sum((V[i][k] - V[j][k]) ** 2 for k in range(4))
            if d2 > 1e-9 and (d2min is None or d2 < d2min):
                d2min = d2
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            d2 = sum((V[i][k] - V[j][k]) ** 2 for k in range(4))
            if abs(d2 - d2min) < 1e-6:
                edges.append((i, j))
    return edges


COUNTS = {'CELL5': (5, 10), 'CELL8': (16, 32), 'CELL16': (8, 24),
          'CELL24': (24, 96), 'CELL600': (120, 720),
          'CELL120': (600, 1200)}

DUAL_KIND = {'CELL5': 'CELL5', 'CELL8': 'CELL16', 'CELL16': 'CELL8',
             'CELL24': 'CELL24', 'CELL120': 'CELL600',
             'CELL600': 'CELL120'}

_CELL_SIZE = {'CELL5': 4, 'CELL8': 8, 'CELL16': 4, 'CELL24': 6,
              'CELL120': 20, 'CELL600': 4}

HALF_TOL = 1e-6


def dual_vertices(kind):
    """Unit-sphere vertices of the dual polytope. Returns (verts,
    dual_kind, key)."""
    if kind == 'CELL5':
        V = [tuple(-x for x in v) for v in polytope_vertices('CELL5')]
        return V, 'CELL5', 'CELL5_DUAL'
    if kind == 'CELL24':
        S = set()
        S |= _signs(_perms_all((1, 0, 0, 0)))
        S |= _signs({(0.5, 0.5, 0.5, 0.5)})
        V = []
        for v in S:
            n = sqrt(sum(x * x for x in v))
            V.append(tuple(x / n for x in v))
        return V, 'CELL24', 'CELL24_DUAL'
    dk = DUAL_KIND[kind]
    return polytope_vertices(dk), dk, dk


def _cell_inradius(kind):
    V = polytope_vertices(kind)
    d = dual_vertices(kind)[0][0]
    k = _CELL_SIZE[kind]
    idx = sorted(range(len(V)),
                 key=lambda i: -sum(V[i][t] * d[t] for t in range(4)))[:k]
    c = [sum(V[i][t] for i in idx) / k for t in range(4)]
    return sqrt(sum(x * x for x in c))


def _half_filter(V, E, F2, tol=HALF_TOL):
    keep = [i for i, v in enumerate(V) if v[3] <= tol]
    remap = {i: n for n, i in enumerate(keep)}
    V2 = [V[i] for i in keep]
    E2 = [(remap[i], remap[j]) for (i, j) in E if i in remap and j in remap]
    F22 = None
    if F2 is not None:
        F22 = [[remap[i] for i in cyc] for cyc in F2
               if all(i in remap for i in cyc)]
    return V2, E2, F22


def _qmul(a, b):
    return (a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
            a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
            a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
            a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0])


def _qorder(q, maxn=12):
    p = q
    for k in range(1, maxn + 1):
        if (abs(p[0] - 1.0) < 1e-9 and abs(p[1]) < 1e-9
                and abs(p[2]) < 1e-9 and abs(p[3]) < 1e-9):
            return k
        p = _qmul(p, q)
    return 0


def _order10_generator(centers):
    g = (PHI / 2, 0.5, 1 / (2 * PHI), 0.0)
    if _qorder(g) == 10:
        return g
    for c in sorted(centers):
        if _qorder(c) == 10:
            return c
    raise RuntimeError("no order-10 element found in 2I")


def hopf_ring_cosets():
    """The 120 cell centers of the 120-cell partitioned into 12 Hopf
    rings of 10 cells. Returns (centers, rings)."""
    centers = polytope_vertices('CELL600')
    g = _order10_generator(centers)
    gpow = [(1.0, 0.0, 0.0, 0.0)]
    for _ in range(9):
        gpow.append(_qmul(gpow[-1], g))

    def nearest(q):
        bi, bd = 0, -2.0
        for i, c in enumerate(centers):
            d = sum(q[t] * c[t] for t in range(4))
            if d > bd:
                bd, bi = d, i
        return bi, bd

    assigned = [False] * len(centers)
    rings = []
    seeds = sorted(range(len(centers)),
                   key=lambda i: tuple(-x for x in centers[i]))
    for s in seeds:
        if assigned[s]:
            continue
        ring = []
        for gp in gpow:
            i, d = nearest(_qmul(centers[s], gp))
            if d < 1.0 - 1e-9 or assigned[i]:
                raise RuntimeError("coset walk left the group")
            assigned[i] = True
            ring.append(i)
        rings.append(ring)
    ring0 = rings[0]

    def mind2(ring):
        return min(sum((centers[i][t] - centers[j][t]) ** 2
                       for t in range(4))
                   for i in ring for j in ring0)

    rest = sorted(rings[1:], key=lambda r: (-mind2(r), tuple(sorted(r))))
    return centers, [ring0] + rest


def ring_cell_points(n_rings, cell_scale, half=False):
    V = polytope_vertices('CELL120')
    centers, rings = hopf_ring_cosets()
    out = []
    for ri in range(min(n_rings, len(rings))):
        for ci in rings[ri]:
            c = centers[ci]
            if half and c[3] > HALF_TOL:
                continue
            idx = sorted(range(len(V)),
                         key=lambda i: -sum(V[i][t] * c[t]
                                            for t in range(4)))[:20]
            cen = [sum(V[i][t] for i in idx) / 20 for t in range(4)]
            out.append((ri, [tuple(cen[t] + (V[i][t] - cen[t]) * cell_scale
                                   for t in range(4)) for i in idx]))
    return out


def rotate4(V, xw, yw, zw, xy):
    """Rotations in the XW, YW, ZW and XY planes (degrees)."""
    out = [list(v) for v in V]
    for (a, b, ang) in ((0, 3, xw), (1, 3, yw), (2, 3, zw), (0, 1, xy)):
        t = math.radians(ang)
        if abs(t) < 1e-12:
            continue
        c, s = cos(t), sin(t)
        for v in out:
            va, vb = v[a], v[b]
            v[a] = va * c - vb * s
            v[b] = va * s + vb * c
    return [tuple(v) for v in out]


def _slerp4(a, b, t):
    d = max(-1.0, min(1.0, sum(a[k] * b[k] for k in range(4))))
    om = math.acos(d)
    if om < 1e-9:
        return a
    sa = sin((1 - t) * om) / sin(om)
    sb = sin(t * om) / sin(om)
    return tuple(a[k] * sa + b[k] * sb for k in range(4))


def project_point(v, dist):
    denom = max(dist - v[3], 0.02)
    s = 1.0 / denom
    return (v[0] * s, v[1] * s, v[2] * s), s


def _pole_angles(V, target=0.3):
    """A small extra 4D rotation (in the XW and YW planes) that moves
    every vertex away from the stereographic pole (w ~ 1), or (0, 0)
    if none is needed. Uses a dense numpy grid search over both
    angles rather than a closed-form guess: a full, symmetric
    polytope only ever has ONE vertex near the pole, so almost any
    escape direction works -- but a reduced, asymmetric point set
    (e.g. after Half Cutaway followed by a manual rotation) can have
    several vertices to dodge at once, which a single guessed
    direction is not guaranteed to clear."""
    def max_w(vs):
        return max(v[3] for v in vs)
    raw = max_w(V)
    if raw < 0.8:
        return 0.0, 0.0
    if np is None:
        base = math.degrees(math.acos(max(min(target, 0.999), -0.999)))
        for ang in ((base, base * 0.35), (base, 0.0), (0.0, 0.0)):
            if max_w(rotate4(V, ang[0], ang[1], 0.0, 0.0)) < raw:
                return ang
        return (0.0, 0.0)
    A = np.asarray(V, dtype=float)
    x, y, w = A[:, 0], A[:, 1], A[:, 3]
    angles = np.deg2rad(np.arange(0.0, 180.0, 2.0))
    ca, sa = np.cos(angles), np.sin(angles)
    X1 = np.outer(ca, x) - np.outer(sa, w)
    W1 = np.outer(sa, x) + np.outer(ca, w)
    best_w, best = raw, (0.0, 0.0)
    for j in range(len(angles)):
        W2 = sa[j] * y[None, :] + ca[j] * W1
        per_i = W2.max(axis=1)
        i = int(np.argmin(per_i))
        if per_i[i] < best_w:
            best_w = float(per_i[i])
            best = (math.degrees(angles[i]), math.degrees(angles[j]))
        if best_w <= target:
            break
    return best


def _sub3(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross3(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _unit3(v):
    l = math.sqrt(sum(x * x for x in v)) or 1.0
    return (v[0] / l, v[1] / l, v[2] / l)


def add_strut(verts, faces, pts, radii, sides):
    n = len(pts)
    rings = []
    prev_n = None
    for i in range(n):
        if i == 0:
            t = _unit3(_sub3(pts[1], pts[0]))
        elif i == n - 1:
            t = _unit3(_sub3(pts[-1], pts[-2]))
        else:
            t = _unit3(_sub3(pts[i + 1], pts[i - 1]))
        if prev_n is None:
            ref = (0.0, 0.0, 1.0) if abs(t[2]) < 0.9 else (1.0, 0.0, 0.0)
            u = _unit3(_cross3(t, ref))
        else:
            u = _unit3(_cross3(t, _cross3(prev_n, t)))
            d = sum(prev_n[k] * t[k] for k in range(3))
            u = _unit3(tuple(prev_n[k] - d * t[k] for k in range(3)))
        w = _cross3(t, u)
        prev_n = u
        ring = []
        for s in range(sides):
            a = 2 * pi * s / sides
            ring.append(len(verts))
            verts.append(tuple(pts[i][k]
                               + radii[i] * (cos(a) * u[k] + sin(a) * w[k])
                               for k in range(3)))
        rings.append(ring)
    for i in range(n - 1):
        r0, r1 = rings[i], rings[i + 1]
        for s in range(sides):
            s2 = (s + 1) % sides
            faces.append([r0[s], r0[s2], r1[s2], r1[s]])
    faces.append(list(reversed(rings[0])))
    faces.append(list(rings[-1]))


def add_sphere(verts, faces, center, radius, seg=8, rings=6):
    base = len(verts)
    verts.append((center[0], center[1], center[2] + radius))
    for r in range(1, rings):
        th = pi * r / rings
        for s in range(seg):
            a = 2 * pi * s / seg
            verts.append((center[0] + radius * sin(th) * cos(a),
                          center[1] + radius * sin(th) * sin(a),
                          center[2] + radius * cos(th)))
    verts.append((center[0], center[1], center[2] - radius))
    last = len(verts) - 1
    ring0 = lambda r: base + 1 + (r - 1) * seg
    for s in range(seg):
        s2 = (s + 1) % seg
        faces.append([base, ring0(1) + s2, ring0(1) + s])
    for r in range(1, rings - 1):
        for s in range(seg):
            s2 = (s + 1) % seg
            faces.append([ring0(r) + s, ring0(r) + s2,
                          ring0(r + 1) + s2, ring0(r + 1) + s])
    for s in range(seg):
        s2 = (s + 1) % seg
        faces.append([last, ring0(rings - 1) + s, ring0(rings - 1) + s2])


def add_cube(verts, faces, center, size):
    """A simple 8-vertex, 6-face cube vertex marker -- far cheaper
    than add_sphere's several dozen faces (a default sphere is 8
    poles + 5*8 ring vertices, ~48 faces; this is 6). Useful as a
    vertex proxy on polytopes with hundreds of vertices, especially
    during animation playback where the mesh rebuilds every frame."""
    base = len(verts)
    for x in (-1, 1):
        for y in (-1, 1):
            for z in (-1, 1):
                verts.append((center[0] + x * size, center[1] + y * size,
                             center[2] + z * size))
    for f in ((0, 1, 3, 2), (4, 6, 7, 5), (0, 4, 5, 1),
             (2, 3, 7, 6), (0, 2, 6, 4), (1, 5, 7, 3)):
        faces.append([base + i for i in f])


def _newell(poly):
    n = [0.0, 0.0, 0.0]
    m = len(poly)
    for i in range(m):
        p, q = poly[i], poly[(i + 1) % m]
        n[0] += (p[1] - q[1]) * (p[2] + q[2])
        n[1] += (p[2] - q[2]) * (p[0] + q[0])
        n[2] += (p[0] - q[0]) * (p[1] + q[1])
    ln = math.sqrt(sum(t * t for t in n)) or 1.0
    return [t / ln for t in n]


def _leonardo_panels(F2, proj, origin, border, panel_thickness, taper,
                     scale):
    faces_o = []
    vsum = {}
    vcnt = {}
    for cyc in F2:
        poly = [proj[i][0] for i in cyc]
        n = _newell(poly)
        c = [sum(p[k] for p in poly) / len(poly) for k in range(3)]
        if sum(n[k] * (c[k] - origin[k]) for k in range(3)) < 0:
            cyc = list(reversed(cyc))
            poly = list(reversed(poly))
            n = [-t for t in n]
        faces_o.append((cyc, poly, n))
        for i in cyc:
            s = vsum.setdefault(i, [0.0, 0.0, 0.0])
            for k in range(3):
                s[k] += n[k]
            vcnt[i] = vcnt.get(i, 0) + 1
    vnormals = {}
    for (cyc, poly, n) in faces_o:
        for i in cyc:
            vnormals.setdefault(i, []).append(n)
    voff = {}
    for i, ns in vnormals.items():
        th = panel_thickness * scale * (proj[i][1] if taper else 1.0)
        M = [[0.0] * 3 for _ in range(3)]
        b = [0.0, 0.0, 0.0]
        for n in ns:
            for r in range(3):
                b[r] += n[r] * th
                for c in range(3):
                    M[r][c] += n[r] * n[c]
        lam = 1e-6 * (M[0][0] + M[1][1] + M[2][2] + 1e-12)
        for r in range(3):
            M[r][r] += lam
        det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
               - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
               + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
        if abs(det) < 1e-12:
            s = vsum[i]
            ln = math.sqrt(sum(t * t for t in s)) or 1.0
            voff[i] = [s[k] / ln * th for k in range(3)]
            continue
        m = []
        for c in range(3):
            Mc = [row[:] for row in M]
            for r in range(3):
                Mc[r][c] = b[r]
            dc = (Mc[0][0] * (Mc[1][1] * Mc[2][2] - Mc[1][2] * Mc[2][1])
                  - Mc[0][1] * (Mc[1][0] * Mc[2][2] - Mc[1][2] * Mc[2][0])
                  + Mc[0][2] * (Mc[1][0] * Mc[2][1] - Mc[1][1] * Mc[2][0]))
            m.append(dc / det)
        ln = math.sqrt(sum(t * t for t in m))
        cap = 3.0 * th
        if ln > cap:
            m = [t * cap / ln for t in m]
        voff[i] = m
    verts = []
    faces = []

    def _tri_n(a, b, c):
        u = [verts[b][k] - verts[a][k] for k in range(3)]
        v = [verts[c][k] - verts[a][k] for k in range(3)]
        n = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
             u[0] * v[1] - u[1] * v[0])
        ln = math.sqrt(sum(t * t for t in n)) or 1.0
        return [t / ln for t in n]

    def quad_tri(i0, i1, i2, i3):
        fold_a = sum(x * y for x, y in zip(_tri_n(i0, i1, i2),
                                           _tri_n(i0, i2, i3)))
        fold_b = sum(x * y for x, y in zip(_tri_n(i0, i1, i3),
                                           _tri_n(i1, i2, i3)))
        if fold_a >= fold_b:
            faces.append([i0, i1, i2])
            faces.append([i0, i2, i3])
        else:
            faces.append([i0, i1, i3])
            faces.append([i1, i2, i3])

    OUT = {}
    INN = {}
    for i in vsum:
        OUT[i] = len(verts)
        verts.append(proj[i][0])
        INN[i] = len(verts)
        verts.append(tuple(proj[i][0][k] - voff[i][k] for k in range(3)))
    rim_done = set()
    for (cyc, poly, n) in faces_o:
        m = len(cyc)
        c = [sum(p[k] for p in poly) / m for k in range(3)]
        avg = sum(proj[i][1] for i in cyc) / m
        th_f = panel_thickness * scale * (avg if taper else 1.0)
        HO = len(verts)
        for p in poly:
            verts.append(tuple(c[k] + (p[k] - c[k]) * (1 - border)
                               for k in range(3)))
        HI = len(verts)
        for p in poly:
            verts.append(tuple(c[k] + (p[k] - c[k]) * (1 - border)
                               - n[k] * th_f for k in range(3)))
        for i in range(m):
            j = (i + 1) % m
            a, b = cyc[i], cyc[j]
            faces.append([OUT[a], OUT[b], HO + j, HO + i])
            quad_tri(HI + i, HI + j, INN[b], INN[a])
            faces.append([HO + j, HO + i, HI + i, HI + j])
            key = (min(a, b), max(a, b))
            if key not in rim_done:
                rim_done.add(key)
                quad_tri(OUT[b], OUT[a], INN[a], INN[b])
    return verts, faces


FAMILIES = {
    'A4': ("[3,3,3] -- 5-cell", (3, 3, 3)),
    'B4': ("[4,3,3] -- tesseract / 16-cell", (4, 3, 3)),
    'F4': ("[3,4,3] -- 24-cell", (3, 4, 3)),
    'H4': ("[5,3,3] -- 120-cell / 600-cell", (5, 3, 3)),
    'D4': ("D4 / Jenn3D's \"Y-family\" -- demitesseract group (order 192)",
           (3, 2, 2, 3, 3, 2)),
}

GROUP_ORDER = {'A4': 120, 'B4': 384, 'F4': 1152, 'H4': 14400, 'D4': 192}


def _full_gram_marks(marks):
    """Normalise either a 3-mark linear chain or a 6-mark general
    diagram into the full (m12, m13, m14, m23, m24, m34) form."""
    if len(marks) == 3:
        p, q, r = marks
        return (p, 2, 2, q, 2, r)
    if len(marks) == 6:
        return tuple(marks)
    raise ValueError("marks must have length 3 (linear chain) or "
                     "6 (full Coxeter diagram)")


def _coxeter_mirrors(m):
    """Unit mirror normals for an n x n Coxeter exponent matrix m
    (m[i][i] = 1, m[i][j] = dihedral order between mirrors i and j),
    via a Cholesky factor of the implied Gram matrix. Rank-agnostic --
    used for both the rank-4 polychora and the rank-3 polyhedra."""
    if np is None:
        raise RuntimeError("the kaleidoscope needs NumPy")
    n = len(m)
    G = np.array([[1.0 if i == j else -math.cos(math.pi / m[i][j])
                   for j in range(n)] for i in range(n)])
    return np.linalg.cholesky(G)


def _coxeter_point(N, rings):
    """The generating point: off the ringed mirrors, on the others."""
    P = np.linalg.solve(N, np.array([1.0 if b else 0.0 for b in rings]))
    return P / np.linalg.norm(P)


def _coxeter_orbit(N, P, cap=200000):
    """Reflect the point through the mirrors until nothing new
    appears."""
    def key(v):
        return tuple(int(round(x * 1e7)) for x in v)
    n = len(N)
    refl = [np.eye(n) - 2.0 * np.outer(v, v) for v in N]
    seen = {key(P): P}
    frontier = [P]
    while frontier:
        nxt = []
        for v in frontier:
            for R in refl:
                w = R @ v
                k = key(w)
                if k not in seen:
                    seen[k] = w
                    nxt.append(w)
                    if len(seen) > cap:
                        raise RuntimeError("orbit too large")
        frontier = nxt
    return [tuple(float(x) for x in v) for v in seen.values()]


def wythoff_edges_exact(N, V, rings, v0_idx=0):
    """Multi-class edge finder: the base vertex v0 (index 0, the point
    the whole orbit was generated from) has one direct neighbour per
    RINGED mirror i, at distance d_i = |s_i(v0) - v0|. Since the group
    acts transitively on vertices, every vertex has the same set of
    edge lengths {d_i} realised among its own neighbours -- so once
    those lengths are known from v0, a spatial-hash search for each
    length across the whole vertex set finds every edge of every
    class. (A single "take the globally shortest pairs" filter only
    ever finds one class and silently drops the others whenever a
    polytope has more than one ringed node with unequal edge lengths
    -- this is what broke the p != q duoprisms earlier, and would
    just as quietly break any multi-ring uniform polytope.)"""
    n = len(N)
    A = np.asarray(V, dtype=float)
    v0 = A[v0_idx]
    refl = [np.eye(n) - 2.0 * np.outer(N[i], N[i]) for i in range(n)]
    lengths = []
    for i in range(n):
        if not rings[i]:
            continue
        d2 = float(np.dot(refl[i] @ v0 - v0, refl[i] @ v0 - v0))
        if d2 > 1e-9:
            lengths.append(d2)
    E = set()
    for d2 in lengths:
        cell = math.sqrt(d2) * 1.0000001
        tol = d2 * 1e-6 + 1e-9
        buckets = {}
        for idx, v in enumerate(A):
            buckets.setdefault(tuple(int(math.floor(x / cell)) for x in v),
                               []).append(idx)
        for c, members in buckets.items():
            near = []
            for off in itertools.product((-1, 0, 1), repeat=n):
                near.extend(buckets.get(tuple(c[k] + off[k]
                                              for k in range(n)), ()))
            for i2 in members:
                for j2 in near:
                    if j2 <= i2:
                        continue
                    dv = A[i2] - A[j2]
                    if abs(float(dv @ dv) - d2) < tol:
                        E.add((i2, j2))
    return sorted(E)


def wythoff_mirrors(marks):
    """Four unit normals realising a rank-4 Coxeter diagram, linear or
    branching."""
    m12, m13, m14, m23, m24, m34 = _full_gram_marks(marks)
    m = [[1, m12, m13, m14], [m12, 1, m23, m24],
         [m13, m23, 1, m34], [m14, m24, m34, 1]]
    return _coxeter_mirrors(m)


def wythoff_point(N, rings):
    return _coxeter_point(N, rings)


def wythoff_orbit(N, P, cap=200000):
    return _coxeter_orbit(N, P, cap)


def wythoff_build(family, rings):
    """(vertices, edges) of the uniform polychoron for this ringing."""
    marks = FAMILIES[family][1]
    N = wythoff_mirrors(marks)
    V = wythoff_orbit(N, wythoff_point(N, rings))
    E = wythoff_edges_exact(N, V, rings)
    return V, E


FAMILIES3 = {
    'A3': ("[3,3] -- tetrahedral (order 24)", (3, 3)),
    'B3': ("[4,3] -- cube / octahedral (order 48)", (4, 3)),
    'H3': ("[5,3] -- dodeca / icosahedral (order 120)", (5, 3)),
}


def _mirrors3(marks):
    p, q = marks
    return _coxeter_mirrors([[1, p, 2], [p, 1, q], [2, q, 1]])


def wythoff_build3(family, rings):
    """(vertices, edges) of the uniform 3D polyhedron for this
    ringing of a rank-3 Coxeter group."""
    N = _mirrors3(FAMILIES3[family][1])
    P = _coxeter_point(N, rings)
    V = _coxeter_orbit(N, P)
    E = wythoff_edges_exact(N, V, rings)
    return V, E


def _bisect_root(f, lo, hi, iters=200):
    """Real root of f on [lo, hi] via bisection (no scipy available)."""
    flo = f(lo)
    for _ in range(iters):
        mid = (lo + hi) / 2.0
        fm = f(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2.0


def _perm_parity(p):
    """0 for an even permutation, 1 for an odd one."""
    p = list(p)
    n = len(p)
    visited = [False] * n
    parity = 0
    for i in range(n):
        if visited[i]:
            continue
        j, clen = i, 0
        while not visited[j]:
            visited[j] = True
            j = p[j]
            clen += 1
        parity += clen - 1
    return parity % 2


def _dedupe_points(V, ndigits=5):
    seen = set()
    out = []
    for v in V:
        k = tuple(round(x, ndigits) for x in v)
        if k not in seen:
            seen.add(k)
            out.append(v)
    return out


def _min_distance_edges(V, tol=1e-3):
    """Edges as the vertex pairs at the minimum distance -- valid here
    because both snub solids are genuinely uniform (single edge
    length), unlike the multi-ring Wythoff cases elsewhere in this
    file. Works in any dimension (unlike polytope_edges, which is
    hardcoded to 4D)."""
    n = len(V)
    dim = len(V[0])
    d2min = None
    for i in range(n):
        for j in range(i + 1, n):
            d2 = sum((V[i][k] - V[j][k]) ** 2 for k in range(dim))
            if d2 > 1e-9 and (d2min is None or d2 < d2min):
                d2min = d2
    E = []
    for i in range(n):
        for j in range(i + 1, n):
            d2 = sum((V[i][k] - V[j][k]) ** 2 for k in range(dim))
            if abs(d2 - d2min) < d2min * tol:
                E.append((i, j))
    return E


def snub_cube_vertices():
    """The 24 vertices of the (chiral) snub cube: all permutations of
    (1, xi, 1/xi), with a sign pattern on each permutation whose
    parity is tied to the permutation's own parity -- where xi is the
    real root of xi^3 + xi^2 + xi = 1. Not reachable via simple
    Wythoff ringing (that only ever gives the non-chiral omnitruncated
    cuboctahedron already covered by TRUNC_CUBOCTA); this uses the
    solid's known closed-form coordinates instead."""
    xi = _bisect_root(lambda x: x ** 3 + x ** 2 + x - 1.0, 0.0, 1.0)
    base = (1.0, xi, 1.0 / xi)
    V = []
    for perm in itertools.permutations(range(3)):
        pp = _perm_parity(perm)
        v0 = [base[perm[i]] for i in range(3)]
        for signs in itertools.product((1, -1), repeat=3):
            nminus = sum(1 for s in signs if s < 0)
            if (nminus % 2) == pp:
                V.append(tuple(v0[i] * signs[i] for i in range(3)))
    return _dedupe_points(V)


def snub_dodecahedron_vertices():
    """The 60 vertices of the (chiral) snub dodecahedron: even
    permutations (with an even number of plus signs) of 5 vectors
    built from xi, the real root of xi^3 - 2*xi = phi. Closed-form
    coordinates, for the same reason as the snub cube above."""
    phi = PHI
    xi = _bisect_root(lambda x: x ** 3 - 2.0 * x - phi, 1.0, 2.0)
    a = xi - 1.0 / xi
    b = xi * phi + phi ** 2 + phi / xi
    base_vectors = [
        (2 * a, 2.0, 2 * b),
        (a + b / phi + phi, -a * phi + b + 1 / phi, a / phi + b * phi - 1),
        (-a / phi + b * phi + 1, -a + b / phi - phi, a * phi + b - 1 / phi),
        (-a / phi + b * phi - 1, a - b / phi - phi, a * phi + b + 1 / phi),
        (a + b / phi - phi, a * phi - b + 1 / phi, a / phi + b * phi + 1),
    ]
    V = []
    for base in base_vectors:
        for perm in itertools.permutations(range(3)):
            if _perm_parity(perm) != 0:
                continue
            v0 = [base[perm[i]] for i in range(3)]
            for signs in itertools.product((1, -1), repeat=3):
                nplus = sum(1 for s in signs if s > 0)
                if nplus % 2 == 0:
                    V.append(tuple(v0[i] * signs[i] for i in range(3)))
    return _dedupe_points(V)


ARCHIMEDEAN = {
    'TRUNC_TETRA': ('A3', (1, 1, 0)),
    'CUBOCTA': ('B3', (0, 1, 0)),
    'TRUNC_CUBE': ('B3', (1, 1, 0)),
    'TRUNC_OCTA': ('B3', (0, 1, 1)),
    'RHOMBICUBOCTA': ('B3', (1, 0, 1)),
    'TRUNC_CUBOCTA': ('B3', (1, 1, 1)),
    'ICOSIDODECA': ('H3', (0, 1, 0)),
    'TRUNC_DODECA': ('H3', (1, 1, 0)),
    'TRUNC_ICOSA': ('H3', (0, 1, 1)),
    'RHOMBICOSIDODECA': ('H3', (1, 0, 1)),
    'TRUNC_ICOSIDODECA': ('H3', (1, 1, 1)),
}

SNUB_SOLIDS = {
    'SNUB_CUBE': snub_cube_vertices,
    'SNUB_DODECA': snub_dodecahedron_vertices,
}


SEED_KINDS = ('TETRA', 'CUBE', 'OCTA', 'DODECA', 'ICOSA')


def hyperprism_from_edges(V3, E3, height=1.0):
    """Hyperprism extrusion directly from a 3D vertex/edge list --
    used for the Archimedean seeds, which mix more than one face size
    (e.g. triangles and octagons on a truncated cube), so they can't
    go through seed_poly's single-face-size cycle walker."""
    n = len(V3)
    V = ([tuple(list(v) + [-0.5 * height]) for v in V3]
         + [tuple(list(v) + [0.5 * height]) for v in V3])
    E = set()
    for a, b in E3:
        E.add((min(a, b), max(a, b)))
        E.add((min(a, b) + n, max(a, b) + n))
    for i in range(n):
        E.add((i, i + n))
    return V, sorted(E)


def hyperprism(seed='CUBE', height=1.0):
    if seed in SNUB_SOLIDS:
        V3 = SNUB_SOLIDS[seed]()
        E3 = _min_distance_edges(V3)
        r = math.sqrt(sum(x * x for x in V3[0]))
        V3 = [tuple(x / r for x in v) for v in V3]
        return hyperprism_from_edges(V3, E3, height)
    if seed in ARCHIMEDEAN:
        fam, rings = ARCHIMEDEAN[seed]
        V3, E3 = wythoff_build3(fam, rings)
        return hyperprism_from_edges(V3, E3, height)
    V3, F3 = seed_poly(seed, unit=True)
    n = len(V3)
    V = ([tuple(list(v) + [-0.5 * height]) for v in V3]
         + [tuple(list(v) + [0.5 * height]) for v in V3])
    E = set()
    for f in F3:
        for i in range(len(f)):
            a, b = f[i], f[(i + 1) % len(f)]
            E.add((min(a, b), max(a, b)))
            E.add((min(a, b) + n, max(a, b) + n))
    for i in range(n):
        E.add((i, i + n))
    return V, sorted(E)


def hyperpyramid(seed='TETRA', height=1.0):
    V3, F3 = seed_poly(seed, unit=True)
    n = len(V3)
    V = [tuple(list(v) + [0.0]) for v in V3] + [(0.0, 0.0, 0.0, height)]
    E = set()
    for f in F3:
        for i in range(len(f)):
            a, b = f[i], f[(i + 1) % len(f)]
            E.add((min(a, b), max(a, b)))
    for i in range(n):
        E.add((i, n))
    return V, sorted(E)


def duoprism_vertices(p, q):
    """The p,q-duoprism: the Cartesian product of a regular p-gon and
    a regular q-gon, embedded as (cos a, sin a, cos b, sin b) / sqrt(2)
    so every vertex lands on the unit 3-sphere -- the standard Clifford
    torus embedding, which is why duoprisms render as nested rings of
    circles under the curved (stereographic) style. Vertex (i, j) is
    at index i * q + j, matched by duoprism_edges below."""
    s = 1.0 / math.sqrt(2.0)
    V = []
    for i in range(p):
        a = 2.0 * math.pi * i / p
        ca, sa = math.cos(a), math.sin(a)
        for j in range(q):
            b = 2.0 * math.pi * j / q
            V.append((ca * s, sa * s, math.cos(b) * s, math.sin(b) * s))
    return V


def duoprism_edges(p, q):
    """Edges of the p,q-duoprism: each (i, j) connects to its neighbour
    around the p-gon and around the q-gon. NOT the same as the
    minimal-distance rule used elsewhere -- when p != q the two
    polygons have different edge lengths, so a single "shortest pair"
    filter would only find one of the two edge families."""
    def idx(i, j):
        return i * q + j
    E = set()
    for i in range(p):
        for j in range(q):
            E.add(tuple(sorted((idx(i, j), idx((i + 1) % p, j)))))
            E.add(tuple(sorted((idx(i, j), idx(i, (j + 1) % q)))))
    return sorted(E)


_RING_MASKS = [tuple(int(c) for c in format(m, '04b')) for m in range(1, 16)]


def wythoff_kind(family, rings):
    return "W:%s:%s" % (family, ''.join(str(b) for b in rings))


def build_polytope_ex(kind='CELL8', style='CURVED', proj_dist=1.05,
                      rot_xw=0.0, rot_yw=0.0, rot_zw=0.0, rot_xy=0.0,
                      arc_segments=12, radius=0.03, sides=6,
                      taper=True, vertex_spheres=True,
                      sphere_factor=1.6, scale=1.0, render='EDGES',
                      border=0.35, panel_thickness=0.03, half=False,
                      dual_compound=False, rings=0,
                      ring_cell_scale=0.9, rings_only=False,
                      vertex_shape='SPHERE'):
    """Builds the mesh data for one 4D-polytope framework, plus
    optional equatorial cutaway, dual compound and Hopf rings. Returns
    (verts, faces, face_mat, stats)."""
    systems = []
    if kind.startswith('P:'):
        _, which, seed, hgt = kind.split(':')
        V4, E = (hyperprism(seed, float(hgt)) if which == 'PRISM'
                 else hyperpyramid(seed, float(hgt)))
        F2 = None
        dual_compound = False
        rings = 0
    elif kind.startswith('D:'):
        _, ps, qs = kind.split(':')
        p, q = int(ps), int(qs)
        V4 = duoprism_vertices(p, q)
        E = duoprism_edges(p, q)
        F2 = (polytope_faces('DUOPRISM', V4, E,
                             cache_key='DUOPRISM:%d:%d' % (p, q))
              if render == 'LEONARDO' else None)
        dual_compound = False
        rings = 0
    elif kind.startswith('W:'):
        _, fam, bits = kind.split(':')
        V4, E = wythoff_build(fam, tuple(int(c) for c in bits))
        F2 = None
        dual_compound = False
        rings = 0
    else:
        V4 = polytope_vertices(kind)
        E = polytope_edges(V4)
        F2 = (polytope_faces(kind, V4, E) if render == 'LEONARDO' else None)
    if half:
        V4, E, F2 = _half_filter(V4, E, F2)
    systems.append((V4, E, F2))
    if dual_compound:
        D4, dkind, dkey = dual_vertices(kind)
        ED = polytope_edges(D4)
        FD = (polytope_faces(dkind, D4, ED, cache_key=dkey)
              if render == 'LEONARDO' else None)
        if style != 'CURVED':
            r = _cell_inradius(kind)
            D4 = [tuple(x * r for x in v) for v in D4]
        if half:
            D4, ED, FD = _half_filter(D4, ED, FD)
        systems.append((D4, ED, FD))
    cells = []
    if rings > 0 and kind == 'CELL120':
        cells = ring_cell_points(rings, ring_cell_scale, half)
    if rings_only and cells:
        systems = []
    systems = [(rotate4(V, rot_xw, rot_yw, rot_zw, rot_xy), Es, Fs)
               for (V, Es, Fs) in systems]
    cells = [(ri, rotate4(P, rot_xw, rot_yw, rot_zw, rot_xy))
             for (ri, P) in cells]
    if style == 'CURVED':
        dist = 1.0
        allv = ([v for (V, _e, _f) in systems for v in V]
                or [p for (_ri, P) in cells for p in P])
        pa, pb = _pole_angles(allv)
        if pa != 0.0 or pb != 0.0:
            systems = [(rotate4(V, pa, pb, 0.0, 0.0), Es, Fs)
                       for (V, Es, Fs) in systems]
            cells = [(ri, rotate4(P, pa, pb, 0.0, 0.0))
                     for (ri, P) in cells]
    else:
        dist = max(proj_dist, 1.001)
    verts = []
    faces = []
    edges = []
    face_mat = []
    curve_polylines = []
    for mi, (V, Es, Fs) in enumerate(systems):
        proj = {}
        for i, v in enumerate(V):
            p, s = project_point(v, dist)
            proj[i] = (tuple(c * scale for c in p), s)
        nf0 = len(faces)
        if render == 'LEONARDO':
            nvp = len(V)
            O = tuple(sum(proj[i][0][k] for i in range(nvp)) / nvp
                      for k in range(3))
            pv, pf = _leonardo_panels(Fs, proj, O, border,
                                      panel_thickness, taper, scale)
            base = len(verts)
            verts.extend(pv)
            faces.extend([[base + i for i in f] for f in pf])
        else:
            for (i, j) in Es:
                if style == 'CURVED':
                    pts = []
                    scls = []
                    for k in range(arc_segments + 1):
                        t = k / arc_segments
                        q = _slerp4(V[i], V[j], t)
                        p, s = project_point(q, dist)
                        pts.append(tuple(c * scale for c in p))
                        scls.append(s)
                else:
                    pts = [proj[i][0], proj[j][0]]
                    scls = [proj[i][1], proj[j][1]]
                    if arc_segments > 1:
                        a, b = pts
                        sa, sb = scls
                        pts = [tuple(a[k]
                                     + (b[k] - a[k]) * t / arc_segments
                                     for k in range(3))
                               for t in range(arc_segments + 1)]
                        scls = [sa + (sb - sa) * t / arc_segments
                                for t in range(arc_segments + 1)]
                if render == 'WIREFRAME':
                    base = len(verts)
                    verts.extend(pts)
                    for a in range(len(pts) - 1):
                        edges.append((base + a, base + a + 1))
                    continue
                if render == 'CURVE':
                    curve_polylines.append((pts, scls))
                    continue
                if taper:
                    radii = [radius * s * scale for s in scls]
                else:
                    radii = [radius * scale] * len(pts)
                strut_sides = 4 if render == 'EDGES' else sides
                add_strut(verts, faces, pts, radii, strut_sides)
            if render == 'BALLSTICK':
                nf_sticks = len(faces)
                face_mat.extend([-1] * (nf_sticks - nf0))
                for i in range(len(V)):
                    p, s = proj[i]
                    r = (radius * sphere_factor
                         * (s if taper else 1.0) * scale)
                    if vertex_shape == 'CUBE':
                        add_cube(verts, faces, p, r)
                    else:
                        add_sphere(verts, faces, p, r)
                face_mat.extend([-2] * (len(faces) - nf_sticks))
                continue
        face_mat.extend([mi] * (len(faces) - nf0))
    ring_base = len(systems)
    n_cells = 0
    for (ri, P) in cells:
        pts = [tuple(c * scale for c in project_point(p, dist)[0])
               for p in P]
        hv, hf = _hull_geometry(pts)
        base = len(verts)
        verts.extend(hv)
        faces.extend([[base + i for i in f] for f in hf])
        face_mat.extend([ring_base + ri] * len(hf))
        n_cells += 1
    stats = {'nv': (len(systems[0][0]) if systems else 0),
             'ne': (len(systems[0][1]) if systems else 0),
             'dual_nv': (len(systems[1][0])
                         if dual_compound and len(systems) > 1 else 0),
             'dual_ne': (len(systems[1][1])
                         if dual_compound and len(systems) > 1 else 0),
             'n_systems': len(systems), 'n_cells': n_cells,
             'n_rings': (min(rings, 12) if cells else 0),
             'wire_edges': edges, 'curve_polylines': curve_polylines}
    return verts, faces, face_mat, stats


def _hull_geometry(pts):
    bm = bmesh.new()
    for p in pts:
        bm.verts.new(p)
    bmesh.ops.convex_hull(bm, input=bm.verts[:])
    bm.verts.index_update()
    used = sorted({v.index for f in bm.faces for v in f.verts})
    remap = {vi: n for n, vi in enumerate(used)}
    bm.verts.ensure_lookup_table()
    hv = [tuple(bm.verts[vi].co) for vi in used]
    hf = [[remap[v.index] for v in f.verts] for f in bm.faces]
    bm.free()
    return hv, hf


def _make_material(name, rgba):
    """Fetch-or-create: a live engine rebuilds constantly, so
    materials are reused by name rather than minted fresh every time
    (which would otherwise pile up as Name.001, Name.002, ...)."""
    mat = bpy.data.materials.get(name)
    if mat is None:
        mat = bpy.data.materials.new(name)
        mat.use_nodes = True
    mat.diffuse_color = rgba
    node = mat.node_tree.nodes.get('Principled BSDF')
    if node is not None:
        node.inputs['Base Color'].default_value = rgba
    return mat


def _ring_color(ri):
    r, g, b = colorsys.hsv_to_rgb((ri / 12.0) % 1.0, 0.75, 0.9)
    return (r, g, b, 1.0)


def _get_or_create_curve_child(obj):
    """The managed Curve Child object for Style = Curve. Kept as a
    persistent reference on the PropertyGroup (rather than found by
    name each time) so a rename doesn't orphan it. Left in place,
    parented to the engine object, and simply emptied whenever the
    style switches away from Curve -- not deleted, since deleting it
    here on every non-Curve rebuild would mean losing it the moment
    you experiment with another style and come back."""
    props = obj.polytope4d
    child = props.curve_child
    if child is not None and child.name in bpy.data.objects:
        return child
    cdata = bpy.data.curves.new("Polytope4D Curve", type='CURVE')
    cdata.dimensions = '3D'
    child = bpy.data.objects.new("Curve Child", cdata)
    child.parent = obj
    for coll in obj.users_collection:
        coll.objects.link(child)
    props.curve_child = child
    return child


def _mesh_rebuild(obj):
    """Rebuild obj's mesh data in place from its polytope4d settings."""
    if obj is None or obj.type != 'MESH':
        return
    props = obj.polytope4d
    t0 = time.perf_counter()
    if props.form == 'UNIFORM':
        kind = wythoff_kind(props.family, props.wythoff_rings)
    elif props.form == 'PRISM':
        kind = "P:PRISM:%s:%.6f" % (props.seed, props.seed_height)
    elif props.form == 'DUOPRISM':
        kind = "D:%d:%d" % (props.duo_p, props.duo_q)
    else:
        kind = props.kind
    verts, faces, face_mat, st = build_polytope_ex(
        kind, props.style, props.proj_dist, props.rot_xw, props.rot_yw,
        props.rot_zw, props.rot_xy, props.arc_segments, props.radius,
        props.sides, props.taper, props.vertex_spheres,
        props.sphere_factor, props.scale, props.render, props.border,
        props.panel_thickness, props.half, props.dual_compound,
        props.rings, props.ring_cell_scale, props.rings_only,
        props.vertex_shape)
    curve_polylines = st.get('curve_polylines', [])
    bbox_pts = list(verts)
    for pts, _scls in curve_polylines:
        bbox_pts.extend(pts)
    if bbox_pts:
        xs = [v[0] for v in bbox_pts]
        ys = [v[1] for v in bbox_pts]
        zs = [v[2] for v in bbox_pts]
        cen = (0.5 * (min(xs) + max(xs)), 0.5 * (min(ys) + max(ys)),
               0.5 * (min(zs) + max(zs)))
        ext = max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))
        s = (2.0 * props.scale / ext) if ext > 1e-9 else 1.0
        verts = [((v[0] - cen[0]) * s, (v[1] - cen[1]) * s,
                  (v[2] - cen[2]) * s) for v in verts]
        curve_polylines = [
            ([((p[0] - cen[0]) * s, (p[1] - cen[1]) * s,
               (p[2] - cen[2]) * s) for p in pts], scls)
            for (pts, scls) in curve_polylines]
    me = obj.data
    me.clear_geometry()
    me.from_pydata(verts, st.get('wire_edges', []), faces)
    me.materials.clear()
    n_rings = st['n_rings']
    n_systems = st['n_systems']
    if props.render == 'BALLSTICK':
        me.materials.append(_make_material(
            "Polytope4D Sticks", (0.75, 0.78, 0.85, 1.0)))
        me.materials.append(_make_material(
            "Polytope4D Spheres", (0.25, 0.55, 0.85, 1.0)))
        for ri in range(n_rings):
            me.materials.append(_make_material(
                "Polytope4D Ring %d" % (ri + 1), _ring_color(ri)))
        mat_idx = [0 if m == -1 else 1 if m == -2
                  else 2 + (m - n_systems) for m in face_mat]
        me.polygons.foreach_set('material_index', mat_idx)
    elif n_systems > 1 or n_rings > 0:
        if n_systems >= 1:
            me.materials.append(_make_material(
                "Polytope4D Primal", (0.75, 0.78, 0.85, 1.0)))
        if n_systems > 1:
            me.materials.append(_make_material(
                "Polytope4D Dual", (0.9, 0.45, 0.15, 1.0)))
        for ri in range(n_rings):
            me.materials.append(_make_material(
                "Polytope4D Ring %d" % (ri + 1), _ring_color(ri)))
        me.polygons.foreach_set('material_index', face_mat)
    me.validate(clean_customdata=True)
    if props.render == 'BALLSTICK':
        if props.vertex_shape == 'CUBE':
            flags = [m == -1 for m in face_mat]
        else:
            flags = [m < n_systems for m in face_mat]
        me.polygons.foreach_set('use_smooth', flags[:len(me.polygons)])
    else:
        me.polygons.foreach_set('use_smooth', [False] * len(me.polygons))
    me.update()

    if props.render == 'CURVE':
        child = _get_or_create_curve_child(obj)
        cdata = child.data
        cdata.splines.clear()
        for pts, scls in curve_polylines:
            spline = cdata.splines.new('POLY')
            spline.points.add(len(pts) - 1)
            for i, pt in enumerate(pts):
                spline.points[i].co = (pt[0], pt[1], pt[2], 1.0)
                spline.points[i].radius = scls[i] if props.taper else 1.0
        cdata.bevel_depth = props.curve_bevel_depth
        cdata.bevel_resolution = props.curve_bevel_resolution
    elif props.curve_child is not None:
        try:
            props.curve_child.data.splines.clear()
        except (AttributeError, ReferenceError):
            pass

    props.last_build_time = time.perf_counter() - t0
    props.last_build_nv = st['nv']
    props.last_build_ne = st['ne']


def _on_prop_update(self, context):
    if self.auto_update:
        _mesh_rebuild(context.object)


def _on_auto_update_toggle(self, context):
    if self.auto_update:
        _mesh_rebuild(context.object)


@persistent
def _frame_change_handler(scene, _depsgraph=None):
    """Keyframed properties are written by the animation system on a
    path that does NOT fire a property's `update` callback (that only
    fires for direct UI edits or `setattr` in a script) -- a
    long-standing Blender limitation, not specific to this add-on. So
    every engine object with Auto Update on is rebuilt explicitly
    on every frame change (playback or manual scrub), rather than
    relying on the rotation sliders' own update callbacks."""
    for obj in scene.objects:
        if obj.type != 'MESH':
            continue
        props = getattr(obj, 'polytope4d', None)
        if props is not None and props.is_engine and props.auto_update:
            _mesh_rebuild(obj)


class POLYTOPE4D_PG_settings(bpy.types.PropertyGroup):
    is_engine: BoolProperty(default=False, options={'HIDDEN'})
    last_build_time: FloatProperty(default=0.0, options={'HIDDEN'})
    last_build_nv: IntProperty(default=0, options={'HIDDEN'})
    last_build_ne: IntProperty(default=0, options={'HIDDEN'})

    auto_update: BoolProperty(
        name="Auto Update", default=True,
        description="Rebuild the mesh immediately whenever a parameter "
                    "changes (including from keyframed animation)",
        update=_on_auto_update_toggle)

    kind: EnumProperty(
        name="Polytope",
        items=[('CELL5', "5-cell", "4-simplex: 5 vertices, 10 edges"),
               ('CELL8', "Tesseract (8-cell)", "16 vertices, 32 edges"),
               ('CELL16', "16-cell", "8 vertices, 24 edges"),
               ('CELL24', "24-cell", "24 vertices, 96 edges"),
               ('CELL120', "120-cell", "600 vertices, 1200 edges"),
               ('CELL600', "600-cell", "120 vertices, 720 edges")],
        default='CELL8',
        description="Which of the six regular convex 4-polytopes to "
                    "project into 3D", update=_on_prop_update)
    form: EnumProperty(
        name="Form",
        items=[('REGULAR', "Regular",
                "One of the six regular convex 4-polytopes"),
               ('UNIFORM', "Uniform (semi-regular)",
                "A vertex-transitive polytope from Wythoff's "
                "kaleidoscope: the truncations, rectifications and "
                "expansions of the regular ones"),
               ('PRISM', "Hyperprism",
                "A polyhedron translated along the fourth axis; the "
                "cube's hyperprism is the tesseract"),
               ('DUOPRISM', "Duoprism",
                "The Cartesian product of two regular polygons, "
                "sitting on a Clifford torus in the 3-sphere")],
        default='REGULAR',
        description="Build a regular polytope, a uniform one, a "
                    "prismatic family, or a duoprism",
        update=_on_prop_update)
    seed: EnumProperty(
        name="Base",
        items=[('TETRA', "Tetrahedron", ""), ('CUBE', "Cube", ""),
               ('OCTA', "Octahedron", ""),
               ('DODECA', "Dodecahedron", ""),
               ('ICOSA', "Icosahedron", ""),
               ('TRUNC_TETRA', "Truncated Tetrahedron", ""),
               ('CUBOCTA', "Cuboctahedron", ""),
               ('TRUNC_CUBE', "Truncated Cube", ""),
               ('TRUNC_OCTA', "Truncated Octahedron", ""),
               ('RHOMBICUBOCTA', "Rhombicuboctahedron", ""),
               ('TRUNC_CUBOCTA', "Truncated Cuboctahedron", ""),
               ('ICOSIDODECA', "Icosidodecahedron", ""),
               ('TRUNC_DODECA', "Truncated Dodecahedron", ""),
               ('TRUNC_ICOSA', "Truncated Icosahedron", ""),
               ('RHOMBICOSIDODECA', "Rhombicosidodecahedron", ""),
               ('TRUNC_ICOSIDODECA', "Truncated Icosidodecahedron", ""),
               ('SNUB_CUBE', "Snub Cube", ""),
               ('SNUB_DODECA', "Snub Dodecahedron", "")],
        default='CUBE',
        description="The 3-dimensional polyhedron the hyperprism is "
                    "raised on -- the 5 Platonic solids, all 13 "
                    "Archimedean solids (the 2 chiral snub solids use "
                    "closed-form coordinates rather than Wythoff "
                    "ringing, which can't reach them). Note: only "
                    "Struts / Ball and Stick / Wireframe support "
                    "these seeds -- Leonardo panels need this "
                    "add-on's face detector, which doesn't yet handle "
                    "a solid with more than one face size",
        update=_on_prop_update)
    seed_height: FloatProperty(
        name="Fourth-Axis Extent", default=1.0, min=0.05, max=10.0,
        description="How far the second copy is translated, or how "
                    "far the apex stands off, along the fourth axis",
        update=_on_prop_update)
    duo_p: IntProperty(
        name="Sides P", default=4, min=3, max=24,
        description="Number of sides of the first factor polygon",
        update=_on_prop_update)
    duo_q: IntProperty(
        name="Sides Q", default=3, min=3, max=24,
        description="Number of sides of the second factor polygon",
        update=_on_prop_update)
    family: EnumProperty(
        name="Symmetry",
        items=[(k, FAMILIES[k][0], "order %d" % GROUP_ORDER[k])
               for k in ('A4', 'B4', 'F4', 'H4', 'D4')],
        default='H4',
        description="Which rank-4 reflection group provides the "
                    "mirrors (uniform only)", update=_on_prop_update)
    wythoff_rings: EnumProperty(
        name="Ringed Nodes",
        items=[(''.join(str(b) for b in bits),
                ''.join(str(b) for b in bits),
                "Hold the generating point off mirror(s) %s"
                % ', '.join(str(i + 1) for i in range(4) if bits[i]))
               for bits in _RING_MASKS],
        default='1000',
        description="Which nodes of the Coxeter diagram are ringed",
        update=_on_prop_update)
    style: EnumProperty(
        name="Edges",
        items=[('CURVED', "Curved (stereographic)",
                "Vertices on the 3-sphere, edges as great-circle "
                "arcs, stereographically projected: circular arcs"),
               ('STRAIGHT', "Straight (perspective)",
                "Direct 4D perspective projection; small distance "
                "approaches a Schlegel diagram")],
        default='CURVED',
        description="Straight 4D perspective edges, or great-circle "
                    "arcs from stereographic projection",
        update=_on_prop_update)
    proj_dist: FloatProperty(
        name="Projection Distance", default=1.05, min=1.001, max=10.0,
        description="Eye distance along w for STRAIGHT edges",
        update=_on_prop_update)
    rot_xw: FloatProperty(name="Rotate XW", default=0.0,
                          min=-180.0, max=180.0,
                          description="Rotation in the XW plane before "
                          "projecting from 4D, in degrees",
                          update=_on_prop_update)
    rot_yw: FloatProperty(name="Rotate YW", default=0.0,
                          min=-180.0, max=180.0,
                          description="Rotation in the YW plane before "
                          "projecting from 4D, in degrees",
                          update=_on_prop_update)
    rot_zw: FloatProperty(name="Rotate ZW", default=0.0,
                          min=-180.0, max=180.0,
                          description="Rotation in the ZW plane before "
                          "projecting from 4D, in degrees",
                          update=_on_prop_update)
    rot_xy: FloatProperty(name="Rotate XY", default=0.0,
                          min=-180.0, max=180.0,
                          description="Rotation in the XY plane before "
                          "projecting from 4D, in degrees",
                          update=_on_prop_update)
    arc_segments: IntProperty(
        name="Arc Segments", default=12, min=1, max=100,
        description="Samples per edge (curved edges and tapering)",
        update=_on_prop_update)
    radius: FloatProperty(name="Strut Radius", default=0.03,
                          min=0.002, max=0.5, step=1, precision=3,
                          description="Radius of the edge struts",
                          update=_on_prop_update)
    sides: IntProperty(name="Strut Sides", default=6, min=3, max=60,
                       description="Cross-section sides of each round "
                       "strut (Ball and Stick style)",
                       update=_on_prop_update)
    taper: BoolProperty(
        name="Taper With Projection", default=True,
        description="Scale strut thickness by the local projection "
                    "factor (near-the-pole features fatter)",
        update=_on_prop_update)
    vertex_spheres: BoolProperty(name="Vertex Spheres", default=True,
                                 description="Place a sphere at each "
                                 "vertex", update=_on_prop_update)
    sphere_factor: FloatProperty(name="Sphere Size", default=1.6,
                                 min=1.0, max=4.0,
                                 description="Vertex marker size "
                                 "relative to the strut radius",
                                 update=_on_prop_update)
    vertex_shape: EnumProperty(
        name="Vertex Shape",
        items=[('SPHERE', "Sphere", "Round vertex marker"),
               ('CUBE', "Cube", "Low-poly cube vertex marker (6 faces "
                "vs. a sphere's several dozen) -- much cheaper on "
                "polytopes with hundreds of vertices, especially "
                "during animation playback where the mesh rebuilds "
                "every frame")],
        default='SPHERE',
        description="Shape of the vertex marker in Ball and Stick "
                    "style", update=_on_prop_update)
    render: EnumProperty(
        name="Style",
        items=[('EDGES', "Struts",
                "Solid tubes along the projected edges (no vertex "
                "spheres), following the curved stereographic arcs"),
               ('BALLSTICK', "Ball and Stick",
                "Struts along the projected edges with a sphere at "
                "every vertex"),
               ('WIREFRAME', "Wireframe",
                "Projected edges as a bare wireframe of polylines"),
               ('LEONARDO', "Leonardo (da Vinci)",
                "A flat open panel per 2D face of the polytope"),
               ('CURVE', "Curve",
                "Same edges as Wireframe, but output as a separate "
                "Curve object (a managed child named 'Curve Child') "
                "instead of mesh geometry -- gives you Blender's "
                "native, non-destructive Bevel Depth/Resolution "
                "instead of hand-built tube meshes. No vertex "
                "markers, same as Wireframe")],
        default='EDGES',
        description="How the projected framework is built",
        update=_on_prop_update)
    curve_bevel_depth: FloatProperty(
        name="Bevel Depth", default=0.03, min=0.0, max=0.5,
        step=1, precision=3,
        description="Radius of the Curve Child's rendered tube (0 = "
                    "a flat line, same look as Wireframe)",
        update=_on_prop_update)
    curve_bevel_resolution: IntProperty(
        name="Bevel Resolution", default=4, min=0, max=32,
        description="Roundness of the tube's cross-section -- higher "
                    "is rounder but adds more geometry",
        update=_on_prop_update)
    curve_child: PointerProperty(
        type=bpy.types.Object, options={'HIDDEN'},
        description="Internal reference to this engine's managed "
                    "Curve Child object -- not meant to be set by hand")
    border: FloatProperty(
        name="Border", default=0.06, min=0.005, max=1.0,
        description="Leonardo panel frame width (fraction of face "
                    "whatever its size)", update=_on_prop_update)
    panel_thickness: FloatProperty(
        name="Panel Thickness", default=0.03, min=0.002, max=0.5,
        step=1, precision=3,
        description="Thickness of the Leonardo face panels",
        update=_on_prop_update)
    scale: FloatProperty(name="Scale", default=1.0, min=0.01, max=100.0,
                         description="Overall size of the framework",
                         update=_on_prop_update)
    half: BoolProperty(
        name="Half (Cutaway)", default=False,
        description="Keep only the elements on one side of the "
                    "equatorial hyperplane (w <= 0, equator included) "
                    "before projection", update=_on_prop_update)
    dual_compound: BoolProperty(
        name="Dual Compound", default=False,
        description="Also build the dual polytope, its vertices at "
                    "this polytope's cell centers, in a second "
                    "material slot", update=_on_prop_update)
    rings: IntProperty(
        name="Hopf Rings", default=0, min=0, max=12,
        description="120-cell only: render N of the 12 rings of 10 "
                    "dodecahedral cells that partition the 120-cell "
                    "along Hopf fibers", update=_on_prop_update)
    ring_cell_scale: FloatProperty(
        name="Ring Cell Scale", default=0.9, min=0.1, max=1.0,
        description="Shrink factor of each ring cell toward its 4D "
                    "centroid (gaps make the rings legible)",
        update=_on_prop_update)
    rings_only: BoolProperty(
        name="Rings Only", default=False,
        description="Drop the edge framework and show just the "
                    "rings of solid cells", update=_on_prop_update)


class MESH_OT_polytope4d_add(bpy.types.Operator):
    """Add a new 4D Polytope engine (edit it afterwards from the "
    "N-panel sidebar)"""
    bl_idname = "mesh.polytope4d_add"
    bl_label = "4D Polytope"
    bl_options = {'UNDO'}

    def execute(self, context):
        me = bpy.data.meshes.new("Polytope4D")
        obj = bpy.data.objects.new("Polytope4D", me)
        context.collection.objects.link(obj)
        obj.location = context.scene.cursor.location
        for o in context.selected_objects:
            o.select_set(False)
        obj.select_set(True)
        context.view_layer.objects.active = obj
        obj.polytope4d.is_engine = True
        _mesh_rebuild(obj)
        self.report({'INFO'}, "%d vertices, %d edges"
                    % (obj.polytope4d.last_build_nv,
                       obj.polytope4d.last_build_ne))
        return {'FINISHED'}


class OBJECT_OT_polytope4d_rebuild(bpy.types.Operator):
    """Rebuild the mesh now (useful if Auto Update is off)"""
    bl_idname = "object.polytope4d_rebuild"
    bl_label = "Rebuild"
    bl_options = {'UNDO'}

    @classmethod
    def poll(cls, context):
        obj = context.object
        return (obj is not None and obj.type == 'MESH'
                and obj.polytope4d.is_engine)

    def execute(self, context):
        _mesh_rebuild(context.object)
        return {'FINISHED'}


class OBJECT_OT_polytope4d_reset_defaults(bpy.types.Operator):
    """Reset every parameter back to its default and rebuild"""
    bl_idname = "object.polytope4d_reset_defaults"
    bl_label = "Reset to Defaults"
    bl_options = {'UNDO'}

    @classmethod
    def poll(cls, context):
        obj = context.object
        return (obj is not None and obj.type == 'MESH'
                and obj.polytope4d.is_engine)

    def execute(self, context):
        props = context.object.polytope4d
        for key in props.bl_rna.properties.keys():
            if key in ('rna_type', 'is_engine', 'curve_child'):
                continue
            try:
                props.property_unset(key)
            except Exception:
                pass
        _mesh_rebuild(context.object)
        return {'FINISHED'}


class OBJECT_OT_polytope4d_detach(bpy.types.Operator):
    """Freeze the mesh as-is and hide the engine panel (the mesh "
    "itself is untouched, only its live link is removed)"""
    bl_idname = "object.polytope4d_detach"
    bl_label = "Detach from Engine"
    bl_options = {'UNDO'}

    @classmethod
    def poll(cls, context):
        obj = context.object
        return (obj is not None and obj.type == 'MESH'
                and obj.polytope4d.is_engine)

    def execute(self, context):
        context.object.polytope4d.is_engine = False
        return {'FINISHED'}


class VIEW3D_PT_polytope4d(bpy.types.Panel):
    bl_label = "4D Polytope"
    bl_idname = "VIEW3D_PT_polytope4d"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "4D Polytope"


    def draw(self, context):
        obj = context.object
        lay = self.layout
        if (obj is None or obj.type != 'MESH'
                or not obj.polytope4d.is_engine):
            lay.label(text="No 4D Polytope engine selected")
            lay.operator('mesh.polytope4d_add', icon='MESH_CUBE',
                        text="Add 4D Polytope")
            return
        props = obj.polytope4d
        lay.use_property_split = True

        lay.prop(props, 'form')
        if props.form == 'UNIFORM':
            lay.prop(props, 'family')
            lay.prop(props, 'wythoff_rings')
        elif props.form == 'PRISM':
            lay.prop(props, 'seed')
            lay.prop(props, 'seed_height')
        elif props.form == 'DUOPRISM':
            lay.prop(props, 'duo_p')
            lay.prop(props, 'duo_q')
        else:
            lay.prop(props, 'kind')
        lay.prop(props, 'style')
        if props.style == 'STRAIGHT':
            lay.prop(props, 'proj_dist')

        lay.prop(props, 'render')
        if props.render == 'LEONARDO':
            lay.prop(props, 'border')
            lay.prop(props, 'panel_thickness')
            lay.prop(props, 'taper')
            row = lay.row(align=True)
            row.label(text="", icon='TRANSFORM_ORIGINS')
            row.prop(props, 'scale')
        elif props.render == 'WIREFRAME':
            row = lay.row(align=True)
            row.label(text="", icon='PARTICLE_PATH')
            row.prop(props, 'arc_segments')
            row = lay.row(align=True)
            row.label(text="", icon='TRANSFORM_ORIGINS')
            row.prop(props, 'scale')
        elif props.render == 'CURVE':
            row = lay.row(align=True)
            row.label(text="", icon='PARTICLE_PATH')
            row.prop(props, 'arc_segments')
            row = lay.row(align=True)
            row.label(text="", icon='GIZMO')
            row.prop(props, 'curve_bevel_depth')
            row = lay.row(align=True)
            row.label(text="", icon='ALIGN_JUSTIFY')
            row.prop(props, 'curve_bevel_resolution')
            lay.prop(props, 'taper')
            row = lay.row(align=True)
            row.label(text="", icon='TRANSFORM_ORIGINS')
            row.prop(props, 'scale')
        else:
            row = lay.row(align=True)
            row.label(text="", icon='PARTICLE_PATH')
            row.prop(props, 'arc_segments')
            row = lay.row(align=True)
            row.label(text="", icon='GIZMO')
            row.prop(props, 'radius')
            if props.render == 'BALLSTICK':
                row = lay.row(align=True)
                row.label(text="", icon='ALIGN_JUSTIFY')
                row.prop(props, 'sides')
            lay.prop(props, 'taper')
            if props.render == 'BALLSTICK':
                row = lay.row(align=True)
                row.label(text="", icon='PRESET_NEW')
                row.prop(props, 'sphere_factor')
                row = lay.row(align=True)
                row.label(text="", icon='SURFACE_NCURVE')
                row.prop(props, 'vertex_shape')
            row = lay.row(align=True)
            row.label(text="", icon='TRANSFORM_ORIGINS')
            row.prop(props, 'scale')

        col = lay.column(align=True)
        for k in ('rot_xw', 'rot_yw', 'rot_zw', 'rot_xy'):
            row = col.row(align=True)
            row.label(text="", icon='ORIENTATION_GIMBAL')
            row.prop(props, k)

        lay.separator()
        col = lay.column(align=True)
        col.prop(props, 'half')
        sub = col.column(align=True)
        sub.enabled = (props.form == 'REGULAR')
        sub.prop(props, 'dual_compound')
        rowr = sub.row(align=True)
        is_120 = (props.form == 'REGULAR' and props.kind == 'CELL120')
        rowr.enabled = is_120
        rowr.prop(props, 'rings')
        if is_120 and props.rings > 0:
            col.prop(props, 'ring_cell_scale')
            col.prop(props, 'rings_only')

        lay.separator()
        row = lay.row(align=True)
        row.prop(props, 'auto_update', toggle=True)
        row.operator('object.polytope4d_rebuild', icon='FILE_REFRESH')
        lay.label(text="Last build %.2fs -- %d vertices, %d edges"
                  % (props.last_build_time, props.last_build_nv,
                     props.last_build_ne))

        lay.separator()
        lay.operator('object.polytope4d_reset_defaults', icon='LOOP_BACK')
        lay.operator('object.polytope4d_detach', icon='UNLINKED')


def _menu_func(self, context):
    self.layout.operator("mesh.polytope4d_add", icon='MESH_CUBE')


classes = (
    POLYTOPE4D_PG_settings,
    MESH_OT_polytope4d_add,
    OBJECT_OT_polytope4d_rebuild,
    OBJECT_OT_polytope4d_reset_defaults,
    OBJECT_OT_polytope4d_detach,
    VIEW3D_PT_polytope4d,
)


def register():
    for cls in classes:
        bpy.utils.register_class(cls)
    bpy.types.Object.polytope4d = PointerProperty(type=POLYTOPE4D_PG_settings)
    bpy.types.VIEW3D_MT_mesh_add.append(_menu_func)
    if _frame_change_handler not in bpy.app.handlers.frame_change_post:
        bpy.app.handlers.frame_change_post.append(_frame_change_handler)


def unregister():
    if _frame_change_handler in bpy.app.handlers.frame_change_post:
        bpy.app.handlers.frame_change_post.remove(_frame_change_handler)
    bpy.types.VIEW3D_MT_mesh_add.remove(_menu_func)
    del bpy.types.Object.polytope4d
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)


if __name__ == "__main__":
    register()
