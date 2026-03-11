"""Pure NumPy implementation of pyvdwsurface for RESP fallback."""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Sequence

import numpy as np

TORAD = math.pi / 180.0
DP_TOL = 0.001
N_REFINE_STEPS = 25
R_H = math.sqrt(1.0 - 2.0 * math.cos(72.0 * TORAD)) / (1.0 - math.cos(72.0 * TORAD))
DOWNSAMPLE_SEED = 0

BONDI_RADII = {
    "H": 1.1,
    "Li": 2.2,
    "Be": 1.9,
    "B": 1.8,
    "C": 1.7,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "Na": 2.4,
    "Mg": 2.2,
    "Al": 2.1,
    "Si": 2.1,
    "P": 1.8,
    "S": 1.8,
    "Cl": 1.75,
    "K": 2.8,
    "Ca": 2.31,
    "Sc": 2.3,
    "Ti": 2.15,
    "V": 2.05,
    "Cr": 2.05,
    "Mn": 2.05,
    "Fe": 2.05,
    "Co": 2.0,
    "Ni": 2.0,
    "Cu": 2.0,
    "Zn": 2.1,
    "Ga": 2.1,
    "As": 2.05,
    "Se": 1.9,
    "Br": 1.83,
    "Rb": 2.9,
    "Sr": 2.55,
    "Y": 2.4,
    "Zr": 2.3,
    "Nb": 2.15,
    "Mo": 2.1,
    "Tc": 2.05,
    "Ru": 2.05,
    "Rh": 2.0,
    "Pd": 2.05,
    "Ag": 2.1,
    "Cd": 2.2,
    "In": 2.2,
    "I": 1.98,
    "W": 2.1,
    "Os": 2.1,
    "Pb": 2.02,
    "Ir": 2.1,
    "La": 1.87,
}


def _vnorm(v: np.ndarray) -> float:
    return float(np.sqrt(np.dot(v, v)))


def _divarc(xyz1: np.ndarray, xyz2: np.ndarray, div1: int, div2: int) -> np.ndarray:
    xd = xyz1[1] * xyz2[2] - xyz2[1] * xyz1[2]
    yd = xyz1[2] * xyz2[0] - xyz2[2] * xyz1[0]
    zd = xyz1[0] * xyz2[1] - xyz2[0] * xyz1[1]
    dd = math.sqrt(xd * xd + yd * yd + zd * zd)
    if dd < DP_TOL:
        raise RuntimeError("_divarc: rotation axis of length")

    d1 = _vnorm(xyz1)
    if d1 < math.sqrt(0.5):
        raise RuntimeError("_divarc: vector 1 of sq.length too small")

    d2 = _vnorm(xyz2)
    if d2 < math.sqrt(0.5):
        raise RuntimeError("_divarc: vector 2 of sq.length too small")

    phi = math.sin(dd / math.sqrt(d1 * d2))
    phi = phi * float(div1) / float(div2)
    sphi = math.sin(phi)
    cphi = math.cos(phi)

    s = (xyz1[0] * xd + xyz1[0] * yd + xyz1[2] * zd) / dd

    x = xd * s * (1.0 - cphi) / dd + xyz1[0] * cphi + (yd * xyz1[2] - xyz1[1] * zd) * sphi / dd
    y = yd * s * (1.0 - cphi) / dd + xyz1[1] * cphi + (zd * xyz1[0] - xyz1[2] * xd) * sphi / dd
    z = zd * s * (1.0 - cphi) / dd + xyz1[2] * cphi + (xd * xyz1[1] - xyz1[0] * yd) * sphi / dd
    out = np.array([x, y, z], dtype=np.float64)
    return out / _vnorm(out)


def _icosahedron_vertices() -> list[np.ndarray]:
    rg = math.cos(72.0 * TORAD) / (1.0 - math.cos(72.0 * TORAD))
    return [
        np.array([0.0, 0.0, 1.0], dtype=np.float64),
        np.array([R_H * math.cos(72.0 * TORAD), R_H * math.sin(72.0 * TORAD), rg], dtype=np.float64),
        np.array([R_H * math.cos(144.0 * TORAD), R_H * math.sin(144.0 * TORAD), rg], dtype=np.float64),
        np.array([R_H * math.cos(216.0 * TORAD), R_H * math.sin(216.0 * TORAD), rg], dtype=np.float64),
        np.array([R_H * math.cos(288.0 * TORAD), R_H * math.sin(288.0 * TORAD), rg], dtype=np.float64),
        np.array([R_H, 0.0, rg], dtype=np.float64),
        np.array([R_H * math.cos(36.0 * TORAD), R_H * math.sin(36.0 * TORAD), -rg], dtype=np.float64),
        np.array([R_H * math.cos(108.0 * TORAD), R_H * math.sin(108.0 * TORAD), -rg], dtype=np.float64),
        np.array([-R_H, 0.0, -rg], dtype=np.float64),
        np.array([R_H * math.cos(252.0 * TORAD), R_H * math.sin(252.0 * TORAD), -rg], dtype=np.float64),
        np.array([R_H * math.cos(324.0 * TORAD), R_H * math.sin(324.0 * TORAD), -rg], dtype=np.float64),
        np.array([0.0, 0.0, -1.0], dtype=np.float64),
    ]


def _dotsphere1(density: int) -> list[np.ndarray]:
    a = math.sqrt((density - 2.0) / 10.0)
    tess = int(math.ceil(a))

    vertices = _icosahedron_vertices()

    if tess > 1:
        a = R_H * R_H * 2.0 * (1.0 - math.cos(72.0 * TORAD))
        for i in range(11):
            for j in range(i + 1, 12):
                d = _vnorm(vertices[i] - vertices[j])
                if abs(a - d * d) > DP_TOL:
                    continue
                for tl in range(tess):
                    vertices.append(_divarc(vertices[i], vertices[j], tl, tess))

    for i in range(10):
        for j in range(i + 1, 11):
            d = _vnorm(vertices[i] - vertices[j])
            if abs(a - d * d) > DP_TOL:
                continue

            for k in range(j + 1, 12):
                d_ik = _vnorm(vertices[i] - vertices[k])
                d_jk = _vnorm(vertices[j] - vertices[k])
                if abs(a - d_ik * d_ik) > DP_TOL or abs(a - d_jk * d_jk) > DP_TOL:
                    continue
                for tl in range(1, tess - 1):
                    ji = _divarc(vertices[j], vertices[i], tl, tess)
                    ki = _divarc(vertices[k], vertices[i], tl, tess)

                    for tl2 in range(1, tess - tl):
                        ij = _divarc(vertices[i], vertices[j], tl2, tess)
                        kj = _divarc(vertices[k], vertices[j], tl2, tess)
                        ik = _divarc(vertices[i], vertices[k], tess - tl - tl2, tess)
                        jk = _divarc(vertices[j], vertices[k], tess - tl - tl2, tess)

                        xyz1 = _divarc(ki, ji, tl2, tess - tl)
                        xyz2 = _divarc(kj, ij, tl, tess - tl2)
                        xyz3 = _divarc(jk, ik, tl, tl + tl2)

                        x = xyz1 + xyz2 + xyz3
                        vertices.append(x / _vnorm(x))

    return vertices


def _dotsphere2(density: int) -> list[np.ndarray]:
    a = math.sqrt((density - 2.0) / 30.0)
    tess = max(int(math.ceil(a)), 1)
    vertices = _icosahedron_vertices()

    a = R_H * R_H * 2.0 * (1.0 - math.cos(72.0 * TORAD))

    for i in range(10):
        for j in range(i + 1, 11):
            d = _vnorm(vertices[i] - vertices[j])
            if abs(a - d * d) > DP_TOL:
                continue
            for k in range(j + 1, 12):
                d_ik = _vnorm(vertices[i] - vertices[k])
                d_jk = _vnorm(vertices[j] - vertices[k])
                if abs(a - d_ik * d_ik) > DP_TOL or abs(a - d_jk * d_jk) > DP_TOL:
                    continue

                x = vertices[i] + vertices[j] + vertices[k]
                vertices.append(x / _vnorm(x))

    if tess > 1:
        adod = 4.0 * (math.cos(108.0 * TORAD) - math.cos(120.0 * TORAD)) / (1.0 - math.cos(120.0 * TORAD))
        ai_d = 2.0 * (1.0 - math.sqrt(1.0 - a / 3.0))

        for i in range(31):
            j1 = 12
            j2 = 32
            aa = ai_d
            if i > 12:
                j1 = i + 1
                aa = adod
            for j in range(j1, j2):
                d = _vnorm(vertices[i] - vertices[j])
                if abs(aa - d * d) > DP_TOL:
                    continue
                for tl in range(1, tess):
                    vertices.append(_divarc(vertices[i], vertices[j], tl, tess))

        for i in range(12):
            for j in range(12, 31):
                d = _vnorm(vertices[i] - vertices[j])
                if abs(ai_d - d * d) > DP_TOL:
                    continue

                for k in range(j + 1, 32):
                    d_ik = _vnorm(vertices[i] - vertices[k])
                    d_jk = _vnorm(vertices[j] - vertices[k])
                    if abs(ai_d - d_ik * d_ik) > DP_TOL or abs(adod - d_jk * d_jk) > DP_TOL:
                        continue
                    for tl in range(1, tess - 1):
                        ji = _divarc(vertices[j], vertices[i], tl, tess)
                        ki = _divarc(vertices[k], vertices[i], tl, tess)
                        for tl2 in range(1, tess - tl):
                            ij = _divarc(vertices[i], vertices[j], tl2, tess)
                            kj = _divarc(vertices[k], vertices[j], tl2, tess)
                            ik = _divarc(vertices[i], vertices[k], tess - tl - tl2, tess)
                            jk = _divarc(vertices[j], vertices[k], tess - tl - tl2, tess)

                            xyz1 = _divarc(ki, ji, tl2, tess - tl)
                            xyz2 = _divarc(kj, ij, tl, tess - tl2)
                            xyz3 = _divarc(jk, ik, tl, tl + tl2)

                            x = xyz1 + xyz2 + xyz3
                            vertices.append(x / _vnorm(x))
    return vertices


def _get_coulomb_energy(points: np.ndarray) -> float:
    e = 0.0
    for i in range(points.shape[0]):
        d = points[i + 1 :] - points[i]
        if d.size == 0:
            continue
        l = np.linalg.norm(d, axis=1)
        mask = l >= 1e-8
        if not np.any(mask):
            continue
        e += float(np.sum(1.0 / l[mask]))
    return e


def _get_coulomb_forces(points: np.ndarray) -> np.ndarray:
    n = points.shape[0]
    forces = np.zeros_like(points)
    for i in range(n):
        r = points[i] - points[i + 1 :]
        if r.size == 0:
            continue
        l = np.linalg.norm(r, axis=1)
        mask = l >= 1e-8
        if not np.any(mask):
            continue
        rr = r[mask]
        ll = l[mask]
        ff = rr / (ll[:, None] ** 3)
        forces[i] += np.sum(ff, axis=0)
        idx = np.where(mask)[0] + (i + 1)
        forces[idx] -= ff
    return forces


def _rand_r(seed: int) -> tuple[int, int]:
    next_ = int(seed) & 0xFFFFFFFF
    next_ = (next_ * 1103515245 + 12345) & 0xFFFFFFFF
    result = (next_ // 65536) % 2048

    next_ = (next_ * 1103515245 + 12345) & 0xFFFFFFFF
    result = (result << 10) ^ ((next_ // 65536) % 1024)

    next_ = (next_ * 1103515245 + 12345) & 0xFFFFFFFF
    result = (result << 10) ^ ((next_ // 65536) % 1024)

    return next_, int(result)


def _sample_indices_like_c_rand_r(total: int, n_keep: int, seed: int = DOWNSAMPLE_SEED) -> np.ndarray:
    keep: set[int] = set()
    s = int(seed)
    while len(keep) < n_keep:
        s, rv = _rand_r(s)
        keep.add(rv % total)
    return np.array(sorted(keep), dtype=np.int64)


def _refine_dotsphere(points: np.ndarray) -> np.ndarray:
    step = 0.005
    if points.shape[0] > 100:
        step /= 50.0

    e0 = _get_coulomb_energy(points)
    for _ in range(N_REFINE_STEPS):
        forces = _get_coulomb_forces(points)
        points = points + forces * step
        points = points / np.linalg.norm(points, axis=1, keepdims=True)
        e = _get_coulomb_energy(points)
        if e0 < e:
            step /= 2.0
        e0 = e
        if step < 1e-8:
            break
    return points


@lru_cache(maxsize=512)
def _dotsphere(density: int) -> np.ndarray:
    i1 = 1
    i2 = 1

    while 10 * i1 * i1 + 2 < density:
        i1 += 1

    while 30 * i2 * i2 + 2 < density:
        i2 += 1

    if 10 * i1 * i1 - 2 < 30 * i2 * i2 - 2:
        vertices = _dotsphere1(density)
    else:
        vertices = _dotsphere2(density)

    if density < len(vertices):
        keep_idx = _sample_indices_like_c_rand_r(len(vertices), density, seed=DOWNSAMPLE_SEED)
        vertices = [vertices[int(i)] for i in keep_idx]

    arr = np.array(vertices, dtype=np.float64)
    arr = _refine_dotsphere(arr)
    return arr


def _normalize_element(el: object) -> str:
    if isinstance(el, bytes):
        return el.decode("utf-8")
    return str(el)


def vdwsurface(
    coordinates: np.ndarray,
    elements: Sequence[object],
    scale_factor: float = 1.0,
    density: float = 1.0,
) -> np.ndarray:
    coords = np.asarray(coordinates, dtype=np.float64)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("coordinates must have shape (n_atoms, 3)")
    if len(elements) != coords.shape[0]:
        raise ValueError("len(elements) must match number of coordinates")

    symbols = [_normalize_element(e) for e in elements]
    radii = np.empty(coords.shape[0], dtype=np.float64)
    for i, sym in enumerate(symbols):
        if sym not in BONDI_RADII:
            raise ValueError(f"{sym} is not a supported element")
        radii[i] = BONDI_RADII[sym] * scale_factor

    surface = []
    for i in range(coords.shape[0]):
        target_n = int(density * ((4.0 / 3.0) * math.pi * radii[i] ** 3))
        if target_n < 1:
            target_n = 1
        dots = _dotsphere(target_n) * radii[i] + coords[i]

        diff = coords - coords[i]
        d = np.linalg.norm(diff, axis=1)
        mask = np.arange(coords.shape[0]) != i
        neigh = np.where(mask & (d < (radii[i] + radii)))[0]
        if neigh.size:
            neigh = neigh[np.argsort(d[neigh])]

        accessible = np.ones(dots.shape[0], dtype=bool)
        for j in neigh:
            dd = np.linalg.norm(coords[j] - dots, axis=1)
            accessible &= dd >= radii[j]
            if not np.any(accessible):
                break

        if np.any(accessible):
            surface.append(dots[accessible])

    if not surface:
        return np.zeros((0, 3), dtype=np.float64)
    return np.vstack(surface)


__all__ = ["vdwsurface", "BONDI_RADII", "DOWNSAMPLE_SEED"]
