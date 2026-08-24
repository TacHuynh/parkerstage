#!/usr/bin/env python3
"""Slice meshes to understand the base/carriage cross-sections, in link frames.

Prints ASCII occupancy maps (even-odd ray fill, same method as mesh_inertia)
for any link's mesh, sliced at a chosen coordinate of a chosen frame.
"""
import sys
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import rebuild_urdf as R  # noqa: E402
from mesh_inertia import read_triangles  # noqa: E402

MESH_DIR = os.path.join(ROOT, "meshes")

RES = 0.5e-3


def tris_in(link_name, frame_link):
    """Triangles of link's visual mesh, in `frame_link`'s frame."""
    import xml.etree.ElementTree as ET
    t = ET.parse(os.path.join(ROOT, "urdf", "parkerstage.urdf"))
    by = {l.get("name"): l for l in t.getroot().findall("link")}
    l = by[link_name]
    vo = l.find("visual/origin")
    oxyz = [float(x) for x in vo.get("xyz").split()]
    fn = l.find("visual/geometry/mesh").get("filename").split("/")[-1]
    tris = read_triangles(os.path.join(MESH_DIR, fn))
    # mesh frame -> link frame (visual origin, no rotation in this export)
    T_link = R.t4([[1, 0, 0], [0, 1, 0], [0, 0, 1]], oxyz)
    # link frame -> frame_link frame
    T_fr = R.mul(R.inv(R.poses[frame_link]), R.poses[link_name])
    out = []
    for a, b, c in tris:
        out.append(tuple(R.vmul(T_fr, R.vmul(T_link, p)) for p in (a, b, c)))
    return out


def slice_map(tris, axis, value, res=RES):
    """Occupancy map in the plane perpendicular to `axis` at `value`.

    axis: 0=x, 1=y, 2=z. Returns (rows, lo0, hi0, lo2, hi2, res) where
    rows[i][j] is True if cell (j, i) is inside the solid; rows are indexed
    top-down along the i2 axis. Each cell is filled by ray-casting along
    `axis` through the cell centre and checking the even-odd span covers 0.
    """
    i0, i2 = [a for a in range(3) if a != axis]
    pad = 2 * res
    lo0 = min(p[i0] for t in tris for p in t) - pad
    hi0 = max(p[i0] for t in tris for p in t) + pad
    lo2 = min(p[i2] for t in tris for p in t) - pad
    hi2 = max(p[i2] for t in tris for p in t) + pad
    n0 = int((hi0 - lo0) / res) + 1
    n2 = int((hi2 - lo2) / res) + 1
    # bucket every triangle by its (y, z) footprint so x_span is O(bucket)
    lo1 = min(p[1] for t in tris for p in t) - pad
    hi1 = max(p[1] for t in tris for p in t) + pad
    ny = int((hi1 - lo1) / res) + 1
    buckets = {}
    for t in tris:
        y0 = int((min(p[1] for p in t) - lo1) / res)
        y1 = int((max(p[1] for p in t) - lo1) / res)
        z0 = int((min(p[2] for p in t) - lo2) / res)
        z1 = int((max(p[2] for p in t) - lo2) / res)
        for yc in range(max(0, y0), min(ny - 1, y1) + 1):
            for zc in range(max(0, z0), min(n2 - 1, z1) + 1):
                buckets.setdefault((yc, zc), []).append(t)

    # Rays are always cast along +x (same as mesh_inertia, validated against
    # the cube and the real meshes): for an extrusion along y the surface
    # triangles are parallel to a y-ray, so a y-cast finds nothing.
    def x_span(y, z):
        """Even-odd x-intervals of the solid at (y, z), along the +x ray."""
        yc = max(0, min(ny - 1, int((y - lo1) / res)))
        zc = max(0, min(n2 - 1, int((z - lo2) / res)))
        ts = []
        for (a, b, c) in buckets.get((yc, zc), []):
            e1 = tuple(b[k] - a[k] for k in range(3))
            e2 = tuple(c[k] - a[k] for k in range(3))
            det = e1[2] * e2[1] - e1[1] * e2[2]  # (e2 x e1) . D with D=(1,0,0)
            if abs(det) < 1e-12:
                continue
            inv = 1.0 / det
            T = (-a[0], y - a[1], z - a[2])
            P = (0.0, -e2[2], e2[1])
            u = (T[0] * P[0] + T[1] * P[1] + T[2] * P[2]) * inv
            if u < 0.0 or u > 1.0:
                continue
            Q = (T[1] * e1[2] - T[2] * e1[1],
                 T[2] * e1[0] - T[0] * e1[2],
                 T[0] * e1[1] - T[1] * e1[0])
            vv = Q[0] * inv
            if vv < 0.0 or u + vv > 1.0:
                continue
            ts.append((e2[0] * Q[0] + e2[1] * Q[1] + e2[2] * Q[2]) * inv)
        ts.sort()
        out = []
        k = 0
        while k + 1 < len(ts):
            if ts[k + 1] - ts[k] > 0.5 * res:
                out.append((ts[k], ts[k + 1]))
            k += 2
        return out

    def inside(jc, zc):
        """True if the cell centre lies inside the solid."""
        v0 = lo0 + (jc + 0.5) * res  # coordinate along i0
        v2 = lo2 + (zc + 0.5) * res  # coordinate along i2
        # ray through (0, y, z): y/z depend on which plane we sliced
        if axis == 0:      # x = const: column is (y, z); test x_slice in span
            y, z = v0, v2
            probe = value
        elif axis == 1:    # y = const: column is (x, z); test x_cell in span
            y, z = value, v2
            probe = v0
        else:              # z = const: column is (x, y); test x_cell in span
            y, z = v0, value
            probe = v0
        for (x0, x1) in x_span(y, z):
            if x0 <= probe <= x1:
                return True
        return False

    rows = []
    for zc in range(n2 - 1, -1, -1):
        rows.append([inside(jc, zc) for jc in range(n0)])
    return rows, lo0, hi0, lo2, hi2, res


def print_slice(tris, axis, value, label=""):
    names = "xyz"
    rows, lo0, hi0, lo2, hi2, res = slice_map(tris, axis, value)
    i0, i2 = [a for a in range(3) if a != axis]
    print("== %s  slice %s=%+.1f mm (plane %s-%s)" % (
        label, names[axis], value * 1000, names[i0], names[i2]))
    for r, row in enumerate(rows):
        v2 = hi2 - (r + 0.5) * res
        print("%+7.1f |%s|" % (v2 * 1000, "".join("#" if c else " " for c in row)))
    # axis0 ruler
    ruler = ""
    for j in range(len(rows[0])):
        v0 = lo0 + (j + 0.5) * res
        if j % 10 == 0:
            ruler += "+"
        else:
            ruler += "-"
    print("        " + ruler)
    print("        %+.1f mm ... %+.1f mm along %s" % (lo0 * 1000, hi0 * 1000, names[i0]))


def voxelize(tris, res=RES):
    """Set of (xc, yc, zc) cell indices inside the solid (x-ray fill, like mesh_inertia)."""
    pad = 2 * res
    lo = [min(p[k] for t in tris for p in t) - pad for k in range(3)]
    hi = [max(p[k] for t in tris for p in t) + pad for k in range(3)]
    n = [int((hi[k] - lo[k]) / res) + 1 for k in range(3)]
    buckets = {}
    for t in tris:
        y0 = int((min(p[1] for p in t) - lo[1]) / res)
        y1 = int((max(p[1] for p in t) - lo[1]) / res)
        z0 = int((min(p[2] for p in t) - lo[2]) / res)
        z1 = int((max(p[2] for p in t) - lo[2]) / res)
        for yc in range(max(0, y0), min(n[1] - 1, y1) + 1):
            for zc in range(max(0, z0), min(n[2] - 1, z1) + 1):
                buckets.setdefault((yc, zc), []).append(t)
    cells = set()
    for yc in range(n[1]):
        y = lo[1] + (yc + 0.5) * res
        for zc in range(n[2]):
            z = lo[2] + (zc + 0.5) * res
            ts = []
            for (a, b, c) in buckets.get((yc, zc), []):
                e1 = tuple(b[k] - a[k] for k in range(3))
                e2 = tuple(c[k] - a[k] for k in range(3))
                det = e1[2] * e2[1] - e1[1] * e2[2]
                if abs(det) < 1e-12:
                    continue
                inv = 1.0 / det
                T = (-a[0], y - a[1], z - a[2])
                P = (0.0, -e2[2], e2[1])
                u = (T[0] * P[0] + T[1] * P[1] + T[2] * P[2]) * inv
                if u < 0.0 or u > 1.0:
                    continue
                Q = (T[1] * e1[2] - T[2] * e1[1],
                     T[2] * e1[0] - T[0] * e1[2],
                     T[0] * e1[1] - T[1] * e1[0])
                vv = Q[0] * inv
                if vv < 0.0 or u + vv > 1.0:
                    continue
                ts.append((e2[0] * Q[0] + e2[1] * Q[1] + e2[2] * Q[2]) * inv)
            ts.sort()
            k = 0
            while k + 1 < len(ts):
                if ts[k + 1] - ts[k] > 0.5 * res:
                    x0 = max(lo[0], ts[k]); x1 = min(hi[0], ts[k + 1])
                    i0 = max(0, int((x0 - lo[0]) / res))
                    i1 = min(n[0] - 1, int((x1 - lo[0]) / res))
                    for xc in range(i0, i1 + 1):
                        cells.add((xc, yc, zc))
                k += 2
    return cells, lo, n, res


def voxel_overlap(name_a, frame_a, name_b, frame_b, res=0.001):
    """Overlap (in mm^3) between two links' meshes, each expressed in frame_a."""
    ta = tris_in(name_a, frame_a)
    tb = tris_in(name_b, frame_a)
    cells_a, lo_a, n_a, res = voxelize(ta, res)
    # voxelize b on the same grid: transform b's cells into a's grid
    cells_b, lo_b, n_b, _ = voxelize(tb, res)
    # snap b cells into a's indexing (same physical grid only if lo/n match)
    common = 0
    # instead: rebuild b's cells in a's grid by cell-centre coordinates
    bc = set()
    for (xc, yc, zc) in cells_b:
        x = lo_b[0] + (xc + 0.5) * res
        y = lo_b[1] + (yc + 0.5) * res
        z = lo_b[2] + (zc + 0.5) * res
        bc.add((int(round((x - lo_a[0]) / res - 0.5)),
                int(round((y - lo_a[1]) / res - 0.5)),
                int(round((z - lo_a[2]) / res - 0.5))))
    common = len(cells_a & bc)
    return common * res ** 3 * 1e6  # mm^3


if __name__ == "__main__":
    # base cross-section in its own link frame, mid-length
    print_slice(tris_in("401200xr__1_", "401200xr__1_"), 1, 0.0, label="base1")
    # carriage cross-section in base1 frame, mid-length
    print_slice(tris_in("401xr___carriage__401xr___carriage_1", "401200xr__1_"),
                1, -0.033, label="carriage(bottom) in base1 frame")
    # same slice but at the two rails' z to see the saddle: slice along y at rail height
    print_slice(tris_in("401xr___carriage__401xr___carriage_1", "401200xr__1_"),
                2, 0.005, label="carriage(bottom) horizontal cut at z=5mm in base1 frame")
    print()
    for (a, fa, b, fb, tag) in [
        ("401200xr__1_", "401200xr__1_", "401xr___carriage__401xr___carriage_1", "401200xr__1_", "base1 vs y-carriage"),
        ("401200xr__1_", "401200xr__1_", "plate", "401200xr__1_", "base1 vs plate"),
        ("401200xr__1_", "401200xr__1_", "401xr___encoder__401xr___encoder_1", "401200xr__1_", "base1 vs y-encoder"),
        ("401200xr__1_", "401200xr__1_", "401xr___carriage_end_caps__401xr___carriage_end_caps_1_1", "401200xr__1_", "base1 vs y-endcap"),
    ]:
        print("%s overlap: %.1f cm^3" % (tag, voxel_overlap(a, fa, b, fb)))
