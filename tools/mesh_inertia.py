#!/usr/bin/env python3
"""Volume, center of mass and inertia tensor of a binary STL mesh (voxel fill).

The Onshape STL exports are nested-shell surfaces (inner + outer), so the
classic signed-tetrahedron decomposition cancels wrongly. Instead we voxelize:
triangles are bucketed by their y-z footprint, each (y, z) column is ray-cast
with an even-odd parity fill along x, and volume / COM / inertia are summed
over the filled voxels. Robust for shell and slightly non-manifold meshes;
accurate to ~1 voxel (~1 mm at the default resolution).

Works for binary STL only (all meshes in this package are binary).
"""
import struct


def read_triangles(path):
    with open(path, "rb") as fh:
        data = fh.read()
    n = struct.unpack("<I", data[80:84])[0]
    tris = []
    for i in range(n):
        off = 84 + i * 50
        tris.append((struct.unpack("<3f", data[off + 12: off + 24]),
                     struct.unpack("<3f", data[off + 24: off + 36]),
                     struct.unpack("<3f", data[off + 36: off + 48])))
    return tris


def mesh_inertia_voxel(path, density=2700.0, res=0.001):
    """Volume / COM / inertia via voxel fill (robust for the nested-shell STLs
    Onshape exports, where the signed-tetrahedron method cancels wrongly).

    Voxels are classified inside/outside with an even-odd ray cast per
    (y, z) column; triangles are bucketed by their y-z footprint first.
    Returns (volume_m3, com_mesh, I_about_com_mesh).
    """
    tris = read_triangles(path)
    mn = [min(t[i][d] for t in tris for i in range(3)) for d in range(3)]
    mx = [max(t[i][d] for t in tris for i in range(3)) for d in range(3)]
    pad = 2 * res
    mn = [m - pad for m in mn]
    mx = [m + pad for m in mx]
    nx = int((mx[0] - mn[0]) / res) + 1
    ny = int((mx[1] - mn[1]) / res) + 1
    nz = int((mx[2] - mn[2]) / res) + 1

    # bucket triangles by (yc, zc) cell range they may cross
    buckets = [[] for _ in range(ny * nz)]
    for ti, (a, b, c) in enumerate(tris):
        y0 = int((min(a[1], b[1], c[1]) - mn[1]) / res)
        y1 = int((max(a[1], b[1], c[1]) - mn[1]) / res)
        z0 = int((min(a[2], b[2], c[2]) - mn[2]) / res)
        z1 = int((max(a[2], b[2], c[2]) - mn[2]) / res)
        for yc in range(max(0, y0), min(ny - 1, y1) + 1):
            for zc in range(max(0, z0), min(nz - 1, z1) + 1):
                buckets[yc * nz + zc].append(ti)

    inside_count = 0
    sum_x = sum_y = sum_z = 0.0
    I_origin = [0.0] * 6  # xx, yy, zz, xy, xz, yz
    m_vox = density * res ** 3

    for yc in range(ny):
        y = mn[1] + (yc + 0.5) * res
        for zc in range(nz):
            z = mn[2] + (zc + 0.5) * res
            xs = []
            for ti in buckets[yc * nz + zc]:
                a, b, c = tris[ti]
                e1 = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
                e2 = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
                # ray O=(0,y,z), D=(1,0,0); Moller-Trumbore, D fixed
                det = e1[2] * e2[1] - e1[1] * e2[2]  # (e1 x e2)_x = n_x
                if abs(det) < 1e-12:
                    continue
                inv = 1.0 / det
                T = (-a[0], y - a[1], z - a[2])
                P = (0.0, -e2[2], e2[1])  # D x e2, D=(1,0,0)
                u = (T[0] * P[0] + T[1] * P[1] + T[2] * P[2]) * inv
                if u < 0.0 or u > 1.0:
                    continue
                Q = (T[1] * e1[2] - T[2] * e1[1],
                     T[2] * e1[0] - T[0] * e1[2],
                     T[0] * e1[1] - T[1] * e1[0])  # T x e1
                vv = Q[0] * inv  # D . Q = Q.x
                if vv < 0.0 or u + vv > 1.0:
                    continue
                t = (e2[0] * Q[0] + e2[1] * Q[1] + e2[2] * Q[2]) * inv
                xs.append(t)  # P = O + t*D with O.x = 0, D = (1,0,0) -> x = t
            xs.sort()
            k = 0
            while k + 1 < len(xs):
                x0, x1 = xs[k], xs[k + 1]
                if x1 - x0 > 0.5 * res:  # inside span
                    i0 = max(0, int((x0 - mn[0]) / res))
                    i1 = min(nx - 1, int((x1 - mn[0]) / res))
                    for xc in range(i0, i1 + 1):
                        x = mn[0] + (xc + 0.5) * res
                        inside_count += 1
                        sum_x += x
                        sum_y += y
                        sum_z += z
                        I_origin[0] += m_vox * (y * y + z * z)   # xx
                        I_origin[1] += m_vox * (x * x + z * z)   # yy
                        I_origin[2] += m_vox * (x * x + y * y)   # zz
                        I_origin[3] -= m_vox * x * y             # xy
                        I_origin[4] -= m_vox * x * z             # xz
                        I_origin[5] -= m_vox * y * z             # yz
                k += 2

    if inside_count == 0:
        raise ValueError("voxel fill found no inside voxels (open mesh?): %s" % path)
    v = inside_count * res ** 3
    com = (sum_x / inside_count, sum_y / inside_count, sum_z / inside_count)
    m = density * v
    # shift to COM: I_com = I_origin + m(com com^T - |com|^2 I)
    rr = sum(c * c for c in com)
    I = [[0.0] * 3 for _ in range(3)]
    I[0][0] = I_origin[0] + m * (com[0] ** 2 - rr)
    I[1][1] = I_origin[1] + m * (com[1] ** 2 - rr)
    I[2][2] = I_origin[2] + m * (com[2] ** 2 - rr)
    I[0][1] = I[1][0] = I_origin[3] + m * com[0] * com[1]
    I[0][2] = I[2][0] = I_origin[4] + m * com[0] * com[2]
    I[1][2] = I[2][1] = I_origin[5] + m * com[1] * com[2]
    return v, com, I


inertia_tuple = lambda I: (I[0][0], I[0][1], I[0][2], I[1][1], I[1][2], I[2][2])


if __name__ == "__main__":
    import sys
    v, com, I = mesh_inertia_voxel(sys.argv[1], density=float(sys.argv[2]) if len(sys.argv) > 2 else 2700.0)
    print("volume m^3: %.6g  (cm^3: %.3f)" % (v, v * 1e6))
    print("mass kg:    %.6g  (g: %.2f)" % (v * 2700.0, v * 2700.0 * 1000))
    print("COM (mesh): (%.6f, %.6f, %.6f)" % com)
    print("inertia about COM (mesh frame):")
    for row in I:
        print("  [% .4e  % .4e  % .4e]" % tuple(row))
