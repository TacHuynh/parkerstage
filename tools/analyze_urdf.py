#!/usr/bin/env python3
"""Analyze parkerstage.urdf: kinematic tree, world poses at zero config,
and per-link world bounding boxes (visual mesh extents incl. origins)."""
import math
import os
import struct
import sys
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
URDF = os.path.join(ROOT, "urdf", "parkerstage.urdf")
MESH_DIR = os.path.join(ROOT, "meshes")


def rx(a):
    c, s = math.cos(a), math.sin(a)
    return [[1, 0, 0], [0, c, -s], [0, s, c]]


def ry(a):
    c, s = math.cos(a), math.sin(a)
    return [[c, 0, s], [0, 1, 0], [-s, 0, c]]


def rz(a):
    c, s = math.cos(a), math.sin(a)
    return [[c, -s, 0], [s, c, 0], [0, 0, 1]]


def mmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def rpy2m(r, p, y):
    return mmul(rz(y), mmul(ry(p), rx(r)))


def t4(R, t):
    return [
        [R[0][0], R[0][1], R[0][2], t[0]],
        [R[1][0], R[1][1], R[1][2], t[1]],
        [R[2][0], R[2][1], R[2][2], t[2]],
        [0, 0, 0, 1],
    ]


def mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def vmul(M, v):
    return [sum(M[i][j] * v[j] for j in range(3)) + M[i][3] for i in range(3)]


def rpy_of(M):
    sy = math.sqrt(M[0][0] ** 2 + M[1][0] ** 2)
    if sy > 1e-9:
        r = math.atan2(M[2][1], M[2][2])
        p = math.atan2(-M[2][0], sy)
        y = math.atan2(M[1][0], M[0][0])
    else:
        r = math.atan2(-M[1][2], M[1][1])
        p = math.atan2(-M[2][0], sy)
        y = 0.0
    return r, p, y


def stl_bbox(path):
    with open(path, "rb") as fh:
        data = fh.read()
    n = struct.unpack("<I", data[80:84])[0]
    mn = [1e9] * 3
    mx = [-1e9] * 3
    for i in range(n):
        off = 84 + i * 50
        for j in range(3):
            vert = struct.unpack("<3f", data[off + 12 + j * 12: off + 24 + j * 12])
            for k in range(3):
                mn[k] = min(mn[k], vert[k])
                mx[k] = max(mx[k], vert[k])
    return mn, mx


def main():
    tree = ET.parse(URDF)
    root = tree.getroot()

    bboxes = {}
    for fn in sorted(os.listdir(MESH_DIR)):
        if fn.endswith(".stl"):
            bboxes[fn] = stl_bbox(os.path.join(MESH_DIR, fn))

    # kinematic tree (parent -> children)
    children = {}
    jtype = {}
    for j in root.findall("joint"):
        p = j.find("parent").get("link")
        c = j.find("child").get("link")
        children.setdefault(p, []).append(c)
        jtype[(p, c)] = j.get("type")

    # world poses at zero config
    poses = {"root": t4([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0])}
    order = ["root"]
    while order:
        p = order.pop(0)
        for c in children.get(p, []):
            j = next(jj for jj in root.findall("joint") if jj.find("parent").get("link") == p and jj.find("child").get("link") == c)
            o = j.find("origin")
            xyz = list(map(float, o.get("xyz").split())) if o is not None else [0, 0, 0]
            rpy = list(map(float, o.get("rpy").split())) if o is not None else [0, 0, 0]
            poses[c] = mul(poses[p], t4(rpy2m(*rpy), xyz))
            order.append(c)

    def world_bbox(linkname):
        l = root.find('link[@name="%s"]' % linkname)
        if l is None:
            return None
        mn = [1e9] * 3
        mx = [-1e9] * 3
        for v in l.findall("visual"):
            o = v.find("origin")
            oxyz = list(map(float, o.get("xyz").split())) if o is not None else [0, 0, 0]
            orpy = list(map(float, o.get("rpy").split())) if o is not None else [0, 0, 0]
            msh = v.find("geometry/mesh")
            if msh is None:
                continue
            fn = msh.get("filename").split("/")[-1]
            bb = bboxes.get(fn)
            if bb is None:
                continue
            bbmin, bbmax = bb
            Rv = rpy2m(*orpy)
            M = mul(poses[linkname], t4(Rv, oxyz))
            for k in range(8):
                vv = [bbmin[i] if not (k >> i) & 1 else bbmax[i] for i in range(3)]
                w = vmul(M, vv)
                for i in range(3):
                    mn[i] = min(mn[i], w[i])
                    mx[i] = max(mx[i], w[i])
        return mn, mx

    print("== kinematic tree ==")
    def print_tree(p, depth):
        for c in children.get(p, []):
            print("  " * depth + "%s  [%s]" % (c, jtype[(p, c)]))
            print_tree(c, depth + 1)
    print_tree("root", 0)

    print("\n== world poses & bboxes (zero config, meters) ==")
    for name in sorted(poses):
        M = poses[name]
        t = [M[i][3] for i in range(3)]
        bb = world_bbox(name)
        if bb is None:
            print("%-55s t=(%8.4f %8.4f %8.4f)  <no visual>" % (name, t[0], t[1], t[2]))
        else:
            mn, mx = bb
            print(
                "%-55s t=(%8.4f %8.4f %8.4f)  bbox x[%7.4f %7.4f] y[%7.4f %7.4f] z[%7.4f %7.4f]"
                % (name, t[0], t[1], t[2], mn[0], mx[0], mn[1], mx[1], mn[2], mx[2])
            )


if __name__ == "__main__":
    sys.exit(main())
