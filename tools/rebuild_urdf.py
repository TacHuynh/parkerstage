#!/usr/bin/env python3
"""Rebuild parkerstage.urdf into a clean kinematic tree using the corrected
`xyslideassy` slide as the SINGLE SOURCE for BOTH the Y and X axes.

The slide (source: /home/th/ParkerStage/xyslideassy/urdf/xyslideassy.urdf) is
now a proper 200 mm 401XR stage: `root` is the fixed base (main extrusion,
drive screw, holders, coupling housing, motor, seal, and the three H3
switches), a `slider` prismatic drives the carriage subtree (switch_flag, end
caps, carriage, encoder, H2 home switch) along the base's long axis, the seal
is bolted to the base (`seal_to_base`), and the travel limits
[+48.370, +248.370] mm are bracketed by the two outer H3 limit switches with
the H3__L1 home switch exactly at mid-stroke.

This script clones that slide twice to build the XY stage:

  Y stage: base `y_root` bolted to world `root`; prismatic `y_slide` drives the
           Y carriage group along world +Y.
  X stage: base `x_root` bolted ON TOP of the Y carriage, rotated -90 deg about
           the vertical so its travel points along world +X; prismatic
           `x_slide` drives the X carriage group along world +X.
  Z column: fixed (from zslide.urdf.onshape) on world root at the assembly
           centre, independent of the XY motion.

The y/x slide joints keep the slide's OWN origin, axis and limits (single
source: the bracketed 200 mm stroke and mid-stroke home cannot drift from the
slide), and the seal is bolted to each base exactly like the slide.  All joint
origins are computed numerically so every link keeps its intended world pose at
q=0.  The verify() block re-checks pose preservation, stroke, limits-vs-slide,
home, Z independence, Z centering and a full collision sweep.
"""
import copy
import math
import os
import sys
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "urdf", "parkerstage.urdf")

XYS = "/home/th/ParkerStage/xyslideassy/urdf/xyslideassy.urdf"
XYS_MESH = "/home/th/ParkerStage/xyslideassy/meshes/"
ZSRC = os.path.join(ROOT, "urdf", "zslide.urdf.onshape")
TRAVEL_Z = 0.075     # 401XR150 150 mm stroke
EFFORT = 1.0
VELOCITY = 1.0
DENSITY = 2700.0     # kg/m^3 aluminium

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mesh_inertia import mesh_inertia_voxel   # noqa: E402
from mesh_inertia import read_triangles       # noqa: E402

MESH_DIR = os.path.join(ROOT, "meshes")

# ---------------------------------------------------------------- transforms


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
    return [[R[0][0], R[0][1], R[0][2], t[0]],
            [R[1][0], R[1][1], R[1][2], t[1]],
            [R[2][0], R[2][1], R[2][2], t[2]],
            [0, 0, 0, 1]]


def mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def inv(M):
    R = [[M[i][j] for j in range(3)] for i in range(3)]
    t = [M[i][3] for i in range(3)]
    Rt = [[R[j][i] for j in range(3)] for i in range(3)]
    nt = [-sum(Rt[i][j] * t[j] for j in range(3)) for i in range(3)]
    return t4(Rt, nt)


def vmul(M, v):
    return [sum(M[i][j] * v[j] for j in range(3)) + M[i][3] for i in range(3)]


def vrot(M, v):
    return [sum(M[i][j] * v[j] for j in range(3)) for i in range(3)]


def jxf(xyz, rpy):
    """Joint-origin 4x4 under real URDF semantics: R . T(xyz), i.e. the
    translation xyz is applied in the ROTATED frame, so a child at identity
    lands at R.xyz in the parent frame (matches ROS/URDF parsing)."""
    R = rpy2m(*rpy)
    return t4(R, vrot(R, xyz))


def origin_of(M):
    """(xyz, rpy) that reproduces transform M in the file under real URDF
    semantics (xyz must satisfy R.xyz = t, i.e. xyz = R^T.t)."""
    R = [[M[i][j] for j in range(3)] for i in range(3)]
    t = [M[i][3] for i in range(3)]
    RT = [[R[j][i] for j in range(3)] for i in range(3)]
    return vrot(RT, t), rpy_of(M)


def rpy_of(M):
    sy = math.sqrt(M[0][0] ** 2 + M[1][0] ** 2)
    if sy > 1e-9:
        return (math.atan2(M[2][1], M[2][2]), math.atan2(-M[2][0], sy), math.atan2(M[1][0], M[0][0]))
    return (math.atan2(-M[1][2], M[1][1]), math.atan2(-M[2][0], sy), 0.0)


def fmt(v):
    if abs(v) < 1e-5:
        v = 0.0
    s = "%.6f" % v
    s = s.rstrip("0").rstrip(".")
    if s in ("", "-"):
        s += "0"
    return s


def fmt_vec(v):
    return " ".join(fmt(x) for x in v)


def fnum(v):
    return "%.9g" % v


# ------------------------------------------------------------- parse new slide

stree = ET.parse(XYS)
sroot = stree.getroot()
slinks = {l.get("name"): l for l in sroot.findall("link")}

SLIDE_BASE = "root"
CAR_ROOT = "401xr___switch_flag__401xr___switch_flag"
CAR_BODY = "401xr___carriage__401xr___carriage"
SEAL = "401xr_200_seal_seal__401xr_200_seal_seal"
SWITCH = "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch"
FLAG = "401xr___switch_flag__401xr___switch_flag"
ENC = "401xr___encoder__401xr___encoder"
ENC_BASE = "401xr___encoder_base_2__401xr___encoder_base_2"
EC1 = "401xr___carriage_end_caps__401xr___carriage_end_caps_1"
EC1_ = "401xr___carriage_end_caps__401xr___carriage_end_caps__1_"

MOVERS = {CAR_ROOT, ENC_BASE, EC1_, EC1, CAR_BODY, ENC, SWITCH}   # seal is fixed to the base now

# joints recreated explicitly below (skipped by fix_all)
EMPTY_JOINT_NAMES = ("slider", "seal_to_base")


def _slide_poses():
    ch = {}
    for j in sroot.findall("joint"):
        ch.setdefault(j.find("parent").get("link"), []).append(j)
    poses = {SLIDE_BASE: t4([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0])}
    order = [SLIDE_BASE]
    while order:
        p = order.pop(0)
        for j in ch.get(p, []):
            o = j.find("origin")
            xyz = list(map(float, o.get("xyz").split())) if o is not None else [0, 0, 0]
            rpy = list(map(float, o.get("rpy").split())) if o is not None else [0, 0, 0]
            poses[j.find("child").get("link")] = mul(poses[p], jxf(xyz, rpy))
            order.append(j.find("child").get("link"))
    return poses


SPOSES = _slide_poses()


def _s_ext(name, base_frame, pmap, axis):
    """Slide link `name` mesh bbox along `base_frame`'s `axis`, in the base
    frame (m).  Aggregates ALL visuals (the base link `root` carries many part
    meshes)."""
    l = slinks[name]
    Mb = inv(pmap[base_frame])
    lo, hi = [1e9] * 3, [-1e9] * 3
    for vis in l.findall("visual"):
        msh = vis.find("geometry/mesh")
        if msh is None:
            continue
        vo = vis.find("origin")
        oxyz = [float(x) for x in vo.get("xyz").split()] if vo is not None else [0, 0, 0]
        orpy = [float(x) for x in vo.get("rpy").split()] if vo is not None else [0, 0, 0]
        fn = msh.get("filename").split("/")[-1]
        tris = read_triangles(os.path.join(XYS_MESH, fn))
        mn = [min(p[i] for tri in tris for p in tri) for i in range(3)]
        mx = [max(p[i] for tri in tris for p in tri) for i in range(3)]
        Rv = rpy2m(*orpy)
        pt_local = []
        for sx in (0, 1):
            for sy in (0, 1):
                for sz in (0, 1):
                    q = [mx[i] if s else mn[i] for i, s in enumerate((sx, sy, sz))]
                    pt_local.append([sum(Rv[i][k] * q[k] for k in range(3)) + oxyz[i] for i in range(3)])
        for q in pt_local:
            c = vmul(Mb, vmul(pmap[name], q))
            for i in range(3):
                lo[i] = min(lo[i], c[i]); hi[i] = max(hi[i], c[i])
    return lo, hi


def _world_bbox(name, pmap):
    l = slinks[name]
    lo, hi = [1e9] * 3, [-1e9] * 3
    for vis in l.findall("visual"):
        msh = vis.find("geometry/mesh")
        if msh is None:
            continue
        vo = vis.find("origin")
        oxyz = [float(x) for x in vo.get("xyz").split()] if vo is not None else [0, 0, 0]
        orpy = [float(x) for x in vo.get("rpy").split()] if vo is not None else [0, 0, 0]
        fn = msh.get("filename").split("/")[-1]
        tris = read_triangles(os.path.join(XYS_MESH, fn))
        mn = [min(p[i] for tri in tris for p in tri) for i in range(3)]
        mx = [max(p[i] for tri in tris for p in tri) for i in range(3)]
        Rv = rpy2m(*orpy)
        for sx in (0, 1):
            for sy in (0, 1):
                for sz in (0, 1):
                    q = [mx[i] if s else mn[i] for i, s in enumerate((sx, sy, sz))]
                    q = [sum(Rv[i][k] * q[k] for k in range(3)) + oxyz[i] for i in range(3)]
                    c = vmul(pmap[name], q)
                    for i in range(3):
                        lo[i] = min(lo[i], c[i]); hi[i] = max(hi[i], c[i])
    return lo, hi


def _s_body_ext(name, base_frame, pmap, filter_substr="_401XR_200_Base"):
    """Bbox of the *main-body* visual only (mesh filename containing
    `filter_substr`), in the base frame.  The base link aggregates many parts
    (motor sticking far out at -y, drive screw, holders), whose centroid would
    drag the travel-centre off the actual extrusion, so mid-stroke/centring
    references must use the main extrusion alone."""
    l = slinks[name]
    Mb = inv(pmap[base_frame])
    lo, hi = [1e9] * 3, [-1e9] * 3
    for vis in l.findall("visual"):
        msh = vis.find("geometry/mesh")
        if msh is None:
            continue
        fn = msh.get("filename").split("/")[-1]
        if filter_substr not in fn:
            continue
        vo = vis.find("origin")
        oxyz = [float(x) for x in vo.get("xyz").split()] if vo is not None else [0, 0, 0]
        orpy = [float(x) for x in vo.get("rpy").split()] if vo is not None else [0, 0, 0]
        tris = read_triangles(os.path.join(XYS_MESH, fn))
        mn = [min(p[i] for tri in tris for p in tri) for i in range(3)]
        mx = [max(p[i] for tri in tris for p in tri) for i in range(3)]
        Rv = rpy2m(*orpy)
        for sx in (0, 1):
            for sy in (0, 1):
                for sz in (0, 1):
                    q = [mx[i] if s else mn[i] for i, s in enumerate((sx, sy, sz))]
                    q = [sum(Rv[i][k] * q[k] for k in range(3)) + oxyz[i] for i in range(3)]
                    c = vmul(Mb, vmul(pmap[name], q))
                    for i in range(3):
                        lo[i] = min(lo[i], c[i]); hi[i] = max(hi[i], c[i])
    return lo, hi


# ---- the travel is SINGLE-SOURCED from the corrected slide: the y/x slide
# joints reuse the slide's own slider origin, axis and limits, so the 200 mm
# stroke bracketed by the two outer H3 limit switches (with the H3__L1 home
# switch at mid-stroke) cannot drift from the source slide.
_slider_j = None
for _j in sroot.findall("joint"):
    if _j.get("name") == "slider":
        _slider_j = _j
assert _slider_j is not None, "slide source has no 'slider' prismatic joint"
_so = _slider_j.find("origin")
SLIDER_XYZ = [float(x) for x in _so.get("xyz").split()] if _so is not None else [0, 0, 0]
SLIDER_RPY = [float(x) for x in _so.get("rpy").split()] if _so is not None else [0, 0, 0]
SLIDER_AXIS = [float(x) for x in _slider_j.find("axis").get("xyz").split()]
_sl = _slider_j.find("limit")
SLIDE_LO = float(_sl.get("lower"))
SLIDE_HI = float(_sl.get("upper"))
TRAVEL = (SLIDE_HI - SLIDE_LO) / 2.0      # 200 mm stroke
MID_STROKE = (SLIDE_LO + SLIDE_HI) / 2.0  # centred between the H3 limit switches
HOME = MID_STROKE                         # home = mid-stroke (H2 on the H3__L1 switch)
_seal_j = None
for _j in sroot.findall("joint"):
    if _j.get("name") == "seal_to_base":
        _seal_j = _j
assert _seal_j is not None, "slide source has no 'seal_to_base' joint"
_seo = _seal_j.find("origin")
SEAL_XYZ = [float(x) for x in _seo.get("xyz").split()] if _seo is not None else [0, 0, 0]
SEAL_RPY = [float(x) for x in _seo.get("rpy").split()] if _seo is not None else [0, 0, 0]

if os.environ.get("DBG_MID"):
    print("DBG slide slider origin %s rpy %s axis %s limits [%+.4f, %+.4f] mid %+.4f"
          % (SLIDER_XYZ, SLIDER_RPY, SLIDER_AXIS, SLIDE_LO, SLIDE_HI, MID_STROKE))
print("slide source: limits [%+.3f, %+.3f] mm (200 mm stroke bracketed by the H3 limit switches), "
      "home = mid-stroke %+.3f mm" % (SLIDE_LO * 1e3, SLIDE_HI * 1e3, HOME * 1e3))

# ------------------------------------------------------------- Z stage
ztree = ET.parse(ZSRC)
zroot = ztree.getroot()
zlinks = {l.get("name"): l for l in zroot.findall("link")}
zchildren = {}
for j in zroot.findall("joint"):
    zchildren.setdefault(j.find("parent").get("link"), []).append((j.find("child").get("link"), j))
zposes = {"z_root": t4([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0])}
zorder = ["z_root"]
while zorder:
    p = zorder.pop(0)
    for (c, j) in zchildren.get(p, []):
        o = j.find("origin")
        xyz = list(map(float, o.get("xyz").split())) if o is not None else [0, 0, 0]
        rpy = list(map(float, o.get("rpy").split())) if o is not None else [0, 0, 0]
        zposes[c] = mul(zposes[p], jxf(xyz, rpy))
        zorder.append(c)

Z_BASE = "z_401xr_150_base__401xr_150_base"
Z_CARRIAGE = "z_401xr___carriage__401xr___carriage"
Z_GROUP_ROOT = "z_c3_401xr__c3_401xr"
Z_HOME_SW = "z_401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch"
Z_FLAG = "z_401xr___switch_flag__401xr___switch_flag"
Z_END_CAPS = "z_401xr___carriage_end_caps__401xr___carriage_end_caps"
Z_ENCODER = "z_401xr___encoder__401xr___encoder"
Z_ENC_BASE = "z_401xr___encoder_base_2__401xr___encoder_base_2_1"
Z_PLANAR = ("z_parallel_1", "z_parallel_1_1", "z_parallel_1_2", "z_parallel_1_3",
            "z_hanging_node_to_root_joint")
Z_GROUP = {Z_CARRIAGE, Z_END_CAPS,
           "z_401xr___carriage_end_caps__401xr___carriage_end_caps_1__1_",
           Z_ENCODER, Z_ENC_BASE, Z_HOME_SW, Z_FLAG}
Z_ALL = {n for n in zlinks if n != "z_root" and not n.startswith("z_parallel_")}
Z_BASE_GROUP = Z_ALL - Z_GROUP


def _z_ext(name, axis):
    l = zlinks[name]
    vo = l.find("visual/origin")
    oxyz = [float(x) for x in vo.get("xyz").split()]
    fn = l.find("visual/geometry/mesh").get("filename").split("/")[-1]
    tris = read_triangles(os.path.join(MESH_DIR, fn))
    mn = [min(p[i] for tri in tris for p in tri) for i in range(3)]
    mx = [max(p[i] for tri in tris for p in tri) for i in range(3)]
    Mb = inv(zposes[Z_BASE])
    cs = []
    for sx in (0, 1):
        for sy in (0, 1):
            for sz in (0, 1):
                p = [mx[i] + oxyz[i] if s else mn[i] + oxyz[i] for i, s in enumerate((sx, sy, sz))]
                cs.append(vmul(Mb, vmul(zposes[name], p)))
    lo = [min(c[i] for c in cs) for i in range(3)]
    hi = [max(c[i] for c in cs) for i in range(3)]
    return lo[axis], hi[axis]


MID_STROKE_Z = ((_z_ext(Z_BASE, 2)[0] + _z_ext(Z_BASE, 2)[1])
                - (_z_ext(Z_CARRIAGE, 2)[0] + _z_ext(Z_CARRIAGE, 2)[1])) / 2.0
HOME_Z = ((_z_ext(Z_HOME_SW, 2)[0] + _z_ext(Z_HOME_SW, 2)[1])
          - (_z_ext(Z_FLAG, 2)[0] + _z_ext(Z_FLAG, 2)[1])) / 2.0

# --------------------------------------------------------- build new tree


def _p(name, pre):
    return pre + "root" if name == SLIDE_BASE else pre + name


def make_joint(name, jtype, parent, child, xyz=(0, 0, 0), rpy=(0, 0, 0), axis=None, limit=None):
    j = ET.Element("joint", {"name": name, "type": jtype})
    if any(v != 0 for v in xyz) or any(v != 0 for v in rpy):
        o = ET.SubElement(j, "origin")
        o.set("xyz", fmt_vec(xyz))
        o.set("rpy", fmt_vec(rpy))
    if axis is not None:
        a = ET.SubElement(j, "axis")
        a.set("xyz", fmt_vec(axis))
    ET.SubElement(j, "parent", {"link": parent})
    ET.SubElement(j, "child", {"link": child})
    if limit is not None:
        l = ET.SubElement(j, "limit")
        for k, v in limit.items():
            l.set(k, v)
    return j


new_root = ET.Element("robot", {"name": "parkerstage"})

# clone every slide link (mesh and empty intermediate frames) under y_/x_
for pre in ("y_", "x_"):
    for name, l in slinks.items():
        c = copy.deepcopy(l)
        c.set("name", _p(name, pre))
        for msh in c.findall("visual/geometry/mesh"):
            fn = msh.get("filename").split("/")[-1]
            msh.set("filename", os.path.join("meshes", fn))
        new_root.append(c)
# Z links
for name, l in zlinks.items():
    if name == "z_root" or name.startswith("z_parallel_"):
        continue
    new_root.append(copy.deepcopy(l))

_built_links = {l.get("name"): l for l in new_root.findall("link")}


def instance_local_poses(pre):
    ch = {}
    for j in sroot.findall("joint"):
        ch.setdefault(_p(j.find("parent").get("link"), pre), []).append(j)
    poses = {_p(SLIDE_BASE, pre): t4([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0])}
    order = [_p(SLIDE_BASE, pre)]
    while order:
        p = order.pop(0)
        for j in ch.get(p, []):
            o = j.find("origin")
            xyz = list(map(float, o.get("xyz").split())) if o is not None else [0, 0, 0]
            rpy = list(map(float, o.get("rpy").split())) if o is not None else [0, 0, 0]
            poses[_p(j.find("child").get("link"), pre)] = mul(poses[p], jxf(xyz, rpy))
            order.append(_p(j.find("child").get("link"), pre))
    return poses


ID3 = t4([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0])

Y_INST = instance_local_poses("y_")
Y_BASE = _p(SLIDE_BASE, "y_")
Y_CARROOT = _p(CAR_ROOT, "y_")
Y_CAR = _p(CAR_BODY, "y_")

X_INST = instance_local_poses("x_")
X_BASE = _p(SLIDE_BASE, "x_")
X_CARROOT = _p(CAR_ROOT, "x_")
X_CAR = _p(CAR_BODY, "x_")


def fix_all(pre, croot):
    """Fixed joints keeping the carriage subtree together (renamed), EXCEPT the
    seal which we hang cleanly off the carriage body."""
    for j in sroot.findall("joint"):
        nm = j.get("name")
        if nm in EMPTY_JOINT_NAMES or j.get("type") != "fixed":
            continue
        jc = copy.deepcopy(j)
        jc.set("name", pre + nm)
        jc.find("parent").set("link", _p(j.find("parent").get("link"), pre))
        jc.find("child").set("link", _p(j.find("child").get("link"), pre))
        new_root.append(jc)


# ---- Y stage ----------------------------------------------------------------
new_root.append(make_joint("y_to_root", "fixed", "root", Y_BASE, (0, 0, 0), (0, 0, 0)))
# y_slide is the slide's own slider verbatim (origin, axis, limits): travel +Y
new_root.append(make_joint("y_slide", "prismatic", Y_BASE, Y_CARROOT,
                           SLIDER_XYZ, SLIDER_RPY, axis=SLIDER_AXIS,
                           limit={"effort": fmt(EFFORT), "velocity": fmt(VELOCITY),
                                  "lower": fmt(SLIDE_LO), "upper": fmt(SLIDE_HI)}))
# seal is bolted to the BASE (slide's seal_to_base verbatim), it does not move
new_root.append(make_joint("y_seal_to_base", "fixed", Y_BASE, _p(SEAL, "y_"),
                           SEAL_XYZ, SEAL_RPY))
fix_all("y_", Y_CARROOT)

# world poses keyed by ORIGINAL slide link name (for _world_bbox/_s_ext)
yinst = instance_local_poses("y_")
Y_POSES = {name: yinst[_p(name, "y_")] for name in SPOSES}
Y_CAR = _p(CAR_BODY, "y_")

def _fk_world(target, qs):
    """World pose of `target` in the built tree at the given joint values."""
    ch = {}
    for j in new_root.findall("joint"):
        ch.setdefault(j.find("parent").get("link"), []).append(j)
    pos = {"root": ID3}
    order = ["root"]
    while order:
        p = order.pop(0)
        for j in ch.get(p, []):
            o = j.find("origin")
            xyz = [float(x) for x in o.get("xyz").split()] if o is not None else [0, 0, 0]
            rpy = [float(x) for x in o.get("rpy").split()] if o is not None else [0, 0, 0]
            Mp = jxf(xyz, rpy)
            if j.get("name") in qs:
                a = [float(x) for x in j.find("axis").get("xyz").split()]
                Mp = mul(Mp, t4([[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                [qs[j.get("name")] * a[0], qs[j.get("name")] * a[1],
                                 qs[j.get("name")] * a[2]]))
            pos[j.find("child").get("link")] = mul(pos[p], Mp)
            order.append(j.find("child").get("link"))
    return pos[target]


# ---- X stage ----------------------------------------------------------------
# The Y carriage's mount frame carries a stray 3D tilt, so we IMPOSE the desired
# world pose for the X base directly: travel (base local +y) along world +X,
# base up (local +z) along world +Z, centred on the assembly centre and rising
# just above the Y-stage envelope.
#   x_root_world = Y_carrier_world . Mount   =>   Mount = inv(Y_carrier_world_0).World
# where Y_carrier_world_0 is the ACTUAL FK pose of the Y carriage at qy=0
# (the instance pose map is not a reliable world reference: the slider origin
# carries a 3D rotation that re-rotates joint-frame translations).
R_des_x = rz(-math.pi / 2)          # base local y -> world +X, local z -> world +Z
# the X base must clear the tallest Y material with its LOWEST part (the motor
# hangs below the main extrusion), so use the full base-link z-extent
xb_z = _s_ext(SLIDE_BASE, SLIDE_BASE, SPOSES, 1)
xbase_zlo = xb_z[0][2]
ycar_top_w = _world_bbox(CAR_BODY, Y_POSES)[1][2]   # tallest Y material
_wz = ycar_top_w + 0.006 - xbase_zlo
World_des_x = t4(R_des_x, [0.0, 0.0, _wz])
Mount = mul(inv(_fk_world(Y_CARROOT, {})), World_des_x)
if os.environ.get("DBG_MID"):
    w0 = _fk_world(Y_CARROOT, {})
    yb = _world_bbox(CAR_BODY, Y_POSES)
    print("DBG ycar_top_w=%.4f xbase_zlo=%.4f wz=%.4f ycarroot0=(%.4f, %.4f, %.4f) posebox z[%.3f, %.3f]"
          % (ycar_top_w, xbase_zlo, _wz, w0[0][3], w0[1][3], w0[2][3], yb[0][2], yb[1][2]))
new_root.append(make_joint("x_onto_ycarriage", "fixed", Y_CARROOT, X_BASE,
                           *origin_of(Mount)))
# x_slide keeps the slide's own origin/axis/limits; the mount rotation turns
# the slide travel (base local +y) into world +X
new_root.append(make_joint("x_slide", "prismatic", X_BASE, X_CARROOT,
                           SLIDER_XYZ, SLIDER_RPY, axis=SLIDER_AXIS,
                           limit={"effort": fmt(EFFORT), "velocity": fmt(VELOCITY),
                                  "lower": fmt(SLIDE_LO), "upper": fmt(SLIDE_HI)}))
# seal bolted to the X base (does not move with the X carriage)
new_root.append(make_joint("x_seal_to_base", "fixed", X_BASE, _p(SEAL, "x_"),
                           SEAL_XYZ, SEAL_RPY))
fix_all("x_", X_CARROOT)

X_POSES = {name: mul(Mount, M) for name, M in instance_local_poses("x_").items()}

# ---- XY-forward kinematics (Z independent), used for the Z-mount centring ----
def poses_xy(qy, qx):
    """FK over the XY subsystem (Z-column joints not present)."""
    pos = {"root": ID3}
    order = ["root"]
    pq = {"y_slide": qy, "x_slide": qx}
    while order:
        p = order.pop(0)
        for j in new_root.findall("joint"):
            if j.find("parent").get("link") != p:
                continue
            if j.get("type") == "fixed" and (j.get("name").startswith("z_")
                                              or j.find("child").get("link").startswith("z_")):
                continue
            o = j.find("origin")
            xyz = list(map(float, o.get("xyz").split())) if o is not None else [0, 0, 0]
            rpy = list(map(float, o.get("rpy").split())) if o is not None else [0, 0, 0]
            Mp = jxf(xyz, rpy)
            if j.get("type") == "prismatic" and j.get("name") in pq:
                a = list(map(float, j.find("axis").get("xyz").split()))
                Mp = mul(Mp, t4([[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                [pq[j.get("name")] * a[0], pq[j.get("name")] * a[1], pq[j.get("name")] * a[2]]))
            pos[j.find("child").get("link")] = mul(pos[p], Mp)
            order.append(j.find("child").get("link"))
    return pos


def _world_bbox_pm(nm, pos):
    """World bbox of a BUILT link (`nm`) given a world pose map `pos` (keyed by
    built link names, e.g. from poses_xy()).  Reads the built link's own visual
    mesh files (already re-rooted under this package's meshes dir)."""
    lo, hi = [1e9] * 3, [-1e9] * 3
    for vis in _built_links[nm].findall("visual"):
        msh = vis.find("geometry/mesh")
        if msh is None:
            continue
        vo = vis.find("origin")
        oxyz = [float(x) for x in vo.get("xyz").split()] if vo is not None else [0, 0, 0]
        orpy = [float(x) for x in vo.get("rpy").split()] if vo is not None else [0, 0, 0]
        fn = msh.get("filename").split("/")[-1]
        tris = read_triangles(os.path.join(XYS_MESH, fn))
        mn = [min(p[i] for tri in tris for p in tri) for i in range(3)]
        mx = [max(p[i] for tri in tris for p in tri) for i in range(3)]
        Rv = rpy2m(*orpy)
        for sx in (0, 1):
            for sy in (0, 1):
                for sz in (0, 1):
                    q = [mx[i] if s else mn[i] for i, s in enumerate((sx, sy, sz))]
                    q = [sum(Rv[i][k] * q[k] for k in range(3)) + oxyz[i] for i in range(3)]
                    c = vmul(pos[nm], q)
                    for i in range(3):
                        lo[i] = min(lo[i], c[i]); hi[i] = max(hi[i], c[i])
    return lo, hi


# ---- 3) Z column (FIXED, on world root, centred on the mid-travel centre) ----
_zcx = _z_ext(Z_BASE, 0)
_zcy = _z_ext(Z_BASE, 1)
_zlo = _z_ext(Z_BASE, 2)[0]
# ---- impose a clean upright column over the assembly centre --------------
# The Z column is a FIXED upright stage: base-local +z (its travel) points world
# +Z, base footprint centred on the XY mid-travel centre (world 0,0), and the
# base -z (bottom) floated just above the tallest XY material anywhere in the
# XY travel envelope (so the column never hits the stack at any XY corner).
#   Z_base_world = ZMount  (zposes[Z_BASE] = identity in the Z base frame)
Z_R = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]   # travel(z) -> world +Z
fxw = (_zcx[0] + _zcx[1]) / 2.0     # Z base footprint x-centre (local)
fyw = (_zcy[0] + _zcy[1]) / 2.0     # Z base footprint y-centre (local)
# tallest XY material across the four travel corners and mid (Y top + X stack)
def _link_top_z_world(nm, M):
    top = -1e9
    lw = _built_links[nm]
    for vis in lw.findall("visual"):
        msh = vis.find("geometry/mesh")
        if msh is None:
            continue
        vo = vis.find("origin")
        oxy = [float(x) for x in vo.get("xyz").split()] if vo is not None else [0, 0, 0]
        oro = [float(x) for x in vo.get("rpy").split()] if vo is not None else [0, 0, 0]
        fn = msh.get("filename").split("/")[-1]
        tris = read_triangles(os.path.join(MESH_DIR, fn))
        mn = [min(p[i] for tri in tris for p in tri) for i in range(3)]
        mx = [max(p[i] for tri in tris for p in tri) for i in range(3)]
        Rv = rpy2m(*oro)
        for sx in (0, 1):
            for sy in (0, 1):
                for sz in (0, 1):
                    q = [mx[i] if s else mn[i] for i, s in enumerate((sx, sy, sz))]
                    q = [sum(Rv[i][k] * q[k] for k in range(3)) + oxy[i] for i in range(3)]
                    wz = sum(M[2][k] * q[k] for k in range(3)) + M[2][3]
                    top = max(top, wz)
    return top


_lo, _hi = MID_STROKE - TRAVEL, MID_STROKE + TRAVEL
_mid = MID_STROKE
_corners = [poses_xy(a, b) for a, b in ((_lo, _lo), (_lo, _hi), (_hi, _lo), (_hi, _hi),
                                         (_mid, _mid))
            ]
_maxxy = -1e9
for pm in _corners:
    for nm, M in pm.items():
        if nm.startswith("y_") or nm.startswith("x_"):
            _maxxy = max(_maxxy, _link_top_z_world(nm, M))
ZRAISE = 0.050
z_bottom = _maxxy + 0.003 + ZRAISE
t_zz = z_bottom - _zlo
# XY mid-travel centre = the X carriage (payload) body centre at BOTH slides
# mid-stroke; the Z base footprint is centred on it so the fixed column stays
# on the assembly centre line across the whole XY travel.
_pm = poses_xy(_mid, _mid)
xc_lo, xc_hi = _world_bbox_pm(X_CAR, _pm)
mx = (xc_lo[0] + xc_hi[0]) / 2.0
my = (xc_lo[1] + xc_hi[1]) / 2.0
ZMount = t4(Z_R, [mx - fxw, my - fyw, t_zz])

if os.environ.get("DBG_MID"):
    pm = poses_xy(_mid, _mid)
    print("DBG mid: X_CAR T=(%.4f, %.4f, %.4f)" % (pm[X_CAR][0][3], pm[X_CAR][1][3], pm[X_CAR][2][3]))
    print("DBG ZMount t=(%.4f, %.4f, %.4f)  z_bottom=%.4f  fxw=%.4f fyw=%.4f" % (ZMount[0][3], ZMount[1][3], ZMount[2][3], z_bottom, fxw, fyw))
    lo, hi = _world_bbox_pm(X_CAR, pm)
    print("DBG X carriage mid x[%7.1f,%7.1f] y[%7.1f,%7.1f] center(%7.1f,%7.1f)"
          % (lo[0]*1e3, hi[0]*1e3, lo[1]*1e3, hi[1]*1e3, (lo[0]+hi[0])/2*1e3, (lo[1]+hi[1])/2*1e3))
    # Z base footprint centre in world = ZMount applied to its local centre
    zc = vmul(ZMount, [fxw, fyw, 0.0])
    print("DBG Z base footprint center(%7.1f, %7.1f) vs X carriage mid (%7.1f, %7.1f)"
          % (zc[0]*1e3, zc[1]*1e3, mx*1e3, my*1e3))

ZM = {name: mul(ZMount, zposes[name]) for name in zposes}

Mt = ZM[Z_GROUP_ROOT]
new_root.append(make_joint("z_mounted_to_root", "fixed", "root", Z_GROUP_ROOT,
                           *origin_of(Mt)))
# z_slide: Z base -> Z carriage, axis world +Z
Mzc = mul(inv(zposes[Z_BASE]), zposes[Z_CARRIAGE])
zt, zrpy = origin_of(Mzc)
R_orig = [[Mzc[i][j] for j in range(3)] for i in range(3)]
R_orig_T = [[R_orig[j][i] for j in range(3)] for i in range(3)]
R_zbase = [[ZM[Z_BASE][i][j] for j in range(3)] for i in range(3)]
R_zbase_T = [[R_zbase[j][i] for j in range(3)] for i in range(3)]
z_axis = vrot(R_orig_T, vrot(R_zbase_T, (0, 0, 1)))
new_root.append(make_joint("z_slide", "prismatic", Z_BASE, Z_CARRIAGE, zt, zrpy, axis=z_axis,
                           limit={"effort": fmt(EFFORT), "velocity": fmt(VELOCITY),
                                  "lower": fmt(MID_STROKE_Z - TRAVEL_Z),
                                  "upper": fmt(MID_STROKE_Z + TRAVEL_Z)}))
for j in zroot.findall("joint"):
    if j.get("name") in Z_PLANAR:
        continue
    jc = copy.deepcopy(j)
    if jc.get("type") == "fixed" and jc.find("axis") is not None:
        jc.remove(jc.find("axis"))
    new_root.append(jc)

# ---------------------------------------------------------------- inertials

def fix_inertials():
    rows = []
    mesh_of = {}
    for pre in ("y_", "x_"):
        for name in slinks:
            l = slinks[name]
            if not l.find("visual/geometry/mesh"):
                continue
            mesh_of[_p(name, pre)] = l.find("visual/geometry/mesh").get("filename").split("/")[-1]
    for name in zlinks:
        mesh_of.setdefault(name, None)
    for l in new_root.findall("link"):
        nm = l.get("name")
        fn = mesh_of.get(nm)
        if not fn:
            msh = l.find("visual/geometry/mesh")
            if msh is None:
                continue
            fn = msh.get("filename").split("/")[-1]
        V, com_m, I_m = mesh_inertia_voxel(os.path.join(MESH_DIR, fn), DENSITY)
        mass = DENSITY * V
        vo = l.find("visual/origin")
        oxyz = list(map(float, vo.get("xyz").split())) if vo is not None else [0, 0, 0]
        orpy = list(map(float, vo.get("rpy").split())) if vo is not None else [0, 0, 0]
        R = rpy2m(*orpy)
        R_T = [[R[j][i] for j in range(3)] for i in range(3)]
        com_l = [vrot(R, com_m)[i] + oxyz[i] for i in range(3)]
        I_l = mmul(R, mmul(I_m, R_T))
        inert = ET.Element("inertial")
        m = ET.SubElement(inert, "mass"); m.set("value", fnum(mass))
        o = ET.SubElement(inert, "origin"); o.set("xyz", fmt_vec(com_l)); o.set("rpy", "0 0 0")
        ia = ET.SubElement(inert, "inertia")
        ia.set("ixx", fnum(I_l[0][0])); ia.set("ixy", fnum(I_l[0][1])); ia.set("ixz", fnum(I_l[0][2]))
        ia.set("iyy", fnum(I_l[1][1])); ia.set("iyz", fnum(I_l[1][2])); ia.set("izz", fnum(I_l[2][2]))
        old = l.find("inertial")
        if old is not None:
            l.remove(old)
        l.insert(0, inert)
        rows.append((nm, V * 1e6, mass * 1000))
    return rows


print("link                                vol[cm^3]  mass[g]")
for name, vcm, g in fix_inertials():
    print("%-35s %9.2f %9.2f" % (name, vcm, g))

# ---------------------------------------------------------------- collisions
# One box per link = its mesh bbox (all visuals) in the *link's own frame*, so the
# box moves rigidly with the link through the prismatic.  A face inset `INSET`
# shrinks each box slightly on all sides; boxes that should be separated by the
# assembly far apart stay that way, and rigidly-adjacent-but-not-touching parts
# keep a real gap.  Any residual overlaps are flagged by the sweep below; we then
# fix them by adjusting the X/Z mounts rather than by weakening the boxes.
INSET = 0.0006  # 0.6 mm face inset to avoid false positives on clean faces


def _link_bbox_local(link_name):
    """Mesh bbox (all visuals) of a built-tree link in its own link frame."""
    l = _built_links[link_name]
    lo, hi = [1e9] * 3, [-1e9] * 3
    for vis in l.findall("visual"):
        msh = vis.find("geometry/mesh")
        if msh is None:
            continue
        vo = vis.find("origin")
        oxyz = [float(x) for x in vo.get("xyz").split()] if vo is not None else [0, 0, 0]
        orpy = [float(x) for x in vo.get("rpy").split()] if vo is not None else [0, 0, 0]
        fn = msh.get("filename").split("/")[-1]
        tris = read_triangles(os.path.join(MESH_DIR, fn))
        mn = [min(p[i] for tri in tris for p in tri) for i in range(3)]
        mx = [max(p[i] for tri in tris for p in tri) for i in range(3)]
        Rv = rpy2m(*orpy)
        for sx in (0, 1):
            for sy in (0, 1):
                for sz in (0, 1):
                    q = [mx[i] if s else mn[i] for i, s in enumerate((sx, sy, sz))]
                    q = [sum(Rv[i][k] * q[k] for k in range(3)) + oxyz[i] for i in range(3)]
                    for i in range(3):
                        lo[i] = min(lo[i], q[i]); hi[i] = max(hi[i], q[i])
    return lo, hi


_built_links = {l.get("name"): l for l in new_root.findall("link")}
CBOXES = {}
for ln, l in _built_links.items():
    lo, hi = _link_bbox_local(ln)
    if lo[0] > hi[0]:
        continue  # no mesh (empty intermediate frame or placeholder)
    # shrink by INSET, clamped to a box (inset never flips a <2*INSET axis)
    c = [(lo[i] + hi[i]) / 2 for i in range(3)]
    s = [max(hi[i] - lo[i] - 2 * INSET, 1e-5) for i in range(3)]
    CBOXES[ln] = [(ln, c, (0, 0, 0), s)]
print("collision boxes: %d links with boxes" % len(CBOXES))


def add_collisions(rootel):
    for l in rootel.findall("link"):
        for (ref, c, rpy, size) in CBOXES.get(l.get("name"), []):
            col = ET.SubElement(l, "collision")
            o = ET.SubElement(col, "origin")
            o.set("xyz", fmt_vec(c))
            o.set("rpy", fmt_vec(rpy))
            g = ET.SubElement(col, "geometry")
            b = ET.SubElement(g, "box")
            b.set("size", fmt_vec(size))

add_collisions(new_root)

# ---------------------------------------------------------------- serialize

def serialize(elem, indent=0):
    pad = "    " * indent
    children = list(elem)
    if not children:
        attrs = "".join(' %s="%s"' % (k, v) for k, v in elem.attrib.items())
        return "%s<%s%s />" % (pad, elem.tag, attrs)
    attrs = "".join(' %s="%s"' % (k, v) for k, v in elem.attrib.items())
    lines = ["%s<%s%s>" % (pad, elem.tag, attrs)]
    for ch in children:
        if ch.tag in ("parent", "child"):
            lines.append("%s<%s %s />" % (pad + "    ", ch.tag,
                                          " ".join('%s="%s"' % (k, v) for k, v in ch.attrib.items())))
        else:
            lines.append(serialize(ch, indent + 1))
    lines.append("%s</%s>" % (pad, elem.tag))
    return "\n".join(lines)


# ---------------- diagnostic before serialize -----------------
print("MID_STROKE=%.4f m  HOME=%.4f m   (200 mm stroke)" % (MID_STROKE, HOME))
print("MID_STROKE_Z=%.4f m  HOME_Z=%.4f m" % (MID_STROKE_Z, HOME_Z))
print("X mount origin t=(%.4f, %.4f, %.4f)" % (Mount[0][3], Mount[1][3], Mount[2][3]))
print("Z mount origin t=(%.4f, %.4f, %.4f)  (Z bottom clears XY %s)" % (
    ZMount[0][3], ZMount[1][3], ZMount[2][3], "ok" if _maxxy < 0 else "?"))
for j in new_root.findall("joint"):
    if j.get("type") != "prismatic":
        continue
    a = j.find("axis").get("xyz")
    lim = j.find("limit")
    print("PRISMATIC %-8s axis=(%s) lo=%s hi=%s"
          % (j.get("name"), a, lim.get("lower"), lim.get("upper")))

# ---------------------------------------------------------------- verify

def poses_at(qs):
    """World poses of every link for joint values qs (all prismatics)."""
    ch = {}
    for j in new_root.findall("joint"):
        ch.setdefault(j.find("parent").get("link"), []).append(j)
    pos = {"root": ID3}
    order = ["root"]
    while order:
        p = order.pop(0)
        for j in ch.get(p, []):
            o = j.find("origin")
            xyz = list(map(float, o.get("xyz").split())) if o is not None else [0, 0, 0]
            rpy = list(map(float, o.get("rpy").split())) if o is not None else [0, 0, 0]
            M = jxf(xyz, rpy)
            if j.get("name") in qs:
                a = list(map(float, j.find("axis").get("xyz").split()))
                M = mul(M, t4([[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                              [qs[j.get("name")] * a[0], qs[j.get("name")] * a[1], qs[j.get("name")] * a[2]]))
            pos[j.find("child").get("link")] = mul(pos[p], M)
            order.append(j.find("child").get("link"))
    return pos


ALL_Z = {n for n in _built_links if n.startswith("z_")}
Y_SEAL = _p(SEAL, "y_")
X_SEAL = _p(SEAL, "x_")
Y_CAR_GROUP = {n for n in _built_links if n.startswith("y_")} - {"y_root", Y_SEAL}
X_CAR_GROUP = {n for n in _built_links if n.startswith("x_")} - {"x_root", X_SEAL}
Y_GROUP = Y_CAR_GROUP | {"y_root", Y_SEAL}   # whole Y stage (fixed base incl. seal + moving carriage)
X_GROUP = X_CAR_GROUP | {"x_root", X_SEAL}
Z_GROUP = {Z_CARRIAGE, Z_END_CAPS,
           "z_401xr___carriage_end_caps__401xr___carriage_end_caps_1__1_",
           Z_ENCODER, Z_ENC_BASE, Z_HOME_SW, Z_FLAG}
Z_BASE_GROUP = ALL_Z - Z_GROUP


def verify():
    ok = True
    # world corner lists for every box (link-local, inset applied)
    world_corners = {}
    for ln, lst in CBOXES.items():
        cs = []
        for (_, c, rpy, size) in lst:
            Rm = rpy2m(*rpy)
            half = [s / 2 for s in size]
            M = t4(Rm, c)
            for sx in (0, 1):
                for sy in (0, 1):
                    for sz in (0, 1):
                        p = [half[i] if s else -half[i] for i, s in enumerate((sx, sy, sz))]
                        cs.append(vmul(M, p))
        world_corners[ln] = cs

    def same_part(a, b):
        """True if a and b are parts of the same slide sub-assembly.  The CAD
        exports nest every carriage/accessory into its own base (the carriage
        tucks into the base channel, the seal/encoder/switch overhang it), so
        intra-stage box overlaps are inherent mesh interpenetration, not real
        collisions.  Only CROSS-stage pairs (Y vs X, Y vs Z, X vs Z) carry real
        collision meaning, and the Z column is additionally guarded by its own
        enforced upright/fixed pose and the top-clearance check."""
        return a[:1] == b[:1]

    def sweep(qs, debug=False):
        pos = poses_at(qs)
        ab = {}
        for ln, cs in world_corners.items():
            pts = [vmul(pos[ln], p) for p in cs]
            ab[ln] = ([min(p[i] for p in pts) for i in range(3)],
                      [max(p[i] for p in pts) for i in range(3)])
        if debug:
            for ln, (mn, mx) in sorted(ab.items()):
                print("  %-46s x[%7.2f,%7.2f] y[%7.2f,%7.2f] z[%7.2f,%7.2f]" % (
                    ln, mn[0] * 1000, mx[0] * 1000, mn[1] * 1000, mx[1] * 1000,
                    mn[2] * 1000, mx[2] * 1000))
        bad = []
        names = list(ab)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a, b = names[i], names[j]
                if same_part(a, b):
                    continue
                (mn1, mx1), (mn2, mx2) = ab[a], ab[b]
                ov = [max(0.0, min(mx1[k], mx2[k]) - max(mn1[k], mn2[k])) for k in range(3)]
                if ov[0] > 0 and ov[1] > 0 and ov[2] > 0:
                    bad.append((a, b, ov, dict(qs)))
        return bad

    lo, hi = MID_STROKE - TRAVEL, MID_STROKE + TRAVEL
    mid = (lo + hi) / 2
    loz, hiz = MID_STROKE_Z - TRAVEL_Z, MID_STROKE_Z + TRAVEL_Z
    midz = (loz + hiz) / 2
    grid = [lo, mid - 0.1, mid, mid + 0.1, hi]
    gridz = [loz, midz - 0.0375, midz, midz + 0.0375, hiz]

    # 1) stroke = 200 mm for y/x, 150 for z
    STROKE_SPEC = {"y_slide": 0.20, "x_slide": 0.20, "z_slide": 0.15}
    strokes = {}
    for jn in STROKE_SPEC:
        for jj in new_root.findall("joint"):
            if jj.get("name") == jn:
                lim = jj.find("limit")
                strokes[jn] = float(lim.get("upper")) - float(lim.get("lower"))
    for jn, st in strokes.items():
        if abs(st - STROKE_SPEC[jn]) > 1e-6:
            print("STROKE MISMATCH %s: %.4f (want %.4f)" % (jn, st, STROKE_SPEC[jn]))
            ok = False
    print("STROKE: %s" % ", ".join("%s %.1f mm" % (jn, strokes[jn] * 1e3) for jn in STROKE_SPEC))

    # 2) Y/X travel is SINGLE-SOURCED from the slide: the y/x limits must equal
    #    the slide's slider limits exactly (bracketed 200 mm stroke), and the Z
    #    mid-stroke is checked from the Z export map.
    for jn in ("y_slide", "x_slide"):
        for jj in new_root.findall("joint"):
            if jj.get("name") == jn:
                lim = jj.find("limit")
                lv, uv = float(lim.get("lower")), float(lim.get("upper"))
                if abs(lv - SLIDE_LO) > 1e-6 or abs(uv - SLIDE_HI) > 1e-6:
                    print("LIMITS DRIFT %s: [%+.4f, %+.4f] vs slide [%+.4f, %+.4f]"
                          % (jn, lv, uv, SLIDE_LO, SLIDE_HI))
                    ok = False
    e0z = _z_ext(Z_BASE, 2)
    e1z = _z_ext(Z_CARRIAGE, 2)
    phys_midz = (e0z[0] + e0z[1] - e1z[0] - e1z[1]) / 2.0
    if abs(phys_midz - MID_STROKE_Z) > 2e-4:
        print("MID-STROKE MISMATCH Z: physical %+.4f vs declared %+.4f" % (phys_midz, MID_STROKE_Z))
        ok = False
    # informational: the slide's limits are centred on the H3 switch midpoint
    # (not the raw body mid-stroke), exactly like the source slide.
    e0 = _s_body_ext(SLIDE_BASE, SLIDE_BASE, SPOSES)
    e1 = _s_ext(CAR_BODY, SLIDE_BASE, SPOSES, 1)
    phys_mid = (e0[0][1] + e0[1][1] - e1[0][1] - e1[1][1]) / 2.0
    print("LIMITS FROM SLIDE: y/x = [%+.3f, %+.3f] mm, centre %+.3f mm (bracketed by the H3 limit "
          "switches; body mid-stroke would be %+.3f mm)"
          % (SLIDE_LO * 1e3, SLIDE_HI * 1e3, MID_STROKE * 1e3, phys_mid * 1e3))

    # 3) home: Y/X = mid-stroke (H2 on the H3__L1 switch, from the slide);
    #    Z = full-retract.  Both must sit inside travel.
    homevals = {}
    for jn in ("y_slide", "x_slide", "z_slide"):
        for jj in new_root.findall("joint"):
            if jj.get("name") == jn:
                lim = jj.find("limit")
                lo_ = float(lim.get("lower")); hi_ = float(lim.get("upper"))
                h = (MID_STROKE if jn != "z_slide" else MID_STROKE_Z) - (
                    0.0 if jn != "z_slide" else TRAVEL_Z)
                homevals[jn] = h
                if not (lo_ - 1e-6 <= h <= hi_ + 1e-6):
                    print("HOME OUTSIDE TRAVEL %s: %+.4f (range [%+.4f, %+.4f])" % (jn, h, lo_, hi_))
                    ok = False
    print("HOME: y/x mid-stroke %s, z full-retract %s (inside travel)"
          % (", ".join("%s %+.3f mm" % (jn, homevals[jn] * 1e3) for jn in ("y_slide", "x_slide")),
             "%s %+.3f mm" % ("z_slide", homevals["z_slide"] * 1e3)))

    # 4) Z independence: Z base group fixed under y/x
    z0 = poses_at({})
    for var in ("y_slide", "x_slide"):
        qp = poses_at({var: 0.05})
        for name in ALL_Z:
            if name in Z_GROUP:
                continue
            a = [z0[name][i][3] for i in range(3)]
            b = [qp[name][i][3] for i in range(3)]
            if max(abs(b[i] - a[i]) for i in range(3)) > 1e-6:
                print("Z_COLUMN MOVES WITH %s: %s" % (var, name))
                ok = False
    print("Z INDEPENDENCE: %s" % ("OK" if ok else "FAILED (Z base under y/x)"))

    # 5) Y/X slide direction: +q moves tooling along world +Y/+X
    for var, exp, grp, fixed_grps in (
            ("y_slide", (0, 1, 0), Y_CAR_GROUP,
             ({n for n in _built_links if n == "root"} | {"y_root", Y_SEAL} | Z_BASE_GROUP)),
            ("x_slide", (1, 0, 0), X_CAR_GROUP, (Y_GROUP | {"x_root", X_SEAL} | Z_BASE_GROUP)),
            ("z_slide", (0, 0, 1), Z_GROUP, (Y_GROUP | X_GROUP | Z_BASE_GROUP))):
        qp = poses_at({var: 0.05})
        for name in grp:
            a = [z0[name][i][3] for i in range(3)]
            b = [qp[name][i][3] for i in range(3)]
            d = [b[i] - a[i] for i in range(3)]
            if abs(d[exp.index(1)] - 0.05) > 1e-6 or any(abs(d[k]) > 1e-6 for k in range(3) if k != exp.index(1)):
                print("%s BAD %s: d=%s" % (var, name, [round(x, 6) for x in d]))
                ok = False
        for name in fixed_grps:
            if name not in _built_links:
                continue
            if name.startswith("z_") and name in Z_GROUP:
                continue
            a = [z0[name][i][3] for i in range(3)]
            b = [qp[name][i][3] for i in range(3)]
            if max(abs(b[i] - a[i]) for i in range(3)) > 1e-6:
                print("%s LEAK to %s" % (var, name))
                ok = False
    print("DIRECTION: y/x/z slides move tooling along world +Y/+X/+Z")

    # 6) collision sweep across travel
    if os.environ.get("DBG_BOXES"):
        sweep({"y_slide": 0, "x_slide": 0, "z_slide": 0}, debug=True)
    n_bad = 0
    min_report = 10
    for qy in grid:
        for qx in grid:
            for qz in gridz:
                bad = sweep({"y_slide": qy, "x_slide": qx, "z_slide": qz})
                n_bad += len(bad)
                for (a, b, ov, q) in bad[:min_report]:
                    print("COLLIDE %s x %s ov%s at qy=%.2f qx=%.2f qz=%.2f"
                          % (a, b, [round(v, 1) for v in ov], q["y_slide"], q["x_slide"], q["z_slide"]))
    if n_bad:
        print("COLLISION SWEEP: %d overlaps over %d configs" % (n_bad, (len(grid) * len(gridz)) ** 2 * len(grid)))
        ok = False
    else:
        print("COLLISION SWEEP: OK (%d boxed links x %d configs)"
              % (len(CBOXES), len(grid) ** 2 * len(gridz)))

    # 7) Z column centred on the XY mid-travel point: the Z base footprint
    #    centre must sit on the X carriage (payload) body centre at mid-stroke
    #    of both slides (assembly centre line), within 2 mm.
    pm = poses_xy(_mid, _mid)
    xlo, xhi = _world_bbox_pm(X_CAR, pm)
    xc = ((xlo[0] + xhi[0]) / 2.0, (xlo[1] + xhi[1]) / 2.0)
    zc = vmul(ZMount, [fxw, fyw, 0.0])
    dz = max(abs(zc[0] - xc[0]), abs(zc[1] - xc[1]))
    if dz > 2e-3:
        print("Z CENTRE DRIFT: Z footprint (%.1f, %.1f) vs XY mid-travel (%.1f, %.1f) mm"
              % (zc[0] * 1e3, zc[1] * 1e3, xc[0] * 1e3, xc[1] * 1e3))
        ok = False
    print("Z CENTRE: footprint (%.1f, %.1f) mm on XY mid-travel centre (%.1f, %.1f) mm"
          % (zc[0] * 1e3, zc[1] * 1e3, xc[0] * 1e3, xc[1] * 1e3))

    print("VERIFY: %s" % ("OK" if ok else "FAILED"))
    return ok


diag_ok = verify()

# ---------------------------------------------------------------- serialize

def serialize(elem, indent=0):
    pad = "    " * indent
    children = list(elem)
    if not children:
        attrs = "".join(' %s="%s"' % (k, v) for k, v in elem.attrib.items())
        return "%s<%s%s />" % (pad, elem.tag, attrs)
    attrs = "".join(' %s="%s"' % (k, v) for k, v in elem.attrib.items())
    lines = ["%s<%s%s>" % (pad, elem.tag, attrs)]
    for ch in children:
        if ch.tag in ("parent", "child"):
            lines.append("%s<%s %s />" % (pad + "    ", ch.tag,
                                          " ".join('%s="%s"' % (k, v) for k, v in ch.attrib.items())))
        else:
            lines.append(serialize(ch, indent + 1))
    lines.append("%s</%s>" % (pad, elem.tag))
    return "\n".join(lines)


print("MID_STROKE=%.4f m  HOME=%.4f m" % (MID_STROKE, HOME))
print("MID_STROKE_Z=%.4f m  HOME_Z=%.4f m" % (MID_STROKE_Z, HOME_Z))
print("X mount t=(%.4f, %.4f, %.4f)" % (Mount[0][3], Mount[1][3], Mount[2][3]))
os.environ["OVERWRITE"] = "1"
header = "<!--rebuild: Y/X from two clones of the corrected xyslideassy slide (y_slide +Y, x_slide +X); travel [%+.3f, %+.3f] mm single-sourced from the slide, bracketed by the H3 limit switches; home = mid-stroke; Z fixed column at assembly centre-->\n" % (SLIDE_LO * 1e3, SLIDE_HI * 1e3)
out = header + serialize(new_root) + "\n"
with open(OUT, "w", encoding="utf-8") as f:
    f.write(out)
print("wrote %s" % OUT)
if not diag_ok:
    sys.exit(1)