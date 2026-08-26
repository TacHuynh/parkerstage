#!/usr/bin/env python3
"""Rebuild parkerstage.urdf into a clean kinematic tree.

The Onshape export decomposed each stage's slide into 2 prismatic + 1
continuous joint (planar mate) with +/-10000 limits and arbitrarily cut the
tree, so nothing is a usable "slide". This script:

  * keeps every link (visuals/inertials verbatim), dropping the empty
    planar_* placeholder links,
  * roots the tree at `root` with base1 (401200xr__1_) rigidly attached,
  * adds a single prismatic joint `y_slide` between base1 and the bottom
    carriage group (plate + carriage_1 + accessories), axis = world Y,
  * bolts the upper stage base (401200xr__3_) to the plate,
  * adds a single prismatic joint `x_slide` between base2 and the top
    carriage group, axis = world X,
  * sets travel limits to +/-0.1 m (Parker 401200XR = 200 mm stroke).

All joint origins are computed numerically so every link keeps its exact
world pose at the zero configuration.
"""
import copy
import math
import os
import sys
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# read the pristine Onshape export if available (backed up before first run)
SRC = os.path.join(ROOT, "urdf", "parkerstage.urdf.onshape")
if not os.path.exists(SRC):
    SRC = os.path.join(ROOT, "urdf", "parkerstage.urdf")
OUT = os.path.join(ROOT, "urdf", "parkerstage.urdf")

TRAVEL = 0.10  # 200 mm stroke -> +/- 0.10 m about mid-stroke
# Physical mid-stroke: carriage center aligned with the base center.  Measured
# in the stage frame (base1/base2 local y): base center y=-13.017 mm, carriage
# center y=-33.750 mm -> q_mid = +20.733 mm.  The CAD export has the carriage
# 20.7 mm off mid-stroke, so the limits are shifted by this offset to keep the
# travel range centered on the physical mid-stroke instead of the CAD zero.
MID_STROKE = 0.020733
# Home-limit switch: the carriage flag trips the switch when the flag centre
# aligns with the switch centre.  Switch at base-y -21.838 mm, flag at
# carriage-y -33.750 mm -> q_home = +11.9 mm (Y stage).  After the X carriage
# flip the X stage carries its switch/flag posed exactly like the Y stage
# (switch at base2-y -21.837 mm, flag at carriage-y -33.750 mm), so the X home
# trips at the same +11.9 mm.
HOME_Y = 0.011912
HOME_X = 0.011913
EFFORT = 1.0
VELOCITY = 1.0
DENSITY = 2700.0  # kg/m^3, aluminum (Parker 401XR bases/carriages are aluminum)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mesh_inertia import mesh_inertia_voxel  # noqa: E402

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
    return [
        [R[0][0], R[0][1], R[0][2], t[0]],
        [R[1][0], R[1][1], R[1][2], t[1]],
        [R[2][0], R[2][1], R[2][2], t[2]],
        [0, 0, 0, 1],
    ]


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


# ---------------------------------------------------------------- parse old

tree = ET.parse(SRC)
root = tree.getroot()

links = {}
for l in root.findall("link"):
    links[l.get("name")] = copy.deepcopy(l)

children = {}
joints_by_name = {}
for j in root.findall("joint"):
    joints_by_name[j.get("name")] = j
    p = j.find("parent").get("link")
    c = j.find("child").get("link")
    children.setdefault(p, []).append((c, j))

# world poses at zero configuration
poses = {"root": t4([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0])}
order = ["root"]
while order:
    p = order.pop(0)
    for (c, j) in children.get(p, []):
        o = j.find("origin")
        xyz = list(map(float, o.get("xyz").split())) if o is not None else [0, 0, 0]
        rpy = list(map(float, o.get("rpy").split())) if o is not None else [0, 0, 0]
        poses[c] = mul(poses[p], t4(rpy2m(*rpy), xyz))
        order.append(c)


def rel(parent, child):
    """Origin (xyz, rpy) of `child` frame relative to `parent` frame at q=0."""
    M = mul(inv(poses[parent]), poses[child])
    t = [M[i][3] for i in range(3)]
    return t, rpy_of(M), M


def joint_axis(parent, child, world_axis):
    """Prismatic axis in the joint frame that maps to `world_axis`.

    axis_world = R_parent_world . R_origin . axis_joint
    with R_origin the joint-origin rotation (= rotation part of
    rel(parent, child)), so
    axis_joint = R_origin^T . R_parent^T . axis_world.
    """
    _, _, M = rel(parent, child)
    return joint_axis_from_origin(M, parent, world_axis)


def joint_axis_from_origin(M, parent, world_axis):
    """Like joint_axis, but taking an explicit parent->child origin transform M
    (e.g. one that has been flipped) instead of the untouched geometry."""
    R_origin = [[M[i][j] for j in range(3)] for i in range(3)]
    R_origin_T = [[R_origin[j][i] for j in range(3)] for i in range(3)]
    R_parent = [[poses[parent][i][j] for j in range(3)] for i in range(3)]
    R_parent_T = [[R_parent[j][i] for j in range(3)] for i in range(3)]
    return vrot(R_origin_T, vrot(R_parent_T, world_axis))


BASE1 = "401200xr__1_"
BASE2 = "401200xr__3_"
Y_CHILD = "plate"
X_CHILD = "401xr___encoder__401xr___encoder"

# ---- X-stage carriage flip -------------------------------------------------
# The Onshape export mounted the X-stage carriage 180 deg mirrored vs the
# Y-stage carriage (encoder/switch on opposite base sides).  We rotate the
# whole X carriage subtree 180 deg about the vertical (base2 z) axis so it is
# posed exactly like the Y carriage.  Carriage stays upright; the horizontal
# encoder/switch sides swap; travel stays world +X.
X_FLIP_GROUP = {
    "401xr___encoder__401xr___encoder",
    "401xr___carriage__401xr___carriage",
    "401xr___carriage_end_caps__401xr___carriage_end_caps_1",
    "401xr___encoder_base_2__401xr___encoder_base_2",
    "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch",
    "401xr___switch_flag__401xr___switch_flag",
}


def _mesh_center_y(name, base, pmap=None):
    """Mesh-bbox centre of link `name` along `base`'s y axis, in base coords.
    (Mirrors tools/home_offset.py: the carriage mesh centre is the sliding-body
    reference, and mid-stroke / home are measured from it.)  Pass the flipped
    pose map (poses_xflipped()) to measure the X stage as the rebuilt URDF
    actually poses it."""
    if pmap is None:
        pmap = poses
    l = links[name]
    vo = l.find("visual/origin")
    oxyz = [float(x) for x in vo.get("xyz").split()]
    fn = l.find("visual/geometry/mesh").get("filename").split("/")[-1]
    from mesh_inertia import read_triangles  # noqa: E402
    tris = read_triangles(os.path.join(MESH_DIR, fn))
    mn = [min(p[i] for tri in tris for p in tri) for i in range(3)]
    mx = [max(p[i] for tri in tris for p in tri) for i in range(3)]
    Mb = inv(pmap[base])
    cs = []
    for sx in (0, 1):
        for sy in (0, 1):
            for sz in (0, 1):
                p = [mn[i] + oxyz[i] if s == 0 else mx[i] + oxyz[i] for i, s in enumerate((sx, sy, sz))]
                cs.append(vmul(Mb, vmul(pmap[name], p)))
    lo = [min(c[i] for c in cs) for i in range(3)]
    hi = [max(c[i] for c in cs) for i in range(3)]
    return (lo[1] + hi[1]) / 2.0


# Rotate about the X carriage's own centre (base2-y through the carriage), not
# the mount, so the flipped carriage keeps the same along-slide position as the
# Y carriage (carriage centre y = -33.75 mm in each stage frame -> same q_mid).
X_CARRIAGE = "401xr___carriage__401xr___carriage"
X_CENTER_Y = _mesh_center_y(X_CARRIAGE, BASE2)
XFLIP4 = t4(rz(math.pi), [0.0, 2.0 * X_CENTER_Y, 0.0])  # 180 deg about base2 z axis thru (0, X_CENTER_Y, z)


def poses_xflipped():
    """World pose map identical to `poses`, but X-stage subtree links are given
    their post-flip world poses (used for collision-box placement and the
    pose-preservation re-check)."""
    # flipped world transform: root -> base2 -> (flip about base2 z) -> base2 -> world
    world_flip = mul(poses[BASE2], mul(XFLIP4, inv(poses[BASE2])))
    return {name: (mul(world_flip, M) if name in X_FLIP_GROUP else M)
            for name, M in poses.items()}


def _mesh_extent(name, base, pmap, axis, links_map=None):
    """(lo, hi) of link `name`'s mesh bbox along `base`'s `axis`, in base coords.
    Generalisation of _mesh_center_y for the Z stage (which slides along the
    base's z axis).  Pass links_map=zlinks for Z-export links."""
    if links_map is None:
        links_map = links
    l = links_map[name]
    vo = l.find("visual/origin")
    oxyz = [float(x) for x in vo.get("xyz").split()]
    fn = l.find("visual/geometry/mesh").get("filename").split("/")[-1]
    from mesh_inertia import read_triangles  # noqa: E402
    tris = read_triangles(os.path.join(MESH_DIR, fn))
    mn = [min(p[i] for tri in tris for p in tri) for i in range(3)]
    mx = [max(p[i] for tri in tris for p in tri) for i in range(3)]
    Mb = inv(pmap[base])
    cs = []
    for sx in (0, 1):
        for sy in (0, 1):
            for sz in (0, 1):
                p = [mn[i] + oxyz[i] if s == 0 else mx[i] + oxyz[i] for i, s in enumerate((sx, sy, sz))]
                cs.append(vmul(Mb, vmul(pmap[name], p)))
    lo = [min(c[i] for c in cs) for i in range(3)]
    hi = [max(c[i] for c in cs) for i in range(3)]
    return lo[axis], hi[axis]


# ---------------------------------------------------------------- Z stage
# Parker 401XR150 vertical slide mounted on the X carriage (the third axis of
# the XYZ compound).  The Onshape export (kept as urdf/zslide.urdf.onshape,
# links/joints namespaced z_*) has the slide axis along the base's local +z
# and the carriage exposed on the -y side (channel); the stage stands on its
# -z end (back-ballscrew-holder end) on the X carriage so the slide travels
# world +Z with the motor on top.  The in-plane rotation (Z +x -> base2 +y,
# Z +y -> base2 -x) centres the 41 x 33 mm base footprint on the X carriage
# and clears the X-stage encoder/switch accessories; the whole Z stage then
# sits above the X carriage top, so the full 150 mm stroke is collision-free.
ZSRC = os.path.join(ROOT, "urdf", "zslide.urdf.onshape")
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
        zposes[c] = mul(zposes[p], t4(rpy2m(*rpy), xyz))
        zorder.append(c)

Z_BASE = "z_401xr_150_base__401xr_150_base"
Z_CARRIAGE = "z_401xr___carriage__401xr___carriage"
Z_GROUP_ROOT = "z_c3_401xr__c3_401xr"  # motor: root of the fixed base group
Z_HOME_SW = "z_401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch"
Z_FLAG = "z_401xr___switch_flag__401xr___switch_flag"
Z_END_CAPS = "z_401xr___carriage_end_caps__401xr___carriage_end_caps"
Z_ENCODER = "z_401xr___encoder__401xr___encoder"
Z_ENC_BASE = "z_401xr___encoder_base_2__401xr___encoder_base_2_1"
TRAVEL_Z = 0.075  # Parker 401XR150 = 150 mm stroke -> +/- 0.075 m about mid-stroke

# ---- Z mount transform (base2 frame -> Z base frame) ----------------------
# Rotation: Z +x -> base2 +y, Z +y -> base2 -x, Z +z -> base2 +z (slide up).
Z_R = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
_zcx = _mesh_extent(Z_BASE, Z_BASE, zposes, 0, zlinks)
_zcy = _mesh_extent(Z_BASE, Z_BASE, zposes, 1, zlinks)
_zlo = _mesh_extent(Z_BASE, Z_BASE, zposes, 2, zlinks)[0]
_xcar_top = _mesh_extent(X_CARRIAGE, BASE2, poses_xflipped(), 2)[1]
Z_MOUNT_T = [(_zcy[0] + _zcy[1]) / 2.0 - MID_STROKE,          # centre the column over the XY travel centre (world -y; +mid via the y slide)
             _mesh_center_y(X_CARRIAGE, BASE2, poses_xflipped()) + MID_STROKE - (_zcx[0] + _zcx[1]) / 2.0,  # on X carriage centre at mid-travel (world +x; +mid via the x slide)
             _xcar_top + 0.0030 + 0.050 - _zlo]              # base2 z: base -z end above the XY sweep
# (the column is FIXED, so its base end must float clear of the XY travel: the
# X-base rails sweep beneath it, tallest at z = rail top; +0.0030 leaves a
# ~2.1 mm gap over the rails, and +0.050 raises the whole Z assembly 50 mm.
# Both are verified by the collision sweep below)
ZMOUNT = t4(Z_R, Z_MOUNT_T)  # base2 coords -> Z base coords


def z_poses_mounted():
    """World poses of the Z links: export poses re-rooted at the mounted Z base."""
    W = mul(poses[BASE2], ZMOUNT)  # Z base frame -> world
    return {name: mul(W, M) for name, M in zposes.items()}


ZM = z_poses_mounted()

# Z slide constants, measured like the Y/X stages (base centre - carriage
# centre, and L1 home switch - flag, along the slide axis); verify() re-checks
# them against the meshes.
MID_STROKE_Z = ((_mesh_extent(Z_BASE, Z_BASE, zposes, 2, zlinks)[0] + _mesh_extent(Z_BASE, Z_BASE, zposes, 2, zlinks)[1])
                - (_mesh_extent(Z_CARRIAGE, Z_BASE, zposes, 2, zlinks)[0] + _mesh_extent(Z_CARRIAGE, Z_BASE, zposes, 2, zlinks)[1])) / 2.0
HOME_Z = ((_mesh_extent(Z_HOME_SW, Z_BASE, zposes, 2, zlinks)[0] + _mesh_extent(Z_HOME_SW, Z_BASE, zposes, 2, zlinks)[1])
          - (_mesh_extent(Z_FLAG, Z_BASE, zposes, 2, zlinks)[0] + _mesh_extent(Z_FLAG, Z_BASE, zposes, 2, zlinks)[1])) / 2.0

# ---------------------------------------------------------------- build new

new_root = ET.Element("robot", {"name": "parkerstage"})

# links: all except empty planar_* placeholders
for name, l in links.items():
    if name.startswith("planar_"):
        continue
    new_root.append(l)

# Z-stage links (namespaced z_*; drop the export's planar placeholders/root)
for name, l in zlinks.items():
    if name.startswith("z_parallel_") or name == "z_root":
        continue
    new_root.append(l)

PLANAR_JOINTS = {
    "planar_1", "planar_1_1", "planar_1_2",
    "planar_2", "planar_2_1", "planar_2_2",
    "planar_3", "planar_3_1", "planar_3_2",
    "hanging_node_to_root_joint",
}


def make_joint(name, jtype, parent, child, xyz=(0, 0, 0), rpy=(0, 0, 0), axis=None, limit=None):
    j = ET.Element("joint", {"name": name, "type": jtype})
    if xyz != (0, 0, 0) or rpy != (0, 0, 0):
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


# 1) base1 rigidly fixed to root
t, rpy, _ = rel("root", BASE1)
new_root.append(make_joint("base1_fixed_to_root", "fixed", "root", BASE1, t, rpy))

# 2) Y slide: base1 -> plate (bottom stage carriage group)
t, rpy, _ = rel(BASE1, Y_CHILD)
y_axis = joint_axis(BASE1, Y_CHILD, (0, 1, 0))  # world +Y
new_root.append(make_joint(
    "y_slide", "prismatic", BASE1, Y_CHILD, t, rpy, axis=y_axis,
    limit={"effort": fmt(EFFORT), "velocity": fmt(VELOCITY),
           "lower": fmt(MID_STROKE - TRAVEL), "upper": fmt(MID_STROKE + TRAVEL)},
))

# 3) upper stage base bolted to the plate
t, rpy, _ = rel(Y_CHILD, BASE2)
new_root.append(make_joint("base2_mounted_to_plate", "fixed", Y_CHILD, BASE2, t, rpy))

# 4) X slide: base2 -> encoder (top stage carriage group)
#
# The Onshape export mounted the X-stage carriage as a 180 deg mirror of the
# Y-stage carriage relative to its (identical) base: encoder toward base -x
# and home-limit switch/flag toward base +x, whereas the Y carriage has the
# encoder toward base +x and switch/flag toward base -x.  Flip the X-stage
# carriage subtree 180 deg about the vertical (base z) axis through the
# carriage centre so it is posed exactly like the Y carriage (encoder toward
# base +x, switch/flag toward base -x).  Carriage stays upright and keeps its
# along-slide position; only the horizontal (encoder/switch) sides swap, so
# the travel axis still points along world +X.
t, rpy, M = rel(BASE2, X_CHILD)
Mx = mul(XFLIP4, M)  # base2 -> X child, flipped 180 deg about the vertical through the carriage centre
# the mount point itself moves through the flip -> update the joint origin
# translation (rotation about (x=0, y=X_CENTER_Y) maps y -> 2*X_CENTER_Y - y)
t = [t[0], 2.0 * X_CENTER_Y - t[1], t[2]]
rpy = rpy_of(Mx)
x_axis = joint_axis_from_origin(Mx, BASE2, (1, 0, 0))  # world +X
new_root.append(make_joint(
    "x_slide", "prismatic", BASE2, X_CHILD, t, rpy, axis=x_axis,
    limit={"effort": fmt(EFFORT), "velocity": fmt(VELOCITY),
           "lower": fmt(MID_STROKE - TRAVEL), "upper": fmt(MID_STROKE + TRAVEL)},
))

# 5) all remaining (fixed) joints, verbatim (axis is meaningless on fixed joints)
for j in root.findall("joint"):
    if j.get("name") in PLANAR_JOINTS:
        continue
    jc = copy.deepcopy(j)
    if jc.get("type") == "fixed":
        a = jc.find("axis")
        if a is not None:
            jc.remove(a)
    new_root.append(jc)

# 5b) Z stage: the motor-rooted base group is a FIXED column bolted to root
# (the world) at the centre of the assembly, independent of the XY motion; the
# z_slide prismatic then drives the Z carriage group along world +Z.  The
# mount pose is the Z base group's world pose at q=0 (ZMOUNT), which centres
# the 41 x 33 mm footprint on the X carriage mid position -- the assembly
# centre -- so the XY stages sweep beneath the column with the carriage clear.
Mt = ZM[Z_GROUP_ROOT]
new_root.append(make_joint("z_mounted_to_root", "fixed", "root", Z_GROUP_ROOT,
                           [Mt[i][3] for i in range(3)], rpy_of(Mt)))

# z_slide: Z base -> Z carriage.  Origin = rel(Z_BASE, Z_CARRIAGE) from the
# export; axis = world +Z expressed in the joint frame.
Mzc = mul(inv(zposes[Z_BASE]), zposes[Z_CARRIAGE])
zt, zrpy = [Mzc[i][3] for i in range(3)], rpy_of(Mzc)
R_orig = [[Mzc[i][j] for j in range(3)] for i in range(3)]
R_orig_T = [[R_orig[j][i] for j in range(3)] for i in range(3)]
R_zbase = [[ZM[Z_BASE][i][j] for j in range(3)] for i in range(3)]
R_zbase_T = [[R_zbase[j][i] for j in range(3)] for i in range(3)]
z_axis = vrot(R_orig_T, vrot(R_zbase_T, (0, 0, 1)))  # world +Z in the joint frame
new_root.append(make_joint(
    "z_slide", "prismatic", Z_BASE, Z_CARRIAGE, zt, zrpy, axis=z_axis,
    limit={"effort": fmt(EFFORT), "velocity": fmt(VELOCITY),
           "lower": fmt(MID_STROKE_Z - TRAVEL_Z), "upper": fmt(MID_STROKE_Z + TRAVEL_Z)},
))

# 5c) remaining Z fixed joints, verbatim (planar chain + root joint dropped)
Z_PLANAR = {"z_parallel_1", "z_parallel_1_1", "z_parallel_1_2", "z_parallel_1_3",
            "z_hanging_node_to_root_joint"}
for j in zroot.findall("joint"):
    if j.get("name") in Z_PLANAR:
        continue
    jc = copy.deepcopy(j)
    if jc.get("type") == "fixed":
        a = jc.find("axis")
        if a is not None:
            jc.remove(a)
    new_root.append(jc)

# 6) realistic inertials from mesh volumes (aluminum density)

def fnum(v):
    return "%.9g" % v


def fix_inertials(new_root):
    rows = []
    for l in new_root.findall("link"):
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
        I_l = mmul(R, mmul(I_m, R_T))  # inertia about COM, aligned with link frame

        inert = ET.Element("inertial")
        m = ET.SubElement(inert, "mass")
        m.set("value", fnum(mass))
        o = ET.SubElement(inert, "origin")
        o.set("xyz", fmt_vec(com_l))
        o.set("rpy", "0 0 0")
        ia = ET.SubElement(inert, "inertia")
        ia.set("ixx", fnum(I_l[0][0]))
        ia.set("ixy", fnum(I_l[0][1]))
        ia.set("ixz", fnum(I_l[0][2]))
        ia.set("iyy", fnum(I_l[1][1]))
        ia.set("iyz", fnum(I_l[1][2]))
        ia.set("izz", fnum(I_l[2][2]))

        old = l.find("inertial")
        if old is not None:
            l.remove(old)
        l.insert(0, inert)
        rows.append((l.get("name"), V * 1e6, mass * 1000))
    return rows


print("link                                vol[cm^3]  mass[g]")
for name, vcm, g in fix_inertials(new_root):
    print("%-35s %9.2f %9.2f" % (name, vcm, g))



# ---------------------------------------------------------------- collisions

# Collision boxes, defined as (ref_frame, xyz_min, xyz_max) in the STAGE
# reference frame (base1 frame for the Y-stage group, base2 frame for the
# X-stage group).  The CAD export interpenetrates in a few places (the
# carriage top plate sits inside the base rail zone, the plate's top edge
# overlaps the upper-stage bottom flange, accessories overlap the base
# walls), so the boxes are built from the *clean* interfaces with explicit
# gaps instead of the raw mesh envelopes:
#   * bases are hollow: a lower body + two outer rail/wall strips, leaving
#     the channel open for the carriage,
#   * the carriage / plate / accessories get central or outward boxes that
#     clear the base geometry by >= 0.5 mm,
#   * the base body stops at the flange (-32.5 mm); the feet below that
#     only exist near the y-ends and nothing moves under them.
COLLISION = {
    # --- bases (identical part, own frame) ---
    "401200xr__1_": [
        ("401200xr__1_", [-0.0205, -0.2052, -0.0325], [0.0205, 0.1791, -0.0080]),
        ("401200xr__1_", [-0.0205, -0.2052, -0.0080], [-0.0185, 0.1791, 0.0100]),
        ("401200xr__1_", [0.0185, -0.2052, -0.0080], [0.0205, 0.1791, 0.0100]),
    ],
    "401200xr__3_": [
        ("401200xr__3_", [-0.0205, -0.2052, -0.0325], [0.0205, 0.1791, -0.0080]),
        ("401200xr__3_", [-0.0205, -0.2052, -0.0080], [-0.0185, 0.1791, 0.0100]),
        ("401200xr__3_", [0.0185, -0.2052, -0.0080], [0.0205, 0.1791, 0.0100]),
    ],
    # --- Y-stage group (ref = base1 frame) ---
    "plate": [
        ("401200xr__1_", [-0.0205, -0.0563, 0.0105], [0.0205, -0.0113, 0.0185]),
    ],
    "401xr___carriage__401xr___carriage_1": [
        ("401200xr__1_", [-0.0175, -0.0578, -0.0070], [0.0175, -0.0135, 0.0092]),
    ],
    "401xr___encoder__401xr___encoder_1": [
        ("401200xr__1_", [0.0207, -0.0587, -0.0239], [0.0365, -0.0127, 0.0092]),
    ],
    "401xr___carriage_end_caps__401xr___carriage_end_caps_1_1": [
        ("401200xr__1_", [-0.0175, -0.0814, 0.0000], [0.0175, -0.0583, 0.0085]),
    ],
    "401xr___carriage_end_caps__401xr___carriage_end_caps_2": [
        ("401200xr__1_", [-0.0175, -0.0123, 0.0000], [0.0175, 0.0139, 0.0085]),
    ],
    "401xr___encoder_base_2__401xr___encoder_base_2_1": [
        ("401200xr__1_", [0.0207, -0.0456, -0.0320], [0.0285, -0.0306, -0.0255]),
    ],
    "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch_1": [
        ("401200xr__1_", [-0.0265, -0.0396, -0.0126], [-0.0207, -0.0041, -0.0026]),
    ],
    "401xr___switch_flag__401xr___switch_flag_1": [
        ("401200xr__1_", [-0.0285, -0.0415, 0.0000], [-0.0207, -0.0260, 0.0072]),
    ],
    # --- X-stage group (ref = base2 frame) ---
    "401xr___carriage__401xr___carriage": [
        ("401200xr__3_", [-0.0175, -0.0578, -0.0070], [0.0175, -0.0135, 0.0092]),
    ],
    "401xr___encoder__401xr___encoder": [
        ("401200xr__3_", [-0.0365, -0.0549, -0.0239], [-0.0207, -0.0089, 0.0092]),
    ],
    "401xr___carriage_end_caps__401xr___carriage_end_caps_1": [
        ("401200xr__3_", [-0.0175, -0.0123, 0.0000], [0.0175, 0.0139, 0.0085]),
    ],
    "401xr___encoder_base_2__401xr___encoder_base_2": [
        ("401200xr__3_", [-0.0285, -0.0369, -0.0320], [-0.0207, -0.0219, -0.0255]),
    ],
    "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch": [
        ("401200xr__3_", [0.0207, -0.0634, -0.0126], [0.0265, -0.0280, -0.0026]),
    ],
    "401xr___switch_flag__401xr___switch_flag": [
        ("401200xr__3_", [0.0207, -0.0415, 0.0000], [0.0285, -0.0260, 0.0072]),
    ],
}

# --- Z-stage group (ref = Z base frame).  The 401XR150 base is a hollow
# extrusion with a central drive-screw slot, so it is boxed as two side walls
# + a top wall; the carriage / end caps / accessories get narrow central or
# outboard boxes that clear the base geometry by >= 0.5 mm (the raw CAD
# interpenetrates: the carriage flanges tuck under the side walls).
COLLISION_Z = {
    Z_BASE: [
        (Z_BASE, [-0.0665, -0.04529, -0.10348], [-0.05173, -0.0120, 0.16352]),  # side wall -x
        (Z_BASE, [0.0775, -0.04529, -0.10348], [0.09263, -0.0120, 0.16352]),   # side wall +x
        (Z_BASE, [-0.0665, -0.0165, -0.10348], [0.0775, -0.0120, 0.16352]),    # top wall
    ],
    Z_CARRIAGE: [
        (Z_BASE, [0.0680, -0.0445, -0.0135], [0.0760, -0.0355, 0.0315]),
    ],
    Z_END_CAPS: [
        (Z_BASE, [0.0680, -0.0510, -0.0387], [0.0760, -0.0440, -0.0150]),  # -z end
        (Z_BASE, [0.0680, -0.0510, 0.0330], [0.0760, -0.0440, 0.0566]),    # +z end
    ],
    Z_ENCODER: [
        (Z_BASE, [0.0935, -0.0510, -0.0150], [0.1075, -0.0200, 0.0295]),
    ],
    Z_ENC_BASE: [
        (Z_BASE, [0.0935, -0.0175, -0.0020], [0.1000, -0.0100, 0.0115]),
    ],
    Z_HOME_SW: [
        (Z_BASE, [0.0455, -0.0405, 0.0035], [0.0510, -0.0310, 0.0380]),
    ],
    Z_FLAG: [
        (Z_BASE, [0.0425, -0.0500, 0.0015], [0.0500, -0.0410, 0.0165]),
    ],
}

# link-frame description of every collision box: (link, ref, center, rpy, size)
def collision_boxes():
    pf = poses_xflipped()  # X-group links given their post-flip world pose
    pf.update(ZM)          # Z links at their mounted poses
    out = {}
    for link, boxes in {**COLLISION, **COLLISION_Z}.items():
        lst = []
        for (ref, mn, mx) in boxes:
            if link in X_FLIP_GROUP:
                # the X-subtree was flipped 180 deg about the vertical, so the
                # box (defined in base2 frame) moves with it: map its 8 corners
                # through XFLIP4 and retake the AABB.
                cps = [vmul(XFLIP4, [mx[i] if s else mn[i] for i, s in enumerate(sw)]) for sw in
                       [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
                        [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]]
                mn = [min(c[i] for c in cps) for i in range(3)]
                mx = [max(c[i] for c in cps) for i in range(3)]
            M = mul(inv(pf[link]), pf[ref])  # ref -> link (flipped/mounted pose map)
            c = vmul(M, [(mn[i] + mx[i]) / 2 for i in range(3)])
            Rm = [[M[i][j] for j in range(3)] for i in range(3)]
            lst.append((ref, c, rpy_of(Rm), [mx[i] - mn[i] for i in range(3)]))
        out[link] = lst
    return out

CBOXES = collision_boxes()


def add_collisions(new_root):
    for l in new_root.findall("link"):
        for (ref, c, rpy, size) in CBOXES.get(l.get("name"), []):
            col = ET.SubElement(l, "collision")
            o = ET.SubElement(col, "origin")
            o.set("xyz", fmt_vec(c))
            o.set("rpy", fmt_vec(rpy))
            g = ET.SubElement(col, "geometry")
            b = ET.SubElement(g, "box")
            b.set("size", fmt_vec(size))


add_collisions(new_root)
print("collision boxes: %d total" % sum(len(v) for v in CBOXES.values()))

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


header = "<!--URDF generated by ONSHAPE BY PTC INC, 1.219-->\n"
header += "<!--rebuild: planar joints removed; y_slide (Y) + x_slide (X) prismatic, 200 mm stroke; z_slide (Z) 150 mm; Z column fixed to root at assembly centre, independent of XY-->\n"
header += "<!--travel centered on physical mid-stroke (q=+%g m Y/X, +%g m Z); home trips at q=%+g m (Y) / %+g m (X) / %+g m (Z)-->\n" % (
    MID_STROKE, MID_STROKE_Z, HOME_Y, HOME_X, HOME_Z)
out = header + serialize(new_root) + "\n"
with open(OUT, "w", encoding="utf-8") as f:
    f.write(out)
print("wrote %s" % OUT)

# ---------------------------------------------------------------- verify

def poses_at(qs):
    """World poses of every link for joint values qs = {'y_slide': v, 'x_slide': v}."""
    ch = {}
    for j in new_root.findall("joint"):
        ch.setdefault(j.find("parent").get("link"), []).append(j)
    pos = {"root": t4([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0])}
    order = ["root"]
    while order:
        p = order.pop(0)
        for j in ch.get(p, []):
            o = j.find("origin")
            xyz = list(map(float, o.get("xyz").split())) if o is not None else [0, 0, 0]
            rpy = list(map(float, o.get("rpy").split())) if o is not None else [0, 0, 0]
            M = t4(rpy2m(*rpy), xyz)
            if j.get("name") in qs:
                a = j.find("axis")
                ax = list(map(float, a.get("xyz").split()))
                M = mul(M, t4([[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                              [qs[j.get("name")] * ax[0], qs[j.get("name")] * ax[1], qs[j.get("name")] * ax[2]]))
            c = j.find("child").get("link")
            pos[c] = mul(pos[p], M)
            order.append(c)
    return pos

Y_GROUP = {"plate", "401xr___carriage__401xr___carriage_1",
           "401xr___encoder_base_2__401xr___encoder_base_2_1",
           "401xr___encoder__401xr___encoder_1",
           "401xr___carriage_end_caps__401xr___carriage_end_caps_1_1",
           "401xr___carriage_end_caps__401xr___carriage_end_caps_2",
           "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch_1",
           "401xr___switch_flag__401xr___switch_flag_1", BASE2,
           "401xr___encoder__401xr___encoder",
           "401xr___carriage__401xr___carriage",
           "401xr___carriage_end_caps__401xr___carriage_end_caps_1",
           "401xr___encoder_base_2__401xr___encoder_base_2",
           "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch",
           "401xr___switch_flag__401xr___switch_flag"}
X_GROUP = {"401xr___encoder__401xr___encoder",
           "401xr___carriage__401xr___carriage",
           "401xr___carriage_end_caps__401xr___carriage_end_caps_1",
           "401xr___encoder_base_2__401xr___encoder_base_2",
           "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch",
           "401xr___switch_flag__401xr___switch_flag"}
# The Z stage is a FIXED column at the assembly centre (bolted to root), so
# it moves with neither the Y nor the X slide; only the Z carriage group moves
# with the z_slide.  (The Z export's planar placeholders are dropped from the
# tree, so they are excluded from the group checks.)
Z_ALL = {n for n in zlinks if not n.startswith("z_parallel_") and n != "z_root"}
Z_GROUP = {Z_CARRIAGE, Z_END_CAPS,
           "z_401xr___carriage_end_caps__401xr___carriage_end_caps_1__1_",
           Z_ENCODER, Z_ENC_BASE, Z_HOME_SW, Z_FLAG}
Z_BASE_GROUP = Z_ALL - Z_GROUP  # fixed links: the column itself

def verify():
    ok = True
    # --- collision sweep: every box pair of different links, across travel ---
    # all boxes are axis-aligned in base1 or base2 frames, and the two stages
    # differ by a 90 deg rotation about z, so transforming corners to the
    # base1 frame leaves every box axis-aligned -> exact AABB overlap test.
    Mb1 = inv(poses["401200xr__1_"])
    box_world_corners = {}
    for link, lst in CBOXES.items():
        corners = []
        for (ref, c, rpy, size) in lst:
            Rm = rpy2m(*rpy)
            half = [s / 2 for s in size]
            M = t4(Rm, c)  # link frame -> collision box frame
            cs = []
            for sx in (0, 1):
                for sy in (0, 1):
                    for sz in (0, 1):
                        p = [half[i] if s == 0 else -half[i] for i, s in enumerate((sx, sy, sz))]
                        cs.append(vmul(M, p))  # in link frame
            corners.append(cs)
        box_world_corners[link] = corners

    def sweep(qs, debug=False):
        pos = poses_at(qs)
        aabbs = {}  # link -> list of (mn, mx) in base1 frame
        for link, lst in box_world_corners.items():
            aa = []
            for cs in lst:
                pts = [vmul(Mb1, vmul(pos[link], p)) for p in cs]
                mn = [min(p[i] for p in pts) for i in range(3)]
                mx = [max(p[i] for p in pts) for i in range(3)]
                aa.append((mn, mx))
            aabbs[link] = aa
        if debug:
            for link, aa in sorted(aabbs.items()):
                for (mn, mx) in aa:
                    print("  AABB %-50s x [%7.2f, %7.2f] y [%7.2f, %7.2f] z [%7.2f, %7.2f]" % (
                        link, mn[0] * 1000, mx[0] * 1000, mn[1] * 1000, mx[1] * 1000,
                        mn[2] * 1000, mx[2] * 1000))
        bad = []
        names = list(aabbs)
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                for (mn1, mx1) in aabbs[names[i]]:
                    for (mn2, mx2) in aabbs[names[j]]:
                        ov = [max(0.0, min(mx1[k], mx2[k]) - max(mn1[k], mn2[k])) for k in range(3)]
                        if ov[0] > 0 and ov[1] > 0 and ov[2] > 0:
                            bad.append((names[i], names[j],
                                        [round(1e3 * v, 1) for v in ov], qs.copy()))
        return bad

    if os.environ.get("DBG_BOXES"):
        print("collision box AABBs in base1 frame at q=0:")
        sweep({"y_slide": 0.0, "x_slide": 0.0}, debug=True)

    lo, hi = MID_STROKE - TRAVEL, MID_STROKE + TRAVEL
    mid = (lo + hi) / 2
    grid = [lo, mid - 0.05, mid, mid + 0.05, hi]
    loz, hiz = MID_STROKE_Z - TRAVEL_Z, MID_STROKE_Z + TRAVEL_Z
    midz = (loz + hiz) / 2
    gridz = [loz, midz - 0.0375, midz, midz + 0.0375, hiz]
    # home positions must sit inside the travel range
    if not (lo <= HOME_Y <= hi and lo <= HOME_X <= hi and loz <= HOME_Z <= hiz):
        print("HOME OUTSIDE TRAVEL: Y %+.4f X %+.4f Z %+.4f vs [%+.4f, %+.4f]/[%+.4f, %+.4f]"
              % (HOME_Y, HOME_X, HOME_Z, lo, hi, loz, hiz))
        ok = False
    else:
        print("HOME in travel: Y %+.4f m, X %+.4f m, Z %+.4f m (ranges [%+.4f, %+.4f] / [%+.4f, %+.4f])"
              % (HOME_Y, HOME_X, HOME_Z, lo, hi, loz, hiz))

    # every slide joint must declare exactly its Parker stroke spec (200 mm for
    # the 401200XR Y/X stages, 150 mm for the 401XR150 Z stage; the two 200 mm
    # slides are also identical to each other).  The spec values are literals so
    # the check still bites if TRAVEL / TRAVEL_Z (and therefore the generated
    # limits) ever drift from the hardware.
    STROKE_SPEC = {"y_slide": 0.20, "x_slide": 0.20, "z_slide": 0.15}  # m, Parker specs
    strokes = {}
    for jn in ("y_slide", "x_slide", "z_slide"):
        for jj in new_root.findall("joint"):
            if jj.get("name") == jn:
                lim = jj.find("limit")
                strokes[jn] = float(lim.get("upper")) - float(lim.get("lower"))
    if any(jn not in strokes for jn in ("y_slide", "x_slide", "z_slide")):
        print("STROKE MISSING joint: %s" % [jn for jn in ("y_slide", "x_slide", "z_slide") if jn not in strokes])
        ok = False
    else:
        for jn, st in strokes.items():
            if abs(st - STROKE_SPEC[jn]) > 1e-6:
                print("STROKE MISMATCH %s: %.4f m (expected %.4f m = %.0f mm)"
                      % (jn, st, STROKE_SPEC[jn], STROKE_SPEC[jn] * 1e3))
                ok = False
        if abs(strokes["y_slide"] - strokes["x_slide"]) > 1e-6:
            print("STROKE MISMATCH y_slide vs x_slide: %.4f m vs %.4f m"
                  % (strokes["y_slide"], strokes["x_slide"]))
            ok = False
        print("STROKE: %s" % ", ".join("%s %.1f mm" % (jn, st * 1e3) for jn, st in strokes.items()))

    # both slides must be centred on the *physical* mid-stroke: the declared
    # MID_STROKE constant and the generated limits must both agree with the
    # geometry (base centre - carriage centre, measured like home_offset.py),
    # so the travel range can't shift away from the physical mid-stroke.
    MID_TOL = 2e-4  # m, absorbs mesh-bbox rounding (~0.002 mm); catches real drift
    y_car = "401xr___carriage__401xr___carriage_1"
    phys = {}
    for tag, base, car, pmap, axis, lmap in (("Y", BASE1, y_car, poses, 1, links),
                                             ("X", BASE2, X_CARRIAGE, poses_xflipped(), 1, links),
                                             ("Z", Z_BASE, Z_CARRIAGE, zposes, 2, zlinks)):
        e0 = _mesh_extent(base, base, pmap, axis, lmap)
        e1 = _mesh_extent(car, base, pmap, axis, lmap)
        phys[tag] = (e0[0] + e0[1] - e1[0] - e1[1]) / 2.0
        declared = MID_STROKE if tag != "Z" else MID_STROKE_Z
        if abs(phys[tag] - declared) > MID_TOL:
            print("MID-STROKE MISMATCH %s: physical %+.4f m vs declared %+.4f m"
                  % (tag, phys[tag], declared))
            ok = False
    for jn, st in strokes.items():
        cent = None
        for jj in new_root.findall("joint"):
            if jj.get("name") == jn:
                lim = jj.find("limit")
                cent = (float(lim.get("lower")) + float(lim.get("upper"))) / 2.0
        declared = MID_STROKE if jn != "z_slide" else MID_STROKE_Z
        if cent is None or abs(cent - declared) > 1e-6:
            print("MID-STROKE MISMATCH %s limits: centre %+.4f m vs declared %+.4f m" % (jn, cent, declared))
            ok = False
    print("MID-STROKE: %s" % ", ".join(
        "%s physical %+.3f mm / limits centred %+.3f mm"
        % (t, phys[t] * 1e3, (MID_STROKE if t != "Z" else MID_STROKE_Z) * 1e3)
        for t in ("Y", "X", "Z")))

    # both home switches must trip at the geometry-derived position and stay
    # inside travel: q_home = switch centre - flag centre along the slide axis,
    # measured from the meshes in each stage frame (like home_offset.py).  This
    # guards both the declared HOME_* constants and the switch/flag geometry.
    homes = {}
    for tag, base, pmap, sw, fl, hconst, axis, lmap in (
            ("Y", BASE1, poses,
             "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch_1",
             "401xr___switch_flag__401xr___switch_flag_1", HOME_Y, 1, links),
            ("X", BASE2, poses_xflipped(),
             "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch",
             "401xr___switch_flag__401xr___switch_flag", HOME_X, 1, links),
            ("Z", Z_BASE, zposes, Z_HOME_SW, Z_FLAG, HOME_Z, 2, zlinks)):
        es = _mesh_extent(sw, base, pmap, axis, lmap)
        ef = _mesh_extent(fl, base, pmap, axis, lmap)
        qh = (es[0] + es[1] - ef[0] - ef[1]) / 2.0
        homes[tag] = qh
        if abs(qh - hconst) > MID_TOL:
            print("HOME MISMATCH %s: physical %+.4f m vs declared %+.4f m" % (tag, qh, hconst))
            ok = False
        if not (lo <= qh <= hi if tag != "Z" else loz <= qh <= hiz):
            print("HOME OUTSIDE TRAVEL %s: %+.4f m (range [%+.4f, %+.4f])"
                  % (tag, qh, (lo, hi) if tag != "Z" else (loz, hiz)))
            ok = False
    print("HOME: %s" % ", ".join(
        "%s physical %+.3f mm (declared %+.3f mm, inside travel)"
        % (t, homes[t] * 1e3, {"Y": HOME_Y, "X": HOME_X, "Z": HOME_Z}[t] * 1e3)
        for t in ("Y", "X", "Z")))
    n_bad = 0
    min_report = 5
    for qy in grid:
        for qx in grid:
            for qz in gridz:
                bad = sweep({"y_slide": qy, "x_slide": qx, "z_slide": qz})
                n_bad += len(bad)
                for (a, b, ov, q) in bad[:min_report]:
                    print("COLLIDE %s x %s overlap(mm) %s at qy=%.2f qx=%.2f qz=%.2f"
                          % (a, b, ov, q["y_slide"], q["x_slide"], q["z_slide"]))
    if n_bad:
        print("COLLISION SWEEP: %d overlaps over %d configs" % (n_bad, len(grid) ** 3))
        ok = False
    else:
        print("COLLISION SWEEP: OK (%d box pairs x %d configs)" % (
            sum(len(v) for v in box_world_corners.values()) ** 2, len(grid) ** 3))

    z0 = poses_at({"y_slide": 0.0, "x_slide": 0.0, "z_slide": 0.0})
    # pose preservation at zero config vs original export.  The X-stage carriage
    # subtree is deliberately flipped (X_FLIP_GROUP), so its links no longer match
    # the raw export; those are re-checked against the intended flipped pose below.
    for name, M in poses.items():
        if name.startswith("planar_") or name in X_FLIP_GROUP:
            continue
        a = [z0[name][i][3] for i in range(3)]
        b = [M[i][3] for i in range(3)]
        if max(abs(x - y) for x, y in zip(a, b)) > 1e-6:
            print("POSE MISMATCH at q=0: %s %s vs %s" % (name, a, b))
            ok = False
    # the flipped X subtree must sit exactly at its post-flip world pose
    pf = poses_xflipped()
    for name in X_FLIP_GROUP:
        a = [z0[name][i][3] for i in range(3)]
        b = [pf[name][i][3] for i in range(3)]
        if max(abs(x - y) for x, y in zip(a, b)) > 1e-6:
            print("X-FLIP POSE MISMATCH: %s %s vs %s" % (name, a, b))
            ok = False
    # Z links must sit at their mounted poses (export poses re-rooted at the
    # mounted Z base)
    for name, M in ZM.items():
        if name == "z_root" or name.startswith("z_parallel_"):
            continue
        a = [z0[name][i][3] for i in range(3)]
        b = [M[i][3] for i in range(3)]
        if max(abs(x - y) for x, y in zip(a, b)) > 1e-6:
            print("Z POSE MISMATCH at q=0: %s" % name)
            ok = False

    # flipped X collision boxes must sit on the *orientation* the Y carriage has in
    # its own base frame: encoder/readhead toward base +x, home-limit switch/flag
    # toward base -x, same z band.  (The y extent differs because the two stages
    # are stacked at different y offsets; the flip preserves the assembly's own
    # geometry, so only the x/z side placement is checked.)
    def _xbox(mn, mx):
        cps = [vmul(XFLIP4, [mx[i] if s else mn[i] for i, s in enumerate(sw)]) for sw in
               [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1],
                [1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]]]
        return [min(c[i] for c in cps) for i in range(3)], \
               [max(c[i] for c in cps) for i in range(3)]
    XVS_Y = {"401xr___encoder__401xr___encoder": "401xr___encoder__401xr___encoder_1",
             "401xr___carriage_end_caps__401xr___carriage_end_caps_1": "401xr___carriage_end_caps__401xr___carriage_end_caps_1_1",
             "401xr___encoder_base_2__401xr___encoder_base_2": "401xr___encoder_base_2__401xr___encoder_base_2_1",
             "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch": "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch_1",
             "401xr___switch_flag__401xr___switch_flag": "401xr___switch_flag__401xr___switch_flag_1"}
    for xl, yl in XVS_Y.items():
        for (_, fx_mn, fx_mx), (_, fy_mn, fy_mx) in zip(COLLISION[xl], COLLISION[yl]):
            xmn, xmx = _xbox(fx_mn, fx_mx)
            if max(abs(xmn[0] - fy_mn[0]), abs(xmx[0] - fy_mx[0]),
                   abs(xmn[2] - fy_mn[2]), abs(xmx[2] - fy_mx[2])) > 0.002:
                print("X-FLIP BOX MISMATCH %s vs %s: X%s Y%s" % (xl, yl, [xmn[0], xmn[2], xmx[0], xmx[2]],
                                                                 [fy_mn[0], fy_mn[2], fy_mx[0], fy_mx[2]]))
                ok = False
    # y_slide moves Y group along world Y only
    zy = poses_at({"y_slide": 0.05, "x_slide": 0.0, "z_slide": 0.0})
    for name in Y_GROUP:
        a = [z0[name][i][3] for i in range(3)]
        b = [zy[name][i][3] for i in range(3)]
        d = [b[i] - a[i] for i in range(3)]
        if abs(d[1] - 0.05) > 1e-6 or abs(d[0]) > 1e-6 or abs(d[2]) > 1e-6:
            print("Y_SLIDE BAD for %s: d=%s" % (name, [round(x, 6) for x in d]))
            ok = False
    for name in set(poses) - Y_GROUP - {"root"}:
        if name.startswith("planar_"):
            continue
        a = [z0[name][i][3] for i in range(3)]
        b = [zy[name][i][3] for i in range(3)]
        if max(abs(b[i] - a[i]) for i in range(3)) > 1e-6:
            print("Y_SLIDE LEAK to %s" % name)
            ok = False
    # x_slide moves X group along world X only
    zx = poses_at({"y_slide": 0.0, "x_slide": 0.05, "z_slide": 0.0})
    for name in X_GROUP:
        a = [z0[name][i][3] for i in range(3)]
        b = [zx[name][i][3] for i in range(3)]
        d = [b[i] - a[i] for i in range(3)]
        if abs(d[0] - 0.05) > 1e-6 or abs(d[1]) > 1e-6 or abs(d[2]) > 1e-6:
            print("X_SLIDE BAD for %s: d=%s" % (name, [round(x, 6) for x in d]))
            ok = False
    for name in set(poses) - X_GROUP - {"root"}:
        if name.startswith("planar_"):
            continue
        a = [z0[name][i][3] for i in range(3)]
        b = [zx[name][i][3] for i in range(3)]
        if max(abs(b[i] - a[i]) for i in range(3)) > 1e-6:
            print("X_SLIDE LEAK to %s" % name)
            ok = False
    # z_slide moves the Z carriage group along world Z only
    zz = poses_at({"y_slide": 0.0, "x_slide": 0.0, "z_slide": 0.05})
    for name in Z_GROUP:
        a = [z0[name][i][3] for i in range(3)]
        b = [zz[name][i][3] for i in range(3)]
        d = [b[i] - a[i] for i in range(3)]
        if abs(d[2] - 0.05) > 1e-6 or abs(d[0]) > 1e-6 or abs(d[1]) > 1e-6:
            print("Z_SLIDE BAD for %s: d=%s" % (name, [round(x, 6) for x in d]))
            ok = False
    for name in (set(poses) | Z_ALL) - Z_GROUP - {"root"}:
        if name.startswith("planar_") or name.startswith("z_parallel_"):
            continue
        a = [z0[name][i][3] for i in range(3)]
        b = [zz[name][i][3] for i in range(3)]
        if max(abs(b[i] - a[i]) for i in range(3)) > 1e-6:
            print("Z_SLIDE LEAK to %s" % name)
            ok = False
    # the Z column is a FIXED frame at the assembly centre: no Z link may move
    # with y_slide or x_slide (the Y/X leak loops only scan the XY export's
    # link set, so a re-mount to a moving link would otherwise go unnoticed)
    z_indep = True
    for jn, qp in (("y_slide", zy), ("x_slide", zx)):
        for name in Z_ALL:
            a = [z0[name][i][3] for i in range(3)]
            b = [qp[name][i][3] for i in range(3)]
            if max(abs(b[i] - a[i]) for i in range(3)) > 1e-6:
                print("Z_COLUMN MOVES WITH %s: %s d=%s"
                      % (jn, name, [round(b[i] - a[i], 6) for i in range(3)]))
                ok = z_indep = False
    print("Z INDEPENDENCE: %s (%d Z links fixed under y_slide/x_slide)"
          % ("OK" if z_indep else "FAILED", len(Z_ALL)))

    # the Z column is a FIXED frame, so it must be centred over the CENTRE of
    # the XY travel envelope (both slides at MID_STROKE), not over any single
    # posed position: only then does the column sit on the assembly centre line
    # across the whole travel.  Check its base footprint centre coincides (in
    # the x-y plane) with the X carriage centre at mid-travel.  A regression
    # that shifts the column off this line is caught here even though the
    # stroke/mid-stroke/home guards are invariant to a planar mount shift.
    z_cent = True
    ZCENT_TOL = 0.0015  # ~1.5 mm: absorbs mesh-bbox rounding
    _ze0 = _mesh_extent(Z_BASE, Z_BASE, zposes, 0, zlinks)
    _ze1 = _mesh_extent(Z_BASE, Z_BASE, zposes, 1, zlinks)
    _cc0 = _mesh_extent(X_CARRIAGE, BASE2, poses_xflipped(), 0, links)
    _cc1 = _mesh_extent(X_CARRIAGE, BASE2, poses_xflipped(), 1, links)
    _zmid = poses_at({"y_slide": MID_STROKE, "x_slide": MID_STROKE, "z_slide": 0.0})
    Mzb = ZM[Z_BASE]                       # world pose of the Z base frame at q=0
    Mcx = _zmid[X_CARRIAGE]                # world pose of the X carriage at mid-travel

    def _wcent(M, lo, hi):
        c = [(lo[k] + hi[k]) / 2.0 for k in range(3)]
        return [sum(M[i][k] * c[k] for k in range(3)) + M[i][3] for i in range(3)]

    wz = _wcent(Mzb, [_ze0[0], _ze1[0], 0], [_ze0[1], _ze1[1], 0])
    wx = _wcent(Mcx, [_cc0[0], _cc1[0], 0], [_cc0[1], _cc1[1], 0])
    for ax, sym in ((0, "x"), (1, "y")):
        if abs(wz[ax] - wx[ax]) > ZCENT_TOL:
            print("Z OFFSET %s: column centre %+.4f m vs X carriage mid-travel %+.4f m"
                  % (sym, wz[ax], wx[ax]))
            ok = z_cent = False
    print("Z CENTERED: %s (base footprint centre (%+.4f, %+.4f) on the mid-travel centre)"
          % ("OK" if z_cent else "FAILED", wz[0], wz[1]))
    print("VERIFY: %s" % ("OK" if ok else "FAILED"))
    return ok

if not verify():
    sys.exit(1)
