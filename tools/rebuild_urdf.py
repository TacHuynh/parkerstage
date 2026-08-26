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

# ---------------------------------------------------------------- build new

new_root = ET.Element("robot", {"name": "parkerstage"})

# links: all except empty planar_* placeholders
for name, l in links.items():
    if name.startswith("planar_"):
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

# link-frame description of every collision box: (link, ref, center, rpy, size)
def collision_boxes():
    pf = poses_xflipped()  # X-group links given their post-flip world pose
    out = {}
    for link, boxes in COLLISION.items():
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
            M = mul(inv(pf[link]), poses[ref])  # ref -> link (flipped link pose for X group)
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
header += "<!--rebuild: planar joints removed; y_slide (Y) + x_slide (X) prismatic, 200 mm stroke-->\n"
header += "<!--travel centered on physical mid-stroke (q=+%g m); home switch trips at q=%+g m (Y) / %+g m (X)-->\n" % (
    MID_STROKE, HOME_Y, HOME_X)
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
    # home positions must sit inside the travel range
    if not (lo <= HOME_Y <= hi and lo <= HOME_X <= hi):
        print("HOME OUTSIDE TRAVEL: Y %+.4f X %+.4f vs [%+.4f, %+.4f]" % (HOME_Y, HOME_X, lo, hi))
        ok = False
    else:
        print("HOME in travel: Y %+.4f m, X %+.4f m (range [%+.4f, %+.4f])" % (HOME_Y, HOME_X, lo, hi))

    # both slide joints must declare exactly the Parker 401200XR 200 mm stroke
    # (and be identical to each other).  The 0.20 m is the *spec* value, kept
    # as a literal here so the check still bites if TRAVEL (and therefore the
    # generated limits) is ever changed away from the 200 mm spec.
    STROKE_SPEC = 0.20  # m, Parker 401200XR
    strokes = {}
    for jn in ("y_slide", "x_slide"):
        for jj in new_root.findall("joint"):
            if jj.get("name") == jn:
                lim = jj.find("limit")
                strokes[jn] = float(lim.get("upper")) - float(lim.get("lower"))
    if any(jn not in strokes for jn in ("y_slide", "x_slide")):
        print("STROKE MISSING joint: %s" % [jn for jn in ("y_slide", "x_slide") if jn not in strokes])
        ok = False
    else:
        for jn, st in strokes.items():
            if abs(st - STROKE_SPEC) > 1e-6:
                print("STROKE MISMATCH %s: %.4f m (expected %.4f m = 200 mm)" % (jn, st, STROKE_SPEC))
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
    for tag, base, car, pmap in (("Y", BASE1, y_car, poses),
                                 ("X", BASE2, X_CARRIAGE, poses_xflipped())):
        phys[tag] = _mesh_center_y(base, base, pmap) - _mesh_center_y(car, base, pmap)
        if abs(phys[tag] - MID_STROKE) > MID_TOL:
            print("MID-STROKE MISMATCH %s: physical %+.4f m vs declared %+.4f m"
                  % (tag, phys[tag], MID_STROKE))
            ok = False
    for jn, st in strokes.items():
        lo = None
        for jj in new_root.findall("joint"):
            if jj.get("name") == jn:
                lim = jj.find("limit")
                lo = (float(lim.get("lower")) + float(lim.get("upper"))) / 2.0
        if lo is None or abs(lo - MID_STROKE) > 1e-6:
            print("MID-STROKE MISMATCH %s limits: centre %+.4f m vs declared %+.4f m" % (jn, lo, MID_STROKE))
            ok = False
    print("MID-STROKE: %s" % ", ".join("%s physical %+.3f mm / limits centred %+.3f mm" % (t, phys[t] * 1e3, MID_STROKE * 1e3)
                                        for t in ("Y", "X")))
    n_bad = 0
    min_report = 5
    for qy in grid:
        for qx in grid:
            bad = sweep({"y_slide": qy, "x_slide": qx})
            n_bad += len(bad)
            for (a, b, ov, q) in bad[:min_report]:
                print("COLLIDE %s x %s overlap(mm) %s at qy=%.2f qx=%.2f" % (a, b, ov, q["y_slide"], q["x_slide"]))
    if n_bad:
        print("COLLISION SWEEP: %d overlaps over %d configs" % (n_bad, len(grid) ** 2))
        ok = False
    else:
        print("COLLISION SWEEP: OK (%d box pairs x %d configs)" % (
            sum(len(v) for v in box_world_corners.values()) ** 2, len(grid) ** 2))

    z0 = poses_at({"y_slide": 0.0, "x_slide": 0.0})
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
    zy = poses_at({"y_slide": 0.05, "x_slide": 0.0})
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
    zx = poses_at({"y_slide": 0.0, "x_slide": 0.05})
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
    print("VERIFY: %s" % ("OK" if ok else "FAILED"))
    return ok

if not verify():
    sys.exit(1)
