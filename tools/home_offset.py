#!/usr/bin/env python3
"""Compute the physical mid-stroke and the home-limit-switch home position.

For each stage (base1 -> Y group, base2 -> X group) in the stage's own frame:

  q_mid  = base_center - carriage_assembly_center
           (the prismatic displacement that centers the carriage on the base)

  q_home = switch_sense_y - flag_y
           (the displacement that puts the carriage flag at the home switch)

All measurements come from the mesh bounding boxes expressed in the stage
reference frame (base1 frame for the Y stage, base2 frame for the X stage).
"""
import sys
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import xml.etree.ElementTree as ET  # noqa: E402
import rebuild_urdf as R  # noqa: E402
from mesh_inertia import read_triangles  # noqa: E402

URDF = os.path.join(ROOT, "urdf", "parkerstage.urdf")
MESHES = os.path.join(ROOT, "meshes")

t = ET.parse(URDF)
by = {l.get("name"): l for l in t.getroot().findall("link")}


def link_bbox(name):
    l = by[name]
    vo = l.find("visual/origin")
    oxyz = [float(x) for x in vo.get("xyz").split()]
    fn = l.find("visual/geometry/mesh").get("filename").split("/")[-1]
    tris = read_triangles(os.path.join(MESHES, fn))
    mn = [min(p[i] for tri in tris for p in tri) for i in range(3)]
    mx = [max(p[i] for tri in tris for p in tri) for i in range(3)]
    return [mn[i] + oxyz[i] for i in range(3)], [mx[i] + oxyz[i] for i in range(3)]


def bbox_in(name, base, poses=None):
    """Link visual bbox expressed in `base`'s frame.

    `poses` defaults to the raw Onshape world poses; pass the X-flipped pose
    map (R.poses_xflipped()) for the X stage so measurements reflect the
    flipped carriage that the rebuilt URDF actually contains."""
    if poses is None:
        poses = R.poses
    mn, mx = link_bbox(name)
    Mb = poses[base]
    corners = []
    for sx in (0, 1):
        for sy in (0, 1):
            for sz in (0, 1):
                p = [mn[i] if s == 0 else mx[i] for i, s in enumerate((sx, sy, sz))]
                corners.append(R.vmul(R.inv(Mb), R.vmul(poses[name], p)))
    return ([min(p[i] for p in corners) for i in range(3)],
            [max(p[i] for p in corners) for i in range(3)])


def center(lo, hi):
    return [(lo[i] + hi[i]) / 2 for i in range(3)]


def stage(axis, base, carriage, endcap_a, endcap_b, switch, flag, tag, poses=None, travel=0.1):
    """axis: local index (0=x,1=y,2=z) of the slide direction.
    `poses` defaults to the raw Onshape poses; pass the X-flipped pose map for
    the X stage (see bbox_in) and the Z-export pose map (R.zposes) for the Z
    stage, whose slide runs along the base local z.  `travel` is the half
    stroke (0.1 m for the 401200XR Y/X stages, 0.075 m for the 401XR150 Z)."""
    if poses is None:
        poses = R.poses
    bmn, bmx = bbox_in(base, base, poses)
    base_c = center(bmn, bmx)[axis]

    # carriage center is the sliding-body reference (the end caps on the Y
    # stage are symmetric about the carriage centre; the X-stage export is
    # missing one end cap, so the carriage itself is the reference there)
    cmn, cmx = bbox_in(carriage, base, poses)
    car_c = center(cmn, cmx)[axis]

    # end-cap symmetry check (Y stage): assembly centre == carriage centre?
    lo, hi = bbox_in(carriage, base, poses)
    for ec in (endcap_a, endcap_b):
        elo, ehi = bbox_in(ec, base, poses)
        lo = [min(lo[i], elo[i]) for i in range(3)]
        hi = [max(hi[i], ehi[i]) for i in range(3)]
    asm_c = center(lo, hi)[axis]

    smn, smx = bbox_in(switch, base, poses)
    fmn, fmx = bbox_in(flag, base, poses)
    sw_c = center(smn, smx)[axis]
    fl_c = center(fmn, fmx)[axis]

    q_mid = base_c - car_c
    q_home = sw_c - fl_c
    print("%s  (slide local axis %s)" % (tag, "xyz"[axis]))
    print("  base center          : %9.3f mm" % (1e3 * base_c))
    print("  carriage center      : %9.3f mm" % (1e3 * car_c))
    print("  +endcaps center      : %9.3f mm (delta %.2f mm)" % (1e3 * asm_c, 1e3 * (asm_c - car_c)))
    print("  switch center        : %9.3f mm" % (1e3 * sw_c))
    print("  flag center          : %9.3f mm" % (1e3 * fl_c))
    print("  q_mid  (carriage onto base center)  : %+8.3f mm" % (1e3 * q_mid))
    print("  q_home (flag at switch)             : %+8.3f mm" % (1e3 * q_home))
    print("  limits centered on mid-stroke: [%+.4f, %+.4f]" % (q_mid - travel, q_mid + travel))
    print()


# Y stage (base1): slide along base1 local y
stage(1, "401200xr__1_",
      "401xr___carriage__401xr___carriage_1",
      "401xr___carriage_end_caps__401xr___carriage_end_caps_1_1",
      "401xr___carriage_end_caps__401xr___carriage_end_caps_2",
      "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch_1",
      "401xr___switch_flag__401xr___switch_flag_1",
      "Y stage (bottom)")

# X stage (base2): slide along base2 local y.  The rebuilt URDF flips the X
# carriage to match the Y stage, so measure with the flipped pose map; the
# export has only one end cap on the X carriage, so the carriage is the
# reference there.
stage(1, "401200xr__3_",
      "401xr___carriage__401xr___carriage",
      "401xr___carriage_end_caps__401xr___carriage_end_caps_1",
      "401xr___carriage_end_caps__401xr___carriage_end_caps_1",
      "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch",
      "401xr___switch_flag__401xr___switch_flag",
      "X stage (top)",
      poses=R.poses_xflipped())

# Z stage (mounted on the X carriage): slide along the Z base local z (the
# stage stands on its -z end, so the slide is world +Z).  Measured from the
# Z-export poses in the Z base frame; the along-slide offsets are invariant
# under the mount.  Half stroke 75 mm (401XR150).
stage(2, "z_401xr_150_base__401xr_150_base",
      "z_401xr___carriage__401xr___carriage",
      "z_401xr___carriage_end_caps__401xr___carriage_end_caps",
      "z_401xr___carriage_end_caps__401xr___carriage_end_caps",
      "z_401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch",
      "z_401xr___switch_flag__401xr___switch_flag",
      "Z stage (vertical)",
      poses=R.zposes,
      travel=0.075)
