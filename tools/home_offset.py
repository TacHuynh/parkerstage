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


def bbox_in(name, base):
    """Link visual bbox expressed in `base`'s frame."""
    mn, mx = link_bbox(name)
    Mb = R.poses[base]
    corners = []
    for sx in (0, 1):
        for sy in (0, 1):
            for sz in (0, 1):
                p = [mn[i] if s == 0 else mx[i] for i, s in enumerate((sx, sy, sz))]
                corners.append(R.vmul(R.inv(Mb), R.vmul(R.poses[name], p)))
    return ([min(p[i] for p in corners) for i in range(3)],
            [max(p[i] for p in corners) for i in range(3)])


def center(lo, hi):
    return [(lo[i] + hi[i]) / 2 for i in range(3)]


def stage(axis, base, carriage, endcap_a, endcap_b, switch, flag, tag):
    """axis: local index (0=x,1=y,2=z) of the slide direction."""
    bmn, bmx = bbox_in(base, base)
    base_c = center(bmn, bmx)[axis]

    # carriage center is the sliding-body reference (the end caps on the Y
    # stage are symmetric about the carriage centre; the X-stage export is
    # missing one end cap, so the carriage itself is the reference there)
    cmn, cmx = bbox_in(carriage, base)
    car_c = center(cmn, cmx)[axis]

    # end-cap symmetry check (Y stage): assembly centre == carriage centre?
    lo, hi = bbox_in(carriage, base)
    for ec in (endcap_a, endcap_b):
        elo, ehi = bbox_in(ec, base)
        lo = [min(lo[i], elo[i]) for i in range(3)]
        hi = [max(hi[i], ehi[i]) for i in range(3)]
    asm_c = center(lo, hi)[axis]

    smn, smx = bbox_in(switch, base)
    fmn, fmx = bbox_in(flag, base)
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
    print("  limits centered on mid-stroke: [%+.4f, %+.4f]" % (q_mid - 0.1, q_mid + 0.1))
    print()


# Y stage (base1): slide along base1 local y
stage(1, "401200xr__1_",
      "401xr___carriage__401xr___carriage_1",
      "401xr___carriage_end_caps__401xr___carriage_end_caps_1_1",
      "401xr___carriage_end_caps__401xr___carriage_end_caps_2",
      "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch_1",
      "401xr___switch_flag__401xr___switch_flag_1",
      "Y stage (bottom)")

# X stage (base2): slide along base2 local y (mirrored part; the export has
# only one end cap on the X carriage, so the carriage itself is the reference)
stage(1, "401200xr__3_",
      "401xr___carriage__401xr___carriage",
      "401xr___carriage_end_caps__401xr___carriage_end_caps_1",
      "401xr___carriage_end_caps__401xr___carriage_end_caps_1",
      "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch",
      "401xr___switch_flag__401xr___switch_flag",
      "X stage (top)")
