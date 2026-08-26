#!/usr/bin/env python3
"""Generate docs/collision_geometry.svg from urdf/parkerstage.urdf.

Renders the collision model of the 401200XR compound XYZ stage:
  * top view (XY)  -- every collision box, color-coded by kinematic group,
    with the slide travel sweep (envelope + extreme positions);
  * sections A-A (XZ) and B-B (YZ) -- channel cross-sections cut through each
    XY carriage at q = 0, showing the clearances (Z column excluded so the
    XY-carriage detail stays readable; the Z stage is shown in the top view);
  * a per-link inventory of all 30 boxes.

Run:  python3 tools/make_collision_svg.py
"""
import math
import os
import xml.etree.ElementTree as ET

import make_tree_svg as TREE  # reuse SHORT names + FONT (import is side-effect free)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
URDF = os.path.join(ROOT, "urdf", "parkerstage.urdf")
OUT = os.path.join(ROOT, "docs", "collision_geometry.svg")

FONT = TREE.FONT
esc = TREE.esc
measure = TREE.measure

# link -> kinematic group (mirrors tools/make_tree_svg.py)
X_ASSEMBLY = {"401xr___encoder__401xr___encoder",
              "401xr___carriage__401xr___carriage",
              "401xr___carriage_end_caps__401xr___carriage_end_caps_1",
              "401xr___encoder_base_2__401xr___encoder_base_2",
              "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch",
              "401xr___switch_flag__401xr___switch_flag"}
Y_ASSEMBLY = {"401xr___carriage__401xr___carriage_1",
              "401xr___encoder_base_2__401xr___encoder_base_2_1",
              "401xr___encoder__401xr___encoder_1",
              "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch_1",
              "401xr___switch_flag__401xr___switch_flag_1",
              "401xr___carriage_end_caps__401xr___carriage_end_caps_1_1",
              "401xr___carriage_end_caps__401xr___carriage_end_caps_2"}

# group -> (fill, stroke)   [light fill, dark stroke -- same family as the tree diagram]
GROUP = {
    "Y base":     ("#e0e7ea", "#78909c"),
    "X base":     ("#cfd8dc", "#546e7a"),
    "plate":      ("#f0f0f0", "#9e9e9e"),
    "Y assembly": ("#e8eaf6", "#7986cb"),
    "X assembly": ("#e3f2fd", "#1e88e5"),
    "Z stage":    ("#fff3e0", "#fb8c00"),
}


def rpy2m(r, p, y):
    cr, sr = math.cos(r), math.sin(r)
    cp, sp = math.cos(p), math.sin(p)
    cy, sy = math.cos(y), math.sin(y)
    return [[cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr],
            [sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr],
            [-sp,     cp * sr,               cp * cr]]


def mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)]
            + [A[i][3] + sum(A[i][k] * B[k][3] for k in range(3))] for i in range(3)]


def t4(R, xyz):
    return [R[i] + [xyz[i]] for i in range(3)]


def parse():
    tree = ET.parse(URDF)
    root = tree.getroot()
    joints = {j.get("name"): j for j in root.findall("joint")}
    links = {l.get("name"): l for l in root.findall("link")}
    return joints, links


def poses_at(joints, q):
    pos = {"root": t4([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [0, 0, 0])}
    order = ["root"]
    while order:
        p = order.pop(0)
        for j in joints.values():
            if j.find("parent").get("link") != p:
                continue
            c = j.find("child").get("link")
            o = j.find("origin")
            xyz = list(map(float, o.get("xyz").split()))
            rpy = list(map(float, o.get("rpy").split()))
            M = mul(pos[p], t4(rpy2m(*rpy), xyz))
            if j.get("type") == "prismatic":
                ax = list(map(float, j.find("axis").get("xyz").split()))
                d = q.get(j.get("name"), 0.0)
                R = [[M[i][k] for k in range(3)] for i in range(3)]
                w = [sum(R[i][k] * ax[k] for k in range(3)) for i in range(3)]
                for i in range(3):
                    M[i][3] += d * w[i]
            pos[c] = M
            order.append(c)
    return pos


def world_aabbs(links, pos):
    """link -> list of world-space (min_xyz, max_xyz) per collision box."""
    out = {}
    for name, l in links.items():
        for c in l.findall("collision"):
            o = c.find("origin")
            xyz = list(map(float, o.get("xyz").split()))
            rpy = list(map(float, o.get("rpy").split()))
            sz = list(map(float, c.find("geometry/box").get("size").split()))
            R = rpy2m(*rpy)
            M = mul(pos[name], t4(R, xyz))
            hs = [s / 2 for s in sz]
            corners = []
            for k in range(8):
                v = [hs[i] if not (k >> i) & 1 else -hs[i] for i in range(3)]
                corners.append([M[i][3] + sum(M[i][j] * v[j] for j in range(3)) for i in range(3)])
            mn = [min(c[i] for c in corners) for i in range(3)]
            mx = [max(c[i] for c in corners) for i in range(3)]
            out.setdefault(name, []).append((mn, mx))
    return out


def group_of(link):
    if link in X_ASSEMBLY:
        return "X assembly"
    if link in Y_ASSEMBLY:
        return "Y assembly"
    if link.startswith("z_"):
        return "Z stage"
    if link == "plate":
        return "plate"
    if link == "401200xr__1_":
        return "Y base"
    if link == "401200xr__3_":
        return "X base"
    raise ValueError(link)


MM = 1e3


def main():
    joints, links = parse()
    lo = float(joints["y_slide"].find("limit").get("lower"))   # meters
    hi = float(joints["y_slide"].find("limit").get("upper"))

    q0 = {"y_slide": 0.0, "x_slide": 0.0, "z_slide": 0.0}
    pos0 = poses_at(joints, q0)
    boxes0 = world_aabbs(links, pos0)
    assert len([b for v in boxes0.values() for b in v]) == 30, "expected 30 boxes (20 XY + 10 Z)"

    # travel envelope (axis-aligned prismatic motion preserves orientation)
    def sweep_rect(mn, mx, moves_x, moves_y):
        r = [list(mn), list(mx)]
        if moves_y:  # q in [lo, hi] along world +Y
            r[0][1] += lo
            r[1][1] += hi
        if moves_x:  # q in [lo, hi] along world +X
            r[0][0] += lo
            r[1][0] += hi
        return r

    # the Z stage is a FIXED column at the assembly centre: it moves with
    # neither slide (only its carriage rides z_slide, which is vertical)
    moving_x = set(X_ASSEMBLY)
    moving_y = set(Y_ASSEMBLY) | {"plate", "401200xr__3_"}

    # extreme configs for dashed outlines
    ext = []
    for qy in (lo, hi):
        for qx in (lo, hi):
            ext.append(({"y_slide": qy / MM, "x_slide": qx / MM}))
    ext_boxes = [world_aabbs(links, poses_at(joints, q)) for q in ext]

    # --- section cut planes (through each carriage centre at q = 0) --------
    y0 = (boxes0["401xr___carriage__401xr___carriage_1"][0][0][1]
          + boxes0["401xr___carriage__401xr___carriage_1"][0][1][1]) / 2
    x0 = (boxes0["401xr___carriage__401xr___carriage"][0][0][0]
          + boxes0["401xr___carriage__401xr___carriage"][0][1][0]) / 2

    def section(boxes, axis, cut):
        """Rectangles (a0,a1) x (b0,b1) of boxes intersecting the cut plane."""
        rects = []
        for name, lst in boxes.items():
            for (mn, mx) in lst:
                if mn[axis] <= cut <= mx[axis]:
                    if axis == 1:      # cut in y -> (x, z)
                        rects.append((mn[0], mn[2], mx[0], mx[2], group_of(name), name))
                    else:              # cut in x -> (y, z)
                        rects.append((mn[1], mn[2], mx[1], mx[2], group_of(name), name))
        return rects

    secA = section(boxes0, 1, y0)   # XZ through the Y carriage
    secB = section(boxes0, 0, x0)   # YZ through the X carriage
    # keep the sections XY-carriage focused: the tall Z column would collapse
    # the panel scale, so drop Z boxes from the cut planes (they are in the top view)
    secA = [r for r in secA if not r[5].startswith("z_")]
    secB = [r for r in secB if not r[5].startswith("z_")]

    # --- clearances (from q=0 AABBs) ----------------------------------------
    yslab, yrail1, yrail2 = (boxes0["401200xr__1_"][0], boxes0["401200xr__1_"][1],
                             boxes0["401200xr__1_"][2])
    ycar = boxes0["401xr___carriage__401xr___carriage_1"][0]
    plate = boxes0["plate"][0]
    xslab, xrail1, xrail2 = (boxes0["401200xr__3_"][0], boxes0["401200xr__3_"][1],
                             boxes0["401200xr__3_"][2])
    xcar = boxes0["401xr___carriage__401xr___carriage"][0]
    gap_side_y = (ycar[0][0] - yrail1[1][0], yrail2[0][0] - ycar[1][0])  # to the rail inner faces
    gap_floor_y = ycar[0][2] - yslab[1][2]
    gap_top_y = yrail1[1][2] - ycar[1][2]
    gap_plate = plate[0][2] - yrail1[1][2]
    gap_side_x = (xcar[0][1] - xrail1[1][1], xrail2[0][1] - xcar[1][1])  # to the rail inner faces
    gap_floor_x = xcar[0][2] - xslab[1][2]
    gap_top_x = xrail1[1][2] - xcar[1][2]

    # --- SVG layout ----------------------------------------------------------
    S = []
    add = S.append

    def new_panel(title):
        return {"title": title, "ox0": 0.0, "oy0": 0.0, "scale": 1.0}

    # --- panel geometry (mm extents) ----------------------------------------
    def extent(rects):
        mn = [1e9, 1e9]
        mx = [-1e9, -1e9]
        for (x0, y0, x1, y1, g, l) in rects:
            mn[0] = min(mn[0], x0)
            mn[1] = min(mn[1], y0)
            mx[0] = max(mx[0], x1)
            mx[1] = max(mx[1], y1)
        return mn, mx

    # top view: q0 boxes + sweep envelopes + extreme outlines
    top_rects = []
    env_rects = []
    for name, lst in boxes0.items():
        grp = group_of(name)
        mx_, my_ = name in moving_x, name in moving_y
        for (mn, mx) in lst:
            top_rects.append((mn[0], mn[1], mx[0], mx[1], grp, name))
            if mx_ or my_:
                e = sweep_rect(mn, mx, mx_, my_)
                env_rects.append((e[0][0], e[0][1], e[1][0], e[1][1], grp, name))
    tmn, tmx = extent(top_rects)
    emn, emx = extent(env_rects)
    allmn = [min(tmn[i], emn[i]) for i in range(2)]
    allmx = [max(tmx[i], emx[i]) for i in range(2)]
    PAD = 0.022   # metres (22 mm of margin around the sweep)
    for i in range(2):
        allmn[i] -= PAD
        allmx[i] += PAD
    top = new_panel("Top view (XY)  ·  q = 0, with travel sweep")
    top["scale"] = min(640 / (allmx[0] - allmn[0]), 560 / (allmx[1] - allmn[1]))
    top["ox0"] = allmn[0]
    top["oy0"] = allmx[1]   # mm max-y maps to panel top

    amn, amx = extent(secA)
    bmn, bmx = extent(secB)
    for i in range(2):
        amn[i] -= 0.015   # 15 mm of margin (metres)
        amx[i] += 0.015
        bmn[i] -= 0.015
        bmx[i] += 0.015
    secA_p = new_panel("Section A-A (XZ, through Y carriage at q = 0)")
    secA_p["scale"] = min(660 / (amx[0] - amn[0]), 240 / (amx[1] - amn[1]))
    secA_p["ox0"] = amn[0]
    secA_p["oy0"] = amx[1]   # mm max-z maps to panel top
    secB_p = new_panel("Section B-B (YZ, through X carriage at q = 0)")
    secB_p["scale"] = min(660 / (bmx[0] - bmn[0]), 240 / (bmx[1] - bmn[1]))
    secB_p["ox0"] = bmn[0]
    secB_p["oy0"] = bmx[1]   # mm max-z maps to panel top

    W = 1500
    # place panels: top full width, sections side by side
    place = {}
    place["top"] = (60, 150)
    place["secA"] = (60, 150 + (allmx[1] - allmn[1]) * top["scale"] + 80)
    place["secB"] = (60 + 700, 150 + (allmx[1] - allmn[1]) * top["scale"] + 80)
    H = place["secA"][1] + (amx[1] - amn[1]) * secA_p["scale"] + 420

    add('<?xml version="1.0" encoding="UTF-8"?>')
    add('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
        'viewBox="0 0 %d %d" font-family="%s">' % (W, H, W, H, FONT))
    add('<rect x="0" y="0" width="%d" height="%d" fill="#ffffff"/>' % (W, H))
    add("""<defs>
  <marker id="arr" markerWidth="9" markerHeight="9" refX="7" refY="4" orient="auto">
    <path d="M0 0 L8 4 L0 8 z" fill="#37474f"/>
  </marker>
</defs>""")

    def panel_origin(p, name):
        ox, oy = place[name]
        p["ox_px"], p["oy_px"] = ox, oy
        add('<text x="%d" y="%d" font-size="14" font-weight="bold" fill="#1a1a1a">%s</text>'
            % (ox, oy - 6, esc(p["title"])))

    def px2(p, x, y):
        return (p["ox_px"] + (x - p["ox0"]) * p["scale"],
                p["oy_px"] + (p["oy0"] - y) * p["scale"])

    nboxes = sum(len(v) for v in boxes0.values())
    nlinks = len(boxes0)

    # ---- header ----
    add('<text x="60" y="46" font-size="17" font-weight="bold" fill="#1a1a1a">'
        'parkerstage  ·  401200XR compound XYZ stage  ·  collision geometry</text>')
    add('<text x="60" y="68" font-size="11" fill="#546e7a">'
        '%d boxes across %d links, color-coded by kinematic group (same groups as the kinematic tree).</text>'
        % (nboxes, nlinks))
    add('<text x="60" y="84" font-size="11" fill="#546e7a">'
        'Boxes are built from the clean interfaces with explicit gaps; the sweep (y_slide, x_slide, z_slide) '
        'was verified zero-overlap over 125 travel configurations.</text>')

    # ---- top view ----
    panel_origin(top, "top")
    # sweep envelopes first (light), then extreme outlines, then q0 boxes
    for (x0, y0, x1, y1, grp, name) in env_rects:
        (ax, ay) = px2(top, x0, y1)
        (bx, by) = px2(top, x1, y0)
        fill = GROUP[grp][1]
        add('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" '
            'fill-opacity="0.10" stroke="none"/>' % (ax, by, bx - ax, ay - by, fill))
    for eb in ext_boxes:
        for name, lst in eb.items():
            if group_of(name) in ("Y base",):
                continue
            grp = group_of(name)
            for (mn, mx) in lst:
                (ax, ay) = px2(top, mn[0], mx[1])
                (bx, by) = px2(top, mx[0], mn[1])
                add('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" '
                    'stroke="%s" stroke-width="0.9" stroke-opacity="0.55" stroke-dasharray="4 3"/>'
                    % (ax, by, bx - ax, ay - by, GROUP[grp][1]))
    for (x0, y0, x1, y1, grp, name) in top_rects:
        (ax, ay) = px2(top, x0, y1)
        (bx, by) = px2(top, x1, y0)
        fill, stroke = GROUP[grp]
        add('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" stroke="%s" '
            'stroke-width="1.1" fill-opacity="0.6"/>' % (ax, by, bx - ax, ay - by, fill, stroke))

    # cut-plane indicators (label at the start of each line)
    for (axis, cut, label) in ((1, y0, "A"), (0, x0, "B")):
        if axis == 1:
            (p1x, p1y) = px2(top, allmn[0], cut)
            (p2x, p2y) = px2(top, allmx[0], cut)
            (lx, ly) = (p1x + 5, p1y - 5)
        else:
            (p1x, p1y) = px2(top, cut, allmn[1])
            (p2x, p2y) = px2(top, cut, allmx[1])
            (lx, ly) = (p1x + 5, p1y + 14)
        add('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#b71c1c" '
            'stroke-width="1" stroke-dasharray="8 4" stroke-opacity="0.6"/>' % (p1x, p1y, p2x, p2y))
        add('<text x="%.1f" y="%.1f" font-size="11" font-weight="bold" fill="#b71c1c">%s</text>'
            % (lx, ly, label))

    # travel arrows + labels (along the sweep envelope edges)
    def arrow(x0, y0, x1, y1, label):
        (ax, ay) = px2(top, x0, y0)
        (bx, by) = px2(top, x1, y1)
        add('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#37474f" '
            'stroke-width="1.4" marker-end="url(#arr)"/>' % (ax, ay, bx, by))
        add('<text x="%.1f" y="%.1f" font-size="11" fill="#37474f" text-anchor="middle">%s</text>'
            % ((ax + bx) / 2, (ay + by) / 2 - 6, esc(label)))
    arrow(allmx[0] - 6, emn[1], allmx[0] - 6, emx[1],
          "y_slide  q in [%+.1f, %+.1f] mm" % (lo * MM, hi * MM))
    arrow(emn[0], allmn[1] + 6, emx[0], allmn[1] + 6,
          "x_slide  q in [%+.1f, %+.1f] mm" % (lo * MM, hi * MM))

    # axis triad at a panel corner (pixel coords), outside the data
    def triad(p, name, pw, ph, h_label, v_label, h_len=40, v_len=30):
        ox = place[name][0] + pw - 90
        oy = place[name][1] + ph + 16
        add('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#455a64" '
            'stroke-width="1.6" marker-end="url(#arr)"/>' % (ox, oy, ox + h_len, oy))
        add('<text x="%.1f" y="%.1f" font-size="11" fill="#455a64">%s</text>'
            % (ox + h_len + 3, oy + 4, esc(h_label)))
        add('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#455a64" '
            'stroke-width="1.6" marker-end="url(#arr)"/>' % (ox, oy, ox, oy - v_len))
        add('<text x="%.1f" y="%.1f" font-size="11" fill="#455a64">%s</text>'
            % (ox + 3, oy - v_len - 4, esc(v_label)))
    triad(top, "top", (allmx[0] - allmn[0]) * top["scale"],
          (allmx[1] - allmn[1]) * top["scale"], "X →", "Y ↑")

    # ---- sections ----
    for (p, name, rects) in ((secA_p, "secA", secA), (secB_p, "secB", secB)):
        panel_origin(p, name)
        for (x0, y0, x1, y1, grp, link) in rects:
            (ax, ay) = px2(p, x0, y1)
            (bx, by) = px2(p, x1, y0)
            fill, stroke = GROUP[grp]
            add('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" stroke="%s" '
                'stroke-width="1.0" fill-opacity="0.85"/>' % (ax, by, bx - ax, ay - by, fill, stroke))
        pw = (amx[0] - amn[0]) * p["scale"] if name == "secA" else (bmx[0] - bmn[0]) * p["scale"]
        ph = (amx[1] - amn[1]) * p["scale"] if name == "secA" else (bmx[1] - bmn[1]) * p["scale"]
        triad(p, name, pw, ph, "X →" if name == "secA" else "Y →", "Z ↑")

    # clearance annotations (leader lines to the small gaps)
    def annotate(p, name, text, tail, head):
        (tx, ty) = px2(p, tail[0], tail[1])
        (hx, hy) = px2(p, head[0], head[1])
        add('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#37474f" stroke-width="0.8"/>'
            % (tx, ty, hx, hy))
        add('<text x="%.1f" y="%.1f" font-size="11" fill="#37474f">%s</text>'
            % (tx - 4, ty + 3, esc(text)))

    # --- legend + inventory ---
    ly = H - 330
    add('<text x="60" y="%d" font-size="14" font-weight="bold" fill="#1a1a1a">Legend</text>' % ly)
    lx = 60
    order = [("Y base", "fixed to root"), ("X base", "bolted to plate"), ("plate", "Y carriage + table"),
             ("Y assembly", "moves with y_slide"), ("X assembly", "moves with y_slide + x_slide"),
             ("Z stage", "fixed column at assembly centre · only z_slide moves")]
    for grp, note in order:
        fill, stroke = GROUP[grp]
        add('<rect x="%d" y="%d" width="18" height="12" rx="2" fill="%s" '
            'stroke="%s" stroke-width="1.3"/>' % (lx, ly + 16, fill, stroke))
        label = "%s (%s)" % (grp, note)
        add('<text x="%d" y="%d" font-size="11" fill="#37474f">%s</text>'
            % (lx + 24, ly + 26, esc(label)))
        lx += 24 + measure(label, 11) + 26
    add('<text x="60" y="%d" font-size="11" fill="#78909c">dashed red lines: section cut planes A-A / B-B</text>'
        % (ly + 44))

    # inventory table (two columns)
    add('<text x="60" y="%d" font-size="14" font-weight="bold" fill="#1a1a1a">'
        'Collision boxes by link (%d boxes / %d links)</text>' % (ly + 76, nboxes, nlinks))
    rows = []
    for name in sorted(boxes0):
        grp = group_of(name)
        short = TREE.SHORT.get(name, name)
        for (mn, mx) in boxes0[name]:
            rows.append((grp, short, (mx[0] - mn[0]) * MM, (mx[1] - mn[1]) * MM, (mx[2] - mn[2]) * MM))
    # column 0: link+dims at x[0],x[1]; column 1 at x[2],x[3]
    col_x = (60, 205, 520, 665)
    col_y0 = ly + 102
    row_h = 20
    half = (len(rows) + 1) // 2
    for c in (0, 1):
        add('<text x="%d" y="%d" font-size="10" fill="#90a4ae" font-weight="bold">LINK</text>'
            % (col_x[c * 2], col_y0 - 4))
        add('<text x="%d" y="%d" font-size="10" fill="#90a4ae" font-weight="bold">BOX (mm)</text>'
            % (col_x[c * 2 + 1], col_y0 - 4))
    for i, (grp, short, dx, dy, dz) in enumerate(rows):
        col = 0 if i < half else 1
        rr = i if col == 0 else i - half
        y = col_y0 + rr * row_h
        fill, stroke = GROUP[grp]
        add('<rect x="%d" y="%d" width="8" height="8" rx="1" fill="%s" stroke="%s" stroke-width="1"/>'
            % (col_x[col * 2] - 14, y - 7, fill, stroke))
        add('<text x="%d" y="%d" font-size="10" fill="#37474f">%s</text>'
            % (col_x[col * 2], y, esc(short)))
        add('<text x="%d" y="%d" font-size="10" fill="#546e7a">%d x %d x %d</text>'
            % (col_x[col * 2 + 1], y, round(dx), round(dy), round(dz)))
    add('<text x="60" y="%d" font-size="10" fill="#78909c">Clearances: carriage-to-wall %.1f mm (both sides), '
        'floor %.1f mm, top %.1f mm, plate-to-rail %.1f mm.</text>'
        % (col_y0 + half * row_h + 8, min(gap_side_y) * MM, gap_floor_y * MM,
           gap_top_y * MM, gap_plate * MM))
    add('</svg>')

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("\n".join(S))
    print("wrote %s (%d x %d)" % (OUT, W, H))


if __name__ == "__main__":
    main()
