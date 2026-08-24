#!/usr/bin/env python3
"""Generate docs/kinematic_tree.svg from urdf/parkerstage.urdf.

Draws the kinematic skeleton of the 401200XR compound XY stage as a top-down
SVG tree: links are boxes, joints are labeled edges.  Chains of fixed joints
are grouped into assembly boxes (every link is still listed, in tree order) so
the two prismatic slide joints -- the only movable DOF -- stand out.

Run:  python3 tools/make_tree_svg.py
"""
import os
import re
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
URDF = os.path.join(ROOT, "urdf", "parkerstage.urdf")
OUT = os.path.join(ROOT, "docs", "kinematic_tree.svg")
REBUILD = os.path.join(ROOT, "tools", "rebuild_urdf.py")

FONT = "DejaVu Sans Mono, Menlo, Consolas, monospace"

# Onshape part names -> readable labels
SHORT = {
    "401xr___encoder__401xr___encoder": "encoder",
    "401xr___encoder__401xr___encoder_1": "encoder",
    "401xr___carriage__401xr___carriage": "carriage",
    "401xr___carriage__401xr___carriage_1": "carriage",
    "401xr___carriage_end_caps__401xr___carriage_end_caps_1": "end caps",
    "401xr___carriage_end_caps__401xr___carriage_end_caps_1_1": "end caps",
    "401xr___carriage_end_caps__401xr___carriage_end_caps_2": "end caps",
    "401xr___encoder_base_2__401xr___encoder_base_2": "encoder base",
    "401xr___encoder_base_2__401xr___encoder_base_2_1": "encoder base",
    "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch": "home-limit switch",
    "401xr___h2__l1___homelimit_switch__401xr___h2__l1___homelimit_switch_1": "home-limit switch",
    "401xr___switch_flag__401xr___switch_flag": "switch flag",
    "401xr___switch_flag__401xr___switch_flag_1": "switch flag",
    "401200xr__1_": "401200xr__1_ (Y base)",
    "401200xr__3_": "401200xr__3_ (X base)",
    "base1_fixed_to_root": "base1_fixed_to_root",
    "base2_mounted_to_plate": "base2_mounted_to_plate",
}

# box id -> (title, fill, stroke)
BOX_STYLE = {
    "root":           ("root (world)",              "#e8f5e9", "#2e7d32"),
    "Y base":         ("401200xr__1_  ·  Y base",   "#eceff1", "#546e7a"),
    "X base":         ("401200xr__3_  ·  X base",   "#eceff1", "#546e7a"),
    "plate":          ("plate  ·  Y carriage + table", "#fafafa", "#9e9e9e"),
    "X assembly":     ("X carriage assembly",       "#e3f2fd", "#1e88e5"),
    "Y assembly":     ("Y carriage assembly",       "#e3f2fd", "#1e88e5"),
}

GAP_COL = 140          # horizontal gap between columns
MARGIN = 48            # outer margin
VGAP = 96              # vertical gap between parent bottom and child top
FS_TITLE, FS_MEMBER, FS_EDGE, FS_META = 13, 10, 10, 11
PAD_X, PAD_Y, LINE_H = 14, 10, 14


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def measure(s, size):
    return len(s) * 0.6 * size


def wrap(s, size, max_w):
    lines, cur = [], ""
    for w in s.split():
        t = (cur + " " + w).strip()
        if measure(t, size) <= max_w or not cur:
            cur = t
        else:
            lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)
    return lines


def preorder(children, link):
    """Links under `link` (inclusive), depth-first pre-order."""
    out = [link]
    for c, _j in children.get(link, []):
        out.extend(preorder(children, c))
    return out


def main():
    tree = ET.parse(URDF)
    root_el = tree.getroot()
    links = {l.get("name") for l in root_el.findall("link")}
    joints = {j.get("name"): j for j in root_el.findall("joint")}

    children = {}
    for name, j in joints.items():
        children.setdefault(j.find("parent").get("link"), []).append(
            (j.find("child").get("link"), name))

    # --- collapse fixed-joint chains into assembly boxes -------------------
    x_slide = joints["x_slide"]
    y_slide = joints["y_slide"]
    fastened = joints["fastened_1"]
    x_root = x_slide.find("child").get("link")
    y_root = fastened.find("child").get("link")
    x_tree = preorder(children, x_root)
    y_tree = preorder(children, y_root)
    assert not (set(x_tree) & set(y_tree)), "X/Y assembly subtrees overlap"

    box_of = {}
    for l in x_tree:
        box_of[l] = "X assembly"
    for l in y_tree:
        box_of[l] = "Y assembly"
    for l in ("plate",):
        box_of[l] = "plate"
    for l in ("401200xr__1_",):
        box_of[l] = "Y base"
    for l in ("401200xr__3_",):
        box_of[l] = "X base"
    box_of["root"] = "root"
    assert set(box_of) == links, "unassigned links: %s" % (links - set(box_of))

    # box member lists (tree order); representative first
    def members(box):
        mem = [l for l, b in box_of.items() if b == box]
        return mem
    order = {
        "root": [("root", [])],
    }
    rep = {
        "Y base": "401200xr__1_",
        "X base": "401200xr__3_",
        "plate": "plate",
        "X assembly": x_root,
        "Y assembly": y_root,
    }
    box_members = {}
    for box, r in rep.items():
        mem = members(box)
        mem.remove(r)
        # keep tree order for the rest (preorder already gives it)
        rank = {l: i for i, l in enumerate(preorder(children, "root"))}
        mem.sort(key=lambda l: rank[l])
        box_members[box] = (r, mem)
    box_members["root"] = ("root", [])

    # --- collapsed tree edges (joints crossing boxes) ----------------------
    edge_order = ["base1_fixed_to_root", "y_slide", "base2_mounted_to_plate",
                  "x_slide", "fastened_1"]
    edges = []
    for name in edge_order:
        j = joints[name]
        p = box_of[j.find("parent").get("link")]
        c = box_of[j.find("child").get("link")]
        if p != c:
            edges.append((p, c, name))

    # --- geometry ----------------------------------------------------------
    depth = {"root": 0, "Y base": 1, "plate": 2, "X base": 3,
             "Y assembly": 3, "X assembly": 4}
    column = {"root": 0, "Y base": 1, "plate": 2, "X base": 3,
              "X assembly": 4, "Y assembly": 5}
    boxes = {}
    for b, (title, fill, stroke) in BOX_STYLE.items():
        r, mem = box_members[b]
        title = r if r == b or b == "root" else title
        title_lines = wrap(title, FS_TITLE, 400)
        h = PAD_Y * 2 + len(title_lines) * (FS_TITLE + 4)
        if mem:
            h += LINE_H * len(mem)
        w = max([measure(t, FS_TITLE) for t in title_lines] +
                [measure(SHORT.get(m, m), FS_MEMBER) for m in mem] +
                [measure("· " + SHORT.get(m, m), FS_MEMBER) for m in mem]) \
            + 2 * PAD_X
        boxes[b] = dict(fill=fill, stroke=stroke, w=w, h=h, title=title_lines,
                        mem=mem, rep=r)

    max_h = max(b["h"] for b in boxes.values())
    row_h = max_h + VGAP
    y_of = lambda b: MARGIN + depth[b] * row_h
    x_of = {}
    x = MARGIN
    for b in sorted(column, key=lambda k: column[k]):
        x_of[b] = x + boxes[b]["w"] / 2
        x += boxes[b]["w"] + GAP_COL
    W = x - GAP_COL + MARGIN
    H = MARGIN + 4 * row_h + max_h + 150  # + legend block

    # --- joint labels -------------------------------------------------------
    def joint_label(name):
        j = joints[name]
        if j.get("type") == "prismatic":
            axis = j.find("axis").get("xyz")
            lim = j.find("limit")
            world = "world +Y" if name == "y_slide" else "world +X"
            return True, ["%s  [prismatic]" % name,
                          "axis %s  ->  %s" % (axis, world),
                          "q in [%s, %s] m" % (lim.get("lower"), lim.get("upper"))]
        return False, ["%s  [fixed]" % name]

    # home offsets (single source of truth: rebuild_urdf.py)
    src = open(REBUILD, encoding="utf-8").read()
    home_y = re.search(r"HOME_Y = ([-\d.]+)", src).group(1)
    home_x = re.search(r"HOME_X = ([-\d.]+)", src).group(1)

    # --- emit ----------------------------------------------------------------
    S = []
    add = S.append
    add('<?xml version="1.0" encoding="UTF-8"?>')
    add('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
        'viewBox="0 0 %d %d" font-family="%s">' % (W, H, W, H, FONT))
    add("""<defs>
  <marker id="arrow" markerWidth="9" markerHeight="9" refX="8" refY="4" orient="auto">
    <path d="M0 0 L8 4 L0 8 z" fill="#1e88e5"/>
  </marker>
</defs>""")
    add('<rect x="0" y="0" width="%d" height="%d" fill="#ffffff"/>' % (W, H))
    add('<text x="%d" y="%d" font-size="17" font-weight="bold" fill="#1a1a1a">'
        % (MARGIN, MARGIN - 22))
    add("parkerstage  ·  401200XR compound XY stage  ·  kinematic tree</text>")
    add('<text x="%d" y="%d" font-size="%d" fill="#546e7a">' % (MARGIN, MARGIN - 4, FS_META))
    add("Two 401XR stages stacked: the bottom slides along world +Y, the top along world +X "
        "(200 mm stroke, centered on mid-stroke).  "
        "17 links / 16 joints / 2 prismatic DOF.  Home: y_slide +%s mm, x_slide %s mm."
        % (home_y, home_x))
    add("</text>")

    # edges first (under boxes)
    for (p, c, jname) in edges:
        p_ = boxes[p]
        c_ = boxes[c]
        x1, y1 = x_of[p], y_of(p) + p_["h"]
        x2, y2 = x_of[c], y_of(c)
        mid = (y1 + y2) / 2
        prism, label = joint_label(jname)
        if prism:
            add('<path d="M %g %g V %g H %g V %g" fill="none" stroke="#1e88e5" '
                'stroke-width="2.5" marker-end="url(#arrow)"/>' % (x1, y1, mid, x2, y2))
        else:
            add('<path d="M %g %g V %g H %g V %g" fill="none" stroke="#90a4ae" '
                'stroke-width="1.5" stroke-dasharray="6 4"/>' % (x1, y1, mid, x2, y2))
        ly = mid - 8 - (len(label) - 1) * 11
        col = "#1565c0" if prism else "#607d8b"
        weight = "bold" if prism else "normal"
        add('<text x="%g" y="%g" font-size="%d" fill="%s" font-weight="%s" text-anchor="middle">'
            % (x2, ly, FS_EDGE, col, weight))
        for k, line in enumerate(label):
            add('<tspan x="%g" dy="%s">%s</tspan>' % (x2, "0" if k == 0 else 11, esc(line)))
        add("</text>")

    # boxes
    for b in BOX_STYLE:
        bb = boxes[b]
        x = x_of[b] - bb["w"] / 2
        y = y_of(b)
        add('<rect x="%g" y="%g" width="%g" height="%g" rx="7" fill="%s" '
            'stroke="%s" stroke-width="1.8"/>' % (x, y, bb["w"], bb["h"], bb["fill"], bb["stroke"]))
        ty = y + PAD_Y + FS_TITLE
        for line in bb["title"]:
            add('<text x="%g" y="%g" font-size="%d" font-weight="bold" fill="#1a1a1a">%s</text>'
                % (x + PAD_X, ty, FS_TITLE, esc(line)))
            ty += FS_TITLE + 4
        for m in bb["mem"]:
            add('<text x="%g" y="%g" font-size="%d" fill="#455a64">%s</text>'
                % (x + PAD_X, ty, FS_MEMBER, esc("· " + SHORT.get(m, m))))
            ty += LINE_H

    # legend
    ly = MARGIN + 4 * row_h + max_h + 28
    items = [("#e8f5e9", "#2e7d32", "root (world)"),
             ("#eceff1", "#546e7a", "fixed base"),
             ("#fafafa", "#9e9e9e", "plate / table"),
             ("#e3f2fd", "#1e88e5", "moving assembly"),
             (None, "#1e88e5", "prismatic joint"),
             (None, "#90a4ae", "fixed joint")]
    lx = MARGIN
    for fill, stroke, label in items:
        if fill:
            add('<rect x="%d" y="%d" width="18" height="12" rx="2" fill="%s" '
                'stroke="%s" stroke-width="1.5"/>' % (lx, ly, fill, stroke))
        else:
            dash = ' stroke-dasharray="4 3"' if stroke == "#90a4ae" else ""
            add('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" '
                'stroke-width="2.5"%s/>' % (lx, ly + 6, lx + 18, ly + 6, stroke, dash))
        add('<text x="%d" y="%d" font-size="%d" fill="#455a64">%s</text>'
            % (lx + 24, ly + 11, FS_EDGE, esc(label)))
        lx += 24 + measure(label, FS_EDGE) + 26
    add('<text x="%d" y="%d" font-size="%d" fill="#78909c">Fixed-joint chains are grouped '
        'into assembly boxes; every link is listed in tree order.</text>'
        % (MARGIN, ly + 34, FS_EDGE))
    add("</svg>")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("\n".join(S))
    print("wrote %s (%d x %d)" % (OUT, W, H))


if __name__ == "__main__":
    main()
