#!/usr/bin/env python3
"""Generate viewer.html: self-contained three.js render of parkerstage.urdf
with all STL meshes embedded as base64 data URIs (the preview server only
serves the html file itself)."""
import base64
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
URDF = os.path.join(ROOT, "urdf", "parkerstage.urdf")
MESHES = os.path.join(ROOT, "meshes")
OUT = os.path.join(ROOT, "viewer.html")

JS = r"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>parkerstage URDF viewer</title>
<style>
  body { margin: 0; overflow: hidden; background: #1a1d23; font-family: monospace; }
  #info { position: absolute; top: 8px; left: 8px; color: #ccc; font-size: 12px; background: rgba(0,0,0,.5); padding: 6px 10px; border-radius: 4px; z-index: 10; white-space: pre; }
  #grid { position: absolute; bottom: 8px; right: 8px; color: #888; font-size: 11px; z-index: 10; }
</style>
</head>
<body>
<div id="info">loading...</div>
<div id="grid">drag to orbit · scroll to zoom</div>
<script type="importmap">
{ "imports": {
    "three": "https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js",
    "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/"
} }
</script>
<script type="module">
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { STLLoader } from 'three/addons/loaders/STLLoader.js';

const URDF_TEXT = __URDF__;
const MESH_DATA = __MESHES__;

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x1a1d23);
const camera = new THREE.PerspectiveCamera(55, innerWidth/innerHeight, 0.001, 100);
camera.position.set(0.38, 0.32, 0.55);
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(innerWidth, innerHeight);
renderer.setPixelRatio(devicePixelRatio);
document.body.appendChild(renderer.domElement);
const controls = new OrbitControls(camera, renderer.domElement);
controls.target.set(0.05, 0.12, 0.02);
controls.update();

scene.add(new THREE.HemisphereLight(0xffffff, 0x334, 1.1));
scene.add(new THREE.DirectionalLight(0xffffff, 1.4));
scene.add(new THREE.DirectionalLight(0xffffff, 0.6));
scene.add(new THREE.AxesHelper(0.05));
const g = new THREE.GridHelper(0.5, 25, 0x444, 0x333);
g.position.y = -0.05;
scene.add(g);

function mat4FromRpyXYZ(rpy, xyz) {
  const [r,p,y] = rpy, [x,yy,z] = xyz;
  const cx=Math.cos(r), sx=Math.sin(r), cy=Math.cos(p), sy=Math.sin(p), cz=Math.cos(y), sz=Math.sin(y);
  const R = [
    [cz*cy, cz*sy*sx - sz*cx, cz*sy*cx + sz*sx],
    [sz*cy, sz*sy*sx + cz*cx, sz*sy*cx - cz*sx],
    [-sy,   cy*sx,            cy*cx          ]
  ];
  const m = new THREE.Matrix4();
  m.set(R[0][0], R[1][0], R[2][0], 0,
        R[0][1], R[1][1], R[2][1], 0,
        R[0][2], R[1][2], R[2][2], 0,
        x, yy, z, 1);
  return m;
}

const stageFrame = new THREE.Group();
stageFrame.rotation.x = -Math.PI / 2; // stage Z-up -> three.js Y-up
scene.add(stageFrame);
const origin = new THREE.Object3D();
stageFrame.add(origin);

const colorFor = (function(){
  const cols = [0xcccccc, 0x66aaff, 0xff9966, 0x99ff66, 0xff6666, 0xcc66ff, 0x66ffdd, 0xffcc66, 0x99ccff, 0xcc9999, 0x99cc99, 0xcc99ff];
  let i = 0; const m = new Map();
  return (name) => { if (!m.has(name)) { m.set(name, cols[i % cols.length]); i++; } return m.get(name); };
})();

async function main() {
  const info = document.getElementById('info');
  const doc = new DOMParser().parseFromString(URDF_TEXT, 'text/xml');
  const links = {};
  for (const l of doc.querySelectorAll('link')) links[l.getAttribute('name')] = l;

  // joint list: parent -> child with origin matrix + optional prismatic axis
  const joints = [];
  const childrenOf = {};
  for (const j of doc.querySelectorAll('joint')) {
    const p = j.querySelector('parent').getAttribute('link');
    const c = j.querySelector('child').getAttribute('link');
    const o = j.querySelector('origin');
    const xyz = o ? o.getAttribute('xyz').split(' ').map(Number) : [0,0,0];
    const rpy = o ? o.getAttribute('rpy').split(' ').map(Number) : [0,0,0];
    const axEl = j.querySelector('axis');
    const axis = axEl ? axEl.getAttribute('xyz').split(' ').map(Number) : null;
    const jt = { name: j.getAttribute('name'), parent: p, child: c,
                 origin: mat4FromRpyXYZ(rpy, xyz), axis };
    joints.push(jt);
    childrenOf[p] = childrenOf[p] || [];
    childrenOf[p].push(jt);
  }

  // group per link; matrices recomputed on slider change
  const node = { root: new THREE.Group() };
  const allNames = new Set(['root']);
  for (const jt of joints) { node[jt.child] = new THREE.Group(); allNames.add(jt.child); }
  for (const n of allNames) origin.add(node[n]);

  const q = {};
  function updatePoses() {
    const pose = { root: new THREE.Matrix4().identity() };
    const order = ['root'];
    while (order.length) {
      const p = order.shift();
      for (const jt of (childrenOf[p] || [])) {
        const M = jt.origin.clone();
        if (jt.axis) {
          const v = q[jt.name] || 0;
          M.multiply(new THREE.Matrix4().makeTranslation(v * jt.axis[0], v * jt.axis[1], v * jt.axis[2]));
        }
        pose[jt.child] = new THREE.Matrix4().multiplyMatrices(pose[p], M);
        order.push(jt.child);
      }
    }
    for (const n of allNames) node[n].matrix.copy(pose[n]);
    for (const n of allNames) node[n].matrixAutoUpdate = false;
  }

  // view presets (stage frame: X right, Y into screen/left, Z up)
  const views = { iso: [[0.40, 0.34, 0.60], [0.05, 0.14, 0.04]],
                  top: [[0.05, 0.75, 0.02], [0.05, 0.14, 0.04]],
                  sideX: [[0.85, 0.18, 0.06], [0.05, 0.14, 0.04]],
                  front: [[0.05, -0.65, 0.30], [0.05, 0.14, 0.04]],
                  end: [[0.05, 0.14, 0.75], [0.05, 0.14, 0.04]] };
  const vbar = document.createElement('div');
  vbar.style.cssText = 'position:absolute;top:8px;left:8px;z-index:10;background:rgba(0,0,0,.6);padding:4px 8px;border-radius:6px;';
  for (const [k, [pos, tgt]] of Object.entries(views)) {
    const b = document.createElement('button');
    b.textContent = k; b.style.cssText = 'font:11px monospace;color:#eee;background:#333;border:1px solid #555;border-radius:4px;margin-right:4px;cursor:pointer;';
    b.onclick = () => { camera.position.set(...pos); controls.target.set(...tgt); controls.update(); };
    vbar.appendChild(b);
  }
  document.body.appendChild(vbar);

  // (view presets injected)
  const panel = document.createElement('div');
  panel.style.cssText = 'position:absolute;top:8px;right:8px;z-index:10;background:rgba(0,0,0,.6);color:#eee;font:12px monospace;padding:8px 12px;border-radius:6px;';
  const labelFor = { y_slide: 'Y slide (bottom stage)', x_slide: 'X slide (middle stage)', z_slide: 'Z slide (vertical stage)' };
  for (const jt of joints) {
    if (!jt.axis) continue;
    const lim = [...doc.querySelectorAll('joint')].find(j => j.getAttribute('name') === jt.name).querySelector('limit');
    const lo = lim ? parseFloat(lim.getAttribute('lower')) : -0.1;
    const hi = lim ? parseFloat(lim.getAttribute('upper')) : 0.1;
    q[jt.name] = 0;
    const row = document.createElement('div');
    const lab = document.createElement('div'); lab.textContent = (labelFor[jt.name] || jt.name) + ': 0.000 m';
    const inp = document.createElement('input');
    inp.type = 'range'; inp.min = lo; inp.max = hi; inp.step = 0.001; inp.value = 0;
    inp.style.width = '170px';
    inp.oninput = () => { q[jt.name] = parseFloat(inp.value); lab.textContent = (labelFor[jt.name] || jt.name) + ': ' + parseFloat(inp.value).toFixed(3) + ' m'; updatePoses(); };
    row.appendChild(lab); row.appendChild(document.createElement('br')); row.appendChild(inp);
    panel.appendChild(row);
  }
  document.body.appendChild(panel);

  const loader = new STLLoader();
  const cache = {};
  function getMesh(filename) {
    const name = filename.split('/').pop();
    if (!cache[name]) cache[name] = loader.loadAsync(MESH_DATA[name]);
    return cache[name];
  }

  function makeLabel(text, color) {
    const cv = document.createElement('canvas');
    cv.width = 512; cv.height = 128;
    const ctx = cv.getContext('2d');
    ctx.fillStyle = 'rgba(0,0,0,0.65)';
    ctx.fillRect(0, 0, 512, 128);
    ctx.font = 'bold 44px monospace';
    ctx.fillStyle = color;
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    ctx.fillText(text, 256, 64);
    const tex = new THREE.CanvasTexture(cv);
    const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, depthTest: false }));
    sp.scale.set(0.09, 0.0225, 1);
    return sp;
  }

  const shortName = (n) => n
    .replace(/__401XR___[A-Za-z0-9_]+$/, '')
    .replace(/^401xr___/, '')
    .replace(/^401200xr__/, 'BASE ');

  let done = 0, total = 0;
  let li = 0;
  const pose = {};
  {
    const o2 = ['root'];
    pose.root = new THREE.Matrix4().identity();
    while (o2.length) {
      const p = o2.shift();
      for (const jt of (childrenOf[p] || [])) {
        pose[jt.child] = new THREE.Matrix4().multiplyMatrices(pose[p], jt.origin);
        o2.push(jt.child);
      }
    }
  }
  const filter = location.hash ? location.hash.slice(1).split(',').filter(Boolean) : null;
  const matchLink = (name, l) => !filter || filter.some(f =>
    name.toLowerCase().includes(f.toLowerCase()) ||
    [...l.querySelectorAll('mesh')].some(m => m.getAttribute('filename').toLowerCase().includes(f.toLowerCase())));
  for (const name in pose) {
    const l = links[name];
    if (!l) continue;
    if (!matchLink(name, l)) continue;
    const geoms = l.querySelectorAll('visual geometry mesh');
    total += geoms.length;
    for (const msh of geoms) {
      const geom = await getMesh(msh.getAttribute('filename'));
      const v = msh.closest('visual');
      const o = v.querySelector('origin');
      const xyz = o ? o.getAttribute('xyz').split(' ').map(Number) : [0,0,0];
      const rpy = o ? o.getAttribute('rpy').split(' ').map(Number) : [0,0,0];
      const mesh = new THREE.Mesh(geom, new THREE.MeshLambertMaterial({ color: colorFor(name) }));
      mesh.applyMatrix4(mat4FromRpyXYZ(rpy, xyz));
      node[name].add(mesh);
      // collision boxes (semi-transparent wireframe overlay)
      for (const col of l.querySelectorAll('collision')) {
        const co = col.querySelector('origin');
        const b = col.querySelector('geometry box');
        if (!b) continue;
        const cxyz = co ? co.getAttribute('xyz').split(' ').map(Number) : [0,0,0];
        const crpy = co ? co.getAttribute('rpy').split(' ').map(Number) : [0,0,0];
        const sz = b.getAttribute('size').split(' ').map(Number);
        const box = new THREE.Mesh(new THREE.BoxGeometry(sz[0], sz[1], sz[2]), new THREE.MeshBasicMaterial({
          color: 0xffcc00, wireframe: true, transparent: true, opacity: 0.35, depthTest: false }));
        box.applyMatrix4(mat4FromRpyXYZ(crpy, cxyz));
        node[name].add(box);
      }
      const structural = name.startsWith('401200') || name.includes('carriage') || name === 'plate';
      if (geoms.length && l.querySelectorAll('inertial').length && structural) {
        const lab = makeLabel(shortName(name), '#' + new THREE.Color(colorFor(name)).getHexString());
        const p = new THREE.Vector3().setFromMatrixPosition(pose[name]);
        lab.position.set(p.x + ((li % 3) - 1) * 0.05, p.y + 0.04, p.z - ((li % 3) - 1) * 0.05);
        node[name].add(lab);
      }
      done++;
      info.textContent = `loaded ${done}/${total}`;
    }
    li++;
  }
  updatePoses();
  info.textContent = `loaded ${done} meshes · drag sliders to slide · orbit to inspect`;
  animate();
}

function animate() {
  requestAnimationFrame(animate);
  controls.update();
  renderer.render(scene, camera);
}
main();
</script>
</body>
</html>
"""


def main():
    with open(URDF, encoding="utf-8") as f:
        urdf = f.read()

    meshes = {}
    for fn in sorted(os.listdir(MESHES)):
        if fn.endswith(".stl"):
            with open(os.path.join(MESHES, fn), "rb") as f:
                meshes[fn] = "data:application/octet-stream;base64," + base64.b64encode(f.read()).decode()

    js = JS.replace("__URDF__", json_dumps(urdf)).replace("__MESHES__", json_dumps(meshes))
    with open(OUT, "w", encoding="utf-8") as f:
        f.write(js)
    print(f"wrote {OUT} ({os.path.getsize(OUT) / 1e6:.1f} MB)")


def json_dumps(obj):
    import json
    return json.dumps(obj)


if __name__ == "__main__":
    main()
