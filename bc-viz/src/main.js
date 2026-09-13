import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import "./style.css";

const canvas = document.getElementById("stage");
const metaLine = document.getElementById("meta-line");
const statsEl = document.getElementById("stats");
const condEl = document.getElementById("conditions");
const scrub = document.getElementById("scrub");
const tLabel = document.getElementById("t-label");
const playBtn = document.getElementById("play");
const showEdges = document.getElementById("show-edges");
const dimIdle = document.getElementById("dim-idle");
const autoRotate = document.getElementById("auto-rotate");

const COLORS = {
  void: new THREE.Color("#070a10"),
  other: new THREE.Color("#3d4a5c"),
  otherDim: new THREE.Color("#1a2230"),
  bc: new THREE.Color("#5eebff"),
  p1: new THREE.Color("#ffd166"),
  mal: new THREE.Color("#ff6b9d"),
  lg6: new THREE.Color("#8b9dff"),
  lg5: new THREE.Color("#c9a0ff"),
  vab3: new THREE.Color("#7cffb2"),
  active: new THREE.Color("#7cffb2"),
};

const ROLE_COLOR = {
  bc: COLORS.bc,
  p1: COLORS.p1,
  mal: COLORS.mal,
  lg6: COLORS.lg6,
  lg5: COLORS.lg5,
  vab3: COLORS.vab3,
  other: COLORS.other,
};

let data = null;
let condition = "bc_exc";
let binIndex = 5;
let playing = true;
const reducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
if (reducedMotion) playing = false;

// --- three setup ---
const renderer = new THREE.WebGLRenderer({
  canvas,
  antialias: true,
  alpha: false,
  powerPreference: "high-performance",
});
renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.setClearColor(COLORS.void, 1);

const scene = new THREE.Scene();
scene.fog = new THREE.FogExp2(COLORS.void.getHex(), 0.0045);

const camera = new THREE.PerspectiveCamera(
  42,
  window.innerWidth / window.innerHeight,
  0.1,
  2000,
);
camera.position.set(140, 90, 220);

const controls = new OrbitControls(camera, canvas);
controls.enableDamping = true;
controls.dampingFactor = 0.06;
controls.minDistance = 40;
controls.maxDistance = 600;
controls.autoRotate = !reducedMotion;
controls.autoRotateSpeed = 0.35;

const root = new THREE.Group();
scene.add(root);

// soft ambient + key light for depth cues on points (points ignore lights; keep for future meshes)
scene.add(new THREE.AmbientLight(0x6aa0c0, 0.35));
const key = new THREE.DirectionalLight(0xb8d4ff, 0.55);
key.position.set(80, 120, 40);
scene.add(key);

// --- state buffers ---
let points = null;
let geometry = null;
let material = null;
let baseColors = null;
let liveColors = null;
let activityStrength = null;
let nVis = 0;
let edgeLines = null;
let stimMarkers = null;
let p1Markers = null;

function roleColor(role) {
  return ROLE_COLOR[role] || COLORS.other;
}

function buildScene(payload) {
  data = payload;
  nVis = payload.meta.n_vis;
  const posArr = new Float32Array(payload.positions);
  baseColors = new Float32Array(nVis * 3);
  liveColors = new Float32Array(nVis * 3);
  activityStrength = new Float32Array(nVis);

  for (let i = 0; i < nVis; i++) {
    const c = roleColor(payload.roles[i]);
    baseColors[i * 3] = c.r;
    baseColors[i * 3 + 1] = c.g;
    baseColors[i * 3 + 2] = c.b;
    liveColors[i * 3] = c.r;
    liveColors[i * 3 + 1] = c.g;
    liveColors[i * 3 + 2] = c.b;
  }

  geometry = new THREE.BufferGeometry();
  geometry.setAttribute("position", new THREE.BufferAttribute(posArr, 3));
  geometry.setAttribute("color", new THREE.BufferAttribute(liveColors, 3));

  material = new THREE.PointsMaterial({
    size: 1.15,
    vertexColors: true,
    sizeAttenuation: true,
    transparent: true,
    opacity: 0.92,
    depthWrite: false,
    blending: THREE.AdditiveBlending,
  });

  points = new THREE.Points(geometry, material);
  root.add(points);

  // b/c markers — larger, always bright
  const bcPos = [];
  for (const vi of payload.bc_vis) {
    bcPos.push(posArr[vi * 3], posArr[vi * 3 + 1], posArr[vi * 3 + 2]);
  }
  const bcGeo = new THREE.BufferGeometry();
  bcGeo.setAttribute("position", new THREE.Float32BufferAttribute(bcPos, 3));
  stimMarkers = new THREE.Points(
    bcGeo,
    new THREE.PointsMaterial({
      color: COLORS.bc,
      size: 6,
      sizeAttenuation: true,
      transparent: true,
      opacity: 1,
      depthWrite: false,
      blending: THREE.AdditiveBlending,
    }),
  );
  root.add(stimMarkers);

  const p1Pos = [];
  for (const vi of payload.p1_vis) {
    p1Pos.push(posArr[vi * 3], posArr[vi * 3 + 1], posArr[vi * 3 + 2]);
  }
  const p1Geo = new THREE.BufferGeometry();
  p1Geo.setAttribute("position", new THREE.Float32BufferAttribute(p1Pos, 3));
  p1Markers = new THREE.Points(
    p1Geo,
    new THREE.PointsMaterial({
      color: COLORS.p1,
      size: 4.5,
      sizeAttenuation: true,
      transparent: true,
      opacity: 0.95,
      depthWrite: false,
      blending: THREE.AdditiveBlending,
    }),
  );
  root.add(p1Markers);

  // edge rays from each b/c cell toward top targets
  const edgePos = [];
  const edgeCol = [];
  const bcCentroid = new THREE.Vector3();
  for (const vi of payload.bc_vis) {
    bcCentroid.x += posArr[vi * 3];
    bcCentroid.y += posArr[vi * 3 + 1];
    bcCentroid.z += posArr[vi * 3 + 2];
  }
  bcCentroid.divideScalar(Math.max(payload.bc_vis.length, 1));

  for (const e of payload.edges) {
    const t = e.target_vis;
    edgePos.push(
      bcCentroid.x,
      bcCentroid.y,
      bcCentroid.z,
      posArr[t * 3],
      posArr[t * 3 + 1],
      posArr[t * 3 + 2],
    );
    const w = Math.min(1, e.weight / 400);
    edgeCol.push(0.2, 0.85 + 0.15 * w, 1.0, 0.15, 0.55, 0.7);
  }
  const eGeo = new THREE.BufferGeometry();
  eGeo.setAttribute("position", new THREE.Float32BufferAttribute(edgePos, 3));
  eGeo.setAttribute("color", new THREE.Float32BufferAttribute(edgeCol, 3));
  edgeLines = new THREE.LineSegments(
    eGeo,
    new THREE.LineBasicMaterial({
      vertexColors: true,
      transparent: true,
      opacity: 0.22,
      blending: THREE.AdditiveBlending,
      depthWrite: false,
    }),
  );
  root.add(edgeLines);

  scrub.max = String(payload.activity[condition].length - 1);
  applyFrame(binIndex, true);
  updateMeta();
}

function currentBins() {
  return data.activity[condition];
}

function applyFrame(bin, force = false) {
  if (!data) return;
  const bins = currentBins();
  bin = Math.max(0, Math.min(bins.length - 1, bin));
  binIndex = bin;
  scrub.value = String(bin);
  tLabel.textContent = `${bin * data.meta.bin_ms} ms`;

  activityStrength.fill(0);
  const pairs = bins[bin];
  for (const [vi, cnt] of pairs) {
    activityStrength[vi] = Math.min(1.5, activityStrength[vi] + cnt * 0.85);
  }
  // temporal bleed: previous bin ghost
  if (bin > 0) {
    for (const [vi, cnt] of bins[bin - 1]) {
      activityStrength[vi] = Math.max(activityStrength[vi], Math.min(0.55, cnt * 0.35));
    }
  }

  const dim = dimIdle.checked;
  for (let i = 0; i < nVis; i++) {
    const s = activityStrength[i];
    const role = data.roles[i];
    const base = roleColor(role);
    let r, g, b;
    if (s > 0.05) {
      const k = Math.min(1, s);
      r = base.r + (COLORS.active.r - base.r) * k;
      g = base.g + (COLORS.active.g - base.g) * k;
      b = base.b + (COLORS.active.b - base.b) * k;
      const boost = 1 + s * 0.8;
      r *= boost;
      g *= boost;
      b *= boost;
    } else if (dim && role === "other") {
      r = COLORS.otherDim.r;
      g = COLORS.otherDim.g;
      b = COLORS.otherDim.b;
    } else {
      r = base.r;
      g = base.g;
      b = base.b;
    }
    liveColors[i * 3] = r;
    liveColors[i * 3 + 1] = g;
    liveColors[i * 3 + 2] = b;
  }
  geometry.attributes.color.needsUpdate = true;

  if (edgeLines) edgeLines.visible = showEdges.checked;
  if (force) updateStats();
}

function activeCount(bin) {
  return currentBins()[bin]?.length || 0;
}

function updateStats() {
  if (!data) return;
  const m = data.meta;
  const cond = m.conditions[condition];
  statsEl.textContent =
    `${cond.label}  ·  t=${binIndex * m.bin_ms}ms  ·  active ${activeCount(binIndex)}  ·  ` +
    `P1 ${m.P1_spikes[condition]}  ·  net ${m.network_spikes[condition]}  ·  b/c ${m.bc_spikes[condition]}`;
}

function updateMeta() {
  if (!data) return;
  const m = data.meta;
  metaLine.innerHTML =
    `${m.n_vis.toLocaleString()} cells with coordinates<br/>` +
    `seed ${m.seed} · ${m.poisson_hz} Hz on b/c · ${m.duration_ms} ms`;
}

function buildConditionUI() {
  condEl.innerHTML = "";
  for (const [key, c] of Object.entries(data.meta.conditions)) {
    const btn = document.createElement("button");
    btn.type = "button";
    btn.textContent = c.label;
    btn.setAttribute("role", "tab");
    btn.setAttribute("aria-selected", key === condition ? "true" : "false");
    btn.addEventListener("click", () => {
      condition = key;
      for (const el of condEl.querySelectorAll("button")) {
        el.setAttribute("aria-selected", "false");
      }
      btn.setAttribute("aria-selected", "true");
      scrub.max = String(currentBins().length - 1);
      applyFrame(binIndex, true);
      updateStats();
    });
    condEl.appendChild(btn);
  }
}

playBtn.addEventListener("click", () => {
  playing = !playing;
  playBtn.textContent = playing ? "Pause" : "Play";
  playBtn.setAttribute("aria-pressed", playing ? "true" : "false");
});

scrub.addEventListener("input", () => {
  playing = false;
  playBtn.textContent = "Play";
  playBtn.setAttribute("aria-pressed", "false");
  applyFrame(Number(scrub.value), true);
});

showEdges.addEventListener("change", () => applyFrame(binIndex, true));
dimIdle.addEventListener("change", () => applyFrame(binIndex, true));
autoRotate.addEventListener("change", () => {
  controls.autoRotate = autoRotate.checked && !reducedMotion;
});

window.addEventListener("resize", () => {
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
});

async function boot() {
  const res = await fetch("./data/trial.json");
  if (!res.ok) {
    metaLine.textContent = "failed to load trial.json";
    throw new Error("trial.json missing");
  }
  data = await res.json();
  buildConditionUI();
  buildScene(data);
  if (reducedMotion) {
    playBtn.textContent = "Play";
    playBtn.setAttribute("aria-pressed", "false");
  }

  const clock = new THREE.Clock();
  let acc = 0;
  const BIN_SEC = data.meta.bin_ms / 1000;

  function tick() {
    requestAnimationFrame(tick);
    const dt = clock.getDelta();
    if (playing && data) {
      acc += dt;
      if (acc >= BIN_SEC) {
        acc = 0;
        let next = binIndex + 1;
        if (next >= currentBins().length) next = 5;
        applyFrame(next, true);
        updateStats();
      }
    }
    if (stimMarkers) {
      const t = clock.elapsedTime;
      stimMarkers.material.size = 5.2 + Math.sin(t * 6) * 1.4;
    }
    if (edgeLines && showEdges.checked) {
      edgeLines.material.opacity = 0.12 + 0.12 * (0.5 + 0.5 * Math.sin(clock.elapsedTime * 3));
    }
    controls.update();
    renderer.render(scene, camera);
  }
  tick();
}

boot().catch((err) => {
  console.error(err);
  metaLine.textContent = String(err.message || err);
});
