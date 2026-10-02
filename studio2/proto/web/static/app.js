// Prototype C (spec 007, jalon J1) : écran Résultats d'un triaxial 2D.
//
// Principe de performance : TOUTES les frames sont chargées une fois dans des
// textures GPU (positions, champs, modes de joints). Changer de frame revient à
// changer un entier uniforme ; aucun tampon n'est recopié. Le CPU ne touche jamais
// un élément individuel après le chargement.

const LARG = 4096;                       // largeur des textures de données
const T0 = performance.now();
const MESURE = new URLSearchParams(location.search).has("mesure");

const CHAMPS = [
  { id: "sigmaYY", nom: "σyy", unite: "MPa", k: 1e-6 },
  { id: "sigmaXX", nom: "σxx", unite: "MPa", k: 1e-6 },
  { id: "sigmaXY", nom: "σxy", unite: "MPa", k: 1e-6 },
  { id: "vonMises", nom: "von Mises", unite: "MPa", k: 1e-6 },
  { id: "epsXX", nom: "εxx", unite: "%", k: 100 },
  { id: "phase", nom: "Phases", unite: "", k: 1, cat: true },
];
const BANDES = 12;                        // bandes de couleur, comme Abaqus
const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const hex = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16) / 255);

// ---------------------------------------------------------------- état
const E = {
  meta: null, champ: CHAMPS[0], frame: 0, masque: 1 | 2 | 4,
  cx: 0, cy: 0, k: 1, sale: true, lecture: false, hist: null, comptes: null,
};

// ---------------------------------------------------------------- WebGL
const cv = document.getElementById("gl");
const gl = cv.getContext("webgl2", { antialias: true, alpha: true });
if (!gl) document.body.textContent = "WebGL2 indisponible.";

function compiler(vs, fs) {
  const p = gl.createProgram();
  for (const [type, src] of [[gl.VERTEX_SHADER, vs], [gl.FRAGMENT_SHADER, fs]]) {
    const s = gl.createShader(type);
    gl.shaderSource(s, src); gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
    gl.attachShader(p, s);
  }
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
  const u = {};
  for (let i = 0; i < gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS); i++) {
    const n = gl.getActiveUniform(p, i).name;
    u[n] = gl.getUniformLocation(p, n);
  }
  return { p, u };
}

const ENTETE = `#version 300 es
precision highp float; precision highp int; precision highp sampler2D; precision highp usampler2D;
uniform sampler2D uPos; uniform int uFrame, uNVert; uniform vec2 uCentre, uEchelle;
ivec2 tc(int i) { return ivec2(i % ${LARG}, i / ${LARG}); }
vec2 pos(int v) { return texelFetch(uPos, tc(uFrame * uNVert + v), 0).rg; }
vec2 ndc(vec2 p) { return (p - uCentre) * uEchelle; }
`;

// Triangles : sommet k -> triangle k/3 (le solveur écrit des nœuds propres à chaque élément).
const progTri = compiler(ENTETE + `
uniform sampler2D uVal; uniform int uFrameVal, uNTri;
flat out float vVal;
void main() {
  int v = gl_VertexID;
  vVal = texelFetch(uVal, tc(uFrameVal * uNTri + v / 3), 0).r;
  gl_Position = vec4(ndc(pos(v)), 0.0, 1.0);
}`, `#version 300 es
precision highp float;
flat in float vVal; uniform float uMin, uMax; uniform int uCat, uBandes; uniform vec3 uCat3[6];
out vec4 o;
vec3 arcenciel(float t) {          // bleu -> cyan -> vert -> jaune -> rouge
  vec3 c[5] = vec3[5](vec3(.16,.27,.82), vec3(.12,.72,.86), vec3(.33,.79,.35), vec3(.96,.82,.24), vec3(.85,.21,.19));
  float x = clamp(t, 0., 1.) * 4.; int i = min(int(x), 3);
  return mix(c[i], c[i + 1], x - float(i));
}
void main() {
  if (uCat == 1) { o = vec4(uCat3[int(vVal) % 6], 1.); return; }
  float t = (vVal - uMin) / (uMax - uMin);
  t = (floor(clamp(t, 0., .9999) * float(uBandes)) + .5) / float(uBandes);
  o = vec4(arcenciel(t), 1.);
}`);

// Joints : un quadrilatère de largeur fixe à l'écran par segment, dessiné par instanciation.
const progJoint = compiler(ENTETE + `
layout(location = 0) in uvec2 aSeg;
uniform usampler2D uMode; uniform int uNJoint, uMasque; uniform vec2 uPixel; uniform float uDemiLarg;
uniform vec3 uCoul[3];
out vec3 vCoul;
void main() {
  int m = int(texelFetch(uMode, tc(uFrame * uNJoint + gl_InstanceID), 0).r);
  if (m == 0 || (m & uMasque) == 0) { gl_Position = vec4(2., 2., 2., 1.); return; }
  vCoul = uCoul[m == 1 ? 0 : (m == 2 ? 1 : 2)];
  int c = gl_VertexID;                                  // 0..5 : deux triangles
  float bout = (c == 1 || c == 2 || c == 4) ? 1. : 0.;
  float cote = (c == 2 || c == 4 || c == 5) ? 1. : -1.;
  vec2 a = ndc(pos(int(aSeg.x))), b = ndc(pos(int(aSeg.y)));
  vec2 d = (b - a) / uPixel; float l = length(d);
  vec2 n = l > 0. ? vec2(-d.y, d.x) / l : vec2(0.);
  gl_Position = vec4(mix(a, b, bout) + n * cote * uDemiLarg * uPixel, 0., 1.);
}`, `#version 300 es
precision highp float; in vec3 vCoul; out vec4 o;
void main() { o = vec4(vCoul, 1.); }`);

const T_SHADERS = performance.now();

function texture(interne, format, type, donnees, n) {
  const h = Math.ceil(n / LARG), comp = donnees.length / n;
  const plein = new donnees.constructor(LARG * h * comp);
  plein.set(donnees);
  const t = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, t);
  gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
  gl.texImage2D(gl.TEXTURE_2D, 0, interne, LARG, h, 0, format, type, plein);
  for (const p of [gl.TEXTURE_MIN_FILTER, gl.TEXTURE_MAG_FILTER]) gl.texParameteri(gl.TEXTURE_2D, p, gl.NEAREST);
  return t;
}

const GPU = { pos: null, val: {}, mode: null, vao: null, nSeg: 0 };

async function bin(run, nom, Type) {
  const r = await fetch(`/cache/${run}/${nom}.bin`);
  return new Type(await r.arrayBuffer());
}

async function chargerRun(run) {
  const t0 = performance.now();
  const meta = await (await fetch(`/cache/${run}/meta.json`)).json();
  const noms = ["pts", ...meta.champs, "phase", "jseg", "jmode", "hist"];
  const types = { jseg: Uint32Array, jmode: Uint8Array };
  const d = Object.fromEntries(await Promise.all(noms.map(async (n) => [n, await bin(run, n, types[n] || Float32Array)])));

  for (const t of [GPU.pos, GPU.mode, ...Object.values(GPU.val)]) if (t) gl.deleteTexture(t);
  // Positions centrées sur l'origine de la boîte : float32 garde ainsi sa précision au micron.
  const [x0, y0] = meta.bbox;
  const p = d.pts;
  for (let i = 0; i < p.length; i += 2) { p[i] -= x0; p[i + 1] -= y0; }
  GPU.pos = texture(gl.RG32F, gl.RG, gl.FLOAT, p, p.length / 2);
  GPU.val = {};
  for (const c of meta.champs) GPU.val[c] = texture(gl.R32F, gl.RED, gl.FLOAT, d[c], d[c].length);
  GPU.val.phase = texture(gl.R32F, gl.RED, gl.FLOAT, d.phase, d.phase.length);
  GPU.mode = texture(gl.R8UI, gl.RED_INTEGER, gl.UNSIGNED_BYTE, d.jmode, d.jmode.length);

  if (GPU.vao) gl.deleteVertexArray(GPU.vao);
  GPU.vao = gl.createVertexArray();
  gl.bindVertexArray(GPU.vao);
  const b = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, b);
  gl.bufferData(gl.ARRAY_BUFFER, d.jseg, gl.STATIC_DRAW);
  gl.enableVertexAttribArray(0);
  gl.vertexAttribIPointer(0, 2, gl.UNSIGNED_INT, 0, 0);
  gl.vertexAttribDivisor(0, 1);
  gl.bindVertexArray(null);

  // Décompte des fissures par mode et par frame, une fois pour toutes.
  const nJ = meta.nJoint;
  E.comptes = [];
  for (let f = 0; f < meta.nFrames; f++) {
    const c = { 1: 0, 2: 0, 4: 0 };
    const tranche = d.jmode.subarray(f * nJ, (f + 1) * nJ);
    for (let j = 0; j < nJ; j++) if (tranche[j]) c[tranche[j]]++;
    E.comptes.push(c);
  }
  E.hist = d.hist;
  E.meta = meta;
  E.bw = meta.bbox[2] - x0; E.bh = meta.bbox[3] - y0;
  E.frame = Math.min(E.frame, meta.nFrames - 1);
  curseur.max = meta.nFrames - 1;
  recadrer();
  majChamp(E.champ);
  majFrame(E.frame);
  await image();
  return performance.now() - t0;
}

function dessiner() {
  const dpr = devicePixelRatio || 1;
  const w = cv.clientWidth, h = cv.clientHeight;
  if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) {
    cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
  }
  gl.viewport(0, 0, cv.width, cv.height);
  gl.clearColor(0, 0, 0, 0);
  gl.clear(gl.COLOR_BUFFER_BIT);
  if (!E.meta) return;
  const m = E.meta, ech = [2 * E.k / w, 2 * E.k / h];

  const champ = E.champ, P = progTri;
  gl.useProgram(P.p);
  gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, GPU.pos);
  gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, GPU.val[champ.id]);
  gl.uniform1i(P.u.uPos, 0); gl.uniform1i(P.u.uVal, 1);
  gl.uniform1i(P.u.uFrame, E.frame); gl.uniform1i(P.u.uFrameVal, champ.cat ? 0 : E.frame);
  gl.uniform1i(P.u.uNVert, m.nVert); gl.uniform1i(P.u.uNTri, m.nTri);
  gl.uniform2f(P.u.uCentre, E.cx, E.cy); gl.uniform2f(P.u.uEchelle, ech[0], ech[1]);
  gl.uniform1i(P.u.uCat, champ.cat ? 1 : 0); gl.uniform1i(P.u.uBandes, BANDES);
  if (!champ.cat) { const [a, b] = m.bornes[champ.id]; gl.uniform1f(P.u.uMin, a); gl.uniform1f(P.u.uMax, b); }
  gl.uniform3fv(P.u["uCat3[0]"], ["#2FBECB", "#E5A54B", "#A78BFA", "#6BCB77", "#E77C8D", "#5A9CF8"].flatMap(hex));
  gl.bindVertexArray(null);
  gl.drawArrays(gl.TRIANGLES, 0, m.nVert);

  const J = progJoint;
  gl.useProgram(J.p);
  gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, GPU.pos);
  gl.activeTexture(gl.TEXTURE2); gl.bindTexture(gl.TEXTURE_2D, GPU.mode);
  gl.uniform1i(J.u.uPos, 0); gl.uniform1i(J.u.uMode, 2);
  gl.uniform1i(J.u.uFrame, E.frame); gl.uniform1i(J.u.uNVert, m.nVert);
  gl.uniform1i(J.u.uNJoint, m.nJoint); gl.uniform1i(J.u.uMasque, E.masque);
  gl.uniform2f(J.u.uCentre, E.cx, E.cy); gl.uniform2f(J.u.uEchelle, ech[0], ech[1]);
  gl.uniform2f(J.u.uPixel, 2 / w, 2 / h);
  gl.uniform1f(J.u.uDemiLarg, Math.max(0.6, Math.min(1.6, E.k * 2.5e-4)));
  gl.uniform3fv(J.u["uCoul[0]"], [css("--trac"), css("--cis"), css("--pre")].flatMap(hex));
  gl.bindVertexArray(GPU.vao);
  gl.drawArraysInstanced(gl.TRIANGLES, 0, 6, m.nJoint);
  gl.bindVertexArray(null);
}

// Boucle de rendu à la demande : on ne redessine que si quelque chose a changé.
let attenteImage = [];
function boucle() {
  if (E.lecture && performance.now() - (E.dernierPas || 0) > 120) {
    E.dernierPas = performance.now();
    majFrame((E.frame + 1) % E.meta.nFrames);
  }
  if (E.sale) { dessiner(); dessinerCourbe(); E.sale = false; }
  const a = attenteImage; attenteImage = [];
  a.forEach((r) => r());
  requestAnimationFrame(boucle);
}
const image = () => new Promise((r) => attenteImage.push(r));
const marquer = () => { E.sale = true; };

// ---------------------------------------------------------------- navigation 2D
function recadrer() {
  const w = cv.clientWidth, h = cv.clientHeight;
  E.k = 0.86 * Math.min(w / E.bw, h / E.bh);
  E.cx = E.bw / 2; E.cy = E.bh / 2;
  marquer();
}
function versMonde(px, py) {
  const r = cv.getBoundingClientRect();
  return [E.cx + (px - r.left - r.width / 2) / E.k, E.cy - (py - r.top - r.height / 2) / E.k];
}
cv.addEventListener("wheel", (ev) => {
  ev.preventDefault();
  const [mx, my] = versMonde(ev.clientX, ev.clientY);
  const f = Math.exp(-ev.deltaY * 0.0015);
  E.k *= f;
  E.cx = mx - (mx - E.cx) / f; E.cy = my - (my - E.cy) / f;
  marquer();
}, { passive: false });
let glisse = null;
cv.addEventListener("pointerdown", (ev) => { glisse = [ev.clientX, ev.clientY]; cv.setPointerCapture(ev.pointerId); });
cv.addEventListener("pointermove", (ev) => {
  if (!glisse) return;
  E.cx -= (ev.clientX - glisse[0]) / E.k; E.cy += (ev.clientY - glisse[1]) / E.k;
  glisse = [ev.clientX, ev.clientY]; marquer();
});
cv.addEventListener("pointerup", () => { glisse = null; });
cv.addEventListener("dblclick", recadrer);
new ResizeObserver(marquer).observe(cv);

// ---------------------------------------------------------------- courbe q-epsilon
const cc = document.getElementById("courbe");
const ctx = cc.getContext("2d");
const COURBE = { marge: { g: 52, d: 16, h: 14, b: 38 } };

function graduations(a, b, n) {
  const pas0 = (b - a) / n, p10 = 10 ** Math.floor(Math.log10(pas0));
  const pas = [1, 2, 2.5, 5, 10].map((m) => m * p10).find((s) => s >= pas0);
  const t = [];
  for (let v = Math.ceil(a / pas) * pas; v <= b + 1e-9; v += pas) t.push(+v.toFixed(10));
  return t;
}

function dessinerCourbe() {
  if (!E.hist) return;
  const dpr = devicePixelRatio || 1, w = cc.clientWidth, h = cc.clientHeight;
  if (cc.width !== Math.round(w * dpr) || cc.height !== Math.round(h * dpr)) {
    cc.width = Math.round(w * dpr); cc.height = Math.round(h * dpr);
  }
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  const H = E.hist, n = H.length / 2, M = COURBE.marge;
  let xa = Infinity, xb = -Infinity, ya = 0, yb = -Infinity;
  for (let i = 0; i < n; i++) {
    xa = Math.min(xa, H[2 * i]); xb = Math.max(xb, H[2 * i]);
    ya = Math.min(ya, H[2 * i + 1]); yb = Math.max(yb, H[2 * i + 1]);
  }
  xa = Math.min(xa, 0); yb *= 1.08;
  const X = (v) => M.g + (v - xa) / (xb - xa) * (w - M.g - M.d);
  const Y = (v) => h - M.b - (v - ya) / (yb - ya) * (h - M.h - M.b);
  COURBE.X = X; COURBE.Y = Y;

  ctx.font = `10.5px ${css("--mono")}`;
  ctx.strokeStyle = css("--border"); ctx.fillStyle = css("--texte3"); ctx.lineWidth = 1;
  ctx.textAlign = "center";
  for (const v of graduations(xa, xb, 5)) {
    ctx.beginPath(); ctx.moveTo(X(v) + .5, M.h); ctx.lineTo(X(v) + .5, h - M.b); ctx.stroke();
    ctx.fillText(v.toLocaleString("fr-FR"), X(v), h - M.b + 15);
  }
  ctx.textAlign = "right";
  for (const v of graduations(ya, yb, 5)) {
    ctx.beginPath(); ctx.moveTo(M.g, Y(v) + .5); ctx.lineTo(w - M.d, Y(v) + .5); ctx.stroke();
    ctx.fillText(v.toLocaleString("fr-FR"), M.g - 6, Y(v) + 3.5);
  }
  ctx.fillStyle = css("--texte2"); ctx.font = `11px ${css("--police")}`;
  ctx.textAlign = "center"; ctx.fillText("ε axial (%)", M.g + (w - M.g - M.d) / 2, h - 6);
  ctx.save(); ctx.translate(13, M.h + (h - M.h - M.b) / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillText("q (MPa)", 0, 0); ctx.restore();

  const k = E.meta.frameVersHist[E.frame];
  ctx.lineJoin = "round";
  ctx.strokeStyle = css("--border-fort"); ctx.lineWidth = 1.6;   // partie à venir, discrète
  ctx.beginPath();
  for (let i = k; i < n; i++) (i === k ? ctx.moveTo : ctx.lineTo).call(ctx, X(H[2 * i]), Y(H[2 * i + 1]));
  ctx.stroke();
  ctx.strokeStyle = css("--accent"); ctx.lineWidth = 2;          // partie parcourue
  ctx.beginPath();
  for (let i = 0; i <= k; i++) (i === 0 ? ctx.moveTo : ctx.lineTo).call(ctx, X(H[2 * i]), Y(H[2 * i + 1]));
  ctx.stroke();
  ctx.fillStyle = css("--bg2"); ctx.strokeStyle = css("--texte"); ctx.lineWidth = 2;
  ctx.beginPath(); ctx.arc(X(H[2 * k]), Y(H[2 * k + 1]), 5, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
}

cc.addEventListener("click", (ev) => {
  if (!COURBE.X) return;
  const r = cc.getBoundingClientRect(), px = ev.clientX - r.left, py = ev.clientY - r.top;
  const H = E.hist, L = E.meta.frameVersHist;
  let meilleur = 0, dmin = Infinity;
  L.forEach((k, f) => {
    const d = (COURBE.X(H[2 * k]) - px) ** 2 + (COURBE.Y(H[2 * k + 1]) - py) ** 2;
    if (d < dmin) { dmin = d; meilleur = f; }
  });
  majFrame(meilleur);
});
new ResizeObserver(marquer).observe(cc);

// ---------------------------------------------------------------- panneaux
const curseur = document.getElementById("curseur");
const fmt = (v, n = 1) => v.toLocaleString("fr-FR", { minimumFractionDigits: n, maximumFractionDigits: n });

function majFrame(f) {
  E.frame = f;
  curseur.value = f;
  const m = E.meta, t = m.temps[f] * 1e3, c = E.comptes[f];
  const k = m.frameVersHist[f], q = E.hist[2 * k + 1], eps = E.hist[2 * k];
  document.getElementById("lectureTemps").textContent =
    `frame ${String(f).padStart(2, "0")} / ${m.nFrames - 1}   ·   t = ${fmt(t, 3)} ms`;
  document.getElementById("infoVue").innerHTML =
    `<div class="grand">${fmt(q, 1)} MPa</div><div class="sous">q à ε = ${fmt(eps, 3)} %  ·  σ₃ = ${fmt(m.sigma3_MPa, 0)} MPa</div>`;
  let qp = -Infinity, ep = 0;
  for (let i = 0; i < E.hist.length / 2; i++) if (E.hist[2 * i + 1] > qp) { qp = E.hist[2 * i + 1]; ep = E.hist[2 * i]; }
  const lignes = [
    ["Temps", `${fmt(t, 3)} ms`], ["ε axial", `${fmt(eps, 3)} %`], ["q", `${fmt(q, 1)} MPa`],
    ["q pic du run", `${fmt(qp, 1)} MPa à ${fmt(ep, 3)} %`],
    ["Rompus en traction", c[1].toLocaleString("fr-FR")], ["Rompus en cisaillement", c[2].toLocaleString("fr-FR")],
    ["Pré-rompus", c[4].toLocaleString("fr-FR")], ["Éléments · joints", `${m.nTri.toLocaleString("fr-FR")} · ${m.nJoint.toLocaleString("fr-FR")}`],
  ];
  document.getElementById("tabStats").innerHTML = lignes.map(([a, b]) => `<tr><td>${a}</td><td>${b}</td></tr>`).join("");
  marquer();
}

function majChamp(ch) {
  E.champ = ch;
  document.querySelectorAll("#segChamp button").forEach((b) => b.classList.toggle("actif", b.dataset.id === ch.id));
  const barre = document.getElementById("echelleBarre"), ticks = document.getElementById("echelleTicks");
  if (ch.cat) {
    document.getElementById("echelleTitre").textContent = "Phase";
    barre.style.background = "linear-gradient(#2FBECB 0 33%, #E5A54B 33% 66%, #A78BFA 66%)";
    ticks.innerHTML = ["quartz", "feldspath", "biotite"].map((n, i) => `<span style="top:${(i + .5) * 33.3}%">${n}</span>`).join("");
  } else {
    const [a, b] = E.meta.bornes[ch.id];
    document.getElementById("echelleTitre").textContent = `${ch.nom} (${ch.unite})`;
    const c = [[.16, .27, .82], [.12, .72, .86], [.33, .79, .35], [.96, .82, .24], [.85, .21, .19]];
    const coul = (t) => { const x = t * 4, i = Math.min(Math.floor(x), 3), u = x - i;
      return `rgb(${c[i].map((v, j) => Math.round(255 * (v + (c[i + 1][j] - v) * u))).join(",")})`; };
    const arrets = [];
    for (let i = 0; i < BANDES; i++) {
      const col = coul((BANDES - 1 - i + .5) / BANDES);
      arrets.push(`${col} ${(i / BANDES * 100).toFixed(2)}%`, `${col} ${((i + 1) / BANDES * 100).toFixed(2)}%`);
    }
    barre.style.background = `linear-gradient(${arrets.join(",")})`;
    ticks.innerHTML = Array.from({ length: BANDES + 1 }, (_, i) =>
      `<span style="top:${i / BANDES * 100}%">${fmt((b - (b - a) * i / BANDES) * ch.k, ch.unite === "%" ? 3 : 0)}</span>`).join("");
  }
  marquer();
}

const seg = document.getElementById("segChamp");
for (const ch of CHAMPS) {
  const b = document.createElement("button");
  b.textContent = ch.nom; b.dataset.id = ch.id;
  b.onclick = () => majChamp(ch);
  seg.appendChild(b);
}
document.querySelectorAll(".puce").forEach((b) => b.onclick = () => {
  E.masque ^= +b.dataset.mode; b.classList.toggle("actif"); marquer();
});
curseur.oninput = () => majFrame(+curseur.value);
const btn = document.getElementById("btnLecture");
const basculer = () => { E.lecture = !E.lecture; btn.textContent = E.lecture ? "❚❚" : "▶"; };
btn.onclick = basculer;
addEventListener("keydown", (ev) => {
  if (!E.meta || ev.target.tagName === "SELECT") return;
  if (ev.key === "ArrowRight") majFrame(Math.min(E.frame + 1, E.meta.nFrames - 1));
  else if (ev.key === "ArrowLeft") majFrame(Math.max(E.frame - 1, 0));
  else if (ev.key === " ") { ev.preventDefault(); basculer(); }
});
const selRun = document.getElementById("selRun");
selRun.onchange = () => chargerRun(selRun.value);

// ---------------------------------------------------------------- démarrage et mesures
let tacheMax = 0, tacheMaxDebut = 0;
try {
  new PerformanceObserver((l) => l.getEntries().forEach((e) => {
    if (e.duration > tacheMax) { tacheMax = e.duration; tacheMaxDebut = e.startTime; }
  }))
    .observe({ type: "longtask", buffered: true });
} catch (e) { /* navigateur sans longtask */ }

async function demarrer() {
  requestAnimationFrame(boucle);
  const runs = await (await fetch("/api/runs")).json();
  selRun.innerHTML = runs.map((r) => `<option>${r}</option>`).join("");
  const voulu = new URLSearchParams(location.search).get("run") || "F7_disc_gbm_P020";
  selRun.value = runs.includes(voulu) ? voulu : runs[0];
  const n3 = await chargerRun(selRun.value);
  const n1 = performance.now();                 // depuis le début de la navigation
  if (MESURE) await mesurer(n1, n3);
  else {
    const q = new URLSearchParams(location.search);
    const ch = CHAMPS.find((c) => c.id === q.get("champ"));
    if (ch) majChamp(ch);
    majFrame(Math.min(+(q.get("frame") ?? 20), E.meta.nFrames - 1));
  }
}

async function mesurer(n1, n3) {
  // N4 : changement de frame mesuré jusqu'à l'image suivante (inclut l'attente de vsync).
  const n4 = [];
  for (let f = 0; f < E.meta.nFrames; f++) {
    const t = performance.now(); majFrame(f); await image(); n4.push(performance.now() - t);
  }
  // Coût GPU+CPU pur d'un rendu, sans vsync : rendu forcé puis gl.finish().
  const rendu = [];
  for (let i = 0; i < 30; i++) { const t = performance.now(); dessiner(); gl.finish(); rendu.push(performance.now() - t); }
  // N5 : 3 s de zoom et déplacement animés.
  const k0 = E.k, cx0 = E.cx, cy0 = E.cy, d0 = performance.now();
  let images = 0, ecartMax = 0, prec = d0;
  while (performance.now() - d0 < 3000) {
    const s = (performance.now() - d0) / 1000;
    E.k = k0 * (1 + 1.5 * (1 - Math.cos(s * 2.1)));
    E.cx = cx0 + 0.004 * Math.sin(s * 1.7); E.cy = cy0 + 0.008 * Math.sin(s * 1.3);
    marquer(); await image(); images++;
    const n = performance.now(); ecartMax = Math.max(ecartMax, n - prec); prec = n;
  }
  E.k = k0; E.cx = cx0; E.cy = cy0;
  majFrame(Math.min(20, E.meta.nFrames - 1));
  await image();
  const moy = (a) => a.reduce((s, v) => s + v, 0) / a.length;
  const dbg = gl.getExtension("WEBGL_debug_renderer_info");
  const res = {
    pile: "C web local (WebGL2)", navigateur: navigator.userAgent, run: E.meta.run,
    gpu: dbg ? gl.getParameter(dbg.UNMASKED_RENDERER_WEBGL) : "inconnu",
    N1_demarrage_ms: +n1.toFixed(0), N3_ouverture_run_ms: +n3.toFixed(0),
    N4_frame_moy_ms: +moy(n4).toFixed(1), N4_frame_max_ms: +Math.max(...n4).toFixed(1),
    rendu_seul_moy_ms: +moy(rendu).toFixed(2), rendu_seul_max_ms: +Math.max(...rendu).toFixed(2),
    N5_images_par_s: +(images / 3).toFixed(1), N5_ecart_max_entre_images_ms: +ecartMax.toFixed(1),
    N6_tache_longue_max_ms: +tacheMax.toFixed(0), N6_debut_ms: +tacheMaxDebut.toFixed(0),
    shaders_compiles_a_ms: +T_SHADERS.toFixed(0),
    taille_vue_px: [cv.width, cv.height],
  };
  document.getElementById("carteMesures").hidden = false;
  document.getElementById("mesures").textContent = JSON.stringify(res, null, 1);
  await fetch("/api/mesures", { method: "POST", body: JSON.stringify(res, null, 1) });
  window.__mesures = res;
}

demarrer();
