// Vue 2D des champs et des fissures, en WebGL2 (moteur validé au jalon J1).
//
// Toutes les frames d'un run sont chargées UNE fois dans des textures GPU (positions, champs,
// modes de joints). Changer de frame = changer un entier uniforme : rien n'est recopié, et le
// coût ne dépend pas de la taille du maillage. Sommet k -> triangle k/3 (nœuds propres par
// élément, comme les écrit le solveur).
import { binaire, css } from "./commun.js";

const LARG = 4096;
export const BANDES = 12;

// Champs connus : facteur d'affichage et unité. Les champs strain* exigent writeStrainFields.
export const CHAMPS = {
  sigmaYY: { nom: "σyy", unite: "MPa", k: 1e-6 }, sigmaXX: { nom: "σxx", unite: "MPa", k: 1e-6 },
  sigmaXY: { nom: "σxy", unite: "MPa", k: 1e-6 }, vonMises: { nom: "von Mises", unite: "MPa", k: 1e-6 },
  strainYY: { nom: "εyy", unite: "%", k: 100 }, strainXX: { nom: "εxx", unite: "%", k: 100 },
  strainXY: { nom: "εxy", unite: "%", k: 100 }, epsXX: { nom: "εxx co-roté", unite: "%", k: 100 },
  phase: { nom: "Phases", unite: "", k: 1, cat: true },
};
export const ARC = [[.16, .27, .82], [.12, .72, .86], [.33, .79, .35], [.96, .82, .24], [.85, .21, .19]];
export function couleurArc(t) {
  const x = Math.min(.9999, Math.max(0, t)) * 4, i = Math.min(Math.floor(x), 3), u = x - i;
  return ARC[i].map((v, j) => v + (ARC[i + 1][j] - v) * u);
}
const hex = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16) / 255);

const ENTETE = `#version 300 es
precision highp float; precision highp int; precision highp sampler2D; precision highp usampler2D;
uniform sampler2D uPos; uniform int uFrame, uNVert; uniform vec2 uCentre, uEchelle;
ivec2 tc(int i) { return ivec2(i % ${LARG}, i / ${LARG}); }
vec2 pos(int v) { return texelFetch(uPos, tc(uFrame * uNVert + v), 0).rg; }
vec2 ndc(vec2 p) { return (p - uCentre) * uEchelle; }
`;
const VS_TRI = ENTETE + `
uniform sampler2D uVal; uniform int uFrameVal, uNTri;
flat out float vVal;
void main() {
  int v = gl_VertexID;
  vVal = texelFetch(uVal, tc(uFrameVal * uNTri + v / 3), 0).r;
  gl_Position = vec4(ndc(pos(v)), 0.0, 1.0);
}`;
const FS_TRI = `#version 300 es
precision highp float;
flat in float vVal; uniform float uMin, uMax; uniform int uCat, uBandes; uniform vec3 uCat3[6];
out vec4 o;
vec3 arc(float t) {
  vec3 c[5] = vec3[5](vec3(.16,.27,.82), vec3(.12,.72,.86), vec3(.33,.79,.35), vec3(.96,.82,.24), vec3(.85,.21,.19));
  float x = clamp(t, 0., 1.) * 4.; int i = min(int(x), 3);
  return mix(c[i], c[i + 1], x - float(i));
}
void main() {
  if (uCat == 1) { o = vec4(uCat3[int(vVal) % 6], 1.); return; }
  float t = (vVal - uMin) / max(uMax - uMin, 1e-30);
  t = (floor(clamp(t, 0., .9999) * float(uBandes)) + .5) / float(uBandes);
  o = vec4(arc(t), 1.);
}`;
const VS_JOINT = ENTETE + `
layout(location = 0) in uvec2 aSeg;
uniform usampler2D uMode; uniform int uNJoint, uMasque; uniform vec2 uPixel; uniform float uDemiLarg;
uniform vec3 uCoul[3];
out vec3 vCoul;
void main() {
  int m = int(texelFetch(uMode, tc(uFrame * uNJoint + gl_InstanceID), 0).r);
  if (m == 0 || (m & uMasque) == 0) { gl_Position = vec4(2., 2., 2., 1.); return; }
  vCoul = uCoul[m == 1 ? 0 : (m == 2 ? 1 : 2)];
  int c = gl_VertexID;
  float bout = (c == 1 || c == 2 || c == 4) ? 1. : 0.;
  float cote = (c == 2 || c == 4 || c == 5) ? 1. : -1.;
  vec2 a = ndc(pos(int(aSeg.x))), b = ndc(pos(int(aSeg.y)));
  vec2 d = (b - a) / uPixel; float l = length(d);
  vec2 n = l > 0. ? vec2(-d.y, d.x) / l : vec2(0.);
  gl_Position = vec4(mix(a, b, bout) + n * cote * uDemiLarg * uPixel, 0., 1.);
}`;
const FS_JOINT = `#version 300 es
precision highp float; in vec3 vCoul; out vec4 o; void main() { o = vec4(vCoul, 1.); }`;

export class VueChamps {
  constructor(canvas) {
    this.cv = canvas;
    this.gl = canvas.getContext("webgl2", { antialias: true, alpha: true });
    if (!this.gl) throw new Error("WebGL2 indisponible");
    this.tri = this._prog(VS_TRI, FS_TRI);
    this.jt = this._prog(VS_JOINT, FS_JOINT);
    this.meta = null; this.tex = {}; this.frame = 0; this.champ = "sigmaYY"; this.masque = 7;
    this.k = 1; this.cx = 0; this.cy = 0; this.sale = true; this.bornes = null;
    this.couleursPhases = [css("--quartz"), css("--feldspath"), css("--biotite"), "#6BCB77", "#E77C8D", "#5A9CF8"];
    this._souris();
    new ResizeObserver(() => this.marquer()).observe(canvas);
    this.apresDessin = null;                 // rappel après chaque image (calques SVG)
    const boucle = () => {
      if (this.sale) { this.dessiner(); this.sale = false; this.apresDessin?.(); }
      requestAnimationFrame(boucle);
    };
    requestAnimationFrame(boucle);
  }

  _prog(vs, fs) {
    const gl = this.gl, p = gl.createProgram();
    for (const [t, src] of [[gl.VERTEX_SHADER, vs], [gl.FRAGMENT_SHADER, fs]]) {
      const s = gl.createShader(t);
      gl.shaderSource(s, src); gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
      gl.attachShader(p, s);
    }
    gl.linkProgram(p);
    const u = {};
    for (let i = 0; i < gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS); i++) {
      const n = gl.getActiveUniform(p, i).name;
      u[n] = gl.getUniformLocation(p, n);
    }
    return { p, u };
  }

  _texture(interne, format, type, donnees, n) {
    const gl = this.gl, h = Math.ceil(n / LARG), comp = donnees.length / n;
    const plein = new donnees.constructor(LARG * h * comp);
    plein.set(donnees);
    const t = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, t);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    gl.texImage2D(gl.TEXTURE_2D, 0, interne, LARG, h, 0, format, type, plein);
    for (const q of [gl.TEXTURE_MIN_FILTER, gl.TEXTURE_MAG_FILTER]) gl.texParameteri(gl.TEXTURE_2D, q, gl.NEAREST);
    return t;
  }

  // base = "/cache/<id>" ; meta = meta.json déjà lu
  async charger(base, meta) {
    const gl = this.gl;
    const noms = ["pts", ...meta.champs, "phase", "jseg", "jmode"];
    const types = { jseg: Uint32Array, jmode: Uint8Array };
    const d = Object.fromEntries(await Promise.all(noms.map(async (n) => [n, await binaire(`${base}/${n}.bin`, types[n] || Float32Array)])));
    for (const t of Object.values(this.tex)) gl.deleteTexture(t);
    const [x0, y0] = meta.bbox, p = d.pts;
    for (let i = 0; i < p.length; i += 2) { p[i] -= x0; p[i + 1] -= y0; }
    this.tex = { pos: this._texture(gl.RG32F, gl.RG, gl.FLOAT, p, p.length / 2),
                 mode: this._texture(gl.R8UI, gl.RED_INTEGER, gl.UNSIGNED_BYTE, d.jmode, d.jmode.length),
                 phase: this._texture(gl.R32F, gl.RED, gl.FLOAT, d.phase, d.phase.length) };
    for (const c of meta.champs) this.tex[c] = this._texture(gl.R32F, gl.RED, gl.FLOAT, d[c], d[c].length);
    if (this.vao) gl.deleteVertexArray(this.vao);
    this.vao = gl.createVertexArray();
    gl.bindVertexArray(this.vao);
    gl.bindBuffer(gl.ARRAY_BUFFER, gl.createBuffer());
    gl.bufferData(gl.ARRAY_BUFFER, d.jseg, gl.STATIC_DRAW);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribIPointer(0, 2, gl.UNSIGNED_INT, 0, 0);
    gl.vertexAttribDivisor(0, 1);
    gl.bindVertexArray(null);
    // décompte des fissures par mode et par frame, une fois pour toutes
    this.comptes = [];
    for (let f = 0; f < meta.nFrames; f++) {
      const c = { 1: 0, 2: 0, 4: 0 }, tr = d.jmode.subarray(f * meta.nJoint, (f + 1) * meta.nJoint);
      for (let j = 0; j < tr.length; j++) if (tr[j]) c[tr[j]]++;
      this.comptes.push(c);
    }
    this.meta = meta;
    this.bw = meta.bbox[2] - x0; this.bh = meta.bbox[3] - y0;
    this.frame = Math.min(this.frame, meta.nFrames - 1);
    if (!meta.champs.includes(this.champ) && this.champ !== "phase") this.champ = meta.champs.includes("sigmaYY") ? "sigmaYY" : meta.champs[0];
    this.recadrer();
  }

  champsDisponibles() { return this.meta ? [...this.meta.champs, "phase"].filter((c) => CHAMPS[c]) : []; }
  bornesChamp() { return this.bornes || this.meta.bornes[this.champ]; }
  marquer() { this.sale = true; }
  setFrame(f) { this.frame = f; this.marquer(); }
  setChamp(c) { this.champ = c; this.bornes = null; this.marquer(); }

  recadrer() {
    if (!this.meta) return;
    const w = this.cv.clientWidth, h = this.cv.clientHeight;
    this.k = 0.88 * Math.min(w / this.bw, h / this.bh);
    this.cx = this.bw / 2; this.cy = this.bh / 2;
    this.marquer();
  }

  // Coordonnées du SOLVEUR (m) <-> pixels CSS du canvas. La vue travaille en coordonnées
  // décalées de l'origine de la boîte (précision float32) : on rajoute bbox[0..1].
  versMonde(clientX, clientY) {
    const r = this.cv.getBoundingClientRect(), [x0, y0] = this.meta.bbox;
    return [x0 + this.cx + (clientX - r.left - r.width / 2) / this.k, y0 + this.cy - (clientY - r.top - r.height / 2) / this.k];
  }
  versEcran(x, y) {
    const [x0, y0] = this.meta.bbox, w = this.cv.clientWidth, h = this.cv.clientHeight;
    return [w / 2 + (x - x0 - this.cx) * this.k, h / 2 - (y - y0 - this.cy) * this.k];
  }

  _souris() {
    const cv = this.cv;
    // Un outil (dessin de fissure) prend la main sur le déplacement : { down, move, up } en coordonnées solveur.
    this.outil = null;
    const monde = (x, y) => { const r = cv.getBoundingClientRect(); return [this.cx + (x - r.left - r.width / 2) / this.k, this.cy - (y - r.top - r.height / 2) / this.k]; };
    cv.addEventListener("wheel", (ev) => {
      ev.preventDefault();
      const [mx, my] = monde(ev.clientX, ev.clientY), f = Math.exp(-ev.deltaY * 0.0015);
      this.k *= f; this.cx = mx - (mx - this.cx) / f; this.cy = my - (my - this.cy) / f;
      this.marquer();
    }, { passive: false });
    let g = null;
    cv.addEventListener("pointerdown", (ev) => {
      cv.setPointerCapture(ev.pointerId);
      if (this.outil && this.meta) return this.outil.down(...this.versMonde(ev.clientX, ev.clientY));
      g = [ev.clientX, ev.clientY];
    });
    cv.addEventListener("pointermove", (ev) => {
      if (this.outil && this.meta) return this.outil.move(...this.versMonde(ev.clientX, ev.clientY));
      if (!g) return;
      this.cx -= (ev.clientX - g[0]) / this.k; this.cy += (ev.clientY - g[1]) / this.k;
      g = [ev.clientX, ev.clientY]; this.marquer();
    });
    cv.addEventListener("pointerup", (ev) => {
      if (this.outil && this.meta) this.outil.up(...this.versMonde(ev.clientX, ev.clientY));
      g = null;
    });
    cv.addEventListener("dblclick", () => this.recadrer());
  }

  dessiner() {
    const gl = this.gl, cv = this.cv, dpr = devicePixelRatio || 1, w = cv.clientWidth, h = cv.clientHeight;
    if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) { cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr); }
    gl.viewport(0, 0, cv.width, cv.height);
    gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT);
    if (!this.meta || !w) return;
    const m = this.meta, ech = [2 * this.k / w, 2 * this.k / h], cat = this.champ === "phase";
    const P = this.tri;
    gl.useProgram(P.p);
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, this.tex.pos);
    gl.activeTexture(gl.TEXTURE1); gl.bindTexture(gl.TEXTURE_2D, this.tex[this.champ]);
    gl.uniform1i(P.u.uPos, 0); gl.uniform1i(P.u.uVal, 1);
    gl.uniform1i(P.u.uFrame, this.frame); gl.uniform1i(P.u.uFrameVal, cat ? 0 : this.frame);
    gl.uniform1i(P.u.uNVert, m.nVert); gl.uniform1i(P.u.uNTri, m.nTri);
    gl.uniform2f(P.u.uCentre, this.cx, this.cy); gl.uniform2f(P.u.uEchelle, ech[0], ech[1]);
    gl.uniform1i(P.u.uCat, cat ? 1 : 0); gl.uniform1i(P.u.uBandes, BANDES);
    if (!cat) { const [a, b] = this.bornesChamp(); gl.uniform1f(P.u.uMin, a); gl.uniform1f(P.u.uMax, b); }
    gl.uniform3fv(P.u["uCat3[0]"], this.couleursPhases.flatMap(hex));
    gl.bindVertexArray(null);
    gl.drawArrays(gl.TRIANGLES, 0, m.nVert);
    const J = this.jt;
    gl.useProgram(J.p);
    gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, this.tex.pos);
    gl.activeTexture(gl.TEXTURE2); gl.bindTexture(gl.TEXTURE_2D, this.tex.mode);
    gl.uniform1i(J.u.uPos, 0); gl.uniform1i(J.u.uMode, 2);
    gl.uniform1i(J.u.uFrame, this.frame); gl.uniform1i(J.u.uNVert, m.nVert);
    gl.uniform1i(J.u.uNJoint, m.nJoint); gl.uniform1i(J.u.uMasque, this.masque);
    gl.uniform2f(J.u.uCentre, this.cx, this.cy); gl.uniform2f(J.u.uEchelle, ech[0], ech[1]);
    gl.uniform2f(J.u.uPixel, 2 / w, 2 / h);
    gl.uniform1f(J.u.uDemiLarg, Math.max(0.6, Math.min(1.6, this.k * 2.5e-4)));
    gl.uniform3fv(J.u["uCoul[0]"], [css("--trac"), css("--cis"), css("--pre")].flatMap(hex));
    gl.bindVertexArray(this.vao);
    gl.drawArraysInstanced(gl.TRIANGLES, 0, 6, m.nJoint);
    gl.bindVertexArray(null);
  }
}
