// Maquettes J4 (jetables) : navigation, schéma, aperçu du maillage, file, comparaison.
// Les données sont réelles (caches de la campagne, deck et synthèse générés par le noyau) ;
// seuls les états de la file sont simulés.
const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const fr = (v, n = 1) => v.toLocaleString("fr-FR", { minimumFractionDigits: n, maximumFractionDigits: n });
const bin = async (run, nom, T = Float32Array) => new T(await (await fetch(`/cache/${run}/${nom}.bin`)).arrayBuffer());
const meta = async (run) => (await fetch(`/cache/${run}/meta.json`)).json();

// ------------------------------------------------------------------ onglets
const dessins = {};
function ouvrir(e) {
  document.querySelectorAll(".ecran").forEach((x) => x.classList.toggle("actif", x.id === "e-" + e));
  document.querySelectorAll("#onglets button").forEach((b) => b.classList.toggle("actif", b.dataset.e === e));
  requestAnimationFrame(() => dessins[e] && dessins[e]());
}
document.querySelectorAll("#onglets button").forEach((b) => (b.onclick = () => ouvrir(b.dataset.e)));
document.querySelectorAll("[data-aller]").forEach((b) => (b.onclick = () => ouvrir(b.dataset.aller)));

function toile(cv) {
  const dpr = devicePixelRatio || 1, w = cv.clientWidth, h = cv.clientHeight;
  cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
  const c = cv.getContext("2d");
  c.setTransform(dpr, 0, 0, dpr, 0, 0);
  return [c, w, h];
}

// ------------------------------------------------------------------ Essai : schéma
function schema() {
  const s = document.getElementById("schema"), k = 6.2;          // px par mm
  const W = 36 * k, H = 72 * k, x0 = 260 - W / 2, y0 = 310 - H / 2;
  const fl = (x1, y1, x2, y2) => `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${css("--accent")}" stroke-width="2" marker-end="url(#pointe)"/>`;
  let g = `<defs><marker id="pointe" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="${css("--accent")}"/></marker>
    <pattern id="grains" width="36" height="36" patternUnits="userSpaceOnUse"><path d="M0 12 L14 6 L26 14 L36 8 M14 6 L12 0 M26 14 L24 26 L36 30 M24 26 L10 30 L0 24 M10 30 L12 36" stroke="${css("--border-fort")}" fill="none"/></pattern></defs>`;
  g += `<rect x="${x0}" y="${y0}" width="${W}" height="${H}" fill="${css("--bg3")}" stroke="${css("--texte2")}"/>`;
  g += `<rect x="${x0}" y="${y0}" width="${W}" height="${H}" fill="url(#grains)"/>`;
  // quelques pré-rompus diffus, dessinés comme des tirets
  let a = 7;
  for (let i = 0; i < 60; i++) {
    a = (a * 9301 + 49297) % 233280;
    const u = a / 233280; a = (a * 9301 + 49297) % 233280; const v = a / 233280;
    a = (a * 9301 + 49297) % 233280; const t = a / 233280 * Math.PI;
    const cx = x0 + 8 + u * (W - 16), cy = y0 + 8 + v * (H - 16);
    g += `<line x1="${cx}" y1="${cy}" x2="${cx + 7 * Math.cos(t)}" y2="${cy + 7 * Math.sin(t)}" stroke="${css("--pre")}" stroke-width="2"/>`;
  }
  // plateaux et chargement axial
  for (const [y, sens] of [[y0 - 16, 1], [y0 + H, -1]]) {
    g += `<rect x="${x0 - 20}" y="${y}" width="${W + 40}" height="16" rx="3" fill="${css("--border-fort")}"/>`;
    const ya = sens > 0 ? y - 44 : y + 60;
    g += fl(260, ya, 260, sens > 0 ? y - 4 : y + 20);
  }
  g += `<text x="${260 + 14}" y="${y0 - 40}" fill="${css("--texte2")}" font-size="13">plateaux, 0,1 m/s</text>`;
  // confinement latéral
  for (let i = 0; i < 7; i++) {
    const y = y0 + (i + 0.5) * H / 7;
    g += fl(x0 - 58, y, x0 - 6, y) + fl(x0 + W + 58, y, x0 + W + 6, y);
  }
  g += `<text x="${x0 - 60}" y="${y0 - 6}" fill="${css("--texte2")}" font-size="13" text-anchor="end">σ₃ = 20 MPa</text>`;
  // cotes
  g += `<text x="260" y="${y0 + H + 52}" fill="${css("--texte3")}" font-size="12" text-anchor="middle" font-family="monospace">36 mm</text>`;
  g += `<text x="${x0 + W + 74}" y="310" fill="${css("--texte3")}" font-size="12" font-family="monospace" transform="rotate(90 ${x0 + W + 74} 310)" text-anchor="middle">72 mm</text>`;
  g += `<text x="260" y="608" fill="${css("--texte3")}" font-size="11.5" text-anchor="middle">GBM 3 phases · grains 3 mm · éléments 0,7 mm · 5 % de joints pré-rompus</text>`;
  s.innerHTML = g;
}

async function deck() {
  const t = (await (await fetch("maquette_deck.txt")).text()).split("\r").join("");
  document.getElementById("deck").innerHTML = t.split("\n").map((l) =>
    l.startsWith("#") ? `<span class="c">${l}</span>` :
    l.includes("=") ? `<span class="k">${l.split("=")[0]}</span>=${l.split("=").slice(1).join("=")}` : l).join("\n");
}

// ------------------------------------------------------------------ Maillage : triangles réels
let MAILLAGE = null;
async function chargerMaillage() {
  const run = "F7_apercu_maillage", m = await meta(run);
  MAILLAGE = { m, pts: await bin(run, "pts"), ph: await bin(run, "phase"),
               seg: await bin(run, "jseg", Uint32Array), mode: await bin(run, "jmode", Uint8Array) };
}
function dessinerMaillage() {
  if (!MAILLAGE) return;
  const { m, pts, ph, seg, mode } = MAILLAGE;
  const [c, w, h] = toile(document.getElementById("cvMaillage"));
  const [x0, y0, x1, y1] = m.bbox, k = 0.9 * Math.min(w / (x1 - x0), h / (y1 - y0));
  const ox = (w - k * (x1 - x0)) / 2, oy = (h + k * (y1 - y0)) / 2;
  const X = (x) => ox + (x - x0) * k, Y = (y) => oy - (y - y0) * k;
  const coul = [css("--quartz"), css("--feldspath"), css("--biotite")];
  for (let p = 0; p < 3; p++) {
    c.fillStyle = coul[p]; c.strokeStyle = coul[p]; c.lineWidth = 0.6;
    c.beginPath();
    for (let t = 0; t < m.nTri; t++) {
      if (ph[t] !== p) continue;
      const i = 6 * t;
      c.moveTo(X(pts[i]), Y(pts[i + 1])); c.lineTo(X(pts[i + 2]), Y(pts[i + 3])); c.lineTo(X(pts[i + 4]), Y(pts[i + 5]));
      c.closePath();
    }
    c.fill(); c.stroke();
  }
  c.strokeStyle = css("--pre"); c.lineWidth = 1.6; c.beginPath();
  for (let j = 0; j < m.nJoint; j++) {
    if (!mode[j]) continue;
    const a = 2 * seg[2 * j], b = 2 * seg[2 * j + 1];
    c.moveTo(X(pts[a]), Y(pts[a + 1])); c.lineTo(X(pts[b]), Y(pts[b + 1]));
  }
  c.stroke();
}
dessins.maillage = dessinerMaillage;

// ------------------------------------------------------------------ File
const FILE = [
  ["lot", "LF_fragile · 20 MPa · 3 runs"],
  ["F1_homog_P020", "Voronoï, une phase", "fini", 1, "1 h 28", ""],
  ["F2_gbm_P020", "GBM 3 phases, grains uniformes", "fini", 1, "1 h 31", ""],
  ["F3_tailles_P020", "GBM 3 phases, tailles par phase", "cours", 0.61, "54 min", "≈ 35 min"],
  ["lot", "Variation σ₃ · F7 GBM + 5 % pré-rompus · 4 runs"],
  ["F7_disc_gbm_P000", "σ₃ = 0 MPa (UCS)", "cours", 0.23, "21 min", "≈ 1 h 10"],
  ["F7_disc_gbm_P010", "σ₃ = 10 MPa", "attente", 0, "", ""],
  ["F7_disc_gbm_P020", "σ₃ = 20 MPa", "fini", 1, "1 h 33", ""],
  ["F7_disc_gbm_P040", "σ₃ = 40 MPa", "attente", 0, "", ""],
  ["lot", "Hors lot"],
  ["A1_elast_P040", "contraste élastique seul, 40 MPa", "echec", 0.02, "3 min", "code 3 · NaN"],
  ["F9_apercu", "aperçu de maillage", "attente", 0, "", "≈ 2 s"],
];
const LIB = { fini: "terminé", cours: "en cours", attente: "en attente", echec: "échec" };
function file() {
  let h = `<tr><th>Run</th><th>État</th><th>Avancement</th><th>Écoulé</th><th>Reste</th><th></th></tr>`;
  for (const r of FILE) {
    if (r[0] === "lot") { h += `<tr class="groupe-lot"><td colspan="6">${r[1]}</td></tr>`; continue; }
    const [nom, desc, et, av, ec, reste] = r;
    const act = et === "cours" ? "■" : et === "attente" ? "▲▼" : "↻";
    h += `<tr class="${nom === "F3_tailles_P020" ? "sel" : ""}"><td><div class="nom">${nom}</div><div class="desc">${desc}</div></td>
      <td><span class="badge ${et}">${LIB[et]}</span></td>
      <td><div class="barre-av"><i style="width:${av * 100}%"></i></div></td>
      <td class="mono2">${ec}</td><td class="mono2">${reste}</td>
      <td><div class="actions-ligne">${[...act.matchAll(/./gu)].map((x) => `<button class="ico" title="${{ "■": "arrêter", "▲": "monter", "▼": "descendre", "↻": "relancer (l'ancien dossier est renommé)" }[x[0]]}">${x[0]}</button>`).join("")}<button class="ico" title="ouvrir les résultats">→</button></div></td></tr>`;
  }
  document.getElementById("tabFile").innerHTML = h;
}

function axes(c, w, h, xa, xb, ya, yb, M, lx, ly) {
  const X = (v) => M.g + (v - xa) / (xb - xa) * (w - M.g - M.d), Y = (v) => h - M.b - (v - ya) / (yb - ya) * (h - M.h - M.b);
  c.strokeStyle = css("--border"); c.fillStyle = css("--texte3"); c.font = `10.5px ${css("--mono")}`; c.lineWidth = 1;
  const pas = (a, b) => { const p = (b - a) / 5, e = 10 ** Math.floor(Math.log10(p)); return [1, 2, 2.5, 5, 10].map((m) => m * e).find((s) => s >= p); };
  c.textAlign = "center";
  for (let v = Math.ceil(xa / pas(xa, xb)) * pas(xa, xb); v <= xb; v += pas(xa, xb)) {
    c.beginPath(); c.moveTo(X(v), M.h); c.lineTo(X(v), h - M.b); c.stroke(); c.fillText(fr(v, v % 1 ? 1 : 0), X(v), h - M.b + 14);
  }
  c.textAlign = "right";
  for (let v = Math.ceil(ya / pas(ya, yb)) * pas(ya, yb); v <= yb; v += pas(ya, yb)) {
    c.beginPath(); c.moveTo(M.g, Y(v)); c.lineTo(w - M.d, Y(v)); c.stroke(); c.fillText(fr(v, 0), M.g - 6, Y(v) + 3);
  }
  c.fillStyle = css("--texte2"); c.font = `11px ${css("--police")}`; c.textAlign = "center";
  c.fillText(lx, M.g + (w - M.g - M.d) / 2, h - 4);
  c.save(); c.translate(12, M.h + (h - M.h - M.b) / 2); c.rotate(-Math.PI / 2); c.fillText(ly, 0, 0); c.restore();
  return [X, Y];
}
let HIST = {};
async function hist(run) { return HIST[run] ||= await bin(run, "hist"); }
async function miniCourbe() {
  const H = await hist("F3_tailles_P020"), n = Math.floor(H.length / 2 * 0.61);
  const [c, w, h] = toile(document.getElementById("miniCourbe"));
  const [X, Y] = axes(c, w, h, 0, 1, 0, 100, { g: 40, d: 10, h: 8, b: 30 }, "ε axial (%)", "q (MPa)");
  c.strokeStyle = css("--accent"); c.lineWidth = 2; c.beginPath();
  for (let i = 0; i < n; i++) (i ? c.lineTo : c.moveTo).call(c, X(H[2 * i]), Y(H[2 * i + 1]));
  c.stroke();
}
dessins.file = miniCourbe;

// ------------------------------------------------------------------ Résultats : comparaison
const SERIES = [["F1_homog_P020", "#2FBECB"], ["F3_tailles_P020", "#E5A54B"], ["F7_disc_gbm_P020", "#A78BFA"]];
let SYN = [];
async function resultats() {
  SYN = await (await fetch("maquette_synthese.json")).json();
  const coul = Object.fromEntries(SERIES);
  document.getElementById("listeRuns").innerHTML = SYN.map((r) =>
    `<label><input type="checkbox" ${coul[r.run] ? "checked" : ""}><span class="trait" style="background:${coul[r.run] || "transparent"}"></span><span>${r.run}</span><span class="mono2">${fr(r.q_pic_MPa)}</span></label>`).join("");
  const cols = [["run", "Run"], ["sigma3_atteint_MPa", "σ₃ atteint"], ["q_pic_MPa", "q pic MPa"], ["E_secante_GPa", "E séc. GPa"],
    ["eps_pic_pct", "ε pic %"], ["chute_post_pic", "chute post-pic"], ["part_intergranulaire", "intergranulaire"],
    ["n_tension", "traction"], ["n_cisaillement", "cisaillement"], ["prerompus_libres", "pré-rompus libres"]];
  const f = (k, v) => v == null ? "—" : k === "run" ? v : k.startsWith("n_") || k === "prerompus_libres" ? v.toLocaleString("fr-FR")
    : k === "eps_pic_pct" ? fr(v, 3) : k === "chute_post_pic" || k === "part_intergranulaire" ? fr(v, 2) : fr(v, 1);
  document.getElementById("tabSynthese").innerHTML = `<tr>${cols.map((c) => `<th>${c[1]}</th>`).join("")}</tr>` +
    SYN.map((r) => `<tr>${cols.map(([k]) => `<td>${f(k, r[k])}</td>`).join("")}</tr>`).join("");
}
async function comparaison() {
  const [c, w, h] = toile(document.getElementById("cvComp"));
  const [X, Y] = axes(c, w, h, 0, 1.4, 0, 100, { g: 44, d: 12, h: 8, b: 32 }, "ε axial (%)", "q (MPa)");
  for (const [run, col] of SERIES) {
    const H = await hist(run);
    c.strokeStyle = col; c.lineWidth = 1.8; c.beginPath();
    for (let i = 0; i < H.length / 2; i++) (i ? c.lineTo : c.moveTo).call(c, X(H[2 * i]), Y(H[2 * i + 1]));
    c.stroke();
  }
  c.font = `11px ${css("--police")}`; c.textAlign = "left";
  SERIES.forEach(([run, col], i) => { c.fillStyle = col; c.fillRect(w - 190, 18 + 16 * i, 12, 3); c.fillStyle = css("--texte2"); c.fillText(run, w - 172, 22 + 16 * i); });
  // champ σyy de la dernière frame de F7, fissures par mode
  const run = "F7_disc_gbm_P020", m = await meta(run);
  const f = m.nFrames - 1, P = await bin(run, "pts"), S = await bin(run, "sigmaYY"), seg = await bin(run, "jseg", Uint32Array), md = await bin(run, "jmode", Uint8Array);
  const [d, w2, h2] = toile(document.getElementById("cvChamp"));
  const [x0, y0, x1, y1] = m.bbox, k = 0.92 * Math.min(w2 / (x1 - x0), h2 / (y1 - y0));
  const ox = (w2 - k * (x1 - x0)) / 2, oy = (h2 + k * (y1 - y0)) / 2, off = f * m.nVert * 2;
  const XX = (x) => ox + (x - x0) * k, YY = (y) => oy - (y - y0) * k;
  const [a, b] = m.bornes.sigmaYY, st = [[.16, .27, .82], [.12, .72, .86], [.33, .79, .35], [.96, .82, .24], [.85, .21, .19]];
  for (let t = 0; t < m.nTri; t++) {
    let u = Math.min(.9999, Math.max(0, (S[f * m.nTri + t] - a) / (b - a)));
    u = (Math.floor(u * 12) + .5) / 12;
    const x = u * 4, i = Math.min(3, Math.floor(x)), r = st[i].map((v, j) => Math.round(255 * (v + (st[i + 1][j] - v) * (x - i))));
    d.fillStyle = d.strokeStyle = `rgb(${r})`; d.lineWidth = .5;
    const q = off + 6 * t;
    d.beginPath(); d.moveTo(XX(P[q]), YY(P[q + 1])); d.lineTo(XX(P[q + 2]), YY(P[q + 3])); d.lineTo(XX(P[q + 4]), YY(P[q + 5])); d.closePath(); d.fill(); d.stroke();
  }
  for (const [mo, col] of [[1, "--trac"], [2, "--cis"]]) {
    d.strokeStyle = css(col); d.lineWidth = 1.4; d.beginPath();
    for (let j = 0; j < m.nJoint; j++) {
      if (md[f * m.nJoint + j] !== mo) continue;
      const p1 = off + 2 * seg[2 * j], p2 = off + 2 * seg[2 * j + 1];
      d.moveTo(XX(P[p1]), YY(P[p1 + 1])); d.lineTo(XX(P[p2]), YY(P[p2 + 1]));
    }
    d.stroke();
  }
}
dessins.resultats = comparaison;

// ------------------------------------------------------------------ démarrage
schema(); deck(); file(); resultats();
chargerMaillage().then(() => dessins.maillage && document.getElementById("e-maillage").classList.contains("actif") && dessinerMaillage());
const voulu = new URLSearchParams(location.search).get("onglet");
if (voulu) ouvrir(voulu);
addEventListener("resize", () => { const a = document.querySelector("#onglets .actif").dataset.e; dessins[a] && dessins[a](); });
