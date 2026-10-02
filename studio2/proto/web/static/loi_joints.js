// Maquette J4 : éditeur de la loi des joints (spec 007 §2.6).
// Les formules sont celles du solveur g1 :
//   setJointLengths (FdemSolver.cpp)  dnE = ft / pj ; dnF = dnE + Gf / (ft I) ; slipF = GfII / (c I)
//   YanSoftening.hpp                  f(D) de Yan (Munjiza), I = intégrale de f sur [0, 1]
//   montée                            pj = k E / h (insertionPenaltyFactor), ou parabole ft (2r - r²)
// Le mode II est tracé sous une hypothèse de forme (montée à pente pj puis pic * f(D)) qui reste
// à vérifier contre tools/yan_point.cpp avant l'implémentation.
const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const fr = (v, n = 2) => v.toLocaleString("fr-FR", { minimumFractionDigits: n, maximumFractionDigits: n });
const H_ELEM = 0.0007;                                     // m, taille d'élément de la campagne

const PRESETS = {
  fragile: { nom: "fragile", E: 15e9, ft: 1.3e6, Gf: 3.8, k: 4, rise: "lin", soft: "yan", a: 0.63, b: 1.8, c: 6, coh: 16.4e6, phi: 23, gfs: 22.105 },
  bohus:   { nom: "bohus",   E: 52e9, ft: 34e6,  Gf: 70,  k: 4, rise: "lin", soft: "yan", a: 0.63, b: 1.8, c: 6, coh: 13.6e6, phi: 13.4, gfs: 10 },
};
let REF = { ...PRESETS.fragile }, S = { ...REF }, SN = 20e6;

// ------------------------------------------------------------------ la loi
function fD(D, a, b, c) {
  if (D <= 0) return 1;
  if (D >= 1) return 0;
  const s = a + b;
  if (Math.abs(s) < 1e-12 || Math.abs(1 - s) < 1e-12) return a * (1 - D) + b * (1 - D) ** c;
  const g = 1 - (s - 1) / s * Math.exp((a + c * b) / (s * (1 - s)) * D);
  return Math.min(1, Math.max(0, g * (a * (1 - D) + b * (1 - D) ** c)));
}
function integrale(a, b, c, n = 1024) {                      // Simpson, comme integralFD
  let s = fD(0, a, b, c) + fD(1, a, b, c);
  for (let i = 1; i < n; i++) s += fD(i / n, a, b, c) * (i % 2 ? 4 : 2);
  return s / n / 3;
}
function derive(s) {
  const pj = s.k * s.E / H_ELEM;
  const I = s.soft === "yan" ? integrale(s.a, s.b, s.c) : 0.5;
  const dnE = s.ft / pj, dnF = dnE + s.Gf / (s.ft * I);
  const GfII = s.gfs * s.Gf, slipF = GfII / (s.coh * I);
  return { pj, I, dnE, dnF, GfII, slipF };
}
const adou = (s, D) => (s.soft === "yan" ? fD(D, s.a, s.b, s.c) : 1 - Math.min(1, Math.max(0, D)));
function sigma(s, d, x) {
  if (x <= d.dnE) {
    const r = x / d.dnE;
    return s.rise === "para" ? s.ft * (2 * r - r * r) : d.pj * x;
  }
  return s.ft * adou(s, (x - d.dnE) / (d.dnF - d.dnE));
}
function tau(s, d, x, sn) {
  const tp = s.coh + Math.tan(s.phi * Math.PI / 180) * sn, sE = tp / d.pj;
  if (x <= sE) return d.pj * x;
  return tp * adou(s, (x - sE) / d.slipF);
}

// ------------------------------------------------------------------ tracés
function toile(cv) {
  const dpr = devicePixelRatio || 1, w = cv.clientWidth, h = cv.clientHeight;
  if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) { cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr); }
  const c = cv.getContext("2d");
  c.setTransform(dpr, 0, 0, dpr, 0, 0);
  c.clearRect(0, 0, w, h);
  return [c, w, h];
}
function pasGrad(a, b) { const p = (b - a) / 5, e = 10 ** Math.floor(Math.log10(p)); return [1, 2, 2.5, 5, 10].map((m) => m * e).find((v) => v >= p); }
function repere(c, w, h, R, M, lx, ly, petit) {
  const X = (v) => M.g + (v - R.x0) / (R.x1 - R.x0) * (w - M.g - M.d);
  const Y = (v) => h - M.b - (v - R.y0) / (R.y1 - R.y0) * (h - M.h - M.b);
  const Xi = (p) => R.x0 + (p - M.g) / (w - M.g - M.d) * (R.x1 - R.x0);
  const Yi = (p) => R.y0 + (h - M.b - p) / (h - M.h - M.b) * (R.y1 - R.y0);
  c.strokeStyle = css("--border"); c.fillStyle = css("--texte3"); c.lineWidth = 1;
  c.font = `${petit ? 9.5 : 10.5}px ${css("--mono")}`;
  const px = pasGrad(R.x0, R.x1), py = pasGrad(R.y0, R.y1), dec = (p) => Math.max(0, -Math.floor(Math.log10(p)));
  c.textAlign = "center";
  for (let v = Math.ceil(R.x0 / px) * px; v <= R.x1 + 1e-12; v += px) {
    c.beginPath(); c.moveTo(X(v), M.h); c.lineTo(X(v), h - M.b); c.stroke(); c.fillText(fr(v, dec(px)), X(v), h - M.b + 13);
  }
  c.textAlign = "right";
  for (let v = Math.ceil(R.y0 / py) * py; v <= R.y1 + 1e-12; v += py) {
    c.beginPath(); c.moveTo(M.g, Y(v)); c.lineTo(w - M.d, Y(v)); c.stroke(); c.fillText(fr(v, dec(py)), M.g - 5, Y(v) + 3);
  }
  if (!petit) {
    c.fillStyle = css("--texte2"); c.font = `11px ${css("--police")}`; c.textAlign = "center";
    c.fillText(lx, M.g + (w - M.g - M.d) / 2, h - 3);
    c.save(); c.translate(11, M.h + (h - M.h - M.b) / 2); c.rotate(-Math.PI / 2); c.fillText(ly, 0, 0); c.restore();
  }
  return { X, Y, Xi, Yi };
}
function courbe(c, T, f, x0, x1, col, larg, tirets) {
  c.strokeStyle = col; c.lineWidth = larg; c.setLineDash(tirets ? [5, 4] : []);
  c.beginPath();
  const n = 600;
  for (let i = 0; i <= n; i++) { const x = x0 + (x1 - x0) * i / n; (i ? c.lineTo : c.moveTo).call(c, T.X(x), T.Y(f(x))); }
  c.stroke(); c.setLineDash([]);
}
function poignee(c, x, y, forme, actif) {
  c.fillStyle = actif ? css("--accent") : css("--bg2"); c.strokeStyle = css("--accent"); c.lineWidth = 2;
  c.beginPath();
  if (forme === "losange") { c.moveTo(x, y - 7); c.lineTo(x + 7, y); c.lineTo(x, y + 7); c.lineTo(x - 7, y); c.closePath(); }
  else c.arc(x, y, 6, 0, 2 * Math.PI);
  c.fill(); c.stroke();
}

// Échelles figées pendant le glissement (sinon le repère fuit sous la souris).
let RI = null, RL = null, RII = null;
function echelles() {
  const d = derive(S), dr = derive(REF);
  RI = { x0: 0, x1: 1.6 * Math.max(d.dnF, dr.dnF) * 1e6, y0: 0, y1: 1.6 * Math.max(S.ft, REF.ft) / 1e6 };
  RL = { x0: 0, x1: 3 * Math.max(d.dnE, dr.dnE) * 1e9, y0: 0, y1: RI.y1 };
  const tp = (s) => s.coh + Math.tan(s.phi * Math.PI / 180) * 60e6;
  RII = { x0: 0, x1: 1.6 * Math.max(d.slipF, dr.slipF) * 1e6, y0: 0, y1: 1.15 * Math.max(tp(S), tp(REF)) / 1e6 };
}
const POIG = {};
function dessiner() {
  const d = derive(S), dr = derive(REF);
  // ---- mode I
  let [c, w, h] = toile(document.getElementById("cvI"));
  let T = repere(c, w, h, RI, { g: 50, d: 280, h: 10, b: 34 }, "ouverture δ (µm)", "σ (MPa)");
  courbe(c, T, (x) => sigma(REF, dr, x * 1e-6) / 1e6, 0, Math.min(RI.x1, dr.dnF * 1e6 * 1.02), css("--texte3"), 1.5, true);
  courbe(c, T, (x) => sigma(S, d, x * 1e-6) / 1e6, 0, Math.min(RI.x1, d.dnF * 1e6 * 1.02), css("--accent"), 2.4);
  // aire Gf hachurée en léger
  c.fillStyle = css("--accent-doux"); c.globalAlpha = 0.55; c.beginPath(); c.moveTo(T.X(d.dnE * 1e6), T.Y(0));
  for (let i = 0; i <= 200; i++) { const x = d.dnE + (d.dnF - d.dnE) * i / 200; c.lineTo(T.X(x * 1e6), T.Y(sigma(S, d, x) / 1e6)); }
  c.lineTo(T.X(d.dnF * 1e6), T.Y(0)); c.fill(); c.globalAlpha = 1;
  c.fillStyle = css("--texte2"); c.font = `12px ${css("--police")}`; c.textAlign = "left";
  c.fillText(`Gf = ${fr(S.Gf, 1)} J/m²`, T.X(d.dnE * 1e6 + 0.12 * (d.dnF - d.dnE) * 1e6), T.Y(0.18 * S.ft / 1e6));
  const xs = d.dnE + 0.5 * (d.dnF - d.dnE);
  POIG.pic = [T.X(d.dnE * 1e6), T.Y(S.ft / 1e6), T];
  POIG.fin = [T.X(d.dnF * 1e6), T.Y(0), T];
  POIG.forme = [T.X(xs * 1e6), T.Y(sigma(S, d, xs) / 1e6), T];
  poignee(c, ...POIG.pic.slice(0, 2), "rond", glisse === "pic");
  poignee(c, ...POIG.fin.slice(0, 2), "rond", glisse === "fin");
  if (S.soft === "yan") poignee(c, ...POIG.forme.slice(0, 2), "losange", glisse === "forme");
  // ---- loupe sur la montée
  [c, w, h] = toile(document.getElementById("cvLoupe"));
  T = repere(c, w, h, RL, { g: 34, d: 8, h: 20, b: 18 }, "", "", true);
  courbe(c, T, (x) => sigma(REF, dr, x * 1e-9) / 1e6, 0, RL.x1, css("--texte3"), 1.3, true);
  courbe(c, T, (x) => sigma(S, d, x * 1e-9) / 1e6, 0, RL.x1, css("--accent"), 2);
  c.fillStyle = css("--texte3"); c.font = `9.5px ${css("--mono")}`; c.textAlign = "right"; c.fillText("δ (nm)", w - 8, h - 22);
  POIG.montee = [T.X(d.dnE * 1e9), T.Y(S.ft / 1e6), T];
  poignee(c, ...POIG.montee.slice(0, 2), "rond", glisse === "montee");
  // ---- mode II
  [c, w, h] = toile(document.getElementById("cvII"));
  T = repere(c, w, h, RII, { g: 50, d: 14, h: 30, b: 34 }, "glissement s (µm)", "τ (MPa)");
  courbe(c, T, (x) => tau(REF, dr, x * 1e-6, SN) / 1e6, 0, RII.x1, css("--texte3"), 1.5, true);
  courbe(c, T, (x) => tau(S, d, x * 1e-6, SN) / 1e6, 0, RII.x1, css("--cis"), 2.4);
  const tp = S.coh + Math.tan(S.phi * Math.PI / 180) * SN;
  POIG.coh = [T.X(tp / d.pj * 1e6), T.Y(tp / 1e6), T];
  POIG.finII = [T.X((tp / d.pj + d.slipF) * 1e6), T.Y(0), T];
  poignee(c, ...POIG.coh.slice(0, 2), "rond", glisse === "coh");
  poignee(c, ...POIG.finII.slice(0, 2), "rond", glisse === "finII");
  c.fillStyle = css("--texte3"); c.font = `10.5px ${css("--police")}`; c.textAlign = "left";
  c.fillText("forme du mode II à vérifier contre tools/yan_point.cpp", T.X(0) + 6, 22);
  panneaux(d);
}

// ------------------------------------------------------------------ panneaux
const champs = { inFt: ["ft", 1e6, 2], inGf: ["Gf", 1, 1], inK: ["k", 1, 2], inC: ["c", 1, 2], inCoh: ["coh", 1e6, 1], inPhi: ["phi", 1, 1], inGfs: ["gfs", 1, 2] };
function panneaux(d) {
  for (const [id, [k, u, n]] of Object.entries(champs)) {
    const el = document.getElementById(id);
    if (document.activeElement !== el) el.value = fr(S[k] / u, n);
  }
  document.getElementById("inC").disabled = S.soft !== "yan";
  document.querySelectorAll("#segMontee button").forEach((b) => b.classList.toggle("actif", b.dataset.v === S.rise));
  document.querySelectorAll("#segAdou button").forEach((b) => b.classList.toggle("actif", b.dataset.v === S.soft));
  const lcz = S.E * S.Gf / S.ft ** 2, dr = derive(REF), lr = REF.E * REF.Gf / REF.ft ** 2;
  const lignes = [
    ["", "Réglage", "Préréglage"],
    ["Ouverture au pic δe", `${fr(d.dnE * 1e9, 1)} nm`, `${fr(dr.dnE * 1e9, 1)} nm`],
    ["Ouverture de rupture δF", `${fr(d.dnF * 1e6, 2)} µm`, `${fr(dr.dnF * 1e6, 2)} µm`],
    ["Intégrale I de f(D)", fr(d.I, 4), fr(dr.I, 4)],
    ["Glissement de rupture", `${fr(d.slipF * 1e6, 1)} µm`, `${fr(dr.slipF * 1e6, 1)} µm`],
    ["Zone cohésive E Gf / ft²", lcz > 1 ? `${fr(lcz, 1)} m` : `${fr(lcz * 1e3, 2)} mm`, lr > 1 ? `${fr(lr, 1)} m` : `${fr(lr * 1e3, 2)} mm`],
    ["Rapport à l'élément (> 2)", fr(lcz / H_ELEM, 1), fr(lr / H_ELEM, 1)],
    ["Pic τ à σn", `${fr((S.coh + Math.tan(S.phi * Math.PI / 180) * SN) / 1e6, 1)} MPa`, `${fr((REF.coh + Math.tan(REF.phi * Math.PI / 180) * SN) / 1e6, 1)} MPa`],
  ];
  document.getElementById("tabLect").innerHTML = lignes.map((l, i) =>
    i ? `<tr><td>${l[0]}</td><td>${l[1]}</td><td style="color:var(--texte3)">${l[2]}</td></tr>` : `<tr><th></th><th>${l[1]}</th><th>${l[2]}</th></tr>`).join("") +
    (lcz / H_ELEM < 2 ? `<tr><td colspan="3" style="color:var(--attention);font-family:var(--police)">Zone cohésive sous 2 éléments : la fissure n'est plus résolue par le maillage.</td></tr>` : "");
  const cles = [];
  const ec = (cle, v, r, f = (x) => x) => { if (Math.abs(v - r) > 1e-9 * Math.max(1, Math.abs(r))) cles.push(`<b>${cle} = ${f(v)}</b>   # préréglage ${f(r)}`); };
  ec("ft", S.ft, REF.ft, (x) => x.toPrecision(4));
  ec("Gf", S.Gf, REF.Gf, (x) => x.toPrecision(4));
  ec("insertionPenaltyFactor", S.k, REF.k, (x) => x.toPrecision(3));
  ec("cohesion", S.coh, REF.coh, (x) => x.toPrecision(4));
  ec("frictionDeg", S.phi, REF.phi, (x) => x.toPrecision(3));
  ec("gfShearFactor", S.gfs, REF.gfs, (x) => x.toPrecision(4));
  if (S.soft === "yan") ec("yanC", S.c, REF.c, (x) => x.toPrecision(3));
  if (S.rise !== REF.rise) cles.push(`<b>jointElastic = ${S.rise === "para" ? "parabolic" : "linear"}</b>`);
  if (S.soft !== REF.soft) cles.push(`<b>jointSoftening = ${S.soft === "yan" ? "yan" : "linear"}</b>`);
  document.getElementById("cles").innerHTML = cles.length ? cles.join("<br>") : "aucune : loi du préréglage";
}

// ------------------------------------------------------------------ cliquer-glisser
let glisse = null;
function attacher(cvId, noms) {
  const cv = document.getElementById(cvId);
  const local = (ev) => { const r = cv.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; };
  const proche = (p) => noms.find((n) => POIG[n] && (n !== "forme" || S.soft === "yan") && Math.hypot(POIG[n][0] - p[0], POIG[n][1] - p[1]) < 12);
  cv.addEventListener("pointerdown", (ev) => { glisse = proche(local(ev)) || null; if (glisse) { cv.setPointerCapture(ev.pointerId); dessiner(); } });
  cv.addEventListener("pointermove", (ev) => {
    const p = local(ev);
    if (!glisse) { cv.style.cursor = proche(p) ? "grab" : "default"; return; }
    cv.style.cursor = "grabbing";
    const T = POIG[glisse][2], d = derive(S);
    if (glisse === "pic") S.ft = Math.max(0.01 * REF.ft, T.Yi(p[1]) * 1e6);
    else if (glisse === "fin") S.Gf = Math.max(0.02 * REF.Gf, (T.Xi(p[0]) * 1e-6 - d.dnE) * S.ft * d.I);
    else if (glisse === "montee") S.k = Math.max(0.05, S.ft / Math.max(1e-12, T.Xi(p[0]) * 1e-9) * H_ELEM / S.E);
    else if (glisse === "forme") {
      // f(0,5) visée -> yanC : recherche sur une grille, f n'étant pas monotone en c partout
      const cible = Math.min(0.97, Math.max(0.01, T.Yi(p[1]) * 1e6 / S.ft));
      let best = S.c, e = Infinity;
      for (let c = 0.3; c <= 40; c += 0.05) { const v = Math.abs(fD(0.5, S.a, S.b, c) - cible); if (v < e) { e = v; best = c; } }
      const Gf = S.Gf; S.c = +best.toFixed(2); S.Gf = Gf;      // Gf garde sa valeur : δF s'ajuste
    } else if (glisse === "coh") S.coh = Math.max(0.01 * REF.coh, T.Yi(p[1]) * 1e6 - Math.tan(S.phi * Math.PI / 180) * SN);
    else if (glisse === "finII") {
      const tp = S.coh + Math.tan(S.phi * Math.PI / 180) * SN;
      S.gfs = Math.max(0.05, (T.Xi(p[0]) * 1e-6 - tp / d.pj) * S.coh * d.I / S.Gf);
    }
    dessiner();
  });
  cv.addEventListener("pointerup", () => { glisse = null; dessiner(); });
}
attacher("cvI", ["pic", "fin", "forme"]);
attacher("cvLoupe", ["montee"]);
attacher("cvII", ["coh", "finII"]);

for (const [id, [k, u]] of Object.entries(champs)) {
  document.getElementById(id).addEventListener("change", (ev) => {
    const v = parseFloat(ev.target.value.replace(/\s/g, "").replace(",", "."));
    if (isFinite(v) && v > 0) { S[k] = v * u; echelles(); }
    dessiner();
  });
}
document.querySelectorAll("#segMontee button").forEach((b) => (b.onclick = () => { S.rise = b.dataset.v; dessiner(); }));
document.querySelectorAll("#segAdou button").forEach((b) => (b.onclick = () => { S.soft = b.dataset.v; echelles(); dessiner(); }));
document.getElementById("inSn").oninput = (ev) => { SN = +ev.target.value * 1e6; document.getElementById("lblSn").textContent = `${ev.target.value} MPa`; dessiner(); };
document.getElementById("selPreset").onchange = (ev) => { REF = { ...PRESETS[ev.target.value] }; S = { ...REF }; echelles(); dessiner(); };
document.getElementById("btnReset").onclick = () => { S = { ...REF }; echelles(); dessiner(); };
addEventListener("resize", dessiner);

// ?demo=1 : un réglage modifié, pour la capture de la maquette
if (new URLSearchParams(location.search).has("demo")) { echelles(); S.ft *= 1.35; S.c = 3; S.Gf *= 0.8; S.rise = "para"; }
else echelles();
dessiner();
