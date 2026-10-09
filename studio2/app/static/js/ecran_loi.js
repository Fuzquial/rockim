// Écran Loi des joints : modeler la loi cohésive à la souris (loi.js = la loi du solveur), voir ce
// que ça change, et le mesurer par un essai éclair (GBM de 10 grains, 1 à 17 s sous g1).
import { api, css, enAttente, fr, SERIES, toile } from "./commun.js";
import { Courbe } from "./courbe.js";
import { abonner, initialiser, lireChoix, lireFormulaire, modifier } from "./etat_essai.js";
import { fD, Loi } from "./loi.js";

let R, nav, SN = 20e6, glisse = null, echelles = null, POIG = {}, eclairs = [], chargementVu = "tx20", courbeEcl;

// ------------------------------------------------------------------ paramètres effectifs
function params(avecModifs = true) {
  const c = lireChoix(), niv = lireFormulaire().niveaux[c.materiau].materiau;
  const s = avecModifs ? c.surcharges || {} : {}, l = avecModifs ? c.loi || {} : {};
  const E = s.E ?? niv.E, h = c.taille_element_mm / 1000, k = l.penalite ?? 4;
  const ft = s.ft ?? niv.ft, coh = s.cohesion ?? niv.cohesion, Gf = s.Gf ?? niv.Gf, gfs = s.gfShearFactor ?? niv.gfShearFactor;
  return { E, h, k, ft, coh, phi: s.frictionDeg ?? niv.frictionDeg, GfI: Gf, gfs, GfII: gfs * Gf, pj: k * E / h,
    c: l.yan_c ?? 6, adoucissement: l.adoucissement ?? "yan", montee: l.montee ?? "linear" };
}
const loiDe = (p) => new Loi(p);

// Écrit une modification : surcharges de matériau ou forme de la loi, jamais les deux mêlées.
function changer({ surcharges = {}, loi = {} }) {
  const c = lireChoix(), ref = params(false);
  // 4 chiffres significatifs : la souris ne mérite pas plus, et le deck reste lisible.
  const arrondi = (o) => Object.fromEntries(Object.entries(o).map(([k, v]) => [k, typeof v === "number" ? +v.toPrecision(4) : v]));
  const s = { ...(c.surcharges || {}), ...arrondi(surcharges) }, l = { ...(c.loi || {}), ...arrondi(loi) };
  const defautL = { penalite: 4, yan_c: 6, adoucissement: "yan", montee: "linear" };
  const defautS = { ft: ref.ft, Gf: ref.GfI, cohesion: ref.coh, frictionDeg: ref.phi, gfShearFactor: ref.gfs, E: ref.E };
  for (const k of Object.keys(s)) if (Math.abs(s[k] - defautS[k]) <= 1e-9 * Math.abs(defautS[k])) delete s[k];
  for (const k of Object.keys(l)) if (l[k] === defautL[k]) delete l[k];
  modifier({ surcharges: s, loi: l });
}

export async function monter(noeud, ctx) {
  R = noeud; nav = ctx;
  await initialiser();
  R.innerHTML = `
  <section class="carte trace loi-I">
    <div class="carte-tete"><span class="titre">Mode I : contrainte normale en fonction de l'ouverture</span><span class="overline">pointillés = préréglage</span></div>
    <canvas class="principal" id="l-I"></canvas>
    <div class="loupe"><span class="t">Montée · loupe</span><canvas id="l-loupe"></canvas></div>
    <div class="poignee-aide">● pic : hauteur = ft · ● fin : aire = Gf · ◆ forme · loupe : raideur de montée</div>
  </section>
  <div class="colonne-droite">
    <section class="carte"><div class="carte-tete"><span class="titre">Réglages</span><button class="bouton petit" id="l-raz">Revenir au préréglage</button></div>
      <div class="corps-carte">
        <div class="ligne"><label>Montée</label><div class="segments plein" id="l-montee"><button data-v="linear">Linéaire</button><button data-v="parabolic">Parabolique</button></div></div>
        <div class="ligne"><label>Adoucissement</label><div class="segments plein" id="l-adou"><button data-v="linear">Linéaire</button><button data-v="yan">Yan</button></div></div>
        <div class="ligne"><label>ft</label><div class="nombre"><input data-p="ft" data-f="1e6" data-d="2"><span>MPa</span></div></div>
        <div class="ligne"><label>Gf</label><div class="nombre"><input data-p="Gf" data-f="1" data-d="2"><span>J/m²</span></div></div>
        <div class="ligne"><label>Pénalité p</label><div class="nombre"><input data-p="penalite" data-f="1" data-d="2"><span>× E / h</span></div></div>
        <div class="ligne"><label>Forme (yanC)</label><div class="nombre"><input data-p="yan_c" data-f="1" data-d="2"><span>—</span></div></div>
        <div class="ligne"><label>Cohésion · φ</label><div class="deux"><div class="nombre"><input data-p="cohesion" data-f="1e6" data-d="1"><span>MPa</span></div><div class="nombre"><input data-p="frictionDeg" data-f="1" data-d="1"><span>°</span></div></div></div>
        <div class="ligne"><label>GfII / GfI</label><div class="nombre"><input data-p="gfShearFactor" data-f="1" data-d="2"><span>—</span></div></div>
      </div></section>
    <section class="carte"><div class="carte-tete"><span class="titre">Ce que ça change</span></div><table class="tab-phases" id="l-lect" style="margin:8px 0 4px"></table></section>
    <section class="carte"><div class="carte-tete"><span class="titre">Clés du deck</span><span class="overline">écart au préréglage</span></div><div class="cles" id="l-cles"></div></section>
  </div>
  <div class="bas">
    <section class="carte trace">
      <div class="carte-tete"><span class="titre">Mode II : cisaillement en fonction du glissement</span></div>
      <div class="curseur-sn">σn <input type="range" id="l-sn" min="0" max="60" value="20"><span class="mono2" id="l-snv">20 MPa</span></div>
      <canvas class="principal" id="l-II"></canvas>
    </section>
    <section class="carte eclair">
      <div class="carte-tete"><span class="titre">Essai éclair : GBM de 10 grains</span><span class="overline">lance g1 · 1 à 17 s · ±5 % de la pleine taille</span></div>
      <div class="eclair-boutons">
        <button class="bouton" data-ch="traction">Traction</button><button class="bouton" data-ch="ucs">Compression simple</button><button class="bouton primaire" data-ch="tx20">Triaxial 20 MPa</button>
        <button class="bouton petit" id="l-effacer" title="Retire les courbes de la superposition (les calculs restent en cache)">Effacer</button></div>
      <div class="eclair-corps"><canvas id="l-ecl"></canvas><table class="tab-phases" id="l-ecl-tab"></table></div>
      <div class="message" id="l-msg" style="padding:0 12px 8px"></div>
    </section>
  </div>`;
  courbeEcl = new Courbe(R.querySelector("#l-ecl"));
  brancher();
  tracerEclairs();
  abonner(() => { if (!glisse) echelles = null; dessiner(); });
  new ResizeObserver(() => dessiner()).observe(R);
}

export function activer() { echelles = null; dessiner(); }

function brancher() {
  R.querySelector("#l-raz").onclick = () => modifier({ surcharges: {}, loi: {} });
  R.querySelectorAll("#l-montee button").forEach((b) => (b.onclick = () => changer({ loi: { montee: b.dataset.v } })));
  R.querySelectorAll("#l-adou button").forEach((b) => (b.onclick = () => changer({ loi: { adoucissement: b.dataset.v } })));
  R.querySelectorAll("input[data-p]").forEach((i) => (i.onchange = () => {
    const v = parseFloat(i.value.replace(/\s/g, "").replace(",", ".")) * +i.dataset.f;
    if (!(v > 0)) return dessiner();
    echelles = null;
    if (i.dataset.p === "penalite" || i.dataset.p === "yan_c") changer({ loi: { [i.dataset.p]: v } });
    else changer({ surcharges: { [i.dataset.p]: v } });
  }));
  R.querySelector("#l-sn").oninput = (ev) => { SN = +ev.target.value * 1e6; R.querySelector("#l-snv").textContent = `${ev.target.value} MPa`; echelles = null; dessiner(); };
  R.querySelectorAll("[data-ch]").forEach((b) => (b.onclick = () => eclair(b.dataset.ch)));
  R.querySelector("#l-effacer").onclick = () => { eclairs = []; tracerEclairs(); };
  attacher("#l-I", ["pic", "fin", "forme"]);
  attacher("#l-loupe", ["montee"]);
  attacher("#l-II", ["coh", "finII"]);
}

// ------------------------------------------------------------------ tracés
function pasGrad(a, b) { const p = (b - a) / 5, e = 10 ** Math.floor(Math.log10(p)); return [1, 2, 2.5, 5, 10].map((m) => m * e).find((v) => v >= p); }
function repere(c, w, h, Rg, M, lx, ly, petit) {
  const X = (v) => M.g + (v - Rg.x0) / (Rg.x1 - Rg.x0) * (w - M.g - M.d), Y = (v) => h - M.b - (v - Rg.y0) / (Rg.y1 - Rg.y0) * (h - M.h - M.b);
  const Xi = (p) => Rg.x0 + (p - M.g) / (w - M.g - M.d) * (Rg.x1 - Rg.x0), Yi = (p) => Rg.y0 + (h - M.b - p) / (h - M.h - M.b) * (Rg.y1 - Rg.y0);
  c.strokeStyle = css("--border"); c.fillStyle = css("--texte3"); c.lineWidth = 1; c.font = `${petit ? 9.5 : 10.5}px ${css("--mono")}`;
  const px = pasGrad(Rg.x0, Rg.x1), py = pasGrad(Rg.y0, Rg.y1), dec = (p) => Math.max(0, -Math.floor(Math.log10(p)));
  c.textAlign = "center";
  for (let v = Math.ceil(Rg.x0 / px) * px; v <= Rg.x1 + 1e-12; v += px) { c.beginPath(); c.moveTo(X(v), M.h); c.lineTo(X(v), h - M.b); c.stroke(); c.fillText(fr(v, dec(px)), X(v), h - M.b + 13); }
  c.textAlign = "right";
  for (let v = Math.ceil(Rg.y0 / py) * py; v <= Rg.y1 + 1e-12; v += py) { c.beginPath(); c.moveTo(M.g, Y(v)); c.lineTo(w - M.d, Y(v)); c.stroke(); c.fillText(fr(v, dec(py)), M.g - 5, Y(v) + 3); }
  if (!petit) {
    c.fillStyle = css("--texte2"); c.font = `11px ${css("--police")}`; c.textAlign = "center";
    c.fillText(lx, M.g + (w - M.g - M.d) / 2, h - 3);
    c.save(); c.translate(11, M.h + (h - M.h - M.b) / 2); c.rotate(-Math.PI / 2); c.fillText(ly, 0, 0); c.restore();
  }
  return { X, Y, Xi, Yi };
}
function trait(c, T, xs, ys, col, ep, tirets) {
  c.strokeStyle = col; c.lineWidth = ep; c.setLineDash(tirets ? [5, 4] : []);
  c.beginPath(); xs.forEach((x, i) => (i ? c.lineTo : c.moveTo).call(c, T.X(x), T.Y(ys[i]))); c.stroke(); c.setLineDash([]);
}
function poignee(c, x, y, forme, actif) {
  c.fillStyle = actif ? css("--accent") : css("--bg2"); c.strokeStyle = css("--accent"); c.lineWidth = 2; c.beginPath();
  if (forme === "losange") { c.moveTo(x, y - 7); c.lineTo(x + 7, y); c.lineTo(x, y + 7); c.lineTo(x - 7, y); c.closePath(); } else c.arc(x, y, 6, 0, 2 * Math.PI);
  c.fill(); c.stroke();
}
const lin = (a, b, n) => Array.from({ length: n + 1 }, (_, i) => a + (b - a) * i / n);

function dessiner() {
  if (!lireChoix() || !R.classList.contains("actif")) return;
  const p = params(), pr = params(false), L = loiDe(p), Lr = loiDe(pr);
  if (!echelles) {
    const tp = (q) => q.coh + Math.tan(q.phi * Math.PI / 180) * Math.max(SN, 20e6);
    echelles = {
      I: { x0: 0, x1: 1.6 * Math.max(L.dnE + L.ot, Lr.dnE + Lr.ot) * 1e6, y0: 0, y1: 1.6 * Math.max(p.ft, pr.ft) / 1e6 },
      L: { x0: 0, x1: 3 * Math.max(L.dnE, Lr.dnE) * 1e9, y0: 0, y1: 1.6 * Math.max(p.ft, pr.ft) / 1e6 },
      II: { x0: 0, x1: 1.6 * Math.max(L.st, Lr.st) * 1e6, y0: 0, y1: 1.15 * Math.max(tp(p), tp(pr)) / 1e6 },
    };
  }
  const E = echelles;
  // mode I
  let [c, w, h] = toile(R.querySelector("#l-I"));
  let T = repere(c, w, h, E.I, { g: 50, d: 280, h: 10, b: 34 }, "ouverture δ (µm)", "σ (MPa)");
  const xs = lin(0, E.I.x1, 600);
  trait(c, T, xs, xs.map((x) => Lr.sigma(x * 1e-6) / 1e6), css("--texte3"), 1.5, true);
  trait(c, T, xs, xs.map((x) => L.sigma(x * 1e-6) / 1e6), css("--accent"), 2.4);
  c.fillStyle = css("--accent-doux"); c.globalAlpha = .55; c.beginPath(); c.moveTo(T.X(L.dnE * 1e6), T.Y(0));
  lin(L.dnE, L.dnE + L.ot, 200).forEach((x) => c.lineTo(T.X(x * 1e6), T.Y(L.sigma(x) / 1e6)));
  c.lineTo(T.X((L.dnE + L.ot) * 1e6), T.Y(0)); c.fill(); c.globalAlpha = 1;
  c.fillStyle = css("--texte2"); c.font = `12px ${css("--police")}`; c.textAlign = "left";
  c.fillText(`Gf = ${fr(p.GfI, 2)} J/m²`, T.X((L.dnE + .12 * L.ot) * 1e6), T.Y(.18 * p.ft / 1e6));
  const xm = L.dnE + .5 * L.ot;
  POIG.pic = [T.X(L.dnE * 1e6), T.Y(p.ft / 1e6), T]; POIG.fin = [T.X((L.dnE + L.ot) * 1e6), T.Y(0), T]; POIG.forme = [T.X(xm * 1e6), T.Y(L.sigma(xm) / 1e6), T];
  poignee(c, POIG.pic[0], POIG.pic[1], "rond", glisse === "pic");
  poignee(c, POIG.fin[0], POIG.fin[1], "rond", glisse === "fin");
  if (L.yan) poignee(c, POIG.forme[0], POIG.forme[1], "losange", glisse === "forme");
  // loupe
  [c, w, h] = toile(R.querySelector("#l-loupe"));
  T = repere(c, w, h, E.L, { g: 34, d: 8, h: 20, b: 18 }, "", "", true);
  const xl = lin(0, E.L.x1, 300);
  trait(c, T, xl, xl.map((x) => Lr.sigma(x * 1e-9) / 1e6), css("--texte3"), 1.3, true);
  trait(c, T, xl, xl.map((x) => L.sigma(x * 1e-9) / 1e6), css("--accent"), 2);
  c.fillStyle = css("--texte3"); c.font = `9.5px ${css("--mono")}`; c.textAlign = "right"; c.fillText("δ (nm)", w - 8, h - 22);
  POIG.montee = [T.X(L.dnE * 1e9), T.Y(p.ft / 1e6), T];
  poignee(c, POIG.montee[0], POIG.montee[1], "rond", glisse === "montee");
  // mode II
  [c, w, h] = toile(R.querySelector("#l-II"));
  T = repere(c, w, h, E.II, { g: 50, d: 14, h: 30, b: 34 }, "glissement s (µm)", "τ (MPa)");
  const ss = lin(0, E.II.x1, 900);
  trait(c, T, ss, Lr.mode2(ss.map((x) => x * 1e-6), -SN).map((t) => t / 1e6), css("--texte3"), 1.5, true);
  const t2 = L.mode2(ss.map((x) => x * 1e-6), -SN).map((t) => t / 1e6);
  trait(c, T, ss, t2, css("--cis"), 2.4);
  const tp = p.coh + L.tanPhi * SN, sE = tp / p.pj;
  POIG.coh = [T.X(sE * 1e6), T.Y(tp / 1e6), T]; POIG.finII = [T.X((sE + L.st) * 1e6), T.Y(L.tanPhi * SN / 1e6), T];
  poignee(c, POIG.coh[0], POIG.coh[1], "rond", glisse === "coh");
  poignee(c, POIG.finII[0], POIG.finII[1], "rond", glisse === "finII");
  if (SN > 0) {
    c.strokeStyle = css("--texte3"); c.setLineDash([2, 3]); c.beginPath(); c.moveTo(T.X(0), T.Y(L.tanPhi * SN / 1e6)); c.lineTo(w - 14, T.Y(L.tanPhi * SN / 1e6)); c.stroke(); c.setLineDash([]);
    c.fillStyle = css("--texte3"); c.font = `10.5px ${css("--police")}`; c.textAlign = "right";
    c.fillText(`frottement résiduel tan φ · σn = ${fr(L.tanPhi * SN / 1e6, 1)} MPa`, w - 18, T.Y(L.tanPhi * SN / 1e6) - 5);
  }
  panneaux(p, pr, L, Lr);
}

function panneaux(p, pr, L, Lr) {
  const vals = { ft: p.ft, Gf: p.GfI, penalite: p.k, yan_c: p.c, cohesion: p.coh, frictionDeg: p.phi, gfShearFactor: p.gfs };
  R.querySelectorAll("input[data-p]").forEach((i) => { if (document.activeElement !== i) i.value = fr(vals[i.dataset.p] / +i.dataset.f, +i.dataset.d); });
  R.querySelector('input[data-p="yan_c"]').disabled = !L.yan;
  R.querySelectorAll("#l-montee button").forEach((b) => b.classList.toggle("actif", b.dataset.v === p.montee));
  R.querySelectorAll("#l-adou button").forEach((b) => b.classList.toggle("actif", b.dataset.v === p.adoucissement));
  const lcz = (q) => q.E * q.GfI / q.ft ** 2, mm = (v) => (v > 1 ? `${fr(v, 1)} m` : `${fr(v * 1e3, 2)} mm`);
  const lignes = [["Ouverture au pic δe", `${fr(L.dnE * 1e9, 1)} nm`, `${fr(Lr.dnE * 1e9, 1)} nm`],
    ["Ouverture de rupture δF", `${fr((L.dnE + L.ot) * 1e6, 2)} µm`, `${fr((Lr.dnE + Lr.ot) * 1e6, 2)} µm`],
    ["Intégrale I de f(D)", fr(L.I, 4), fr(Lr.I, 4)],
    ["Glissement de rupture", `${fr(L.st * 1e6, 1)} µm`, `${fr(Lr.st * 1e6, 1)} µm`],
    ["Zone cohésive E Gf / ft²", mm(lcz(p)), mm(lcz(pr))],
    ["Rapport à l'élément (> 2)", fr(lcz(p) / p.h, 1), fr(lcz(pr) / pr.h, 1)],
    ["Pic τ à σn", `${fr((p.coh + L.tanPhi * SN) / 1e6, 1)} MPa`, `${fr((pr.coh + Lr.tanPhi * SN) / 1e6, 1)} MPa`]];
  R.querySelector("#l-lect").innerHTML = `<tr><th></th><th>Réglage</th><th>Préréglage</th></tr>` +
    lignes.map((l) => `<tr><td>${l[0]}</td><td>${l[1]}</td><td style="color:var(--texte3)">${l[2]}</td></tr>`).join("") +
    (lcz(p) / p.h < 2 ? `<tr><td colspan="3" style="color:var(--attention);font-family:var(--police)">Zone cohésive sous 2 éléments : la fissure n'est plus résolue par le maillage.</td></tr>` : "");
  const c = lireChoix(), mods = [...Object.entries(c.surcharges || {}), ...Object.entries(c.loi || {})];
  const nomCle = { penalite: "insertionPenaltyFactor", yan_c: "yanC", montee: "jointElastic", adoucissement: "jointSoftening" };
  R.querySelector("#l-cles").innerHTML = mods.length ? mods.map(([k, v]) => `<b>${nomCle[k] || k} = ${typeof v === "number" ? +v.toPrecision(5) : v}</b>`).join("<br>") +
    `<br><span style="color:var(--texte3)">appliquées à l'essai en cours (écran Essai)</span>` : "aucune : loi du préréglage";
}

// ------------------------------------------------------------------ cliquer-glisser
function attacher(sel, noms) {
  const cv = R.querySelector(sel);
  const local = (ev) => { const r = cv.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; };
  const proche = (q) => noms.find((n) => POIG[n] && (n !== "forme" || params().adoucissement === "yan") && Math.hypot(POIG[n][0] - q[0], POIG[n][1] - q[1]) < 12);
  cv.addEventListener("pointerdown", (ev) => { glisse = proche(local(ev)) || null; if (glisse) { cv.setPointerCapture(ev.pointerId); dessiner(); } });
  cv.addEventListener("pointermove", (ev) => {
    const q = local(ev);
    if (!glisse) { cv.style.cursor = proche(q) ? "grab" : "default"; return; }
    cv.style.cursor = "grabbing";
    const T = POIG[glisse][2], p = params(), L = loiDe(p), pr = params(false);
    if (glisse === "pic") changer({ surcharges: { ft: Math.max(0.01 * pr.ft, T.Yi(q[1]) * 1e6) } });
    else if (glisse === "fin") changer({ surcharges: { Gf: Math.max(0.02 * pr.GfI, (T.Xi(q[0]) * 1e-6 - L.dnE) * p.ft * L.I) } });
    else if (glisse === "montee") changer({ loi: { penalite: Math.max(0.05, p.ft / Math.max(1e-12, T.Xi(q[0]) * 1e-9) * p.h / p.E) } });
    else if (glisse === "forme") {
      const cible = Math.min(0.97, Math.max(0.01, T.Yi(q[1]) * 1e6 / p.ft));
      let best = p.c, e = Infinity;
      for (let cc = 0.3; cc <= 40; cc += 0.05) { const v = Math.abs(fD(0.5, 0.63, 1.8, cc) - cible); if (v < e) { e = v; best = cc; } }
      changer({ loi: { yan_c: +best.toFixed(2) } });
    } else if (glisse === "coh") changer({ surcharges: { cohesion: Math.max(0.01 * pr.coh, T.Yi(q[1]) * 1e6 - L.tanPhi * SN) } });
    else if (glisse === "finII") {
      const sE = (p.coh + L.tanPhi * SN) / p.pj;
      changer({ surcharges: { gfShearFactor: Math.max(0.05, (T.Xi(q[0]) * 1e-6 - sE) * p.coh * L.I / p.GfI) } });
    }
  });
  cv.addEventListener("pointerup", () => { glisse = null; dessiner(); });
}

// ------------------------------------------------------------------ essai éclair
const NOMS_CH = { traction: "Traction", ucs: "Compression simple", tx20: "Triaxial 20 MPa" };
async function eclair(ch) {
  const msg = R.querySelector("#l-msg"), c = lireChoix();
  const libelle = [...Object.entries(c.surcharges || {}), ...Object.entries(c.loi || {})].map(([k, v]) => `${k} ${typeof v === "number" ? +v.toPrecision(3) : v}`).join(", ") || "préréglage";
  chargementVu = ch;
  try {
    const t0 = performance.now();
    let e = await api("/api/court/eclair", { choix: c, chargement: ch });
    while (e.etat === "en_cours" || e.etat === "absent") {
      msg.textContent = `${NOMS_CH[ch]} : g1 calcule… ${fr((performance.now() - t0) / 1000, 0)} s`; msg.className = "message";
      await enAttente(300);
      e = await api(`/api/court/${e.cle}`);
    }
    if (e.etat !== "fini") throw new Error(e.etat + "\n" + (e.journal || []).join("\n"));
    const h = await api(`/api/runs/${encodeURIComponent(e.id)}/historique`);
    const dejà = eclairs.find((x) => x.cle === e.cle);
    if (!dejà) eclairs.push({ cle: e.cle, ch, libelle, h, s: e.synthese, duree: (performance.now() - t0) / 1000 });
    msg.textContent = dejà ? "Ce réglage était déjà calculé : courbe existante." : `${NOMS_CH[ch]} : fini en ${fr((performance.now() - t0) / 1000, 1)} s.`;
    msg.className = "message ok";
    tracerEclairs();
  } catch (err) { msg.textContent = "Échec : " + err.message; msg.className = "message err"; }
}

function tracerEclairs() {
  const l = eclairs.filter((x) => x.ch === chargementVu);
  courbeEcl.definir(l.map((x, i) => ({ x: x.h.eps, y: x.h.q, couleur: SERIES[i % 6], nom: x.libelle.slice(0, 28) })));
  R.querySelector("#l-ecl-tab").innerHTML = l.length ? `<tr><th>${NOMS_CH[chargementVu]}</th><th>Pic filtré</th><th>ε pic</th></tr>` +
    l.map((x, i) => `<tr><td><span class="trait" style="display:inline-block;margin-right:6px;background:${SERIES[i % 6]}"></span>${x.libelle}</td><td>${fr(x.s.q_pic_filtre_MPa, chargementVu === "traction" ? 2 : 1)} MPa</td><td>${fr(x.s.eps_pic_filtre_pct, 3)} %</td></tr>`).join("")
    : `<tr><td style="color:var(--texte3);padding:16px">Chaque clic lance g1 sur le banc de 10 grains avec la loi affichée, et superpose la courbe aux réglages précédents.</td></tr>`;
}
