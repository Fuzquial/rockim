// Écran Essai : choix physiques -> essai (noyau/formulaire.py), schéma, vérifications, deck, ajout à la file.
import { api, css, fr } from "./commun.js";
import { abonner, initialiser, lireChoix, lireFormulaire, modifier, reinitialiser } from "./etat_essai.js";

let R, nav;

// ------------------------------------------------------------------ briques de formulaire
const seg = (cle, opts) => `<div class="segments plein" data-cle="${cle}">${opts.map(([v, l]) => `<button data-v="${v}">${l}</button>`).join("")}</div>`;
const num = (cle, unite, f = 1, d = 2) => `<div class="nombre"><input data-cle="${cle}" data-f="${f}" data-d="${d}"><span>${unite}</span></div>`;
const ligne = (label, ctrl, si = "") => `<div class="ligne" ${si ? `data-si="${si}"` : ""}><label>${label}</label>${ctrl}</div>`;
const inter = (cle, label, si = "") => `<div class="interrupteur" data-cle="${cle}" ${si ? `data-si="${si}"` : ""}><i></i>${label}</div>`;
const aide = (id, si = "") => `<div class="aide-ligne" id="${id}" ${si ? `data-si="${si}"` : ""}></div>`;

function condition(expr, c) {
  return expr.split("&").every((t) => {
    const [k, op, v] = t.match(/^([\w]+)(!=|=)?(.*)$/).slice(1);
    if (!op) return !!c[k];
    return op === "=" ? String(c[k]) === v : String(c[k]) !== v;
  });
}

export async function monter(noeud, ctx) {
  R = noeud; nav = ctx;
  await initialiser();
  R.innerHTML = `
  <div class="formulaire">
    <section class="carte"><div class="carte-tete"><span class="titre">Essai</span><div style="display:flex;gap:6px"><button class="bouton petit" id="e-lien" title="Copie un lien qui recrée cet essai">Lien</button><button class="bouton petit" id="e-raz">Repartir du défaut</button></div></div>
      <div class="corps-carte">
        ${ligne("Nom", `<div class="nombre"><input data-cle="nom" data-texte="1" style="text-align:left"></div>`)}
        ${ligne("Préréglage d'essai", `<select id="e-preset"><option value="">—</option></select>`)}
        ${seg("type_essai", [["triaxial", "Triaxial"], ["bresilien", "Brésilien"], ["traction", "Traction directe"]])}
        ${ligne("Confinement σ₃", num("sigma3_MPa", "MPa", 1, 1), "type_essai=triaxial")}
        ${ligne("Vitesse", num("vitesse", "m/s", 1, 2))}
        ${ligne("Arrêt après le pic", num("chute_arret", "% de chute", 100, 0), "type_essai=triaxial")}
        ${aide("e-aide-charge", "type_essai!=traction")}
      </div></section>
    <section class="carte"><div class="carte-tete"><span class="titre">Éprouvette et matériau</span></div>
      <div class="corps-carte">
        ${ligne("Largeur × hauteur", `<div class="deux">${num("W_mm", "mm", 1, 1)}${num("H_mm", "mm", 1, 1)}</div>`, "type_essai!=bresilien")}
        ${ligne("Diamètre du disque", num("D_mm", "mm", 1, 1), "type_essai=bresilien")}
        ${ligne("Méplat · plateau", `<div class="deux">${num("aplatissement_deg", "° total", 1, 0)}${num("plateau_mm", "mm ½", 1, 2)}</div>`, "type_essai=bresilien")}
        ${ligne("Préréglage", `<select data-cle="materiau"><option value="fragile">Red Bohus, fragile (Yan 2023)</option><option value="bohus">Red Bohus, calibré</option></select>`)}
        ${aide("e-aide-mat")}
        <div class="loi-modifiee" id="e-loi" hidden></div>
      </div></section>
    <section class="carte"><div class="carte-tete"><span class="titre">Microstructure</span></div>
      <div class="corps-carte">
        ${ligne("Maillage", seg("maillage", [["voronoi", "Voronoï"], ["gmsh", "Gmsh homogène"]]))}
        ${ligne("Phases", seg("phases", [["une", "Une"], ["trois", "Trois minéraux"]]), "maillage=voronoi")}
        ${ligne("Taille de grain", num("taille_grain_mm", "mm", 1, 2), "maillage=voronoi")}
        ${ligne("Taille d'élément", num("taille_element_mm", "mm", 1, 2))}
        ${aide("e-aide-maillage")}
        ${ligne("Tailles de grain", seg("tailles", [["uniformes", "Uniformes"], ["dispersees", "Dispersées"], ["par_phase", "Par phase"]]), "maillage=voronoi")}
        ${ligne("Dispersion", num("dispersion", "—", 1, 2), "maillage=voronoi&tailles!=uniformes")}
        ${ligne("Contraste minéral", seg("contraste", [["elastique", "Élastique seul"], ["complet", "Complet"]]), "maillage=voronoi&phases=trois")}
        <div data-si="maillage=voronoi&phases=trois"><table class="tab-phases" id="e-phases"></table></div>
        ${ligne("Joints de grain", `<div class="segments plein" data-cle="joints_grain"><button data-v="niveau" id="e-jg-niveau">Niveau</button><button data-v="identiques">Identiques</button><button data-v="paires" data-si="phases=trois">Par paire</button><button data-v="inter_intra">Inter/intra</button></div>`, "maillage=voronoi")}
        <div class="deux-alphas" data-si="maillage=voronoi&joints_grain=inter_intra" id="e-alphas"></div>
      </div></section>
    <section class="carte"><div class="carte-tete"><span class="titre">Discontinuités</span><span class="overline">combinables</span></div>
      <div class="corps-carte">
        ${inter("diffus", "Joints pré-rompus diffus")}
        ${ligne("Fraction", num("fraction_diffuse", "%", 100, 1), "diffus")}
        ${inter("plans", "Famille de plans affaiblis")}
        ${ligne("Pendage · espacement", `<div class="deux">${num("pendage", "°", 1, 0)}${num("espacement_mm", "mm", 1, 1)}</div>`, "plans")}
        ${ligne("Facteur · part rompue", `<div class="deux">${num("facteur_plans", "× ft", 1, 2)}${num("fraction_rompue", "%", 100, 0)}</div>`, "plans")}
        ${aide("e-aide-segments")}
        ${aide("e-aide-mu")}
      </div></section>
    <section class="carte"><div class="carte-tete"><span class="titre">Calcul</span></div>
      <div class="corps-carte">
        ${ligne("Durée simulée max.", num("T", "ms", 1000, 2))}
        ${ligne("Frames écrites", num("frames", "", 1, 0))}
        ${ligne("Graine", num("graine", "", 1, 0))}
      </div></section>
  </div>
  <section class="carte schema"><div class="carte-tete"><span class="titre">Schéma de l'essai</span><span class="overline">à l'échelle</span></div>
    <svg id="e-schema" viewBox="0 0 520 640"></svg></section>
  <div class="colonne-droite">
    <section class="carte"><div class="carte-tete"><span class="titre">Estimation</span><span class="overline">avant tout calcul</span></div>
      <div class="chiffres" id="e-chiffres"></div></section>
    <section class="carte"><div class="carte-tete"><span class="titre">Vérifications</span><span class="overline" id="e-bilan"></span></div>
      <ul class="avis" id="e-avis"></ul>
      <div class="boutons">
        <button class="bouton" id="e-apercu">Aperçu du maillage</button>
        <button class="bouton primaire" id="e-ajouter">Ajouter à la file</button>
        <div class="ligne" data-si="type_essai=triaxial" style="grid-template-columns:1fr auto"><div class="nombre"><input id="e-variation" value="0 10 20 40" style="text-align:left"><span>σ₃ MPa</span></div>
          <button class="bouton" id="e-varier" style="padding:0 12px">Varier σ₃</button></div>
        <div class="message" id="e-message"></div>
      </div></section>
    <section class="carte"><div class="carte-tete"><span class="titre">Deck g1</span><span class="overline">généré, lecture seule</span></div>
      <div class="deck" id="e-deck"></div></section>
  </div>`;
  brancher();
  abonner(afficher);
}

export async function activer() {
  const { lireResultat } = await import("./etat_essai.js");
  if (lireResultat()) afficher(lireResultat(), lireChoix());
}

// ------------------------------------------------------------------ saisie
function brancher() {
  R.querySelectorAll(".segments[data-cle]").forEach((s) => s.querySelectorAll("button").forEach((b) => (b.onclick = () => modifier({ [s.dataset.cle]: b.dataset.v }))));
  R.querySelectorAll(".interrupteur[data-cle]").forEach((t) => (t.onclick = () => modifier({ [t.dataset.cle]: !lireChoix()[t.dataset.cle] })));
  R.querySelectorAll("select[data-cle]").forEach((s) => (s.onchange = () => modifier({ [s.dataset.cle]: s.value })));
  R.querySelectorAll("input[data-cle]").forEach((i) => (i.onchange = () => {
    if (i.dataset.texte) return modifier({ [i.dataset.cle]: i.value.trim().replace(/[^\w.-]+/g, "_") || "essai" });
    const v = parseFloat(i.value.replace(/\s/g, "").replace(",", "."));
    if (Number.isFinite(v)) modifier({ [i.dataset.cle]: v / +i.dataset.f });
    else if (i.dataset.cle === "T" || i.dataset.cle === "vitesse") modifier({ [i.dataset.cle]: null });
  }));
  R.querySelector("#e-raz").onclick = () => reinitialiser();
  const pr = R.querySelector("#e-preset"), presets = lireFormulaire().presets || {};
  pr.innerHTML += Object.entries(presets).map(([k, v]) => `<option value="${k}">${v.description}</option>`).join("");
  pr.onchange = () => { if (pr.value) { modifier({ ...lireFormulaire().defaut, ...presets[pr.value].choix }, true); pr.value = ""; } };
  R.querySelector("#e-lien").onclick = async () => {
    const { lienDePartage } = await import("./etat_essai.js");
    const m = R.querySelector("#e-message");
    try { await navigator.clipboard.writeText(lienDePartage()); m.textContent = "Lien de l'essai copié."; m.className = "message ok"; }
    catch (e) { m.textContent = lienDePartage(); m.className = "message"; }
  };
  R.querySelector("#e-apercu").onclick = () => nav.ouvrir("maillage");
  R.querySelector("#e-ajouter").onclick = () => ajouter();
  R.querySelector("#e-varier").onclick = () => {
    const l = R.querySelector("#e-variation").value.split(/[\s;,]+/).map(Number).filter((x) => Number.isFinite(x) && x >= 0);
    if (l.length) ajouter(l);
  };
}

async function ajouter(sigma3) {
  const r = (await import("./etat_essai.js")).lireResultat(), m = R.querySelector("#e-message");
  try {
    const rep = await api("/api/essai/ajouter", { essai: r.essai, sigma3_liste: sigma3 || null, lot: sigma3 ? `Variation σ₃ · ${r.essai.nom}` : null });
    m.innerHTML = `${rep.ajoutes.length} run${rep.ajoutes.length > 1 ? "s" : ""} ajouté${rep.ajoutes.length > 1 ? "s" : ""} à la file. <a href="#file" id="e-vers-file">Voir la file →</a>`;
    m.className = "message ok";
    R.querySelector("#e-vers-file").onclick = (ev) => { ev.preventDefault(); nav.ouvrir("file"); };
  } catch (e) { m.textContent = e.message; m.className = "message err"; }
}

// ------------------------------------------------------------------ affichage
function synchroniser(c) {
  R.querySelectorAll("[data-si]").forEach((x) => (x.hidden = !condition(x.dataset.si, c)));
  R.querySelectorAll(".segments[data-cle]").forEach((s) => s.querySelectorAll("button").forEach((b) => b.classList.toggle("actif", String(c[s.dataset.cle]) === b.dataset.v)));
  R.querySelectorAll(".interrupteur[data-cle]").forEach((t) => t.classList.toggle("on", !!c[t.dataset.cle]));
  R.querySelectorAll("select[data-cle]").forEach((s) => (s.value = c[s.dataset.cle]));
  R.querySelectorAll("input[data-cle]").forEach((i) => {
    if (document.activeElement === i) return;
    const v = c[i.dataset.cle];
    i.value = i.dataset.texte ? v : v == null ? "" : fr(v * +i.dataset.f, +i.dataset.d);
  });
}

function afficher(r, c) {
  synchroniser(c);
  // champs laissés vides : on montre la valeur effective, en gris
  for (const [cle, v, f, d] of [["vitesse", r.essai.chargement.vitesse, 1, 2], ["T", r.essai.sorties.T, 1000, 2]]) {
    const i = R.querySelector(`input[data-cle="${cle}"]`);
    if (c[cle] == null && document.activeElement !== i) { i.value = fr(v * f, d); i.classList.add("implicite"); }
    else i.classList.remove("implicite");
  }
  const f = lireFormulaire(), niv = f.niveaux[c.materiau], m = niv.materiau, e = r.essai, est = r.estimation;
  const bn = R.querySelector("#e-jg-niveau");
  bn.textContent = `Niveau ${fr(niv.alpha, 1)}`;
  bn.title = `Atténuation des joints de grain du préréglage : gbAlpha = ${fr(niv.alpha, 1)}`;
  R.querySelector("#e-aide-mat").textContent = `E ${fr(m.E / 1e9, 0)} GPa · ν ${fr(m.nu, 2)} · ft ${fr(m.ft / 1e6, 1)} MPa · c ${fr(m.cohesion / 1e6, 1)} MPa · φ ${fr(m.frictionDeg, 1)}° · Gf ${fr(m.Gf, 1)} J/m²`;
  R.querySelector("#e-aide-charge").textContent = c.type_essai === "bresilien"
    ? `Arrêt ${fr(e.chargement.delai_arret_bresilien * 1e6, 0)} µs après la chute post-pic ; jauge élastique lue pour σt entre ${fr(e.chargement.jauge_elastique[0], 1)} et ${fr(e.chargement.jauge_elastique[1], 1)} ft.`
    : `Confinement établi en ${fr(e.chargement.rampe_confinement * 1e3, 1)} ms, charge axiale à partir de ${fr(e.chargement.delai_axial * 1e3, 1)} ms.`;
  const mod = Object.keys(c.surcharges || {}).length + Object.keys(c.loi || {}).length;
  const zl = R.querySelector("#e-loi");
  zl.hidden = !mod;
  if (mod) {
    zl.innerHTML = `Loi des joints modifiée : ${[...Object.entries(c.surcharges), ...Object.entries(c.loi)].map(([k, v]) => `${k} = ${typeof v === "number" ? v.toPrecision(4) : v}`).join(" · ")}
      <button class="lien" id="e-loi-ed">modifier</button><button class="lien" id="e-loi-raz">revenir au préréglage</button>`;
    R.querySelector("#e-loi-ed").onclick = () => nav.ouvrir("loi");
    R.querySelector("#e-loi-raz").onclick = () => modifier({ surcharges: {}, loi: {} }, true);
  }
  // phases
  R.querySelector("#e-phases").innerHTML = `<tr><th>Phase</th><th>Aire</th><th>E GPa</th><th>ft MPa</th><th>c MPa</th><th>φ °</th></tr>` +
    niv.phases.map((p, i) => `<tr><td><i class="pastille-phase" style="background:${[css("--quartz"), css("--feldspath"), css("--biotite")][i]}"></i>${{ quartz: "Quartz", feldspar: "Feldspath", biotite: "Biotite" }[p.nom] || p.nom}</td>
      <td>${fr(100 * p.fraction, 1)} %</td><td>${fr(p.E / 1e9, 1)}</td>${c.contraste === "complet" ? `<td>${fr(p.ft / 1e6, 2)}</td><td>${fr(p.cohesion / 1e6, 1)}</td><td>${fr(p.frictionDeg, 1)}</td>` : `<td colspan="3" style="text-align:center;color:var(--texte3)">ceux du volume</td>`}</tr>`).join("");
  const al = R.querySelector("#e-alphas");
  if (!al.dataset.fait) {
    al.dataset.fait = 1;
    al.innerHTML = `<div class="aide-ligne" style="margin:0 0 4px">Rapport joint de grain / joint intragranulaire (Table 2 du papier granite par défaut)</div>` +
      ["Ten", "Coh", "Gf", "E", "Fric"].map((k) => ligne({ Ten: "Traction", Coh: "Cohésion", Gf: "Énergie", E: "Raideur", Fric: "Frottement" }[k], `<div class="nombre"><input data-alpha="${k}"><span>× intra</span></div>`)).join("");
    al.querySelectorAll("input").forEach((i) => (i.onchange = () => {
      const v = parseFloat(i.value.replace(",", "."));
      if (Number.isFinite(v) && v > 0) modifier({ alphas: { ...lireChoix().alphas, [i.dataset.alpha]: v } });
    }));
  }
  al.querySelectorAll("input").forEach((i) => { if (document.activeElement !== i) i.value = fr(c.alphas[i.dataset.alpha], 2); });
  // aides
  const epg = est.elements_par_grain;
  const am = R.querySelector("#e-aide-maillage");
  am.textContent = c.maillage === "voronoi" ? `≈ ${fr(epg, 0)} éléments par grain (plage 35-90), ${est.grains} grains` : `≈ ${est.elements.toLocaleString("fr-FR")} triangles Gmsh`;
  am.className = "aide-ligne " + (c.maillage !== "voronoi" || (epg >= 35 && epg <= 90) ? "ok" : "att");
  const segs = r.segments_prerompus.length, dess = (c.segments || []).length;
  R.querySelector("#e-aide-segments").textContent = c.plans || dess ? `${segs} segment${segs > 1 ? "s" : ""} pré-rompu${segs > 1 ? "s" : ""} (${dess} dessiné${dess > 1 ? "s" : ""} sur l'aperçu)` : "";
  R.querySelector("#e-aide-mu").textContent = c.diffus || c.plans || dess ? `Frottement résiduel des pré-rompus : tan ${fr(e.materiau.frictionDeg, 1)}° = ${fr(Math.tan(e.materiau.frictionDeg * Math.PI / 180), 3)} (frottement de pic)` : "";
  // estimation
  const duree = (s) => s < 90 ? `${fr(s, 0)} s` : s < 5400 ? `${fr(s / 60, 0)} min` : `${fr(s / 3600, 1)} h`;
  R.querySelector("#e-chiffres").innerHTML = `
    <div><small>Éléments</small><b>${est.elements.toLocaleString("fr-FR")}</b></div>
    <div><small>${c.maillage === "voronoi" ? "Grains" : "Maillage"}</small><b>${c.maillage === "voronoi" ? est.grains : "Gmsh"}</b></div>
    <div><small>Pas de temps</small><b>${fr(est.dt * 1e9, 1)} ns</b></div>
    <div><small>Durée au plus</small><b>≈ ${duree(est.duree_max_s)}</b></div>`;
  // vérifications
  const ordre = { erreur: 0, alerte: 1, info: 2 }, cls = { erreur: "err", alerte: "att", info: "info" };
  const av = [...r.avis].sort((a, b) => ordre[a.niveau] - ordre[b.niveau]);
  const ne = av.filter((a) => a.niveau === "erreur").length, na = av.filter((a) => a.niveau === "alerte").length;
  R.querySelector("#e-bilan").textContent = `${ne} erreur${ne > 1 ? "s" : ""} · ${na} alerte${na > 1 ? "s" : ""}`;
  const recette = { triaxial: "triaxiale : plateaux, confinement latéral, arrêt après le pic",
    traction: "de traction directe par mors", bresilien: "brésilienne : disque à méplats entre plateaux, jauge élastique de bande, arrêt après le pic" };
  R.querySelector("#e-avis").innerHTML = (ne + na ? "" : `<li class="ok">Recette ${recette[c.type_essai]}.</li>`) +
    av.map((a) => `<li class="${cls[a.niveau]}">${a.message}</li>`).join("");
  R.querySelector("#e-ajouter").disabled = R.querySelector("#e-varier").disabled = ne > 0;
  R.querySelector("#e-apercu").disabled = ne > 0 || c.maillage !== "voronoi";
  // deck
  R.querySelector("#e-deck").innerHTML = r.deck.split("\n").map((l) => l.startsWith("#") ? `<span class="c">${l}</span>`
    : l.includes("=") ? `<span class="k">${l.split("=")[0]}</span>=${l.split("=").slice(1).join("=")}` : l).join("\n");
  schema(c, r);
}

// ------------------------------------------------------------------ schéma à l'échelle
function schema(c, r) {
  const svg = R.querySelector("#e-schema"), bd = c.type_essai === "bresilien";
  const Wm = bd ? c.D_mm : c.W_mm, Hm = bd ? c.D_mm : c.H_mm;
  const k = Math.min(300 / Wm, 420 / Hm), W = Wm * k, H = Hm * k, x0 = 260 - W / 2, y0 = 320 - H / 2;
  const X = (mm) => x0 + mm * k, Y = (mm) => y0 + H - mm * k;            // repère du solveur : origine en bas à gauche
  const acc = css("--accent");
  // disque à méplats : demi-angle alpha, cordes horizontales en haut et en bas
  function disque(attr = "") {
    const rr = W / 2, a = (c.aplatissement_deg / 2) * Math.PI / 180, cx = x0 + rr, cy = y0 + rr;
    if (a <= 0) return `<circle cx="${cx}" cy="${cy}" r="${rr}" ${attr}/>`;
    const sx = rr * Math.sin(a), sy = rr * Math.cos(a);
    return `<path d="M${cx - sx},${cy - sy} L${cx + sx},${cy - sy} A${rr},${rr} 0 0 1 ${cx + sx},${cy + sy} L${cx - sx},${cy + sy} A${rr},${rr} 0 0 1 ${cx - sx},${cy - sy} Z" ${attr}/>`;
  }
  const fl = (x1, y1, x2, y2) => `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${acc}" stroke-width="2" marker-end="url(#pointe)"/>`;
  let g = `<defs><marker id="pointe" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="${acc}"/></marker>
    <clipPath id="eprouvette">${bd ? disque() : `<rect x="${x0}" y="${y0}" width="${W}" height="${H}"/>`}</clipPath></defs>`;
  g += bd ? disque(`fill="${css("--bg3")}" stroke="${css("--texte2")}"`)
    : `<rect x="${x0}" y="${y0}" width="${W}" height="${H}" fill="${css("--bg3")}" stroke="${css("--texte2")}"/>`;
  // grains suggérés (pas le maillage réel : celui-ci est dans l'écran Maillage)
  if (c.maillage === "voronoi") {
    const d = c.taille_grain_mm * k;
    g += `<g clip-path="url(#eprouvette)" stroke="${css("--border-fort")}" fill="none">`;
    let a = 11;
    for (let yy = y0 - d; yy < y0 + H + d; yy += d * 0.87) for (let xx = x0 - d + ((yy / d) % 2) * d / 2; xx < x0 + W + d; xx += d) {
      a = (a * 9301 + 49297) % 233280; const j = (a / 233280 - .5) * d * .5;
      g += `<polygon points="${[0, 1, 2, 3, 4, 5].map((i) => `${xx + j + d * .55 * Math.cos(i * Math.PI / 3)},${yy + d * .55 * Math.sin(i * Math.PI / 3)}`).join(" ")}"/>`;
    }
    g += `</g>`;
  }
  // famille de plans (convention du solveur : n . X = k s, ancrée à l'origine)
  if (c.plans) {
    const b = c.pendage * Math.PI / 180, nx = -Math.sin(b), ny = Math.cos(b), s = c.espacement_mm;
    const kmax = Math.ceil((Math.abs(nx) * Wm + Math.abs(ny) * Hm) / s) + 2;
    g += `<g clip-path="url(#eprouvette)" stroke="${css("--attention")}" stroke-width="1.4" stroke-dasharray="6 4">`;
    for (let i = -kmax; i <= kmax; i++) {
      const t = 2 * (Wm + Hm), px = nx * i * s, py = ny * i * s;
      g += `<line x1="${X(px - ny * t)}" y1="${Y(py + nx * t)}" x2="${X(px + ny * t)}" y2="${Y(py - nx * t)}"/>`;
    }
    g += `</g>`;
  }
  for (const [x1, y1, x2, y2] of r.segments_prerompus)
    g += `<line x1="${X(x1 * 1e3)}" y1="${Y(y1 * 1e3)}" x2="${X(x2 * 1e3)}" y2="${Y(y2 * 1e3)}" stroke="${css("--pre")}" stroke-width="3"/>`;
  if (c.diffus) {
    let a = 7;
    const n = Math.round(Math.min(160, 25 + 900 * c.fraction_diffuse));
    for (let i = 0; i < n; i++) {
      a = (a * 9301 + 49297) % 233280; const u = a / 233280; a = (a * 9301 + 49297) % 233280; const v = a / 233280;
      a = (a * 9301 + 49297) % 233280; const t = a / 233280 * Math.PI;
      const cx = x0 + 6 + u * (W - 12), cy = y0 + 6 + v * (H - 12);
      g += `<line x1="${cx}" y1="${cy}" x2="${cx + 6 * Math.cos(t)}" y2="${cy + 6 * Math.sin(t)}" stroke="${css("--pre")}" stroke-width="2"/>`;
    }
  }
  // chargement
  if (bd) {
    const a = (c.aplatissement_deg / 2) * Math.PI / 180, rr = W / 2, yh = y0 + rr - rr * Math.cos(a), yb = y0 + rr + rr * Math.cos(a), lp = c.plateau_mm * k;
    g += `<rect x="${260 - lp}" y="${yh - 12}" width="${2 * lp}" height="12" rx="2" fill="${css("--border-fort")}"/>`;
    g += `<rect x="${260 - lp}" y="${yb}" width="${2 * lp}" height="12" rx="2" fill="${css("--border-fort")}"/>`;
    g += fl(260, yh - 56, 260, yh - 16) + fl(260, yb + 52, 260, yb + 16);
    g += `<text x="274" y="${yh - 30}" fill="${css("--texte2")}" font-size="13">plateaux, ${fr(r.essai.chargement.vitesse, 2)} m/s au total</text>`;
    g += `<text x="${x0 - 10}" y="${y0 + rr}" fill="${css("--texte3")}" font-size="12" text-anchor="end">σt = 2P / (π D t)</text>`;
  } else if (c.type_essai === "triaxial") {
    for (const [y, sens] of [[y0 - 16, 1], [y0 + H, -1]]) {
      g += `<rect x="${x0 - 20}" y="${y}" width="${W + 40}" height="16" rx="3" fill="${css("--border-fort")}"/>`;
      g += fl(260, sens > 0 ? y - 44 : y + 60, 260, sens > 0 ? y - 4 : y + 20);
    }
    g += `<text x="274" y="${y0 - 40}" fill="${css("--texte2")}" font-size="13">plateaux, ${fr(r.essai.chargement.vitesse, 2)} m/s</text>`;
    if (c.sigma3_MPa > 0) {
      for (let i = 0; i < 7; i++) { const y = y0 + (i + 0.5) * H / 7; g += fl(x0 - 58, y, x0 - 6, y) + fl(x0 + W + 58, y, x0 + W + 6, y); }
      g += `<text x="${x0 - 60}" y="${y0 - 6}" fill="${css("--texte2")}" font-size="13" text-anchor="end">σ₃ = ${fr(c.sigma3_MPa, 0)} MPa</text>`;
    }
  } else {
    for (const [y, sens] of [[y0 - 16, -1], [y0 + H, 1]]) {
      g += `<rect x="${x0}" y="${y}" width="${W}" height="16" fill="${css("--border-fort")}"/>`;
      g += fl(260, sens < 0 ? y - 4 : y + 20, 260, sens < 0 ? y - 44 : y + 60);
    }
    g += `<text x="274" y="${y0 - 40}" fill="${css("--texte2")}" font-size="13">mors, ${fr(r.essai.chargement.vitesse, 2)} m/s</text>`;
  }
  if (!bd) g += `<text x="${x0 + W / 2 + 60}" y="${y0 + H + 40}" fill="${css("--texte3")}" font-size="12" text-anchor="middle" font-family="monospace">${fr(Wm, 0)} mm</text>`;
  if (bd) g += `<text x="260" y="${y0 + H + 90}" fill="${css("--texte3")}" font-size="12" text-anchor="middle" font-family="monospace">D = ${fr(Wm, 1)} mm</text>`;
  else g += `<text x="${x0 + W + 74}" y="320" fill="${css("--texte3")}" font-size="12" font-family="monospace" transform="rotate(90 ${x0 + W + 74} 320)" text-anchor="middle">${fr(Hm, 0)} mm</text>`;
  svg.innerHTML = g;
}
