// Écran Résultats : comparaison de runs (courbes + synthèse) et détail d'un run (champs, fissures).
import { api, css, el, enAttente, fr, SERIES } from "./commun.js";
import { Courbe } from "./courbe.js";
import { BANDES, CHAMPS, couleurArc, VueChamps } from "./vue_champs.js";

let R, nav, runs = [], coches = [], courbeComp, syntheses = {}, hist = {};
let D = null;                                         // état du mode détail

const COLONNES = [
  ["q_pic_filtre_MPa", "q pic MPa", 1], ["q_pic_MPa", "max brut", 1], ["eps_pic_filtre_pct", "ε pic %", 3],
  ["E_secante_GPa", "E séc. GPa", 1], ["sigma3_atteint_MPa", "σ₃ atteint", 1], ["chute_post_pic", "chute", 2],
  ["part_intergranulaire", "intergr.", 2], ["n_tension", "traction", 0], ["n_cisaillement", "cisaill.", 0],
  ["prerompus_libres", "pré-r. libres", 0],
];

export async function monter(noeud, ctx) {
  R = noeud; nav = ctx;
  // Lien direct vers une comparaison : #resultats?c=id1,id2
  const c = new URLSearchParams(location.hash.split("?")[1] || "").get("c");
  if (c) coches = c.split(",");
  R.innerHTML = `
    <section class="carte liste-runs"><div class="carte-tete"><span class="titre">Runs</span><span class="overline">cocher pour comparer</span></div>
      <div id="r-liste" style="margin-top:6px"></div></section>
    <div id="r-droite" style="min-height:0;display:grid"></div>`;
  await rafraichirListe();
  setInterval(() => { if (R.classList.contains("actif")) rafraichirListe(); }, 5000);
}

export function activer(params = {}) {
  if (params.run) ouvrirDetail(params.run);
  else if (!D) vueComparaison();
}

// ------------------------------------------------------------------ liste
async function rafraichirListe() {
  const n = await api("/api/runs");
  const cle = (l) => l.map((r) => r.id + r.etat).join();
  if (cle(n) === cle(runs)) return;
  runs = n;
  const groupes = {};
  for (const r of runs) (groupes[r.source] ||= []).push(r);
  R.querySelector("#r-liste").innerHTML = Object.entries(groupes).map(([src, l]) => `
    <div class="overline" style="padding:8px 12px 4px">${src === "espace" ? "Espace de travail" : src}</div>
    ${l.map((r) => `<label data-id="${r.id}">
        <input type="checkbox" ${coches.includes(r.id) ? "checked" : ""} ${r.etat !== "fini" && r.etat !== "deja_fait" ? "disabled" : ""}>
        <span class="trait" style="background:${coches.includes(r.id) ? SERIES[coches.indexOf(r.id) % 6] : "transparent"}"></span>
        <a class="lien-run" href="#resultats/${encodeURIComponent(r.id)}">${r.nom}</a>
        <span class="mono2">${r.etat === "fini" || r.etat === "deja_fait" ? "" : r.etat.replace("_", " ")}</span></label>`).join("")}`).join("");
  R.querySelectorAll("#r-liste input").forEach((cb) => (cb.onchange = () => {
    const id = cb.closest("label").dataset.id;
    coches = cb.checked ? [...coches, id] : coches.filter((x) => x !== id);
    runs = []; rafraichirListe();
    vueComparaison();
  }));
  R.querySelectorAll(".lien-run").forEach((a) => (a.onclick = (ev) => {
    ev.preventDefault(); nav.ouvrir("resultats", { run: a.closest("label").dataset.id });
  }));
}

// ------------------------------------------------------------------ comparaison
async function vueComparaison() {
  D = null;
  const d = R.querySelector("#r-droite");
  d.className = "comparaison-runs";
  d.innerHTML = `
    <section class="carte trace-comp"><div class="carte-tete"><span class="titre">q en fonction de ε axial</span>
      <span class="overline">${coches.length ? "clic sur un nom : détail du run" : "cochez des runs à gauche"}</span></div>
      <canvas id="r-courbe"></canvas></section>
    <section class="carte synthese"><div class="carte-tete"><span class="titre">Synthèse</span>
      <span class="overline">pic filtré (médiane glissante ±7) · max brut pour mémoire · définitions de depouille.py</span></div>
      <table id="r-synth"></table></section>`;
  courbeComp = new Courbe(d.querySelector("#r-courbe"));
  for (const id of coches) {
    hist[id] ||= await api(`/api/runs/${encodeURIComponent(id)}/historique`);
    syntheses[id] ||= await api(`/api/runs/${encodeURIComponent(id)}/synthese`);
  }
  const nom = (id) => runs.find((r) => r.id === id)?.nom || id.split("~")[1];
  courbeComp.definir(coches.map((id, i) => ({ x: hist[id].eps, y: hist[id].q, couleur: SERIES[i % 6], nom: nom(id) })));
  d.querySelector("#r-synth").innerHTML = `<tr><th>Run</th>${COLONNES.map((c) => `<th>${c[1]}</th>`).join("")}</tr>` +
    coches.map((id, i) => `<tr data-id="${id}"><td><span class="trait" style="display:inline-block;margin-right:8px;background:${SERIES[i % 6]}"></span><a href="#resultats/${encodeURIComponent(id)}" class="lien-synth">${nom(id)}</a></td>` +
      COLONNES.map(([k, , n]) => `<td>${n === 0 && syntheses[id][k] != null ? syntheses[id][k].toLocaleString("fr-FR") : fr(syntheses[id][k], n)}</td>`).join("") + "</tr>").join("");
  d.querySelectorAll(".lien-synth").forEach((a) => (a.onclick = (ev) => { ev.preventDefault(); nav.ouvrir("resultats", { run: a.closest("tr").dataset.id }); }));
}

// ------------------------------------------------------------------ détail d'un run
async function ouvrirDetail(id) {
  const d = R.querySelector("#r-droite");
  d.className = "detail-run";
  d.innerHTML = `
    <div class="detail-tete">
      <button class="bouton" id="d-retour" style="height:30px;padding:0 12px">← Comparaison</button>
      <span class="titre" id="d-nom"></span>
      <div class="segments" id="d-champs"></div>
      <div class="puces" id="d-puces">
        <button class="puce actif" data-mode="1"><i class="pastille trac"></i>Traction</button>
        <button class="puce actif" data-mode="2"><i class="pastille cis"></i>Cisaillement</button>
        <button class="puce actif" data-mode="4"><i class="pastille pre"></i>Pré-rompus</button></div>
    </div>
    <section class="carte vue"><canvas id="d-gl"></canvas>
      <div class="incrust haut-gauche" id="d-info"></div>
      <div class="echelle"><div class="echelle-titre" id="d-et"></div><div class="echelle-corps"><div class="echelle-barre" id="d-eb"></div><div class="echelle-ticks" id="d-tk"></div></div></div>
      <div class="incrust bas-gauche aide">molette : zoom · glisser : déplacer · double-clic : recadrer</div>
      <div class="voile-attente" id="d-attente">Préparation de l'affichage…</div></section>
    <div class="colonne-detail">
      <section class="carte courbe"><div class="carte-tete"><span class="titre">q en fonction de ε axial</span><span class="overline">clic : aller à la frame</span></div>
        <div class="courbe-zone"><canvas id="d-courbe"></canvas></div></section>
      <section class="carte stats"><div class="carte-tete"><span class="titre">Frame courante</span></div><table id="d-stats"></table></section>
    </div>
    <footer class="temps detail-temps">
      <button id="d-lecture" class="bouton-rond" title="Lecture / pause (espace)">▶</button>
      <input id="d-curseur" type="range" min="0" max="0" value="0"><div class="temps-lecture mono" id="d-temps"></div></footer>`;
  d.querySelector("#d-retour").onclick = () => { history.replaceState(null, "", "#resultats"); vueComparaison(); };
  const r = runs.find((x) => x.id === id) || { nom: id.split("~")[1] };
  d.querySelector("#d-nom").textContent = r.nom;
  const vue = new VueChamps(d.querySelector("#d-gl"));
  const courbe = new Courbe(d.querySelector("#d-courbe"));
  D = { id, vue, courbe, lecture: false };

  // cache d'affichage : converti une fois, côté serveur
  let e = await api(`/api/runs/${encodeURIComponent(id)}/cache`, {});
  const t0 = performance.now();
  while (e.etat === "en_cours") {
    d.querySelector("#d-attente").textContent = `Préparation de l'affichage (une seule fois) : ${fr((performance.now() - t0) / 1000, 0)} s`;
    await enAttente(400);
    if (D?.id !== id) return;
    e = await api(`/api/runs/${encodeURIComponent(id)}/cache`);
  }
  if (e.etat !== "pret") { d.querySelector("#d-attente").textContent = e.etat; return; }
  d.querySelector("#d-attente").remove();
  const meta = e.meta;
  hist[id] ||= await api(`/api/runs/${encodeURIComponent(id)}/historique`);
  syntheses[id] ||= await api(`/api/runs/${encodeURIComponent(id)}/synthese`);
  await vue.charger(`/cache/${encodeURIComponent(id)}`, meta);
  D.meta = meta;

  const seg = d.querySelector("#d-champs");
  seg.innerHTML = vue.champsDisponibles().map((c) => `<button data-c="${c}">${CHAMPS[c].nom}</button>`).join("");
  seg.querySelectorAll("button").forEach((b) => (b.onclick = () => { vue.setChamp(b.dataset.c); echelle(); }));
  d.querySelectorAll("#d-puces .puce").forEach((b) => (b.onclick = () => { vue.masque ^= +b.dataset.mode; b.classList.toggle("actif"); vue.marquer(); }));
  const curseur = d.querySelector("#d-curseur");
  curseur.max = meta.nFrames - 1;
  curseur.oninput = () => frame(+curseur.value);
  courbe.surClic = ({ index }) => {
    let best = 0, dm = Infinity;
    meta.frameVersHist.forEach((k, f) => { const x = Math.abs(k - index); if (x < dm) { dm = x; best = f; } });
    frame(best);
  };
  const btn = d.querySelector("#d-lecture");
  btn.onclick = () => basculer();
  echelle();
  frame(meta.nFrames - 1);
}

function basculer() {
  if (!D) return;
  D.lecture = !D.lecture;
  R.querySelector("#d-lecture").textContent = D.lecture ? "❚❚" : "▶";
  const pas = () => {
    if (!D?.lecture) return;
    frame((D.vue.frame + 1) % D.meta.nFrames);
    setTimeout(pas, 140);
  };
  pas();
}

function frame(f) {
  const { vue, courbe, meta, id } = D;
  vue.setFrame(f);
  R.querySelector("#d-curseur").value = f;
  const k = meta.frameVersHist[f], h = hist[id], s = syntheses[id], c = vue.comptes[f];
  courbe.definir([{ x: h.eps, y: h.q, couleur: css("--accent"), nom: "", jusqua: k }], { serie: 0, index: k });
  R.querySelector("#d-temps").textContent = `frame ${String(f).padStart(2, "0")} / ${meta.nFrames - 1}   ·   t = ${fr(meta.temps[f] * 1e3, 3)} ms`;
  R.querySelector("#d-info").innerHTML = `<div class="grand">${fr(h.q[k], 1)} MPa</div><div class="sous">q à ε = ${fr(h.eps[k], 3)} % · σ₃ = ${fr(meta.sigma3_MPa, 0)} MPa</div>`;
  const lignes = [["Temps", `${fr(meta.temps[f] * 1e3, 3)} ms`], ["ε axial", `${fr(h.eps[k], 3)} %`], ["q", `${fr(h.q[k], 1)} MPa`],
    ["q pic (filtré)", `${fr(s.q_pic_filtre_MPa, 1)} MPa à ${fr(s.eps_pic_filtre_pct, 3)} %`],
    ["Rompus en traction", c[1].toLocaleString("fr-FR")], ["Rompus en cisaillement", c[2].toLocaleString("fr-FR")],
    ["Pré-rompus", c[4].toLocaleString("fr-FR")], ["Éléments · joints", `${meta.nTri.toLocaleString("fr-FR")} · ${meta.nJoint.toLocaleString("fr-FR")}`]];
  R.querySelector("#d-stats").innerHTML = lignes.map(([a, b]) => `<tr><td>${a}</td><td>${b}</td></tr>`).join("");
}

function echelle() {
  const { vue } = D, ch = CHAMPS[vue.champ];
  R.querySelectorAll("#d-champs button").forEach((b) => b.classList.toggle("actif", b.dataset.c === vue.champ));
  const barre = R.querySelector("#d-eb"), ticks = R.querySelector("#d-tk");
  if (ch.cat) {
    R.querySelector("#d-et").textContent = "Phase";
    const c = vue.couleursPhases;
    barre.style.background = `linear-gradient(${c[0]} 0 33%, ${c[1]} 33% 66%, ${c[2]} 66%)`;
    ticks.innerHTML = ["phase 0", "phase 1", "phase 2"].map((n, i) => `<span style="top:${(i + .5) * 33.3}%">${n}</span>`).join("");
    return;
  }
  const [a, b] = vue.bornesChamp();
  R.querySelector("#d-et").textContent = `${ch.nom} (${ch.unite})`;
  const arrets = [];
  for (let i = 0; i < BANDES; i++) {
    const col = `rgb(${couleurArc((BANDES - 1 - i + .5) / BANDES).map((v) => Math.round(255 * v))})`;
    arrets.push(`${col} ${(i / BANDES * 100).toFixed(2)}%`, `${col} ${((i + 1) / BANDES * 100).toFixed(2)}%`);
  }
  barre.style.background = `linear-gradient(${arrets.join(",")})`;
  ticks.innerHTML = Array.from({ length: BANDES + 1 }, (_, i) =>
    `<span style="top:${i / BANDES * 100}%">${fr((b - (b - a) * i / BANDES) * ch.k, ch.unite === "%" ? 3 : 0)}</span>`).join("");
}

addEventListener("keydown", (ev) => {
  if (!D?.meta || !R.classList.contains("actif") || ev.target.tagName === "INPUT" && ev.target.type !== "range") return;
  if (ev.key === "ArrowRight") frame(Math.min(D.vue.frame + 1, D.meta.nFrames - 1));
  else if (ev.key === "ArrowLeft") frame(Math.max(D.vue.frame - 1, 0));
  else if (ev.key === " ") { ev.preventDefault(); basculer(); }
});
