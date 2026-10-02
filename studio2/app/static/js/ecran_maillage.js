// Écran Maillage : aperçu du maillage RÉEL (run g1 tronqué à une frame), confronté à l'estimation,
// et dessin de fissures pré-rompues à la souris.
import { api, css, enAttente, fr } from "./commun.js";
import { abonner, initialiser, lireChoix, lireResultat, modifier } from "./etat_essai.js";
import { VueChamps } from "./vue_champs.js";

let R, nav, vue = null, apercu = null, empreinte = null, dessin = false, trace = null;

// Ce qui change le maillage ou les discontinuités : le reste (nom, durée, frames...) ne périme pas l'aperçu.
const SANS_EFFET = ["nom", "description", "T", "frames", "vitesse", "chute_arret", "sigma3_MPa", "type_essai", "deformations_historique", "champs_deformation", "delai_axial"];
const empreinteDe = (c) => JSON.stringify(Object.fromEntries(Object.entries(c).filter(([k]) => !SANS_EFFET.includes(k)).sort()));

export async function monter(noeud, ctx) {
  R = noeud; nav = ctx;
  await initialiser();
  R.innerHTML = `
  <section class="carte vue-maillage">
    <canvas id="m-gl"></canvas>
    <svg id="m-calque" class="calque"></svg>
    <div class="incrust haut-gauche"><div style="font-size:15px;font-weight:600">Aperçu du maillage réel</div>
      <div class="sous" id="m-sous">généré par g1, run tronqué à une frame</div></div>
    <div class="outils">
      <button class="bouton" id="m-dessin">✎ Dessiner une fissure</button>
      <button class="bouton" id="m-effacer">Effacer les fissures dessinées</button>
    </div>
    <div class="legende">
      <span><i class="pastille-phase" style="background:var(--quartz)"></i>quartz</span>
      <span><i class="pastille-phase" style="background:var(--feldspath)"></i>feldspath</span>
      <span><i class="pastille-phase" style="background:var(--biotite)"></i>biotite</span>
      <span><i class="pastille" style="background:var(--pre);opacity:1;width:14px"></i> pré-rompus</span>
      <span><i class="pastille" style="background:var(--attention);opacity:1;width:14px"></i> dessinées</span>
    </div>
    <div class="voile-attente" id="m-voile">
      <div style="text-align:center;max-width:420px;line-height:1.6">
        <div style="font-size:15px;color:var(--texte);margin-bottom:6px">Pas encore d'aperçu pour cet essai</div>
        Le maillage d'un GBM est construit par le solveur lui-même. L'aperçu lance g1 sur une seule frame
        (quelques secondes, un fil) et montre exactement le maillage du calcul.<br><br>
        <button class="bouton primaire" id="m-generer" style="padding:0 18px">Générer l'aperçu</button>
        <div class="message" id="m-message" style="margin-top:10px"></div>
      </div></div>
  </section>
  <div class="colonne-droite">
    <section class="carte" id="m-perime" hidden><div class="corps-carte" style="color:var(--attention)">
      L'essai a changé depuis cet aperçu. <button class="lien" id="m-regenerer">Régénérer</button></div></section>
    <section class="carte"><div class="carte-tete"><span class="titre">Ce que le solveur a construit</span></div>
      <table class="tab-phases" id="m-construit" style="margin:8px 0 4px"></table></section>
    <section class="carte"><div class="carte-tete"><span class="titre">Phases</span><span class="overline">aire réalisée / visée</span></div>
      <table class="tab-phases" id="m-phases" style="margin:8px 0 4px"></table></section>
    <section class="carte"><div class="carte-tete"><span class="titre">Discontinuités posées</span></div>
      <ul class="avis" id="m-disc"></ul>
      <div class="boutons"><button class="bouton" id="m-essai">← Modifier l'essai</button>
        <button class="bouton primaire" id="m-ajouter">Maillage validé : ajouter à la file</button>
        <div class="message" id="m-msg-ajout"></div></div></section>
  </div>`;
  vue = new VueChamps(R.querySelector("#m-gl"));
  vue.champ = "phase";
  vue.apresDessin = calque;
  R.querySelector("#m-generer").onclick = generer;
  R.querySelector("#m-regenerer").onclick = generer;
  R.querySelector("#m-essai").onclick = () => nav.ouvrir("essai");
  R.querySelector("#m-ajouter").onclick = ajouter;
  R.querySelector("#m-dessin").onclick = () => basculerDessin();
  R.querySelector("#m-effacer").onclick = () => modifier({ segments: [] }, true);
  abonner((r, c) => { tableaux(); R.querySelector("#m-perime").hidden = !apercu || empreinteDe(c) === empreinte; calque(); });
}

export async function activer() {
  tableaux();
  // Un aperçu déjà calculé pour cet essai s'affiche tout seul, sans relancer g1.
  if (apercu && empreinteDe(lireChoix()) === empreinte) return;
  try {
    const e = await api("/api/court/apercu/existant", lireChoix());
    if (e.etat === "fini") await generer();
  } catch (err) { /* pas d'aperçu : le voile reste */ }
}

async function generer() {
  const msg = R.querySelector("#m-message");
  R.querySelector("#m-voile").hidden = false;
  msg.textContent = "Lancement de g1…"; msg.className = "message";
  const choix = lireChoix();
  try {
    let e = await api("/api/court/apercu", choix);
    const t0 = performance.now();
    while (e.etat === "en_cours" || e.etat === "absent" || (e.etat === "fini" && e.cache?.etat !== "pret")) {
      msg.textContent = `${e.etat === "fini" ? "Préparation de l'affichage" : "g1 construit le maillage"} : ${fr((performance.now() - t0) / 1000, 0)} s`;
      await enAttente(400);
      e = await api(`/api/court/${e.cle}`);
    }
    if (e.etat !== "fini") throw new Error(`${e.etat}\n${(e.journal || []).join("\n")}`);
    apercu = e;
    empreinte = empreinteDe(choix);
    await vue.charger(`/cache/${encodeURIComponent(e.id)}`, e.cache.meta);
    vue.setChamp("phase");
    vue.setFrame(0);
    R.querySelector("#m-voile").hidden = true;
    R.querySelector("#m-perime").hidden = true;
    R.querySelector("#m-sous").textContent = `généré par g1 en une frame · clé ${e.cle}`;
    tableaux();
  } catch (err) {
    msg.textContent = "Échec : " + err.message; msg.className = "message err";
  }
}

function tableaux() {
  const r = lireResultat();
  if (!r) return;
  const est = r.estimation, d = apercu?.diagnostics || {};
  const val = (v, f = (x) => x.toLocaleString("fr-FR")) => (v == null ? "—" : f(v));
  const lignes = [["Éléments", val(est.elements), val(d.elements)], ["Grains", val(est.grains), val(d.grains)],
    ["Éléments par grain", fr(est.elements_par_grain, 1), d.grains ? fr(d.elements / d.grains, 1) : "—"],
    ["Joints", "—", val(d.joints)], ["dont joints de grain", "—", val(d.joints_de_grain)],
    ["Hauteur min. d'élément", "—", d.hmin_m ? `${fr(d.hmin_m * 1e3, 2)} mm` : "—"]];
  R.querySelector("#m-construit").innerHTML = `<tr><th></th><th>Estimé</th><th>Réalisé</th></tr>` +
    lignes.map((l) => `<tr><td>${l[0]}</td><td>${l[1]}</td><td>${l[2]}</td></tr>`).join("");
  const noms = { quartz: "Quartz", feldspar: "Feldspath", biotite: "Biotite" }, couleurs = [css("--quartz"), css("--feldspath"), css("--biotite")];
  const ph = d.phases_realisees || [];
  R.querySelector("#m-phases").innerHTML = ph.length ? `<tr><th>Phase</th><th>Visée</th><th>Réalisée</th><th>Écart</th></tr>` +
    ph.map((p, i) => `<tr><td><i class="pastille-phase" style="background:${couleurs[i]}"></i>${noms[p.nom] || p.nom}</td><td>${fr(p.cible_pct, 2)} %</td><td>${fr(p.realise_pct, 2)} %</td><td>${(p.realise_pct >= p.cible_pct ? "+" : "") + fr(p.realise_pct - p.cible_pct, 2)}</td></tr>`).join("")
    : `<tr><td style="color:var(--texte3)">${apercu ? "une seule phase" : "après l'aperçu"}</td></tr>`;
  const c = lireChoix(), av = [];
  if (d.prerompus != null) av.push(`<li class="info">${d.prerompus.toLocaleString("fr-FR")} joints pré-rompus sur ${d.joints.toLocaleString("fr-FR")} (${fr(d.prerompus_pct, 1)} %).</li>`);
  if (d.prerompus_libres != null) av.push(`<li class="att">${d.prerompus_libres.toLocaleString("fr-FR")} seulement ont une extrémité déjà scindée : fraction glissante effective ${fr(100 * d.prerompus_libres / d.joints, 1)} %.</li>`);
  if (d.joints_affaiblis != null) av.push(`<li class="info">${d.joints_affaiblis.toLocaleString("fr-FR")} joints affaiblis sur les plans.</li>`);
  if ((c.segments || []).length) av.push(`<li class="info">${c.segments.length} fissure${c.segments.length > 1 ? "s" : ""} dessinée${c.segments.length > 1 ? "s" : ""}${apercu && empreinteDe(c) !== empreinte ? " (régénérer pour voir les joints retenus)" : ""}.</li>`);
  if (!av.length) av.push(`<li class="ok">${apercu ? "Aucune discontinuité préexistante." : "Générer l'aperçu pour voir les discontinuités posées par le solveur."}</li>`);
  R.querySelector("#m-disc").innerHTML = av.join("");
  const err = r.avis.some((a) => a.niveau === "erreur");
  R.querySelector("#m-ajouter").disabled = err;
  R.querySelector("#m-generer").disabled = err || c.maillage !== "voronoi";
  if (c.maillage !== "voronoi") R.querySelector("#m-message").textContent = "En Gmsh homogène, le maillage vient du fichier .msh : pas d'aperçu par le solveur.";
}

// ------------------------------------------------------------------ fissures dessinées
function basculerDessin() {
  dessin = !dessin;
  R.querySelector("#m-dessin").classList.toggle("primaire", dessin);
  R.querySelector("#m-gl").style.cursor = dessin ? "crosshair" : "";
  vue.outil = dessin ? {
    down: (x, y) => { trace = [x, y, x, y]; calque(); },
    move: (x, y) => { if (trace) { trace[2] = x; trace[3] = y; calque(); } },
    up: (x, y) => {
      if (!trace) return;
      const [x1, y1] = trace, [W, H] = [vue.meta.bbox[2], vue.meta.bbox[3]];
      const borne = (v, m) => Math.min(m, Math.max(0, v));
      const s = [borne(x1, W), borne(y1, H), borne(x, W), borne(y, H)].map((v) => +v.toFixed(6));
      trace = null;
      if (Math.hypot(s[2] - s[0], s[3] - s[1]) > 5e-4) modifier({ segments: [...(lireChoix().segments || []), s] }, true);
      calque();
    },
  } : null;
}

function calque() {
  const svg = R.querySelector("#m-calque");
  if (!vue?.meta) { svg.innerHTML = ""; return; }
  const segs = [...(lireChoix().segments || []), ...(trace ? [trace] : [])];
  svg.innerHTML = segs.map(([x1, y1, x2, y2]) => {
    const [a, b] = vue.versEcran(x1, y1), [c, d] = vue.versEcran(x2, y2);
    return `<line x1="${a}" y1="${b}" x2="${c}" y2="${d}" stroke="${css("--attention")}" stroke-width="2.5" stroke-linecap="round"/>`;
  }).join("");
}

async function ajouter() {
  const r = lireResultat(), m = R.querySelector("#m-msg-ajout");
  try {
    await api("/api/essai/ajouter", { essai: r.essai });
    m.innerHTML = `Ajouté à la file. <a href="#file">Voir la file →</a>`; m.className = "message ok";
  } catch (e) { m.textContent = e.message; m.className = "message err"; }
}
