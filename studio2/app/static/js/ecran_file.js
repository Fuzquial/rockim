// Écran File : état en direct, actions, réglages d'exécution ; courbe et journal du run choisi.
import { api, css, fr } from "./commun.js";
import { Courbe } from "./courbe.js";

let R, nav, choisi = null, courbe, derniereCle = "";
const LIB = { en_attente: ["attente", "en attente"], en_cours: ["cours", "en cours"], fini: ["fini", "terminé"],
  deja_fait: ["fini", "déjà fait"], echec: ["echec", "échec"], arrete: ["attente", "arrêté"], interrompu: ["echec", "interrompu"] };

const duree = (s) => s == null ? "" : s < 90 ? `${Math.round(s)} s` : s < 5400 ? `${Math.round(s / 60)} min`
  : `${Math.floor(s / 3600)} h ${String(Math.round(s % 3600 / 60)).padStart(2, "0")}`;

export async function monter(noeud, ctx) {
  R = noeud; nav = ctx;
  R.innerHTML = `
    <section class="carte" style="overflow:auto;min-height:0">
      <div class="barre-file">
        <span class="titre">File de calculs</span>
        <span class="mono2">exe</span><select id="f-exe" style="min-width:180px"></select>
        <span class="mono2">jobs</span><div class="nombre" style="width:64px"><input id="f-jobs"></div>
        <span class="mono2">fils</span><div class="nombre" style="width:64px"><input id="f-fils"></div>
        <span style="margin-left:auto" class="mono2">les calculs continuent si l'interface est fermée</span>
        <button class="bouton" id="f-pause" style="padding:0 14px"></button>
      </div>
      <table class="tab-file" id="f-tab"></table>
      <div class="vide-file" id="f-vide" hidden>La file est vide. Les essais s'ajoutent depuis l'écran Essai.</div>
    </section>
    <div class="colonne-droite">
      <section class="carte"><div class="carte-tete"><span class="titre" id="f-nom">Aucun run choisi</span><span class="overline" id="f-etat"></span></div>
        <div style="position:relative;height:230px"><canvas id="f-courbe" style="position:absolute;inset:4px 8px;width:calc(100% - 16px);height:calc(100% - 8px)"></canvas></div>
        <div class="chiffres" id="f-chiffres"></div></section>
      <section class="carte"><div class="carte-tete"><span class="titre">Journal</span><span class="overline">sortie du solveur</span></div>
        <div class="journal" id="f-journal"></div></section>
    </div>`;
  courbe = new Courbe(R.querySelector("#f-courbe"));
  const etat = await api("/api/etat");
  R.querySelector("#f-exe").innerHTML = etat.exes.map((e) => `<option ${e === etat.reglages.exe ? "selected" : ""}>${e}</option>`).join("");
  R.querySelector("#f-jobs").value = etat.reglages.jobs;
  R.querySelector("#f-fils").value = etat.reglages.fils;
  const regler = (cle, v) => api("/api/reglages", { [cle]: v }).then(rafraichir);
  R.querySelector("#f-exe").onchange = (ev) => regler("exe", ev.target.value);
  R.querySelector("#f-jobs").onchange = (ev) => regler("jobs", Math.max(1, parseInt(ev.target.value) || 1));
  R.querySelector("#f-fils").onchange = (ev) => regler("fils", Math.max(1, parseInt(ev.target.value) || 1));
  R.querySelector("#f-pause").onclick = async () => {
    const f = await api("/api/file");
    regler("pause", !f.pause);
  };
  await rafraichir();
  setInterval(() => { if (R.classList.contains("actif")) rafraichir(); }, 1500);
}

export function activer() { rafraichir(); }

async function rafraichir() {
  const f = await api("/api/file");
  R.querySelector("#f-pause").textContent = f.pause ? "Reprendre les lancements" : "Suspendre les lancements";
  R.querySelector("#f-pause").classList.toggle("primaire", f.pause);
  R.querySelector("#f-vide").hidden = f.travaux.length > 0;
  const maintenant = Date.now() / 1000;
  // regroupement par lot, dans l'ordre de la file
  let html = `<tr><th>Run</th><th>État</th><th>Avancement</th><th>Écoulé</th><th>Reste</th><th></th></tr>`, lot = undefined;
  for (const t of f.travaux) {
    if (t.lot !== lot) { lot = t.lot; html += `<tr class="groupe-lot"><td colspan="6">${lot || "Hors lot"}</td></tr>`; }
    const [cls, lib] = LIB[t.etat] || ["attente", t.etat];
    const ecoule = t.debut ? (t.fin || maintenant) - t.debut : null;
    const reste = t.etat === "en_cours" && t.avancement > 0.02 ? ecoule * (1 - t.avancement) / t.avancement : null;
    const actions = t.etat === "en_cours" ? [["arreter", "■", "arrêter (le dossier reste exploitable)"]]
      : t.etat === "en_attente" ? [["monter", "▲", "plus tôt"], ["descendre", "▼", "plus tard"], ["retirer", "✕", "retirer de la file"]]
      : [["relancer", "↻", "relancer (l'ancien dossier est renommé, pas supprimé)"], ["retirer", "✕", "retirer de la file (les fichiers restent)"]];
    html += `<tr data-id="${t.id}" class="${t.id === choisi ? "sel" : ""}">
      <td><div class="nom">${t.nom}</div><div class="desc">${t.resume}</div></td>
      <td><span class="badge ${cls}">${lib}${t.code && t.etat === "echec" ? " · code " + t.code : ""}</span></td>
      <td><div class="barre-av"><i style="width:${(100 * t.avancement).toFixed(1)}%"></i></div></td>
      <td class="mono2">${duree(ecoule)}</td><td class="mono2">${reste ? "≈ " + duree(reste) : ""}</td>
      <td><div class="actions-ligne">${actions.map(([a, s, titre]) => `<button class="ico" data-a="${a}" title="${titre}">${s}</button>`).join("")}
        ${t.etat === "fini" || t.etat === "deja_fait" ? `<button class="ico" data-a="ouvrir" title="ouvrir les résultats">→</button>` : ""}</div></td></tr>`;
  }
  const cle = html.replace(/style="width:[^"]*"|<td class="mono2">[^<]*<\/td>/g, "");
  if (cle !== derniereCle) {                     // ne reconstruit le tableau que si sa structure change
    R.querySelector("#f-tab").innerHTML = html;
    derniereCle = cle;
    R.querySelectorAll("#f-tab tr[data-id]").forEach((tr) => (tr.onclick = (ev) => {
      const a = ev.target.closest("button")?.dataset.a, id = tr.dataset.id;
      if (a === "ouvrir") return nav.ouvrir("resultats", { run: "espace~" + tr.querySelector(".nom").textContent });
      if (a) return api(`/api/file/${id}/${a}`, {}).then(rafraichir).catch((e) => alert(e.message));
      choisi = id; derniereCle = ""; rafraichir();
    }));
  } else {                                        // sinon, seulement les barres et les durées
    f.travaux.forEach((t, i) => {
      const tr = R.querySelector(`#f-tab tr[data-id="${t.id}"]`);
      if (!tr) return;
      tr.querySelector(".barre-av i").style.width = `${(100 * t.avancement).toFixed(1)}%`;
      const ecoule = t.debut ? (t.fin || maintenant) - t.debut : null;
      const reste = t.etat === "en_cours" && t.avancement > 0.02 ? ecoule * (1 - t.avancement) / t.avancement : null;
      const tds = tr.querySelectorAll("td.mono2");
      tds[0].textContent = duree(ecoule); tds[1].textContent = reste ? "≈ " + duree(reste) : "";
    });
  }
  if (!choisi && f.travaux.length) choisi = (f.travaux.find((t) => t.etat === "en_cours") || f.travaux[0]).id;
  const t = f.travaux.find((x) => x.id === choisi);
  if (t) await detail(t);
}

async function detail(t) {
  R.querySelector("#f-nom").textContent = t.nom;
  R.querySelector("#f-etat").textContent = `${(LIB[t.etat] || [, t.etat])[1]} · ${fr(100 * t.avancement, 0)} %`;
  try {
    const h = await api(`/api/runs/${encodeURIComponent("espace~" + t.nom)}/historique`);
    courbe.definir(h.eps.length ? [{ x: h.eps, y: h.q, couleur: css("--accent"), nom: t.nom }] : []);
    const n = h.q.length - 1;
    const pic = n >= 0 ? Math.max(...h.q.filter((v) => v != null)) : null;
    R.querySelector("#f-chiffres").innerHTML = `
      <div><small>t simulé</small><b>${n >= 0 ? fr(h.t[n], 2) + " ms" : "—"}</b></div>
      <div><small>q actuel</small><b>${n >= 0 ? fr(h.q[n], 1) + " MPa" : "—"}</b></div>
      <div><small>q max (brut)</small><b>${pic != null ? fr(pic, 1) + " MPa" : "—"}</b></div>
      <div><small>Lignes d'historique</small><b>${(n + 1).toLocaleString("fr-FR")}</b></div>`;
  } catch (e) { courbe.definir([]); }
  const j = await api(`/api/file/${t.id}/journal?n=80`);
  const zone = R.querySelector("#f-journal"), bas = zone.scrollTop + zone.clientHeight >= zone.scrollHeight - 4;
  zone.textContent = j.lignes.join("\n") || "(pas encore de sortie)";
  if (bas) zone.scrollTop = zone.scrollHeight;
}
