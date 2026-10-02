// Coquille : onglets, chargement paresseux des écrans, état global (file) en bandeau.
import { api } from "./commun.js";

const ECRANS = {
  essai: () => import("./ecran_essai.js"),
  maillage: () => import("./ecran_maillage.js"),
  loi: () => import("./ecran_loi.js"),
  file: () => import("./ecran_file.js"),
  resultats: () => import("./ecran_resultats.js"),
};
const montes = {};
const conteneur = document.getElementById("ecrans");

export async function ouvrir(nom, params = {}) {
  if (!ECRANS[nom]) nom = "resultats";
  document.querySelectorAll("#onglets button").forEach((b) => b.classList.toggle("actif", b.dataset.e === nom));
  for (const [n, m] of Object.entries(montes)) m.noeud.classList.toggle("actif", n === nom);
  if (!montes[nom]) {
    const noeud = document.createElement("main");
    noeud.className = "ecran actif";
    noeud.id = "e-" + nom;
    conteneur.appendChild(noeud);
    try {
      const mod = await ECRANS[nom]();
      montes[nom] = { noeud, mod };
      await mod.monter(noeud, { ouvrir });
    } catch (e) {
      noeud.innerHTML = `<section class="carte" style="padding:20px">Écran « ${nom} » indisponible : ${e.message}</section>`;
      montes[nom] = { noeud, mod: {} };
    }
  }
  montes[nom].mod.activer?.(params);
  const h = "#" + nom + (params.run ? "/" + encodeURIComponent(params.run) : "");
  if (location.hash.split("?")[0] !== h) history.replaceState(null, "", h);
}

document.querySelectorAll("#onglets button").forEach((b) => (b.onclick = () => ouvrir(b.dataset.e)));

// Bandeau : exécutable, jobs, et nombre de runs en cours / en attente.
async function etat() {
  try {
    const f = await api("/api/file");
    const enCours = f.travaux.filter((t) => t.etat === "en_cours").length;
    const attente = f.travaux.filter((t) => t.etat === "en_attente").length;
    document.getElementById("resumeFile").textContent = enCours || attente ? `${enCours} en cours · ${attente}` : "";
    document.getElementById("texteEtat").textContent = `${f.exe} · ${f.jobs} jobs × ${f.fils} fils`;
    document.getElementById("pointEtat").style.background = "";
  } catch (e) {
    document.getElementById("texteEtat").textContent = "serveur injoignable";
    document.getElementById("pointEtat").style.background = "var(--erreur)";
  }
}
etat();
setInterval(etat, 2000);

function suivreAncre() {
  const [nom, run] = location.hash.slice(1).split("?")[0].split("/");
  ouvrir(nom || "resultats", run ? { run: decodeURIComponent(run) } : {});
}
addEventListener("hashchange", suivreAncre);
suivreAncre();
