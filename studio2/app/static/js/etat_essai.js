// État partagé de l'essai en préparation (écrans Essai, Maillage, Loi des joints).
// Le brouillon est gardé dans le navigateur : fermer l'onglet ne perd pas la saisie.
import { api } from "./commun.js";

const CLE = "rockim.essai.brouillon";
let choix = null, defaut = null, formulaire = null, resultat = null, abonnes = [], minuterie = null, version = 0;

export async function initialiser() {
  if (formulaire) return;
  formulaire = await api("/api/formulaire");
  defaut = formulaire.defaut;
  try { choix = JSON.parse(localStorage.getItem(CLE)) || null; } catch (e) { choix = null; }
  // Lien de partage : #<écran>?brouillon=<choix JSON encodé> remplace le brouillon local.
  const partage = new URLSearchParams(location.hash.split("?")[1] || "").get("brouillon");
  if (partage) {
    try { choix = JSON.parse(partage); localStorage.setItem(CLE, partage); } catch (e) { /* lien abîmé : on l'ignore */ }
  }
  choix = { ...defaut, ...(choix || {}) };
  await reconstruire();
}

export const lireChoix = () => choix;
export const lireFormulaire = () => formulaire;
export const lireResultat = () => resultat;
export function abonner(f) { abonnes.push(f); if (resultat) f(resultat, choix); }

export function modifier(changements, immediat = false) {
  choix = { ...choix, ...changements };
  try { localStorage.setItem(CLE, JSON.stringify(choix)); } catch (e) { /* stockage indisponible */ }
  clearTimeout(minuterie);
  minuterie = setTimeout(reconstruire, immediat ? 0 : 120);
}

export function reinitialiser() {
  choix = { ...defaut };
  modifier({}, true);
}

async function reconstruire() {
  const v = ++version;
  const r = await api("/api/essai/construire", choix);
  if (v !== version) return;                       // une saisie plus récente est déjà partie
  resultat = r;
  abonnes.forEach((f) => f(resultat, choix));
}

// Lien qui recrée l'essai courant chez un collègue (même serveur ou un autre).
export function lienDePartage() {
  const modifie = Object.fromEntries(Object.entries(choix).filter(([k, v]) => JSON.stringify(v) !== JSON.stringify(defaut[k])));
  return `${location.origin}${location.pathname}#essai?brouillon=${encodeURIComponent(JSON.stringify(modifie))}`;
}
