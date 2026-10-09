// Outils partagés par les écrans : API, format français, jetons de thème, toile 2D.

export const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

export function fr(v, n = 1) {
  if (v === null || v === undefined || Number.isNaN(v)) return "—";
  const x = Math.abs(v) < 0.5 * 10 ** -n ? 0 : Number(v);      // pas de « -0 »
  return x.toLocaleString("fr-FR", { minimumFractionDigits: n, maximumFractionDigits: n });
}

export async function api(chemin, corps) {
  const r = await fetch(chemin, corps === undefined ? {} : {
    method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(corps),
  });
  const d = r.headers.get("Content-Type")?.startsWith("application/json") ? await r.json() : await r.arrayBuffer();
  if (!r.ok) throw new Error(d.erreur || `HTTP ${r.status}`);
  return d;
}

export async function binaire(url, Type = Float32Array) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`${url} : HTTP ${r.status}`);
  return new Type(await r.arrayBuffer());
}

export function el(html) {
  const t = document.createElement("template");
  t.innerHTML = html.trim();
  return t.content.firstElementChild;
}

// Canvas à la taille de son conteneur, en pixels physiques ; rend le contexte et la taille CSS.
export function toile(cv) {
  const dpr = devicePixelRatio || 1, w = cv.clientWidth, h = cv.clientHeight;
  if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) {
    cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
  }
  const c = cv.getContext("2d");
  c.setTransform(dpr, 0, 0, dpr, 0, 0);
  c.clearRect(0, 0, w, h);
  return [c, w, h];
}

// Palette des séries de données (rockim-design), dans l'ordre.
export const SERIES = ["#2FBECB", "#E5A54B", "#A78BFA", "#6BCB77", "#E77C8D", "#5A9CF8"];

export const enAttente = (ms) => new Promise((r) => setTimeout(r, ms));
