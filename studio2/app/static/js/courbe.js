// Traceur de courbes sur canvas 2D : plusieurs séries, marqueur, clic vers le point le plus proche.
// Fait maison parce que les courbes q-epsilon ne sont pas monotones (retour après le pic).
import { css, fr, toile } from "./commun.js";

function pas(a, b, n = 5) {
  const p = (b - a) / n, e = 10 ** Math.floor(Math.log10(p || 1));
  return [1, 2, 2.5, 5, 10].map((m) => m * e).find((s) => s >= p);
}

export class Courbe {
  constructor(canvas, { titreX = "ε axial (%)", titreY = "q (MPa)" } = {}) {
    this.cv = canvas; this.titreX = titreX; this.titreY = titreY;
    this.series = []; this.marqueur = null; this.surClic = null; this.M = { g: 52, d: 16, h: 14, b: 38 };
    canvas.addEventListener("click", (ev) => this._clic(ev));
    new ResizeObserver(() => this.dessiner()).observe(canvas);
  }

  // series : [{ x, y, couleur, nom, epais, jusqua }] ; marqueur : { serie, index }
  definir(series, marqueur = null) { this.series = series; this.marqueur = marqueur; this.dessiner(); }

  dessiner() {
    const [c, w, h] = toile(this.cv), M = this.M;
    if (!this.series.length || !w) return;
    let xa = 0, xb = -Infinity, ya = 0, yb = -Infinity;
    for (const s of this.series) for (let i = 0; i < s.x.length; i++) {
      if (s.x[i] == null || s.y[i] == null) continue;
      xa = Math.min(xa, s.x[i]); xb = Math.max(xb, s.x[i]); ya = Math.min(ya, s.y[i]); yb = Math.max(yb, s.y[i]);
    }
    if (!(xb > xa)) xb = xa + 1;
    yb = yb > ya ? yb * 1.06 : ya + 1;
    const X = (v) => M.g + (v - xa) / (xb - xa) * (w - M.g - M.d), Y = (v) => h - M.b - (v - ya) / (yb - ya) * (h - M.h - M.b);
    this.X = X; this.Y = Y;
    c.font = `10.5px ${css("--mono")}`; c.strokeStyle = css("--border"); c.fillStyle = css("--texte3"); c.lineWidth = 1;
    const px = pas(xa, xb), py = pas(ya, yb), dx = Math.max(0, -Math.floor(Math.log10(px))), dy = Math.max(0, -Math.floor(Math.log10(py)));
    c.textAlign = "center";
    for (let v = Math.ceil(xa / px) * px; v <= xb + 1e-12; v += px) {
      c.beginPath(); c.moveTo(X(v) + .5, M.h); c.lineTo(X(v) + .5, h - M.b); c.stroke(); c.fillText(fr(v, dx), X(v), h - M.b + 15);
    }
    c.textAlign = "right";
    for (let v = Math.ceil(ya / py) * py; v <= yb + 1e-12; v += py) {
      c.beginPath(); c.moveTo(M.g, Y(v) + .5); c.lineTo(w - M.d, Y(v) + .5); c.stroke(); c.fillText(fr(v, dy), M.g - 6, Y(v) + 3.5);
    }
    c.fillStyle = css("--texte2"); c.font = `11px ${css("--police")}`; c.textAlign = "center";
    c.fillText(this.titreX, M.g + (w - M.g - M.d) / 2, h - 6);
    c.save(); c.translate(13, M.h + (h - M.h - M.b) / 2); c.rotate(-Math.PI / 2); c.fillText(this.titreY, 0, 0); c.restore();
    c.lineJoin = "round";
    for (const s of this.series) {
      const fin = s.jusqua ?? s.x.length - 1;
      const trace = (de, a, col, ep) => {
        c.strokeStyle = col; c.lineWidth = ep; c.beginPath(); let lev = true;
        for (let i = de; i <= a; i++) {
          if (s.x[i] == null || s.y[i] == null) { lev = true; continue; }
          if (lev) c.moveTo(X(s.x[i]), Y(s.y[i])); else c.lineTo(X(s.x[i]), Y(s.y[i]));
          lev = false;
        }
        c.stroke();
      };
      if (fin < s.x.length - 1) trace(fin, s.x.length - 1, css("--border-fort"), 1.4);
      trace(0, fin, s.couleur, s.epais || 1.8);
    }
    if (this.series.length > 1) {
      c.font = `11px ${css("--police")}`; c.textAlign = "left";
      this.series.forEach((s, i) => {
        c.fillStyle = s.couleur; c.fillRect(w - M.d - 190, M.h + 6 + 16 * i, 12, 3);
        c.fillStyle = css("--texte2"); c.fillText(s.nom, w - M.d - 172, M.h + 10 + 16 * i);
      });
    }
    if (this.marqueur) {
      const s = this.series[this.marqueur.serie], i = this.marqueur.index;
      if (s && s.x[i] != null) {
        c.fillStyle = css("--bg2"); c.strokeStyle = css("--texte"); c.lineWidth = 2;
        c.beginPath(); c.arc(X(s.x[i]), Y(s.y[i]), 5, 0, 2 * Math.PI); c.fill(); c.stroke();
      }
    }
  }

  _clic(ev) {
    if (!this.surClic || !this.X) return;
    const r = this.cv.getBoundingClientRect(), px = ev.clientX - r.left, py = ev.clientY - r.top;
    let best = null, dmin = Infinity;
    this.series.forEach((s, k) => {
      for (let i = 0; i < s.x.length; i++) {
        if (s.x[i] == null || s.y[i] == null) continue;
        const d = (this.X(s.x[i]) - px) ** 2 + (this.Y(s.y[i]) - py) ** 2;
        if (d < dmin) { dmin = d; best = { serie: k, index: i }; }
      }
    });
    if (best) this.surClic(best);
  }
}
