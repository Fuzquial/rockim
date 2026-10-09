// Loi cohésive des joints de g1, côté navigateur : copie de studio2/noyau/loi.py (elle-même
// confrontée à tools/yan_point.cpp à 1e-9 près). Module pur, sans DOM : tests/test_loi_js.py
// l'exécute sous Node et le compare à la version Python.

export function fD(D, a = 0.63, b = 1.8, c = 6.0) {
  if (D <= 0) return 1;
  if (D >= 1) return 0;
  const s = a + b;
  if (Math.abs(s) < 1e-12 || Math.abs(1 - s) < 1e-12) return a * (1 - D) + b * (1 - D) ** c;
  const K = (s - 1) / s, alpha = (a + c * b) / (s * (1 - s));
  const f = (1 - K * Math.exp(alpha * D)) * (a * (1 - D) + b * (1 - D) ** c);
  return Math.min(1, Math.max(0, f));
}

export function integrale(a = 0.63, b = 1.8, c = 6.0, n = 4096) {
  if (n % 2) n += 1;
  const h = 1 / n;
  let s = fD(0, a, b, c) + fD(1, a, b, c);
  for (let i = 1; i < n; i++) s += fD(i * h, a, b, c) * (i % 2 ? 4 : 2);
  return s * h / 3;
}

export class Loi {
  // p : { pj, ft, coh, phi, GfI, GfII, a, b, c, adoucissement: "yan"|"linear", montee: "linear"|"parabolic" }
  constructor(p) {
    Object.assign(this, { a: 0.63, b: 1.8, c: 6.0, adoucissement: "yan", montee: "linear" }, p);
    this.tanPhi = Math.tan(this.phi * Math.PI / 180);
    this.yan = this.adoucissement === "yan";
    this.I = this.yan ? integrale(this.a, this.b, this.c) : 0.5;
    this.dnE = this.ft / this.pj;
    this.ot = this.GfI / (this.ft * this.I);
    this.st = this.GfII / (this.coh * this.I);
  }
  f(D) { return this.yan ? fD(D, this.a, this.b, this.c) : Math.min(1, Math.max(0, 1 - D)); }
  sigma(dn) {
    if (dn <= this.dnE) {
      const r = dn / this.dnE;
      return this.montee === "parabolic" ? this.ft * (2 * r - r * r) : this.pj * dn;
    }
    return Math.min(this.pj * dn, this.f(Math.min(1, (dn - this.dnE) / this.ot)) * this.ft);
  }
  // Mode II, retour radial pas à pas : seule la cohésion s'adoucit, le frottement reste.
  mode2(dtg, sigmaN = 0) {
    let slip = 0, D = 0;
    return dtg.map((x) => {
      const tr = this.pj * (x - slip);
      const lim = this.f(D) * this.coh + this.tanPhi * Math.max(0, -sigmaN);
      const t = Math.max(-lim, Math.min(lim, tr));
      if (t !== tr) { slip += (tr - t) / this.pj; D = Math.max(D, Math.min(1, Math.abs(slip) / this.st)); }
      return t;
    });
  }
}
