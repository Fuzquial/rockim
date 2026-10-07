# -*- coding: utf-8 -*-
"""check_s2.py — verification FALSIFIANTE de contactForcePairs (S2, campagne du 13/09).

usage : python tests_f2/campagne13/S2/check_s2.py <out_fc> <out_ref> <log_fc> [<log_ref>]

Criteres (chacun imprime OK / ECHEC) :
  C1  Fc_insert_rock_* == 0 sur TOUTES les lignes (a 30 us l onde n a pas atteint la roche)
      et Fc_rock_insert_* == 0 aussi.
  C2  Fc_bit_piston_z != 0 apres le contact piston/bit ; instant de premiere valeur non nulle,
      signe attendu POSITIF (le bit repousse le piston vers +z), Fc_piston_bit_z negatif.
  C3  antisymetrie EXACTE : Fc_piston_bit + Fc_bit_piston == 0 (meme sommes, signes opposes).
  C4  grpFz (trackGroup = bit, force de contact nette sur le bit) == somme des Fc_<X>_bit_z
      a 1e-9 relatif (les deux accumulateurs lisent les memes forces, ordre different).
  C5  quantite de mouvement du piston : integrale trapezes de Fc_bit_piston_z + m_p g_z
      contre m_p (vz_p(T) - vz_p(0)) — comparaison a l estimateur de tools/fig_fp.py
      (m_p lue dans le journal). Ecart imprime, seuil 5 % (echantillonnage des lignes).
  C6  bit-identite : history.csv du temoin == history.csv du deck fc PRIVE des colonnes Fc_*
      (texte identique), VTU byte-identiques, journaux identiques hors la ligne d annonce
      et la ligne de temps mur.
"""
import csv
import hashlib
import os
import re
import sys


def load(run):
    with open(os.path.join(run, "history.csv"), newline="") as f:
        rows = [r for r in csv.DictReader(f) if r]
    return rows


def mass(log, body):
    m = re.search(r"corps '%s':.*?masse = ([\d.eE+-]+) kg" % body, open(log, errors="replace").read())
    return float(m.group(1)) if m else None


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def main():
    out_fc, out_ref, log_fc = sys.argv[1:4]
    log_ref = sys.argv[4] if len(sys.argv) > 4 else None
    R = load(out_fc)
    cols = list(R[0].keys())
    fc_cols = [c for c in cols if c.startswith("Fc_")]
    print("lignes history :", len(R), "; colonnes Fc :", fc_cols)
    ok_all = True

    def verdict(name, ok, msg):
        nonlocal ok_all
        ok_all = ok_all and ok
        print("%-4s %s : %s" % (name, "OK   " if ok else "ECHEC", msg))

    # C1
    mx = 0.0
    for r in R:
        for c in ("Fc_insert_rock_x", "Fc_insert_rock_y", "Fc_insert_rock_z",
                  "Fc_rock_insert_x", "Fc_rock_insert_y", "Fc_rock_insert_z"):
            mx = max(mx, abs(float(r[c])))
    verdict("C1", mx == 0.0, "max |Fc_insert_rock|, |Fc_rock_insert| = %g N (attendu 0 exact)" % mx)

    # C2
    t = [float(r["t"]) for r in R]
    fbp = [float(r["Fc_bit_piston_z"]) for r in R]
    fpb = [float(r["Fc_piston_bit_z"]) for r in R]
    first = next((i for i, v in enumerate(fbp) if v != 0.0), None)
    if first is None:
        verdict("C2", False, "Fc_bit_piston_z reste nul sur %d lignes" % len(R))
    else:
        imax = max(range(len(fbp)), key=lambda i: abs(fbp[i]))
        pos = sum(1 for v in fbp if v > 0)
        neg = sum(1 for v in fbp if v < 0)
        verdict("C2", fbp[first] > 0 and fpb[first] < 0,
                "premiere force non nulle a t = %.4g us (ligne %d) : Fc_bit_piston_z = %+.4g N, "
                "Fc_piston_bit_z = %+.4g N ; max |Fc_bit_piston_z| = %.4g N a t = %.4g us ; "
                "lignes > 0 : %d, < 0 : %d" % (t[first] * 1e6, first, fbp[first], fpb[first],
                                               abs(fbp[imax]), t[imax] * 1e6, pos, neg))

    # C3
    mx = 0.0
    for r in R:
        for a in "xyz":
            mx = max(mx, abs(float(r["Fc_piston_bit_" + a]) + float(r["Fc_bit_piston_" + a])))
    verdict("C3", mx == 0.0, "max |Fc_piston_bit + Fc_bit_piston| = %g N (attendu 0 exact)" % mx)

    # C4 — history.csv est ecrit a 6 chiffres significatifs (ostream par defaut) :
    # la somme de termes ARRONDIS (1e5 N -> resolution 0,1 N) ne peut egaler le total
    # arrondi qu a la demi-unite du 6e chiffre de chaque terme pres. Tolerance =
    # somme de ces demi-unites (terme par terme + total) ; en deca, les deux
    # accumulateurs lisent les memes forces. Un terme manquant (une paire oubliee,
    # un point d accumulation absent) donnerait un ecart de l ordre du terme lui-meme.
    import math

    def ulp6(x):
        x = abs(x)
        return 0.5 * 10.0 ** (math.floor(math.log10(x)) - 5) if x > 0 else 0.0

    partners = [c for c in fc_cols if c.endswith("_bit_z")]
    mxabs = 0.0
    mxratio = 0.0
    worst = ""
    for r in R:
        g = float(r["grpFz"])
        terms = [float(r[c]) for c in partners]
        s = sum(terms)
        tol = ulp6(g) + sum(ulp6(v) for v in terms)
        d = abs(g - s)
        mxabs = max(mxabs, d)
        ratio = d / tol if tol > 0 else (0.0 if d == 0 else float("inf"))
        if ratio > mxratio:
            mxratio = ratio
            worst = "t = %.4g us : grpFz = %.10g, somme = %.10g, ecart %.3g N, tolerance d arrondi %.3g N" % (
                float(r["t"]) * 1e6, g, s, d, tol)
    verdict("C4", mxratio <= 1.0, "colonnes %s ; ecart absolu max %.3g N ; ecart / tolerance d arrondi 6 chiffres max %.3g (seuil 1) ; pire : %s"
            % (partners, mxabs, mxratio, worst))

    # C5
    mp = mass(log_fc, "piston")
    vp = [float(r["vz_piston"]) for r in R]
    g = 9.81
    imp = 0.0
    for i in range(1, len(R)):
        imp += 0.5 * (fbp[i] + fbp[i - 1]) * (t[i] - t[i - 1])
    imp_g = -mp * g * (t[-1] - t[0])
    dP = mp * (vp[-1] - vp[0])
    rel = abs(imp + imp_g - dP) / max(abs(dP), 1e-300)
    verdict("C5", rel <= 0.05,
            "m_p = %.6g kg ; impulsion trapezes de Fc_bit_piston_z = %.6g N.s, pesanteur %.3g N.s ; "
            "m_p dv_p = %.6g N.s ; ecart %.2f %% (seuil 5 %% ; %d lignes, pas moyen %.3g us)"
            % (mp, imp, imp_g, dP, 100 * rel, len(R), 1e6 * (t[-1] - t[0]) / max(1, len(R) - 1)))

    # C6
    with open(os.path.join(out_fc, "history.csv"), newline="") as f:
        hdr = f.readline().rstrip("\n").split(",")
        keep = [i for i, c in enumerate(hdr) if not c.startswith("Fc_")]
        fc_txt = [",".join(hdr[i] for i in keep)]
        for line in f:
            p = line.rstrip("\n").split(",")
            if len(p) < 2:
                continue
            fc_txt.append(",".join(p[i] for i in keep))
    with open(os.path.join(out_ref, "history.csv"), newline="") as f:
        ref_txt = [l.rstrip("\n") for l in f if l.strip()]
    same_hist = fc_txt == ref_txt
    ndiff = sum(1 for a, b in zip(fc_txt, ref_txt) if a != b) + abs(len(fc_txt) - len(ref_txt))
    verdict("C6a", same_hist, "history.csv temoin == fc prive des Fc_* : %d lignes, %d differentes"
            % (len(ref_txt), ndiff))
    vt = sorted(x for x in os.listdir(out_fc) if x.endswith(".vtu"))
    vr = sorted(x for x in os.listdir(out_ref) if x.endswith(".vtu"))
    same_vtu = vt == vr and all(sha(os.path.join(out_fc, x)) == sha(os.path.join(out_ref, x)) for x in vt)
    verdict("C6b", same_vtu, "%d VTU, sha256 identiques : %s" % (len(vt), same_vtu))
    if log_ref:
        def strip(p):
            out = []
            keys_line = ""
            for l in open(p, errors="replace"):
                if "contactForcePairs" in l or "wall time" in l or "Fc_" in l:
                    continue
                if l.startswith("[rockim] cles :"):      # comptage des cles du deck :
                    keys_line = l.strip()[:60]           # differe d une cle par construction
                    continue
                if "potential stats" in l:               # compteurs gardes, temps (tGrid, tLoop) masques
                    l = l.split(", tGrid")[0] + "\n"
                out.append(l.rstrip("\n"))
            return out, keys_line
        (a, ka), (b, kb) = strip(log_fc), strip(log_ref)
        nd = sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))
        verdict("C6c", nd == 0, "journaux hors annonce, temps mur et comptage des cles : %d lignes differentes "
                "(comptage fc « %s », ref « %s »)" % (nd, ka, kb))
        if nd:
            for x, y in zip(a, b):
                if x != y:
                    print("   fc  :", x[:160])
                    print("   ref :", y[:160])
                    break

    print("BILAN S2 :", "TOUS OK" if ok_all else "AU MOINS UN ECHEC")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
