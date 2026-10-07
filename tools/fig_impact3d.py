#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_impact3d.py — depouillement de l'impact 3D d'insert unique : GIF de coupe
# et planche de courbes, comparees a Yang et al. 2026 (IJRMMS 206, 106660).
#
#   python tools/fig_impact3d.py out_3d_fin [out_3d_gros]
#
# Le premier dossier est le run PRINCIPAL (trace epais, coupe de droite), le
# second la REFERENCE optionnelle (trait fin, coupe de gauche).
#
# COUPE EXACTE, pas un nuage de centroides. Chaque tetraedre est intersecte
# analytiquement avec le plan y = D/2 : le polygone rendu est la section reelle.
# La version « tricontourf sur les centroides » fabriquait une fausse surface
# bombee, parce que le centroide d'un element de 15 mm est 7 mm sous la peau
# alors que celui d'un element de 2 mm n'est qu'a 1 mm (bug du 2026-09-03).
#
# CADRAGE SERRE (28 x 15 mm) : le cratere fait moins d'1 mm de profond et les
# ruptures tiennent dans r < 12 mm. Le cadrage large d'origine (64 x 54 mm)
# rendait le sujet illisible — l'utilisateur l'a signale, a juste titre.
#
# La position de l'outil est le CENTRE de la sphere, lu dans history.csv :
# z_centre = toolZ, donc R - delta au-dessus de la surface. La poser a -delta
# enfonce le cercle de R = 8,51 mm dans la roche (autre bug du meme jour).
#
# CIBLES PUBLIEES tracees d'office (leur §5 et fig. 9b, cas piston 9 m/s) :
#   enfoncement 1,07 mm | sortie 4,65 m/s | e = 0,83 | 23,84 -> 16,31 J
# ---------------------------------------------------------------------------
import csv
import io
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from PIL import Image

FIELD = "vonMises"    # --field : vonMises | sigma1 | pMean | tauMax
CLIM = (0.0, 140.0)   # --clim MIN MAX [MPa]
CMAP = "YlGnBu_r"     # --cmap ; sigma1 et pMean sont SIGNES -> divergent
R = 0.00851           # rayon de l'insert [m]
H0 = 0.100            # hauteur du bloc [m]
XC = YC = 0.075       # axe d'impact [m]
V0 = 5.62             # vitesse d'entree [m/s]
KE0 = 23.83           # energie injectee [J]
KE_OUT = 16.31        # leur energie de sortie [J]
V_OUT = 4.65          # leur vitesse de sortie [m/s]
D_PUB = 1.07          # leur enfoncement [mm]
TEAL, GREY, OCHRE = "#12808C", "#5a6570", "#B8752A"
ERR, TOOL, PUR = "#8E1B10", "#F09A3E", "#7B4EA8"
XL, XR, YB, YT = -14.0, 14.0, -11.0, 4.5
E6 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def rd(p):
    return io.open(p, encoding="utf-8", errors="ignore").read()


def pts(s):
    m = re.search(r"<Points>.*?<DataArray[^>]*>(.*?)</DataArray>", s, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)


def cn(s, n):
    m = re.search(r'Name="connectivity"[^>]*>(.*?)</DataArray>', s, re.S)
    return np.fromstring(m.group(1), sep=" ").astype(int).reshape(-1, n)


def ar(s, nm):
    m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % nm, s, re.S)
    return np.fromstring(m.group(1), sep=" ") if m else None


def crushed(s):
    """Champ d'ENDOMMAGEMENT DE BROYAGE de l'element, quelle que soit la loi.

    L'outil lisait `bulkD` en dur. Or `bulkDamage` est NEUTRALISE des qu'une
    loi de volume existe (Fdem3dSolver.cpp:2736) : sous la loi de la note de
    septembre 2026 (`law = saksala` + `compDamage = crackband`) le champ
    n'est pas ecrit et `bd[ei]` levait un TypeError sur None.

    Le porteur equivalent y est `omegaC`, l'endommagement de COMPRESSION de
    l'eq. 6 de la note — exactement la meme grandeur physique (le broyage),
    ecrit au VTU par le lot du 2026-09-11. On prend donc bulkD s'il existe
    (comportement inchange, bit-identique sur tous les runs anterieurs),
    omegaC sinon, et un vecteur nul si aucun des deux n'est present.
    """
    bd = ar(s, "bulkD")
    if bd is not None:
        return bd
    wc = ar(s, "omegaC")
    if wc is not None:
        return wc
    n = ar(s, "vonMises")
    return np.zeros(0 if n is None else len(n))


def slice_tets(P, C, y0):
    """Intersection EXACTE des tetraedres avec le plan y = y0."""
    s = P[C][:, :, 1] - y0
    po, kp = [], []
    for i in np.where(~((s > 0).all(1) | (s < 0).all(1)))[0]:
        v = P[C[i]]
        d = v[:, 1] - y0
        q = []
        for a, b in E6:
            if d[a] * d[b] < 0:
                w = d[a] / (d[a] - d[b])
                q.append(v[a] + w * (v[b] - v[a]))
            elif d[a] == 0.0:
                q.append(v[a])
        if len(q) < 3:
            continue
        q = np.unique(np.round(np.array(q), 12), axis=0)
        if len(q) < 3:
            continue
        xz = np.c_[(q[:, 0] - XC) * 1e3, (q[:, 2] - H0) * 1e3]
        ang = np.arctan2(xz[:, 1] - xz[:, 1].mean(), xz[:, 0] - xz[:, 0].mean())
        po.append(xz[np.argsort(ang)])
        kp.append(i)
    return po, np.array(kp, dtype=int)


def slice_tris(P, C, y0):
    """Trace des facettes de joint dans le plan : segments, pas des points."""
    s = P[C][:, :, 1] - y0
    sg = []
    for i in np.where(~((s > 0).all(1) | (s < 0).all(1)))[0]:
        v = P[C[i]]
        d = v[:, 1] - y0
        q = []
        for a, b in [(0, 1), (1, 2), (2, 0)]:
            if d[a] * d[b] < 0:
                w = d[a] / (d[a] - d[b])
                q.append(v[a] + w * (v[b] - v[a]))
            elif d[a] == 0.0:
                q.append(v[a])
        if len(q) >= 2:
            q = np.array(q[:2])
            sg.append(np.c_[(q[:, 0] - XC) * 1e3, (q[:, 2] - H0) * 1e3])
    return sg


def hist(D):
    h = list(csv.DictReader(io.open(D + "/history.csv", encoding="utf-8")))
    g = lambda k: np.array([float(a[k]) for a in h])
    return dict(t=g("t") * 1e6, F=np.abs(g("toolFz")) / 1e3,
                d=(H0 + R + 1e-4 - g("toolZ")) * 1e3, tz=(g("toolZ") - H0) * 1e3,
                vz=g("toolVz"), nb=g("nBroken"), nf=g("nFrag"), ke=g("toolKE"))


def frames(D):
    return {int(l.split(",")[0]): float(l.split(",")[1])
            for l in rd(D + "/frames.csv").strip().split("\n")[1:]}


def state(h):
    """Renvoie (fini, e, t_separation) — e n'a de sens qu'apres separation."""
    k = int(np.argmax(h["d"]))
    c = np.where(h["F"] > 0.05)[0]
    done = (h["t"][-1] - h["t"][c[-1]]) > 5.0
    return done, h["vz"][k:].max() / V0, h["t"][c[-1]]


# --------------------------------------------------------------------------
def curves(runs, out):
    fig, ax = plt.subplots(2, 2, figsize=(13.2, 8.6), dpi=125)
    (A, B), (C, D) = ax
    for i, (Dn, lab, col, ls) in enumerate(runs):
        h = hist(Dn)
        lw = 2.6 if i == 0 else 1.7
        k = int(np.argmax(h["d"]))
        A.plot(h["d"], h["F"], color=col, ls=ls, lw=lw, label=lab)
        A.plot(h["d"][k], h["F"][k], "o", ms=6, color=col)
        B.plot(h["t"], h["vz"], color=col, ls=ls, lw=lw, label=lab)
        C.plot(h["t"], h["F"], color=col, ls=ls, lw=lw, label=lab)
        D.plot(h["t"], h["ke"], color=col, ls=ls, lw=lw, label=lab)
    h = hist(runs[0][0])
    done, e, tsep = state(h)
    tag = "SEPARE a %.0f µs" % tsep if done else "contact ENCORE actif"
    A.plot(h["d"][-1], h["F"][-1], "s", ms=9, color=TEAL, zorder=6)
    if ON_YANG_BENCH:                      # cibles du banc de Yang seulement
        A.axvline(D_PUB, color=ERR, lw=1.3, ls="-.")
        A.annotate("publié %.2f mm" % D_PUB, (D_PUB, 5), textcoords="offset points",
                   xytext=(7, 0), fontsize=9, color=ERR)
    A.set_xlabel("pénétration δ (mm)", fontsize=10.5)
    A.set_ylabel("force F (kN)", fontsize=10.5)
    A.set_title("Force–pénétration", fontsize=11.5)
    if ON_YANG_BENCH:
        A.set_xlim(0, 1.25)          # echelle du banc ; ailleurs, auto

    B.axhline(0, color="k", lw=.8, ls=":")
    B.axhline(V_OUT, color=ERR, lw=1.4, ls="-.")
    B.annotate("sortie publiée %.2f m/s  →  e = 0,83" % V_OUT, (14, V_OUT + .2),
               fontsize=9.5, color=ERR)
    B.plot(h["t"][-1], h["vz"][-1], "s", ms=9, color=TEAL, zorder=6)
    B.annotate("%.2f m/s\n%s" % (h["vz"][-1], tag), (h["t"][-1], h["vz"][-1]),
               textcoords="offset points", xytext=(-124, 22), fontsize=9,
               color=TEAL, arrowprops=dict(arrowstyle="->", lw=1.2, color=TEAL))
    B.set_xlabel("t (µs)", fontsize=10.5)
    B.set_ylabel("vitesse de l'outil (m/s)", fontsize=10.5)
    B.set_title("Vitesse : c'est ici que se lit e", fontsize=11.5)
    B.set_ylim(-6, 5.6)

    C.axvline(300, color=GREY, lw=1, ls=":")
    C.annotate("fin des 2 premiers runs", (300, 3), textcoords="offset points",
               xytext=(-115, 0), fontsize=8.5, color=GREY)
    C.axvspan(373, 406, color=ERR, alpha=.10)
    C.annotate("séparation\nextrapolée", (390, 50), ha="center", fontsize=8.5, color=ERR)
    C.set_xlabel("t (µs)", fontsize=10.5)
    C.set_ylabel("force F (kN)", fontsize=10.5)
    C.set_title("La séparation, c'est F = 0", fontsize=11.5)

    D.axhline(KE0, color="k", lw=1, ls=":")
    D.annotate("énergie injectée %.2f J" % KE0, (12, KE0 + .6), fontsize=9)
    D.axhline(KE_OUT, color=ERR, lw=1.4, ls="-.")
    D.annotate("KE de sortie publiée %.2f J" % KE_OUT, (12, KE_OUT + .6),
               fontsize=9.5, color=ERR)
    D.set_xlabel("t (µs)", fontsize=10.5)
    D.set_ylabel("énergie cinétique de l'outil (J)", fontsize=10.5)
    D.set_title("Le budget : %.2f J entrent, combien ressortent ?" % KE0, fontsize=11.5)
    D.set_ylim(0, 26)

    for a in (A, B, C, D):
        a.grid(alpha=.22, lw=.5)
        a.legend(fontsize=9, frameon=False)
        for s in ["top", "right"]:
            a.spines[s].set_visible(False)
    for a in (B, C, D):
        a.set_xlim(0, 450)
    pct = min(100.0, h["t"][-1] / 4.5)
    fig.suptitle("Impact 3D Kuru — 1,509 kg à 5,62 m/s — %s"
                 % ("run TERMINÉ, e = %.3f" % e if done else "run à %.0f %%" % pct),
                 fontsize=12.5, y=.985)
    fig.tight_layout(rect=[0, 0, 1, .955])
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return done, e


def gif(runs, out, tmp="_gif3d"):
    os.makedirs(tmp, exist_ok=True)
    HH = {d: hist(d) for d, _, _, _ in runs}
    FR = {d: frames(d) for d, _, _, _ in runs}
    main = runs[0][0]
    ks = [k for k in sorted(FR[main]) if k <= max(FR[main]) - 1]
    imgs = []
    for k in ks:
        t = FR[main][k] * 1e6
        fig = plt.figure(figsize=(12.8, 6.6), dpi=115)
        gs = fig.add_gridspec(2, len(runs), height_ratios=[.55, 1.35],
                              hspace=.32, wspace=.09)
        a0 = fig.add_subplot(gs[0, :])
        for i, (D, lab, col, ls) in enumerate(runs):
            h = HH[D]
            m = h["t"] <= t
            if i == 0:
                a0.plot(h["t"], h["vz"], color="#9aa5ad", lw=.7, alpha=.30)
            a0.plot(h["t"][m], h["vz"][m], color=col, ls=ls,
                    lw=2.2 if i == 0 else 1.4, label=lab)
        a0.axhline(0, color="k", lw=.8, ls=":")
        a0.axhline(V_OUT, color=ERR, lw=1.3, ls="-.")
        a0.annotate("publié %.2f m/s (e = 0,83)" % V_OUT, (9, V_OUT + .25),
                    fontsize=9, color=ERR)
        a0.axvline(t, color="k", lw=1, alpha=.45)
        a0.set_xlim(0, 450)
        a0.set_ylim(-6.2, 5.8)
        a0.set_xlabel("t (µs)", fontsize=9)
        a0.set_ylabel("vitesse (m/s)", fontsize=9)
        a0.grid(alpha=.18, lw=.5)
        a0.legend(fontsize=9, frameon=False, loc="lower right")
        for s_ in ["top", "right"]:
            a0.spines[s_].set_visible(False)
        a0.set_title("Vitesse de l'outil — le rebond", fontsize=10)
        for j, (D, lab, col, ls) in enumerate(runs):
            axx = fig.add_subplot(gs[1, j])
            kk = k if k in FR[D] else max(FR[D])
            s = rd("%s/fdem3d_%04d.vtu" % (D, kk))
            js = rd("%s/fdem3d_joints_%04d.vtu" % (D, kk))
            P, C = pts(s), cn(s, 4)
            # --field (14/09) : quel champ colorie la coupe. Defaut
            # vonMises, echelle 0-140 MPa = comportement d origine.
            # sigma1 reproduit la fig. 16 de Yang (contrainte principale
            # maximale, traction > 0) ; pMean la pression moyenne ;
            # tauMax le cisaillement maximal.
            vm, bd = ar(s, FIELD) / 1e6, crushed(s)
            po, ei = slice_tets(P, C, YC)
            pc = PolyCollection(po, array=np.clip(vm[ei], CLIM[0], CLIM[1]),
                                cmap=CMAP, edgecolors="#5c686f", linewidths=.32)
            pc.set_clim(*CLIM)
            axx.add_collection(pc)
            pv = bd[ei] > 0.5
            if pv.any():
                axx.add_collection(PolyCollection(
                    [po[q] for q in np.where(pv)[0]], facecolors="none",
                    edgecolors=PUR, linewidths=1.9))
            jp, jc, dm = pts(js), cn(js, 3), ar(js, "damage")
            sg = slice_tris(jp, jc[dm >= 0.999], YC)
            if sg:
                axx.add_collection(LineCollection(sg, colors=ERR, linewidths=2.4))
            h = HH[D]
            tt = min(t, h["t"][-1])
            axx.add_patch(plt.Circle((0, np.interp(tt, h["t"], h["tz"])), R * 1e3,
                                     facecolor=TOOL, alpha=.45, edgecolor=TOOL,
                                     lw=2, zorder=6))
            axx.axhline(0, color="k", lw=.7, ls=":", alpha=.55)
            axx.set_xlim(XL, XR)
            axx.set_ylim(YB, YT)
            axx.set_aspect("equal")
            axx.set_xlabel("x − axe (mm)", fontsize=8.5)
            axx.tick_params(labelsize=7.5)
            if j == 0:
                axx.set_ylabel("z − surface (mm)", fontsize=8.5)
            else:
                axx.set_yticklabels([])
            sfx = "" if k in FR[D] else "  — FIN DU RUN"
            axx.set_title("%s · %d rompus · %d frag%s"
                          % (lab, np.interp(tt, h["t"], h["nb"]),
                             np.interp(tt, h["t"], h["nf"]), sfx),
                          fontsize=10.5, color=col)
        hf = HH[main]
        fig.suptitle("t = %5.1f µs  ·  delta %5.3f mm  ·  F %5.1f kN  ·  vz %+5.2f m/s"
                     % (t, R * 1e3 + 0.1 - np.interp(t, hf["t"], hf["tz"]),
                        np.interp(t, hf["t"], hf["F"]), np.interp(t, hf["t"], hf["vz"])),
                     fontsize=11, y=.985, family="monospace")
        p = "%s/f%03d.png" % (tmp, k)
        fig.savefig(p, bbox_inches="tight")
        plt.close(fig)
        imgs.append(p)
    ims = [Image.open(p).convert("P", palette=Image.ADAPTIVE) for p in imgs]
    ims[0].save(out, save_all=True, append_images=ims[1:], duration=460, loop=0)
    return len(ims), os.path.getsize(out) // 1024


def geometry_override(argv):
    """Geometrie du banc, surchargeable depuis la ligne de commande.

    Les constantes en tete de ce fichier decrivent le banc de Yang (bloc de
    100 mm, insert R = 8,51 mm, axe en 75 mm, 5,62 m/s). Elles etaient CODEES
    EN DUR : le script rendait donc une coupe faussement cadree — et des
    annotations fausses — sur tout autre deck. Ces drapeaux sont OPTIONNELS et
    leurs defauts sont les valeurs d'origine : sans eux, comportement inchange.

      --field vonMises|sigma1|pMean|tauMax   champ colorant la coupe (14/09)
      --clim MIN MAX [MPa]           --cmap  nom matplotlib
      --R  rayon de l'outil [m]      --H   hauteur du bloc [m]
      --XC axe d'impact x [m]        --YC  axe d'impact y [m]
      --V0 vitesse d'entree [m/s]    --KE0 energie injectee [J]
      --win XL XR YB YT  fenetre de la coupe [mm]
      --label "..."      legende du run principal
    Les cibles publiees (D_PUB, V_OUT, KE_OUT, e = 0,83) ne sont tracees que si
    --V0 n'est pas pose : hors du banc de Yang elles n'ont aucun sens.
    """
    global R, H0, XC, YC, V0, KE0, XL, XR, YB, YT, ON_YANG_BENCH
    global FIELD, CLIM, CMAP
    global D_PUB, V_OUT, KE_OUT
    ON_YANG_BENCH = "--V0" not in argv
    if not ON_YANG_BENCH:
        # Hors du banc de Yang, ses cibles publiees n'ont aucun sens. On les
        # met a NaN plutot que de garder un axvline/axhline a une valeur
        # etrangere au deck : matplotlib ne trace alors rien, et aucune des
        # annotations correspondantes n'atterrit sur la figure. C'est le
        # minimum qui rende l'outil honnete sur un autre bloc.
        D_PUB = V_OUT = KE_OUT = float("nan")
    rest, i = [], 0
    while i < len(argv):
        a = argv[i]
        if a == "--win" and i + 4 < len(argv):
            XL, XR, YB, YT = [float(x) for x in argv[i + 1:i + 5]]
            i += 5
        elif a in ("--R", "--H", "--XC", "--YC", "--V0", "--KE0") and i + 1 < len(argv):
            v = float(argv[i + 1])
            if a == "--R":   R = v
            elif a == "--H": H0 = v
            elif a == "--XC": XC = v
            elif a == "--YC": YC = v
            elif a == "--V0": V0 = v
            else:            KE0 = v
            i += 2
        elif a == "--field" and i + 1 < len(argv):
            FIELD = argv[i + 1]
            if FIELD not in ("vonMises", "sigma1", "pMean", "tauMax"):
                raise SystemExit("--field : vonMises | sigma1 | pMean | tauMax")
            # defauts adaptes au champ, surchargeables par --clim / --cmap
            if FIELD == "sigma1":   CLIM, CMAP = (-20.0, 20.0), "RdBu_r"
            elif FIELD == "pMean":  CLIM, CMAP = (-200.0, 20.0), "RdBu_r"
            elif FIELD == "tauMax": CLIM, CMAP = (0.0, 80.0), "YlGnBu_r"
            i += 2
        elif a == "--clim" and i + 2 < len(argv):
            CLIM = (float(argv[i + 1]), float(argv[i + 2]))
            i += 3
        elif a == "--cmap" and i + 1 < len(argv):
            CMAP = argv[i + 1]
            i += 2
        elif a == "--label" and i + 1 < len(argv):
            rest.append(("label", argv[i + 1]))
            i += 2
        else:
            rest.append(a)
            i += 1
    return rest


ON_YANG_BENCH = True


def main():
    argv = geometry_override(sys.argv[1:])
    lab = "1,3 mm · µ 0,18"
    for e in list(argv):
        if isinstance(e, tuple) and e[0] == "label":
            lab = e[1]
            argv.remove(e)
    if not argv:
        raise SystemExit("usage: fig_impact3d.py <run> [reference] "
                         "[--R m --H m --XC m --YC m --V0 m/s --KE0 J "
                         "--win XL XR YB YT --label txt]")
    runs = [(argv[0], lab, TEAL, "-")]
    if len(argv) > 1:
        runs.append((argv[1], "2,0 mm · µ 0,60", GREY, "--"))
    D = argv[0]
    done, e = curves(runs, D + "/fp_raffine.png")
    n, ko = gif(runs, D + "/raffine_ab.gif")
    h = hist(D)
    print("courbes -> %s/fp_raffine.png" % D)
    print("GIF %d ko, %d trames -> %s/raffine_ab.gif" % (ko, n, D))
    print("t = %.1f µs | delta max %.3f mm | pic %.1f kN | %d rompus | KE %.2f J"
          % (h["t"][-1], h["d"].max(), h["F"].max(), h["nb"][-1], h["ke"][-1]))
    print("e = %.3f  %s   (cible 0,83)" % (e, "MESURE" if done else "provisoire"))


if __name__ == "__main__":
    main()
