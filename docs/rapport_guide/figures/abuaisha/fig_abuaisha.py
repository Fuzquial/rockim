#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_abuaisha.py — figures de synthese du banc AbuAisha et al. 2017 (B4)
# pour le rapport-guide rockim.
#
#   python fig_abuaisha.py [HIST] [RUNS]
#
# HIST : historiques archives du banc (binaire corrige du 20/08, rockim_e1 a
#        rockim_e3, Windows/MSVC) — defaut :
#        /home/user/phd_geothermie/FDEM/rockim/bench_abuaisha/historiques
# RUNS : sorties rejouees le 2026-10-06 avec le binaire f0209ef (Linux,
#        1 fil) — defaut : scratchpad/runs/abuaisha
#
# Aucune donnee n'est inventee : toutes les courbes sont relues dans les
# history.csv / fdem_*.vtu. Les seules valeurs ecrites en dur sont les
# references publiees ou analytiques (citees en commentaire).
#
# Produit :
#   fig_abuaisha_pression.pdf  — p(t) contre les seuils ; pic decompose
#   fig_abuaisha_rejeu.pdf     — Lame (paroi) et Parker (ouverture), rejoues
# ---------------------------------------------------------------------------
import csv
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
HIST = sys.argv[1] if len(sys.argv) > 1 else \
    "/home/user/phd_geothermie/FDEM/rockim/bench_abuaisha/historiques"
RUNS = sys.argv[2] if len(sys.argv) > 2 else \
    "/tmp/claude-0/-home-user/95ba17a5-6321-5243-8157-c0b3e0616268/scratchpad/runs/abuaisha"

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10,
    "legend.fontsize": 8.5,
})
# --- Style commun du rapport : Computer Modern (CMU Serif, a defaut Latin Modern
# Roman 10), virgule decimale sur toutes les graduations, signe moins ASCII.
import glob as _glob
import matplotlib.ticker as _mticker
from matplotlib import font_manager as _fm
for _f in (_glob.glob("/usr/share/fonts/**/cmun*.[ot]tf", recursive=True)
           + _glob.glob("/usr/share/texmf/fonts/opentype/public/lm/lmroman10-*.otf")
           + _glob.glob("/usr/share/texlive/texmf-dist/fonts/opentype/public/lm/lmroman10-*.otf")):
    _fm.fontManager.addfont(_f)
_fm.fontManager.ttflist = [_e for _e in _fm.fontManager.ttflist
                           if _e.name != "Latin Modern Roman" or "lmroman10" in _e.fname]
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "Latin Modern Roman", "Computer Modern Roman", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
    "axes.formatter.use_locale": False,
    "pdf.fonttype": 3,
})
_sf_call = _mticker.ScalarFormatter.__call__


def _sf_virgule(self, x, pos=None):
    s = _sf_call(self, x, pos)
    return s.replace(".", "{,}") if "$" in s else s.replace(".", ",")


_mticker.ScalarFormatter.__call__ = _sf_virgule
# --- fin du style commun


def virg(x, nd=1):
    """nombre a la francaise"""
    return ("%.*f" % (nd, x)).replace(".", ",")


def fmt_axes(ax, nx=1, ny=1):
    from matplotlib.ticker import FuncFormatter
    ax.xaxis.set_major_formatter(FuncFormatter(
        lambda v, _: ("%g" % round(v, 6)).replace(".", ",")))
    ax.yaxis.set_major_formatter(FuncFormatter(
        lambda v, _: ("%g" % round(v, 6)).replace(".", ",")))


def read_hist(path):
    with open(path) as f:
        r = list(csv.DictReader(f))
    return {k: np.array([float(x[k]) for x in r]) for k in r[0]}


def read_vtu_points(path):
    with open(path) as f:
        txt = f.read()
    P = np.array(re.search(r"<Points>.*?<DataArray[^>]*>(.*?)</DataArray>",
                           txt, re.S).group(1).split(), float).reshape(-1, 3)[:, :2]
    C = np.array(re.search(r'Name="connectivity"[^>]*>(.*?)</DataArray>',
                           txt, re.S).group(1).split(), int).reshape(-1, 3)
    return P, C


def complete(path):
    with open(path, "rb") as f:
        f.seek(max(0, os.path.getsize(path) - 200))
        return b"</VTKFile>" in f.read()


def frames(run):
    return [f for f in sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
            if complete(f)]


# ------------------------------------------------------------------ donnees
# leur eq. 10 (Hubbert-Willis / Fjaer) : p_HF = 3 s'_h - s'_H + f_t
FT = 5.0
CIBLE = {"aniso": 3 * 4.6 - 6.8 + FT, "iso46": 3 * 4.6 - 4.6 + FT,
         "iso35": 3 * 3.5 - 3.5 + FT}

RUNS_ARCH = [  # (fichier, etiquette, cible)
    ("hf_aniso", "aniso., pompe dès $t=0$", CIBLE["aniso"]),
    ("e3_aniso", "aniso., protocole article", CIBLE["aniso"]),
    ("e2_aniso", r"aniso., $D_{\mathrm{wet}}=0$", CIBLE["aniso"]),
    ("hf_iso", "iso. 4,6, pompe dès $t=0$", CIBLE["iso46"]),
    ("e1_iso_cible12", "iso. 3,5, pompe dès $t=0$", CIBLE["iso35"]),
    ("e3_iso12", "iso. 3,5, protocole article", CIBLE["iso35"]),
    ("e2_iso12", r"iso. 3,5, $D_{\mathrm{wet}}=0$", CIBLE["iso35"]),
]


def pic(h):
    i = int(np.argmax(h["hydroP"]))
    out = dict(p=h["hydroP"][i] / 1e6, t=h["t"][i] * 1e3, ins=None)
    if "nInserted" in h and h["nInserted"].max() > 0:
        j = int(np.argmax(h["nInserted"] > 0))
        out["ins"] = h["hydroP"][j] / 1e6
    return out


def fig_pression():
    H = {n: read_hist(os.path.join(HIST, n + ".csv")) for n, _, _ in RUNS_ARCH}
    rej = {}
    for n in ("hf_aniso_hydro_c", "hf_iso_hydro_c",
              "hf_aniso_hydro_m6", "hf_iso_hydro_m6"):
        p = os.path.join(RUNS, "out_" + n, "history.csv")
        if os.path.exists(p):
            rej[n] = read_hist(p)

    fig, (a, b) = plt.subplots(1, 2, figsize=(7.4, 3.9),
                               gridspec_kw=dict(width_ratios=[1.15, 1]))
    # (a) p(t), runs de reference + rejeux grossiers eventuels
    sty = {"hf_aniso": ("C0", "-", r"aniso. 6,8/4,6, maille 3 mm (20/08)"),
           "hf_iso": ("C3", "-", r"iso. 4,6, maille 3 mm (20/08)")}
    for n, (c, ls, lab) in sty.items():
        h = H[n]
        a.plot(h["t"] * 1e3, h["hydroP"] / 1e6, color=c, ls=ls, lw=1.4,
               label=lab)
        k = pic(h)
        a.plot(k["t"], k["p"], "o", color=c, ms=4)
        a.annotate(virg(k["p"], 2), (k["t"], k["p"]), xytext=(-30, 4),
                   textcoords="offset points", color=c, fontsize=8.5)
    # rejeux du 2026-10-06 (binaire f0209ef) sur maillages plus grossiers
    lab_rej = {"hf_aniso_hydro_c": ("C0", ":", "aniso., maille 12 mm"),
               "hf_iso_hydro_c": ("C3", ":", "iso., maille 12 mm"),
               "hf_aniso_hydro_m6": ("C0", "-.", "aniso., maille 6 mm"),
               "hf_iso_hydro_m6": ("C3", "-.", "iso., maille 6 mm")}
    for n, h in rej.items():
        c, ls, lab = lab_rej[n]
        a.plot(h["t"] * 1e3, h["hydroP"] / 1e6, color=c, ls=ls, lw=1.0,
               label=lab + " (rejoué 06/10)")
    a.axhline(CIBLE["aniso"], color="C0", ls="--", lw=0.8)
    a.axhline(CIBLE["iso46"], color="C3", ls="--", lw=0.8)
    a.text(0.05, CIBLE["aniso"] + 0.2, "éq. 10 : 12,0 MPa", color="C0",
           fontsize=8)
    a.text(0.05, CIBLE["iso46"] + 0.2, "éq. 10 : 14,2 MPa", color="C3",
           fontsize=8)
    # leur valeur numerique : ~12,5 MPa (texte, apres eq. 10), 11,69 lu fig. 11b
    a.axhspan(11.69, 12.5, color="0.85", zorder=0)
    a.text(0.05, 10.6, "gris : Y-Geo, 11,69 à 12,5",
           fontsize=7.5, color="0.35")
    a.set_xlabel("temps depuis le début du calcul [ms]")
    a.set_ylabel("pression de puits $p$ [MPa]")
    a.set_xlim(0, 4.0)
    a.set_ylim(0, 19.5)
    a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2,
             frameon=False, fontsize=7.0)
    a.set_title("(a) Pression de puits contre les seuils fermés", loc="left")
    fmt_axes(a)

    # (b) pic decompose, en ecart relatif a la cible de l'eq. 10
    y = np.arange(len(RUNS_ARCH))[::-1]
    for yi, (n, lab, cib) in zip(y, RUNS_ARCH):
        k = pic(H[n])
        tot = (k["p"] - cib) / cib * 100
        if k["ins"] is not None:
            sta = (k["ins"] - cib) / cib * 100
            b.barh(yi, sta, color="#e39a63", height=0.6)
            b.barh(yi, tot - sta, left=sta, color="#7fb3d5", height=0.6)
        else:
            b.barh(yi, tot, color="0.7", height=0.6)
        b.text(tot + 0.4, yi, "+" + virg(tot, 1) + " %", va="center",
               fontsize=8)
    b.set_yticks(y)
    b.set_yticklabels([r[1] for r in RUNS_ARCH], fontsize=8)
    b.set_xlim(0, 44)
    b.set_ylim(-0.6, len(RUNS_ARCH) - 0.2)
    b.axvline(4.2, color="0.4", ls=":", lw=0.8)  # 12,5 / 12 - 1 = +4,2 %
    b.text(4.6, y[0] + 0.45, "Y-Geo, +4 %", fontsize=7.5, color="0.35")
    from matplotlib.patches import Patch
    b.legend(bbox_to_anchor=(1.0, -0.17), handles=[Patch(color="#e39a63", label="statique (cible → insertion)"),
                      Patch(color="#7fb3d5", label="incubation (insertion → pic)"),
                      Patch(color="0.7", label="non instrumenté")],
             loc="upper right", frameon=False, fontsize=7.0)
    b.set_xlabel("dépassement du seuil de l'éq. 10 [%]")
    b.set_title("(b) Dépassement décomposé, 7 calculs", loc="left")
    fmt_axes(b)
    b.set_yticklabels([r[1] for r in RUNS_ARCH], fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_abuaisha_pression.pdf"),
                bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ rejeux
CX, CY, RB = 4.0, 4.0, 0.05
E, NU = 35.0e9, 0.27            # leur Table 1
G = E / (2 * (1 + NU))


def wall_history(run):
    fs = frames(run)
    P0, _ = read_vtu_points(fs[0])
    r0 = np.hypot(P0[:, 0] - CX, P0[:, 1] - CY)
    w = np.where(np.abs(r0 - RB) < 2.5e-3)[0]
    ur = []
    for f in fs:
        P, _ = read_vtu_points(f)
        ur.append((np.hypot(P[w, 0] - CX, P[w, 1] - CY) - r0[w]).mean())
    tf = read_hist(os.path.join(run, "frames.csv"))["t"][:len(fs)]
    return tf, np.array(ur)


def parker_profile(run, l=0.75):
    """ouverture SIGNEE (levre haute - levre basse), cote fixe par le
    centroide de l'element porteur — meme methode que vv/B_articles
    (m_parker) : un controle par une norme ne controle pas le signe."""
    fs = frames(run)
    P0, C = read_vtu_points(fs[0])
    P, _ = read_vtu_points(fs[-1])
    side = np.zeros(len(P0))
    cyel = P0[C].mean(axis=1)[:, 1]
    for k in range(3):
        side[C[:, k]] = np.sign(cyel - CY)
    on = (np.abs(P0[:, 1] - CY) < 1e-9) & (np.abs(P0[:, 0] - CX) <= l + 1e-9)
    xs = np.round(P0[on, 0] - CX, 9)
    ys, ss = P[on, 1], side[on]
    xa, wa = [], []
    for xv in np.unique(xs):
        g = xs == xv
        up, lo = ys[g & (ss > 0)], ys[g & (ss < 0)]
        if len(up) and len(lo):
            xa.append(xv)
            wa.append(up.mean() - lo.mean())
    tf = read_hist(os.path.join(run, "frames.csv"))["t"][len(fs) - 1]
    return np.array(xa), np.array(wa), tf, len(fs)


def fig_rejeu():
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    # (a) Lame : u_r = p a / 2G, deux chemins de chargement
    hh = read_hist(os.path.join(RUNS, "out_signe_hydro", "history.csv"))
    for run, c, mk, lab in (("out_signe_conf", "C2", "o",
                             "pression imposée (confinement)"),
                            ("out_signe_hydro", "C0", "x",
                             "module hydro (pression imposée)")):
        t, ur = wall_history(os.path.join(RUNS, run))
        a.plot(t * 1e3, ur * 1e6, mk, color=c, ms=5, mfc="none", label=lab)
    p = hh["hydroP"]
    a.plot(hh["t"] * 1e3, p * RB / (2 * G) * 1e6, "k-", lw=1,
           label=r"Lamé $u_r = p\,a/2G$")
    t, ur = wall_history(os.path.join(RUNS, "out_signe_hydro"))
    err = (ur[-1] - 12e6 * RB / (2 * G)) / (12e6 * RB / (2 * G)) * 100
    a.text(0.03, 0.62, "final : %s µm\nLamé : %s µm (%s %%)"
           % (virg(ur[-1] * 1e6, 3), virg(12e6 * RB / (2 * G) * 1e6, 3),
              virg(err, 1)), transform=a.transAxes, fontsize=8)
    a.set_xlabel("temps [ms]")
    a.set_ylabel(r"déplacement radial moyen de paroi [µm]")
    a.set_title("(a) Paroi du forage contre Lamé", loc="left")
    a.legend(loc="lower right", frameon=False, fontsize=7.2)
    fmt_axes(a)

    # (b) Parker : leur eq. A.1, w = 2 s'(1-nu^2)/E sqrt(l^2-x^2) par levre
    run = os.path.join(RUNS, "out_parker_hydro_c")
    if os.path.isdir(run) and frames(run):
        xa, wa, tf, nf = parker_profile(run)
        Ep, nup, sp, l = 45e9, 0.2, 2e6, 0.75  # annexe A (E = 45 GPa)
        xx = np.linspace(-l, l, 400)
        w1 = 2 * sp * (1 - nup ** 2) / Ep * np.sqrt(l * l - xx ** 2)
        b.plot(xx, 2 * w1 * 1e3, "k-", lw=1,
               label=r"$2w$, éq. A.1 (Parker) : $2w(0)=0{,}128$ mm")
        b.plot(xx, w1 * 1e3, "k--", lw=0.8,
               label=r"$w$, éq. A.1 : $w(0)=0{,}064$ mm")
        b.plot(xa, wa * 1e3, ".", color="C3", ms=3,
               label="rockim, ouverture signée")
        i0 = int(np.argmin(np.abs(xa)))
        w0an = 2 * 2 * sp * (1 - nup ** 2) / Ep * l
        b.text(0.5, 0.03, "centre : %s mm, soit %s %% de $2w(0)$\n"
               "$t$ = %s ms (trame %d)"
               % (virg(wa[i0] * 1e3, 4), virg((wa[i0] / w0an - 1) * 100, 1),
                  virg(tf * 1e3, 1), nf - 1),
               transform=b.transAxes, fontsize=8, ha="center")
        b.legend(loc="upper right", frameon=False, fontsize=7.0)
    b.set_xlabel("abscisse depuis le centre de la fissure [m]")
    b.set_ylabel("ouverture [mm]")
    b.set_ylim(0, 0.17)
    b.set_title("(b) Fissure de Parker sous 2 MPa nets", loc="left")
    fmt_axes(b)
    fig.text(0.5, -0.01, "calculs rejoués sur maillages grossiers "
             "(maille fine de 12 mm)", ha="center",
             fontsize=7.5, color="0.35")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_abuaisha_rejeu.pdf"),
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig_pression()
    fig_rejeu()
    print("ecrit :", os.path.join(HERE, "fig_abuaisha_pression.pdf"),
          os.path.join(HERE, "fig_abuaisha_rejeu.pdf"))
