# -*- coding: utf-8 -*-
"""Le jeu de la base LHS le plus proche de l'experimental : facies des trois
essais + comparaison chiffree aux cibles Dumoulin et al. (2024)."""
import csv, glob, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
RUNS = os.path.join(BASE, "runs")
T = {"ucs": 126.6, "bts": 10.27, "tx20": 424.8}
TAG = sys.argv[1] if len(sys.argv) > 1 else "L001"


def arr(t, n):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n, t, re.S)
    return np.fromstring(m.group(1), sep=" ") if m else None


def pts(t):
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", t, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)[:, :2]


def draw(ax, run, title):
    d = os.path.join(RUNS, run)
    fs = sorted(glob.glob(os.path.join(d, "fdem_[0-9][0-9][0-9][0-9].vtu")))
    if not fs:
        ax.set_axis_off(); return
    t0 = open(fs[0]).read()
    conn = arr(t0, "connectivity").astype(int).reshape(-1, 3)
    grain = arr(t0, "grain")
    P = pts(open(fs[-1]).read())
    rng = np.random.default_rng(3)
    lut = rng.random((int(grain.max()) + 2, 3)) * 0.35 + 0.6
    ax.add_collection(PolyCollection([P[c] * 1e3 for c in conn],
                                     facecolors=[lut[int(g)] for g in grain],
                                     edgecolors="0.5", linewidths=0.1))
    jf = fs[-1].replace("fdem_", "fdem_joints_")
    n = 0
    if os.path.exists(jf):
        jt = open(jf).read()
        jP, jc = pts(jt), arr(jt, "connectivity").astype(int)
        tb, off = arr(jt, "tBreak"), arr(jt, "offsets").astype(int)
        segs, s = [], 0
        for i, e in enumerate(off):
            idx = jc[s:e]; s = e
            if tb is not None and tb[i] >= 0 and len(idx) >= 2:
                segs.append(jP[idx[:2]] * 1e3)
        n = len(segs)
        if segs:
            ax.add_collection(LineCollection(segs, colors="crimson", lw=1.3))
    P0 = pts(t0)
    ax.set_xlim(P0[:, 0].min() * 1e3 - 2, P0[:, 0].max() * 1e3 + 2)
    ax.set_ylim(P0[:, 1].min() * 1e3 - 2, P0[:, 1].max() * 1e3 + 2)
    ax.set_aspect("equal"); ax.set_xlabel("x (mm)")
    ax.set_title("%s\n%d joints rompus" % (title, n), fontsize=10)


def main():
    row = [r for r in csv.DictReader(open(os.path.join(BASE, "lhs_results.csv")))
           if r["tag"] == TAG][0]
    sim = {"ucs": float(row["ucs_peak_MPa"]),
           "bts": float(row["bts_sigma_t_MPa"]),
           "tx20": float(row["tx20_peak_MPa"])}

    fig = plt.figure(figsize=(14, 6.4))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.15])
    draw(fig.add_subplot(gs[0]), "%s_ucs_s4211" % TAG,
         "Compression simple\n%.1f MPa (cible %.1f)" % (sim["ucs"], T["ucs"]))
    draw(fig.add_subplot(gs[1]), "%s_bts_s4211" % TAG,
         "Brésilien\n%.2f MPa (cible %.2f)" % (sim["bts"], T["bts"]))
    draw(fig.add_subplot(gs[2]), "%s_tx20_s4211" % TAG,
         "Triaxial $\\sigma_3$ = 20 MPa\n%.1f MPa (cible %.1f)" % (sim["tx20"], T["tx20"]))

    a = fig.add_subplot(gs[3])
    names = ["UCS", "$\\sigma_t$ (BTS)", "$\\sigma_1$ à $\\sigma_3$=20"]
    y = np.arange(3)
    s = [sim["ucs"], sim["bts"], sim["tx20"]]
    t = [T["ucs"], T["bts"], T["tx20"]]
    a.barh(y + 0.18, [si / ti * 100 for si, ti in zip(s, t)], height=0.34,
           color="C0", label="rockim")
    a.barh(y - 0.18, [100] * 3, height=0.34, color="C2", alpha=0.65,
           label="expérimental")
    for i, (si, ti) in enumerate(zip(s, t)):
        a.text(si / ti * 100 + 3, i + 0.18, "%.1f" % si, va="center", fontsize=9)
        a.text(103, i - 0.18, "%.1f" % ti, va="center", fontsize=9)
    a.set_yticks(y); a.set_yticklabels(names)
    a.set_xlabel("% de la cible expérimentale")
    a.axvline(100, color="k", lw=0.8)
    a.set_xlim(0, 160); a.legend(fontsize=9); a.grid(alpha=0.3, axis="x")
    a.set_title("Écarts aux essais Dumoulin et al. (2024)", fontsize=10)

    p = " ".join("%s=%.3g" % (k, float(row[k])) for k in
                 ("ft", "cohesion", "frictionDeg", "Gf", "gfShearFactor"))
    fig.suptitle("Jeu %s — le plus proche de l'expérimental dans la base LHS\n%s"
                 % (TAG, p), fontsize=11)
    fig.tight_layout()
    out = os.path.join(BASE, "figures", "best_%s.png" % TAG)
    fig.savefig(out, dpi=140)
    print("ecrit", out)


if __name__ == "__main__":
    main()
