# -*- coding: utf-8 -*-
"""GIF de l'impact banc moyen (coupe mediane) — gere le multi-corps
toolShape=none que tools/make_gif.py ne sait pas rendre (sentinelle 1e9).
Roche coloree par von Mises, insert en gris, fragments detaches en rouge."""
import io, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

snap = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "out_banc_mid")
out_gif = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "banc_mid_impact.gif")

def arr(text, name):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name, text, re.S)
    return np.fromstring(m.group(1), sep=" ")

def points(text):
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", text, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)

hist = np.genfromtxt(os.path.join(snap, "history.csv"), delimiter=",", names=True, invalid_raise=False)
frames_meta = np.genfromtxt(os.path.join(snap, "frames.csv"), delimiter=",", names=True, invalid_raise=False)

vtus = sorted(f for f in os.listdir(snap)
              if re.fullmatch(r"fdem3d_\d{4}\.vtu", f))
imgs = []
for k, fn in enumerate(vtus):
    txt = open(os.path.join(snap, fn)).read()
    P = points(txt)
    conn = arr(txt, "connectivity").astype(int).reshape(-1, 4)
    cen = P[conn].mean(axis=1)
    phase = arr(txt, "phase")
    frag = arr(txt, "fragment")
    vm = arr(txt, "vonMises")

    tk = float(np.atleast_1d(frames_meta["t"])[k])
    Fz = float(np.interp(tk, hist["t"], hist["grpFz"]))
    vz = float(np.interp(tk, hist["t"], hist["grpVz"]))

    sl = np.abs(cen[:, 1] - 0.06) < 0.005          # coupe mediane y = D/2
    rock = sl & (phase < 0.5)
    ins = sl & (phase >= 0.5)
    # fragment principal de la roche = id majoritaire ; le reste = detache
    ids, counts = np.unique(frag[phase < 0.5], return_counts=True)
    main = ids[np.argmax(counts)]
    det = rock & (frag != main)

    fig, ax = plt.subplots(figsize=(7, 6.2))
    sc = ax.scatter(cen[rock, 0], cen[rock, 2], c=vm[rock] * 1e-6,
                    cmap="viridis", vmin=0, vmax=40, s=14, marker="s")
    ax.scatter(cen[ins, 0], cen[ins, 2], color="0.45", s=8, marker="o")
    if det.any():
        ax.scatter(cen[det, 0], cen[det, 2], color="red", s=18, marker="s",
                   label="fragment détaché")
        ax.legend(loc="lower left", fontsize=8)
    ax.set_xlim(0, 0.12); ax.set_ylim(0, 0.135)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)"); ax.set_ylabel("z (m)")
    ax.set_title("t = %.0f µs   Fz = %.1f kN   vz = %+.2f m/s" % (tk * 1e6, Fz * 1e-3, vz))
    cb = fig.colorbar(sc, ax=ax, shrink=0.8); cb.set_label("von Mises (MPa)")
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100)
    plt.close(fig)
    buf.seek(0)
    imgs.append(Image.open(buf).convert("P"))

imgs[0].save(out_gif, save_all=True, append_images=imgs[1:],
             duration=700, loop=0)
print("ecrit", out_gif, "(%d frames)" % len(imgs))
