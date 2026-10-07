# -*- coding: utf-8 -*-
"""Depouillement des runs fem3d de l'etude (phase A et suivantes).

usage : python extract_fem3d.py <dossier_runs> [out.csv]
  <dossier_runs> contient des sous-dossiers out_* (history.csv, frames.csv,
  fem3d_NNNN.vtu ASCII). Pour chaque run :
    e_r      = KE_outil(fin) / KE_outil(0)               restitution
    dmax     = z_outil(0) - min z_outil                   enfoncement max [mm]
    dz_fin   = z_outil(fin) - min z_outil                 remontee (rebond) [mm]
    W_abs    = travail de l'outil (fin) [J]
    F_pic    = max |Fz| [kN]  (brut, non filtre : rang 2)
    nEroded, craterVol [mm3]
    V(D>=0,9), V(D>=0,5), V(epvEq>=0,05), V(erodes) [mm3]  (dernier VTU, cellules)
  Les colonnes ajoutees par fieldStats (wPlas, wDamT, wDamC, ...) sont
  reprises telles quelles si presentes (fin de ligne)."""
import io, os, re, sys, csv, glob
import numpy as np


def read_history(path):
    with io.open(path, encoding="utf-8") as f:
        head = f.readline().strip().split(",")
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    if data.ndim == 1:
        data = data[None, :]
    return head, data


def read_vtu_cells(path):
    """VTU ASCII ecrit par VtkWriter : points, tets, champs cellule."""
    txt = io.open(path, encoding="utf-8", errors="ignore").read()
    def arr(name):
        m = re.search(r'Name="%s"[^>]*>\s*(.*?)\s*</DataArray>' % re.escape(name), txt, re.S)
        return np.fromstring(m.group(1), sep=" ") if m else None
    pts = arr("Points") if arr("Points") is not None else None
    if pts is None:
        m = re.search(r"<Points>\s*<DataArray[^>]*>\s*(.*?)\s*</DataArray>", txt, re.S)
        pts = np.fromstring(m.group(1), sep=" ")
    pts = pts.reshape(-1, 3)
    con = arr("connectivity").astype(int).reshape(-1, 4)
    fields = {}
    for name in re.findall(r'<DataArray type="Float64" Name="([^"]+)"', txt):
        if name in ("Points", "velocity"):
            continue
        a = arr(name)
        if a is not None and a.size == con.shape[0]:
            fields[name] = a
    # volumes courants des tets
    a, b, c, d = (pts[con[:, k]] for k in range(4))
    vol = np.abs(np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a)) / 6.0
    return con, vol, fields


def summarize(run_dir):
    head, h = read_history(os.path.join(run_dir, "history.csv"))
    col = {n: i for i, n in enumerate(head)}
    out = {"run": os.path.basename(run_dir)}
    t = h[:, col["t"]]
    if "toolKE" in col:
        ke = h[:, col["toolKE"]]
        z = h[:, col["toolZ"]]
        out["e_r"] = ke[-1] / ke[0] if ke[0] > 0 else float("nan")
        out["dmax_mm"] = 1e3 * (z[0] - z.min())
        out["dz_fin_mm"] = 1e3 * (z[-1] - z.min())
        out["W_abs_J"] = h[-1, col["work"]]
        out["F_pic_kN"] = 1e-3 * np.abs(h[:, col["toolFz"]]).max()
        out["nEroded"] = int(h[-1, col["nEroded"]])
        out["craterVol_mm3"] = 1e9 * h[-1, col["craterVol"]]
        out["t_fin_us"] = 1e6 * t[-1]
        for k in head[10:]:                       # colonnes ajoutees (fieldStats)
            out[k] = h[-1, col[k]]
        # coupe (scenario shear) : fenetre stationnaire x in [4, 12] mm
        if "toolFx" in col and "toolX" in col:
            x = h[:, col["toolX"]]
            fx = h[:, col["toolFx"]]
            fz = h[:, col["toolFz"]]
            if x.max() - x.min() > 5e-3:          # l'outil a avance (coupe)
                sel = (x >= 4e-3) & (x <= 12e-3)
                if sel.sum() > 10:
                    out["Fx_mean_kN"] = -1e-3 * fx[sel].mean()
                    out["Fx_cv"] = fx[sel].std() / max(abs(fx[sel].mean()), 1e-9)
                    out["Fz_mean_kN"] = -1e-3 * fz[sel].mean()
                    out["Fx_max_kN"] = -1e-3 * fx[sel].min()
                # travail de coupe sur la fenetre / volume balaye (d x D x L) = MSE
                w = h[:, col["work"]]
                i0 = np.argmax(x >= 4e-3) if (x >= 4e-3).any() else None
                i1 = np.argmax(x >= 12e-3) if (x >= 12e-3).any() else None
                if i0 is not None and i1 is not None and i1 > i0:
                    out["W_coupe_J"] = w[i1] - w[i0]
                out["x_fin_mm"] = 1e3 * x[-1]
                out["atteint_12mm"] = bool(x[-1] >= 12e-3)
    vtus = sorted(glob.glob(os.path.join(run_dir, "fem3d_*.vtu")))
    if vtus:
        con, vol, f = read_vtu_cells(vtus[-1])
        D = f.get("damage")
        ero = f.get("eroded")
        epv = f.get("epvEq")
        alive = ero < 0.5 if ero is not None else np.ones(len(vol), bool)
        if D is not None:
            out["V_D09_mm3"] = 1e9 * vol[(D >= 0.9) & alive].sum()
            out["V_D05_mm3"] = 1e9 * vol[(D >= 0.5) & alive].sum()
        if epv is not None:
            out["V_epv005_mm3"] = 1e9 * vol[(epv >= 0.05) & alive].sum()
            out["epv_max"] = epv.max()
        if ero is not None:
            out["V_eroded_mm3"] = 1e9 * vol[~alive].sum()
        if "omegaC" in f:
            out["V_wc05_mm3"] = 1e9 * vol[(f["omegaC"] >= 0.5) & alive].sum()
        out["vtu"] = os.path.basename(vtus[-1])
    log = os.path.join(os.path.dirname(run_dir), "logs", out["run"].replace("out_", "") + ".log")
    if os.path.exists(log):
        s = io.open(log, encoding="utf-8", errors="ignore").read()
        m = re.search(r"wall time: ([0-9.]+)", s)
        out["wall_s"] = float(m.group(1)) if m else float("nan")
        m = re.search(r"dt = ([0-9.e+-]+) s, steps = (\d+)", s)
        if m:
            out["dt_s"] = float(m.group(1)); out["steps"] = int(m.group(2))
        out["crash"] = "error" in s.lower()
    return out


def main():
    root = sys.argv[1]
    outcsv = sys.argv[2] if len(sys.argv) > 2 else os.path.join(root, "resultats.csv")
    runs = sorted(d for d in glob.glob(os.path.join(root, "out_*")) if os.path.isdir(d))
    rows = [summarize(d) for d in runs]
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with io.open(outcsv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    # tableau lisible
    show = [k for k in ("run", "e_r", "dmax_mm", "dz_fin_mm", "W_abs_J", "F_pic_kN", "nEroded",
                        "craterVol_mm3", "V_D09_mm3", "V_D05_mm3", "V_epv005_mm3", "epv_max", "wall_s") if k in keys]
    print(" | ".join("%-14s" % k for k in show))
    for r in rows:
        cells = []
        for k in show:
            v = r.get(k, "")
            cells.append("%-14s" % (("%.4g" % v) if isinstance(v, float) else str(v)))
        print(" | ".join(cells))
    print("->", outcsv)


if __name__ == "__main__":
    main()
