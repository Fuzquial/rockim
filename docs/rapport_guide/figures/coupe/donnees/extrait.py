#!/usr/bin/env python3
"""Extraction des observables de coupe d'un run rockim fdem (2D, scenario shear).

usage: extrait.py out_dir log [vOutil=10]
Imprime : ratio d'injection, v max / 2v, joints rompus, fragments, residu B4,
force de coupe (pic |Fx|, pic |F|, moyenne de |Fx| en contact, fraction du temps
en contact depuis le premier contact), comparees a la cible Heilman 3,08 MN/m.
"""
import re
import sys

import numpy as np

out, log = sys.argv[1], sys.argv[2]
v = float(sys.argv[3]) if len(sys.argv) > 3 else 10.0
txt = open(log, encoding="utf-8", errors="replace").read()


def grab(pat, cast=float):
    m = re.search(pat, txt)
    return cast(m.group(1)) if m else None


res = {
    "outil->solide": grab(r"outil->solide:\s*([-\d.e+]+)"),
    "tool work output": grab(r"tool work output\s*:\s*([-\d.e+]+)"),
    "injection ratio (imprime)": grab(r"injection outil\s*:.*ratio\s*([-\d.e+]+)"),
    "injection trapeze ratio": grab(r"injection \(trapeze\):.*ratio\s*([-\d.e+]+)"),
    "v nodale max": grab(r"v nodale max\s*:\s*([-\d.e+]+)"),
    "integration": grab(r"integration\s*:\s*\+?([-\d.e+]+)"),
    "Cundall": grab(r"Cundall\s*:\s*([-\d.e+]+)"),
    "contact general": grab(r"contact\s*:\s*([-\d.e+]+) J/m \(dont"),
    "KE finale": grab(r"KE 0 -> ([-\d.e+]+)"),
    "residu": grab(r"residu\s*:\s*([-\d.e+]+) J/m"),
    "residu pct": grab(r"residu.*\(([-\d.e+]+) % de"),
    "broken": grab(r"broken joints\s*:\s*(\d+)", int),
    "fragments": grab(r"fragments\s*:\s*(\d+)", int),
    "peak tool force (log)": grab(r"peak tool force\s*:\s*([-\d.e+]+)"),
    "noeuds profonds": grab(r"noeuds profonds\s*:\s*(\d+)", int),
    "collees pct": grab(r"evaluations COLLEES \(([-\d.e+]+) %\)"),
    "wall": grab(r"wall time:\s*([-\d.e+]+)"),
}
if res["outil->solide"] and res["tool work output"]:
    res["ratio toolWork/work"] = res["outil->solide"] / res["tool work output"]
if res["v nodale max"]:
    res["vmax/2v"] = res["v nodale max"] / (2 * v)

H = np.genfromtxt(f"{out}/history.csv", delimiter=",", names=True)
t, fx, fy, x = H["t"], H["toolFx"], H["toolFy"], H["toolX"]
F = np.hypot(fx, fy)
on = F > 1.0                                   # N/m : contact present
if on.any():
    i0 = np.argmax(on)
    sl = slice(i0, None)
    res["premier contact t (s)"] = float(t[i0])
    res["premier contact toolX (mm)"] = float(x[i0] * 1e3)
    res["fraction contact apres 1er contact"] = float(on[sl].mean())
    res["pic |Fx| (MN/m)"] = float(np.abs(fx).max() / 1e6)
    res["pic |F| (MN/m)"] = float(F.max() / 1e6)
    res["moyenne |Fx| en contact (MN/m)"] = float(np.abs(fx[on]).mean() / 1e6)
    res["mediane |Fx| en contact (MN/m)"] = float(np.median(np.abs(fx[on])) / 1e6)
    res["moyenne |Fx| depuis 1er contact (MN/m)"] = float(np.abs(fx[sl]).mean() / 1e6)
    res["releves > 1 MN/m"] = f"{int((F > 1e6).sum())} / {len(F)}"
    res["pic |Fx| / 3,08"] = res["pic |Fx| (MN/m)"] / 3.08
for k, val in res.items():
    print(f"{k:40s} {val}")
