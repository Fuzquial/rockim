#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_decks.py — du plan d experiences aux decks rockim + liste de jobs.
#
#   python calib_quick/calib/make_decks.py design.csv --template calib_quick/q1u070_P050.cfg
#          --tag h1 --conf 20 50 [--seeds 1] [--outdir calib_quick/runs_h1]
#
# Chaque ligne du design (colonnes = parametres physiques + id) devient un
# deck par confinement (et par seed) : le template est recopie, les cles
# listees dans MAPPING sont remplacees. Parametres derives :
#   Gf = lcz * ft^2 / E          (l_cz = E Gf / ft^2, tenu explicitement)
#   phi -> frictionDeg, muRes -> jointResidualMu, gII -> gfShearFactor
# Ecrit <outdir>/<tag>_<id>_P<conf>[_s<seed>].cfg et <outdir>/jobs_<tag>.json.
# NE LANCE RIEN.
# ---------------------------------------------------------------------------
import argparse
import json
import os
import re
import sys

import csv

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))


def read_csv(path):
    """design.csv (ecrit par ParamSpace.write_csv) -> (noms des parametres, lignes dict)"""
    with open(path, newline="", encoding="ascii") as fh:
        rows = list(csv.reader(fh))
    header = [h.strip() for h in rows[0]]
    names = [h for h in header if h != "id"]
    out = []
    for r in rows[1:]:
        if not r or all(not c.strip() for c in r):
            continue
        d = {h: (float(v) if h != "id" else v.strip()) for h, v in zip(header, r)}
        out.append(d)
    return names, out

# parametre du design -> cle rockim (valeur formatee). Les cles absentes du
# design gardent la valeur du template.
MAPPING = {
    "ft": ("ft", lambda v: f"{v:.6g}"),
    "c": ("cohesion", lambda v: f"{v:.6g}"),
    "phi": ("frictionDeg", lambda v: f"{v:.4g}"),
    "gII": ("gfShearFactor", lambda v: f"{v:.4g}"),
    "muRes": ("jointResidualMu", lambda v: f"{v:.4g}"),
    "E": ("E", lambda v: f"{v:.6g}"),
    "nu": ("nu", lambda v: f"{v:.4g}"),
    "m": ("jointWeibullM", lambda v: f"{v:.4g}"),
    "corrLen": ("strengthCorrLength", lambda v: f"{v:.6g}"),
    "pen": ("jointPenaltyFactor", lambda v: f"{v:.4g}"),
    "Gf": ("Gf", lambda v: f"{v:.6g}"),            # Gf direct (sinon derive de lcz)
    "muRes": ("jointResidualMu", lambda v: f"{v:.4g}"),
}


def set_key(lines, key, value, comment=""):
    pat = re.compile(r"^\s*" + re.escape(key) + r"\s*=")
    new = f"{key} = {value}" + (f"    # {comment}" if comment else "")
    for i, l in enumerate(lines):
        if pat.match(l):
            lines[i] = new
            return lines
    lines.append(new)
    return lines


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("design")
    ap.add_argument("--template", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--conf", type=float, nargs="+", default=[20, 50])
    ap.add_argument("--seeds", type=int, nargs="*", default=[])
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--extra", nargs="*", default=[], help="cle=valeur imposees a tous les decks")
    ap.add_argument("--bts", default=None, help="template du bresilien (un deck BTS par point, memes joints)")
    a = ap.parse_args()
    outdir = a.outdir or os.path.join(ROOT, "calib_quick", f"runs_{a.tag}")
    os.makedirs(outdir, exist_ok=True)
    tpl = open(a.template, encoding="utf-8", errors="replace").read().splitlines()
    tplB = open(a.bts, encoding="utf-8", errors="replace").read().splitlines() if a.bts else None
    names, rows = read_csv(a.design)
    jobs = []
    for irow, p in enumerate(rows):
        pid = str(p.get("id", f"d{irow:03d}"))
        if tplB is not None:
            lines = list(tplB)
            lines[0] = f"# {a.tag} point {pid} BTS (bresilien) - " + ", ".join(f"{k} = {p[k]:.4g}" for k in names)
            for k, (key, fmt) in MAPPING.items():
                if k in p:
                    lines = set_key(lines, key, fmt(p[k]))
            if "lcz" in p and "ft" in p and "Gf" not in p:
                m = [l for l in lines if re.match(r"^\s*E\s*=", l)]
                E = p.get("E") or (float(m[0].split("=")[1].split("#")[0]) if m else 77.7e9)
                lines = set_key(lines, "Gf", f"{p['lcz'] * p['ft'] ** 2 / E:.6g}", f"l_cz = {p['lcz']*1e3:.3g} mm")
            for kv in a.extra:
                k, v = kv.split("=", 1)
                if k.strip() in ("historyStrains", "gripsStopAfterPeak", "gripsStopDelay", "stopPeakDrop"):
                    continue                      # cles du triaxial, sans objet en bresilien
                lines = set_key(lines, k.strip(), v.strip())
            name = f"{a.tag}_{pid}_BTS"
            cfg = os.path.join(outdir, name + ".cfg")
            open(cfg, "w", encoding="ascii", errors="replace").write("\n".join(lines) + "\n")
            jobs.append({"cfg": os.path.relpath(cfg, ROOT), "out": os.path.relpath(os.path.join(outdir, "out_" + name), ROOT),
                         "id": pid, "conf": "BTS", "seed": None})
        for conf in a.conf:
            seeds = a.seeds or [None]
            for s in seeds:
                lines = list(tpl)
                lines[0] = f"# {a.tag} point {pid} sigma3 = {conf:g} MPa" + (f" seed {s}" if s is not None else "") + \
                           " - " + ", ".join(f"{k} = {p[k]:.4g}" for k in names)
                for k, (key, fmt) in MAPPING.items():
                    if k in p:
                        lines = set_key(lines, key, fmt(p[k]))
                if "lcz" in p and "ft" in p:
                    E = p.get("E")
                    if E is None:
                        m = [l for l in lines if re.match(r"^\s*E\s*=", l)]
                        E = float(m[0].split("=")[1].split("#")[0]) if m else 71e9
                    Gf = p["lcz"] * p["ft"] ** 2 / E
                    lines = set_key(lines, "Gf", f"{Gf:.6g}", f"l_cz = {p['lcz']*1e3:.3g} mm")
                lines = set_key(lines, "confiningPressure", f"{conf:g}e6")
                if s is not None:
                    lines = set_key(lines, "seed", str(s))
                for kv in a.extra:
                    k, v = kv.split("=", 1)
                    lines = set_key(lines, k.strip(), v.strip())
                name = f"{a.tag}_{pid}_P{int(conf):03d}" + (f"_s{s}" if s is not None else "")
                cfg = os.path.join(outdir, name + ".cfg")
                open(cfg, "w", encoding="ascii", errors="replace").write("\n".join(lines) + "\n")
                jobs.append({"cfg": os.path.relpath(cfg, ROOT), "out": os.path.relpath(os.path.join(outdir, "out_" + name), ROOT),
                             "id": pid, "conf": conf, "seed": s})
    jp = os.path.join(outdir, f"jobs_{a.tag}.json")
    json.dump(jobs, open(jp, "w"), indent=1)
    print(f"{len(rows)} points x {len(a.conf)} confinements x {max(1, len(a.seeds))} seeds = {len(jobs)} decks -> {outdir}")
    print("jobs :", jp)


if __name__ == "__main__":
    main()
