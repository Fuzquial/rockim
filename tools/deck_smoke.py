#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# deck_smoke.py — lancement de FUMEE d un ou plusieurs decks (T4, campagne du
# 13/09) : chaque deck est recopie dans un dossier de travail avec T et frames
# SURCHARGES (2 us / 1 trame par defaut), lance avec OMP_NUM_THREADS fixe, et le
# journal est depouille : demarrage sans erreur de cle ?, dt, nombre de pas,
# masses par corps, avertissements, temps mur. Il en deduit le COUT ESTIME du
# deck complet (pas = T_deck / dt) a un cout par pas donne.
#
#   python tools/deck_smoke.py --exe rockim_g1y16.exe --out tests_f2/campagne13/T4 \
#          [--T 2e-6] [--frames 1] [--threads 4] [--mesh meshes/x.msh] \
#          [--drop-key jointBreakModeRef] [--tag nom] [--timeout 600] \
#          [--ms-shared 77.73] [--ms-alone 17.0] deck1.cfg deck2.cfg ...
#
# --mesh      remplace la ligne meshFile (quand le maillage vise n existe pas
#             encore, ex. la serie T3) — la surcharge est ECRITE dans le resume.
# --drop-key  retire une ligne `cle = ...` du deck de travail (pour prouver
#             qu un refus vient de CETTE cle et de rien d autre).
# Le deck original n est jamais modifie ; les runs vont dans <out>/<stem>_<tag>/.
# Le resume est ecrit dans <out>/smoke_<tag>.md (tableau) et imprime.
# Code de retour : 0 si tous les runs ont demarre ET fini (code 0), 1 sinon.
# ---------------------------------------------------------------------------
import argparse
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

RE_DT = re.compile(r"(\d+) tets, (\d+) joints, (\d+) nodes, dt = ([0-9.eE+-]+) s, steps = (\d+)")
RE_T = re.compile(r"^T\s*=\s*([0-9.eE+-]+)", re.M)
RE_MESH = re.compile(r"^meshFile\s*=\s*(\S+)", re.M)
RE_WALL = re.compile(r"wall time: ([0-9.]+) s")
RE_MASS = re.compile(r"corps '(\w+)': .*masse = ([0-9.eE+-]+) kg")
RE_ERR = re.compile(r"^\[rockim\] error.*$|^\s*- .*$", re.M)
RE_RES = re.compile(r"residu\s*: ([0-9.eE+-]+) J")


def rewrite(text, T, frames, mesh, drop):
    n_t = n_f = 0
    out = []
    for line in text.splitlines():
        s = line.strip()
        if re.match(r"^T\s*=", s):
            out.append("T = %g" % T); n_t += 1; continue
        if re.match(r"^frames\s*=", s):
            out.append("frames = %d" % frames); n_f += 1; continue
        if mesh and re.match(r"^meshFile\s*=", s):
            out.append("meshFile = %s" % mesh); continue
        if drop and re.match(r"^%s\s*=" % re.escape(drop), s):
            out.append("# [deck_smoke --drop-key] " + line); continue
        out.append(line)
    if n_t != 1 or n_f != 1:
        raise SystemExit("deck : %d ligne(s) T, %d ligne(s) frames (1 attendue)" % (n_t, n_f))
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("decks", nargs="+")
    ap.add_argument("--exe", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--T", type=float, default=2e-6)
    ap.add_argument("--frames", type=int, default=1)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--mesh", default=None)
    ap.add_argument("--drop-key", default=None)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--timeout", type=float, default=600.0)
    ap.add_argument("--ms-shared", type=float, default=77.73,
                    help="ms/pas machine partagee 14 fils (v3P : 5498 s / 70733 pas)")
    ap.add_argument("--ms-alone", type=float, default=17.0,
                    help="ms/pas machine seule (~20 min / 70733 pas, estimation du cadrage)")
    a = ap.parse_args()

    exe = os.path.abspath(os.path.join(ROOT, a.exe)) if not os.path.isabs(a.exe) else a.exe
    tag = a.tag or os.path.splitext(os.path.basename(a.exe))[0]
    out_root = os.path.abspath(os.path.join(ROOT, a.out))
    os.makedirs(out_root, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(a.threads))

    rows = []
    all_ok = True
    for deck in a.decks:
        src = os.path.abspath(os.path.join(ROOT, deck))
        stem = os.path.splitext(os.path.basename(src))[0]
        with open(src, encoding="utf-8") as f:
            text = f.read()
        m = RE_T.search(text)
        T_deck = float(m.group(1)) if m else float("nan")
        mm = RE_MESH.search(text)
        mesh_deck = mm.group(1) if mm else "?"
        work = os.path.join(out_root, "%s_%s" % (stem, tag))
        os.makedirs(work, exist_ok=True)
        cfg = os.path.join(work, stem + "_smoke.cfg")
        with open(cfg, "w", encoding="utf-8", newline="\n") as f:
            f.write(rewrite(text, a.T, a.frames, a.mesh, a.drop_key))
        run_dir = os.path.join(work, "run")
        log = os.path.join(work, stem + ".log")
        t0 = time.time()
        try:
            p = subprocess.run([exe, cfg, run_dir], cwd=ROOT, env=env,
                               capture_output=True, text=True, timeout=a.timeout,
                               encoding="utf-8", errors="replace")
            rc, so, se = p.returncode, p.stdout, p.stderr
        except subprocess.TimeoutExpired as e:
            rc, so, se = -1, (e.stdout or ""), (e.stderr or "")
            if isinstance(so, bytes): so = so.decode("utf-8", "replace")
            if isinstance(se, bytes): se = se.decode("utf-8", "replace")
        wall = time.time() - t0
        with open(log, "w", encoding="utf-8", newline="\n") as f:
            f.write(so)
            if se:
                f.write("\n---- stderr ----\n" + se)
        both = so + "\n" + se
        d = RE_DT.search(both)
        errs = [l.strip() for l in RE_ERR.findall(both) if l.strip()]
        masses = {k: float(v) for k, v in RE_MASS.findall(both)}
        res = RE_RES.search(both)
        nwarn = both.count("AVERTISSEMENT") + both.count("WARNING")
        if d:
            dt = float(d.group(4)); steps_2us = int(d.group(5))
            steps_full = int(round(T_deck / dt)) if T_deck == T_deck else -1
            cost_sh = steps_full * a.ms_shared / 1000.0
            cost_al = steps_full * a.ms_alone / 1000.0
        else:
            dt = float("nan"); steps_2us = steps_full = -1; cost_sh = cost_al = float("nan")
        started = d is not None and not errs
        ok = started and rc == 0
        all_ok = all_ok and ok
        rows.append(dict(deck=deck, rc=rc, started=started, ok=ok, dt=dt,
                         steps_2us=steps_2us, T_deck=T_deck, steps_full=steps_full,
                         cost_sh=cost_sh, cost_al=cost_al, wall=wall, nwarn=nwarn,
                         masses=masses, errs=errs, mesh=a.mesh or mesh_deck,
                         res=res.group(1) if res else "-"))
        print("%-45s rc=%2d started=%s dt=%.4g steps(%gus)=%d wall=%.0fs warn=%d" %
              (deck, rc, started, dt, a.T * 1e6, steps_2us, wall, nwarn))
        for e in errs[:8]:
            print("    " + e)

    md = []
    md.append("# Fumee `%s` — exe `%s`, T = %g s, frames = %d, OMP_NUM_THREADS = %d%s%s\n" % (
        tag, a.exe, a.T, a.frames, a.threads,
        (", meshFile surcharge -> `%s`" % a.mesh) if a.mesh else "",
        (", cle retiree `%s`" % a.drop_key) if a.drop_key else ""))
    md.append("Cout estime du deck complet = pas(T_deck/dt) x ms/pas ; %.2f ms/pas partage "
              "(results/yang_bench_s25_v3P.log : 5 498 s / 70 733 pas, 14 fils) ; "
              "%.1f ms/pas seul (estimation du cadrage, non mesuree).\n" % (a.ms_shared, a.ms_alone))
    md.append("| deck | maillage | demarre | code | dt (ns) | pas (%g us) | T deck (us) | pas total | "
              "cout partage | cout seul | mur fumee | avert. | residu B4 (J) | masses piston / bit / insert (kg) |" % (a.T * 1e6))
    md.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        ms = r["masses"]
        mtxt = " / ".join("%.4g" % ms[k] if k in ms else "-" for k in ("piston", "bit", "insert"))
        md.append("| `%s` | `%s` | %s | %d | %s | %s | %s | %s | %s | %s | %.0f s | %d | %s | %s |" % (
            os.path.basename(r["deck"]), os.path.basename(r["mesh"]),
            "oui" if r["started"] else "NON", r["rc"],
            ("%.3f" % (r["dt"] * 1e9)) if r["dt"] == r["dt"] else "-",
            r["steps_2us"] if r["steps_2us"] >= 0 else "-",
            ("%g" % (r["T_deck"] * 1e6)) if r["T_deck"] == r["T_deck"] else "-",
            ("%d" % r["steps_full"]) if r["steps_full"] >= 0 else "-",
            ("%.0f s (%.1f h)" % (r["cost_sh"], r["cost_sh"] / 3600)) if r["cost_sh"] == r["cost_sh"] else "-",
            ("%.0f s (%.0f min)" % (r["cost_al"], r["cost_al"] / 60)) if r["cost_al"] == r["cost_al"] else "-",
            r["wall"], r["nwarn"], r["res"], mtxt))
    for r in rows:
        if r["errs"]:
            md.append("\nRefus / erreurs de `%s` :\n" % os.path.basename(r["deck"]))
            md += ["    " + e for e in r["errs"]]
    summary = os.path.join(out_root, "smoke_%s.md" % tag)
    with open(summary, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(md) + "\n")
    print("\n".join(md))
    print("\nresume : %s" % summary)
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
