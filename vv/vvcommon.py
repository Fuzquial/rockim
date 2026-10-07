#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# vvcommon.py : briques communes des bancs de la campagne V&V (vv/<banc>/).
#
# Chaque banc est un dossier vv/<ID>_<nom>/ avec un script <id>.py qui expose :
#   JOBS(args) -> liste de Job (un run rockim chacun, deck deja ecrit)
#   analyse()  -> ecrit resultats.json et les figures, imprime le verdict
# et une fonction main() : `python3 vv/<banc>/<id>.py {prepare,run,analyse,all}`.
# La file commune `python3 vv/run_queue.py <bancs...> --slots 4 --threads 2`
# execute les Job de plusieurs bancs en parallele (un processus par slot).
#
# Regles de la campagne (docs/VV_campagne.md) :
#   - uniquement des cles rockim existantes, aucune modification de src/ ;
#   - criteres d'acceptation ecrits dans le script AVANT le premier calcul ;
#   - un run est « fini » si son journal contient « wall time » ;
#   - sorties brutes dans <banc>/out/ (ignore par git), maillages dans
#     <banc>/meshes/ (*.msh ignores), figures PDF + PNG dans <banc>/.
# ---------------------------------------------------------------------------
import os, re, subprocess, time
from dataclasses import dataclass, field

VV = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(VV)
EXE_DEFAULT = os.path.join(ROOT, "build_vv", "rockim")   # HEAD b80a4ac ; bit-identique a build_nofma sur 2 decks V1


@dataclass
class Job:
    bench: str            # identifiant du banc, ex. "P2.1"
    name: str             # identifiant du run, ex. "potential_v1"
    cfg: str              # chemin du deck
    outdir: str           # dossier de sortie (outputDir du deck)
    weight: float = 1.0   # cout relatif estime (pour ordonner la file)
    meta: dict = field(default_factory=dict)


def write_deck(outdir, lines, header=None):
    """ecrit outdir/deck.cfg (outputDir ajoute en fin) et renvoie son chemin"""
    os.makedirs(outdir, exist_ok=True)
    body = ([f"# {header}"] if header else []) + list(lines) + [f"outputDir = {outdir}"]
    cfg = os.path.join(outdir, "deck.cfg")
    with open(cfg, "w") as f:
        f.write("\n".join(body) + "\n")
    return cfg


def is_done(outdir):
    lg = os.path.join(outdir, "rockim.log")
    return os.path.exists(lg) and "wall time" in open(lg, errors="replace").read()


def running_rockim():
    """nombre de calculs rockim en cours sur la machine (toutes files confondues)"""
    out = subprocess.run(["pgrep", "-f", r"/rockim .*\.cfg"], capture_output=True, text=True).stdout
    return len([l for l in out.split() if l.strip()])


def wait_budget(max_procs):
    """plafond GLOBAL de calculs simultanes, partage par toutes les files
    (VV_MAX_PROCS, defaut 5 : 5 x 2 fils = 10 coeurs)"""
    import random
    while running_rockim() >= max_procs:
        time.sleep(20 + 10 * random.random())


def run_job(job, exe=EXE_DEFAULT, threads=None, force=False):
    """lance rockim sur job.cfg ; journal dans outdir/rockim.log ; saute un run fini"""
    if is_done(job.outdir) and not force:
        return 0, 0.0
    wait_budget(int(os.environ.get("VV_MAX_PROCS", "5")))
    env = dict(os.environ)
    if threads:
        env["OMP_NUM_THREADS"] = str(threads)
    t0 = time.time()
    with open(os.path.join(job.outdir, "rockim.log"), "w") as log:
        rc = subprocess.run([exe, job.cfg], stdout=log, stderr=subprocess.STDOUT,
                            cwd=ROOT, env=env).returncode
    return rc, time.time() - t0


def parse_log(outdir):
    """grandeurs communes du journal fdem3d / fem3d / fdem"""
    txt = open(os.path.join(outdir, "rockim.log"), errors="replace").read()

    def f(pat, cast=float):
        m = re.search(pat, txt)
        return cast(m.group(1).rstrip(",")) if m else None

    return dict(
        done="wall time" in txt,
        wall=f(r"wall time: (\S+) s"),
        dt=f(r"dt = (\S+) s"),
        steps=f(r"steps = (\d+)", int),
        ntet=f(r"(\d+) tets", int),
        budget_pct=f(r"residu\s*:\s*\S+ J \((\S+) % de l'echelle\)"),
        broken=f(r"broken joints\s*:\s*(\d+)", int),
        text=txt,
    )


def read_history(outdir, name="history.csv"):
    import numpy as np
    return np.genfromtxt(os.path.join(outdir, name), delimiter=",", names=True)


def plot_style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIX Two Text"],
                         "mathtext.fontset": "stix", "font.size": 9})
    return plt


def savefig(fig, path_noext):
    fig.savefig(path_noext + ".pdf")
    fig.savefig(path_noext + ".png", dpi=200)


def verdict_line(label, value, crit, fmt="{:+.3%}", absval=True):
    ok = (abs(value) if absval else value) <= crit
    return ok, f"  {label:28s} {fmt.format(value):>12s}  critere {crit:g}  {'PASSE' if ok else 'ECHEC'}"
