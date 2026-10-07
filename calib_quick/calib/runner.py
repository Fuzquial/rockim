#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# runner.py — file d attente de runs rockim avec BUDGET DE THREADS global.
#
#   python calib_quick/calib/runner.py jobs.json --exe rockim_f2l.exe --parallel 3 --threads 4
#
# jobs.json : [{"cfg": "calib_quick/xxx.cfg", "out": "out_xxx"}, ...]
# Regle maison : un seul "job 14 cpus" a la fois -> parallel x threads <= 14.
# Chaque run ecrit out/_log.txt et out/_run.json (mur, code retour) ; un run
# dont out/history.csv existe deja et dont _run.json dit "ok" est SAUTE
# (reprise apres interruption). Ne lance JAMAIS sans validation de la liste.
# ---------------------------------------------------------------------------
import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def run_one(job, exe, threads):
    cfg, out = job["cfg"], job["out"]
    stamp = os.path.join(out, "_run.json")
    if os.path.exists(stamp):
        try:
            if json.load(open(stamp)).get("rc") == 0 and os.path.exists(os.path.join(out, "history.csv")):
                return dict(job, skipped=True, rc=0, wall=0.0)
        except Exception:
            pass
    os.makedirs(out, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    t0 = time.time()
    with open(os.path.join(out, "_log.txt"), "w", encoding="utf-8", errors="replace") as log:
        p = subprocess.run([exe, cfg, out], cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    res = dict(job, rc=p.returncode, wall=time.time() - t0, threads=threads, exe=exe,
               finished=time.strftime("%Y-%m-%d %H:%M:%S"))
    json.dump(res, open(stamp, "w"), indent=1)
    return res


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs")
    ap.add_argument("--exe", default=os.path.join(ROOT, "rockim_f2l.exe"))
    ap.add_argument("--parallel", type=int, default=3)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--budget", type=int, default=14, help="threads totaux autorises")
    a = ap.parse_args()
    if a.parallel * a.threads > a.budget:
        raise SystemExit(f"parallel x threads = {a.parallel * a.threads} > budget {a.budget}")
    if not os.path.isabs(a.exe):                 # subprocess resout l exe depuis le cwd PARENT, pas depuis cwd=ROOT
        a.exe = os.path.join(ROOT, a.exe)
    if not os.path.exists(a.exe):
        raise SystemExit(f"executable introuvable : {a.exe}")
    jobs = json.load(open(a.jobs))
    print(f"[runner] {len(jobs)} runs, {a.parallel} en parallele x {a.threads} threads, exe {a.exe}")
    t0 = time.time()
    done = 0
    with ThreadPoolExecutor(max_workers=a.parallel) as ex:
        futs = {ex.submit(run_one, j, a.exe, a.threads): j for j in jobs}
        for f in as_completed(futs):
            r = f.result(); done += 1
            tag = "saute" if r.get("skipped") else ("ok" if r["rc"] == 0 else f"RC={r['rc']}")
            print(f"[runner] {done:3d}/{len(jobs)}  {r['out']:36s} {tag:6s} {r['wall']:6.0f} s   "
                  f"(total {time.time() - t0:6.0f} s)", flush=True)
    print(f"[runner] fini en {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
