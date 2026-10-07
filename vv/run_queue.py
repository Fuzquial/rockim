#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# run_queue.py : file d'execution commune des bancs V&V.
#
#   python3 vv/run_queue.py P3_fissure_guo P4_essais_guo --slots 4 --threads 2
#   python3 vv/run_queue.py --list
#
# Importe vv/<banc>/<script>.py (le .py du dossier qui definit JOBS), recupere
# ses Job, les trie du plus lourd au plus leger et les execute sur `slots`
# processus paralleles a `threads` fils OpenMP chacun. Un run deja fini (journal
# avec « wall time ») est saute. Un recapitulatif est ecrit dans
# vv/run_queue.log. Les analyses restent a lancer par banc (`<script> analyse`).
# Regle : slots x threads <= nombre de coeurs (10 sur le Mac de reference).
# ---------------------------------------------------------------------------
import argparse, importlib.util, os, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed

VV = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, VV)
import vvcommon as C  # noqa: E402


def load_bench(folder):
    d = os.path.join(VV, folder)
    for fn in sorted(os.listdir(d)):
        if fn.endswith(".py") and not fn.startswith(("viz", "_")):
            spec = importlib.util.spec_from_file_location(folder + "." + fn[:-3], os.path.join(d, fn))
            mod = importlib.util.module_from_spec(spec)
            sys.path.insert(0, d)
            spec.loader.exec_module(mod)
            if hasattr(mod, "JOBS"):
                return mod
    raise SystemExit(f"{folder} : aucun script avec JOBS()")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("benches", nargs="*")
    ap.add_argument("--slots", type=int, default=4)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--exe", default=C.EXE_DEFAULT)
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list:
        for d in sorted(os.listdir(VV)):
            if os.path.isdir(os.path.join(VV, d)) and not d.startswith(("_", ".")):
                print(d)
        return
    jobs = []
    for b in a.benches:
        jobs += load_bench(b).JOBS(None)
    jobs.sort(key=lambda j: -j.weight)
    todo = [j for j in jobs if not C.is_done(j.outdir)]
    print(f"[file] {len(jobs)} runs, {len(todo)} a faire, {a.slots} slots x {a.threads} fils", flush=True)
    log = open(os.path.join(VV, "run_queue.log"), "a")
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=a.slots) as ex:
        fut = {ex.submit(C.run_job, j, os.path.abspath(a.exe), a.threads): j for j in todo}
        for f in as_completed(fut):
            j = fut[f]
            rc, wall = f.result()
            line = f"{time.strftime('%H:%M:%S')}  {j.bench:6s} {j.name:30s} rc = {rc}  {wall:8.1f} s"
            print(line, flush=True)
            log.write(line + "\n"); log.flush()
    print(f"[file] termine en {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
