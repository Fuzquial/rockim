#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# run_daemon.py : lanceur de fond de la campagne V&V.
#
#   nohup python3 vv/run_daemon.py --slots 4 --threads 2 &
#
# Toutes les `poll` secondes, cherche les bancs vv/<banc>/ qui contiennent un
# fichier PRET (pose apres relecture du banc) et pas encore FINI. Pour chacun,
# dans l'ordre de pose des PRET : execute ses Job par la file commune
# (vv/run_queue.py : slots processus x threads fils, runs finis sautes), puis
# lance `<script> analyse`, puis pose FINI (ou ECHEC avec le code de retour).
# Le fichier PRET peut contenir des options, une par ligne :
#   slots = 2
#   threads = 4
# Les bancs prets partent EN PARALLELE, chacun dans sa file ; le plafond
# global de calculs simultanes est VV_MAX_PROCS (vvcommon.wait_budget).
# Etats : PRET -> ENCOURS -> FINI | ECHEC.
# Arret propre : creer vv/STOP (les files en cours se terminent d'abord).
# Journal : vv/run_daemon.log.
# ---------------------------------------------------------------------------
import argparse, os, subprocess, sys, time

VV = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(VV)


def log(msg):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {msg}"
    print(line, flush=True)
    with open(os.path.join(VV, "run_daemon.log"), "a") as f:
        f.write(line + "\n")


def ready_benches():
    out = []
    for d in os.listdir(VV):
        p = os.path.join(VV, d)
        if (os.path.isdir(p) and os.path.exists(os.path.join(p, "PRET"))
                and not os.path.exists(os.path.join(p, "FINI"))
                and not os.path.exists(os.path.join(p, "ECHEC"))):
            out.append((os.path.getmtime(os.path.join(p, "PRET")), d))
    return [d for _, d in sorted(out)]


def options(bench, slots, threads):
    for line in open(os.path.join(VV, bench, "PRET")):
        if "=" in line:
            k, v = (x.strip() for x in line.split("=", 1))
            if k == "slots":
                slots = int(v)
            elif k == "threads":
                threads = int(v)
    return slots, threads


def script_of(bench):
    d = os.path.join(VV, bench)
    for fn in sorted(os.listdir(d)):
        if fn.endswith(".py") and not fn.startswith(("viz", "_")) and "JOBS" in open(os.path.join(d, fn)).read():
            return os.path.join(d, fn)
    return None


def finish(b, script, rc):
    rc2 = subprocess.run([sys.executable, script, "analyse"], cwd=ROOT,
                         stdout=open(os.path.join(VV, b, "analyse.log"), "w"),
                         stderr=subprocess.STDOUT).returncode
    mark = "FINI" if rc == 0 and rc2 == 0 else "ECHEC"
    enc = os.path.join(VV, b, "ENCOURS")
    if os.path.exists(enc):
        os.remove(enc)
    open(os.path.join(VV, b, mark), "w").write(f"runs rc {rc}, analyse rc {rc2}, {time.ctime()}\n")
    log(f"{b} : {mark} (analyse rc = {rc2}, voir {b}/analyse.log)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--slots", type=int, default=4)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--poll", type=float, default=60.0)
    a = ap.parse_args()
    log(f"demarrage (bancs en parallele, plafond global VV_MAX_PROCS = "
        f"{os.environ.get('VV_MAX_PROCS', '5')} calculs), scrutation {a.poll:g} s")
    running = {}                                   # banc -> (Popen, script, t0)
    while True:
        for b, (pr, script, t0) in list(running.items()):
            rc = pr.poll()
            if rc is not None:
                log(f"{b} : runs termines, rc = {rc}, {time.time() - t0:.0f} s ; analyse")
                finish(b, script, rc)
                del running[b]
        if os.path.exists(os.path.join(VV, "STOP")):
            if not running:
                break
            time.sleep(a.poll)
            continue
        for b in ready_benches():
            if b in running or os.path.exists(os.path.join(VV, b, "ENCOURS")):
                continue
            slots, threads = options(b, a.slots, a.threads)
            script = script_of(b)
            if script is None:
                log(f"{b} : aucun script avec JOBS, banc ignore")
                open(os.path.join(VV, b, "ECHEC"), "w").write("pas de JOBS\n")
                continue
            open(os.path.join(VV, b, "ENCOURS"), "w").write(time.ctime() + "\n")
            pr = subprocess.Popen([sys.executable, os.path.join(VV, "run_queue.py"), b,
                                   "--slots", str(slots), "--threads", str(threads)],
                                  cwd=ROOT, stdout=open(os.path.join(VV, b, "queue.log"), "a"),
                                  stderr=subprocess.STDOUT)
            running[b] = (pr, script, time.time())
            log(f"{b} : file lancee ({slots} slots x {threads} fils, sous le plafond global)")
        time.sleep(a.poll)
    log("arret demande (vv/STOP)")


if __name__ == "__main__":
    main()
