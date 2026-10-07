#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# bitid.py — preuve de BIT-IDENTITE generique de rockim (mesure B6 du plan de
# robustesse du 2026-09-05). Remplace les scripts bitid_wNN.sh copies-colles
# par build (exe, decks et nombre de fils codes en dur, ancre changee deux
# fois, hachages jamais compares au meme nombre de fils).
#
#   python tools/bitid.py --exe rockim_f2w18.exe            # compare a l'ancre
#   python tools/bitid.py --exe rockim_f2w18.exe --update   # (re)prend l'ancre
#   python tools/bitid.py --exe build/rockim.exe --only fem3d
#   python tools/bitid.py --list
#
# Ce que le script fait, pour UN exe :
#   1. joue la liste FIXE `DECKS` (8 decks courts de tests_f2/bitid/ : 3 fem3d,
#      3 fdem3d, 2 fdem), chacun < 5 min a OMP_NUM_THREADS = 4, l'exe lance
#      depuis la RACINE rockim_f2 (les meshFile des decks y sont relatifs) ;
#   2. hache (SHA-256 du contenu) history.csv, frames.csv et le(s) .vtu de la
#      DERNIERE trame (fdem : fdem_NNNN.vtu et fdem_joints_NNNN.vtu), plus
#      probes.csv s'il existe ; releve le pic de force et le temps de calcul ;
#   3. compare a tools/bitid_refs.json (l'ANCRE, au depot) au MEME nombre de
#      fils que l'ancre — les reductions OpenMP rendent les hachages
#      incomparables entre nombres de fils differents — et imprime
#      « N/N IDENTIQUE » ou la liste des differences avec les pics de force ;
#   4. --update reecrit l'ancre (exe_sha256, date, fils, hachages). Regle :
#      toute ancre changee = une ligne dans CHANGELOG.md (« ancre changee :
#      raison ») et le json committe dans le meme commit. Voir tools/BITID.md.
#
# Boucle de reprise Apex One integree (l'antivirus refuse par intermittence
# un exe frais : WinError 5 / rc 126) ; les .vtu sont supprimes apres hachage
# (--keep pour les garder). Les sorties (history.csv, frames.csv, journal
# <deck>.log) vont dans un dossier TEMPORAIRE (tempfile.mkdtemp, meme
# convention que verify_suite.py ; le chemin est imprime) sauf --outroot ;
# le rapport durable est --json. Codes de sortie : 0 identique, 1 difference
# ou deck sans reference, 2 echec de run / exe introuvable.
# ---------------------------------------------------------------------------
import argparse, hashlib, json, os, platform, re, shutil, socket, subprocess, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DECKDIR = "tests_f2/bitid"   # separateur posix : la meme chaine sur Windows et sur la CI Linux (B8)
REFS_DEFAULT = os.path.join(HERE, "bitid_refs.json")
THREADS_DEFAULT = 4          # l'ancre est prise a 4 fils (14 + 4 = 18 coeurs logiques)

# ---- liste FIXE des decks (nom, deck relatif a la racine, mode, loi, note) ----
# `threads` : None = valeur commune (--threads, defaut 4) ; un entier force ce
# deck a ce nombre de fils (a n'utiliser que pour un deck NON deterministe a
# 4 fils, avec la raison dans `note`).
DECKS = [
    dict(name="fem3d_cdp_PQ_court", cfg=f"{DECKDIR}/fem3d_cdp_PQ_court.cfg",
         mode="fem3d", law="cdp", threads=None,
         note="quart de bloc Q1_c05 (13 687 tets), cdp inverse, penalite, confinement 20 MPa, 130 us"),
    dict(name="fem3d_dpr_T1_court", cfg=f"{DECKDIR}/fem3d_dpr_T1_court.cfg",
         mode="fem3d", law="dpr", threads=None,
         note="T1_c05_clean (98 342 tets), dpr + briques (apex, rankine stress, power, crackband) + signorini + bulkViscosity + hencky"),
    dict(name="fem3d_sk2011_cyl_court", cfg=f"{DECKDIR}/fem3d_sk2011_cyl_court.cfg",
         mode="fem3d", law="saksala2011", threads=None,
         note="cylindre interne Kuhn + jitter, saksala2011, penalite"),
    dict(name="fdem3d_kuru9_court", cfg=f"{DECKDIR}/fdem3d_kuru9_court.cfg",
         mode="fdem3d", law="(joints cohesifs + bulkDamage yang)", threads=None,
         note="impact insert unique Kuru (6 corps, phases, groupBond), adaptive + potential + DIF yang-fig2"),
    dict(name="fdem3d_cut3d_heilman_court", cfg=f"{DECKDIR}/fdem3d_cut3d_heilman_court.cfg",
         mode="fdem3d", law="(joints cohesifs, volume elastique)", threads=None,
         note="coupe PDC 3D Heilman, toolShape = pdc, signorini, adaptive + potential"),
    dict(name="fdem3d_visc_yan_3d", cfg=f"{DECKDIR}/fdem3d_visc_yan_3d.cfg",
         mode="fdem3d", law="(joints cohesifs + viscosite Yan)", threads=None,
         note="= test visc_yan_3d de verify_suite (cube 6 000 tets, bulkViscosityXi = 2)"),
    dict(name="fdem_toolcontact_signorini", cfg=f"{DECKDIR}/fdem_toolcontact_signorini.cfg",
         mode="fdem", law="(joints cohesifs, volume elastique)", threads=None,
         note="= test t1_toolcontact_signorini de verify_suite (banc T1, cutter PDC 2D)"),
    dict(name="fdem_ucs_yan_adaptive_court", cfg=f"{DECKDIR}/fdem_ucs_yan_adaptive_court.cfg",
         mode="fdem", law="(joints cohesifs yan, insertion adaptative)", threads=None,
         note="UCS de Yan 2023 sur Voronoi/Delaunay, raccourci"),
    # 9e deck (conseil du 12/09, N1) : le chemin de la loi de joint du deck v2
    # Yang — intrinseque + munjiza + parabolic + yang + coulomb + origin +
    # midedge + majority + deltaC + DIF continu — qu aucun des 8 premiers
    # n exerce. Reference LATERALE tools/bitid_refs_jointlaw.json (--refs),
    # pour ne pas rouvrir l ancre principale a chaque lot de loi de joint :
    #   python tools/bitid.py --exe X --only yang_v2 --refs tools/bitid_refs_jointlaw.json
    dict(name="fdem3d_yang_v2_court", cfg=f"{DECKDIR}/fdem3d_yang_v2_court.cfg",
         mode="fdem3d", law="(joints cohesifs, conventions Solidity, DIF continu)", threads=None, side=True,
         note="banc Yang s = 2,5 (10 563 tets), bit lance a 9 m/s, 20 us : la loi de joint v2 sous fracture ; ancre laterale bitid_refs_jointlaw.json"),
]

RX_PEAK = re.compile(r"peak tool force\s*:\s*(-?[\d.eE+]+)")
RX_WALL = re.compile(r"wall time:\s*([\d.eE+]+)\s*s")
RX_DT   = re.compile(r"\bdt = ([\d.eE+-]+) s")
RX_STEP = re.compile(r"steps = (\d+)")
RX_VTU  = re.compile(r"^(.*?)_(\d{4})\.vtu$")
APEX_RC = {126, 5, -1073741790, 3221225506}     # rc 126 (bash), WinError 5, 0xC0000022
# message d'acces refuse (anglais, francais « Accès refusé » / « L'accès est refusé »,
# accents eventuellement remplaces au decodage) — un motif ETROIT : « Acc » seul
# aurait pris `fragBrushAccel` ou tout « Accept » d'un vrai echec pour Apex One
RX_APEX = re.compile(r"Permission denied|Access is denied|Acc\S{0,3}s\s+(?:est\s+)?refus", re.I)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_exe(exe):
    for cand in (exe, os.path.join(ROOT, exe), os.path.join(os.getcwd(), exe)):
        if os.path.isfile(cand):
            return os.path.abspath(cand)
    return None


def last_frame_vtus(outdir):
    """Les .vtu de la DERNIERE trame, par famille (fdem3d_, fdem3d_joints_...)."""
    fam = {}
    for f in os.listdir(outdir):
        m = RX_VTU.match(f)
        if m:
            fam.setdefault(m.group(1), []).append((int(m.group(2)), f))
    return sorted(f for lst in fam.values() for (_, f) in [max(lst)])


def drop_vtus(out):
    """Supprime les .vtu d'un dossier de sortie (gros, deja haches ou inutiles)."""
    if os.path.isdir(out):
        for f in os.listdir(out):
            if f.endswith(".vtu"):
                os.remove(os.path.join(out, f))


def run_deck(exe, deck, threads, outroot, timeout, apex_retries, keep):
    cfg = os.path.join(ROOT, deck["cfg"])
    if not os.path.isfile(cfg):
        return dict(rc=None, error=f"deck introuvable : {cfg}")
    out = os.path.join(outroot, deck["name"])
    log = out + ".log"
    shutil.rmtree(out, ignore_errors=True)
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    tries, rc, text, t0 = 0, None, "", time.time()
    while True:
        tries += 1
        try:
            p = subprocess.run([exe, cfg, out], cwd=ROOT, env=env, timeout=timeout,
                               capture_output=True, text=True, errors="replace")
            rc, text = p.returncode, p.stdout + p.stderr
        except PermissionError as e:                 # Apex One : exe frais refuse
            rc, text = 5, f"PermissionError: {e}"
        except subprocess.TimeoutExpired as e:     # l'enfant est tue ; sortie partielle gardee
            part = (e.stdout or b"") if isinstance(e.stdout, (bytes, bytearray)) else (e.stdout or "")
            if isinstance(part, (bytes, bytearray)):
                part = part.decode("utf-8", "replace")
            with open(log, "w", encoding="utf-8", errors="replace") as f:
                f.write(part + f"\n[bitid] TIMEOUT {timeout} s (deck tue)\n")
            if not keep:
                drop_vtus(out)
            return dict(rc=124, error=f"TIMEOUT {timeout} s", tries=tries)
        apex = rc in APEX_RC or (rc != 0 and RX_APEX.search(text) is not None
                                 and "history.csv" not in text)
        if not apex or tries > apex_retries:
            break
        time.sleep(20)
    wall_total = time.time() - t0
    with open(log, "w", encoding="utf-8", errors="replace") as f:
        f.write(text)
    res = dict(rc=rc, tries=tries, apex_retries=tries - 1, wall_total_s=round(wall_total, 1))
    if rc != 0:
        tail = text.strip().splitlines()[-1] if text.strip() else "?"
        res["error"] = f"rc={rc} : {tail}"
        if not keep:
            drop_vtus(out)
        return res
    m = RX_PEAK.search(text); res["peak_force_N"] = float(m.group(1)) if m else None
    m = RX_WALL.search(text); res["wall_s"] = float(m.group(1)) if m else None
    m = RX_DT.search(text);   res["dt_s"] = float(m.group(1)) if m else None
    m = RX_STEP.search(text); res["steps"] = int(m.group(1)) if m else None
    hashes = {}
    for fn in ("history.csv", "frames.csv", "probes.csv"):
        p = os.path.join(out, fn)
        if os.path.isfile(p):
            hashes[fn] = sha256_file(p)
    for fn in last_frame_vtus(out):
        hashes[fn] = sha256_file(os.path.join(out, fn))
    res["hashes"] = hashes
    if not keep:
        drop_vtus(out)
    return res


def fmt_peak(v):
    return "pic n/a" if v is None else f"pic {v:.6g} N"


def main():
    ap = argparse.ArgumentParser(description="bit-identite rockim : joue les 8 decks de "
                                 "tests_f2/bitid/, hache les sorties, compare a l'ancre")
    ap.add_argument("--exe", help="binaire (chemin, ou nom dans la racine rockim_f2)")
    ap.add_argument("--refs", default=REFS_DEFAULT, help="ancre json (defaut tools/bitid_refs.json)")
    ap.add_argument("--update", action="store_true",
                    help="reecrit l'ancre avec les hachages de CE run (+ CHANGELOG + commit !)")
    ap.add_argument("--threads", type=int, default=None,
                    help="OMP_NUM_THREADS ; defaut = celui de l'ancre par deck (4 si pas d'ancre). "
                         "Un nombre different de l'ancre rend la comparaison NON COMPARABLE")
    ap.add_argument("--only", default=None, help="filtre sur le nom du deck (sous-chaine)")
    ap.add_argument("--timeout", type=float, default=900.0, help="par deck, s (defaut 900)")
    ap.add_argument("--apex-retries", type=int, default=10, dest="apex_retries")
    ap.add_argument("--keep", action="store_true", help="garder les .vtu apres hachage")
    ap.add_argument("--outroot", default=None,
                    help="dossier des sorties (defaut : dossier TEMPORAIRE rockim_bitid_<exe>_omp<N>_* "
                         "cree par tempfile.mkdtemp, comme verify_suite.py ; le chemin est imprime)")
    ap.add_argument("--json", default=None, help="rapport json de ce run")
    ap.add_argument("--list", action="store_true", help="imprime la liste des decks et sort")
    args = ap.parse_args()

    try:   # utf-8 et ligne par ligne (un log redirige se lit pendant le run)
        sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)
    except Exception:
        pass

    if args.list:
        for d in DECKS:
            print(f"{d['name']:<30s} mode={d['mode']:<6s} law={d['law']:<45s} "
                  f"threads={d['threads'] or 'commun'}  {d['cfg']}")
        return 0
    if not args.exe:
        ap.error("--exe requis (ou --list)")
    exe = resolve_exe(args.exe)
    if not exe:
        print(f"[bitid] exe introuvable : {args.exe}")
        return 2
    exe_sha = sha256_file(exe)
    exe_stem = os.path.splitext(os.path.basename(exe))[0]

    refs = {}
    if os.path.isfile(args.refs):
        with open(args.refs, encoding="utf-8") as f:
            refs = json.load(f)
    ref_decks = refs.get("decks", {})
    meta = refs.get("_meta", {})

    # Les decks `side` (ancre LATERALE, ex. fdem3d_yang_v2_court) ne sont joues
    # que sur --only ou sur une ancre --refs autre que la principale : la passe
    # 8/8 par defaut reste ce qu elle etait (ancre principale intacte).
    sel = [d for d in DECKS
           if (args.only is None or args.only in d["name"])
           and (not d.get("side") or args.only is not None
                or os.path.abspath(args.refs) != os.path.abspath(REFS_DEFAULT))]
    if not sel:
        print(f"[bitid] aucun deck ne contient '{args.only}'")
        return 2

    # nombre de fils par deck : --threads > ancre > liste > 4
    def threads_for(d):
        if args.threads:
            return args.threads
        r = ref_decks.get(d["name"])
        if r and r.get("threads"):
            return int(r["threads"])
        return d["threads"] or THREADS_DEFAULT

    thr_set = sorted({threads_for(d) for d in sel})
    thr_tag = "omp" + "-".join(map(str, thr_set))
    if args.outroot:
        outroot = os.path.abspath(args.outroot)
        os.makedirs(outroot, exist_ok=True)
    else:   # dossier temporaire (rien n'est ecrit dans le depot ; le rapport durable = --json)
        outroot = tempfile.mkdtemp(prefix=f"rockim_bitid_{exe_stem}_{thr_tag}_")

    print(f"[bitid] exe = {exe}")
    print(f"[bitid] exe_sha256 = {exe_sha}")
    if meta:
        same = "MEME exe que l'ancre (rejeu)" if meta.get("exe_sha256") == exe_sha else \
               f"ancre = {meta.get('exe')} ({str(meta.get('exe_sha256', ''))[:16]}..., {meta.get('date')})"
        print(f"[bitid] ancre {args.refs} : {len(ref_decks)} decks, threads {meta.get('threads')} — {same}")
    elif not args.update:
        print(f"[bitid] AUCUNE ancre ({args.refs} absent) : rien a comparer, utiliser --update")
    print(f"[bitid] {len(sel)} decks, OMP_NUM_THREADS = {thr_set}, sorties {outroot}")

    results, n_id, n_diff, n_noref, n_fail = {}, 0, 0, 0, 0
    for d in sel:
        thr = threads_for(d)
        t0 = time.time()
        r = run_deck(exe, d, thr, outroot, args.timeout, args.apex_retries, args.keep)
        r.update(cfg=d["cfg"], mode=d["mode"], law=d["law"], threads=thr)
        results[d["name"]] = r
        if r.get("rc") != 0:
            n_fail += 1
            print(f"  ECHEC      {d['name']:<30s} {r.get('error')} (essais {r.get('tries')})")
            continue
        hh = r["hashes"]
        short = " ".join(f"{k}={v[:16]}" for k, v in hh.items())
        ref = ref_decks.get(d["name"])
        wall = r.get("wall_s") or (time.time() - t0)
        if args.update or ref is None:
            tag = "ANCRE     " if args.update else "SANS REF  "
            if not args.update:
                n_noref += 1
            print(f"  {tag} {d['name']:<30s} {wall:7.1f} s  omp{thr}  {fmt_peak(r.get('peak_force_N'))}  {short}")
            continue
        if int(ref.get("threads", 0)) != thr:
            n_diff += 1
            print(f"  NON COMPARABLE {d['name']:<26s} threads ancre {ref.get('threads')} != {thr}  {short}")
            r["verdict"] = "non comparable (threads)"
            continue
        diffs = []
        for k in sorted(set(hh) | set(ref.get("hashes", {}))):
            a, b = ref.get("hashes", {}).get(k), hh.get(k)
            if a != b:
                diffs.append(f"{k} ancre {str(a)[:16] if a else 'absent'} != {str(b)[:16] if b else 'absent'}")
        if diffs:
            n_diff += 1
            r["verdict"] = "DIFFERENT"
            print(f"  DIFFERENT  {d['name']:<30s} {wall:7.1f} s  omp{thr}  "
                  f"pic ancre {ref.get('peak_force_N')} -> {r.get('peak_force_N')} N ; " + " ; ".join(diffs))
        else:
            n_id += 1
            r["verdict"] = "IDENTIQUE"
            print(f"  IDENTIQUE  {d['name']:<30s} {wall:7.1f} s  omp{thr}  {fmt_peak(r.get('peak_force_N'))}")

    n = len(sel)
    if args.update:
        if n_fail:
            print(f"[bitid] {n_fail} echec(s) : ancre NON ecrite")
            return 2
        now = time.strftime("%Y-%m-%d %H:%M:%S")
        new = dict(refs) if refs else {}
        full_meta = dict(exe=os.path.basename(exe), exe_path=exe, exe_sha256=exe_sha,
                         exe_size=os.path.getsize(exe),
                         exe_mtime=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(exe))),
                         date=now, threads=thr_set[0] if len(thr_set) == 1 else thr_set,
                         host=socket.gethostname(), platform=platform.platform(), python=platform.python_version(),
                         hash="sha256 du contenu de history.csv, frames.csv, probes.csv (si present) et des .vtu de la "
                              "derniere trame ; comparable au MEME nombre de fils seulement",
                         regle="toute ancre changee = ligne CHANGELOG.md « ancre bitid changee : raison » + commit du json")
        if args.only and new.get("_meta"):
            # mise a jour PARTIELLE : l'ancre garde son _meta, la retouche est journalisee ;
            # un exe different de celui de l'ancre donne une ancre MIXTE (chaque deck porte son exe_sha256)
            m = new["_meta"]
            m.setdefault("partial_updates", []).append(dict(date=now, decks=sorted(results), exe=os.path.basename(exe),
                                                            exe_sha256=exe_sha, only=args.only))
            if m.get("exe_sha256") != exe_sha:
                m["MIXTE"] = "des decks proviennent d'un exe different de _meta.exe_sha256 (voir decks[*].exe_sha256)"
                print(f"[bitid] ATTENTION : ancre MIXTE — {os.path.basename(exe)} != exe de l'ancre {m.get('exe')}")
        else:
            new["_meta"] = full_meta
        decks = new.setdefault("decks", {})
        for name, r in results.items():
            decks[name] = dict(cfg=r["cfg"], mode=r["mode"], law=r["law"], threads=r["threads"],
                               hashes=r["hashes"], peak_force_N=r.get("peak_force_N"), wall_s=r.get("wall_s"),
                               dt_s=r.get("dt_s"), steps=r.get("steps"), apex_retries=r.get("apex_retries", 0),
                               exe=os.path.basename(exe), exe_sha256=exe_sha, date=now)
        with open(args.refs, "w", encoding="utf-8") as f:
            json.dump(new, f, indent=1, ensure_ascii=False, sort_keys=True)
        print(f"[bitid] ancre ecrite : {args.refs} ({len(decks)} decks"
              f"{', MISE A JOUR PARTIELLE --only' if args.only else ''}) — ajouter la ligne CHANGELOG et committer")
        rc = 0
    else:
        verdict = "IDENTIQUE" if n_id == n else "ECHEC"
        extra = []
        if n_diff:  extra.append(f"{n_diff} different(s)")
        if n_noref: extra.append(f"{n_noref} sans reference")
        if n_fail:  extra.append(f"{n_fail} echec(s) de run")
        print(f"[bitid] {n_id}/{n} {verdict}" + (" — " + ", ".join(extra) if extra else ""))
        rc = 0 if n_id == n else (2 if n_fail else 1)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(dict(exe=exe, exe_sha256=exe_sha, date=time.strftime("%Y-%m-%d %H:%M:%S"),
                           refs=args.refs, update=args.update, results=results), f, indent=1, ensure_ascii=False)
    return rc


if __name__ == "__main__":
    sys.exit(main())
