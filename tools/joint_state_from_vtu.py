# -*- coding: utf-8 -*-
"""
Lit les VTU d'un run fdem3d trame a trame et reconstruit, par facette, la
separation normale dn et le glissement |ds| (paire de noeuds A/B), a partir de
fdem3d_XXXX.vtu (tets + points dedoubles) et fdem3d_joints_XXXX.vtu (triangles
= noeuds cote A, champs damage / bonded / tmIns / dmF / breakMode).

Le cote B n'est pas exporte : on le retrouve par la geometrie de la trame 0
(les deux tetraedres qui partagent les trois positions de la facette).

Usage : python tools/joint_state_from_vtu.py <out_dir> <frame> [<frame> ...] [--pj=<Pa/m>] [--mu=<tanphi>] [--coh=<Pa>]
"""
import sys
import numpy as np
import meshio


def load_topology(out_dir):
    m0 = meshio.read(f"{out_dir}/fdem3d_0000.vtu")
    X0 = m0.points
    tets = None
    for cb in m0.cells:
        if cb.type == "tetra":
            tets = cb.data
    j0 = meshio.read(f"{out_dir}/fdem3d_joints_0000.vtu")
    tris = None
    for cb in j0.cells:
        if cb.type == "triangle":
            tris = cb.data
    elemOf = -np.ones(len(X0), int)
    for e, t in enumerate(tets):
        elemOf[t] = e
    # cle geometrique d'une face : positions arrondies triees
    scale = 1e-9
    key_of_pos = lambda p: tuple(np.round(p / scale).astype(np.int64))
    faces = {}
    for e, t in enumerate(tets):
        for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
            ids = [t[q] for q in f]
            k = tuple(sorted(key_of_pos(X0[i]) for i in ids))
            faces.setdefault(k, []).append((e, ids))
    a = tris.copy()
    b = np.zeros_like(a)
    eA = elemOf[a[:, 0]]
    eB = -np.ones(len(a), int)
    for m, tri in enumerate(a):
        k = tuple(sorted(key_of_pos(X0[i]) for i in tri))
        cands = faces.get(k, [])
        other = [(e, ids) for (e, ids) in cands if e != eA[m]]
        if len(other) != 1:
            b[m] = -1
            continue
        eB[m], ids = other[0]
        pos = {key_of_pos(X0[i]): i for i in ids}
        for q in range(3):
            b[m, q] = pos[key_of_pos(X0[tri[q]])]
    return X0, tets, a, b, eA, eB


def joint_kinematics(P, a, b):
    """dn (moyenne des 3 paires), |ds| (moyenne), normale sortante de A."""
    ok = b[:, 0] >= 0
    M = np.array([0.5 * (P[a[:, q]] + P[b[:, q]]) for q in range(3)])  # 3,M,3
    nr = np.cross(M[1] - M[0], M[2] - M[0])
    nn = np.linalg.norm(nr, axis=1)
    n = nr / np.maximum(nn, 1e-30)[:, None]
    dn = np.zeros((len(a), 3))
    ds = np.zeros((len(a), 3))
    for q in range(3):
        d = P[b[:, q]] - P[a[:, q]]
        dnq = np.einsum("ij,ij->i", d, n)
        dn[:, q] = dnq
        ds[:, q] = np.linalg.norm(d - dnq[:, None] * n, axis=1)
    return dn, ds, ok


if __name__ == "__main__":
    args = [x for x in sys.argv[1:] if not x.startswith("--")]
    opts = dict(x[2:].split("=") for x in sys.argv[1:] if x.startswith("--"))
    out_dir = args[0]
    frames = [int(x) for x in args[1:]]
    pj = float(opts.get("pj", 10 * 60e9 / 1.5e-3))
    mu = float(opts.get("mu", np.tan(np.radians(35))))
    coh = float(opts.get("coh", 30e6))
    X0, tets, a, b, eA, eB = load_topology(out_dir)
    print(f"{out_dir}: {len(tets)} tets, {len(a)} facettes, cote B retrouve pour {(b[:,0]>=0).sum()}")
    for fr in frames:
        m = meshio.read(f"{out_dir}/fdem3d_{fr:04d}.vtu")
        j = meshio.read(f"{out_dir}/fdem3d_joints_{fr:04d}.vtu")
        P = m.points
        cd = j.cell_data
        D = np.array(cd["damage"][0]); bonded = np.array(cd["bonded"][0]) > 0.5
        bm = np.array(cd["breakMode"][0])
        tmax = np.array(m.cell_data["tauMax"][0])
        dn, ds, ok = joint_kinematics(P, a, b)
        dnm = dn.mean(1); dsm = ds.mean(1)
        ins = (~bonded) & ok
        brk = ins & (D >= 1.0)
        comp = brk & (dnm < 0)
        print(f"--- trame {fr}: inserees {ins.sum()}, rompues {brk.sum()}, rompues en COMPRESSION (dn<0) {comp.sum()} "
              f"({100.0*comp.sum()/max(1,brk.sum()):.0f} %)")
        if comp.sum():
            tn = pj * dnm[comp]
            tau = mu * np.abs(tn)
            print(f"    rompues comprimees : |dn| um p50 {np.median(-dnm[comp])*1e6:.3f} p90 {np.percentile(-dnm[comp],90)*1e6:.3f} max {(-dnm[comp]).max()*1e6:.3f}"
                  f" ; |ds| um p50 {np.median(dsm[comp])*1e6:.3f} p90 {np.percentile(dsm[comp],90)*1e6:.3f}")
            print(f"    -> t_n = pj dn : MPa p50 {np.median(tn)/1e6:.1f} max {tn.min()/1e6:.1f} ; tether camacho |tau| = mu|t_n| : MPa p50 {np.median(tau)/1e6:.1f} p90 {np.percentile(tau,90)/1e6:.1f} max {tau.max()/1e6:.1f} (cohesion d insertion c = {coh/1e6:.0f} MPa)")
            small = comp & (dsm < 1e-8)
            print(f"    rompues comprimees a |ds| < 10 nm (tether a pleine valeur sur un glissement quasi nul) : {small.sum()}")
            # contrainte nominale des elements voisins
            nb = np.unique(np.concatenate([eA[comp], eB[comp]]))
            print(f"    tauMax des {len(nb)} elements adjacents : MPa p50 {np.median(tmax[nb])/1e6:.1f} p90 {np.percentile(tmax[nb],90)/1e6:.1f} max {tmax[nb].max()/1e6:.1f} ; "
                  f"tous elements : p50 {np.median(tmax)/1e6:.1f} max {tmax.max()/1e6:.1f}")
        if ins.sum():
            print(f"    toutes inserees : dn um p50 {np.median(dnm[ins])*1e6:+.3f} (dn<0 : {100.0*(dnm[ins]<0).mean():.0f} %) ; |ds| um p50 {np.median(dsm[ins])*1e6:.3f} max {dsm[ins].max()*1e6:.2f}")
