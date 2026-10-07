# -*- coding: utf-8 -*-
"""
Banc UN POINT D'INTEGRATION — port en Python, terme pour terme, de la branche
`jointTSL = camacho` de Fdem3dSolver::jointForces (src/Fdem3dSolver.cpp
~3743-3880) et de la branche `penalty` (retour radial, ~4160-4180), avec les
fonctions pures de include/rockim/JointTsl.hpp.

Question : le travail fourni par le joint aux noeuds, jw = trac . (vA - vB) dt
(la meme somme que le poste eJnt de history.csv), peut-il etre POSITIF sur un
cycle FERME de separation (dn, ds) ?  Un joint cohesif + frottant physique ne
peut que stocker (elastique, restitue) ou dissiper : sur un cycle ferme,
W <= 0.  Si W > 0, le joint CREE de l'energie.

Cycle teste (joint en COMPRESSION, D fige) :
   1. dn = -d1, ds : 0 -> s
   2. ds = s,  dn : -d1 -> -d2   (d2 > d1 : on comprime davantage)
   3. dn = -d2, ds : s -> 0
   4. ds = 0,  dn : -d2 -> -d1
C'est le mouvement d'une facette sous le bouton : l'onde de compression fait
osciller dn, le cisaillement de la matrice fait osciller ds, avec un dephasage.

Aucun fichier du depot n'est modifie. Lancer :  python tools/facet_cycle_camacho.py
"""
import math
import numpy as np


def mac(x):
    return x if x > 0.0 else 0.0


# ---------------------------------------------------------------- JointTsl.hpp
class Stamp:
    def __init__(self, tn, ts, tn0, ts0, GIc, GIIc, eta, rise):
        self.beta = ts0 / tn0 if (tn0 > 0 and ts0 > 0) else 1.0
        a = mac(tn)
        b = ts / self.beta
        self.tmIns = math.sqrt(a * a + b * b)
        if eta > 0.0:
            if a > 0:
                r = abs(ts) / a
                m = r * r / (self.beta ** 2 + r * r)
            else:
                m = 1.0
            self.Gc = GIc + (GIIc - GIc) * m ** eta
        else:
            self.Gc = GIc
        self.dmF = 2.0 * self.Gc / self.tmIns if self.tmIns > 0 else 0.0
        self.en = a / self.tmIns
        self.es = (abs(ts) / self.beta) / self.tmIns
        self.dm0 = rise * self.dmF if rise > 0 else 0.0

    def ok(self):
        return self.tmIns > 0 and self.dmF > 0


def traction(S, dmEff, dmMaxEff):
    dmx = max(dmMaxEff, dmEff)
    g = dmx - S.dm0
    if g >= S.dmF:
        return 0.0
    env = S.tmIns * (1.0 - mac(g) / S.dmF)
    if dmEff >= dmx:
        return env
    return env * (dmEff / dmx) if dmx > 0 else 0.0


def damage(S, dmMaxEff):
    return min(1.0, max(0.0, (dmMaxEff - S.dm0) / S.dmF))


def split(tm, dm, dn, beta):
    if not (dm > 0) or not (tm > 0):
        return 0.0, 0.0
    k = tm / dm
    return k * mac(dn), k * beta * beta


def shearCap(cohEff, mu, tnNeg, D, mobilised):
    fric = mu * mac(-tnNeg)
    lim = cohEff + (D * fric if mobilised else fric)
    return lim if lim > 0 else 0.0


# ------------------------------------------- branche camacho (3743-3880)
class CamachoPoint:
    """Un point d'integration k de la branche camacho. n = ez ; le plan de
    la facette est (ex, ey). tsDir = ex (traction tamponnee en cisaillement
    selon ex)."""

    def __init__(self, S, pj, tanPhi, dmMax0=None):
        self.S, self.pj, self.mu = S, pj, tanPhi
        self.dmMax = S.dm0 if dmMax0 is None else dmMax0
        self.D = 0.0

    def force(self, dn, ds_vec):
        S = self.S
        ds = np.asarray(ds_vec, float)
        tdir = np.array([1.0, 0.0, 0.0])
        offS = S.dm0 * S.es / S.beta if S.beta > 0 else 0.0
        dsEffV = ds + offS * tdir
        dsEffN = float(np.linalg.norm(dsEffV))
        dnEff = mac(dn) + S.dm0 * S.en
        dmEff = math.sqrt(dnEff ** 2 + S.beta ** 2 * dsEffN ** 2)
        if dmEff > self.dmMax:
            self.dmMax = dmEff
        dmx = self.dmMax
        tm = traction(S, dmEff, dmx)
        Dk = damage(S, dmx)
        if Dk > self.D:
            self.D = min(1.0, Dk)
        tnCoh, tsScale = split(tm, dmEff, dnEff, S.beta)
        tn = tnCoh if dn >= 0.0 else tnCoh + self.pj * dn
        tauCohMag = tsScale * dsEffN
        lim = shearCap(tauCohMag, self.mu, tn, Dk, True)      # fricMob = damage
        tau = (lim / dsEffN) * dsEffV if dsEffN > 1e-30 else np.zeros(3)
        return tn, tau


# ------------------------------------------- branche penalty (plastic)
class PenaltyPoint:
    def __init__(self, pj, ft, coh, tanPhi, Gf, GfII, dn0=0.0, D=0.0):
        self.pj, self.ft, self.coh, self.mu = pj, ft, coh, tanPhi
        self.dnE = ft / pj
        self.dnF = self.dnE + 2.0 * Gf / ft          # aire = Gf (forme lineaire)
        self.slipF = 2.0 * GfII / coh
        self.dn0 = dn0
        self.slip = np.zeros(3)
        self.D = D

    def force(self, dn_geo, ds_vec):
        dn = dn_geo + self.dn0
        dt3 = np.asarray(ds_vec, float)
        if dn >= 0.0:
            if dn <= self.dnE:
                env = self.pj * dn
            elif dn >= self.dnF:
                env = 0.0
            else:
                env = self.ft * (self.dnF - dn) / (self.dnF - self.dnE)
            tr2 = (1.0 - self.D) * self.pj * dn
            if tr2 > env:
                sig = env
                if dn > 1e-30:
                    Dn = 1.0 - env / (self.pj * dn)
                    if Dn > self.D:
                        self.D = min(1.0, Dn)
            else:
                sig = tr2
        else:
            sig = self.pj * dn
        coh = (1.0 - self.D) * self.coh
        tauLim = shearCap(coh, self.mu, sig, self.D, True)
        tauTr = self.pj * (dt3 - self.slip)
        tnorm = float(np.linalg.norm(tauTr))
        tau = tauTr.copy()
        if tnorm > tauLim and tnorm > 0:
            tau *= tauLim / tnorm
            self.slip += (tauTr - tau) / self.pj                # retour radial
            Dt2 = float(np.linalg.norm(self.slip)) / self.slipF
            if Dt2 > self.D:
                self.D = min(1.0, Dt2)
        return sig, tau


def cycle_work(point, d1, d2, s, nsub=400):
    """Travail du joint sur les noeuds, W = -integral trac . d(delta)
    (= somme des trac . (vA - vB) dt du solveur), sur le cycle ferme."""
    n = np.array([0, 0, 1.0])
    path = []
    for i in range(nsub + 1):                       # 1
        path.append((-d1, np.array([s * i / nsub, 0, 0])))
    for i in range(1, nsub + 1):                    # 2
        path.append((-d1 - (d2 - d1) * i / nsub, np.array([s, 0, 0])))
    for i in range(1, nsub + 1):                    # 3
        path.append((-d2, np.array([s * (1 - i / nsub), 0, 0])))
    for i in range(1, nsub + 1):                    # 4
        path.append((-d2 + (d2 - d1) * i / nsub, np.array([0.0, 0, 0])))
    W = 0.0
    prev = None
    for dn, ds in path:
        tn, tau = point.force(dn, ds)
        trac = tn * n + tau
        if prev is not None:
            ddelta = (dn - prev[0]) * n + (ds - prev[1])
            # regle du point milieu sur la force
            W += -0.5 * np.dot(trac + prev[2], ddelta)
        prev = (dn, ds, trac)
    return W


if __name__ == "__main__":
    # Valeurs du deck loi_note_2026 (facette moyenne du bloc 40 mm : h ~ 2 mm)
    E, h = 60e9, 2.0e-3
    pj = 10.0 * E / h                                   # insertionPenaltyFactor = 10
    ft, coh, phi = 11.4e6, 30e6, 35.0
    mu = math.tan(math.radians(phi))
    Gf, GfII = 100.0, 1000.0
    A = 1.0e-6                                          # 1 mm^2, pour lire en J

    print("pj = %.3e Pa/m, mu = %.3f" % (pj, mu))
    print()
    print("Cycle ferme en compression : dn -d1 -> -d2, ds 0 -> s -> 0 ; W > 0 = ENERGIE CREEE")
    print("%-34s %10s %10s %10s %14s %14s" % ("cas", "d1 [um]", "d2 [um]", "s [um]",
                                              "W camacho [J]", "W penalty [J]"))
    # joint insere en cisaillement pur (tn <= 0, ts = coh) : es = 1, en = 0
    for D0, label in ((1.0, "rompu D=1"), (0.5, "D=0.5"), (0.1, "D=0.1")):
        for d1, d2, s in ((0.5, 1.0, 0.5), (0.5, 1.0, 0.05), (0.1, 0.2, 0.01), (1.0, 3.0, 1.0)):
            S = Stamp(-1.0e6, coh, ft, coh, Gf, GfII, 0.0, 1e-3)   # jointMixLaw = none
            dmMax0 = S.dm0 + D0 * S.dmF
            cp = CamachoPoint(S, pj, mu, dmMax0=dmMax0)
            pp = PenaltyPoint(pj, ft, coh, mu, Gf, GfII, dn0=0.0, D=D0)
            um = 1e-6
            Wc = cycle_work(cp, d1 * um, d2 * um, s * um) * A
            Wp = cycle_work(pp, d1 * um, d2 * um, s * um) * A
            print("%-34s %10.2f %10.2f %10.3f %14.3e %14.3e" % (label, d1, d2, s, Wc, Wp))
            # formule fermee attendue pour camacho a D fige :
            # W = D mu pj (d2 - d1) s  (le terme frottant D mu <-tn> le long de dsEffV)
        Wth = D0 * mu * pj * (1.0 - 0.5) * 1e-6 * 0.5e-6 * A
        print("   -> formule D mu pj (d2-d1) s pour la 1re ligne : %.3e J" % Wth)
    print()
    print("Meme cycle, mu = 0 (frictionDeg = 0) : le terme frottant disparait")
    S = Stamp(-1.0e6, coh, ft, coh, Gf, GfII, 0.0, 1e-3)
    cp = CamachoPoint(S, pj, 0.0, dmMax0=S.dm0 + 1.0 * S.dmF)
    print("   W camacho (D=1, mu=0) = %.3e J" % (cycle_work(cp, 0.5e-6, 1e-6, 0.5e-6) * A))
    print()
    print("Discontinuite en ds = 0 (joint rompu, dn = -1 um) : tau(+eps) et tau(-eps)")
    cp = CamachoPoint(S, pj, mu, dmMax0=S.dm0 + 1.0 * S.dmF)
    for eps in (1e-9, 1e-12, -1e-12, -1e-9):
        tn, tau = cp.force(-1e-6, [eps, 0, 0])
        print("   ds = %+.0e m : tn = %+.3e Pa, tau_x = %+.3e Pa" % (eps, tn, tau[0]))
    print("   (la loi de penalite donne pj*ds = %+.3e Pa a ds = 1e-9)" % (pj * 1e-9))
