"""Loi cohésive des joints de g1 (jointSoftening = yan), évaluée en un point.

Référence de l'éditeur de loi (spec 007 §2.6). Transcription de tools/yan_point.cpp,
qui reproduit la mise à jour des tractions de FdemSolver::jointForces avec les fonctions
mêmes du solveur (YanSoftening.hpp). Le test tests/test_loi.py compare ce module à ce
programme compilé, point par point.

Mode I (ouverture monotone) : sigma = pj dn jusqu'à dnE = ft / pj, puis f(D) ft avec
D = (dn - dnE) / ot et ot = GfI / (ft I) ; l'aire de la branche adoucissante vaut GfI.
Mode II (glissement à contrainte normale fixée, retour radial) : la limite vaut
f(D) c + tan(phi) max(0, -sigma_n). Seule la COHÉSION s'adoucit, le frottement reste :
sous compression, tau tend vers tan(phi) |sigma_n|, pas vers zéro. D = slip / st avec
st = GfII / (c I) et slip le glissement PLASTIQUE.

Montée parabolique (jointElastic = parabolic, Guo éq. 2.31) : sigma = ft (2r - r^2),
r = dn / dnE, pour 0 <= dn <= dnE. Adoucissement linéaire : f(D) = 1 - D, I = 1/2.
"""
import math

import numpy as np


def fD(D, a=0.63, b=1.8, c=6.0):
    """f(D) de Yan et al. 2023 éq. 11, transcrit de YanSoftening.hpp (bornée à [0, 1])."""
    if D <= 0.0:
        return 1.0
    if D >= 1.0:
        return 0.0
    s = a + b
    if abs(s) < 1e-12 or abs(1.0 - s) < 1e-12:
        return a * (1.0 - D) + b * (1.0 - D) ** c
    K = (s - 1.0) / s
    alpha = (a + c * b) / (s * (1.0 - s))
    f = (1.0 - K * math.exp(alpha * D)) * (a * (1.0 - D) + b * (1.0 - D) ** c)
    return min(1.0, max(0.0, f))


def integrale(a=0.63, b=1.8, c=6.0, n=4096):
    """I = intégrale de f sur [0, 1], Simpson composite (integralFD)."""
    if n % 2:
        n += 1
    h = 1.0 / n
    s = fD(0.0, a, b, c) + fD(1.0, a, b, c)
    for i in range(1, n):
        s += fD(i * h, a, b, c) * (4.0 if i % 2 else 2.0)
    return s * h / 3.0


class Loi:
    def __init__(self, pj, ft, coh, phi_deg, GfI, GfII, a=0.63, b=1.8, c=6.0,
                 adoucissement="yan", montee="linear"):
        self.pj, self.ft, self.coh = pj, ft, coh
        self.tan_phi = math.tan(math.radians(phi_deg))
        self.a, self.b, self.c = a, b, c
        self.yan = adoucissement == "yan"
        self.montee = montee
        self.I = integrale(a, b, c) if self.yan else 0.5
        self.dnE = ft / pj
        self.ot = GfI / (ft * self.I)
        self.st = GfII / (coh * self.I)

    def f(self, D):
        return fD(D, self.a, self.b, self.c) if self.yan else min(1.0, max(0.0, 1.0 - D))

    def sigma(self, dn):
        """Mode I, ouverture monotone."""
        if dn <= self.dnE:
            r = dn / self.dnE
            return self.ft * (2 * r - r * r) if self.montee == "parabolic" else self.pj * dn
        D = min(1.0, (dn - self.dnE) / self.ot)
        return min(self.pj * dn, self.f(D) * self.ft)

    def mode2(self, dtg, sigma_n=0.0):
        """Mode II : tractions le long d'un chemin de glissement croissant dtg (tableau),
        par retour radial pas à pas, exactement comme Point::tau de yan_point."""
        slip, D, out = 0.0, 0.0, []
        for x in dtg:
            tr = self.pj * (x - slip)
            lim = self.f(D) * self.coh + self.tan_phi * max(0.0, -sigma_n)
            t = max(-lim, min(lim, tr))
            if t != tr:
                slip += (tr - t) / self.pj
                D = max(D, min(1.0, abs(slip) / self.st))
            out.append(t)
        return np.array(out)
