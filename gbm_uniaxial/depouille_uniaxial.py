# -*- coding: utf-8 -*-
r"""Depouillement des deux essais uniaxiaux : microstructure GBM contre temoin homogene,
SUR LE MEME MAILLAGE, donc comparables element par element.

Deux questions, deux mesures.
  (1) Le module apparent : ou tombe-t-il entre les bornes de Voigt et de Reuss ?
  (2) LE CONTRASTE DE RAIDEUR CONCENTRE-T-IL LES CONTRAINTES ? On mesure la dispersion de
      la contrainte equivalente (q99 sur mediane), sa repartition par phase, et surtout si
      les maxima se logent AUX FRONTIERES DE GRAIN plutot qu'a l'interieur des grains.

    python depouille_uniaxial.py     ->  UNIAXIAL.md + figures/uniaxial_synthese.pdf/.png
"""
from __future__ import print_function
import os, re, io
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.family': 'serif', 'font.serif': ['CMU Serif', 'DejaVu Serif'],
                     'mathtext.fontset': 'cm', 'font.size': 9, 'pdf.fonttype': 42})
PH = ['feldspath', 'quartz', 'biotite']
COL = ['#d9c8a9', '#8fb8de', '#4a3b2a']
FRAC = np.array([0.62, 0.31, 0.07])
EPH = np.array([70.0, 83.1, 29.3])                      # GPa


def vtu(path):
    s = io.open(path, encoding='utf-8', errors='ignore').read()
    out = {}
    for m in re.finditer(r'<DataArray([^>]*)>(.*?)</DataArray>', s, re.S):
        nm = re.search(r'Name="([^"]+)"', m.group(1))
        if nm:
            out[nm.group(1)] = np.fromstring(m.group(2).strip(), sep=' ')
    p = re.search(r'<Points>\s*<DataArray[^>]*>(.*?)</DataArray>', s, re.S)
    return np.fromstring(p.group(1).strip(), sep=' ').reshape(-1, 3), out


def hist(path):
    rows = [l.split(',') for l in io.open(path, encoding='utf-8').read().splitlines()]
    h = rows[0]
    d = dict((k, []) for k in h)
    for r in rows[1:]:
        if len(r) != len(h):
            continue
        try:
            v = [float(x) for x in r]
        except ValueError:
            break
        for k, x in zip(h, v):
            d[k].append(x)
    return dict((k, np.array(v)) for k, v in d.items())


def distance_au_joint(P, conn, grain):
    """Distance du centre de chaque element au JOINT DE GRAIN le plus proche, en metres.

    La definition naive - partager une face avec un autre grain - attrape 91 % du volume
    quand un grain ne compte qu'une quarantaine de tetraedres : elle ne discrimine rien.
    On mesure donc une vraie distance : centre des faces inter-grains, puis plus proche
    voisin par arbre k-d.
    """
    from scipy.spatial import cKDTree
    fi = [(1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)]
    F = np.concatenate([conn[:, list(f)] for f in fi])
    own = np.tile(np.arange(len(conn)), 4)
    K = np.sort(F, axis=1)
    o = np.lexsort((K[:, 2], K[:, 1], K[:, 0]))
    K, F, own = K[o], F[o], own[o]
    eq = np.all(K[1:] == K[:-1], axis=1)
    a, b = own[:-1][eq], own[1:][eq]
    inter = np.where(grain[a] != grain[b])[0]
    faces = F[:-1][eq][inter]
    cf = P[faces].mean(axis=1)
    cen = P[conn].mean(axis=1)
    return cKDTree(cf).query(cen)[0], len(cf)


def quantiles(x, w):
    o = np.argsort(x)
    xs, ws = x[o], w[o]
    c = np.cumsum(ws) / ws.sum()
    q = lambda p: float(np.interp(p, c, xs))
    return q(0.5), q(0.99), q(0.999), float(xs[-1])


def seuil(x, w, p):
    o = np.argsort(x)
    return float(np.interp(p, np.cumsum(w[o]) / w.sum(), x[o]))


def main():
    last = sorted([f for f in os.listdir('out_U_gbm') if f.endswith('.vtu')])[-1]
    P, Ag = vtu(os.path.join('out_U_gbm', last))
    _, Ah = vtu(os.path.join('out_U_homo', last))
    conn = Ag['connectivity'].astype(np.int64).reshape(-1, 4)
    ph = Ag['phase'].astype(int)
    grain = Ag['grain'].astype(int)
    vol = np.abs(np.einsum('ij,ij->i',
                           np.cross(P[conn[:, 1]] - P[conn[:, 0]], P[conn[:, 2]] - P[conn[:, 0]]),
                           P[conn[:, 3]] - P[conn[:, 0]])) / 6.0
    sg = Ag['vonMises'] / 1e6
    sh = Ah['vonMises'] / 1e6
    dist, nfaces = distance_au_joint(P, conn, grain)
    Hg, Hh = hist('out_U_gbm/history.csv'), hist('out_U_homo/history.csv')

    L = []
    L.append(u"# Essai uniaxial : microstructure GBM contre temoin homogene\n")
    L.append(u"*Meme maillage (1001 grains de Voronoi, 40 724 tetraedres, eprouvette "
             u"16 x 16 x 32 mm), meme chargement, loi elastique. Seul le materiau change : trois "
             u"phases minerales contre un materiau unique aux moyennes de Voigt. Toute difference "
             u"vient donc du CONTRASTE DE RAIDEUR et de rien d'autre. Image depouillee : %s.*\n"
             % last)
    L.append(u"\n> **Reserve emise par le solveur lui-meme.** Le coefficient de Poisson varie de "
             u"0,17 a 0,36 entre les phases. Les tetraedres lineaires a un point d'integration se "
             u"raidissent d'autant plus que ce coefficient est grand : l'artefact est CORRELE A LA "
             u"PHASE et SOUS-ESTIME le contraste. Les ecarts ci-dessous sont des bornes basses.\n")
    L.append(u"\n> **Seconde reserve.** L'eprouvette est un parallelepipede et non un cylindre, et "
             u"le maillage intragranulaire est l'eventail par defaut de la tessellation interne. "
             u"Ces deux points sont a corriger avant tout usage en rapport.\n")

    # 1 -- module apparent
    Ev = float((FRAC * EPH).sum())
    Er = float(1.0 / (FRAC / EPH).sum())
    L.append(u"\n## 1. Module apparent\n")
    L.append(u"| grandeur | valeur |")
    L.append(u"|---|---|")
    L.append(u"| borne de Voigt, moyenne des modules | %.2f GPa |" % Ev)
    L.append(u"| borne de Reuss, moyenne des souplesses | %.2f GPa |" % Er)
    mods = {}
    for nm, H in ((u'GBM', Hg), (u'temoin homogene', Hh)):
        s, e = H['sigZZmid'], H['epsAxMid']
        amp = np.abs(e).max()
        m = (np.abs(e) > 0.25 * amp) & (np.abs(e) < 0.95 * amp)
        if m.sum() > 5:
            E = float(np.polyfit(e[m], s[m], 1)[0]) / 1e9
            mods[nm] = E
            L.append(u"| module mesure, %s | **%.2f GPa** |" % (nm, E))
    if len(mods) == 2:
        a, b = mods[u'GBM'], mods[u'temoin homogene']
        L.append(u"| ecart GBM / temoin | **%+.2f %%** |" % (100.0 * (a - b) / b))
    L.append(u"\nLe temoin doit retomber sur Voigt, puisqu'on lui a donne cette moyenne. Le GBM "
             u"doit tomber EN DESSOUS : un assemblage reel est toujours plus souple que la moyenne "
             u"des modules, parce que les grains souples se deforment davantage.\n")

    # 2 -- concentration
    mg, q99g, q999g, mxg = quantiles(sg, vol)
    mh, q99h, q999h, mxh = quantiles(sh, vol)
    L.append(u"\n## 2. Concentration de contrainte\n")
    L.append(u"Contrainte equivalente de von Mises, ponderee par le volume des elements.\n")
    L.append(u"| | mediane | quantile 99 % | quantile 99,9 % | maximum | q99 / mediane |")
    L.append(u"|---|---|---|---|---|---|")
    L.append(u"| GBM | %.1f MPa | %.1f | %.1f | %.1f | **%.2f** |"
             % (mg, q99g, q999g, mxg, q99g / mg))
    L.append(u"| temoin homogene | %.1f MPa | %.1f | %.1f | %.1f | **%.2f** |"
             % (mh, q99h, q999h, mxh, q99h / mh))
    L.append(u"| rapport GBM / temoin | %.2f | %.2f | %.2f | %.2f | %.2f |"
             % (mg / mh, q99g / q99h, q999g / q999h, mxg / mxh, (q99g / mg) / (q99h / mh)))
    L.append(u"\nLa derniere colonne est la mesure directe de la concentration : de combien la "
             u"queue haute depasse le coeur de la distribution. Le temoin homogene donne la "
             u"dispersion due au seul maillage et aux bords ; tout ce que le GBM a en plus vient "
             u"du contraste.\n")

    # 3 -- par phase
    L.append(u"\n## 3. Repartition par phase\n")
    L.append(u"| phase | fraction volumique | module | sigma_eq moyen, GBM | idem temoin | rapport |")
    L.append(u"|---|---|---|---|---|---|")
    for i, nm in enumerate(PH):
        m = ph == i
        if not m.any():
            continue
        a = float(np.average(sg[m], weights=vol[m]))
        b = float(np.average(sh[m], weights=vol[m]))
        L.append(u"| %s | %.1f %% | %.1f GPa | %.1f MPa | %.1f MPa | **%.2f** |"
                 % (nm, 100 * vol[m].sum() / vol.sum(), EPH[i], a, b, a / b))
    L.append(u"\nUn contraste de raideur charge les grains RAIDES et decharge les souples. Si les "
             u"rapports suivent l'ordre des modules (quartz > feldspath > biotite), l'effet est "
             u"bien elastique et non un artefact.\n")

    # 4 -- frontieres de grain, par la DISTANCE
    L.append(u"\n## 4. Les maxima se logent-ils aux frontieres de grain ?\n")
    L.append(u"La definition naive, *partager une face avec un autre grain*, attrape 91 % du "
             u"volume quand un grain ne compte qu'une quarantaine de tetraedres : elle ne "
             u"discrimine rien. On mesure donc la DISTANCE du centre de chaque element au joint "
             u"de grain le plus proche, rapportee au rayon de grain equivalent.\n")
    R = (3.0 * vol.sum() / (4.0 * np.pi * len(np.unique(grain)))) ** (1.0 / 3.0)
    dn = dist / R
    L.append(u"Rayon de grain equivalent : %.2f mm. Distance mediane au joint : %.2f mm, "
             u"soit %.2f rayon.\n" % (1e3 * R, 1e3 * float(np.median(dist)), float(np.median(dn))))
    L.append(u"| population | distance mediane au joint, GBM | idem temoin | ecart |")
    L.append(u"|---|---|---|---|")
    for p, lab in ((0.90, u'les 10 % les plus charges'), (0.99, u'les 1 % les plus charges'),
                   (0.999, u'le 0,1 % le plus charge')):
        mg2 = sg >= seuil(sg, vol, p)
        mh2 = sh >= seuil(sh, vol, p)
        a = float(np.median(dn[mg2])); b = float(np.median(dn[mh2]))
        L.append(u"| %s | %.3f rayon | %.3f rayon | **%+.0f %%** |" % (lab, a, b, 100 * (a - b) / b))
    L.append(u"\nUne distance PLUS COURTE que celle du temoin veut dire que les maxima se "
             u"rapprochent des joints. Le temoin, materiau uniforme sur le meme maillage, donne "
             u"la reference : ce qu'on lit chez lui vient du seul maillage.\n")
    q = [0.0, 0.15, 0.30, 0.50, 1.0, 10.0]
    L.append(u"| tranche de distance au joint | part du volume | sigma_eq moyen, GBM | idem temoin | rapport |")
    L.append(u"|---|---|---|---|---|")
    for k in range(len(q) - 1):
        m2 = (dn >= q[k]) & (dn < q[k + 1])
        if vol[m2].sum() / vol.sum() < 0.005:
            continue
        a = float(np.average(sg[m2], weights=vol[m2]))
        b = float(np.average(sh[m2], weights=vol[m2]))
        L.append(u"| %.2f a %.2f rayon | %.1f %% | %.2f MPa | %.2f MPa | **%.3f** |"
                 % (q[k], q[k + 1], 100 * vol[m2].sum() / vol.sum(), a, b, a / b))

    io.open('UNIAXIAL.md', 'w', encoding='utf-8').write(u"\n".join(L) + u"\n")
    for l in L:
        print(l)

    # figures
    os.makedirs('figures', exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(9.6, 3.0), constrained_layout=True)
    for nm, H, c in ((u'GBM', Hg, '#c0392b'), (u'homogene', Hh, '#2c6fbb')):
        ax[0].plot(np.abs(H['epsAxMid']) * 100, np.abs(H['sigZZmid']) / 1e6, lw=1.2,
                   color=c, label=nm)
    ax[0].set_xlabel(u'deformation axiale [%]')
    ax[0].set_ylabel(u'contrainte axiale [MPa]')
    ax[0].legend(frameon=False, fontsize=8)
    ax[0].set_title(u'reponse globale', fontsize=9)

    hi = float(np.percentile(sg, 99.8))
    b = np.linspace(0, hi, 70)
    ax[1].hist(sg, bins=b, weights=vol, density=True, histtype='step', color='#c0392b',
               lw=1.3, label=u'GBM')
    ax[1].hist(sh, bins=b, weights=vol, density=True, histtype='step', color='#2c6fbb',
               lw=1.3, label=u'homogene')
    ax[1].set_xlabel(u'$\\sigma_{eq}$ [MPa]')
    ax[1].set_ylabel(u'densite en volume')
    ax[1].legend(frameon=False, fontsize=8)
    ax[1].set_title(u'distribution des contraintes', fontsize=9)

    for i, nm in enumerate(PH):
        m = ph == i
        if m.any():
            ax[2].hist(sg[m], bins=b, weights=vol[m], density=True, histtype='step',
                       color=COL[i], lw=1.4, label=nm)
    ax[2].set_xlabel(u'$\\sigma_{eq}$ [MPa]')
    ax[2].legend(frameon=False, fontsize=8)
    ax[2].set_title(u'GBM, par phase', fontsize=9)
    for e in ('pdf', 'png'):
        fig.savefig('figures/uniaxial_synthese.' + e, dpi=200)
    print('-> UNIAXIAL.md et figures/uniaxial_synthese.pdf et .png')


if __name__ == '__main__':
    main()
