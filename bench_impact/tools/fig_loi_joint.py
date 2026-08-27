# -*- coding: utf-8 -*-
"""Loi de joint cohesif du schema ADAPTATIF de rockim (Yan 2023 / Guo 2014).
Style article : PDF vectoriel, Latin Modern (Computer Modern Unicode) + mathtext cm."""
import numpy as np, glob, matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm, matplotlib.pyplot as plt
for f in glob.glob('/usr/share/texmf/fonts/opentype/public/lm/lmroman10-*.otf'):
    fm.fontManager.addfont(f)
plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['Latin Modern Roman', 'DejaVu Serif'],
    'mathtext.fontset': 'cm', 'axes.unicode_minus': True,
    'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 9.5,
    'legend.fontsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'axes.linewidth': 0.7, 'lines.linewidth': 1.6,
})
RED, AMB, GRY = '#c1121f', '#b5730a', '#5a5a5a'      # palette validee (6 checks)

E, ft, c, tanphi = 57e9, 7.0e6, 18.8e6, 1.0
GI, GII = 12.0, 800.4
A, B, C = 0.63, 1.8, 6.0
fD = lambda D: (1-(A+B-1)/(A+B)*np.exp(D*(A+B*C)/((A+B)*(1-A-B))))*(A*(1-D)+B*(1-D)**C)
gd = np.linspace(0, 1, 200001); I = float(np.trapezoid(fD(gd), gd)); kI = 1.0/I
dnc, dsc0 = kI*GI/ft, kI*GII/c
SIG = -500e6; fs = c + tanphi*abs(SIG); dscC = kI*GII/fs
D = np.linspace(0, 1, 500)

fig, ax = plt.subplots(2, 2, figsize=(10.8, 7.8))
def law(a, xc, pk, col, ls, lab=None, z=3, fill=True):
    x, y = D*xc*1e6, fD(D)*pk
    a.plot(x, y, color=col, ls=ls, lw=1.9, label=lab, zorder=z)
    if fill: a.fill_between(x, 0, y, color=col, alpha=0.13, lw=0, zorder=1)
    return x, y

# (a) mode I : intrinseque vs adaptatif
a = ax[0,0]; el = 0.55
a.plot([0, el], [0, ft/1e6], color=GRY, ls=':', lw=1.7, zorder=2)
a.plot(D*dnc*1e6+el, fD(D)*ft/1e6, color=GRY, ls=':', lw=1.7, zorder=2,
       label="intrinsèque : branche élastique, puis adoucissement")
law(a, dnc, ft/1e6, RED, '-', "adaptatif : inséré AU pic, adoucissement seul")
a.annotate("", xy=(0, 7.55), xytext=(el, 7.55), arrowprops=dict(arrowstyle='<->', color='k', lw=0.7))
a.text(el/2, 7.75, r"$\delta_{np}$", ha='center', fontsize=8)
a.text(2.25, 0.35, "aire $= G_I$  (identique)", color=RED, fontsize=8.5)
a.plot(0, ft/1e6, 'o', color=RED, ms=5, zorder=6)
a.annotate("insertion", xy=(0.02, ft/1e6), xytext=(1.15, 7.9), color=RED, fontsize=7.8,
           arrowprops=dict(arrowstyle="->", color=RED, lw=0.8))
a.set_xlim(-0.15, 5.5); a.set_ylim(0, 8.6)
a.set_xlabel(r"ouverture $\delta_n$  [$\mu$m]"); a.set_ylabel(r"$\sigma$  [MPa]")
a.set_title("(a)  Mode I — le joint adaptatif naît au pic", loc='left')
a.legend(loc='upper right', frameon=False)

# (b) enveloppe de Mohr-Coulomb a coupure
b = ax[0,1]
sn = np.linspace(-42e6, ft, 300)
b.plot(sn/1e6, (c-sn*tanphi)/1e6, color=RED, lw=1.9, label=r"$f_s=c-\sigma_n\tan\phi$")
b.plot([ft/1e6]*2, [0, (c-ft*tanphi)/1e6], color=RED, lw=1.9)
b.axvline(0, color='k', lw=0.6); b.axhline(0, color='k', lw=0.6)
b.plot(0, c/1e6, 'o', color=RED, ms=4.5, zorder=5); b.text(1.4, c/1e6+1.2, "$c$", fontsize=9)
b.text(ft/1e6+0.6, 2.5, "coupure $f_t$", fontsize=7.6, color=RED)
b.annotate("plus c'est comprimé,\nplus le pic est haut",
           xy=(-36, (c+36e6)/1e6), xytext=(-30, 20), fontsize=8, color=GRY,
           arrowprops=dict(arrowstyle='->', lw=0.7, color=GRY))
b.text(-40, 6, "sous l'insert : $\\sigma_n \\sim -0{,}5$ à $-1$ GPa,\n"
               "donc $f_s \\sim 0{,}5$ à $1$ GPa (hors échelle)", fontsize=7.6, color=GRY)
b.set_xlim(-42, 10); b.set_ylim(0, 66)
b.set_xlabel(r"$\sigma_n$  [MPa]      (traction $>0$)"); b.set_ylabel(r"$\tau$  [MPa]")
b.set_title(r"(b)  Le pic du mode II, c'est $f_s(\sigma_n)$", loc='left')
b.legend(loc='upper right', frameon=False)

# (c) mode II sans confinement
cc = ax[1,0]
law(cc, dsc0, c/1e6, RED, '-', r"adaptatif, $\sigma_n=0$   (pic $=c$)")
cc.text(46, 12.4, "aire $= G_{II}$\n"+r"$\delta_{sc}=3G_{II}/c=110\ \mu$m", color=RED, fontsize=8.5)
cc.set_xlim(-1.5, 116); cc.set_ylim(0, 22.6)
cc.set_xlabel(r"glissement $\delta_s$  [$\mu$m]"); cc.set_ylabel(r"$\tau$  [MPa]")
cc.set_title("(c)  Mode II sans confinement", loc='left')
cc.legend(loc='upper right', frameon=False)

# (d) mode II sous confinement : le verrou
d = ax[1,1]
law(d, dsc0, fs/1e6, AMB, '--', z=2)
xb, yb = law(d, dscC, fs/1e6, RED, '-', z=3)
d.annotate("plage $3G_{II}/f_s(\\sigma_n)$ — la loi publiée\n"
           + r"$\delta_{sc}=4{,}0\ \mu$m,   aire $=G_{II}$",
           xy=(4.6, 300), xytext=(11, 470), color=RED, fontsize=8.4,
           arrowprops=dict(arrowstyle='->', color=RED, lw=0.9))
d.text(19, 52, "plage FIGÉE $3G_{II}/c$ — la transcription\n"
       + r"$\delta_{sc}=110\ \mu$m,   aire $=22$ kJ/m$^2 = 28\,G_{II}$",
       color=AMB, fontsize=8.4)
d.set_xlim(-1.5, 116); d.set_ylim(0, 620)
d.set_xlabel(r"glissement $\delta_s$  [$\mu$m]"); d.set_ylabel(r"$\tau$  [MPa]")
d.set_title(r"(d)  Mode II à $\sigma_n=-500$ MPa — sous l'insert", loc='left')
iz = d.inset_axes([0.63, 0.50, 0.34, 0.42])
iz.plot(xb, yb, color=RED, lw=1.5); iz.fill_between(xb, 0, yb, color=RED, alpha=0.15, lw=0)
iz.set_xlim(0, 4.2); iz.set_ylim(0, 560); iz.tick_params(labelsize=6)
iz.set_title(r"zoom $\times 28$", fontsize=7)
for s in iz.spines.values(): s.set_linewidth(0.5)

for a_ in ax.flat:
    a_.grid(alpha=0.16, lw=0.4); a_.set_axisbelow(True)
    for s in ('top', 'right'): a_.spines[s].set_visible(False)
fig.suptitle("Loi de joint cohésif du schéma ADAPTATIF de rockim — calcaire St Anne "
             r"($f_t=7$ MPa, $c=18{,}8$ MPa, $\phi=45^\circ$, $G_I=12$, $G_{II}=800$ J/m$^2$)",
             fontsize=10, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.955])
P = '/tmp/claude-0/-home-user/0263e025-847d-592e-aea7-7fec643bb1d6/scratchpad/fig_loi_joint'
fig.savefig(P+'.pdf'); fig.savefig(P+'.png', dpi=170)
print("kI=%.4f | mode I %.2f um | mode II sigma0 %.1f um | fs=%.0f MPa -> %.2f um" %
      (kI, dnc*1e6, dsc0*1e6, fs/1e6, dscC*1e6))
