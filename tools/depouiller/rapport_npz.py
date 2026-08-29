# -*- coding: utf-8 -*-
"""Rapport P1 (outil MAILLE) depuis un .npz produit par tools/pack_run.py
   usage: rapport_npz.py [chemin.npz] [tag]"""
import sys, glob, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm, matplotlib.pyplot as plt
for f in glob.glob('/usr/share/texmf/fonts/opentype/public/lm/lmroman10-*.otf'):
    fm.fontManager.addfont(f)
plt.rcParams.update({'font.family':'serif','font.serif':['Latin Modern Roman'],
                     'mathtext.fontset':'cm','font.size':9,'axes.linewidth':0.7})
RED,AMB,BLUE,INK='#c1121f','#b5730a','#1f4e9c','#3a3a3a'
SP='/tmp/claude-0/-home-user/0263e025-847d-592e-aea7-7fec643bb1d6/scratchpad'
NPZ=sys.argv[1] if len(sys.argv)>1 else '/home/user/rockim/bench_impact/donnees/P1.npz'
TAG=sys.argv[2] if len(sys.argv)>2 else 'P1'
M_TOOL=1.11813+0.0635586; A_BIT=np.pi*0.015**2

def liss(x,n=21):
    h=n//2; xp=np.concatenate([x[h:0:-1],x,x[-2:-h-2:-1]])   # miroir, jamais zero
    return np.convolve(xp,np.ones(n)/n,mode='valid')

d=np.load(NPZ); h=d['hist']; c=list(d['hist_cols']); I=lambda n:c.index(n)
t=h[:,I('t')]*1e6
# CONVENTION : la penetration se mesure sous la SURFACE LIBRE, pas depuis la
# position initiale de l outil — le maillage laisse un JEU (mesure ici, pas
# supposee) que l outil franchit en vol libre. Meme piege qu au run Crebond.
surf = d['el_xyz'][0][:,:,2].max()          # max des NOEUDS, pas des barycentres
jeu  = float(d['tool_z'][0,0]) - surf
p=((h[0,I('z_insert')]-h[:,I('z_insert')]) - jeu)*1e3
p=np.maximum(p, 0.0)
vb=h[:,I('vz_bit')]; vi=h[:,I('vz_insert')]
Fg=-h[:,I('szz_bit')]*A_BIT/1e3                      # jauge du bit
# route inertielle : derivee de la QUANTITE DE MOUVEMENT des deux corps
# (bit + insert brases). Signe : l outil descend (v < 0) et est freine, donc
# dv/dt > 0 et la force de reaction de la roche vaut +dP/dt (vers le haut).
M_BIT, M_INS = 1.11813, 0.0635586
P_tot = M_BIT*vb + M_INS*vi
Fi = np.gradient(liss(P_tot), h[:,I('t')].astype(np.float64))/1e3
npul=h[:,I('nPulv')]; bdw=-h[:,I('bdWork')]; nb=h[:,I('nBroken')]

fig,ax=plt.subplots(2,2,figsize=(11.6,9.0))
# (a) force-penetration, deux routes
a=ax[0,0]; k=p>1e-4
a.plot(p[k],liss(Fi)[k],color=RED,lw=1.6,label=r'inertie  $\dot P = \mathrm{d}(mv)/\mathrm{d}t$')
a.plot(p[k],liss(Fg)[k],color=BLUE,lw=1.2,ls='--',label='jauge du bit  $\\sigma_{zz}A$')
a.set_xlabel("enfoncement de l'insert  [mm]"); a.set_ylabel('force  [kN]')
a.set_title('(a)  Force–pénétration : deux routes indépendantes',loc='left')
a.legend(fontsize=8,frameon=False,loc='upper left')
_dF = liss(Fi).max()-liss(Fg).max()
a.text(0.03,0.55,f"écart des pics {_dF:+.1f} kN\n(la jauge est à 0,13 m de\nl'insert : transit d'onde\net moyenne sur 29 tets)",
       transform=a.transAxes,fontsize=7.5,color=INK,va='top')
im=int(np.argmax(p)); a.annotate(f"rebroussement\n{p[im]:.2f} mm à {t[im]:.0f} µs",
    xy=(p[im],liss(Fi)[im]),xytext=(0.55,0.72),textcoords='axes fraction',fontsize=8,
    color=INK,arrowprops=dict(arrowstyle='->',color=INK,lw=.8))
a.annotate(f"résiduel {p[-1]:.2f} mm",xy=(p[-1],0),xytext=(0.05,0.12),
    textcoords='axes fraction',fontsize=8,color=INK,
    arrowprops=dict(arrowstyle='->',color=INK,lw=.8))
# (b) vitesses et enfoncement
b=ax[0,1]
b.plot(t,vb,color=BLUE,lw=1.6,label='bit'); b.plot(t,vi,color=INK,lw=.9,ls=':',label='insert')
b.axhline(0,color='0.6',lw=.6)
b2=b.twinx(); b2.plot(t,p,color=RED,lw=1.3,ls='--'); b2.set_ylabel('enfoncement  [mm]',color=RED)
b2.tick_params(axis='y',labelcolor=RED)
b.set_xlabel('temps  [µs]'); b.set_ylabel(r'$v_z$  [m/s]',color=BLUE)
b.tick_params(axis='y',labelcolor=BLUE); b.set_title('(b)  Vitesses des corps et enfoncement',loc='left')
b.legend(fontsize=8,frameon=False,loc='center left')
i0=int(np.argmin(np.abs(vb)))
b.annotate(f"$v=0$ à {t[i0]:.0f} µs",xy=(t[i0],0),xytext=(0.30,0.35),textcoords='axes fraction',
           fontsize=8,color=INK,arrowprops=dict(arrowstyle='->',color=INK,lw=.8))
b.annotate(f"rebond {vb[-1]:+.2f} m/s\n$e = {abs(vb[-1])/9.5:.2f}$",xy=(t[-1],vb[-1]),
    xytext=(0.62,0.80),textcoords='axes fraction',fontsize=8,color=BLUE,
    arrowprops=dict(arrowstyle='->',color=BLUE,lw=.8))
# (c) pulverisation vs force
cc=ax[1,0]
cc.plot(t,liss(Fi),color='0.75',lw=1.0,label='force (inertie)')
cc.set_ylabel('force  [kN]',color='0.45'); cc.tick_params(axis='y',labelcolor='0.45')
c2=cc.twinx(); c2.plot(t,npul,color=AMB,lw=1.8,label="éléments pulvérisés")
c2.plot(t,bdw*10,color=RED,lw=1.3,ls='--',label=r'$10\times$ énergie de pulvérisation [J]')
c2.set_ylabel("éléments à $D=D_{max}$   |   $10\\times$ J",color=AMB)
c2.tick_params(axis='y',labelcolor=AMB)
cc.set_xlabel('temps  [µs]'); cc.set_title('(c)  Pulvérisation — le canal de WP6',loc='left')
nz=np.nonzero(npul>0)[0]
if len(nz):
    cc.axvline(t[nz[0]],color=INK,lw=.7,ls=':')
    cc.annotate(f"1er élément pulvérisé\n{t[nz[0]]:.0f} µs",xy=(t[nz[0]],0),
        xytext=(0.30,0.55),textcoords='axes fraction',fontsize=8,color=INK,
        arrowprops=dict(arrowstyle='->',color=INK,lw=.8))
cc.axvline(62.1,color=AMB,lw=.7,ls=':')
cc.annotate("1er contact à $\\mu$ résiduel\n62,1 µs",xy=(62.1,0),xytext=(0.04,0.78),
    textcoords='axes fraction',fontsize=8,color=AMB,
    arrowprops=dict(arrowstyle='->',color=AMB,lw=.8))
hh,ll=c2.get_legend_handles_labels(); c2.legend(hh,ll,fontsize=8,frameon=False,loc='upper left')
# (d) energies
e=ax[1,1]
for nom,lab,col,ls in [('eFric','frottement',RED,'-'),('eEl','éléments',BLUE,'-'),
                       ('eJnt','joints',AMB,'-'),('eGc','contact (net)',INK,'--'),
                       ('eLys','absorbantes','0.6','-')]:
    e.plot(t,np.abs(h[:,I(nom)]),color=col,ls=ls,lw=1.3,label=lab)
e.plot(t,bdw,color=RED,ls=':',lw=1.6,label='dont pulvérisation')
e.plot(t,0.5*M_TOOL*vb**2,color='0.35',lw=1.0,ls='-.',label="énergie cinétique de l'outil")
e.axhline(0.5*M_TOOL*9.5**2,color='0.8',lw=.7)
e.set_xlabel('temps  [µs]'); e.set_ylabel('énergie cumulée  [J]')
e.set_title("(d)  Canaux d'énergie",loc='left'); e.legend(fontsize=7.5,frameon=False,ncol=2)
for a_ in ax.flat: a_.grid(alpha=.15,lw=.4); a_.set_axisbelow(True)
fig.suptitle(f"{TAG} — pulvérisation × coulomb × contact résiduel  |  "
             f"cycle complet : rebroussement {t[im]:.0f} µs, rebond {vb[-1]:+.2f} m/s "
             f"($e={abs(vb[-1])/9.5:.2f}$), {int(nb[-1])} fissures",fontsize=10.5)
fig.tight_layout(rect=[0,0,1,0.955])
for ext in ('png','pdf'): fig.savefig(f'{SP}/p1/fig/rapport_{TAG}.{ext}',dpi=145)
print(f"F_pic jauge {Fg.max():.1f} kN | F_pic inertie {liss(Fi).max():.1f} kN | "
      f"p_max {p.max():.3f} mm | p_res {p[-1]:.3f} mm | e {abs(vb[-1])/9.5:.3f} | "
      f"nPulv {int(npul[-1])} | bdWork {bdw[-1]:.2f} J")
print(f'{SP}/p1/fig/rapport_{TAG}.png')
