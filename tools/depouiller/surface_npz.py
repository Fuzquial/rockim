# -*- coding: utf-8 -*-
"""Facies de SURFACE d un run a outil maille, depuis un .npz.
   Chaque joint rompu qui DEBOUCHE possede une arete dans le plan de la
   surface libre : c est elle qu on trace. Rendu en TRAITS, comme une photo
   de la surface apres essai (et non des facettes projetees, illisibles).
   usage: surface_npz.py [npz] [tag] [frame|-1] [demi-largeur mm]"""
import sys, glob, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm, matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
for f in glob.glob('/usr/share/texmf/fonts/opentype/public/lm/lmroman10-*.otf'):
    fm.fontManager.addfont(f)
plt.rcParams.update({'font.family':'serif','font.serif':['Latin Modern Roman'],
                     'mathtext.fontset':'cm','font.size':9,'axes.linewidth':0.7})
RED,AMB,BLUE,INK='#c1121f','#b5730a','#1f4e9c','#3a3a3a'
SP='/tmp/claude-0/-home-user/0263e025-847d-592e-aea7-7fec643bb1d6/scratchpad'
NPZ=sys.argv[1] if len(sys.argv)>1 else '/home/user/rockim/bench_impact/donnees/P1.npz'
TAG=sys.argv[2] if len(sys.argv)>2 else 'P1'
FR =int(sys.argv[3]) if len(sys.argv)>3 else -1
DEMI=float(sys.argv[4]) if len(sys.argv)>4 else 30.0
R_INS=0.00851

d=np.load(NPZ); XX=d['el_xyz']; BD=d['el_bulkD']
cx,cy,_=d['axe']
h=d['hist']; c=list(d['hist_cols']); I=lambda n:c.index(n)
N=len(XX); T=float(h[-1,I('t')]); tf=np.linspace(0,T,N); t_fr=tf[FR]
JT=d['jt_xyz']; JM=d['jt_mode']; JTB=d['jt_tBreak']

def apparier(J,nodes,cell=1.5e-3):
    V=J.reshape(-1,3); K=(V/cell).astype(np.int64); KN=(nodes/cell).astype(np.int64)
    tab={}
    for i,k in enumerate(map(tuple,KN)): tab.setdefault(k,[]).append(i)
    best=np.full(len(V),-1,np.int64); dist=np.full(len(V),9e9)
    for j in range(len(V)):
        a0,b0,c0=K[j]; cand=[]
        for a in(-1,0,1):
            for b in(-1,0,1):
                for e in(-1,0,1): cand.extend(tab.get((a0+a,b0+b,c0+e),()))
        if not cand: continue
        cand=np.array(cand); dd=((nodes[cand]-V[j])**2).sum(1); m=dd.argmin()
        best[j]=cand[m]; dist[j]=dd[m]**0.5
    tri=best.reshape(-1,3); dist=dist.reshape(-1,3)
    return tri,(tri>=0).all(1)&(dist<5e-4).all(1)

TRI,OKJ=apparier(JT,XX[-1].reshape(-1,3))
n0=XX[0].reshape(-1,3); zs=float(n0[:,2].max())
# ARETE DEBOUCHANTE : les 2 sommets du triangle qui sont dans le plan de la
# surface libre, jauges dans la configuration de REFERENCE (frame 0).
zref=n0[TRI][:,:,2]; surf=zref>=zs-1e-5
deb=OKJ&(surf.sum(1)==2)
idx2=np.array([TRI[i][surf[i]] for i in np.nonzero(deb)[0]])   # (n, 2)
MOD=JM[deb]; TB=JTB[deb]
nod=XX[FR].reshape(-1,3)
SEG=(nod[idx2][:,:,:2]-np.array([cx,cy]))*1e3                  # (n, 2, 2) en mm
vu=TB<=t_fr
LON=np.linalg.norm(nod[idx2][:,0]-nod[idx2][:,1],axis=1)*1e3
MID=SEG.mean(1); RAY=np.hypot(MID[:,0],MID[:,1]); AZ=np.arctan2(MID[:,1],MID[:,0])

fig=plt.figure(figsize=(11.0,10.4))
ax_a=fig.add_subplot(2,2,1); ax_b=fig.add_subplot(2,2,2,projection='polar')
ax_c=fig.add_subplot(2,2,3); ax_d=fig.add_subplot(2,2,4)
th=np.linspace(0,2*np.pi,200)
def deco(a,titre):
    a.plot(R_INS*1e3*np.cos(th),R_INS*1e3*np.sin(th),color=BLUE,lw=1.0)
    a.set_xlim(-DEMI,DEMI); a.set_ylim(-DEMI,DEMI); a.set_aspect('equal')
    a.set_xlabel('$x$  [mm]'); a.set_ylabel('$y$  [mm]'); a.set_title(titre,loc='left')
    a.grid(alpha=.15,lw=.4); a.set_axisbelow(True)

# (a) LE FACIES : traces debouchantes, en traits
ax_a.add_collection(LineCollection(list(SEG[vu]),colors=np.where(MOD[vu]==2,AMB,RED),
                                   linewidths=1.4,alpha=.9))
ax_a.plot([],[],color=RED,lw=1.4,label=f'traction {int((MOD[vu]!=2).sum())}')
ax_a.plot([],[],color=AMB,lw=1.4,label=f'cisaillement {int((MOD[vu]==2).sum())}')
ax_a.legend(fontsize=7.5,frameon=False,loc='upper right')
deco(ax_a,'(a)  Faciès de surface — arêtes débouchantes')

# (b) densite azimutale des traces, PONDEREE PAR LEUR LONGUEUR
nsec=36; bins=np.linspace(-np.pi,np.pi,nsec+1)
dens,_=np.histogram(AZ[vu],bins=bins,weights=LON[vu])
ctr=(bins[:-1]+bins[1:])/2
ax_b.bar(ctr,dens,width=2*np.pi/nsec,color=RED,alpha=.55,edgecolor=INK,linewidth=.4)
moy=dens.mean()
ax_b.plot(np.linspace(-np.pi,np.pi,200),np.full(200,moy),color=INK,lw=.9,ls='--')
ax_b.set_theta_zero_location('E'); ax_b.set_title('(b)  Longueur de trace par secteur  [mm]',loc='left')
ax_b.tick_params(labelsize=7)
contr=dens.max()/moy if moy>0 else 0
ax_b.text(0.5,-0.10,f"contraste max/moyenne = {contr:.2f}   "
                    f"(isotrope $\\to$ 1)",transform=ax_b.transAxes,ha='center',fontsize=8,color=INK)

# (c) pulverisation
V=BD[FR]; m=V>0.05
if m.any():
    Pp=(XX[FR][m][:,:,:2]-np.array([cx,cy]))*1e3; o=np.argsort(V[m])
    pc=PolyCollection(Pp[o],array=V[m][o],cmap='inferno_r',clim=(0,0.9),edgecolors='none')
    ax_c.add_collection(pc)
    cb=fig.colorbar(pc,ax=ax_c,fraction=.046,pad=.03); cb.set_label('$D$',fontsize=8)
    cb.ax.tick_params(labelsize=7)
deco(ax_c,f'(c)  Pulvérisation  ({int((V>=0.9).sum())} éléments à $D_{{max}}$)')

# (d) cratere
sf=XX[0][:,:,2].max(axis=1)>zs-3e-4
dz=(XX[FR][sf][:,:,2].mean(1)-XX[0][sf][:,:,2].mean(1))*1e3
Ps=(XX[FR][sf][:,:,:2]-np.array([cx,cy]))*1e3
enpl=np.abs(dz)<3.0; bour=dz[enpl].max() if enpl.any() else 0.0
lim=max(-dz.min(),bour,0.5); o=np.argsort(-np.abs(dz))
pcd=PolyCollection(Ps[o],array=np.clip(dz[o],-lim,lim),cmap='RdBu_r',clim=(-lim,lim),
                   edgecolors='none')
ax_d.add_collection(pcd)
cb=fig.colorbar(pcd,ax=ax_d,fraction=.046,pad=.03)
cb.set_label('déplacement vertical  [mm]\n(éjecta saturés)',fontsize=8); cb.ax.tick_params(labelsize=7)
deco(ax_d,f'(d)  Cratère — creux {dz.min():+.2f} mm, bourrelet {bour:+.2f} mm')
ax_d.text(0.02,0.03,f"{int((~enpl).sum())} éclats en vol (jusqu'à {dz.max():+.1f} mm)",
          transform=ax_d.transAxes,fontsize=7.5,color=INK)

fig.suptitle(f"{TAG} — surface libre à $t$ = {t_fr*1e6:.0f} µs   |   "
             f"{int(vu.sum())} fissures débouchantes sur {int(OKJ.sum())} "
             f"({LON[vu].sum():.0f} mm de trace cumulée)",fontsize=10.5)
fig.tight_layout(rect=[0,0,1,0.96])
for ext in ('png','pdf'): fig.savefig(f'{SP}/p1/fig/surface_{TAG}.{ext}',dpi=145)
print(f"debouchantes {int(vu.sum())}/{int(OKJ.sum())} ({int((MOD[vu]==2).sum())} cisaillement) | "
      f"trace cumulee {LON[vu].sum():.0f} mm | rayon max {RAY[vu].max():.1f} mm, "
      f"p95 {np.percentile(RAY[vu],95):.1f} mm | contraste azimutal {contr:.2f}")
print(f'{SP}/p1/fig/surface_{TAG}.png')
