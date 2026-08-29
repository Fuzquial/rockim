# -*- coding: utf-8 -*-
"""Film VU DE DESSUS : le facies de surface qui se construit. Chaque joint
   rompu DEBOUCHANT est trace par son arete dans le plan de la surface libre
   (traits), le noyau pulverise est peint par-dessus (rendu elements).
   usage: filmhaut_npz.py [npz] [tag] [demi mm] [ms/image]"""
import sys, os, glob, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm, matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
from PIL import Image
for f in glob.glob('/usr/share/texmf/fonts/opentype/public/lm/lmroman10-*.otf'):
    fm.fontManager.addfont(f)
plt.rcParams.update({'font.family':'serif','font.serif':['Latin Modern Roman'],
                     'mathtext.fontset':'cm','font.size':9,'axes.linewidth':0.7})
RED,AMB,BLUE,INK='#c1121f','#b5730a','#1f4e9c','#3a3a3a'
SP='/tmp/claude-0/-home-user/0263e025-847d-592e-aea7-7fec643bb1d6/scratchpad'
NPZ=sys.argv[1] if len(sys.argv)>1 else '/home/user/rockim/bench_impact/donnees/P1.npz'
TAG=sys.argv[2] if len(sys.argv)>2 else 'P1'
DEMI=float(sys.argv[3]) if len(sys.argv)>3 else 28.0
MSPF=int(sys.argv[4]) if len(sys.argv)>4 else 260
R_INS=0.00851

d=np.load(NPZ); XX=d['el_xyz']; BD=d['el_bulkD']
cx,cy,_=d['axe']
h=d['hist']; c=list(d['hist_cols']); I=lambda n:c.index(n)
N=len(XX); T=float(h[-1,I('t')]); tf=np.linspace(0,T,N)
JT=d['jt_xyz']; JM=d['jt_mode']; JTB=d['jt_tBreak']
nb=np.interp(tf,h[:,I('t')],h[:,I('nBroken')])

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
surf=n0[TRI][:,:,2]>=zs-1e-5
deb=OKJ&(surf.sum(1)==2)
idx2=np.array([TRI[i][surf[i]] for i in np.nonzero(deb)[0]])
MOD=JM[deb]; TB=JTB[deb]
print(f"aretes debouchantes : {len(idx2)} sur {int(OKJ.sum())} fissures "
      f"({int((MOD==2).sum())} en cisaillement)")

os.makedirs(f'{SP}/p1/frames_h',exist_ok=True); imgs=[]
th=np.linspace(0,2*np.pi,200)
for k in range(N):
    nod=XX[k].reshape(-1,3)
    fig,ax=plt.subplots(figsize=(6.2,6.0),dpi=140)
    vu=TB<=tf[k]
    if vu.any():
        S=(nod[idx2[vu]][:,:,:2]-np.array([cx,cy]))*1e3
        ax.add_collection(LineCollection(list(S),colors=np.where(MOD[vu]==2,AMB,RED),
                                         linewidths=1.4,alpha=.9,zorder=3))
    V=BD[k]; m=V>0.15
    if m.any():
        Pp=(XX[k][m][:,:,:2]-np.array([cx,cy]))*1e3; o=np.argsort(V[m])
        ax.add_collection(PolyCollection(Pp[o],array=V[m][o],cmap='inferno_r',
                                         clim=(0,0.9),edgecolors='none',alpha=.9,zorder=2))
    ax.plot(R_INS*1e3*np.cos(th),R_INS*1e3*np.sin(th),color=BLUE,lw=1.1,zorder=6)
    ax.set_xlim(-DEMI,DEMI); ax.set_ylim(-DEMI,DEMI); ax.set_aspect('equal')
    ax.set_xlabel('$x$  [mm]'); ax.set_ylabel('$y$  [mm]')
    ax.set_title(f"{TAG} — faciès de surface",loc='left',fontsize=10)
    ax.text(0.97,0.97,f"$t$ = {tf[k]*1e6:5.0f} µs",transform=ax.transAxes,ha='right',
            va='top',fontsize=11,color=INK)
    ax.text(0.97,0.905,f"{int(vu.sum())} traces débouchantes\n"
                       f"{int(nb[k])} fissures au total\n"
                       f"{int((V>=0.9).sum())} éléments pulvérisés",
            transform=ax.transAxes,ha='right',va='top',fontsize=8,color=INK)
    ax.add_patch(plt.Rectangle((0.06,0.025),0.88,0.011,transform=ax.transAxes,
                               facecolor='0.88',edgecolor='none',zorder=7))
    ax.add_patch(plt.Rectangle((0.06,0.025),0.88*k/(N-1),0.011,transform=ax.transAxes,
                               facecolor=RED,edgecolor='none',zorder=8))
    ax.grid(alpha=.13,lw=.4); ax.set_axisbelow(True)
    p=f'{SP}/p1/frames_h/f{k:03d}.png'; fig.savefig(p,bbox_inches='tight'); plt.close(fig)
    imgs.append(Image.open(p).convert('P',palette=Image.ADAPTIVE,colors=192))
gif=f'{SP}/p1/fig/filmhaut_{TAG}.gif'
imgs[0].save(gif,save_all=True,append_images=imgs[1:],duration=MSPF,loop=0,optimize=True)
print(f"{N} images | {gif} ({os.path.getsize(gif)/1e6:.1f} Mo)")
