# -*- coding: utf-8 -*-
"""Film de la coupe d un run a outil maille, depuis un .npz.
   Roche coloree par bulkD (echelle FIXE), fissures en traits CUMULEES.
   usage: film_npz.py [npz] [tag] [demi-largeur mm] [ms par image]"""
import sys, os, glob, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm, matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
from PIL import Image
for f in glob.glob('/usr/share/texmf/fonts/opentype/public/lm/lmroman10-*.otf'):
    fm.fontManager.addfont(f)
plt.rcParams.update({'font.family':'serif','font.serif':['Latin Modern Roman'],
                     'mathtext.fontset':'cm','font.size':9,'axes.linewidth':0.7})
RED,AMB,INK='#c1121f','#b5730a','#3a3a3a'
SP='/tmp/claude-0/-home-user/0263e025-847d-592e-aea7-7fec643bb1d6/scratchpad'
NPZ=sys.argv[1] if len(sys.argv)>1 else '/home/user/rockim/bench_impact/donnees/P1.npz'
TAG=sys.argv[2] if len(sys.argv)>2 else 'P1'
DEMI=float(sys.argv[3]) if len(sys.argv)>3 else 24.0     # mm
MSPF=int(sys.argv[4]) if len(sys.argv)>4 else 260        # ms par image
R_INS=0.00851

d=np.load(NPZ); XX=d['el_xyz']; BD=d['el_bulkD']
cx,cy,_=d['axe']; C=float(cy); zs=float(XX[0][:,:,2].max())
h=d['hist']; c=list(d['hist_cols']); I=lambda n:c.index(n)
N=len(XX); T=float(h[-1,I('t')])
tf=np.linspace(0.0,T,N)                                   # frames uniformes (verifie)
nb=np.interp(tf,h[:,I('t')],h[:,I('nBroken')])
npu=np.interp(tf,h[:,I('t')],h[:,I('nPulv')])
JT=d['jt_xyz']; JTB=d['jt_tBreak']; JM=d['jt_mode']
E=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]

def coupe_cellules(P):
    out=[]; idx=[]
    for i,q in enumerate(P):
        pts=[]
        for a,b in E:
            ya,yb=q[a,1]-C,q[b,1]-C
            if ya*yb<0:
                s=ya/(ya-yb); pts.append(q[a]+s*(q[b]-q[a]))
        if len(pts)>=3:
            pts=np.array(pts); m=pts.mean(axis=0)
            o=np.argsort(np.arctan2(pts[:,2]-m[2],pts[:,0]-m[0]))
            out.append(np.stack([(pts[o,0]-cx)*1e3,(pts[o,2]-zs)*1e3],1)); idx.append(i)
    return out,np.array(idx,int)

# LES FISSURES SUIVENT LA MATIERE. Le .npz ne garde les joints qu a la
# DERNIERE frame (tBreak y porte toute la chronologie) : les dessiner a cette
# position sur une roche qui bouge fait flotter des traits dans le vide des
# que la matiere est ejectee. Reparation : le maillage est DISCONTINU (4
# noeuds propres par tetraedre), et les sommets des triangles de joint sont
# EXACTEMENT des noeuds d elements — apparies ici au nanometre pres. On
# retrouve donc la geometrie de chaque fissure a n importe quelle frame.
def apparier(JT, nodes, cell=1.5e-3):
    V=JT.reshape(-1,3); K=(V/cell).astype(np.int64); KN=(nodes/cell).astype(np.int64)
    tab={}
    for i,k in enumerate(map(tuple,KN)): tab.setdefault(k,[]).append(i)
    best=np.full(len(V),-1,np.int64); dist=np.full(len(V),9e9)
    for j in range(len(V)):
        a0,b0,c0=K[j]; cand=[]
        for a in (-1,0,1):
            for b in (-1,0,1):
                for c in (-1,0,1): cand.extend(tab.get((a0+a,b0+b,c0+c),()))
        if not cand: continue
        cand=np.array(cand); dd=((nodes[cand]-V[j])**2).sum(1); m=dd.argmin()
        best[j]=cand[m]; dist[j]=dd[m]**0.5
    tri=best.reshape(-1,3); dist=dist.reshape(-1,3)
    ok=(tri>=0).all(1) & (dist<5e-4).all(1)
    return tri, ok

TRI, OKJ = apparier(JT, XX[-1].reshape(-1,3))
JM_ok, JTB_ok = JM[OKJ], JTB[OKJ]
print(f"fissures suivies : {int(OKJ.sum())} / {len(JT)} "
      f"(erreur de reconstruction a la derniere frame "
      f"{np.abs(XX[-1].reshape(-1,3)[TRI[OKJ]] - JT[OKJ]).max()*1e3:.4f} mm)")

def coupe_fissures(P3):
    """P3 (n,3,3) -> (segments, indices) de l intersection avec y = C. Vectorise."""
    y=P3[:,:,1]-C
    pts=np.full((len(P3),3,2),np.nan)
    for e,(a,b) in enumerate(((0,1),(1,2),(2,0))):
        m=(y[:,a]*y[:,b])<0
        if m.any():
            s=(y[m,a]/(y[m,a]-y[m,b]))[:,None]
            q=P3[m,a,:]+s*(P3[m,b,:]-P3[m,a,:])
            pts[m,e,0]=(q[:,0]-cx)*1e3; pts[m,e,1]=(q[:,2]-zs)*1e3
    n=(~np.isnan(pts[:,:,0])).sum(1)
    idx=np.nonzero(n==2)[0]
    seg=[pts[i][~np.isnan(pts[i,:,0])] for i in idx]
    return seg, idx

os.makedirs(f'{SP}/p1/frames',exist_ok=True)
imgs=[]
for k in range(N):
    poly,ie=coupe_cellules(XX[k])
    fig,ax=plt.subplots(figsize=(6.6,5.2),dpi=140)
    pc=PolyCollection(poly,array=BD[k][ie],cmap='inferno_r',clim=(0,0.9),
                      edgecolors='0.88',linewidths=.12)
    ax.add_collection(pc)
    vu=JTB_ok<=tf[k]
    if vu.any():
        P3=XX[k].reshape(-1,3)[TRI[OKJ][vu]]          # geometrie A CETTE FRAME
        seg,idx=coupe_fissures(P3)
        if len(seg):
            md=JM_ok[vu][idx]
            ax.add_collection(LineCollection(seg,colors=np.where(md==2,AMB,RED),
                                             linewidths=1.1,alpha=.92,zorder=3))
    th=np.linspace(np.pi,2*np.pi,120); zc=(float(d['tool_z'][k,0])+R_INS-zs)*1e3
    ax.plot(R_INS*1e3*np.cos(th),zc+R_INS*1e3*np.sin(th),color=INK,lw=1.5,zorder=4)
    ax.axhline(0,color=INK,lw=.5,ls=':')
    ax.set_xlim(-DEMI,DEMI); ax.set_ylim(-DEMI*0.92,DEMI*0.34); ax.set_aspect('equal')
    ax.set_xlabel('$x$  [mm]'); ax.set_ylabel('$z$ sous la surface  [mm]')
    ax.set_title(f"{TAG} — coupe $y=y_0$",loc='left',fontsize=10)
    ax.text(0.98,0.96,f"$t$ = {tf[k]*1e6:5.0f} µs",transform=ax.transAxes,ha='right',
            va='top',fontsize=11,color=INK)
    ax.text(0.98,0.89,f"{int(nb[k])} fissures\n{int(npu[k])} éléments pulvérisés",
            transform=ax.transAxes,ha='right',va='top',fontsize=8,color=INK)
    ax.add_patch(plt.Rectangle((0.06,0.03),0.88,0.012,transform=ax.transAxes,
                               facecolor='0.88',edgecolor='none',zorder=5))
    ax.add_patch(plt.Rectangle((0.06,0.03),0.88*(k/(N-1)),0.012,transform=ax.transAxes,
                               facecolor=RED,edgecolor='none',zorder=6))
    cb=fig.colorbar(pc,ax=ax,fraction=.042,pad=.02); cb.set_label("$D$ (pulvérisation)",fontsize=8)
    cb.ax.tick_params(labelsize=7)
    p=f'{SP}/p1/frames/f{k:03d}.png'; fig.savefig(p,bbox_inches='tight'); plt.close(fig)
    imgs.append(Image.open(p).convert('P',palette=Image.ADAPTIVE,colors=192))
gif=f'{SP}/p1/fig/film_{TAG}.gif'
imgs[0].save(gif,save_all=True,append_images=imgs[1:],duration=MSPF,loop=0,optimize=True)
print(f"{N} images | {int(OKJ.sum())} fissures suivies | {gif} ({os.path.getsize(gif)/1e6:.1f} Mo)")
