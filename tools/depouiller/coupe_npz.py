# -*- coding: utf-8 -*-
"""Coupe verticale exacte d un run a outil maille, depuis un .npz.
   Roche coloree par bulkD (pulverisation), fissures en TRAITS.
   usage: coupe_npz.py [npz] [tag] [frame|-1] [champ bulkD|vonMises]"""
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
CH =sys.argv[4] if len(sys.argv)>4 else 'bulkD'
R_INS=0.00851

d=np.load(NPZ); X=d['el_xyz'][FR]; V=d['el_'+CH][FR]
cx,cy,_=d['axe']; C=float(cy)
zs=float(d['el_xyz'][0][:,:,2].max())   # surface libre = max des NOEUDS
h=d['hist']; c=list(d['hist_cols']); I=lambda n:c.index(n)
t_fr=h[-1,I('t')] if FR in (-1,len(d['el_xyz'])-1) else h[:,I('t')][int(round((FR)/(len(d['el_xyz'])-1)*(len(h)-1)))]

E=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
def coupe_cellules(P):
    """P (n, k, 3) -> liste de polygones (x,z) de l intersection avec y = C."""
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
            out.append(np.stack([pts[o,0],pts[o,2]],1)); idx.append(i)
    return out,np.array(idx,int)

def coupe_triangles(T):
    """T (n,3,3) -> segments (x,z) de l intersection avec y = C."""
    seg=[]; idx=[]
    for i,q in enumerate(T):
        pts=[]
        for a,b in ((0,1),(1,2),(2,0)):
            ya,yb=q[a,1]-C,q[b,1]-C
            if ya*yb<0:
                s=ya/(ya-yb); p=q[a]+s*(q[b]-q[a]); pts.append([p[0],p[2]])
        if len(pts)==2: seg.append(pts); idx.append(i)
    return seg,np.array(idx,int)

poly,ie=coupe_cellules(X)
seg,ij=coupe_triangles(d['jt_xyz'])
mode=d['jt_mode'][ij] if len(ij) else np.array([])
tb=d['jt_tBreak'][ij] if len(ij) else np.array([])
z_tip=float(d['tool_z'][FR,0])

fig,axs=plt.subplots(1,2,figsize=(12.4,5.6))
for ax,demi,titre in ((axs[0],0.060,'(a)  vue large'),(axs[1],0.020,'(b)  zoom sur le noyau')):
    P=[(p-np.array([cx,zs]))*1e3 for p in poly]
    vals=V[ie]
    if CH=='bulkD':
        pc=PolyCollection(P,array=vals,cmap='inferno_r',clim=(0,0.9),
                          edgecolors='0.85',linewidths=.15)
    else:
        pc=PolyCollection(P,array=vals/1e6,cmap='YlGnBu',
                          clim=(0,np.percentile(vals,98)/1e6),edgecolors='none')
    ax.add_collection(pc)
    if len(seg):
        S=[[( (p[0]-cx)*1e3,(p[1]-zs)*1e3) for p in s] for s in seg]
        ax.add_collection(LineCollection(S,colors=np.where(mode==2,AMB,RED),
                                         linewidths=1.0,alpha=.9,zorder=3))
    th=np.linspace(np.pi,2*np.pi,120)                       # hemisphere de l insert
    zc=(z_tip+R_INS-zs)*1e3
    ax.plot(R_INS*1e3*np.cos(th),zc+R_INS*1e3*np.sin(th),color=INK,lw=1.4,zorder=4)
    ax.axhline(0,color=INK,lw=.6,ls=':')
    ax.set_xlim(-demi*1e3,demi*1e3); ax.set_ylim(-demi*1e3*0.9,demi*1e3*0.35)
    ax.set_aspect('equal'); ax.set_xlabel('$x$  [mm]'); ax.set_ylabel('$z$ sous la surface  [mm]')
    ax.set_title(titre,loc='left')
cb=fig.colorbar(pc,ax=axs,fraction=.028,pad=.02)
cb.set_label("endommagement de pulvérisation  $D$" if CH=='bulkD' else 'von Mises [MPa]',fontsize=9)
axs[1].plot([],[],color=RED,lw=1.2,label='fissure en traction')
axs[1].plot([],[],color=AMB,lw=1.2,label='fissure en cisaillement')
axs[1].legend(fontsize=8,frameon=False,loc='lower right')
npul=int((V>=0.9).sum())
fig.suptitle(f"{TAG} — coupe $y = y_0$ à $t$ = {t_fr*1e6:.0f} µs   |   "
             f"{len(seg)} traces de fissures, {npul} éléments à $D=D_{{max}}$ dans le volume, "
             f"$D_{{max}}$ atteint = {V.max():.2f}",fontsize=10.5)
for ext in ('png','pdf'): fig.savefig(f'{SP}/p1/fig/coupe_{TAG}.{ext}',dpi=145,bbox_inches='tight')
prof=(zs-X[V>=0.5][:,:,2].min())*1e3 if (V>=0.5).any() else 0
ray=np.abs(X[V>=0.5][:,:,0]-cx).max()*1e3 if (V>=0.5).any() else 0
print(f"cellules coupees {len(poly)} | traces {len(seg)} ({int((mode==2).sum())} cisaillement) | "
      f"D max {V.max():.3f} | zone D>0.5 : profondeur {prof:.1f} mm, rayon {ray:.1f} mm | "
      f"pointe outil {(zs-z_tip)*1e3:+.2f} mm sous la surface")
print(f'{SP}/p1/fig/coupe_{TAG}.png')
