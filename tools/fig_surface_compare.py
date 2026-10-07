"""Compare initial free-surface crack edges with a subsurface section.
Read-only frozen simulation frames. No reconstruction of a post-ejection crater.
"""
from pathlib import Path
import argparse
import json
import hashlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from fig_icl_stanne_review import ROOT, SHIFT, RED, YELLOW, GRAY, facecolors
from imp_lib import read_vtu, broken_mask, frame_times
from fig_joints_cuts import traces


def surface_edges(P, tol=1e-6):
    """P in mm relative to initial surface. Return unique edges on z=0.
    A triangle touching the surface at only one vertex has no line trace.
    Coplanar triangles are counted separately; they are not crack-mouth edges.
    """
    on=np.abs(P[:,:,2]) <= tol
    full=np.flatnonzero(on.sum(axis=1)==3)
    pairs=[]; ids=[]; local=[]; keys=set()
    for i in np.flatnonzero(on.sum(axis=1)==2):
        ab=np.flatnonzero(on[i]); s=P[i,ab]
        key=tuple(sorted(tuple(v) for v in np.rint(s/tol).astype(np.int64)))
        if key in keys: continue
        keys.add(key); pairs.append(s); ids.append(i); local.append(ab)
    return (np.array(pairs).reshape(-1,2,3),np.array(ids,int),
            np.array(local,int).reshape(-1,2),len(full),int(np.sum(on.sum(axis=1)==1)))


def checks():
    p=np.array([[[0,0,0],[1,0,0],[0,1,-1]],
                [[0,0,0],[0,1,-1],[1,1,-1]],
                [[1,0,0],[0,0,0],[0,-1,-1]],
                [[0,0,0],[1,0,0],[0,1,0]]],float)
    s,i,ab,nc,nv=surface_edges(p)
    assert len(s)==1 and nc==1 and nv==1
    assert np.allclose(s[0],p[0,:2])
    assert len(traces(p,2,0)[0])==0  # Negative control: old strict filter misses it.
    assert len(traces(p,2,-.5)[0])==3
    assert len(surface_edges(p,tol=1e-7)[0])==1


def panel(ax,segments,colors,title):
    ax.set_facecolor(GRAY)
    ax.add_collection(LineCollection(segments[:,:,:2],colors=colors,linewidths=.85,
                                    capstyle='round'))
    ax.set(xlim=(-21,21),ylim=(-21,21),xlabel='x [mm]',ylabel='y [mm]')
    ax.set_aspect('equal');ax.set_xticks([-20,-10,0,10,20]);ax.set_yticks([-20,-10,0,10,20])
    ax.tick_params(labelsize=8,length=2,color='#777777')
    for spine in ax.spines.values():spine.set_visible(False)
    ax.set_title(title,loc='left',fontsize=10,pad=9)


def metric(s):
    return dict(segments=len(s),total_trace_length_mm=float(np.linalg.norm(s[:,1]-s[:,0],axis=1).sum()),
                max_radius_mm=float(np.linalg.norm(s[:,:,:2],axis=2).max()) if len(s) else 0)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--frame',type=int,default=11)
    ap.add_argument('--run',default='out_stanne2025_rock137')
    ap.add_argument('--out',default='output/pdf/stanne_surface_165us');a=ap.parse_args()
    checks()
    run=ROOT/a.run;out=ROOT/a.out;out.mkdir(parents=True,exist_ok=True)
    f0=run/'fdem3d_joints_0000.vtu';fk=run/f'fdem3d_joints_{a.frame:04d}.vtu'
    p0,c0,_=read_vtu(f0);p,c,f=read_vtu(fk)
    assert np.array_equal(c,c0)
    mask=broken_mask(f);P=(p0[c[mask]]-SHIFT)*1000;Pc=(p[c[mask]]-SHIFT)*1000
    mode=f['breakMode'][mask]
    s0,ids,local,ncop,nvertex=surface_edges(P)
    assert ncop==0, 'Coplanar facets need an explicit boundary treatment'
    for tol in (1e-7,1e-5):assert len(surface_edges(P,tol)[0])==len(s0)
    current=np.array([Pc[i,ab] for i,ab in zip(ids,local)])
    ss,ii=traces(P,2,-.5)
    colors=facecolors(P,mode)
    t=frame_times(run)[a.frame]*1e6
    fig,axes=plt.subplots(2,2,figsize=(10.7,10.6))
    fig.subplots_adjust(left=.085,right=.965,bottom=.13,top=.855,wspace=.20,hspace=.32)
    fig.suptitle(f'Fissures en surface ou à -0,5 mm ?  |  {t:.0f} µs',x=.085,y=.974,ha='left',
                 fontsize=17,fontweight='bold')
    fig.text(.085,.925,'St Anne, rock137 - même instant, mêmes limites et mêmes épaisseurs de trait',fontsize=10,color='#555555')
    panel(axes[0,0],s0,colors[ids],f'(a) Surface initiale : z = 0\n{len(s0)} arêtes uniques débouchantes')
    panel(axes[0,1],ss,colors[ii],f'(b) Coupe sous la surface : z = -0,5 mm\n{len(ss)} segments d’intersection')
    panel(axes[1,0],s0,'#155fa4','(c) Superposition dans la géométrie initiale\nBleu : surface ; orange : coupe à -0,5 mm')
    axes[1,0].add_collection(LineCollection(ss[:,:,:2],colors='#e27721',linewidths=.65,alpha=.9))
    panel(axes[1,1],current,colors[ids],
          '(d) Mêmes arêtes de surface, positions courantes\nProjection XY ; ce n’est pas une coupe à z = 0')
    fig.legend(handles=[Line2D([0],[0],color=RED,label='Traction (mode exporté)'),
                        Line2D([0],[0],color=YELLOW,label='Cisaillement (mode exporté)')],
               loc='lower center',bbox_to_anchor=(.53,.061),frameon=False,ncol=2,fontsize=9)
    fig.text(.085,.035,'Les arêtes de surface sont extraites par leur appartenance au plan initial, y compris les intersections tangentes.',fontsize=8)
    fig.text(.085,.018,'Aucun lissage. Pas de reconstruction des nouvelles surfaces exposées après éjection. Les couleurs de mode restent indicatives.',fontsize=8)
    for ext in ('png','pdf'):fig.savefig(out/f'comparaison_surface_coupe.{ext}',dpi=220)
    plt.close(fig)
    report=dict(run=str(run),frame=a.frame,time_us=t,broken_facets=int(mask.sum()),
                surface=metric(s0),subsurface=metric(ss),surface_current=metric(current),
                point_only_contacts=nvertex,coplanar_triangles=ncop,
                old_strict_filter_surface_segments=len(traces(P,2,0)[0]),
                tolerance_mm=1e-6,tolerance_sensitivity='same count at 1e-7, 1e-6, 1e-5 mm',
                synthetic_tests='PASS: edge hit, vertex-only hit, duplicate edge, coplanar detection, negative strict filter',
                surface_displacement_mm_max=float(np.linalg.norm(current-s0,axis=2).max()),
                source_sha256={q.name:hashlib.sha256(q.read_bytes()).hexdigest() for q in (f0,fk)})
    (out/'mesures.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(report,indent=2,ensure_ascii=True))


if __name__=='__main__':main()
