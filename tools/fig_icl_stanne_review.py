"""Scientific postprocessing only: opaque fracture surfaces and exact sections.

Frozen snapshot, no edits to solver, input decks or live outputs.
Run from repository root with the bundled Python runtime.
"""
from pathlib import Path
import sys
import json
import hashlib
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bench_impact/tools'))
sys.path.insert(0, str(ROOT / 'tools'))
from imp_lib import read_vtu, broken_mask, frame_times
from crack_paths import unify_vertices, analyse, facet_adjacency_components
from fig_joints_cuts import traces

RED = '#d52329'
YELLOW = '#f6bd18'
GRAY = '#cfcec6'
SHIFT = np.array([.125, .125, .15])
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                     'axes.titlesize': 10, 'pdf.fonttype': 42,
                     'svg.fonttype': 'none', 'savefig.facecolor': 'white'})


def facecolors(P, mode, shade=False):
    colors = np.array([matplotlib.colors.to_rgba(RED if m < 1.5 else YELLOW) for m in mode])
    if shade and len(P):
        n = np.cross(P[:, 1]-P[:, 0], P[:, 2]-P[:, 0])
        n /= np.maximum(np.linalg.norm(n, axis=1)[:, None], 1e-20)
        light = np.array([-.35, -.45, .82]); light /= np.linalg.norm(light)
        colors[:, :3] *= (.56+.44*np.abs(n@light))[:, None]
    return colors


def surface(ax, P, mode, elev=90, azim=-90, limit=23, wire=False):
    col = facecolors(P, mode, shade=not wire)
    mesh = Poly3DCollection(P, facecolors='none' if wire else col,
                           edgecolors=col if wire else 'none',
                           linewidths=.30 if wire else 0, antialiaseds=False,
                           zsort='average')
    ax.add_collection3d(mesh)
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit), zlim=(-24, 1))
    ax.set_box_aspect((2*limit, 2*limit, 25))
    ax.set_proj_type('ortho'); ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    # Match the 3D top row spatial scale to the 2D maps (mplot3d adds margins).
    ax.set_position(ax.get_position())


def plane(ax, P, mode, axis=2, value=-.5, lim=23, bg=True, width=.75):
    S, ids = traces(P, axis, value)
    xy = (0, 1) if axis == 2 else ((0, 2) if axis == 1 else (1, 2))
    ax.set_facecolor(GRAY if bg else 'white')
    if len(S):
        ax.add_collection(LineCollection(S[:, :, xy], colors=facecolors(P, mode)[ids],
                                        linewidths=width, capstyle='round'))
    ax.set_xlim(-lim, lim)
    ax.set_ylim((-lim, lim) if axis == 2 else (-24, 2))
    ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    if axis != 2:
        ax.axhline(0, color='#77756e', lw=.65)
    return S, ids


def scale(ax, vertical=False):
    x, y = -20, -21 if not vertical else -21
    ax.plot([x, x+10], [y, y], color='#252525', lw=1.6)
    ax.text(x+5, y+1, '10 mm', ha='center', va='bottom', fontsize=8)


def legend(fig, y=.03):
    fig.legend(handles=[Patch(color=RED, label='Traction (mode exporté)'),
                        Patch(color=YELLOW, label='Cisaillement (mode exporté)')],
               loc='lower center', bbox_to_anchor=(.5,y), ncol=2, frameon=False,
               fontsize=9)


def save(fig, out, name):
    fig.savefig(out / (name+'.png'), dpi=220)
    fig.savefig(out / (name+'.pdf'))
    plt.close(fig)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--run', default='out_stanne2025_rock137')
    ap.add_argument('--frame', type=int, default=11)
    ap.add_argument('--out', default='output/pdf/stanne_ICL_165us')
    a=ap.parse_args()
    run=ROOT/a.run; out=ROOT/a.out; out.mkdir(parents=True, exist_ok=True)
    times=frame_times(run)
    p0, con0, _=read_vtu(run/'fdem3d_joints_0000.vtu')
    uid=unify_vertices(p0)
    ks=[k for k in (6,8,9,10,a.frame) if k in times]
    data=[]; hashes={}
    for k in dict.fromkeys(ks):
        path=run/f'fdem3d_joints_{k:04d}.vtu'
        p,c,f=read_vtu(path)
        assert np.array_equal(c,con0), 'Topology changed'
        assert 'tBreak' in f
        sel=broken_mask(f)
        # All broken non-bonded joints; no morphology or mode filtering.
        P0=p0[c[sel]]
        assert np.max(P0[:,:,2]) <= .15001, 'Unexpected non-rock broken facet'
        D=dict(k=k,t=times[k], P=(P0-SHIFT)*1000,
               current=(p[c[sel]]-SHIFT)*1000,
               tri=uid[c[sel]], mode=f['breakMode'][sel], n=int(sel.sum()))
        data.append(D)
        hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    D=data[-1]; P=D['P']; mode=D['mode']; t=D['t']*1e6

    # Figure 1: explicitly identified operations, common framing across time.
    fig=plt.figure(figsize=(13.2,8.0))
    gs=fig.add_gridspec(3,len(data), left=.105,right=.985,bottom=.14,top=.87,
                        wspace=.055,hspace=.12,height_ratios=[1.0,1.05,.67])
    fig.suptitle('St Anne - développement du réseau de rupture', x=.105, ha='left',
                 fontsize=17, fontweight='bold', y=.965)
    fig.text(.105,.916,'Impact à 10,66 m/s  |  rock137  |  géométrie initiale  |  facettes rompues uniquement',
             fontsize=10,color='#555555')
    for j,d in enumerate(data):
        ax=fig.add_subplot(gs[0,j],projection='3d'); surface(ax,d['P'],d['mode'])
        ax.set_title(f"{d['t']*1e6:.0f} µs\n{d['n']:,} facettes".replace(',',' '),pad=1)
        ax=fig.add_subplot(gs[1,j]); plane(ax,d['P'],d['mode'])
        if j==0: scale(ax)
        ax=fig.add_subplot(gs[2,j]); plane(ax,d['P'],d['mode'],axis=1,value=0,bg=False)
        if j==0: scale(ax,True)
    for y,label in [(.755,'Surfaces 3D\nvues de dessus'),(.49,'Coupe horizontale\nz = -0,5 mm'),(.262,'Coupe verticale\ny = 0 mm')]:
        fig.text(.015,y,label,fontsize=9,va='center')
    legend(fig,.055)
    fig.text(.105,.026,'Présentation inspirée de Yang et al. (2025), fig. 14. Coupes exactes, sans lissage ; aucune valeur Kuru superposée.',
             fontsize=8,color='#555555')
    save(fig,out,'01_evolution_ICL')

    # Figure 2: exactly the same triangles, operations compared transparently.
    fig=plt.figure(figsize=(12.5,7.5))
    gs=fig.add_gridspec(2,3,left=.05,right=.975,top=.85,bottom=.14,hspace=.26,wspace=.12)
    fig.suptitle(f'Mêmes données, impressions différentes - {t:.0f} µs',x=.05,ha='left',
                 fontsize=17,fontweight='bold',y=.965)
    fig.text(.05,.91,f"{D['n']:,} facettes rompues ; aucun changement de loi, aucun retrait de fissures".replace(',',' '),
             fontsize=10,color='#555555')
    ax=fig.add_subplot(gs[0,0]);
    # Wireframe projected in xy: includes occluded and deep edges.
    ax.add_collection(PolyCollection(P[:,:,:2],facecolors='none',edgecolors=facecolors(P,mode),linewidths=.38))
    ax.set(xlim=(-23,23),ylim=(-23,23)); ax.set_aspect('equal'); ax.set_axis_off()
    ax.set_title('(a) Toutes les arêtes projetées\nToutes les profondeurs superposées')
    ax=fig.add_subplot(gs[0,1],projection='3d');surface(ax,P,mode)
    ax.set_title('(b) Surfaces opaques\nFacettes cachées occultées')
    ax=fig.add_subplot(gs[0,2]);plane(ax,P,mode)
    ax.set_title('(c) Intersection avec z = -0,5 mm\nTraces de fissures près de la surface');scale(ax)
    for j,(axis,value,title) in enumerate([(1,0,'(d) Coupe exacte y = 0 mm'),
                                         (0,0,'(e) Coupe exacte x = 0 mm'),
                                         (2,-3,'(f) Coupe exacte z = -3 mm')]):
        ax=fig.add_subplot(gs[1,j]);plane(ax,P,mode,axis,value,bg=axis==2)
        ax.set_title(title);scale(ax,axis!=2)
    legend(fig,.048)
    fig.text(.05,.018,'Une tranche projetée de 10 mm mélange les fissures ; une coupe exacte ne montre que leurs intersections avec un plan.',
             fontsize=8,color='#555555')
    save(fig,out,'02_effet_du_rendu')

    # Figure 3: geometric connectivity, all background facets retained.
    res=analyse(P/1000+SHIFT,D['tri'])
    rad=[r for r in res['arms'] if r['orient']=='radiale' and r['n']>=10]
    rad.sort(key=lambda r:-r['rmax']); chosen=rad[:3]
    fig=plt.figure(figsize=(12,7.7))
    gs=fig.add_gridspec(1,2,left=.04,right=.98,top=.86,bottom=.26,wspace=.12)
    fig.suptitle(f'Branches reliées au noyau - {t:.0f} µs',x=.05,ha='left',
                 fontsize=17,fontweight='bold',y=.96)
    fig.text(.05,.91,'Continuité par arête partagée dans le maillage initial ; noyau défini par r < 6 mm',
             fontsize=10,color='#555555')
    ax=fig.add_subplot(gs[0,0],projection='3d');surface(ax,P,mode,elev=23,azim=-63)
    ax.set_title('Réseau complet : surfaces de rupture opaques',pad=4)
    ax=fig.add_subplot(gs[0,1]);
    ax.add_collection(PolyCollection(P[:,:,:2],facecolors='#e4e3df',edgecolors='none'))
    palette=['#bd2140','#315d9f','#168675']
    stats=[]
    for i,(r,col) in enumerate(zip(chosen,palette)):
        idx=res['arm_lab']==r['id']; pp=P[idx]
        ax.add_collection(PolyCollection(pp[:,:,:2],facecolors=col,edgecolors='none'))
        rv=np.linalg.norm(pp[:,:,:2],axis=2); tip=pp.reshape(-1,3)[np.argmax(rv)]
        az=np.radians(r['az']); labelxy=np.array([19*np.cos(az),19*np.sin(az)])
        ax.annotate(chr(65+i),xy=tip[:2],xytext=labelxy,ha='center',va='center',fontsize=12,
                    fontweight='bold',color=col,arrowprops=dict(arrowstyle='-',color=col,lw=.9),
                    bbox=dict(facecolor='white',edgecolor='none',pad=2))
        stats.append(dict(label=chr(65+i),**r))
    ax.set(xlim=(-23,23),ylim=(-23,23));ax.set_aspect('equal');ax.set_axis_off();scale(ax)
    ax.set_title('Trois grandes branches ; reste du réseau en gris')
    for i,(r,col) in enumerate(zip(chosen,palette)):
        fig.text(.075,.205-i*.034,
                 f"{chr(65+i)}    {r['n']} facettes connectées    |    centre-pointe : {r['rmax']*1000:.2f} mm    |    azimut : {r['az']:.0f}°",
                 color=col,fontsize=10)
    fig.text(.075,.068,'Les couleurs A/B/C identifient des branches géométriques, pas des modes de rupture. Le noyau n’est pas un cratère mesuré.',
             fontsize=8,color='#555555')
    fig.text(.075,.039,'La continuité du réseau est établie à cette résolution ; la convergence et la comparaison quantitative à Yang restent à vérifier.',
             fontsize=8,color='#555555')
    save(fig,out,'03_branches_connectees')

    sensitivity=[]
    for rcore in (6,8,10):
        rr=analyse(P/1000+SHIFT,D['tri'],rcore=rcore/1000)
        rrads=sorted([r for r in rr['arms'] if r['orient']=='radiale' and r['n']>=10],key=lambda r:-r['rmax'])
        sensitivity.append(dict(rcore_mm=rcore,branches=rrads[:5]))
    meta=dict(run=str(run),frame=D['k'],time_us=t,facets=D['n'],
              connected_components=res['ncomp'],largest_component_facets=res['comps'][0]['n'],
              branches=stats,core_sensitivity=sensitivity,source_sha256=hashes,
              reference_geometry=True,all_broken_facets_kept=True,
              mode_warning='jointBreakModeRef=slipF in effective deck; exported mode is indicative',
              rendering='Opaque flat-shaded triangles, no edges, no smoothing. Exact planar sections.')
    (out/'mesures_et_provenance.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps({'out':str(out),'frame':D['k'],'n':D['n'],'branches':stats,
                      'sensitivity':sensitivity},ensure_ascii=True,indent=2))


if __name__=='__main__':main()
