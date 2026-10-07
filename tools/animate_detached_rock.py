"""Actual deformed tetrahedra; fragments reconstructed from effective joint death."""
from pathlib import Path
import sys, json, csv, xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PIL import Image
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'bench_impact/tools'))
from imp_lib import read_vtu, frame_times

def tetread(path):
    tree=ET.parse(path); root=tree.getroot()
    pts=np.fromstring(root.find('.//Points/DataArray').text,sep=' ').reshape(-1,3)
    arr={n.attrib['Name']:np.fromstring(n.text,sep=' ') for n in root.findall('.//DataArray') if 'Name' in n.attrib}
    return pts,arr['connectivity'].astype(int).reshape(-1,4),arr

def main():
    run=ROOT/'out_stanne2025_rock137';out=ROOT/'output/fragments_stanne';out.mkdir(parents=True,exist_ok=True)
    tf=frame_times(run);ks=[k for k in sorted(tf) if k<=12] # frozen complete snapshot
    p0,tets,f0=tetread(run/'fdem3d_0000.vtu')
    rock=f0['phase']==0;ri=np.flatnonzero(rock)
    coords,uid=np.unique(np.round(p0*1e9).astype(np.int64),axis=0,return_inverse=True)
    node_lookup={tuple(p):i for i,p in enumerate(coords)}
    uv=uid[tets[rock]]
    face_local=np.array([[0,1,2],[0,1,3],[0,2,3],[1,2,3]])
    faces=np.sort(uv[:,face_local].reshape(-1,3),axis=1)
    unique,inv,count=np.unique(faces,axis=0,return_inverse=True,return_counts=True)
    order=np.argsort(inv);start=np.r_[0,np.cumsum(count)[:-1]]
    shared=np.flatnonzero(count==2)
    aa=order[start[shared]]//4;bb=order[start[shared]+1]//4
    face_lookup={tuple(unique[fid]):j for j,fid in enumerate(shared)}
    jp,jc,_=read_vtu(run/'fdem3d_joints_0000.vtu')
    ju=np.array([node_lookup[tuple(x)] for x in np.round(jp*1e9).astype(np.int64)])
    jmap=np.array([face_lookup.get(tuple(x),-1) for x in np.sort(ju[jc],axis=1)])
    ref=p0[tets[rock]];vol=np.abs(np.linalg.det(np.stack([ref[:,1]-ref[:,0],ref[:,2]-ref[:,0],ref[:,3]-ref[:,0]],axis=2)))/6
    assert np.all(count<=2)
    images=[];metrics=[]
    for k in ks:
        p,t,a=tetread(run/f'fdem3d_{k:04d}.vtu');assert np.array_equal(t,tets)
        _,_,jf=read_vtu(run/f'fdem3d_joints_{k:04d}.vtu')
        dead=jf['dead']>.5 if 'dead' in jf else jf['tBreak']>=0
        cut=jmap[dead & (jmap>=0)]
        keep=np.ones(len(shared),bool);keep[cut]=False
        graph=coo_matrix((np.ones(keep.sum()),(aa[keep],bb[keep])),shape=(len(ri),len(ri)))
        n,lab=connected_components(graph,directed=False)
        masses=np.bincount(lab,weights=vol);mainid=np.argmax(masses)
        detached=lab!=mainid
        P=(p[tets[rock]]-np.array([.125,.125,.15]))*1000
        vel=a['velocity'].reshape(-1,3)[tets[rock]].mean(axis=1)
        # All exposed faces of detached components (no internal triangulation faces).
        boundary=count==1
        isface=np.repeat(detached,4)
        outward=np.zeros(len(faces),bool)
        outward[order[start[boundary]]]=True
        broken_faceids=shared[~keep]
        outward[order[start[broken_faceids]]]=True
        outward[order[start[broken_faceids]+1]]=True
        selected=np.flatnonzero(outward & isface)
        tid=selected//4;fid=selected%4
        triangles=np.array([P[e,face_local[q]] for e,q in zip(tid,fid)]).reshape(-1,3,3)
        upward=vel[tid,2]>0
        colors=np.where(upward[:,None],np.array([.1,.65,.77,1]),np.array([.87,.4,.12,1]))
        # Only initial top boundary of remaining main body, as a translucent context.
        top=(np.abs((ref[:,face_local,2]-.15)).max(axis=2)<1e-8).reshape(-1)
        topidx=np.flatnonzero(top & ~isface)
        maintri=np.array([P[e//4,face_local[e%4]] for e in topidx]).reshape(-1,3,3)
        if len(maintri):maintri=maintri[np.linalg.norm(maintri.mean(axis=1)[:,:2],axis=1)<26]
        ids=np.flatnonzero(masses>0);ids=ids[ids!=mainid]
        fullyabove=0;upabove=0
        for i in ids:
            mask=lab==i
            if P[mask,:,2].min()>.05:
                fullyabove+=1
                if np.average(vel[mask,2],weights=vol[mask])>0:upabove+=1
        metric=dict(frame=k,time_us=tf[k]*1e6,components_detached=n-1,
                    detached_tets=int(detached.sum()),detached_volume_mm3=float(vol[detached].sum()*1e9),
                    fully_above_005mm=fullyabove,fully_above_and_rising=upabove,
                    max_detached_vertex_z_mm=float(P[detached,:,2].max()) if detached.any() else None)
        metrics.append(metric)
        fig=plt.figure(figsize=(11,5.8));gs=fig.add_gridspec(1,2,left=.03,right=.96,top=.82,bottom=.19,wspace=.12)
        fig.suptitle(f"St Anne - morceaux détachés | {tf[k]*1e6:.0f} µs",x=.05,ha='left',fontsize=17,fontweight='bold')
        fig.text(.05,.875,f"{n-1} composantes séparées du massif | {vol[detached].sum()*1e9:.2f} mm³ | déplacements réels ×1",fontsize=10)
        ax=fig.add_subplot(gs[0,0],projection='3d')
        if len(maintri):ax.add_collection3d(Poly3DCollection(maintri,facecolor='#aaaaaa',edgecolor='none',alpha=.13))
        if len(triangles):ax.add_collection3d(Poly3DCollection(triangles,facecolor=colors,edgecolor='#633e23',linewidth=.08))
        ax.set(xlim=(-23,23),ylim=(-23,23),zlim=(-20,12),xlabel='x [mm]',ylabel='y [mm]',zlabel='z [mm]')
        ax.view_init(elev=18,azim=-62);ax.set_proj_type('ortho');ax.set_box_aspect((46,46,32));ax.tick_params(labelsize=7)
        ax.set_title('Vue oblique ; surface restante translucide',fontsize=10)
        ax=fig.add_subplot(gs[0,1])
        if len(triangles):ax.add_collection(PolyCollection(triangles[:,:,[0,2]],facecolors=colors,edgecolors='#633e23',linewidths=.08))
        ax.axhline(0,color='#666666',linestyle='--',linewidth=1,label='Surface initiale')
        ax.set(xlim=(-23,23),ylim=(-20,12),xlabel='x [mm]',ylabel='z [mm]');ax.set_aspect('equal')
        ax.set_title('Vue latérale : tous les morceaux projetés',fontsize=10)
        ax.text(.02,.96,f"Entièrement au-dessus de z = +0,05 mm : {fullyabove}\nDont vitesse moyenne ascendante : {upabove}",transform=ax.transAxes,va='top',fontsize=8)
        fig.text(.05,.105,'Bleu : vitesse verticale ascendante ; orange : descendante ou nulle. Gris : surface initialement libre du massif.',fontsize=9)
        fig.text(.05,.062,'Détaché = sans liaison de joint au massif ; un morceau peut rester coincé ou en contact. Insert masqué pour voir sous lui.',fontsize=8)
        fig.text(.05,.028,'Connectivité reconstruite avec les joints effectivement morts. Pas de brosse de tri ; aucune amplification des déplacements.',fontsize=8)
        fig.canvas.draw();im=Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy());images.append(im)
        if k==ks[-1]:fig.savefig(out/'fragments_derniere_trame.png',dpi=190)
        plt.close(fig)
        print(metric,flush=True)
    images[0].save(out/'fragments_mouvement.gif',save_all=True,append_images=images[1:],duration=350,loop=0)
    (out/'mesures.json').write_text(json.dumps(metrics,indent=2),encoding='utf-8')
    (out/'LISEZMOI.md').write_text('Animation des tetraedres de roche aux positions exportees (deja deformees), echelle 1.\nLa connectivite est reconstruite en supprimant les faces partagees dont le joint est mort (dead, sinon tBreak), car le champ fragment natif utilise D<1.\nLes couleurs indiquent le signe de la vitesse verticale moyenne par tetraedre, pas une identite stable de fragment.\nUn morceau sans liaison cohesive peut rester coince ou en contact : etre detache ne prouve pas un vol libre.\nLe compteur au-dessus utilise tous les sommets du fragment a z>0,05 mm et le sens de sa vitesse moyenne ponderee par volume. Ce sont des indices cinematiques, pas un test de contact nul.\nLa surface grise est la portion restante de la surface initiale du massif, pas une reconstruction complete des faces nouvellement exposees. Insert masque pour la visibilite.\nTrames figees 0-12, jusqu a 179,994 us ; aucun fichier du run modifie.\nReproduction : python tools/animate_detached_rock.py\n',encoding='utf-8')

if __name__=='__main__':main()
