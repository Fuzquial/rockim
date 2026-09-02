import re, glob, sys
import numpy as np
def pts(f):
    raw=open(f,'rb').read().decode('utf8','ignore')
    m=re.search(r'<Points>.*?<DataArray[^>]*>(.*?)</DataArray>',raw,re.S)
    P=np.fromstring(m.group(1).strip(),sep=' ').reshape(-1,3)[:,:2]
    m=re.search(r'Name="connectivity"[^>]*>(.*?)</DataArray>',raw,re.S)
    T=np.fromstring(m.group(1).strip(),sep=' ',dtype=int).reshape(-1,3)
    m=re.search(r'Name="sigmaYY"[^>]*>(.*?)</DataArray>',raw,re.S)
    syy=np.fromstring(m.group(1).strip(),sep=' ')
    return P,T,syy
for g,cible in ((0,1.5910),(90,4.3307)):
    fr=sorted(glob.glob(f'out_t14b_{g}/fdem_0*.vtu'))
    P0,T,_=pts(fr[0]); P1,_,syy=pts(fr[-1])
    H=P0[:,1].max()-P0[:,1].min(); y0=P0[:,1].min()
    lo=(P0[:,1]>y0+0.28*H)&(P0[:,1]<y0+0.32*H); hi=(P0[:,1]>y0+0.68*H)&(P0[:,1]<y0+0.72*H)
    L0=P0[hi,1].mean()-P0[lo,1].mean(); L1=P1[hi,1].mean()-P1[lo,1].mean(); eps=(L1-L0)/L0
    c=P0[T].mean(1); band=(c[:,1]>y0+0.3*H)&(c[:,1]<y0+0.7*H)
    s=syy[band].mean(); E=s/eps
    print(f'  litage {g:2d} deg : eps_bande = {eps:.3e}, sigmaYY moyen = {s/1e3:.2f} kPa, E_app = {E/1e9:.4f} GPa  (cible {cible} : ecart {100*(E/1e9/cible-1):+.2f} %)')
