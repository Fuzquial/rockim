import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"STIXGeneral","mathtext.fontset":"stix","font.size":9})
a,b,c=0.63,1.8,6.0
def f(D):
    D=np.clip(D,0,1); s=a+b
    return (1-(s-1)/s*np.exp(D*(a+c*b)/(s*(1-s))))*(a*(1-D)+b*(1-D)**c)
I=np.trapezoid(f(np.linspace(0,1,4001)),np.linspace(0,1,4001))
ft=1.0; ot=1.0          # plage d'adoucissement (normalisee)
dnE=0.12                # ouverture au pic exageree pour la lisibilite
out="/home/user/rockim/docs/rapport_guide/figures/insertion/"
fig,ax=plt.subplots(1,3,figsize=(7.2,2.5))
o=np.linspace(-0.15,1.25,1401)
# (a) intrinseque : Munjiza (parabole) et rockim lineaire, abscisse = ouverture geometrique
def munj(o):
    r=o/dnE
    return np.where(o<0,2*ft*r,np.where(o<=dnE,ft*(2*r-r*r),ft*f((o-dnE)/ot)))
def lin(o):
    return np.where(o<=dnE,ft*o/dnE,ft*f((o-dnE)/ot))
ax[0].plot(o,munj(o),"k-",lw=1.2,label="Munjiza / Wang (parabole)")
ax[0].plot(o,lin(o),"k--",lw=1.0,label="rockim intrinsèque linéaire")
ax[0].fill_between(o,0,munj(o),where=(o>=dnE),color="0.85")
ax[0].axvline(dnE,color="0.5",lw=0.6,ls=":")
ax[0].text(dnE+0.02,0.05,r"$o_p$",fontsize=8)
ax[0].set_title("(a) joint intrinsèque",fontsize=9)
# (b) adaptatif : abscisse = ouverture geometrique depuis l'insertion
og=np.linspace(-0.15,1.25,1401)
yan=np.where(og<0,np.nan,ft*f(og/ot))
ax[1].plot(og,yan,"k-",lw=1.2,label="Yan 2023, rockim adaptatif (née en traction)")
# rockim : effective dn = og + dn0, dn0 = dnE ; compression lineaire jusqu'a og = -dnE
oc=np.linspace(-0.15,0,50)
ax[1].plot(oc,ft*(oc+dnE)/dnE,"k-",lw=1.2)
# pointe : nee a ft/k, k = 1.6
k=1.6; dn0=ft/k*dnE/ft
dn=og+dn0
tip=np.where(dn<=dnE,ft*dn/dnE,ft*f((dn-dnE)/ot))
ax[1].plot(og,tip,"-",color="0.45",lw=1.0,label=r"née en pointe, $k_{tip}=1{,}6$")
# camacho : lineaire depuis t_ins = 1.1 ft (depassement compris), dmF = 2 G/t
G=ft*ot*I; tins=1.1*ft; dmF=2*G/tins
cam=np.where(og<0,np.nan,np.clip(tins*(1-og/dmF),0,None))
ax[1].plot(og,cam,":",color="k",lw=1.2,label=r"camacho, $t_{ins}=1{,}1\,f_t$")
ax[1].set_title("(b) joint inséré",fontsize=9)
# (c) decharge depuis o_max = 0.3 ot
om=0.3; Dm=om/ot; sm=ft*f(Dm)
x=np.linspace(0,om,50)
ax[2].plot(og,yan,"-",lw=0.8,color="0.6")
ax[2].plot(x,sm*x/om,"k-",lw=1.2,label="Yan éq. 17 : sécante vers $o=0$")
xe=np.linspace(-dnE,om,50)
ax[2].plot(xe,sm*(xe+dnE)/(om+dnE),"k--",lw=1.2,label=r"rockim : sécante vers $o=-\delta_{n0}$")
ax[2].set_title("(c) décharge",fontsize=9)
for A in ax:
    A.axhline(0,color="0.3",lw=0.5); A.axvline(0,color="0.3",lw=0.5)
    A.set_xlabel(r"ouverture / $o_t$"); A.set_xlim(-0.17,1.25); A.set_ylim(-0.6,1.75)
    A.legend(fontsize=6.3,frameon=False,loc="upper right")
ax[0].set_ylabel(r"$\sigma / f_t$")
fig.tight_layout()
fig.savefig(out+"fig_lois_insertion.pdf"); fig.savefig(out+"fig_lois_insertion.png",dpi=200)
print("I =",I)
