# Benchmark AbuAisha et al. 2017 — la part reproductible SANS hydro

AbuAisha, Eaton, Priest & Wong, *Hydro-mechanically coupled FDEM framework to
investigate near-wellbore hydraulic fracturing in homogeneous and fractured rock
formations*, J. Petrol. Sci. Eng. **154** (2017) 100-113. Code Y-Geo.

L'architecture mecanique est celle de rockim : triangles Delaunay elastiques +
elements cohesifs, mode I de Hillerborg, mode II en slip-weakening, couplage
mixte **elliptique** (leur eq. 3). Ce qui manque a rockim est le fluide :
calculateur de volume de cavite, compressibilite, pompe. Consequence unique mais
decisive, ecrite dans le source de rockim lui-meme (`confineFaces = bore`) :

> *faces born from cracking receive nothing*

**La pression ne suit pas la fissure.** Tout ce qui depend de la propagation
hydraulique est hors de portee ; tout ce qui se joue AVANT l'amorcage, ou sous
pression imposee, est reproductible.

---

## Les cinq essais reproductibles

| | figure | ce qu'on mesure | cible |
|---|---|---|---|
| **B1** | A.20 / A.21 | ouverture d'une discontinuite sous pression **uniforme** | **solution analytique de Parker (1981), eq. A.1** — w(0) = 0,0640 mm |
| **B2** | 11, branche montante + pic | pression de rupture | **12 MPa** anisotrope (leur eq. 10), ~12,5 MPa leur numerique ; **14,2 MPa** isotrope |
| **B3** | 9, panneau « no defects » | champ de contrainte longitudinale avant amorcage | Kirsch avec pression interne (`tunnel_edz/tools/kirsch_check.py`) |
| **B4** | 7, t = 1,2 ms | facies **a l'amorcage** | etoile radiale (isotrope) / bi-aile sur sigma_H (anisotrope) |
| **B5** | 11 avec joints, pic seul | decalage du seuil du au joint preexistant | **+5,2 %** joints longitudinal et oblique, ~0 % transverse |

**B1 est le seul point de tout l'article confronte a une solution fermee**, et
il est sans hydro par construction : pression uniforme imposee, et resistances
mises a des valeurs « irreal high » pour interdire toute propagation. C'est donc
le controle le plus fort disponible, et il ne teste que de l'elasticite sous
charge suiveuse.

**B5 est reproductible parce que le decalage de seuil est un effet PRE-amorcage** :
il vient de la perturbation du champ de contrainte par le joint, pas de
l'interaction fluide-joint.

## Hors de portee sans developper le couplage HM

Le plateau post-pic a 5,5 MPa (= contrainte lointaine effective), le branchement
et l'incurvation des figures 7 et 10, l'interaction fissure-joint des figures 12,
13 et 15, l'essai de terrain Montney de la figure 19. La microsismicite des
figures 16-17 est un pur post-traitement (leurs eq. 11-13 sur masse et vitesse
nodales, que rockim sort) mais le semis d'evenements suit le trajet de fracture,
donc il ne coincidera pas au-dela de l'amorcage.

## Un signal favorable sur la longueur cohesive

l_cz = E·G_Ic/ft² = 35 GPa × 10 / (5 MPa)² = **14 mm**, pour une maille de 3 mm
et un forage de 100 mm :

> 2 dx = 6 mm  <  l_cz = 14 mm  <  100 mm

**La regle maison est respectee des deux cotes.** C'est le contraire du cas
Heilman (coupe PDC), ou l_cz valait 65 mm pour une passe de 1 mm — violation par
le haut qui interdisait toute lecture du facies.

---

## Mode d'emploi

```
# B1 — Parker (maillage de production : hFine = 0,003)
python bench_abuaisha/tools/make_crack_mesh.py 8.0 8.0 1.5 0.003 0.3 \
       meshes/parker_crack.msh 1
rockim_tun.exe bench_abuaisha\configs\parker.cfg out_parker
python bench_abuaisha/tools/parker_check.py out_parker

# B2/B3/B4 — forage (outil de maillage repris tel quel de l'etude tunnel)
python tunnel_edz/tools/make_circle_mesh.py 8.0 8.0 0.05 0.003 0.4 0.3 \
       meshes/hf_bore.msh 1
rockim_tun.exe bench_abuaisha\configs\hf_aniso.cfg out_hf_aniso
rockim_tun.exe bench_abuaisha\configs\hf_iso.cfg   out_hf_iso
```

Versions grossieres pour degrossir (hFine = 0,012) : `parker_crack_c.msh`
(22 444 triangles) et `hf_bore_c.msh` (12 581 triangles, 26 elements sur le
pourtour). Les maillages de production visent ~105 elements sur le pourtour du
forage, comme l'article.

## Points techniques a connaitre

**Le dedoublement des levres.** Une discontinuite d'epaisseur nulle ne peut pas
etre une fente decoupee : l'ouverture attendue vaut 0,065 mm pour une maille de
3 mm. `make_crack_mesh.py` maille donc le bloc plein avec la ligne imposee comme
contrainte interne, puis dedouble ses noeuds avec le greffon `Crack` de gmsh.
Les levres deviennent des faces exterieures confondues a t = 0, sans joint entre
elles, et `confineFaces = bore` les selectionne — et elles seules.

**La contrainte nette en B1.** L'article pose sigma_h = 15, sigma_v = 10 MPa et
p = 12 MPa ; la solution ne depend que de sigma' = p - sigma_n = 2 MPa. On
applique donc 2 MPa dans un milieu non precontraint : meme sigma', meme
ouverture, et pas de levres plaquees l'une contre l'autre a t = 0 qu'il faudrait
faire gerer au contact general.

**Le piege de nomenclature.** Leur sigma_H est la contrainte LONGITUDINALE (x),
leur sigma_h la TRANSVERSALE (y). La cle rockim `insituSh` designe
l'HORIZONTALE (x). Donc `insituSh` <- leur sigma_H et `insituSv` <- leur
sigma_h : les deux « h » ne designent pas la meme chose.

**Une coquille de l'article.** L'annexe A ecrit « E=45 MPa ». C'est 45 GPa :
w(0) = 2 x 2e6 x 0,96 x 0,75 / E, et leur figure A.21 lit 0,065 mm, d'ou
E = 44,3 GPa.

**Une reserve assumee sur B2.** rockim ne sait pas retarder la rampe de
confinement, qui demarre a t = 0 en meme temps que l'excavation. Celle-ci est
bouclee a 0,2 ms, quand la paroi ne porte que 1,2 MPa, soit 10 % du seuil vise.
A verifier au depouillement que le pic n'en depend pas.
