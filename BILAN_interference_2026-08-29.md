# Broyage : qui gouverne le mode de rupture ? — experience du 2026-08-29

## 1. La question

Le run P1 (pulverisation x coulomb x contact residuel) ne casse que **27,8 %**
de ses joints en cisaillement, alors que le correctif `jointShearRange =
coulomb` avait ete adopte sur la dominance INVERSE (Bgrad : 78 % de
cisaillement dans le noyau ; banc polyaxial : 100 %). Le depouillement de P1
et sa contre-expertise ont impute cette chute a `bulkDamage`, sur une chaine
de code exacte : la loi amortit la contrainte (`sig *= bdCd(1-D)`,
Fdem3dSolver.cpp l. 2289), cette contrainte amortie alimente `sigG` (l. 2319),
et le critere d'INSERTION des joints lit `sigG` (l. 1733). Dans le noyau
pulverise la contrainte vue par le critere est donc divisee par ~10 : les
joints n'y naissent pas, donc ne peuvent pas casser.

L'imputation etait plausible. Elle est FAUSSE, et c'est la mesure qui le dit.

## 2. Le protocole

Trois runs, meme maillage (18 186 tets), meme deck, meme binaire, T = 80 us,
une seule variable a la fois. Comparaison A TEMPS EGAL (les fractions derivent
fortement en debut de run : les deux configurations adaptatives passent de
5,3 % et 25,4 % a 30 us pour converger toutes deux vers 28 % a 80 us).

## 3. Resultats (t <= 65 us)

| configuration | joints rompus | cisaillement | % |
|---|---|---|---|
| adaptatif + bulkDamage (= le reglage de P1) | 364 | 105 | **28,8 %** |
| adaptatif SANS bulkDamage | 834 | 217 | **26,0 %** |
| **INTRINSEQUE** + bulkDamage | 1 428 | 592 | **41,5 %** |

Fraction de cisaillement de la configuration intrinseque au fil du temps :
52,0 % (20 us), 48,5 % (30), 42,5 % (40), 40,8 % (50), 41,9 % (60), 41,5 %
(65) — **stabilisee depuis 50 us**, ce n'est pas un transitoire.

## 4. Ce que ca etablit

1. **`bulkDamage` ne gouverne PAS le mode de rupture.** Retire entierement, la
   fraction de cisaillement ne bouge pas : 26,0 % contre 28,8 %. L'imputation
   du depouillement de P1 est refutee sur ce point precis.
2. **`bulkDamage` etouffe la quantite de fissuration** : 364 ruptures contre
   834 sans lui, soit **-56 %**, et 2 137 joints inseres contre 3 290. La
   chaine de code decrite au §1 est reelle ; c'est sa CONSEQUENCE qui avait
   ete mal identifiee. Elle supprime de la casse, elle n'en choisit pas le
   mode.
3. **Le SCHEMA D'INSERTION gouverne le mode.** L'intrinseque donne 41,5 % de
   cisaillement contre 28,8 % pour l'adaptatif a physique identique (facteur
   1,44) et quatre fois plus de ruptures. Explication : le critere adaptatif
   cree un joint quand la contrainte franchit un seuil bati sur la traction ;
   sous l'indenteur, ou tout est comprime, il ne se declenche pas, donc les
   ruptures en cisaillement n'ont jamais de joint sur lequel se produire.
   C'est precisement pourquoi Yang et Solidity inserent tous les joints
   d'emblee : on ne peut pas deviner a l'avance ou la roche cedera en
   cisaillement sous compression.

## 5. Reserves

* Les 41,5 % sont une fraction GLOBALE sur tout le bloc ; les 78 % de Bgrad
  etaient mesures DANS LE NOYAU. Les deux ne se comparent pas directement —
  la valeur au noyau de la configuration intrinseque reste a extraire.
* Le run intrinseque coute cher : tous les joints actifs des le depart au lieu
  de 6-14 % ; ~2 h sur un coeur contre 12 min en adaptatif a quatre coeurs
  pour la meme fenetre de 80 us.
* `strainRateDIFArm = insertion` (le defaut) EXIGE `insertion = adaptive`
  (l. 114-115) : le passage a l'intrinseque impose `strainRateDIFArm =
  envelope`, ce qui change le moment ou le facteur dynamique est fige. Cet
  ecart n'est pas isole dans l'experience et devra l'etre.

## 6. Consequence pratique

Le levier a essayer n'est pas de retirer `bulkDamage` — il ne sert a rien pour
le mode — mais de changer le schema d'insertion, ou de rendre au critere
d'insertion la contrainte NON amortie. Les deux se testent a bas cout sur ce
meme jumeau de 80 us avant d'engager le moindre gros run.

## 7. RESULTATS FINAUX (runs termines, T = 80 us)

| | adaptatif + bulkDamage | adaptatif SANS bulkDamage | INTRINSEQUE + bulkDamage |
|---|---|---|---|
| joints rompus | 605 | 1 498 | **2 240** |
| cisaillement | 28,6 % | 27,5 % | **44,3 %** |
| poste joints (fissuration) | 0,654 J | — | **2,373 J** |
| poste elements | 11,59 J | — | 5,38 J |
| pulverisation | 3,26 J, **25** elements a Dmax | — | 0,114 J, **0** element |
| frottement de contact | 0,510 J | — | 8,645 J |
| contact residuel (WP6) | engage a 42,1 us | — | **JAMAIS engage** |
| energie absorbee par le bloc | 12,87 J | — | 5,36 J |
| vz du bit a 80 us | -8,29 m/s | — | -8,90 m/s |

Le confondant redoute est LEVE par la banniere du solveur : `strainRateDIFArm
= envelope` gele le DIF au franchissement de l enveloppe, c est-a-dire sur LE
MEME critere que l insertion adaptative (contrainte moyenne des deux tetras
contre ft dynamique et Mohr-Coulomb). Les deux schemas partagent le critere et
ne different que par sa consequence : l un cree le joint, l autre s y contente
de geler le DIF. Le passage a l intrinseque n a donc pas change la resistance.

## 8. LES DEUX DECOUVERTES DU RUN INTRINSEQUE

**(a) bulkDamage et l insertion intrinseque s excluent mutuellement.** Sous
l intrinseque la pulverisation devient INERTE : 0,114 J contre 3,26 J, et
ZERO element a Dmax contre 25. Consequence en cascade : `contactResidualMu`
n a JAMAIS ete engage — la banniere de fin l ecrit noir sur blanc, « aucun
element pulverise n a touche un contact ». Le mecanisme n est pas celui qu on
croyait : le critere de pulverisation s ecrit dm = hEl x eps_dev (l. 2274), et
quand les joints peuvent casser librement des t = 0, la deformation
deviatorique est RELACHEE PAR LA FISSURATION et n atteint jamais le seuil.
Les deux modeles ne se genent pas seulement, ils se disputent la meme energie
de deformation : celui qui agit le premier interdit l autre. Le partage
energetique le montre — le poste joints passe de 0,65 a 2,37 J (x3,6) pendant
que le poste elements tombe de 11,59 a 5,38 J. Regime fragile contre regime
de deformation volumique.
NB : Yang et al. font tourner les deux ensemble. L ecart doit donc venir de la
calibration (delta0 = 1,4e-5 m est une valeur GRANITE appliquee au calcaire,
et dm depend de la taille d element) ou de leur loi de joint. A instruire.

**(b) L intrinseque aggrave l injection par la branche NORMALE du contact.**
`contact = -2,471 J (dont frottement 8,645 J)` : la part normale a donc INJECTE
11,1 J, soit 20 % de KE0 — contre 3,66 J (6,9 %) sur P1. L explication est
mecanique : quatre fois plus de joints rompus font quatre fois plus de
surfaces libres, donc de paires de contact qui NAISSENT, et le facteur de
naissance par defaut (`gcBirth = ramp`, rampe par volume) rend la force
normale dependante du chemin. C est un argument fort, et desormais chiffre,
pour poser `gcBirth = penalty` AVANT d adopter l intrinseque.

## 9. CE QU IL FAUT EN RETENIR

L insertion intrinseque restaure le cisaillement (44,3 % contre 28,6 %) mais,
en l etat, elle eteint la pulverisation et fait exploser l injection de
contact. Elle n est donc PAS adoptable telle quelle : il faut d abord
`gcBirth = penalty`, puis re-instruire la calibration de bulkDamage pour le
calcaire. Le run des nombres et le run du facies lances le 29/08 restent
volontairement en insertion ADAPTATIVE, pour ne changer qu une famille de
variables a la fois.

Decks : `J1_intrinseque.cfg`, `J2_coulombseul.cfg` (scratchpad du 29/08).

---

# CORRECTION DU 29/08 (soir) — le code d'Imperial est sur la machine

Le doctorant a demande : « comment font-ils, eux, avec l intrinseque a
Imperial ? ». La reponse a ete cherchee DANS LEUR SOURCE, presente en
`/home/user/solidity`. Elle invalide une partie de ce qui precede. Les
sections 1 a 9 sont conservees telles quelles (on ne reecrit pas
l historique) ; ce qui suit les corrige.

## C1. Ce qu on lit chez eux, verifie de premiere main

* **Solidity est INTEGRALEMENT intrinseque** : zero occurrence de `adaptive`
  ou `extrinsic` dans les 15 448 lignes de `src/`. Leurs joints sont donc
  aussi libres de casser des t = 0 que les notres, et la loi de Yang y vit.
  « Des joints cohesifs partout » NE PEUT PAS etre la cause de notre
  probleme : c est la seule variable sur laquelle les deux codes sont
  identiques.
* **Perimetre des joints** : `I1PEJP  3 0 0 0` (`examples/BST.Y3D`
  l. 117001-117002) — joints dans le PREMIER corps seulement, la roche.
  L acier et le carbure sont des continua sans un seul joint. rockim en pose
  34 507, soit toute face interne des TROIS corps : on amollit l outil et on
  paie le pas de temps du carbure pour rien.
* **Raideur de joint** : `Spring_Stiffness = 900e9 Pa` (`mat.txt` l. 10),
  soit 15 E sur leur granite a 60 GPa — et NON les 52,6 E que
  `DOCUMENTATION_rockim.md` §5.4 avait derives. Le « 3 000 GPa » de cette
  derivation est `D1PEPE`, la penalite d ELEMENT/CONTACT, identique pour le
  granite, l acier et le carbure : elle ne peut pas valoir alpha x E.
* **L INGREDIENT QU ILS ONT ET QUE ROCKIM N A PAS** (`src/Y3Did.c`
  l. 1263-1265, lu verbatim) :
      d_fact = MINIM((R1-d1df[ielem]),(R1-d1df[jelem]));
      if (d_fact < 0.041) d_fact = d_fact/1000.0;
      penalty = penalty*d_fact;
  et `mu = mud*d_fact` (l. 995 et 1044). Chez eux la PENALITE DE CONTACT et
  le FROTTEMENT sont multiplies par (1-D) EN CONTINU des que D > 0, avec un
  effondrement par 1000 au-dela de D = 0,959. C est l effondrement de
  PORTANCE. Chez rockim, `ctcMu` est un echelon binaire a D >= Dmax
  (Fdem3dSolver.hpp l. 497-508) et **la penalite de contact n est jamais
  touchee par D**. WP6 n implemente donc que la moitie tangentielle du
  couplage, et en tout-ou-rien. C est la piste identifiee apres le banc A/B
  du 28/08 (« le canal manquant de l effondrement de portance ») : elle est
  desormais confirmee sur leur source.
* Ils **calibrent avec la complaisance dedans** : E = 60 GPa au deck, aucun
  relevement pour compenser les joints intrinseques.

## C2. Trois resultats du present bilan sont RETIRES

1. **La clause causale du §8(a) est FAUSSE.** « Des joints libres de casser
   relachent la deformation deviatorique AVANT qu elle n atteigne le seuil »
   est contredit par la colonne `bdWork` des history.csv, que personne
   n avait ouverte : a 12 us l intrinseque a deja 6,2e-4 J de pulverisation
   quand l adaptatif est a ZERO ; il MENE la course pendant quatre
   microsecondes, puis GELE (1,37e-2 J a 20 us comme a 30 us) pendant que
   l adaptatif file a 0,760 J. Ce n est pas un non-franchissement, c est un
   emballement ARRETE. La loi est auto-entretenue (`sig *= Cd(1-D)`) : le
   premier mecanisme qui prend verrouille l autre.
2. **Le « facteur 30 » sur le poste pulverisation est trompeur** : il
   contient 2,4 de pur ecart d energie absorbee (12,87 J contre 5,36 J a
   80 us). Normalise par le budget element il vaut x13,3.
3. **Les « 44,3 % contre 28,6 % » melangent deux causes.** A PENALITE EGALE,
   l ecart de schema ne vaut que **+1,5 point** (28,6 % -> 30,1 %) ; les
   ~14 points restants sont de la PENALITE. Mesure a temps egal (30 us),
   normalisee par le budget element : adaptatif pf4 = 29,60 % ; intrinseque
   pf20 = 3,27 % ; intrinseque pf4 = 0,082 %. Levier schema x9,1 ; levier
   penalite x39,7, soit **4,4 fois plus fort**. Le confondant n etait pas
   secondaire, il etait DOMINANT.

## C3. Un defaut de constitution mis au jour

`insertionPenaltyFactor` n est lue NULLE PART sous `insertion = intrinsic`
(l. 1399 et l. 536 sont toutes deux gardees par `adaptive_`), et la penalite
reellement appliquee — `jointPenaltyFactor`, defaut 20 — n est JAMAIS
imprimee. Le deck J1 posait `insertionPenaltyFactor = 4` : la cle etait
MORTE et le run a tourne a 20 E/h sans une ligne de journal. C est exactement
la faute que le depot s interdit (« une capacite active et muette est
indiscernable d une capacite inerte »). CORRECTIF : sortir l impression de la
penalite du bloc `if (adaptive_)` et avertir quand un deck ecrit la cle de
l autre schema.

## C4. Ce qui reste ouvert

* La base unitaire de `delta_m` : leur deck donne `D1PEM0 = 5,0e-3` et
  `D1PEMF = 1,0`, impossibles en METRES sur un tetraedre de 1 mm, plausibles
  comme DEFORMATIONS. rockim en fait une longueur (`dm = hEl x
  sqrt(2/3)||dev eps||`, l. 2274). Si c est une deformation, le facteur `hEl`
  est de trop et la loi devient objective au maillage — ce qui changerait
  aussi la lecture de la sur-pulverisation. A trancher sur l article.
* Le carre 2x2 n est ferme a AUCUN instant ou la pulverisation travaille : il
  manque un run ADAPTATIF a penalite 20 au-dela de 12 us.

---

# CORRECTION 2 (29/08, nuit) — CE QUE `/home/user/solidity` EST VRAIMENT

**AVERTISSEMENT POUR TOUTE SESSION FUTURE.** Le dossier `/home/user/solidity`
n'est PAS le code qui a produit les articles de Yang, Xiang et Latham. F.
Uzquiano l'a signale, et le code le confirme.

PREUVE, `src/Y3Dfd.c` l. 749-760, fonction `CauchyTet4` :

    /*calculate damage factor */
    df=R0;
    *deldam=df;
    ...
    T[i][j]=(R1-df)*(dpemu/detf)*B[i][j]+dpeks*D[i][j];

Le facteur d'endommagement est **cable a zero**. Les parametres `d1pem0`,
`d1pemf`, `d1pedm` sont bien lus au deck (`Y3Drd.c` l. 1171) et transportes
jusqu'a `CauchyTet4` (`Y3Dfd.c` l. 848), mais **le calcul de D a partir de
delta_m — l'equation (4) de l'article — est ABSENT**. Par consequent
`d1df[ielem]` vaut identiquement 0, donc `d_fact = 1`, donc le couplage
`penalty *= d_fact` et `mu = mud*d_fact` (`Y3Did.c` l. 995, 1044, 1263-1265)
est **INERTE dans ce code**.

## Ce qui TOMBE de la CORRECTION 1

La formule « l ingredient qu ils ont et que rockim n a pas, je l ai LU dans
leur source » est trop forte. J ai lu une FORME (l endroit ou le couplage
s appliquerait) et j en ai conclu une IMPLEMENTATION. Ce n est pas la meme
chose, et la distinction est exactement celle que ce depot s impose partout
ailleurs entre ce qui est ETABLI et ce qui est INFERE.

## Ce qui TIENT malgre tout

* La **forme** du couplage est reelle et informative : que `T` porte
  `(1-df)` et que `penalty`/`mu` portent `d_fact` dit ou le couplage se
  branche dans une architecture FDEM de cette famille. C est une piste
  d implementation, PAS une preuve de ce que fait Imperial.
* L **article** (IJRMMS 206, 2026), lui, est une source primaire et il decrit
  explicitement « severe local stiffness degradation, LOSS OF LOAD-BEARING
  CAPACITY, and fragment-support reduction beneath the insert » (p. 4). WP7
  reste donc motive — mais par l ARTICLE, pas par ce code.
* Les valeurs de la **Table 1** viennent de l article : elles sont solides.
* Les faits de deck (`I1PEJP`, l absence d insertion adaptative) decrivent CE
  code, et ne peuvent plus etre attribues a Imperial sans autre source.

## Consequence

Determiner ce qui manque a rockim exige desormais un travail de SOURCES :
articles, theses, codes publies, et non la lecture d un depot local. C est
l objet de [MISSION_etat_de_l_art_2026-08-29.md](MISSION_etat_de_l_art_2026-08-29.md).
