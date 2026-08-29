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

Decks : `J1_intrinseque.cfg`, `J2_coulombseul.cfg` (scratchpad du 29/08).
