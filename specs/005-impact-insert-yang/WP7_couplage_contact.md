# WP7 — frottement par phase et couplage continu endommagement -> contact
*2026-08-29 soir. Deux capacites opt-in, defaut bit-identique. Source : la
lecture DIRECTE du code Solidity (/home/user/solidity) et de l article
Yang et al., IJRMMS 206 (2026) 106660 (fiche
biblio_insertion/yang2026_pulverisation.md).*

## 1. Pourquoi

Le depouillement de P1 et la contre-expertise de la question « comment font-
ils, eux, avec l intrinseque ? » ont etabli que rockim differait d Imperial
par DEUX canaux de contact, pas par le schema d insertion :

1. **Le frottement est chez eux une propriete MATERIAU** (Table 1 : granite
   0,18, carbure 0,6, acier 0,6), la paire prenant le MINIMUM des deux
   (Y3Did.c l. 1292). Le deck rockim posait contactMu = 0,6 globalement :
   le frottement roche/roche etait 3,3x trop eleve sur presque tous les
   contacts — coherent avec un poste frottement a 30 % de KE0 chez nous
   contre 65 % chez eux.
2. **La raideur du contact et le frottement suivent (1-D) EN CONTINU**
   (penalty *= d_fact l. 1265, mu = mud*d_fact l. 995 et 1044, effondrement
   /1000 sous d_fact < 0,041 l. 1264). C est l effondrement de PORTANCE
   (« loss of load-bearing capacity », p. 4 de l article). WP6 n en
   implementait que la moitie tangentielle, en echelon binaire a D = Dmax —
   d ou le paradoxe du banc A/B du 28/08 (outil MOINS freine avec le mu
   residuel : la resistance tangentielle baissait, la portance normale
   jamais).

## 2. Les cles

### contactMu.<nom de phase>
Frottement glissant par phase. Regle de paire = minimum des deux phases en
presence (l outil analytique, sans phase, prend le mu de l element touche).
Gardes : exige des phases nommees (throw sinon) ; banniere listant chaque
phase avec sa valeur et son origine (deck | global).

### contactDamageCoupling = off | solidity
En mode solidity : toute force NORMALE de contact impliquant un element
endommage est multipliee par d_fact = min(1-D_i, 1-D_j) — branche potentiel
(via le facteur sc), relais penalite, outil analytique — et le frottement
suit via ctcMu (d_fact au carre sur Ft, comme chez eux ou le cap mu*fn porte
d_fact des deux cotes). Effondrement d_fact /= 1000 sous 0,041 — NB : a
Dmax = 0,9 (Table 1) le plancher vaut 0,1 et cette branche ne s engage
jamais ; elle n existe que pour Dmax > 0,959, comme chez eux. Exception
documentee : Signorini, frottement seulement (le normal y est une CONTRAINTE,
pas une penalite). Gardes : exige bulkDamage = yang (throw) ; EXCLUSIF de
contactResidualMu (throw — le couplage continu contient l echelon de WP6).
Compteurs : nCplEval_ (evaluations degradees), nCplColl_ (effondrements),
tCpl0_ (premier engagement) ; resume de fin de run avec repli « JAMAIS
engage ».

Implemente en MIROIR 2D/3D : 3D = ctcMu + 4 sites normaux ; 2D = ctcMu + 5
sites normaux (paire, relais, PDC, flat, disc).

## 3. Premiere application : les decks Kuru Grey

bench_impact/configs/impact_kuru9.cfg et impact_kuru11.cfg — replique du cas
vitrine (granite, Table 1 en jeu INTEGRAL, montage allege bit+insert a la
vitesse d indentation mesuree chez eux). Predictions inscrites au deck :
e = 0,827 a 5,62 m/s et e = 0,389 a 6,83 m/s — LE verdict est la BASCULE
non lineaire du rebond entre les deux, le phenomene que leur modele de
pulverisation existe pour reproduire (leur §5 et fig. 9d). Reserves
assumees : masse mobile 1,20 contre 1,509 kg, maillage s = 1,5, pas de
piston — comparer le RAPPORT e d abord.

## 4. Etat

* Suite fast 42/42 avec le binaire patche (build_cpl, OMP = 1) : en cours au
  moment de la redaction, verdict au journal suite_cpl.log.
* Demonstration 2D : couplage engage (29 625 evaluations degradees, 0
  effondrement — attendu a Dmax = 0,9), bannieres et resume conformes.
* Trois gardes verifiees au throw (exclusivite, bulkDamage absent, phases
  absentes).
* Fumee des deux decks Kuru : bannieres PAR PHASE et couplage presentes,
  dt inchange (2,93693e-9 s), vitesses correctes.
* kuru9 en file derriere la certification de la suite ; kuru11 attend la
  liberation d un coeur (fin de P2nombres).
* Revue adverse du diff (3 lentilles : fidelite a Solidity, bit-identite et
  concurrence, coherence physique) : en cours, verdict a appliquer avant
  le commit definitif.
