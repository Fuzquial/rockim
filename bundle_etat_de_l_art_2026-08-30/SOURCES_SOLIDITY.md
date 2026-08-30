# Provenance des citations `Y3D*.c` — la référence unique
*Arrêtée le 2026-08-30 (action **B4a** du plan de contre-audit). **Toute citation
d'une ligne du code de Solidity dans ce dépôt renvoie ici.***

---

## 1. La source, une fois pour toutes

> **Dépôt** : [`ImperialCollegeLondon/solidity-solver-open`](https://github.com/ImperialCollegeLondon/solidity-solver-open)
> **Licence** : LGPL-3.0 — **Volume** : C, ~17 000 lignes, format `.Y3D`
> **Lignée** : Munjiza (Y2D/Y3D → VGW → VGeST → Solidity), celle de la thèse de
> Guo et des articles de Yang *et al.*
> **Cloné et lu le** : **2026-08-26**
> **Commit** : ⚠️ **NON RELEVÉ** — voir §4, c'est la lacune principale.

C'est **bien le code d'Imperial College London**. Ce point a été contesté en
interne pendant trois jours puis rétabli :
[`chantier_imperial_2026-08-29/A03_resourcer_attributions.md`](chantier/A03_resourcer_attributions.md) §2.
Les noms de valeur `solidity` dans les clés de rockim sont donc **exacts** et ne
seront pas renommés.

## 2. ⚠️ Ce qui n'est PAS ce code

**Ce n'est pas la version qui a produit l'article de 2026.** Son facteur
d'endommagement d'élément y est **câblé à zéro** (`df = R0`) et son DIF est
**neutre** (`dpeftdif = R1`), alors que l'article publie les équations (3)-(4)
d'un modèle d'endommagement. Lecture la plus simple : **version ouverte en retard
sur la version interne** — banal pour un code de recherche.

> **Conséquence de méthode, et elle vaut pour chaque ligne du §3 :** lire une
> **FORME** dans ce code et en conclure une **implémentation de ce que décrit
> l'article de 2026** est une faute. Non pas parce que ce serait le code de
> quelqu'un d'autre — c'est bien le leur, même lignée, mêmes auteurs — mais parce
> que **ce n'est pas la version dont l'article parle**.

**Trois statuts, jamais deux.** Pour chaque ligne : ce que disent les **articles
publiés** ; ce que fait le **code public** ; ce que fait la **version interne**
(inconnue, non consultable). La colonne « statut » du §3 les distingue.

## 3. Les treize références, et rien de plus

Le dépôt porte **79 lignes** citant un fichier `Y3D*.c`, dont **72** avec un
numéro de ligne — mais elles ne renvoient qu'à **13 endroits distincts**. Les
voici tous.

> ⚠️ **La colonne « ce que rockim en dit » rapporte ce que le dépôt AFFIRME
> trouver à cette ligne. Je ne l'ai pas vérifié : le code de Solidity n'est pas
> accessible depuis l'environnement où cette table est écrite.** C'est
> exactement ce que le §4 demande de faire.

| # | référence | ce que **rockim en dit** | clé concernée | statut |
|---|---|---|---|---|
| 1 | `Y3Dfd.c` **1098-1099** | ouverture au pic mode I, `op = R2·el·dpeft/dpepe` | `jointDeltaC = solidity` | **code public seul** |
| 2 | `Y3Dfd.c` **1099** | δ_c = `dnE + max(2·dnE, 3·Gf/ft)` | `jointDeltaC = solidity` | **code public seul** — la convention `guo` (éq. 2.30) oublie l'offset **et** le plancher |
| 3 | `Y3Dfd.c` **1110-1126** | plage d'adoucissement mode II, plancher `2·sp` | `jointDeltaC`, plancher | **partiellement publié** : la part Guo (éq. 2.24/2.30) l'est ; **le plancher `2·sp` n'a aucune source d'article** |
| 4 | `Y3Dfd.c` **1126** | plage de cisaillement divisée par `fs = c + tan(φ)·|σ_n|` | `shearRangeCoulomb` | **code public seul** |
| 5 | `Y3Dfd.c` **1175** | `if((nfail>1)&&...)` — deux points d'intégration au moins | `jointFailRule = majority` | ✅ **DEUX sources indépendantes** — manuscrit UCL **p. 14** : « *at least two integration points have zero stress components* » |
| 6 | `Y3Dfd.c` **1448** | taux de l'élément pris **tel quel**, sans lissage | `strainRateFilter = none` | **code public seul** — aucun article ne décrit la mesure du taux |
| 7 | `Y3Dfd.c` **1448-1456** | `dpeftdif` local à la boucle, repris à chaque pas | `strainRateDIFArm = continuous` | **code public seul** (l'architecture) ⚠️ **et le DIF y est NEUTRE** (`= R1`) : c'est la forme qui est reprise, pas un comportement |
| 8 | `Y3Did.c` **915-964** | naissance de contact sur joint rompu, ré-échelonnement de la pénalité de la paire | `gcBirth = penalty` | ✅ **partiellement racheté** — manuscrit UCL **p. 17** pose le problème et **l'éq. (18)** publie une rampe **linéaire sur ~10 pas**. Seul le **ré-échelonnement** est propre à rockim |
| 9 | `Y3Did.c` **995** | pénalité normale du contact × `d_fact` | `contactDamageCoupling = solidity` | **code public** ⚠️ **et INERTE là-bas** (`df = R0` ⇒ `d_fact = 1`). rockim **active une branche morte** |
| 10 | `Y3Did.c` **1017** | `ktss = 2.0/7.0·penalty` — rapport tangent/normal | `potTangentFactor / potPenaltyFactor` | **code public seul** ⚠️ **aucune valeur de k_t n'est publiée**, sur **8 sources** ([lot 2c §4](fiches/2026-08-29_lot2c_frottement_tangentiel.md), B9) |
| 11 | `Y3Did.c` **1264** | effondrement `/1000` sous `d_fact < 0,041` | `contactDamageCoupling` | **code public**, inerte là-bas (idem n° 9) |
| 12 | `Y3Did.c` **1265** | raideur normale de la paire suivant l'endommagement | `contactDamageCoupling` | **code public**, inerte là-bas (idem n° 9) |
| 13 | `Y3Did.c` **1292** | règle de paire = **minimum** (`mu = min`) | `contactMu.<phase>` | **code public seul** ⚠️ et **la question ne se pose pas chez eux** : dans tous leurs essais publiés les corps sont du **même matériau**. Choix propre à rockim, à assumer — **sans leur reprocher un silence** |

**Bilan des statuts** : **2** rachetées par un article (n° 5 et 8), **1**
partiellement (n° 3), **4** relevées dans un code où le mécanisme est **inerte**
(n° 7, 9, 11, 12), **6** conventions de code sans contrepartie publiée.
**Aucune n'est « sans source »** — c'était l'erreur du lot 4 §1 d'origine.

## 4. ⚠️ Ce qui reste à faire, et cela demande le clone (B4b)

**Les numéros de ligne ci-dessus ne sont pas reproductibles en l'état.** Le dépôt
est **activement maintenu** (dernier push relevé au 2026-03-31) : **les lignes
bougent**. Elles valent pour l'état lu le **2026-08-26**, et **aucun commit n'a
été relevé**. Un rapporteur qui reclone aujourd'hui ne les retrouvera pas
nécessairement.

**Trois gestes, dix minutes avec le clone sous la main** — et impossibles sans :

1. **Relever le commit** : `git -C <clone> rev-parse HEAD` et `git log -1 --date=short`.
   L'inscrire au §1 à la place du « NON RELEVÉ ».
2. **Vérifier les treize** d'un coup :

   ```sh
   for r in 1098 1099 1110 1126 1175 1448; do
       printf '%-14s ' "Y3Dfd.c:$r"; sed -n "${r}p" <clone>/Y3Dfd.c
   done
   for r in 915 995 1017 1264 1265 1292; do
       printf '%-14s ' "Y3Did.c:$r"; sed -n "${r}p" <clone>/Y3Did.c
   done
   ```

   Chaque ligne imprimée doit correspondre à la colonne « ce que rockim en dit ».
   **Toute divergence est un défaut d'attribution**, à consigner ici.
3. **Répondre aux trois questions ouvertes** de la fiche A03 §7 : le dépôt a-t-il
   un tag correspondant aux articles ? le facteur d'endommagement est-il absent de
   **toute** l'histoire du dépôt ou seulement de la tête ? les auteurs disent-ils
   quelque part que le code public diffère de celui des articles ?

## 5. Ce que cette note remplace

Elle ne remplace **rien** : les 72 citations restent où elles sont, **inchangées**.
Ce qui change, c'est qu'elles ont désormais **un seul endroit où renvoyer**. Les
quatre fichiers de code qui en portent (`src/FdemSolver.cpp`,
`src/Fdem3dSolver.cpp` et les deux en-têtes) pointent ici depuis leur tête.

**Pourquoi pas 72 éditions** : elles se réduisent à **13 endroits distincts**.
Écrire treize fois la provenance en table vaut mieux que soixante-douze fois en
commentaire — c'est plus court, c'est vérifiable d'un coup d'œil, et cela évite
soixante-douze occasions d'introduire une erreur dans un diff mécanique.
