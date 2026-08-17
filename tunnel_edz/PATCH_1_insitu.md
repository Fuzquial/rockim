# PATCH 1 — contrainte in situ (FDEM 2D)

*À coller tel quel. 3 hunks. Aucun effet tant que `insituSh`/`insituSv` valent
zéro : les configs existantes restent bit-identiques (constitution, principe I).*

**Nouvelles clés**

| clé | défaut | sens |
|---|---|---|
| `insituSh` | 0 | contrainte horizontale, **pression positive** [Pa] (comme `confiningPressure`) |
| `insituSv` | 0 | contrainte verticale, pression positive [Pa] |
| `insituSxy` | 0 | cisaillement, convention tension-positive (nul dans l'article) |

**Le principe.** Le tenseur σ₀ est ajouté à la contrainte de CHAQUE élément :
σ_total = σ₀ + D:ε(u). C'est la méthode classique dite « des contraintes
initiales » ; avec les frontières à rouleaux, un champ σ₀ uniforme est
auto-équilibré, donc l'état initial est un équilibre exact et rien ne bouge
tant que la cavité n'est pas relâchée (patch 2).

**Une seule ligne de physique.** Dans `elementForces()`, la force interne
s'écrit `fe = A0·thk·(P·dN)` avec `P = σ_global·R`. Ajouter σ₀ à la contrainte
globale revient donc à ajouter `σ₀·R` à `P`. Trois conséquences gratuites, car
elles lisent toutes `sigG = P·Rᵀ` :

* la sortie VTU (`sigmaXX/YY/XY`) porte la contrainte **totale** ;
* la jauge `achievedConfinement()` mesure directement σ₀ (contrôle V1) ;
* **le critère d'insertion adaptative** (`insertionSweep()`, qui moyenne
  `A.sxx/B.sxx…`) voit l'in situ — donc les joints s'insèrent sous la
  contrainte totale, et `activateJoint()` stampe la continuité `dn0` et
  `slip` sur cette même contrainte totale. Rien d'autre à faire côté joints.

---

## Hunk 1/3 — `include/rockim/FdemSolver.hpp`

**Chercher** (vers la ligne 666) :

```cpp
    double confP_ = 0.0, confRamp_ = 0.0, confL0_ = 0.0;
    std::vector<BEdge> confEdges_;         // ORIGINAL exterior faces only
```

**Coller juste APRÈS** :

```cpp
    // ---- contrainte in situ (etude tunnel EDZ, Wang et al. 2024) ----------
    // Tenseur de contrainte initial GLOBAL, en convention rockim TENSION
    // POSITIVE (les cles, elles, sont des pressions positives : insituSh =
    // 5e6 signifie 5 MPa de COMPRESSION horizontale). Ajoute a la contrainte
    // de chaque element dans elementForces() : sigma_total = sigma0 + D:eps.
    bool hasInsitu_ = false;
    Eigen::Matrix2d insituS_ = Eigen::Matrix2d::Zero();
```

## Hunk 2/3 — `src/FdemSolver.cpp`, dans `init()`

**Chercher** (vers la ligne 92) :

```cpp
    gravity_ = cfg_.getd("gravity", 0.0);
    if (gravity_ < 0.0)
        throw std::runtime_error("gravity is a magnitude in m/s^2 (it acts along -y): use a positive value");
```

**Coller juste APRÈS** :

```cpp
    // ---- contrainte in situ (etude tunnel EDZ) ---------------------------
    // insituSh / insituSv : contraintes principales horizontale et verticale
    // du massif, en PRESSION POSITIVE (5e6 = 5 MPa de compression), meme
    // convention que confiningPressure. Le coefficient de pression laterale
    // de l'article est lambda = insituSh / insituSv.
    {
        double sh = cfg_.getd("insituSh", 0.0);
        double sv = cfg_.getd("insituSv", 0.0);
        double sxy0 = cfg_.getd("insituSxy", 0.0);
        if (sh < 0.0 || sv < 0.0)
            throw std::runtime_error("insituSh / insituSv sont des PRESSIONS "
                                     "positives (compression) : utiliser 5e6 "
                                     "pour 5 MPa de compression");
        if (sh > 0.0 || sv > 0.0 || sxy0 != 0.0) {
            hasInsitu_ = true;
            insituS_ << -sh, sxy0, sxy0, -sv;
            std::cout << "[FDEM] in situ : sigma_h = " << sh / 1e6
                      << " MPa, sigma_v = " << sv / 1e6
                      << " MPa (compression), lambda = "
                      << (sv > 0.0 ? sh / sv : 0.0) << "\n";
            // NB : adaptive_ n'est lu que plus bas dans init(), on relit donc
            // la cle elle-meme.
            if (cfg_.gets("insertion", "intrinsic") != "adaptive")
                std::cout << "[FDEM] WARNING: contrainte in situ avec "
                             "insertion = intrinsic — les joints intrinseques "
                             "ne transmettent sigma0 qu'apres s'etre fermes de "
                             "sigma0/pj, donc l'etat initial n'est PAS un "
                             "equilibre exact (transitoire parasite, controle "
                             "a charge nulle inexact). insertion = adaptive "
                             "rend l'etat initial exact : les aretes non "
                             "encore inserees sont des noeuds partages.\n";
            if (cfg_.getd("crushCap", 8.0 * cfg_.getd("cohesion", 0.0)) < 1e11)
                std::cout << "[FDEM] NOTE: le cap deviatorique (crushCap) ne "
                             "voit PAS la contrainte in situ (e.svm est "
                             "calcule avant l'ajout). Sous champ anisotrope il "
                             "peut mordre pour une raison purement numerique : "
                             "poser crushCap tres grand si le bulk doit rester "
                             "elastique.\n";
        }
    }
```

> Le dernier `if` relit la clé brute (`cfg_.getd`) et non `crushCapP_` : ce
> vecteur par phase n'est pas encore rempli à cet endroit de `init()`. C'est un
> message d'aide, aucune physique n'en dépend — supprimable si besoin.

## Hunk 3/3 — `src/FdemSolver.cpp`, dans `elementForces()`

**Chercher** (vers la ligne 2439) :

```cpp
        Eigen::Matrix2d P = R * sig;                   // rotate back
```

**Remplacer par** :

```cpp
        Eigen::Matrix2d P = R * sig;                   // rotate back
        // ---- contrainte in situ ---------------------------------------
        // sigma_global_total = R sig R^T + sigma0, et la force interne
        // s'ecrit avec P = sigma_global * R : il suffit donc d'ajouter
        // sigma0 * R ici. sigG (calcule juste en dessous par P * R^T) porte
        // alors la contrainte TOTALE, et c'est elle que lisent la sortie
        // VTU, la jauge achievedConfinement() et — surtout — le critere
        // d'insertion adaptative de insertionSweep().
        if (hasInsitu_) P += insituS_ * R;
```

---

## Contrôle immédiat après compilation

1. `python tools/verify_suite.py --exe rockim_tun.exe --tier fast` → **15/15**,
   valeurs inchangées au bit près (aucune config du parc ne pose `insituS*`).
2. `configs/verif_zeroload_insitu.cfg` → 0 joint cassé, KE ~1e-10 J,
   `achievedConfinement` = **−5,000 MPa** (c'est la mesure directe de σ₀ :
   si elle sort à +5 MPa, le signe de `insituS_` est inversé ; si elle sort à
   0, `hasInsitu_` n'a pas été armé).
