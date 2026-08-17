# PATCH 2 — excavation par relâchement de la paroi (FDEM 2D)

*À coller après le PATCH 1, dont il dépend. 5 hunks. Inerte tant que
`excavRelease` n'est pas posé : configs existantes bit-identiques.*

**Nouvelles clés**

| clé | défaut | sens |
|---|---|---|
| `excavRelease` | false | active le relâchement |
| `excavStart` | 0 | instant du début du relâchement [s] (avant : équilibre) |
| `excavRamp` | — | durée du relâchement [s], **obligatoire** > 0 |
| `boreCX`, `boreCY`, `boreSelectR` | — | sélection des faces de la cavité, **clés déjà existantes** (`confineFaces = bore`) |

**Le principe.** Le massif est maillé avec sa cavité et pré-contraint. Les
faces de la paroi portent donc, à t = 0, la traction que la roche excavée
exerçait : **t = σ₀·n**. On la rétablit intégralement (facteur `rel = 1`,
équilibre exact), puis on la fait décroître en cosinus jusqu'à zéro —
convergence-confinement. Appliquer σ₀·n face par face, et non une pression
scalaire, rend le relâchement **exact pour λ ≠ 1**.

```
rel(t) = 1                                    t <= excavStart
       = 1/2 (1 + cos(pi (t - start)/ramp))   pendant la rampe
       = 0                                    ensuite
```

Le travail du relâchement est versé dans `confWork_`, le canal « confinement »
du bilan d'énergie B4 — donc le résidu reste juste sans toucher à
l'assemblage du budget.

---

## Hunk 1/5 — `include/rockim/FdemSolver.hpp`, déclarations de méthodes

**Chercher** (vers la ligne 247) :

```cpp
    void setupConfinement();
```

**Coller juste APRÈS** :

```cpp
    void setupExcavation();                // etude tunnel : relachement paroi
    void excavationForces();
    double excavRelief() const;            // facteur de relachement rel(t)
```

## Hunk 2/5 — `include/rockim/FdemSolver.hpp`, membres

**Chercher** (le bloc ajouté par le patch 1) :

```cpp
    bool hasInsitu_ = false;
    Eigen::Matrix2d insituS_ = Eigen::Matrix2d::Zero();
```

**Coller juste APRÈS** :

```cpp
    // ---- excavation : relachement de la traction de paroi ----------------
    bool excRelease_ = false;
    double excStart_ = 0.0, excRamp_ = 0.0;
    std::vector<BEdge> excEdges_;          // faces ORIGINALES de la cavite
```

## Hunk 3/5 — `src/FdemSolver.cpp`, dans `init()`

**Chercher** (vers la ligne 385) :

```cpp
    setupConfinement();
```

**Remplacer par** :

```cpp
    setupConfinement();
    setupExcavation();                     // no-op si excavRelease = false
```

## Hunk 4/5 — `src/FdemSolver.cpp`, dans `step()`

**Chercher** (vers la ligne 2230) :

```cpp
    confiningForces();                     // no-op when confiningPressure = 0
```

**Remplacer par** :

```cpp
    confiningForces();                     // no-op when confiningPressure = 0
    excavationForces();                    // no-op si excavRelease = false
```

## Hunk 5/5 — `src/FdemSolver.cpp`, les trois fonctions

**Chercher la fin de `confiningForces()`** (vers la ligne 3825) :

```cpp
        // V2/B4 : travail de la pression sur le solide (compteur en v-)
        confWork_ += dt_ * (half.dot(v_[be.na]) + half.dot(v_[be.nb]));
    }
}
```

**Coller juste APRÈS ce `}`** :

```cpp
// ---------------------------------------------------------------------------
// EXCAVATION (etude tunnel EDZ, Wang et al. 2024). Le massif est maille AVEC
// sa cavite et pre-contraint (insituSh / insituSv, patch 1). Les faces de la
// paroi portent donc a t = 0 une traction DESEQUILIBREE : celle que la roche
// excavee exercait, t = sigma0 . n. On la retablit exactement (rel = 1 :
// equilibre parfait, rien ne bouge), puis on la fait decroitre jusqu'a zero.
// C'est la methode convergence-confinement ; elle remplace le "core modulus
// reduction" de l'article (meme etat initial, meme etat final) sans avoir a
// faire varier un module en cours de calcul.
//
// La traction est evaluee FACE PAR FACE comme sigma0 . n : exact pour un champ
// anisotrope (lambda != 1), ce qu'une pression scalaire ne saurait pas faire.
// Convention de normale et de signe identiques a confiningForces() : n sortant
// du solide, et sigma0 . n avec sigma0 = -p I redonne exactement -p n.
// ---------------------------------------------------------------------------
void FdemSolver::setupExcavation() {
    excRelease_ = cfg_.getb("excavRelease", false);
    if (!excRelease_) return;
    if (!hasInsitu_)
        throw std::runtime_error("excavRelease exige une contrainte in situ "
                                 "(insituSh / insituSv) : sans elle il n'y a "
                                 "rien a relacher");
    excStart_ = cfg_.getd("excavStart", 0.0);
    excRamp_ = cfg_.getd("excavRamp", 0.0);
    if (excRamp_ <= 0.0)
        throw std::runtime_error("excavRamp doit etre > 0 : un relachement "
                                 "instantane lance une onde de choc dans le "
                                 "massif (meme piege que confiningRamp = 0). "
                                 "Compter plusieurs transits d'onde.");
    double bcx = cfg_.getd("boreCX", 0.5 * W_);
    double bcy = cfg_.getd("boreCY", 0.5 * H_);
    double bsr = cfg_.getd("boreSelectR", 0.0);
    if (bsr <= 0.0)
        throw std::runtime_error("excavRelease exige boreSelectR > 0 : rayon "
                                 "de selection des faces de paroi autour de "
                                 "(boreCX, boreCY)");
    for (const auto& be : exterior_) {
        Eigen::Vector2d M = 0.5 * (X0_[be.na] + X0_[be.nb]);
        double dx = M.x() - bcx, dy = M.y() - bcy;
        if (dx * dx + dy * dy <= bsr * bsr) excEdges_.push_back(be);
    }
    if (excEdges_.empty())
        throw std::runtime_error("excavRelease : aucune face de cavite "
                                 "selectionnee — verifier boreCX/boreCY/"
                                 "boreSelectR contre la geometrie du maillage");
    double len = 0.0;
    for (const auto& be : excEdges_)
        len += (X0_[be.nb] - X0_[be.na]).norm();
    std::cout << "[FDEM] excavation : " << excEdges_.size()
              << " faces de paroi, perimetre " << len << " m, relachement de "
              << excStart_ << " s a " << (excStart_ + excRamp_) << " s\n";
}

double FdemSolver::excavRelief() const {
    if (!excRelease_) return 0.0;
    if (t_ <= excStart_) return 1.0;
    double x = (t_ - excStart_) / excRamp_;
    if (x >= 1.0) return 0.0;
    return 0.5 * (1.0 + std::cos(M_PI * x));
}

void FdemSolver::excavationForces() {
    if (!excRelease_) return;
    double rel = excavRelief();
    if (rel <= 0.0) return;
    for (const auto& be : excEdges_) {
        Eigen::Vector2d P = X0_[be.na] + u_[be.na];
        Eigen::Vector2d Q = X0_[be.nb] + u_[be.nb];
        Eigen::Vector2d d = Q - P;
        double L = d.norm();
        if (L < 1e-14) continue;
        Eigen::Vector2d n(d.y() / L, -d.x() / L);      // sortante du solide
        // traction que la roche excavee exercait sur la paroi : t = sigma0 . n
        Eigen::Vector2d half = 0.5 * rel * L * thk_ * (insituS_ * n);
        f_[be.na] += half;
        f_[be.nb] += half;
        // B4 : verse dans le canal confinement (deja cable au budget)
        confWork_ += dt_ * (half.dot(v_[be.na]) + half.dot(v_[be.nb]));
    }
}
```

---

## Contrôles immédiats

1. **Suite fast** inchangée (aucune config du parc ne pose `excavRelease`).
2. **V1, charge nulle in situ** : `excavStart` > T ⇒ `rel` reste à 1 tout du
   long ⇒ le massif est en équilibre exact. Attendu : **0 joint cassé**,
   KE ~1e-10 J, résidu B4 < 0,01 %.
   *Si des joints cassent dès t = 0, c'est un signe : la traction de paroi et
   la pré-contrainte ne se compensent pas — regarder d'abord le signe de
   `insituS_` (patch 1), puis la convention de normale.*
3. **V2, Kirsch** : contrôle quantitatif de tout l'ensemble (voir README §4).

## Extension déjà prévue (à ne PAS faire tout de suite)

Le même mécanisme donne gratuitement le **soutènement** : au lieu de descendre
`rel` jusqu'à 0, l'arrêter à `rel = rel_min` — c'est une pression de
soutènement résiduelle, et la courbe convergence-confinement du tunnel s'en
déduit. Une clé `excavRelMin` (défaut 0) suffirait. L'article ne modélise
aucun soutènement ; à garder pour plus tard.
