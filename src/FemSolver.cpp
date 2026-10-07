#include "rockim/FemSolver.hpp"
#include "rockim/Guards.hpp"
#include "rockim/VtkWriter.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

namespace rockim {

// ===========================================================================
// Construction / setup
// ===========================================================================

FemSolver::FemSolver(const Config& cfg, std::string outDir)
    : cfg_(cfg), out_(std::move(outDir)) {}

// ---------------------------------------------------------------------------
// Cles de JOINT refusees en fem (2026-09-07), calquees sur Fem3dSolver.
// Le mode fem est un continuum a NOEUDS PARTAGES : la frontiere de grain y
// est un saut de proprietes, pas une interface qui peut s'ouvrir. Or
// PhaseSet::from LIT gbAlphaTen/Coh/Gf/E/Fric, gbHeteroFactor et les
// gb.<a>.<b>.* : sans cette garde elles seraient consommees, jugees
// legitimes par l'audit des cles, et parfaitement INERTES. Un deck GBM
// porte depuis fdem qui affaiblit ses joints tournerait a raideur pleine
// sans le moindre message.
// ---------------------------------------------------------------------------
void FemSolver::phaseKeyGuards() const {
    const char* why =
        " : c'est une propriete de JOINT (frontiere cohesive). Le mode fem est "
        "un continuum a NOEUDS PARTAGES — il n'existe aucun joint, la frontiere "
        "de grain y est un simple saut de proprietes (contraste de raideur et "
        "de resistance). La cle serait lue, validee et sans le moindre effet : "
        "utiliser mode = fdem pour une interface cohesive.";
    for (const char* k : {"gbAlphaTen", "gbAlphaCoh", "gbAlphaGf", "gbAlphaE",
                          "gbAlphaFric", "gbHeteroFactor"})
        if (cfg_.has(k))
            throw std::runtime_error(std::string("fem : '") + k + "'" + why);
    for (const char* pref : {"gb.", "groupBond."}) {
        auto ks = cfg_.keysWithPrefix(pref);
        if (!ks.empty())
            throw std::runtime_error("fem : '" + ks.front() + "'" + why);
    }
}

// Audit COMPLET de la famille `phase.<nom>.<propriete>` (miroir fem3d).
// keysWithPrefix MARQUE les cles rendues comme consommees : cette fonction
// doit donc etre exhaustive, c'est elle qui remplace l'audit generique.
static void auditPhaseKeysFem(const Config& cfg, const PhaseSet& ps) {
    static const char* kProps[] = {"fraction", "E", "nu", "rho", "ft",
                                   "cohesion", "frictionDeg", "Gf",
                                   "gfShearFactor", "grainSize", nullptr};
    const bool declared = cfg.has("phases");
    std::string known;
    for (const auto& nm : ps.name) known += " " + nm;
    for (const std::string& k : cfg.keysWithPrefix("phase.")) {
        if (!declared)
            throw std::runtime_error("fem : cle '" + k + "' sans cle `phases` — "
                "aucune phase n'est declaree, donc AUCUNE fiche de phase n'est "
                "lue et celle-ci resterait sans le moindre effet (le materiau "
                "global s'appliquerait partout). Declarer les phases : "
                "`phases = <noms separes par des espaces>`, avec mesh = voronoi "
                "comme source de phase.");
        const auto dot = k.rfind('.');
        if (dot == std::string::npos || dot <= 6)
            throw std::runtime_error("fem : cle '" + k + "' malformee — la "
                "syntaxe est phase.<nom>.<propriete>");
        const std::string nm = k.substr(6, dot - 6);
        const std::string prop = k.substr(dot + 1);
        bool okName = false;
        for (const auto& p : ps.name) if (p == nm) okName = true;
        if (!okName)
            throw std::runtime_error("fem : cle '" + k + "' — la phase '" + nm
                + "' n'est pas declaree. Phases declarees par la cle `phases` :"
                + (known.empty() ? std::string(" (aucune)") : known) + ".");
        bool okProp = false;
        for (const char** p = kProps; *p; ++p) if (prop == *p) okProp = true;
        if (okProp) continue;
        throw std::runtime_error(
            "fem : '" + k + "' — '" + prop + "' n'est pas une propriete de "
            "phase. Proprietes admises : fraction, E, nu, rho, ft, cohesion, "
            "frictionDeg, Gf, gfShearFactor, grainSize.");
    }
}

void FemSolver::init() {
    mat_ = Material::from(cfg_);
    PhaseSet::validate(mat_, "global");
    // ---- materiau PAR PHASE (2026-09-07) : la garde des cles de JOINT AVANT
    // PhaseSet::from (qui les consommerait), puis le jeu de phases. Sans la
    // cle `phases`, phases_ = {mat_} et TOUT ce qui suit est bit-identique.
    phaseKeyGuards();
    phases_ = PhaseSet::from(cfg_);
    phasesDeclared_ = cfg_.has("phases");
    auditPhaseKeysFem(cfg_, phases_);

    W_ = cfg_.getd("W", W_);
    H_ = cfg_.getd("H", H_);
    nx_ = cfg_.geti("nx", nx_);
    ny_ = cfg_.geti("ny", ny_);
    T_  = cfg_.getd("T", 200e-6);

    std::string s = cfg_.gets("scenario", "percussion");
    if      (s == "percussion") scen_ = Scenario::PERCUSSION;
    else if (s == "shear")      scen_ = Scenario::SHEAR;
    else if (s == "bar_wave")   scen_ = Scenario::BAR_WAVE;
    else if (s == "tension")    scen_ = Scenario::TENSION;
    else throw std::runtime_error("FemSolver: unknown scenario '" + s
                                  + "' (percussion | shear | bar_wave | "
                                    "tension)");

    // ---- source de maillage : grille (defaut) ou tessellation de Voronoi ---
    const std::string ms = cfg_.gets("mesh", "grid");
    if (ms != "grid" && ms != "voronoi")
        throw std::runtime_error("fem : mesh must be grid | voronoi (got '"
                                 + ms + "')");
    voronoi_ = (ms == "voronoi");
    if (!voronoi_) {
        // Les cles de tessellation posees hors du chemin voronoi seraient lues
        // par personne : meme garde qu'en fem3d, pour la meme raison (un deck
        // GBM ou l'on oublie `mesh = voronoi` tourne comme un maillage
        // homogene en ayant l'air de decrire une microstructure).
        for (const char* k : {"grainSize", "grainJitter", "grainSeeding",
                              "lloydIters", "refineLevels", "vertexMergeFrac",
                              "grainMesh", "grainElemSize", "grainMeshRandom",
                              "grainSizeSpread"})
            if (cfg_.has(k))
                throw std::runtime_error(std::string("fem : '") + k
                    + "' n'est lue que par mesh = voronoi (tessellation de "
                      "grains) ; avec mesh = " + ms + " elle serait sans le "
                      "moindre effet — le maillage vient de la grille "
                      "structuree. Poser mesh = voronoi, ou retirer la cle.");
        if (phasesDeclared_)
            throw std::runtime_error("fem : la cle `phases` exige mesh = "
                "voronoi — la grille structuree n'a AUCUNE source de phase, "
                "les fiches seraient lues et le materiau global s'appliquerait "
                "partout.");
    }

    // Amortissement local de Cundall (2026-09-07). Le defaut vaut 0 sur les
    // trois scenarios historiques — le FEM 2D n'a jamais eu d'amortissement,
    // et le chemin reste bit-identique — et 0,7 en tension, ne le meme jour,
    // pour s'aligner sur fdem et fem3d en quasi-statique.
    damping_ = cfg_.getd("dampingLocal",
                         scen_ == Scenario::TENSION ? 0.7 : 0.0);
    if (!(damping_ >= 0.0) || damping_ >= 1.0)
        throw std::runtime_error("dampingLocal doit etre dans [0 ; 1[ "
                                 "(coefficient de Cundall, sans dimension ; "
                                 "0 = aucun amortissement)");
    if (damping_ > 0.0)
        std::cout << "[FEM] amortissement local de Cundall : dampingLocal = "
                  << damping_ << " (force non visqueuse -alpha |f| sign(v), "
                     "nulle a l'equilibre)\n";

    damageOn_  = cfg_.getb("damage", scen_ != Scenario::BAR_WAVE);
    erodeD_    = cfg_.getd("erodeD", erodeD_);
    strainCap_ = cfg_.getd("strainCap", strainCap_);
    nanEvery_  = cfg_.geti("nanCheckEvery", 256);       // C4 (w20), 0 = off

    // Damage initiation thresholds expressed as equivalent strains:
    //   tension: kappa0_t = ft / E        (Rankine, strain at first cracking)
    //   shear  : kappa0_s = k_DP / (2G)   (DP stress measure mapped to strain)
    DmP_.clear(); dpAlphaP_.clear(); dpKP_.clear();
    kappa0TP_.clear(); kappa0SP_.clear();
    for (const Material& m : phases_.mat) {
        double a = 0.0, k = 0.0;
        m.dpParams(a, k);
        DmP_.push_back(m.Dmat());
        dpAlphaP_.push_back(a);
        dpKP_.push_back(k);
        kappa0TP_.push_back(m.ft / m.E);
        kappa0SP_.push_back(k / (2.0 * m.G()));
    }

    if (voronoi_) buildMeshVoronoi(); else buildMesh();
    finishMesh();
    applyBoundaryConditions();
    placeTool();

    // Penalty contact stiffness. kp ~ E * t gives a contact compliance of the
    // same order as one element row, i.e. stiff enough not to pollute the
    // response, soft enough not to dominate the stable time step (it is
    // accounted for in computeStableDt() anyway).
    // En traction il n'y a pas d'outil : kp_ = 0 retire la limite de contact
    // du pas de temps au lieu de la laisser mordre pour rien.
    const double kpF = cfg_.getd("kpFactor", 1.0);   // lue dans tous les cas
    kp_  = (scen_ == Scenario::TENSION) ? 0.0 : kpF * phases_.maxE() * thk_;
    xiC_ = cfg_.getd("contactXi", xiC_);
    muC_ = cfg_.getd("contactMu", muC_);

    computeStableDt();
    toolKE0_ = tool_.ke();

    if (scen_ == Scenario::BAR_WAVE) {
        barV0_ = cfg_.getd("barV0", 1.0);
        double frac = cfg_.getd("barGaugeFrac", 0.8);
        gaugeX_ = frac * W_;
        double dx = W_ / nx_;
        for (int i = 0; i < (int)X0_.size(); ++i) {
            if (X0_[i].x() < 0.5 * dx) barBcNodes_.push_back(i);
            if (std::abs(X0_[i].x() - gaugeX_) < 0.51 * dx) gaugeNodes_.push_back(i);
        }
    }

    std::cout << "[FEM] " << el_.size() << " CST elements, " << X0_.size()
              << " nodes, dt = " << dt_ << " s, steps = "
              << (long)std::ceil(T_ / dt_) << "\n";
}

void FemSolver::buildMesh() {
    // Structured grid, each quad split into two CSTs with alternating
    // diagonals ("union-jack"-ish) to reduce directional mesh bias in the
    // crack patterns.
    int nnx = nx_ + 1, nny = ny_ + 1;
    double dx = W_ / nx_, dy = H_ / ny_;

    X0_.resize((std::size_t)nnx * nny);
    for (int j = 0; j < nny; ++j)
        for (int i = 0; i < nnx; ++i)
            X0_[(std::size_t)j * nnx + i] = {i * dx, j * dy};

    u_.assign(X0_.size(), Eigen::Vector2d::Zero());
    v_.assign(X0_.size(), Eigen::Vector2d::Zero());
    f_.assign(X0_.size(), Eigen::Vector2d::Zero());
    m_.assign(X0_.size(), 0.0);
    fix_.assign(X0_.size(), 0);
    active_.assign(X0_.size(), 1);

    auto nid = [nnx](int i, int j) { return j * nnx + i; };

    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int a = nid(i, j), b = nid(i + 1, j), c = nid(i + 1, j + 1), d = nid(i, j + 1);
            if ((i + j) % 2 == 0) {           // diagonal a-c
                el_.push_back({{a, b, c}});
                el_.push_back({{a, c, d}});
            } else {                          // diagonal b-d
                el_.push_back({{a, b, d}});
                el_.push_back({{b, c, d}});
            }
        }
    }

}

// ---------------------------------------------------------------------------
// GBM 2D en FEM (2026-09-07), opt-in `mesh = voronoi`.
//
// POURQUOI. Le mode fdem sait deja mailler une microstructure de Voronoi avec
// un maillage de Delaunay a l'interieur des grains ; le mode fem, lui, n'avait
// que la grille structuree. La question du rapport — FEM ou FEMDEM ? — ne se
// tranche qu'a MICROSTRUCTURE EGALE : meme tessellation, meme germe, memes
// fiches de phase, meme maillage intra-grain, et pour seule difference le
// traitement de la discontinuite (joint cohesif qui s'ouvre contre endommagement
// continu a noeuds partages). Sans ce chemin, toute comparaison confondait
// l'effet de la methode et l'effet du maillage.
//
// CE QUE LE MODE fem FAIT DE LA TESSELLATION, ET CE QU'IL N'EN FAIT PAS. Les
// noeuds sont PARTAGES : une frontiere de grain n'est pas une surface qui peut
// s'ouvrir, c'est un saut de proprietes (E, nu, rho, ft, c, phi, Gf) entre deux
// elements voisins. La fissuration intergranulaire n'est donc PAS imposee par
// la geometrie : elle doit emerger du contraste. C'est exactement l'hypothese
// que la comparaison met a l'epreuve, et c'est pourquoi les cles de joint sont
// refusees plus haut au lieu d'etre acceptees et ignorees.
//
// La tessellation reste maitresse du maillage intra-grain : `grainMesh =
// delaunay` + `grainElemSize` + `grainMeshRandom` donnent le maillage non
// structure de la litterature GBM (h/d_grain ~ 0,15-0,18), au lieu de
// l'eventail depuis le centroide qui force toute fissure transgranulaire a
// passer par le centre du grain.
// ---------------------------------------------------------------------------
void FemSolver::buildMeshVoronoi() {
    const double d   = cfg_.reqd("grainSize");
    const double jit = cfg_.getd("grainJitter", 0.5);
    const int lloyd  = cfg_.geti("lloydIters", 2);
    const double mf  = cfg_.getd("vertexMergeFrac", 0.12);
    const int refine = cfg_.geti("refineLevels", 0);
    if (!(d > 0.0))
        throw std::runtime_error("grainSize doit etre > 0 (diametre moyen de "
                                 "grain vise, en metres)");
    const std::string seeding = cfg_.gets("grainSeeding", "hex");
    if (seeding != "hex" && seeding != "random")
        throw std::runtime_error("grainSeeding must be hex | random (got '"
                                 + seeding + "')");
    const std::string gm = cfg_.gets("grainMesh", "fan");
    if (gm != "fan" && gm != "delaunay")
        throw std::runtime_error("grainMesh must be fan | delaunay (got '"
                                 + gm + "')");
    const double gh = cfg_.getd("grainElemSize", 0.0);
    const double gSpread = cfg_.getd("grainSizeSpread", 0.0);
    if (gSpread < 0.0 || gSpread > 1.5)
        throw std::runtime_error("grainSizeSpread doit etre dans [0 ; 1,5] "
                                 "(ecart-type de ln(taille))");
    if (gSpread > 0.0 && seeding != "random")
        throw std::runtime_error("grainSizeSpread exige grainSeeding = random");
    const bool gRandom = cfg_.getb("grainMeshRandom", false);
    if (gRandom && gm != "delaunay")
        throw std::runtime_error("grainMeshRandom exige grainMesh = delaunay");
    std::vector<double> pSize;
    bool anyPSize = false;
    for (const std::string& nm : phases_.name) {
        const double v = cfg_.getd("phase." + nm + ".grainSize", -1.0);
        pSize.push_back(v);
        if (v > 0.0) anyPSize = true;
    }
    if (anyPSize && phases_.n() < 2)
        throw std::runtime_error("phase.<nom>.grainSize n a de sens qu avec "
                                 "plusieurs phases");
    if (!anyPSize) pSize.clear();

    std::mt19937 rng(cfg_.geti("seed", 12345));
    Tessellation T = Tessellation::build(W_, H_, d, jit, lloyd, mf, refine,
                                         phases_.fraction, rng,
                                         seeding == "random",
                                         gm == "delaunay", gh,
                                         gSpread, pSize, gRandom);
    nGrains_ = T.nGrains;

    // ---- COMPACTAGE DES SOMMETS. La tessellation cree une entree de vtx par
    // composante d'union-find et abandonne les faces qui s'effondrent sous le
    // triangle : des sommets peuvent n'etre references par AUCUN triangle. Le
    // chemin FEMDEM ne le voit jamais (il dedouble les noeuds), mais en noeuds
    // partages ce sont des noeuds de MASSE NULLE. Remap deterministe (ordre
    // croissant de l'ancien indice), jamais d'epinglage.
    std::vector<std::array<int, 3>> tris;
    std::vector<int> triGrain, triPhase;
    tris.reserve(T.tri.size());
    triGrain.reserve(T.tri.size());
    triPhase.reserve(T.tri.size());
    for (const auto& t : T.tri) {
        tris.push_back(t.v);
        triGrain.push_back(t.grain);
        if (t.grain < 0 || t.grain >= (int)T.phaseOfGrain.size())
            throw std::runtime_error("mesh = voronoi : triangle sans grain "
                                     "valide");
        triPhase.push_back(T.phaseOfGrain.empty() ? 0
                                                  : T.phaseOfGrain[t.grain]);
    }
    std::vector<char> used = guards::referencedMask(T.vtx.size(), tris);
    std::vector<int> remap(T.vtx.size(), -1);
    X0_.clear();
    X0_.reserve(T.vtx.size());
    std::size_t nDrop = 0;
    for (std::size_t i = 0; i < T.vtx.size(); ++i) {
        if (!used[i]) { ++nDrop; continue; }
        remap[i] = (int)X0_.size();
        X0_.push_back(T.vtx[i]);
    }
    if (nDrop > 0) {
        const double frac = (double)nDrop / (double)T.vtx.size();
        std::cout << "[FEM] voronoi : " << nDrop << " sommets sur "
                  << T.vtx.size() << " (" << 100.0 * frac << " %) "
                     "n'appartiennent a aucun triangle (arete contractee) — "
                     "retires par remap deterministe\n";
        if (frac > 5e-3)
            throw std::runtime_error("mesh = voronoi : " + std::to_string(nDrop)
                + " sommets orphelins (" + std::to_string(100.0 * frac)
                + " % > 0,5 %) — la contraction d'aretes a effondre trop de "
                  "faces : baisser vertexMergeFrac ou augmenter grainSize.");
        for (auto& tt : tris) for (int& vv : tt) vv = remap[(std::size_t)vv];
    }

    u_.assign(X0_.size(), Eigen::Vector2d::Zero());
    v_.assign(X0_.size(), Eigen::Vector2d::Zero());
    f_.assign(X0_.size(), Eigen::Vector2d::Zero());
    m_.assign(X0_.size(), 0.0);
    fix_.assign(X0_.size(), 0);
    active_.assign(X0_.size(), 1);

    // Orientation CCW : finishMesh calcule A = det/2 et la garde
    // checkDegenerate refuse une aire negative. Meme correctif que le chemin
    // mesh = file du mode fdem.
    for (auto& t : tris) {
        const Eigen::Vector2d& A = X0_[t[0]];
        const Eigen::Vector2d& B = X0_[t[1]];
        const Eigen::Vector2d& C = X0_[t[2]];
        if ((B.x() - A.x()) * (C.y() - A.y())
            - (C.x() - A.x()) * (B.y() - A.y()) < 0.0) std::swap(t[1], t[2]);
    }

    el_.clear();
    el_.reserve(tris.size());
    for (std::size_t k = 0; k < tris.size(); ++k) {
        Elem e;
        e.n = tris[k];
        e.grain = triGrain[k];
        e.phase = triPhase[k];
        el_.push_back(e);
    }

    // ---- CE QUE VAUT LE MAILLAGE, imprime AVANT tout calcul ---------------
    // Trois chiffres decident si le resultat sera lisible : elements par
    // grain (viser 35-90, la litterature GBM-FDEM travaille a h/d = 0,15-0,18),
    // la taille d'element equivalente rapportee au grain, et les fractions
    // d'aire REALISEES (l'affectation des phases est un glouton par aire : sur
    // peu de grains l'ecart a la cible est important, et sans cette ligne la
    // composition modale annoncee serait celle du deck, pas celle du calcul).
    {
        double aTot = 0.0;
        std::vector<double> aPh(phases_.n(), 0.0);
        for (std::size_t k = 0; k < tris.size(); ++k) {
            const Eigen::Vector2d& p1 = X0_[tris[k][0]];
            const Eigen::Vector2d& p2 = X0_[tris[k][1]];
            const Eigen::Vector2d& p3 = X0_[tris[k][2]];
            const double a = 0.5 * std::abs((p2.x() - p1.x()) * (p3.y() - p1.y())
                                          - (p3.x() - p1.x()) * (p2.y() - p1.y()));
            aTot += a;
            aPh[triPhase[k]] += a;
        }
        const double perGrain = nGrains_ > 0
                                ? (double)tris.size() / nGrains_ : 0.0;
        const double hEq = tris.empty() ? 0.0
                           : std::sqrt(4.0 * (aTot / tris.size()) / std::sqrt(3.0));
        std::cout << "[FEM] mesh = voronoi : " << nGrains_ << " grains, "
                  << el_.size() << " triangles, " << X0_.size() << " noeuds, "
                     "boite " << W_ << " x " << H_ << " m\n"
                  << "[FEM]   " << perGrain << " elements par grain, taille "
                     "d'element equivalente " << 1e3 * hEq << " mm, rapport "
                     "h/d_grain = " << (d > 0.0 ? hEq / d : 0.0)
                  << " (litterature GBM : 0,15-0,18 ; l'eventail depuis le "
                     "centroide vaut ~0,6 et force les fissures "
                     "transgranulaires par le centre du grain)\n";
        if (phases_.n() > 1)
            for (int p = 0; p < phases_.n(); ++p)
                std::cout << "[FEM]   " << phases_.name[p]
                          << " : fraction d'aire REALISEE "
                          << 100.0 * aPh[p] / aTot << " % (cible "
                          << 100.0 * phases_.fraction[p] << " %)\n";
    }
}

// ---------------------------------------------------------------------------
// finishMesh : geometrie, matrices B, masses concentrees et gardes. Commun a
// la grille structuree et a la tessellation (2026-09-07). La seule difference
// est la masse volumique, prise PAR PHASE — sans la cle `phases` il n'y a
// qu'une fiche et rho vaut mat_.rho, comme avant.
// ---------------------------------------------------------------------------
void FemSolver::finishMesh() {
    for (auto& e : el_) {
        const auto& p1 = X0_[e.n[0]];
        const auto& p2 = X0_[e.n[1]];
        const auto& p3 = X0_[e.n[2]];
        double x1 = p1.x(), y1 = p1.y(), x2 = p2.x(), y2 = p2.y(), x3 = p3.x(), y3 = p3.y();
        double det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
        e.A = 0.5 * det;                       // CCW ordering guarantees det > 0
        double b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
        double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;
        double inv = 1.0 / det;                // = 1/(2A)
        e.B.setZero();
        e.B(0, 0) = b1 * inv; e.B(0, 2) = b2 * inv; e.B(0, 4) = b3 * inv;
        e.B(1, 1) = c1 * inv; e.B(1, 3) = c2 * inv; e.B(1, 5) = c3 * inv;
        e.B(2, 0) = c1 * inv; e.B(2, 1) = b1 * inv;
        e.B(2, 2) = c2 * inv; e.B(2, 3) = b2 * inv;
        e.B(2, 4) = c3 * inv; e.B(2, 5) = b3 * inv;

        double l12 = std::hypot(x2 - x1, y2 - y1);
        double l23 = std::hypot(x3 - x2, y3 - y2);
        double l31 = std::hypot(x1 - x3, y1 - y3);
        double lmax = std::max({l12, l23, l31});
        e.hMin = 2.0 * e.A / lmax;             // smallest altitude
        e.lc   = std::sqrt(2.0 * e.A);         // crack-band width ~ grid spacing

        double mNode = phases_.mat[e.phase].rho * e.A * thk_ / 3.0;
        for (int k = 0; k < 3; ++k) m_[e.n[k]] += mNode;
    }
    {   // C3 (w20) : orphelin, triangle degenere, masse nulle = erreurs nommees
        std::vector<double> areas;
        std::vector<std::array<int, 3>> conn;
        areas.reserve(el_.size());
        conn.reserve(el_.size());
        for (const auto& e : el_) { areas.push_back(e.A); conn.push_back(e.n); }
        guards::checkOrphans(X0_, conn, guards::kNoIds);
        guards::checkDegenerate("fem", areas, conn, X0_, guards::kNoIds);
        guards::checkMasses("fem", m_, X0_, guards::kNoMask, guards::kNoIds);
    }
}

// ---------------------------------------------------------------------------
// applyBoundaryConditions : appuis, mors de traction, frontieres absorbantes.
// ---------------------------------------------------------------------------
void FemSolver::applyBoundaryConditions() {
    int nnx = nx_ + 1, nny = ny_ + 1;
    double dx = W_ / nx_, dy = H_ / ny_;
    auto nid = [nnx](int i, int j) { return j * nnx + i; };
    (void)dx;
    // Dimensionnes ICI et pas plus bas : integrate() indexe cAbsX_[i] a chaque
    // pas, y compris sur les chemins qui sortent tot (traction).
    cAbsX_.assign(X0_.size(), 0.0);
    cAbsY_.assign(X0_.size(), 0.0);
    kAbsX_.assign(X0_.size(), 0.0);
    kAbsY_.assign(X0_.size(), 0.0);

    // ---- scenario = tension (2026-09-07) : mors bas / mors haut -----------
    // Tolerance geometrique : la grille pose ses noeuds exactement sur y = 0
    // et y = H, la tessellation les y snappe aussi (les cellules sont
    // decoupees par le rectangle), donc 1e-9 suffit dans les deux cas. On
    // MESURE le nombre de noeuds saisis : un mors vide donnerait une
    // eprouvette libre et un pic nul, en silence.
    grip_.assign(X0_.size(), 0);
    if (scen_ == Scenario::TENSION) {
        pullV_    = cfg_.getd("pullV", 0.05);
        pullRamp_ = cfg_.getd("pullRamp", 0.0);
        gripFree_ = cfg_.getb("gripLateralFree", false);
        gripSection_ = cfg_.getd("gripSection", W_ * thk_);
        if (!(gripSection_ > 0.0))
            throw std::runtime_error("gripSection doit etre > 0 [m^2]");
        long nBot = 0, nTop = 0;
        for (int i = 0; i < (int)X0_.size(); ++i) {
            if (X0_[i].y() < 1e-9)           { grip_[i] = 1; ++nBot; }
            else if (X0_[i].y() > H_ - 1e-9) { grip_[i] = 2; ++nTop; }
        }
        if (nBot == 0 || nTop == 0)
            throw std::runtime_error("scenario = tension : mors vide ("
                + std::to_string(nBot) + " noeuds en bas, "
                + std::to_string(nTop) + " en haut) — l'eprouvette serait "
                  "libre et le pic nul");
        for (int e = 0; e < (int)el_.size(); ++e) {
            double yc = 0.0;
            for (int a = 0; a < 3; ++a) yc += X0_[el_[e].n[a]].y();
            yc /= 3.0;
            if (yc > H_ / 3.0 && yc < 2.0 * H_ / 3.0) midEl_.push_back(e);
        }
        std::cout << "[FEM] scenario = tension : " << nBot << " noeuds au mors "
                     "bas, " << nTop << " au mors haut, vitesse " << pullV_
                  << " m/s (" << (pullV_ < 0.0 ? "COMPRESSION" : "traction")
                  << "), rampe " << pullRamp_ << " s, section " << gripSection_
                  << " m^2, " << midEl_.size() << " elements dans la jauge "
                     "centrale\n";
        return;                       // pas d'appui ni de frontiere absorbante
    }

    // Boundary conditions: bottom fully fixed (bedrock support) except for
    // the bar-wave test, which needs a free-free bar with a velocity BC.
    if (scen_ != Scenario::BAR_WAVE) {
        for (int i = 0; i < (int)X0_.size(); ++i)
            if (X0_[i].y() < 1e-12) fix_[i] = 3;
        if (cfg_.getb("fixSides", false))
            for (int i = 0; i < (int)X0_.size(); ++i)
                if (X0_[i].x() < 1e-12 || X0_[i].x() > W_ - 1e-12) fix_[i] |= 1;
    }

    // ------------------------------------------------------------------
    // Lysmer-Kuhlemeyer absorbing (quiet) boundaries.
    // A plane wave hitting a boundary of normal n is absorbed exactly (and
    // oblique incidence approximately) by viscous tractions matching the 1D
    // impedances:  t_n = rho c_p v_n,  t_t = rho c_s v_t.  Lumped on each
    // boundary node over its tributary edge length, this gives one dashpot
    // per direction. 'absorbing = sides' treats the lateral faces (bottom
    // stays fixed), 'absorbing = all' also replaces the bottom support by
    // dashpots. Applied implicitly in the velocity update (see integrate()),
    // so it is unconditionally stable and never controls the time step.
    // ------------------------------------------------------------------
    std::string ab = cfg_.gets("absorbing", "none");
    if (scen_ != Scenario::BAR_WAVE && ab != "none") {
        if (ab != "sides" && ab != "all")
            throw std::runtime_error("absorbing must be none | sides | all");
        // Les dashpots sont poses par INDICE DE GRILLE (nid(i,j)) et leur
        // longueur tributaire vaut dy : sur une tessellation ces deux
        // quantites n'ont aucun sens, et la cle serait lue en posant des
        // amortisseurs sur des noeuds arbitraires. Refus explicite plutot
        // qu'une frontiere absorbante fausse et silencieuse.
        if (voronoi_)
            throw std::runtime_error("fem : absorbing = " + ab + " n'est pas "
                "portee sur mesh = voronoi — les amortisseurs de Lysmer y sont "
                "poses par indice de grille et par longueur tributaire "
                "uniforme, deux notions que la tessellation n'a pas. Utiliser "
                "mesh = grid, ou absorbing = none (eprouvette de laboratoire : "
                "les bords SONT libres).");
        // Pure viscous (Lysmer) boundaries have ZERO static stiffness: under
        // the slow part of the loading the block behaves as if free-floating,
        // which puts spurious bending tension at the bottom centre. The
        // viscous-SPRING boundary (Deeks & Randolph 1994) restores the static
        // support of the truncated half-space with a spring in parallel:
        //   k_n = G/R, k_t = G/(2R) per unit boundary length,
        // R being the spreading distance from the loaded zone to that face.
        // The spring is orders softer than the medium, so wave absorption is
        // essentially unaffected.
        double zP = mat_.rho * mat_.cP() * thk_;   // impedance per unit length
        double zS = mat_.rho * mat_.cS() * thk_;
        double G  = mat_.E / (2.0 * (1.0 + mat_.nu));
        double sF = cfg_.getd("absorbSpringFactor", 1.0);   // 0 = pure Lysmer
        double Rside = cfg_.getd("absorbSpringR", 0.5 * W_);
        double Rbot  = cfg_.getd("absorbSpringR", H_);
        for (int j = 0; j < nny; ++j) {            // lateral faces, n = +-x
            double L = dy * ((j == 0 || j == nny - 1) ? 0.5 : 1.0);
            for (int i : {nid(0, j), nid(nnx - 1, j)}) {
                cAbsX_[i] += zP * L;
                cAbsY_[i] += zS * L;
                kAbsX_[i] += sF * G / Rside * L * thk_;
                kAbsY_[i] += sF * G / (2.0 * Rside) * L * thk_;
            }
        }
        if (ab == "all") {                         // bottom face, n = -y
            for (int i = 0; i < nnx; ++i) {
                double L = dx * ((i == 0 || i == nnx - 1) ? 0.5 : 1.0);
                cAbsY_[nid(i, 0)] += zP * L;
                cAbsX_[nid(i, 0)] += zS * L;
                kAbsY_[nid(i, 0)] += sF * G / Rbot * L * thk_;
                kAbsX_[nid(i, 0)] += sF * G / (2.0 * Rbot) * L * thk_;
            }
            for (int i = 0; i < (int)X0_.size(); ++i)
                if (X0_[i].y() < 1e-12) fix_[i] = 0;   // dashpots replace the support
        }
    }
}

void FemSolver::placeTool() {
    if (scen_ == Scenario::BAR_WAVE || scen_ == Scenario::TENSION) return;
    tool_.mass = cfg_.getd("toolMass", 5.0);
    double gap = cfg_.getd("toolGap", 1e-4);

    std::string sh = cfg_.gets("toolShape", scen_ == Scenario::SHEAR ? "disc" : "flat");
    tool_.shape = (sh == "disc") ? Tool::Shape::DISC : Tool::Shape::FLAT;
    tool_.width  = cfg_.getd("toolWidth", 0.02);
    tool_.radius = cfg_.getd("toolRadius", 0.01);

    if (scen_ == Scenario::PERCUSSION) {
        tool_.motion = Tool::Motion::FREE;
        double vImp = cfg_.getd("impactSpeed", 15.0);
        double xc = cfg_.getd("toolX", 0.5 * W_);
        if (tool_.shape == Tool::Shape::FLAT) tool_.x = {xc, H_ + gap};
        else                                  tool_.x = {xc, H_ + tool_.radius + gap};
        tool_.v = {0.0, -vImp};
    } else {  // SHEAR: prescribed lateral pass at a fixed depth of cut
        tool_.motion = Tool::Motion::PRESCRIBED;
        tool_.shape  = Tool::Shape::DISC;   // lateral cutting needs a 2D normal
        double depth = cfg_.getd("cutDepth", 0.003);
        double vCut  = cfg_.getd("cutSpeed", 10.0);
        tool_.x = {-tool_.radius - gap, H_ - depth + tool_.radius};
        tool_.v = {vCut, 0.0};
    }
}

// ---------------------------------------------------------------------------
// Stable explicit time step.
//
// Central difference / velocity-Verlet is conditionally stable:
//     dt <= dt_crit = 2 / omega_max.
// For a lumped-mass CST mesh, omega_max is bounded element-wise and the usual
// CFL estimate is
//     dt_crit ~ min_e ( h_min(e) / c_p ),
// where h_min is the smallest element altitude and c_p the plane-strain
// dilatational wave speed (the fastest wave). The penalty contact adds a
// spring k_p on single nodes, i.e. an extra frequency sqrt(k_p/m_node), so
// its own limit 2*sqrt(m_min/k_p) is taken into account too. A safety factor
// (default 0.7) covers the softening branch of the damage law and the
// contact nonlinearity.
// ---------------------------------------------------------------------------
void FemSolver::computeStableDt() {
    // CFL prise sur la phase la PLUS RAIDE : une seule phase rapide fixe le
    // pas de temps du bloc entier. Sans la cle `phases` c'est mat_.cP().
    double cP = phases_.maxCp();
    double hMin = 1e30;
    for (const auto& e : el_) hMin = std::min(hMin, e.hMin);
    double dtMesh = hMin / cP;

    double mMin = 1e30;
    for (double mm : m_) mMin = std::min(mMin, mm);
    double dtContact = (kp_ > 0) ? 2.0 * std::sqrt(mMin / kp_) : 1e30;

    double cfl = cfg_.getd("cfl", 0.7);
    dt_ = cfl * std::min(dtMesh, dtContact);
}

// ===========================================================================
// Time stepping
// ===========================================================================

void FemSolver::step() {
    for (auto& fv : f_) fv.setZero();

    if (scen_ != Scenario::BAR_WAVE && scen_ != Scenario::TENSION) applyContact();
    internalForcesAndDamage();
    if (activeDirty_) refreshActiveNodes();
    integrate();

    t_ += dt_;
    ++nanStep_;                          // C4 (w20) : detecteur reel
    if (nanEvery_ > 0 && nanStep_ % nanEvery_ == 0) checkFinite();
}

void FemSolver::checkFinite() {
    static const char* const names[3] = {"u", "v", "f"};
    guards::checkFinite("FEM", nanStep_, t_, X0_.size(), 3, names,
        [&](std::size_t i, int k) -> const Eigen::Vector2d& {
            return k == 0 ? u_[i] : k == 1 ? v_[i] : f_[i];
        },
        [&](std::size_t i) -> const Eigen::Vector2d& { return X0_[i]; },
        [&](std::size_t i) -> long {
            for (std::size_t e = 0; e < el_.size(); ++e)
                for (int a = 0; a < 3; ++a)
                    if ((std::size_t)el_[e].n[a] == i) return (long)e;
            return -1;
        },
        {{"work", work_}, {"|toolF|", tool_.F.norm()}});
}

// ---------------------------------------------------------------------------
// Rigid tool <-> mesh contact, node-to-rigid-surface penalty:
//     f_n = k_p * penetration + c * v_rel_n,   c = 2 xi sqrt(k_p m_node)
// clamped to compression only (no adhesion), plus regularized Coulomb
// friction f_t = -mu |f_n| tanh(v_slip / v_reg).
// ---------------------------------------------------------------------------
void FemSolver::applyContact() {
    tool_.resetForce();
    for (int i = 0; i < (int)X0_.size(); ++i) {
        if (!active_[i]) continue;             // orphan (fully eroded) nodes fly free
        Eigen::Vector2d p = X0_[i] + u_[i];
        Eigen::Vector2d Fc = Eigen::Vector2d::Zero();

        if (tool_.shape == Tool::Shape::FLAT) {
            if (std::abs(p.x() - tool_.x.x()) > 0.5 * tool_.width) continue;
            double pen = p.y() - tool_.x.y();  // node above the bottom face
            if (pen <= 0) continue;
            double vreln = v_[i].y() - tool_.v.y();
            double c = 2.0 * xiC_ * std::sqrt(kp_ * m_[i]);
            double fn = kp_ * pen + c * vreln;
            if (fn < 0) fn = 0;
            double vslip = v_[i].x() - tool_.v.x();
            double ftang = -muC_ * fn * std::tanh(vslip / vReg_);
            Fc = {ftang, -fn};                 // tool pushes the node down
        } else {                               // DISC
            Eigen::Vector2d d = p - tool_.x;
            double dist = d.norm();
            if (dist >= tool_.radius || dist < 1e-14) continue;
            Eigen::Vector2d n = d / dist;      // outward normal (pushes node out)
            Eigen::Vector2d tdir(-n.y(), n.x());
            double pen = tool_.radius - dist;
            Eigen::Vector2d vrel = v_[i] - tool_.v;
            double c = 2.0 * xiC_ * std::sqrt(kp_ * m_[i]);
            double fn = kp_ * pen - c * vrel.dot(n);
            if (fn < 0) fn = 0;
            double ftang = -muC_ * fn * std::tanh(vrel.dot(tdir) / vReg_);
            Fc = fn * n + ftang * tdir;
        }

        f_[i] += Fc;
        tool_.F -= Fc;                         // Newton's third law
    }
    peakF_ = std::max(peakF_, tool_.F.norm());
    // Energy delivered by the tool (thrust/cutting work). Tool-side
    // bookkeeping deliberately includes the interface friction and dashpot
    // dissipation: that is part of the drilling energy input.
    work_ += -tool_.F.dot(tool_.v) * dt_;
}

// ---------------------------------------------------------------------------
// Internal forces + damage.
//
// Small-strain, total-strain format:  eps = B u_e,  effective stress
// s_eff = D eps (with s_zz = nu (s_xx + s_yy) in plane strain), nominal
// stress s = (1 - D) s_eff, and  f_int = B^T s A t.
//
// Damage model (engineering-grade continuum damage, see README limitations):
//  * Tension (Rankine): driver  e_t = max principal effective stress / E,
//    threshold kappa0_t = ft/E.
//  * Shear (Drucker-Prager): driver e_s = (sqrt(J2) + alpha I1)/(2G),
//    threshold kappa0_s = k/(2G). Pure hydrostatic compression never damages
//    (alpha I1 << 0 keeps the driver below threshold), which is the point of
//    using DP instead of a von Mises-type surface for rock.
//  * Both use the exponential softening law
//        D = 1 - (kappa0/kappa) exp( -(kappa - kappa0) / eps_f )
//    with eps_f tied to the fracture energy and the element size through the
//    crack-band scaling (Bazant-Oh):  eps_f = G_f / (f * l_c). This keeps the
//    dissipated energy per unit crack area ~ G_f independent of the mesh.
//  * D = max(D_t, D_s), monotonic. Element erosion at D >= erodeD (or at a
//    strain cap, to kill inverted junk) removes the element => cracks and
//    material removal appear as bands/craters of eroded elements.
// ---------------------------------------------------------------------------
void FemSolver::internalForcesAndDamage() {
    for (auto& e : el_) {
        if (e.eroded) continue;

        Eigen::Matrix<double, 6, 1> ue;
        for (int k = 0; k < 3; ++k) {
            ue(2 * k)     = u_[e.n[k]].x();
            ue(2 * k + 1) = u_[e.n[k]].y();
        }
        Eigen::Vector3d eps = e.B * ue;

        // strain cap: erode grossly distorted elements before they misbehave
        if (std::abs(eps(0)) > strainCap_ || std::abs(eps(1)) > strainCap_ ||
            std::abs(eps(2)) > 2.0 * strainCap_) {
            erode(e);
            continue;
        }

        Eigen::Vector3d sEff = DmP_[e.phase] * eps;
        double szz = phases_.mat[e.phase].nu * (sEff(0) + sEff(1));  // plane strain

        if (damageOn_) {
            updateDamage(e, sEff, szz);
            if (e.D >= erodeD_) { erode(e); continue; }
        }

        Eigen::Vector3d sig = (1.0 - e.D) * sEff;

        // stored for output
        double sm = (sig(0) + sig(1) + (1.0 - e.D) * szz) / 3.0;
        double dxx = sig(0) - sm, dyy = sig(1) - sm, dzz = (1.0 - e.D) * szz - sm;
        e.smean = sm;
        e.syy = sig(1);                      // jauge axiale du scenario tension
        e.svm = std::sqrt(1.5 * (dxx * dxx + dyy * dyy + dzz * dzz) + 3.0 * sig(2) * sig(2));

        Eigen::Matrix<double, 6, 1> fe = e.B.transpose() * sig * (e.A * thk_);
        for (int k = 0; k < 3; ++k) {
            f_[e.n[k]].x() -= fe(2 * k);
            f_[e.n[k]].y() -= fe(2 * k + 1);
        }
    }
}

void FemSolver::updateDamage(Elem& e, const Eigen::Vector3d& sEff, double szz) {
    double sxx = sEff(0), syy = sEff(1), txy = sEff(2);

    // In-plane principal + out-of-plane
    double ctr = 0.5 * (sxx + syy);
    double R   = std::sqrt(0.25 * (sxx - syy) * (sxx - syy) + txy * txy);
    double s1  = std::max(ctr + R, szz);

    // Invariants (tension positive), z included
    double I1 = sxx + syy + szz;
    double sm = I1 / 3.0;
    double dxx = sxx - sm, dyy = syy - sm, dzz = szz - sm;
    double J2 = 0.5 * (dxx * dxx + dyy * dyy + dzz * dzz) + txy * txy;

    // --- tension (Rankine) ---
    // Fiche et constantes derivees de la PHASE de l'element (2026-09-07).
    // Sans la cle `phases` il n'y a qu'une fiche : mp == mat_, k0T == kappa0T_
    // de l'ancien code, etc. — chemin bit-identique.
    const Material& mp = phases_.mat[e.phase];
    const double k0T = kappa0TP_[e.phase], k0S = kappa0SP_[e.phase];

    double et = s1 / mp.E;
    if (et > e.kappaT) e.kappaT = et;
    double Dt = 0.0;
    if (e.kappaT > k0T) {
        double efT = std::max(mp.Gf / (mp.ft * e.lc), 1.5 * k0T);
        Dt = 1.0 - (k0T / e.kappaT) * std::exp(-(e.kappaT - k0T) / efT);
    }

    // --- shear (Drucker-Prager) ---
    double q = std::sqrt(J2) + dpAlphaP_[e.phase] * I1;
    double es = q / (2.0 * mp.G());
    if (es > e.kappaS) e.kappaS = es;
    double Ds = 0.0;
    if (e.kappaS > k0S) {
        double efS = std::max(mp.gfShearFactor * mp.Gf
                              / (dpKP_[e.phase] * e.lc), 1.5 * k0S);
        Ds = 1.0 - (k0S / e.kappaS) * std::exp(-(e.kappaS - k0S) / efS);
    }

    double Dn = std::min(std::max(Dt, Ds), 0.999);
    if (Dn > e.D) e.D = Dn;
}

void FemSolver::erode(Elem& e) {
    e.eroded = true;
    e.D = 1.0;
    erodedVol_ += e.A * thk_;
    ++nEroded_;
    activeDirty_ = true;
}

void FemSolver::refreshActiveNodes() {
    std::fill(active_.begin(), active_.end(), 0);
    for (const auto& e : el_)
        if (!e.eroded)
            for (int k = 0; k < 3; ++k) active_[e.n[k]] = 1;
    activeDirty_ = false;
}

// Central-difference / leapfrog kick-drift:
//   v^{n+1/2} = v^{n-1/2} + dt a^n ;  u^{n+1} = u^n + dt v^{n+1/2}
void FemSolver::integrate() {
    // ---- traction / compression : la REACTION DES MORS se lit sur f_ AVANT
    // que la mise a jour des vitesses ne l'ecrase. La jauge centrale mesure
    // sigma_yy en moyenne d'aire sur le tiers median : c'est elle qui vaut
    // contrainte de l'essai, la reaction des mors incluant l'inertie du mors.
    if (scen_ == Scenario::TENSION) {
        gripF_.setZero();
        for (int i = 0; i < (int)X0_.size(); ++i)
            if (grip_[i] == 2) gripF_ += f_[i];
        sigmaPeak_ = std::max(sigmaPeak_, std::abs(gripF_.y()) / gripSection_);
        double s = 0.0, aa = 0.0;
        for (int e : midEl_)
            if (!el_[e].eroded) { s += el_[e].syy * el_[e].A; aa += el_[e].A; }
        sigMid_ = aa > 0.0 ? s / aa : 0.0;
    }

    for (int i = 0; i < (int)X0_.size(); ++i) {
        if (scen_ == Scenario::TENSION && grip_[i]) {
            // mors : y impose, x libre ou bloque selon gripLateralFree.
            // gripLateralFree = true retire le frettage des mors, qui en
            // compression cree un cone de confinement et gonfle l'UCS.
            if (gripFree_) v_[i].x() += dt_ * f_[i].x() / m_[i];
            else           v_[i].x() = 0.0;
            if (grip_[i] == 1) {
                v_[i].y() = 0.0;                  // mors bas : encastrement
            } else {
                double vg = pullV_;
                if (pullRamp_ > 0.0 && t_ < pullRamp_)
                    vg *= 0.5 * (1.0 - std::cos(M_PI * t_ / pullRamp_));
                v_[i].y() = vg;
            }
            continue;
        }
        // boundary spring (viscous-spring quiet boundary), anchored at X0
        if (active_[i] && (kAbsX_[i] > 0 || kAbsY_[i] > 0)) {
            f_[i].x() -= kAbsX_[i] * u_[i].x();
            f_[i].y() -= kAbsY_[i] * u_[i].y();
        }
        // Amortissement local de Cundall, MEME expression qu'en fdem et fem3d :
        // par composante, on retire alpha |f| dans le sens du mouvement. La
        // vitesse n'intervient que par son SIGNE, donc la force s'annule
        // exactement a l'equilibre et n'introduit aucune echelle de temps.
        // Elle est appliquee APRES les ressorts de frontiere et AVANT
        // l'acceleration, comme dans les deux autres solveurs.
        if (damping_ > 0.0) {
            for (int a2 = 0; a2 < 2; ++a2) {
                const double s = v_[i](a2) > 0 ? 1.0
                               : (v_[i](a2) < 0 ? -1.0 : 0.0);
                f_[i](a2) -= damping_ * std::abs(f_[i](a2)) * s;
            }
        }
        Eigen::Vector2d a = f_[i] / m_[i];
        if (!(fix_[i] & 1)) v_[i].x() += dt_ * a.x(); else v_[i].x() = 0;
        if (!(fix_[i] & 2)) v_[i].y() += dt_ * a.y(); else v_[i].y() = 0;
        // Lysmer dashpot, implicit in v:  m dv/dt = F - c v  discretized as
        //   v <- (v + dt F/m) / (1 + dt c/m)
        // exact impedance at low frequency, unconditionally stable.
        if (cAbsX_[i] > 0) v_[i].x() /= 1.0 + dt_ * cAbsX_[i] / m_[i];
        if (cAbsY_[i] > 0) v_[i].y() /= 1.0 + dt_ * cAbsY_[i] / m_[i];
    }
    if (scen_ == Scenario::BAR_WAVE)
        for (int i : barBcNodes_) v_[i] = {barV0_, 0.0};   // velocity step at x=0

    for (int i = 0; i < (int)X0_.size(); ++i) u_[i] += dt_ * v_[i];

    if (scen_ != Scenario::BAR_WAVE && scen_ != Scenario::TENSION) {
        tool_.integrate(dt_);
        if (tool_.motion == Tool::Motion::FREE && tool_.v.y() > 0 &&
            cfg_.getb("stopOnRebound", false)) { /* optional early stop hook */ }
    }

    if (scen_ == Scenario::BAR_WAVE && tArrive_ < 0)
        for (int i : gaugeNodes_)
            if (std::abs(v_[i].x()) > 0.1 * std::abs(barV0_)) { tArrive_ = t_; break; }
}

// ===========================================================================
// Output
// ===========================================================================

void FemSolver::writeFrame(int frame) {
    std::vector<Eigen::Vector2d> pts(X0_.size());
    for (std::size_t i = 0; i < X0_.size(); ++i) pts[i] = X0_[i] + u_[i];

    std::vector<std::array<int, 3>> tris;
    std::vector<double> dmg, svm, smean, phs, grn;
    for (const auto& e : el_) {
        if (e.eroded) continue;                // eroded elements = open cracks/crater
        tris.push_back(e.n);
        dmg.push_back(e.D);
        svm.push_back(e.svm);
        smean.push_back(e.smean);
        phs.push_back((double)e.phase);
        grn.push_back((double)e.grain);
    }

    char name[64];
    std::snprintf(name, sizeof(name), "/fem_%04d.vtu", frame);
    vtk::ScalarField sf{{"damage", &dmg}, {"vonMises", &svm},
                        {"meanStress", &smean}};
    // phase / grain n'existent que sur le chemin GBM : sur la grille ils
    // vaudraient 0 partout et alourdiraient chaque frame pour rien.
    if (voronoi_) { sf["phase"] = &phs; sf["grain"] = &grn; }
    vtk::writeTriMesh(out_ + name, pts, tris, sf,
                      {{"velocity", &v_}, {"displacement", &u_}});

    // frame index -> (time, tool pose): lets post-processing (make_gif.py)
    // place the rigid tool without re-deriving the frame stride.
    std::ofstream fm(out_ + "/frames.csv",
                     frame == 0 ? std::ios::trunc : std::ios::app);
    if (frame == 0) fm << "frame,t,toolX,toolY\n";
    fm << frame << "," << t_ << "," << tool_.x.x() << "," << tool_.x.y() << "\n";
}

void FemSolver::historyHeader(std::ostream& os) const {
    if (scen_ == Scenario::BAR_WAVE) { os << "t,gaugeVx\n"; return; }
    if (scen_ == Scenario::TENSION) {
        os << "t,sigmaGrip,sigmaMid,strain,nEroded\n";
        return;
    }
    os << "t,toolFx,toolFy,toolX,toolY,toolVx,toolVy,work,toolKE,"
          "erodedVol,nEroded,specificEnergy\n";
}

void FemSolver::historyRow(std::ostream& os) const {
    if (scen_ == Scenario::TENSION) {
        // Deformation nominale des mors : |pullV| t / H (la rampe la surestime
        // un peu au demarrage — la jauge centrale sigmaMid reste la mesure de
        // reference, la reaction des mors portant aussi leur inertie).
        const double eps = (H_ > 0.0) ? std::abs(pullV_) * t_ / H_ : 0.0;
        os << t_ << "," << std::abs(gripF_.y()) / gripSection_ << ","
           << sigMid_ << "," << eps << "," << nEroded_ << "\n";
        return;
    }
    if (scen_ == Scenario::BAR_WAVE) {
        double vg = 0;
        for (int i : gaugeNodes_) vg = std::max(vg, std::abs(v_[i].x()));
        os << t_ << "," << vg << "\n";
        return;
    }
    double Es = erodedVol_ > 0 ? work_ / erodedVol_ : 0.0;
    os << t_ << "," << tool_.F.x() << "," << tool_.F.y() << ","
       << tool_.x.x() << "," << tool_.x.y() << ","
       << tool_.v.x() << "," << tool_.v.y() << ","
       << work_ << "," << tool_.ke() << ","
       << erodedVol_ << "," << nEroded_ << "," << Es << "\n";
}

void FemSolver::finalize() {
    if (nanEvery_ > 0) checkFinite();    // C4 (w20) : dernier controle
    // Convenience CSV snapshot of the final damage field (easy to plot)
    std::ofstream fe(out_ + "/fem_final_elements.csv");
    fe << "cx,cy,damage,eroded\n";
    for (const auto& e : el_) {
        Eigen::Vector2d c = Eigen::Vector2d::Zero();
        for (int k = 0; k < 3; ++k) c += (X0_[e.n[k]] + u_[e.n[k]]) / 3.0;
        fe << c.x() << "," << c.y() << "," << e.D << "," << (e.eroded ? 1 : 0) << "\n";
    }

    std::cout << "\n[FEM] ---- summary ----\n";
    if (scen_ == Scenario::TENSION) {
        const bool comp = pullV_ < 0.0;
        std::cout << "[FEM] " << (comp ? "compression" : "traction")
                  << " uniaxiale ("
                  << (voronoi_ ? "grains de Voronoi" : "grille uniforme")
                  << (phases_.n() > 1
                      ? ", " + std::to_string(phases_.n()) + " phases"
                      : std::string(", monophase"))
                  << ") :\n"
                  << "[FEM]   peak macro stress = " << sigmaPeak_ / 1e6
                  << " MPa (reaction des mors / " << gripSection_ << " m^2)\n"
                  << "[FEM]   ft de la fiche globale = " << mat_.ft / 1e6
                  << " MPa, rapport = "
                  << (mat_.ft > 0.0 ? sigmaPeak_ / mat_.ft : 0.0) << "\n"
                  << "[FEM]   elements erodes = " << nEroded_ << " / "
                  << el_.size() << " ("
                  << 100.0 * (double)nEroded_
                     / (double)std::max<std::size_t>(el_.size(), 1) << " %)\n";
        if (nEroded_ == 0)
            std::cout << "[FEM]   AUCUN element erode : le pic ci-dessus n'est "
                         "PAS une rupture — allonger T ou augmenter |pullV|\n";
        return;
    }
    if (scen_ == Scenario::BAR_WAVE) {
        double cTh = mat_.cBar();
        if (tArrive_ > 0) {
            double cMe = gaugeX_ / tArrive_;
            double err = 100.0 * (cMe - cTh) / cTh;
            std::cout << "[FEM] bar-wave verification: measured c = " << cMe
                      << " m/s, theory sqrt(E/rho) = " << cTh
                      << " m/s, error = " << err << " %  "
                      << (std::abs(err) < 3.0 ? "[PASS]" : "[FAIL]") << "\n";
        } else {
            std::cout << "[FEM] bar-wave verification: wave did not reach the gauge [FAIL]\n";
        }
        return;
    }
    double Es = erodedVol_ > 0 ? work_ / erodedVol_ : 0.0;
    std::cout << "[FEM] peak tool force   : " << peakF_ << " N/m\n"
              << "[FEM] tool work output  : " << work_ << " J/m";
    if (tool_.motion == Tool::Motion::FREE)
        std::cout << "  (tool KE loss: " << toolKE0_ - tool_.ke() << " J/m)";
    std::cout << "\n"
              << "[FEM] eroded volume     : " << erodedVol_ << " m^3/m ("
              << nEroded_ << " elements)\n"
              << "[FEM] specific energy   : " << Es << " J/m^3\n";
}

} // namespace rockim
