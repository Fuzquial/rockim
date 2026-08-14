#pragma once
// ---------------------------------------------------------------------------
// Fdem3dSolver — 3D combined finite-discrete element method, extending the
// 2D FdemSolver to tetrahedra. Every safeguard learned the hard way in 2D is
// built in from the start (see the README debugging story): consistently
// outward-oriented faces, joints that die by clear separation only,
// initial-penetration relief at contact birth, and a quasi-plastic soft
// debris contact.
//
//   * Two mesh front-ends behind one topology builder, as in 2D:
//     mesh = grid    — structured hex grid split into 6 Kuhn tetrahedra per
//                      cell (compatible face diagonals), optional jitter;
//     mesh = voronoi — 3D Voronoi grain structure (Tessellation3) with
//                      per-grain mineral phases: the GBM mode. Crack paths
//                      then follow grain boundaries and intra-grain fans.
//   * Linear tets, CO-ROTATIONAL: R from F by 3 Higham iterations
//     R <- (R + R^-T)/2 (globally convergent for det F > 0), Biot strain in
//     the co-rotated frame, crush cap on the deviator + mean-tension cap
//     (per-phase in the GBM mode).
//   * Triangular 6-node cohesive joints on every interior face, 3 node-pair
//     integration points: mode I softens from ft over 2 Gf/ft of opening,
//     mode II from c over 2 Gf_II/c of frictional slip (vector return
//     mapping in the face plane), friction tan(phi)(-sigma_n) throughout.
//     Properties live PER JOINT: intra-grain joints carry the phase
//     material, grain-boundary joints the attenuated mean of the two phases
//     (type 0 = intra-grain, 1 = homophase boundary, 2 = heterophase), and
//     jointWeibullM statistical strengths (independent draws or one
//     correlated random field, RandomField3) scale ft and cohesion.
//   * General node-triangle penalty contact (cell grid, clipped box,
//     birth-gap, quasi-plastic restitution, co-location exclusion; AABB
//     binning on the voronoi mesh, whose faces outgrow the grid cells).
//   * Rigid tool: sphere (free percussion / prescribed shear) or
//     flat-ended punch (percussion), viscous-spring quiet boundaries on
//     the four lateral faces and the bottom, tension scenario with
//     pullRamp and gripLateralFree as in 2D.
// ---------------------------------------------------------------------------
#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include "rockim/Config.hpp"
#include "rockim/MatLaw.hpp"
#include "rockim/Material.hpp"
#include "rockim/Solver.hpp"
#include "rockim/YanSoftening.hpp"

namespace rockim {

class Fdem3dSolver : public Solver {
public:
    Fdem3dSolver(const Config& cfg, std::string outDir);

    void init() override;
    void step() override;
    void writeFrame(int frame) override;
    void historyHeader(std::ostream&) const override;
    void historyRow(std::ostream&) const override;
    bool finished() const override;        // E2 : moniteur d'energie
    void finalize() override;

private:
    enum Flag { FREE = 0, FIXED = 1, PRESCRIBED = 2 };
    enum class Scenario { PERCUSSION, SHEAR, TENSION };

    struct Elem {
        std::array<int, 4> n;
        Eigen::Matrix<double, 3, 4> dN;    // reference shape-fn gradients
        double V0;
        double svm = 0.0;
        int phase = 0;                     // mineral phase (index in phases_)
        int grain = 0;                     // grain id (voronoi) / 0 (grid)
        MatState st;                       // only used when law_ is set
        // global Cauchy stress R sig R^T, stored for the adaptive-insertion
        // face criterion (and cheap to keep: written once per step)
        Eigen::Matrix3d sigG = Eigen::Matrix3d::Zero();
    };

    // Triangular cohesive joint. Properties live per joint (GBM): type 0 =
    // intra-grain (bulk of the phase), 1 = homophase grain boundary, 2 =
    // heterophase; stat is the Weibull strength factor (output as ftScale).
    struct Joint {
        int eA, eB;
        std::array<int, 3> a, b;            // node pairs (a[k] ~ b[k])
        double A0;                          // face area
        double D = 0.0;
        std::array<Eigen::Vector3d, 3> slip;  // jointShearUnload = plastic : le
                                           // glissement plastique cumule par le
                                           // retour radial vectoriel.
                                           // = origin : l'ORIGINE FIGEE de la
                                           // secante de l'eq. 18, stampee par
                                           // activateJoint() (0 en intrinseque),
                                           // seulement reprojetee dans le plan.
        bool dead = false;
        double tBreak = -1.0;
        int type = 0;
        double pj = 0.0;                    // penalty per area [Pa/m]
        double ft = 0.0, coh = 0.0;         // strengths [Pa]
        double Gf = 0.0, GfII = 0.0;        // fracture energies [J/m^2]
        double dnE = 0.0, dnF = 0.0;        // mode I elastic / final opening
        double slipF = 0.0;                 // mode II softening slip
        double tanPhi = 0.0;                // friction
        double stat = 1.0;                  // Weibull strength factor (output)
        // ---- adaptive insertion (Yan et al. 2023, ported from FdemSolver) --
        // bonded = the joint does not exist yet: the co-located node copies
        // are rigidly bound (exact shared-node FEM) and processJoint skips
        // it. dn0 is the opening offset stamped at activation for stress
        // continuity (the newborn joint transmits at zero geometric opening
        // exactly the traction the bonded face was carrying).
        bool bonded = false;
        double dn0 = 0.0;
        // largest opening ever reached per integration point (jointSoftening
        // = yan: the omax of eq. 17 — origin-secant unloading), as in 2D
        std::array<double, 3> omax{{0.0, 0.0, 0.0}};
        // largest SLIDING ever reached, the s_max of eq. 18, mesure depuis
        // l'origine figee slip[k] (jointShearUnload = origin seulement)
        std::array<double, 3> smax{{0.0, 0.0, 0.0}};
        // ---- failure mode, recorded ONCE when D first reaches 1 -----------
        // Ported from the 2D solver: partition of the eq. 16 damage driver at
        // the breaking instant. bmode = 1 tensile, 2 shear, 0 intact;
        // failMode = rn^2/(rn^2+rs^2) in [0,1], -1 while intact. Output only.
        int bmode = 0;
        double rnB = 0.0, rsB = 0.0;
        double failMode = -1.0;
    };

    struct BFace {                          // active contact face
        int elem;
        std::array<int, 3> n;
    };

    struct Tool3 {
        bool free = true;
        bool flat = false;   // flat-ended cylindrical punch (axis z); x is
                             // then the center of the bottom face — the 3D
                             // lift of the 2D FLAT tool (percussion only)
        double mass = 0.5, radius = 0.015;
        Eigen::Vector3d x{0, 0, 0}, v{0, 0, 0}, F{0, 0, 0};
        void integrate(double dt) { if (free) v += (dt / mass) * F; x += dt * v; }
        double ke() const { return 0.5 * mass * v.squaredNorm(); }
    };

    void buildMesh();                      // dispatch: grid | voronoi | file
    void buildMeshVoronoi();
    // mesh = file — import d'un maillage tetraedrique NON STRUCTURE (Gmsh
    // MSH 2.2 ASCII, elements type 4), la voie "maillage a la Yan et al." :
    // simplexes uniformes sans structure de grains. Le bloc est translate a
    // l'origine et W/D/H sont relus de la boite englobante.
    void buildMeshFile();
    void buildFromTets(const std::vector<Eigen::Vector3d>& vpos,
                       const std::vector<std::array<int, 4>>& tets,
                       const std::vector<int>& tetGrain,
                       const std::vector<int>& grainPhase);
    void assignJointProps();
    void applyJointStatistics();
    // adaptive insertion machinery (insertion = adaptive), as in 2D
    void buildBindingTables();
    void rebindVertex(int v);
    void insertionSweep();
    void activateJoint(int jI, double sig, const Eigen::Vector3d& tauV,
                       double fsNow);
    void placeTool();
    void setupBoundaries();
    // triaxial 3D : pression suiveuse sur les faces exterieures LATERALES
    // d'origine (la membrane de la cellule), rampe cosinus, jauge de la
    // contrainte laterale atteinte dans le coeur — le portage direct du
    // confinement 2D (FdemSolver), scenario = tension + pullV < 0.
    void setupConfinement();
    void confiningForces();
    double achievedConfinement() const;
    void computeStableDt();

    void elementForces();
    void jointForces();
    void rebuildContactFaces();
    void generalContact();
    void toolContact();
    void integrate();
    void computeFragments();

    Config cfg_;
    std::string out_;
    Material mat_;
    PhaseSet phases_;                      // per-phase materials (>= 1)
    Scenario scen_ = Scenario::PERCUSSION;

    double W_ = 0.08, D_ = 0.08, H_ = 0.06, hmin_ = 4e-3;
    int nx_ = 20, ny_ = 20, nz_ = 15;
    bool voronoi_ = false;
    int nGrains_ = 1;

    // joint law: per-joint properties live in Joint; xiJ_ is the shared
    // dashpot ratio (bilateral on intact joints, clipped resultant on broken
    // ones — the 2026-08-05 2D refonte, ported)
    double xiJ_ = 0.05;

    // ---- post-peak softening shape (jointSoftening), ported from 2D --------
    // linear (default, unchanged): ft -> 0 over dnF - dnE = 2 Gf / ft.
    // yan: the exponential reduction factor f(D) of Yan et al. (IJRMMS 169,
    // 2023) eq. 11, critical opening/slip calibrated so the area under the
    // softening branch is still exactly GfI / GfII:
    //     dnF - dnE = Gf / (ft yanI_),  slipF = GfII / (c yanI_).
    bool yanSoft_ = false;
    bool yanFricScaled_ = false;
    yan::Params yanP_;
    double yanI_ = 1.0;                    // int_0^1 f(D) dD

    // ---- decharge en CISAILLEMENT (jointShearUnload), miroir exact du 2D ----
    // plastic (defaut, inchange) : plasticite a retour vectorielle — la
    //   decharge suit la secante de PENALITE, le glissement plastique reste.
    // origin : l'eq. 18 de Yan et al., symetrique de l'eq. 17 du mode I — la
    //   decharge ET la recharge suivent la secante a l'ORIGINE passant par
    //   (s_max, tau_env(s_max)), avec s_p = (c + tan(phi)|sigma_n|)/p le
    //   glissement au pic de Munjiza. Forme litterale de l'article.
    // ⚠️ l'eq. 18 place TOUT le cap dans la secante, frottement de Coulomb
    // compris : combinee a jointFrictionScaled = 0 elle rend le glissement
    // frottant reversible. Forme litterale = origin + jointFrictionScaled = 1.
    bool shearOrigin_ = false;

    // Work done by the joint dashpot. A dashpot can only DISSIPATE, so this
    // must stay <= 0 — the direct detector of the rectifier/anti-damping
    // failure modes (2D lesson, 2026-08-05). One multiply per point.
    double dampWork_ = 0.0;

    // ---- V2/B4 : bilan d'energie par sous-systeme (instrumentation PURE,
    // aucun flottant de la physique ne change ; biais O(dt) des compteurs
    // travail comme gcWork_, documente). Theoreme travail-energie :
    //   KE(t) - KE(0) = elWork_ + jointWork_ + gcWork_ + cundWork_
    //                 + lysWork_ + toolWork_ + bcWork_ + residu(O(dt))
    // Postes lisibles : -elWork_ = elastique stocke + caps ; -(jointWork_ -
    // dampWork_) = cohesif (fissuration + stocke transitoire) ; -cundWork_ =
    // amortissement local Cundall ; -lysWork_ = amortisseurs de Lysmer ;
    // gcFricWork_ = part frottement de gcWork_ ; toolWork_ = travail de
    // l'outil rigide SUR le solide ; bcWork_ = travail des noeuds PRESCRIBED
    // (platines). Les compteurs multi-threads se reduisent en ordre de
    // thread, meme statut de reproductibilite que dampWork_.
    double elWork_ = 0.0;      // forces internes des elements
    double jointWork_ = 0.0;   // tractions des joints (visqueux INCLUS)
    double cundWork_ = 0.0;    // damping local de Cundall (<= 0)
    double lysWork_ = 0.0;     // amortisseurs frontiere (<= 0)
    double gcFricWork_ = 0.0;  // part tangentielle du contact general
    double toolWork_ = 0.0;    // outil rigide -> solide
    double bcWork_ = 0.0;      // platines (PRESCRIBED) -> solide
    double confWork_ = 0.0;    // pression de confinement/bore -> solide
    // E2 : moniteur d'energie runtime (opt-in budgetAbortPct, 0 = off)
    double eAbortPct_ = -1.0;  // -1 = pas encore lu dans la config
    bool eAbort_ = false;
    void checkEnergyAbort();
    // arret post-rupture en tension/compression (opt-in stopPeakDrop)
    double stopDrop_ = -1.0;
    bool peakStop_ = false;
    // V2/B2 : force de contact NETTE sur le corps suivi (trackGroup) au pas
    // courant — remise a zero en tete de generalContact, sommee dans les
    // deux lois (penalite et potentiel). La F-delta se lit alors en direct
    // dans history (grpFx/y/z) au lieu de M dv/dt.
    Eigen::Vector3d grpF_ = Eigen::Vector3d::Zero();
    double biasW_ = 0.0;       // correction leapfrog EXACTE : les compteurs
                               // lisent v- ; le theoreme discret veut
                               // (v- + v+)/2 -> ecart = f_tot^2 dt^2 / 2m par
                               // noeud et par pas, accumule ici (>= 0)
    double keInit_ = -1.0;     // KE du solide au premier pas (< 0 = non prise)

    // ---- adaptive insertion (insertion = adaptive), as in FdemSolver -------
    // No joint exists at t = 0: bonded faces are handled kinematically (node
    // groups = connected components of the tet fan around each original
    // vertex over still-bonded faces), the sweep averages the two elements'
    // global stress on every bonded face and activates the joint when
    // sigma_n >= ft or |tau| >= fs (Mohr-Coulomb). Node splitting falls out
    // of re-running the union-find at the face's three vertices.
    bool adaptive_ = false;
    long nInserted_ = 0;
    std::vector<std::vector<int>> copiesOfVert_;
    std::vector<std::vector<int>> jointsOfVert_;
    std::vector<std::vector<std::vector<int>>> grpsOfVert_;
    int nVert_ = 0;
    // per-phase element tables (indexed by Elem::phase)
    std::vector<double> lamP_, mu2P_, crushCapP_, ftP_, rhoP_;
    // Optional BULK constitutive law (law = dpr | saksala | saksala2011 |
    // dpdfh) — the coupled configuration of the thesis' vumat_fdem_coupled:
    // a dissipative bulk PLUS discrete cohesive joints. Simpler than in 2D:
    // the element loop already works on 3x3 Biot strain, so the law plugs in
    // directly with no plane-strain embedding. Absent by default, so every
    // earlier result stays reproducible.
    std::unique_ptr<MatLaw> law_;
    std::vector<double> hEl_;              // per-element inscribed size 6V/A

    double kp_ = 0.0, muC_ = 0.5, xiC_ = 0.05, vReg_ = 1e-3;
    double kpGC_ = 0.0, xiGC_ = 0.8, gcRest_ = 0.2, gcWork_ = 0.0, relax_ = 1.0;
    std::unordered_map<uint64_t, double> pen0_;

    double damping_ = 0.05, pullV_ = 0.05;
    double pullRamp_ = 0.0;                // grip velocity rise time [s]
    bool gripFree_ = false;                // frictionless tension grips

    std::vector<Eigen::Vector3d> X0_, u_, v_, f_;
    std::vector<double> m_;
    std::vector<int> flag_, elemOf_, vOf_;
    std::vector<Elem> el_;
    std::vector<Joint> jt_;
    std::vector<int> fragId_;
    std::vector<BFace> exterior_;

    std::vector<Eigen::Vector3d> cAbs_, kAbs_;

    std::vector<BFace> act_;
    std::vector<int> actNodes_;
    // OpenMP scratch: per-thread force buffers for the joint scatter,
    // reduced serially in thread order (deterministic for a fixed thread
    // count; 1 thread = bit-identical to the serial build), as in 2D
    std::vector<std::vector<Eigen::Vector3d>> fTL_;
    std::vector<std::vector<char>> seenTL_;
    std::vector<std::vector<int>> touchedTL_;
    long actStamp_ = -1;
    double cell_ = 0.0;
    double cellV_ = 0.0;                   // voronoi contact cell size (2 x
                                           // median hEl: 2 x hmin would put
                                           // ~(L/hmin)^3 cells on the box —
                                           // the dense grid is REALLOCATED
                                           // every step, ~0.5 GB/step on the
                                           // percussion demo)
    Eigen::Vector3d gmin_;
    int gx_ = 1, gy_ = 1, gz_ = 1;
    std::vector<std::vector<int>> grid_;   // dense grid (grid mesh)
    std::unordered_map<uint64_t, std::vector<int>> gridV_;  // sparse (voronoi)

    // ---- activation adaptative du contact (gcActivation = adaptive) --------
    // Miroir exact du 2D (voir FdemSolver.hpp pour les trois regles C/A/B, la
    // cadence par v_max et l'approximation assumee). Pool = tout l'exterieur
    // (pas de gcXwindow en 3D). Defaut full, bit-identique.
    bool gcAdaptive_ = false;
    double gcActMargin_ = 2.0;             // marge d'activation [cellules]
    long gcActEvery_ = 64;                 // cadence max du balayage [pas]
    std::vector<BFace> pool_;              // exterieur (fige a l'init)
    std::vector<char> extOn_;              // drapeau actif, par face du pool
    std::vector<char> elemDam_;            // regle C : element au bord casse
                                           // + un anneau par sommet
    std::vector<std::vector<int>> vElems_; // sommet -> elements (statique)
    std::vector<int> bodyOf_;              // composante connexe par element
    std::vector<long> lastTouch_;          // dernier pas de contact du noeud
    long bodyStamp_ = -1;                  // nBroken_ au dernier union-find
    int nBodies_ = 1;
    long nextSweep_ = 0;
    long sweepBroken_ = -1;
    bool poolBuilt_ = false;
    bool haveDead_ = false;
    long nActivated_ = 0;
    std::vector<BFace> deadList_;          // cache des faces liberees (timing
                                           // du declencheur historique)
    void activationSweep();

    // ---- contact par POTENTIEL de Munjiza (contact = potential) — A3 -------
    // Miroir exact du 2D (FdemSolver.hpp) en tet-tet : force distribuee
    // p ∮ (phi_A - phi_B) n dG sur le bord du POLYEDRE de recouvrement
    // (PotentialContact.hpp, pot3), frottement incremental vectoriel de
    // l'eq. 4-5, detection element-element O(N), exclusion des paires liees
    // par un joint vivant, releve de naissance par VOLUME (le pen0_ du
    // potentiel — voir la lecon 2D : une rampe temporelle INJECTE).
    bool contactPot_ = false;
    double potP_ = 0.0;                    // penalite normale [Pa]
    double potKt_ = 0.0;                   // raideur tangentielle [N/m]
    struct PotHist {
        Eigen::Vector3d Ft{0.0, 0.0, 0.0}; // force tangentielle sur l'el. MIN
        long step = -1000;
        double vRef = -1.0;                // reference de naissance (volume) ;
                                           // < 0 = jamais recouvert (l'entree
                                           // peut preexister via le cache
                                           // d'axe separateur)
        int sepAxis = -1;                  // plan separateur en cache (Baraff)
    };
    std::unordered_map<uint64_t, PotHist> potFt_;
    std::unordered_map<uint64_t, int> jointOfPair_;
    // Compteurs DIAGNOSTIQUES du potentiel (issues par paire + chronos de
    // sections), cumules sur le run et imprimes au resume. Aucune influence
    // sur la physique — pur observatoire pour piloter l'optimisation N1.
    struct PotStats {
        uint64_t pairs = 0;      // candidates AABB (apres dedup, hors joint)
        uint64_t joint = 0;      // ecartees : joint vivant porte la paire
        uint64_t sepHint = 0;    // separees par l'axe en cache (1 test)
        uint64_t sepFace = 0;    // separees par scan, axe de face (0-7)
        uint64_t sepEdge = 0;    // separees par scan, axe d'arete (8-43)
        uint64_t clipMiss = 0;   // clip complet -> pas de recouvrement reel
        uint64_t clipHit = 0;    // clip complet -> force appliquee
        double tGrid = 0.0;      // s : grille + collecte + tri des paires
        double tLoop = 0.0;      // s : boucle des paires (SAT + clips)
    } potStats_;
    void potentialContact();

    Tool3 tool_;
    double toolKE0_ = 0.0;

    // confinement triaxial (0 = inactif)
    double confP_ = 0.0, confRamp_ = 0.0;
    double pullDelay_ = 0.0;               // equilibrage de sigma3 avant l'axial
    std::vector<BFace> confFaces_;         // faces LATERALES d'origine seulement
    bool confLatched_ = false;
    double confAchieved_ = 0.0;

    // ---- physical groups Gmsh (mesh = file) — V1 ---------------------------
    // $PhysicalNames (dim 3) -> un GROUPE par volume physique : materiau par
    // groupe (phase du meme nom, ou groupPhase.<nom>), PAS de joint entre
    // groupes (les deux faces vont a exterior_ : l'interface est un contact,
    // penalite ou potentiel) — c'est le multi-corps piston-taillant-roche.
    // groupVel.<nom> = "vx vy vz" donne une vitesse initiale au groupe :
    // l'outil MAILLE remplace l'outil analytique (toolShape = none).
    std::vector<int> elemGroup_;           // groupe par element (vide = mono)
    std::vector<int> tetGroupTmp_;         // par tet, pour buildFromTets
    std::vector<std::string> groupName_;
    int nGroups_ = 1;
    bool toolNone_ = false;                // toolShape = none : outil maille
    int trackGroup_ = -1;                  // groupe suivi dans history.csv

    long stepCount_ = 0;
    double work_ = 0.0, peakF_ = 0.0;
    long nBroken_ = 0;
    int nFrag_ = 1;
    double detachedVol_ = 0.0;
    Eigen::Vector3d gripF_{0, 0, 0};
    double sigmaPeak_ = 0.0;
};

} // namespace rockim
