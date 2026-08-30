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
        // Taux de deformation FILTRE de l element [1/s], pour le DIF de
        // Yang et al. 2025. Miroir du champ edot du solveur 2D. La mesure
        // est la principale MAXIMALE en valeur absolue du tenseur taux :
        // les courbes DIF de la litterature sont mesurees en essai
        // UNIAXIAL, ou la principale max EST le taux axial, sans facteur
        // de convention. Vaut 0 tant que ni la viscosite ni le DIF ne sont
        // armes (la branche n est pas calculee).
        double edot = 0.0;
        int phase = 0;                     // mineral phase (index in phases_)
        int grain = 0;                     // grain id (voronoi) / 0 (grid)
        MatState st;                       // only used when law_ is set
        // global Cauchy stress R sig R^T, stored for the adaptive-insertion
        // face criterion (and cheap to keep: written once per step)
        Eigen::Matrix3d sigG = Eigen::Matrix3d::Zero();
        // WP1 pulverisation (bulkDamage = yang) : endommagement D de
        // l element et max historique du deplacement effectif h_e*eps_vm.
        // Restent a 0 (et ne sont pas ecrits) quand la cle est absente.
        double bdD = 0.0;
        double bdDm = 0.0;
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
        // ---- DIF de Yang et al. 2025, FIGE a l instant de l insertion --
        // Sortie seulement : les facteurs ont deja ete appliques a ft, coh,
        // Gf et GfII lors de activateJoint().
        double difT = 1.0, difC = 1.0, edotIns = 0.0;
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
    void applyJointSizeEffect(double mWeib);  // eq. 42 : ft *= (Zeff/V_J)^(1/m)
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

    // jointShearRange = cohesion | coulomb (defaut cohesion, bit-identique).
    // Miroir exact du 2D (FdemSolver.hpp) : la plage d'adoucissement de mode
    // II est divisee a chaque pas par fs = c + tan(phi) |sigma_n| (compression
    // seule, plancher 2 sE) — la forme publiee du modele : Guo 2014 p. 65 +
    // eq. 2.24 (fs Mohr-Coulomb), 2.30 (delta_c = 3 Gf/f) et 2.33 (driver),
    // et le code Solidity (Y3Dfd.c l. 1110-1126). La forme cohesion-seule
    // (3 GfII/c) etait une erreur de transcription, active uniquement pour
    // les joints COMPRIMES — le noyau broye d'indentation.
    bool shearRangeCoulomb_ = false;

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
    double eAbortMin_ = 0.0;   // plancher absolu [J]
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
    // --- PORTAGE 2D -> 3D du 2026-08-18 ------------------------------------
    // jointContactPenalty = fixed (defaut) | adaptive : k- = k+(D) = (1-D) pj
    //   Ghesquiere-Dierickx, Molinari & Anciaux, arXiv:2511.14323 sec. 4. La
    //   penalite de COMPRESSION suit la secante d'endommagement, supprimant le
    //   saut de raideur a dn = 0 que leur article designe comme la source
    //   dominante d'instabilite. Mesure 2D du 2026-08-18 : injection outil
    //   /14, vitesse nodale max /11, plateau de fissuration au lieu d'une
    //   inflation. Les auteurs le qualifient de DIAGNOSTIC, pas de remede :
    //   l'interpenetration croit avec D.
    // toolImpulseCap = kappa : |Fc| <= kappa * 2 * |v_outil| * m_i / dt.
    //   L'ecretage historique borne la PENETRATION (0,6 h), pas l'IMPULSION.
    //   DEFAUT CONNU, herite du 2D : kappa borne l'increment PAR PAS alors que
    //   la borne physique porte sur toute la collision — un noeud en contact
    //   soutenu accumule. Il faut une condition sur la VITESSE (CD-Lagrange).
    bool jcAdaptive_ = false;
    double toolVCap_ = 0.0;

    // --- A1 : CONTACT OUTIL EN CONDITION DE VITESSE — miroir exact du 2D ----
    // toolContact = penalty (defaut, bit-identique) | signorini.
    // Voir FdemSolver.hpp pour la formulation complete, les sources
    // (CD-Lagrange de Fekak/Brun/Gravouil ; Dureisseix et al., JTCAM 2024 ;
    // arXiv:2606.01355) et la raison pour laquelle le schema reste EXPLICITE
    // ici : masse concentree + obstacle rigide => operateur de Delassus
    // diagonal et spherique par noeud, H = 1/m_i, impulsion en forme fermee.
    // La geometrie est DUPLIQUEE dans une lambda distincte plutot que
    // factorisee : la voie penalite ne doit pas etre touchee.
    bool toolSig_ = false;
    double toolSigRelax_ = 0.0;

    // --- TRI DES FRAGMENTS (Yang et al. 2025, IJRMMS 191, 106125, sec. 2.3)
    // Porte du 2D le 2026-08-18. Corrige un defaut de computeFragments(), qui
    // declare detache TOUT ce qui n'est pas la plus grosse composante : un bloc
    // topologiquement libre peut etre geometriquement encastre entre des blocs
    // intacts. Le remede reproduit le geste experimental — on passe un pinceau
    // sur le cratere pour ramasser les debris et les peser :
    //
    //   1. le chargement doit etre TERMINE (leur etape 1) ;
    //   2. on donne aux candidats une vitesse et une acceleration OPPOSEES a la
    //      direction d'impact, choisies pour ne creer AUCUNE fissure ;
    //   3. est un debris ce qui se deplace d'au moins beta fois ce que ferait
    //      une PARTICULE LIBRE soumise aux memes conditions : d_ref = v0 t +
    //      a t^2 / 2. Ce ne sont ni un seuil ni un critere qui tranchent, c'est
    //      le contact et le frottement.
    //
    // POURQUOI LE 3D EST LE BON TERRAIN. Leur etape 1 est gratuite pour une
    // PERCUSSION — l'insert frappe, rebondit, s'en va — et impossible pour une
    // coupe continue. Mesure du 2026-08-18 : vitesse residuelle des candidats
    // 2,5 mm/s apres impact 2D, contre 108 m/s en coupe, soit un facteur
    // 43 000. Le tri a donc pleinement son sens sur le banc percussion 3D.
    //
    // RIEN N'EST SUPPRIME : c'est une CLASSIFICATION. Le travail du tri va dans
    // un poste d'energie SEPARE (brushWork_) et n'entre jamais dans le bilan
    // B4 — apres la journee du 2026-08-18 passee a traquer des pompes logees
    // dans des canaux comptabilises, en fabriquer une serait l'erreur a ne pas
    // commettre.
    //
    // Cles : gravity, toolStop, fragBrushStart, fragBrushV0, fragBrushAccel,
    //        fragBrushBeta, fragBrushDirX/Y/Z, fragBrushZeroV.
    double gravity_ = 0.0;                 // force volumique, agit selon -z
    double toolStop_ = 0.0;                // arret de l'outil (0 = jamais)
    bool toolStopped_ = false;
    double brushStart_ = 0.0;              // armement du tri (0 = desarme)
    double brushV0_ = 2.5e-3;
    double brushA_ = 98.1;                 // 10 g, leur valeur
    double brushBeta_ = 0.8;               // leur valeur (plateau de leur fig. 4)
    Eigen::Vector3d brushDir_ = Eigen::Vector3d(0.0, 0.0, 1.0);
    bool brushZeroV_ = false;
    bool brushArmed_ = false;
    double brushT0_ = 0.0;
    double brushWork_ = 0.0;               // poste SEPARE (cf. energyBodyForces)
    // ---- energyBodyForces : les deux forces VOLUMIQUES dans le bilan B4 ----
    // Contre-audit du 2026-08-30. Le theoreme travail-energie ne distingue pas
    // une force « physique » d'une force « numerique » : TOUTE force appliquee
    // aux noeuds fait un travail, et ce travail est soit dans sumW, soit dans
    // le RESIDU. Deux y manquaient :
    //   * la PESANTEUR. Aucun compteur n'existait, alors que gravity = 9.81 est
    //     pose dans 20 des 22 decks de bench_impact/configs, et qu'ARMA 24-0952
    //     pose explicitement l'energie potentielle gravitaire dans ses eq. 3-7.
    //     C'est le SEPTIEME poste de leur bilan, et le seul que rockim n'avait
    //     pas. Magnitude sur un impact : deplacements micrometriques sur 600 us,
    //     donc ~1e-4 J contre ~49 J injectes — INVISIBLE. Le defaut est
    //     STRUCTUREL, pas numerique : il devient reel des qu'un run est long ou
    //     quasi statique.
    //   * le TRI des fragments (brushWork_), tenu hors bilan a dessein — la
    //     raison ecrite plus haut (« ne pas fabriquer une pompe logee dans un
    //     canal comptabilise ») est bonne, mais elle ne protege pas de ce
    //     qu'elle craint : hors de sumW, ce travail tombe entierement dans le
    //     residu, ou budgetAbortPct peut couper un run SAIN sur un artefact
    //     purement numerique. Le mesurer et le montrer protege ; le cacher, non.
    //
    // Le compromis retenu : la MESURE est inconditionnelle et IMPRIMEE (elle ne
    // touche aucune force, donc la physique reste bit-identique dans les deux
    // cas) ; seule l'ENTREE DANS sumW est opt-in.
    //   off (defaut, bit-identique) : sumW inchange, et le resume DIT en toutes
    //                                 lettres combien de J tombent dans le residu.
    //   on                          : gravWork_ et brushWork_ entrent dans sumW
    //                                 et dans l'echelle, donc dans le verdict
    //                                 [OK|CHECK] et dans budgetAbortPct.
    double gravWork_ = 0.0;                // travail de la PESANTEUR (B4)
    bool   eBody_ = false;                 // energyBodyForces = on
    std::vector<char> brushCand_;
    std::vector<int> brushFrag_;
    int brushNFrag_ = 0;
    std::vector<Eigen::Vector3d> brushU0_;
    void bodyForces();
    void armBrush();
    void brushReport();

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

    // ---- viscosite newtonienne de Yan et al. 2023 (leur eq. 6, 2 mu D) --
    // muEl_ porte mu PAR ELEMENT. En mode global il est rempli d une seule
    // valeur ; en mode gradue (bulkViscosityGraded) il suit la taille de
    // maille. La distinction n est pas cosmetique : la borne diffusive de
    // stabilite se cale sur le PLUS PETIT element, donc un mu global sur un
    // maillage gradue 17x paie le pas de temps du plus fin tetra partout.
    // Chiffrage sur p1_ultra (127 147 tetras, h de 0,109 a 2,92 mm) :
    // xi = 2 global divise le pas par 41,8 (2 h -> 3,8 JOURS) contre 2,35
    // en gradue. Le mode gradue est le seul finançable en percussion.
    std::vector<double> muEl_;
    bool viscOn_ = false;                  // un mu > 0 quelque part
    // viscousInInsertion : la contrainte visqueuse entre-t-elle dans sigG,
    // donc dans le critere d insertion ? Defaut 1 = OUI, fidele a Yan (leur
    // eq. 6 EST la contrainte que lit leur eq. 7) et parite avec le 2D.
    // Mis a 0, le visqueux ne sert plus qu aux forces : le critere, la
    // jauge de confinement et le stocke elastique du bilan redeviennent
    // purement elastiques. A considerer avant toute calibration de ft sur
    // un run visqueux, et surtout si le DIF est arme en meme temps — sinon
    // le taux de deformation agit DEUX FOIS, une fois en gonflant la
    // contrainte d essai et une fois en gonflant le seuil.
    bool viscIns_ = true;
    double viscWork_ = 0.0;                // ventilation, DEJA dans elWork_
    // WP1 pulverisation (Yang et al. 2026). bdWork_ = ventilation de la
    // dissipation d endommagement volumique, DEJA comptee dans elWork_.
    bool bdOn_ = false;
    double bdD0_ = 1.4e-5, bdDf_ = 4.0e-4, bdDmax_ = 0.9, bdCd_ = 1.0;
    // WP6 (spec 005, plan WP6_contact_residuel.md, 2026-08-28) : mu de
    // contact RESIDUEL post-pulverisation — l ingredient "mobilite" du
    // modele de Yang et al. 2026 (sliding friction 0,18 sur le granite,
    // contre tan(phi) = 1,85 intact). Quand une interaction de contact
    // implique un element PULVERISE (bulkDamage : D arrive a Dmax), son mu
    // passe de contactMu a cette valeur — bascule BINAIRE au franchissement
    // de Dmax, comme le relais joint mort -> contact, car le papier
    // l applique au materiau ROMPU (une rampe en D degraderait la force de
    // penetration en pleine charge, delta_0 = 0,014 mm etant franchi
    // presque immediatement sous l insert). Sites : contact OUTIL (les
    // decisifs — il n existe pas de joint outil-roche, et sous l insert en
    // jointDeath = separation le relais contact ne s engage jamais) ET
    // contact general (ejection des debris). HORS perimetre : les plateaux
    // de compression (frontiere de machine, pas un support de fragments).
    // muCRes_ < 0 = cle absente = comportement historique, bit-identique.
    double muCRes_ = -1.0;
    unsigned long long nCtcPulv_ = 0;   // evaluations au mu residuel
    double tCtcPulv0_ = -1.0;           // premier engagement (s)
    // ---- frottement PAR PHASE (Table 1 de Yang et al. 2026 : le
    // coefficient glissant est une propriete MATERIAU — 0,18 granite,
    // 0,6 carbure, 0,6 acier). Cle : contactMu.<nom de phase>. Regle de
    // paire = le MINIMUM des deux (Solidity Y3Did.c l. 1292 :
    // if(d1pefr[iprop]>d1pefr[jprop]) iprop=jprop). Vide = contactMu
    // global partout, bit-identique.
    std::vector<double> muPhase_;
    bool muPerPhase_ = false;
    // ---- couplage CONTINU endommagement -> contact (Solidity Y3Did.c :
    // d_fact = min(1-D_i, 1-D_j), effondrement /1000 sous 0,041 l. 1264,
    // penalty *= d_fact l. 1265, mu = mud*d_fact l. 995 et 1044).
    // contactDamageCoupling = solidity ; 0 = off = bit-identique. La
    // raideur NORMALE et le frottement suivent (1-D) des que D > 0 —
    // c est l effondrement de portance que WP6 (echelon binaire sur le
    // seul mu) ne fournissait pas.
    int cplMode_ = 0;
    unsigned long long nCplEval_ = 0;   // evaluations avec d_fact < 1
    unsigned long long nCplColl_ = 0;   // effondrements (d_fact < 0,041)
    double tCpl0_ = -1.0;               // premier engagement (s)
    inline double cplDf(int eA, int eB) const {
        double d = 1.0;
        if (eA >= 0) d = std::min(d, 1.0 - el_[eA].bdD);
        if (eB >= 0) d = std::min(d, 1.0 - el_[eB].bdD);
        return (d < 0.041) ? d * 1e-3 : d;   // Y3Did.c l. 1264
    }
    // mu effectif d une interaction impliquant les elements eA (et eB si
    // >= 0). Compteurs mis a jour au premier engagement et a chaque
    // evaluation degradee (atomic : appele depuis les boucles OMP).
    inline double ctcMu(int eA, int eB = -1) {
        double mu = muC_;
        if (muPerPhase_) {               // le plus FAIBLE gouverne la paire
            double m = 1e300;
            if (eA >= 0) m = std::min(m, muPhase_[el_[eA].phase]);
            if (eB >= 0) m = std::min(m, muPhase_[el_[eB].phase]);
            if (m < 1e299) mu = m;
        }
        if (cplMode_) {
            double dr = 1.0;
            if (eA >= 0) dr = std::min(dr, 1.0 - el_[eA].bdD);
            if (eB >= 0) dr = std::min(dr, 1.0 - el_[eB].bdD);
            if (dr < 1.0) {
#ifdef _OPENMP
#pragma omp atomic
#endif
                ++nCplEval_;
                if (dr < 0.041) {
#ifdef _OPENMP
#pragma omp atomic
#endif
                    ++nCplColl_;
                    dr *= 1e-3;
                }
                if (tCpl0_ < 0.0) tCpl0_ = t_;   // course benigne
                mu *= dr;
            }
            return mu;
        }
        if (muCRes_ < 0.0) return mu;
        bool p = (eA >= 0 && el_[eA].bdD >= bdDmax_)
              || (eB >= 0 && el_[eB].bdD >= bdDmax_);
        if (!p) return mu;
#ifdef _OPENMP
#pragma omp atomic
#endif
        ++nCtcPulv_;
        if (tCtcPulv0_ < 0.0) tCtcPulv0_ = t_;   // course benigne : ~meme t
        return muCRes_;
    }
    double bdWork_ = 0.0;
    long nPulv_ = 0;                       // elements a D = Dmax (lecture)

    // ---- DIF de Yang et al. 2025 (leurs eq. 2-3) ------------------------
    // Voir include/rockim/YangDif.hpp pour les formules et la coquille de
    // l exposant. FIGE a l insertion : le schema extrinseque insere le
    // joint a l instant precis ou le materiau atteint son enveloppe, donc
    // c est l instant dont le taux gouverne la rupture ; et un DIF
    // reversible sur un joint deja endommage demanderait un cliquet.
    // ---- enveloppe de cisaillement : Yan (eq. 8) ou Yang (eq. 1) --------
    // Voir include/rockim/YangDif.hpp. Defaut yan = comportement historique,
    // bit-identique. La forme de Yang affaiblit le cisaillement dans les
    // zones TENDUES (elle y fait decroitre le terme de frottement au lieu de
    // l annuler), ce qui deplace le partage traction / cisaillement — donc le
    // facies. Les deux formes sont rigoureusement identiques en compression.
    bool yangEnv_ = false;
    // Facteur du cap de traction moyenne (pm <= f * ft). Historiquement 3, en
    // dur et sans echappatoire en 3D. Ce cap N EST DANS AUCUN des deux
    // articles : il faut pouvoir le desarmer pour faire tourner le modele de
    // quelqu un d autre. <= 0 le desarme. Defaut 3 = inchange.
    double mtCap_ = 3.0;
    bool difOn_ = false;
    double difExpT_ = 0.07;                // 0,07 litteral | 0,1707 fig. 2b
    double srTau_ = 0.0, srRelax_ = 0.0;   // filtre du taux

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
    // dtBudgetTangential : la raideur tangentielle du contact par potentiel
    // entre-t-elle dans le budget de pas de temps ? Defaut false =
    // bit-identique. Xiang, Munjiza, Latham & Guises, Eng. Comput. 26(6)
    // (2009) 673-687, p. 677 : le calcul des forces TANGENTIELLES exige un pas
    // plus petit que le cas sans frottement. Voir computeStableDt().
    bool dtTangential_ = false;
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
    double toolKEStop_ = -1.0;   // KE a l arret toolStop (cf. miroir 2D)

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
    // WP2 : matrice (nGroups x nGroups) des paires de corps LIES par
    // groupBond.<A>.<B> = joints. Vide ou 0 partout = aucun joint
    // inter-corps, comportement historique.
    std::vector<char> gbond_;
    // WP3 : cinematique par corps (trackGroups) et jauges en tranche
    // (gauge.<nom> = z0 z1). Vides = aucune colonne ajoutee.
    struct Gauge3 { int grp = -1; double z0 = 0, z1 = 0;
                    std::vector<int> elems; };
    std::vector<int> trkGrps_;
    std::vector<Gauge3> gauges_;
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
