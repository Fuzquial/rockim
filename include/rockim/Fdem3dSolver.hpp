#pragma once
//
// ---------------------------------------------------------------------------
// PROVENANCE DES CITATIONS « Y3D*.c l. NNNN » DE CE FICHIER — lire d abord
//   ../../SOURCES_SOLIDITY.md  (a la racine du depot, note B4a du 2026-08-30)
//
// Elles renvoient au code d Imperial College London,
// github.com/ImperialCollegeLondon/solidity-solver-open, LGPL-3.0, LU LE
// 2026-08-26. Trois choses a savoir avant d en citer une :
//   1. c est BIEN leur code — le contraire a ete affirme puis rectifie ;
//   2. ce n est PAS la version qui a produit l article de 2026 (facteur
//      d endommagement cable a zero, DIF neutre) : y lire une FORME et en
//      conclure une implementation de l article est une faute ;
//   3. LES NUMEROS DE LIGNE NE SONT PAS ANCRES SUR UN COMMIT. Le depot est
//      maintenu, donc ils bougent. Ils valent pour le 2026-08-26.
// Les 72 citations du depot ne visent que 13 endroits distincts : la table
// des 13, avec leur statut (article / code public / version interne), est
// dans SOURCES_SOLIDITY.md §3.
// ---------------------------------------------------------------------------
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
#include "rockim/JointTsl.hpp"
#include "rockim/MatLaw.hpp"
#include "rockim/Material.hpp"
#include "rockim/Solver.hpp"
#include "rockim/YanSoftening.hpp"

namespace rockim {

// ---------------------------------------------------------------------------
//  CURSEUR DE COULOMB D UN JOINT INSERE SOUS jointTSL = camacho
//  (RETOUR_v3_2026-09-11 §1.4 quater ; note v3 §3.4)
//
//  Le frottement d un joint insere est une traction d ESSAI de collage sur
//  le glissement GEOMETRIQUE, ecretee au cap de Coulomb, avec retour de
//  glissement plastique — le mecanisme de la branche penalisee de g0
//  (Fdem3dSolver.cpp, retour radial vectoriel sur J.slip[k]) applique a la
//  SEULE part frottante :
//
//      tau_tr   = pj (delta_s - slip_p)
//      f_cap    = [D] mu <-t_n>              (jtsl::shearCap a cohesion nulle)
//      |tau_tr| <= f_cap : tau_fric = tau_tr
//      |tau_tr| >  f_cap : tau_fric = tau_tr f_cap/|tau_tr|,
//                          slip_p  += (tau_tr - tau_fric)/pj
//      f_cap == 0        : tau_fric = 0 ET slip_p = delta_s
//
//  Le dernier cas (traction, ou D = 0 sous jointFrictionMobilised) est
//  essentiel : sans lui le ressort de collage pj ajouterait en traction une
//  raideur tangentielle que la loi cohesive de l eq. 19 ne prevoit pas.
//
//  POURQUOI cette forme et pas l addition. La forme livree le 2026-09-11
//  appliquait [D] mu <-t_n> comme une TRACTION de magnitude constante dirigee
//  par le DEPLACEMENT de glissement effectif (tau = (lim/dsEffN) dsEffV) :
//  un ressort sec a force constante, discontinu a l origine, sur chaque joint
//  insere comprime — et sur chaque joint ROMPU comprime, que jointDeath =
//  separation ne tue jamais. Cinq isolations (RETOUR_v3 §1.4 ter/quater) :
//  explosion a ~91 us quelle que soit la matrice (elastique comprise), quel
//  que soit dt, avec ou sans branche ascendante ; x25 insertions en cascade ;
//  frottement mobilise eteint => explose PLUS TOT ; joints penalises de g0
//  (frottement en CAP avec retour de glissement) => 100 us sains. La v3 §3.4
//  ecrit la bonne forme, t_s^fric = -mu <-t_n> d(delta_s)/dt/||d(delta_s)/dt||
//  regularise : une force qui s oppose a la VITESSE, pas au deplacement. Le
//  curseur elastique-parfaitement-plastique ci-dessus en est la
//  regularisation standard (raideur pj, deja au budget CFL : max(k0, kPara
//  pj) sur toutes les facettes, computeStableDt).
//
//  Proprietes verifiees par tests_f2/check_camacho_friction3d.cpp :
//   (i)   trajet en traction pure : tau_fric = 0 identiquement ;
//   (ii)  cycle de glissement alterne sous compression : dissipation par
//         cycle = 2 f_cap (P - 2 f_cap/pj) >= 0 avec P l amplitude crete a
//         crete (-> 2 f_cap P quand pj -> inf), travail cumule jamais sous
//         -f_cap^2/(2 pj) (l energie de collage stockee) ;
//   (iii) continuite en delta_s = 0 (pente pj, contre un saut de 2 f_cap
//         pour la forme retiree).
//
//  Fonction LIBRE et sans etat pour etre banchee sans le solveur. Le glissement
//  slip est maintenu dans le plan de la facette par l appelant, comme dans la
//  branche historique. Retourne tau_fric.
// ---------------------------------------------------------------------------
inline Eigen::Vector3d camachoFrictionSlider(const Eigen::Vector3d& ds,
                                             double pj, double fCap,
                                             Eigen::Vector3d& slip)
{
    if (!(fCap > 0.0) || !(pj > 0.0)) {
        slip = ds;                          // aucun ressort de collage
        return Eigen::Vector3d::Zero();
    }
    Eigen::Vector3d tauTr = pj * (ds - slip);
    const double tt = tauTr.norm();
    if (tt <= fCap) return tauTr;          // collage : elastique
    Eigen::Vector3d tauF = tauTr * (fCap / tt);
    slip += (tauTr - tauF) / pj;            // retour radial vectoriel
    return tauF;
}

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
        // ---- facetRate = tensor (§2.1 eq. 11 de la note 2026) -------------
        // TENSEUR taux de deformation de l element, dans le repere GLOBAL —
        // c est celui dans lequel la normale de facette est exprimee, donc
        // le seul ou n . eps_point . n ait un sens sans rotation
        // supplementaire. Le champ `edot` ci-dessus n en garde qu un
        // SCALAIRE (la principale maximale en valeur absolue), si bien que
        // eps_point_eq (mode I) et gamma_point (mode II) ne sont pas
        // formables et que DIF_t et DIF_s partagent le meme argument — le
        // defaut de l audit du 2026-09-11, brique 2.1.
        // Reste a ZERO, et n est meme pas calcule, sans `facetRate = tensor`.
        Eigen::Matrix3d edotG = Eigen::Matrix3d::Zero();
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
                                           // jointTSL = camacho : glissement
                                           // plastique du CURSEUR DE COULOMB
                                           // (camachoFrictionSlider), 0 a
                                           // l insertion ; la part cohesive de
                                           // l eq. 19 ne le lit pas.
        bool dead = false;
        // Force NORMALE nette que le joint transmettait a l instant EXACT de
        // sa mort, en N (negatif = compression). Sortie de mesure seulement :
        // c est la charge que le relais au contact doit reprendre. Voir
        // jointDeath dans le header ci-dessous.
        double fDeath = 0.0;
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
        // En schema INTRINSEQUE (difIntrinsic_), le meme gel a lieu au
        // franchissement de l enveloppe et non a l insertion : ce drapeau
        // garantit qu il n a lieu QU UNE FOIS. Inerte en adaptatif.
        bool difStamped = false;
        double tanPhi = 0.0;                // friction
        double stat = 1.0;                  // Weibull strength factor (output)
        // ---- adaptive insertion (Yan et al. 2023, ported from FdemSolver) --
        // bonded = the joint does not exist yet: the co-located node copies
        // are rigidly bound (exact shared-node FEM) and processJoint skips
        // it. dn0 is the opening offset stamped at activation for stress
        // continuity (the newborn joint transmits at zero geometric opening
        // exactly the traction the bonded face was carrying).
        bool bonded = false;
        // groupContinuum.<corps> = true (2026-09-11) : facette interieure a
        // un corps declare CONTINU (acier du train de frappe, carbure de
        // l insert). Liee pour toujours : jamais evaluee par processJoint,
        // jamais inseree par le balayage adaptatif, absente du budget CFL.
        // Le corps est alors un continuum EF exact (noeuds partages via la
        // liaison par groupes), sans la complaisance E/(1 + 1/pf) ni le
        // pas de temps des ressorts de penalite.
        bool perm = false;
        double dn0 = 0.0;
        // largest opening ever reached per integration point (jointSoftening
        // = yan: the omax of eq. 17 — origin-secant unloading), as in 2D
        std::array<double, 3> omax{{0.0, 0.0, 0.0}};
        // largest SLIDING ever reached, the s_max of eq. 18, mesure depuis
        // l'origine figee slip[k] (jointShearUnload = origin seulement)
        std::array<double, 3> smax{{0.0, 0.0, 0.0}};
        // ---- jointSecantRatchet = on (2026-09-12, conseil M1/M7) ----------
        // Raideurs SECANTES non croissantes, par point d integration :
        // knr = sigma/dn du mode I (eq. 17), ksr = tau/s de l eq. 18. Une
        // secante ne peut que baisser : ni le DIF (ft(t) sous continuous) ni
        // la pression courante (tau_lim(sigma_n)) ne peuvent la remonter a
        // ouverture ou glissement fixe — c est la condition Phi >= 0 d une
        // loi d endommagement. -1 = jamais posee. Inertes si la cle est off.
        std::array<double, 3> knr{{-1.0, -1.0, -1.0}};
        std::array<double, 3> ksr{{-1.0, -1.0, -1.0}};
        // ---- failure mode, recorded ONCE when D first reaches 1 -----------
        // Ported from the 2D solver: partition of the eq. 16 damage driver at
        // the breaking instant. bmode = 1 tensile, 2 shear, 0 intact;
        // failMode = rn^2/(rn^2+rs^2) in [0,1], -1 while intact. Output only.
        int bmode = 0;
        double rnB = 0.0, rsB = 0.0;
        double failMode = -1.0;
        // ---- jointFailRule = majority (Solidity Y3Dfd.c, Sigma_tau) -------
        // Solidity ne porte PAS un endommagement par joint : chaque point
        // d integration a son propre z, sa propre traction, et le joint ne
        // meurt qu une fois DEUX points sur trois arrives a z >= 1
        //   (`if((nfail>1)&&(i1elpr[ielem]>=0))`, ligne 1175).
        // Dk[k] est ce z par point. Inerte sous `any` (defaut), ou seul le
        // scalaire partage J.D pilote la loi — comportement historique.
        std::array<double, 3> Dk{{0.0, 0.0, 0.0}};
        // ---- strainRateDIFArm = continuous (Solidity, meme fichier) -------
        // Leur DIF est un multiplicateur RECALCULE a chaque pas dans la
        // boucle element (`dpeftdif` lignes 1448-1456), jamais fige. On ne
        // peut donc pas l appliquer en place sur ft/coh/Gf/GfII sans le
        // composer a l infini : on garde ici les valeurs de BASE, intactes,
        // et la loi reconstruit ft = ftB*DIF a chaque pas.
        // Renseignes une seule fois par assignJointProps/applyJointStatistics
        // (via snapBase()), inertes sous `insertion` et `envelope`.
        double ftB = 0.0, cohB = 0.0, GfB = 0.0, GfIIB = 0.0;
        bool baseSnapped = false;
        // ---- insertionHoldSteps (§2.2 de la note 2026) --------------------
        // Nombre de pas CONSECUTIFS pendant lesquels Phi_F >= 1. Remis a
        // ZERO des que Phi < 1 : c est ce qui distingue une sollicitation
        // reelle d une oscillation numerique, et ce qui evite les insertions
        // en rafale que la note redoute. Toujours tenu a jour (il ne coute
        // qu un entier par facette et par pas) mais sans aucun effet tant
        // que insertionHoldSteps vaut 1, son defaut.
        int nPhi = 0;
        // ---- jointTSL = camacho (§2.4 eq. 16-19) --------------------------
        // Etat FIGE a l instant de l insertion : traction effective
        // REELLEMENT transmise t_m^ins, beta = t_s0/t_n0, G_C (melange de
        // Benzeggagh-Kenane si jointMixLaw = bk) et delta_m^f = 2 G_C/t_m^ins.
        // Reste a son etat par defaut (ok() == false) sous `penalty`, ou la
        // loi de penalite historique tourne inchangee.
        jtsl::Stamp tsl;
        // Plus grande separation EFFECTIVE atteinte, PAR POINT d integration
        // — meme granularite que omax/smax de la loi de penalite. Pilote la
        // decharge secante et D = delta_m^max / delta_m^f (eq. 17).
        // Sous jointTSLRise > 0 c est une separation DECALEE (jtsl::effOffset)
        // : elle vaut dm0 = rise * delta_m^f a l insertion, jamais 0, et
        // l ouverture GEOMETRIQUE maximale s en deduit par dmMax - tsl.dm0.
        // Sous `penalty` ce tableau sert de MESURE seulement : on y range le
        // maximum de l ouverture normale, pour le recensement des facettes
        // inserees mais jamais ouvertes (§3.2 de la note).
        std::array<double, 3> dmMax{{0.0, 0.0, 0.0}};
        // ---- jointTSLRise (§2.4, branche ascendante courte) ----------------
        // DIRECTION TANGENTIELLE UNITAIRE de la traction tamponnee a
        // l insertion, dans le plan de la facette (tauV / |tauV|). C est le
        // long de cette direction que le decalage dm0 es / beta est pose sur
        // le glissement, pour que l eq. 19 rende exactement (t_n^ins, t_s^ins)
        // a separation geometrique NULLE — sans lui split() rendait ZERO a
        // l insertion, un Dirac de -t_ins sur chacune des ~30 000 insertions
        // du run du 2026-09-11. Reprojetee dans le plan courant a chaque pas,
        // comme l origine slip[k] de la loi de penalite. Nulle si la facette
        // s est inseree sans cisaillement (es = 0 : le decalage est purement
        // normal). Inerte sous `penalty` et a rise = 0.
        Eigen::Vector3d tsDir = Eigen::Vector3d::Zero();
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
        // ---- toolShape = pdc (2026-09-03, scenario shear seulement) ------
        // Cutter PDC : DISQUE fini chanfreine, geometrie dans ToolPdc3d.hpp.
        // `x` est l'ARETE DE COUPE (convention 2D), `radius` = cutterDia/2.
        // n, u, w = repere de la face (Tool.hpp:67-74 releve y -> z).
        // pdc = false (defaut) : sphere / plat, bit-identique.
        bool pdc = false;
        double thick = 0.0, rakeDeg = -20.0, cham = 0.0, chamDeg = 45.0;
        bool floorFlat = false;
        Eigen::Vector3d n{1, 0, 0}, u{0, 0, 1}, w{0, 1, 0};
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
    void checkFinite();                    // C4 (w20) : NaN/Inf reel
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
    // `n` = normale de la facette a l instant de l activation. Ajoutee pour
    // §2.1 eq. 11 : sous facetRate = tensor le DIF a besoin de la normale
    // pour projeter le tenseur taux (eps_point_eq et gamma_point), et sous
    // §2.4 la meme normale sert au tampon d insertion.
    void activateJoint(int jI, double sig, const Eigen::Vector3d& tauV,
                       double fsNow, const Eigen::Vector3d& n);
    // ---- DIF de Yang : application des facteurs a UN joint ---------------
    // Extrait de activateJoint() le 2026-08-25 pour etre partage avec le gel
    // a l AMORCAGE du schema intrinseque (difIntrinsic_). Arithmetique
    // inchangee : les reperes dif_yang_* de la suite verrouillent la
    // bit-identite du chemin adaptatif.
    // erT / erS : les DEUX taux de l eq. 13 — eps_point_eq (mode I) et
    // gamma_point (mode II). Sous facetRate = scalar (defaut) l appelant
    // passe DEUX FOIS le meme scalaire filtre, ce qui restitue exactement le
    // comportement de rockim_g0 (DIF_t et DIF_s partagent leur argument).
    void stampDif(Joint& J, double erT, double erS);
    // dnE / dnF / slipF depuis les resistances COURANTES du joint — un seul
    // site, partage par assignJointProps, applyJointStatistics, stampDif et
    // refreshDif (voir jointDeltaC).
    void setJointLengths(Joint& J);
    void snapBase(Joint& J);              // proprietes de base, une fois
    void refreshDif(Joint& J, double erT, double erS);  // DIFArm = continuous
    // ---- §2.1 eq. 10-11 : les DEUX moyennes de facette, en UN seul site --
    // facetWA rend le poids du tetra A ; 0,5 sous `facetAverage = arith`
    // (defaut), V_A/(V_A + V_B) sous `volume`. Toutes les moyennes de
    // facette du solveur passent par ces trois fonctions : c est ce qui
    // garantit que la contrainte et le taux voient la MEME ponderation —
    // l audit en avait denombre six sites independants.
    double facetWA(const Joint& J) const;
    Eigen::Matrix3d facetStress(const Joint& J) const;
    double facetEdot(const Joint& J) const;
    // §2.1 eq. 11 : taux a l interface, sur le MEME patch pondere.
    // epsEq = n . eps_point_F . n (mode I), gam = 2 ||eps_point_F . n -
    // (n . eps_point_F . n) n|| (mode II). Sous `facetRate = scalar`
    // (defaut) les deux recoivent le MEME scalaire filtre facetEdot(J) :
    // c est litteralement le comportement de rockim_g0, rendu visible.
    void facetRates(const Joint& J, const Eigen::Vector3d& n,
                    double& epsEq, double& gam) const;
    // §3.2 eq. 26 : les trois postes de dissipation de la MATIERE, integres
    // sur V0 — D_vp (wPlas), D_omega_c en traction (wDamT) et en compression
    // (wDamC), en J. Ecrits en FIN de ligne de history.csv sous
    // energyBreakdown = on, jamais sans.
    void energyBreakdownRow(std::ostream& os) const;
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
    // ---- jointResidualMu : le frottement RESIDUEL du joint rompu ---------
    // Coefficient de frottement vers lequel le joint GLISSE quand il
    // s endommage, par la meme f(D) que la cohesion : mu_eff = muRes +
    // (tan(frictionDeg) - muRes) f(D). A D = 0 c est le frottement de PIC, a
    // D = 1 c est muRes. C est la distinction de Y-Geo (AbuAisha et al. 2015,
    // eq. 7.5 : angle de frottement de FRACTURE phi_f distinct de l angle
    // interne), et l equivalent de ce que Solidity obtient en remettant le
    // joint rompu au contact et a son glissement (0,6 calcaire / 0,18
    // granite). rockim gardait jusqu ici le frottement de PIC a vie.
    // < 0 = non pose = comportement historique, bit-identique.
    // Generalise jointFrictionScaled : muRes = tan(frictionDeg) reproduit le
    // defaut, muRes = 0 reproduit jointFrictionScaled = 1. Les deux cles ne
    // peuvent donc pas etre posees ensemble (erreur de config).
    double muRes_ = -1.0;
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
    // ⚠️ CONSEIL DU 12/09 (M1) : sous `origin`, tau = tau_lim(sigma_n)/s_max
    // · s des que le cap est actif — la raideur secante SUIT la compression
    // courante. dtau/ddn != dsigma/ds : aucun potentiel, et un cycle a
    // glissement fixe (comprimer, glisser en retour, relacher) rend
    // ½(k2 - k1) s² > 0. Mesure : 2D V0 vs V19, 3D reference vs B1 (16 J
    // crees en 81 us). `plastic` est la forme conforme (retour radial,
    // Phi = tau·ds_p >= 0). `origin` reste disponible (principe VIII) avec
    // un AVERTISSEMENT ; jointSecantRatchet = on le rend dissipatif.
    bool shearOrigin_ = false;

    // jointSecantRatchet = on | off (defaut off, bit-identique) — conseil du
    // 12/09 (M1, M7). Les secantes de decharge des eq. 17 (mode I, toutes
    // branches) et 18 (mode II, branche origin) deviennent NON CROISSANTES
    // dans le temps (Joint::knr, Joint::ksr). C est la condition de
    // dissipation d une loi d endommagement : ni la remontee du DIF (ft(t)
    // sous strainRateDIFArm = continuous, qui reecrit dnE a chaque pas) ni
    // la remontee de tau_lim avec la compression ne peuvent restituer plus
    // d energie que la charge n en a stockee. Le DIF continue d elever
    // l ENVELOPPE atteignable en charge (le critere), pas la raideur.
    bool secRatchet_ = false;

    // jointNormalProxy = penalty | law (defaut penalty = bit-identique) — audit A
    // #1 du 13/09. La contrainte normale « geometrique » qui fixe s_E, la plage
    // coulomb et l amorcage du DIF etait pj·dn, alors que sous jointElastic =
    // parabolic la loi transmet 2·pj·dn en compression : sigma_n divise par
    // deux dans tout ce qui fixe une resistance de mode II (le cap, lui, lit
    // la vraie traction). `law` pose pjN_ = 2 sous parabolic (1 sinon).
    double pjN_ = 1.0;

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
    // insertion = none : continuum pur — miroir EXACT du 2D, voir
    // FdemSolver.hpp pour la mesure qui a motive la cle (porte de
    // insertion-pointe, commit 2ead636 du 2026-08-25).
    bool noJoints_ = false;
    // ---- Insertion preferentielle en POINTE — miroir EXACT du 2D --------
    // (principe III : memes cles, meme loi). Voir FdemSolver.hpp pour la
    // mesure qui la motive. En 3D la pointe est un FRONT, mais le test reste
    // porte par les SOMMETS de la facette : une facette dont un sommet porte
    // deja un joint rompu est en propagation. Defaut 1,0 = bit-identique.
    double tipFactor_ = 1.0;      // insertionTipFactor
    double tipD_ = 0.5;           // insertionTipDamage
    std::vector<char> vertTip_;
    long nNuc_ = 0, nProp_ = 0;
    // =====================================================================
    //  NOTE DE TRAVAIL DE SEPTEMBRE 2026 — « FDEM hybride a insertion
    //  adaptative », sections 2.1 a 2.6. Le noyau physique partage est
    //  include/rockim/JointTsl.hpp ; ici ne vivent que les INTERRUPTEURS et
    //  les compteurs. TOUS les defauts ci-dessous reproduisent rockim_g0
    //  bit pour bit : chaque branche nouvelle est un chemin mort explicite
    //  tant que sa cle n est pas posee (principe VIII, croissance par
    //  addition).
    // =====================================================================
    // ---- facetAverage = arith (defaut) | volume — §2.1 eq. 10 -----------
    // La note pondere la contrainte de facette par les VOLUMES des deux
    // tetras, sigma_F = (V+ sigma+ + V- sigma-)/(V+ + V-). rockim_g0 fait
    // une moyenne ARITHMETIQUE 0,5/0,5 alors que les volumes sont deja la
    // (Elem::V0) : sur un maillage gradue 17x (percussion), un gros tetra
    // et un petit pesent autant, et c est le petit — celui qui voit la
    // singularite — qui est dilue de moitie.
    bool facetVolAvg_ = false;
    // facetAverage = max (2026-09-11 nuit) : le critere d insertion evalue la
    // traction de facette sur CHACUN des deux tetraedres et retient le plus
    // charge, au lieu de leur moyenne. Motif : sous l insert la moyenne
    // dilue de moitie l element qui porte l anneau de traction hertzien ;
    // premier joint a 90 us en adaptatif contre 64 us en intrinseque (banc
    // s = 2,5), dix fois moins de facettes inserees. Opt-in ; la contrainte
    // de facette servie ailleurs (facetStress) reste la moyenne.
    bool facetMaxIns_ = false;
    // facetAverage = nodal (2026-09-12) : le critere d insertion lit la
    // traction REELLEMENT transmise par la facette liee, reconstruite par
    // partition des forces nodales internes (Camacho & Ortiz 1996, Pandolfi &
    // Ortiz 2002) : a chaque sommet, somme des forces internes des copies
    // situees du cote A du plan de la facette, attribuee a la facette au
    // prorata de son aire tributaire, puis t = -F / A_f. Ni moyenne de
    // deux elements, ni maximum : la force que le cote A pousse a travers
    // l interface. Opt-in ; elCenNow_ = centroides courants des elements,
    // remplis au debut de chaque balayage.
    bool facetNodal_ = false;
    std::vector<Eigen::Vector3d> elCenNow_;
    // Acceleration nodale du pas PRECEDENT (copie par copie, ecrite dans
    // integrate() sous facetNodal_) : la partition dynamique retranche
    // l inertie, F_A->B = sum_A (f_i - m_i a_i), sans quoi la reaction
    // interne des elements a une charge de contact sur le sommet (l insert,
    // ~1 GPa) serait lue comme une traction transmise a l interface.
    std::vector<Eigen::Vector3d> accN_;
    bool facetTractionNodal(const Joint& J, const Eigen::Vector3d& n,
                            Eigen::Vector3d& t) const;
    // ---- facetRate = scalar (defaut) | tensor — §2.1 eq. 11 -------------
    bool facetTensor_ = false;
    // ---- insertionCriterion = or (defaut) | elliptic — §2.2 eq. 12 ------
    // Le critere historique est un OU LOGIQUE. La note veut une NORME
    // EFFECTIVE : une facette a 80 % de sa resistance en traction ET a 80 %
    // en cisaillement s insere (Phi = 1,28) la ou le OU la laisse intacte.
    bool ellipticIns_ = false;
    // n_h de la note (2-3 recommandes). 1 = comportement historique.
    int holdSteps_ = 1;
    // Mesures de l hysteresis : plus grand compteur atteint, et nombre de
    // remises a zero (= oscillations effectivement filtrees). Sans elles la
    // cle serait active et muette.
    long nHoldMax_ = 0, nHoldReset_ = 0;
    // ---- difExpS — §2.2 eq. 13 ------------------------------------------
    // a_s, exposant du DIF de CISAILLEMENT, separe de a_t (difExpT_). La
    // note impose a_s < a_t : le cisaillement est moins sensible a la
    // vitesse que la traction. Defaut 0,07 = l exposant LITTERAL de l eq. 2
    // de Yang et al. 2025, celui que porte rockim::difCompressionYang(edot)
    // — donc bit-identique tant que la cle est absente. L exposant est passe
    // a la surcharge rockim::difCompressionYang(edot, n) de YangDif.hpp,
    // l UNIQUE transcription de l article (la copie locale difShearExp du
    // premier lot a ete resorbee le 2026-09-11). Pourquoi ce defaut n est
    // PAS difExpT_ : il le serait sous strainRateDIF = yang-fig2, ou difExpT_
    // vaut 0,1707 — la bit-identite passe avant la lettre du contrat, et
    // l eq. 2 de Yang est de toute facon CONFIRMEE a 0,07 par leur fig. 2a.
    double difExpS_ = 0.07;
    // ---- jointTSL = penalty (defaut) | camacho — §2.4 eq. 16-19 ---------
    // `camacho` : loi INITIALEMENT RIGIDE. Pas de raideur K_n physique (la
    // branche pre-rupture est portee par la matrice) ; depuis le 2026-09-11
    // elle porte la courte branche ascendante jointTSLRise (ci-dessous), et
    // sa raideur de charge t_ins/dm0 ENTRE au budget CFL pour toutes les
    // facettes, liees comprises — le « gain sur dt » du premier lot n existe
    // pas, il etait achete a credit. t_m^ins est la traction
    // REELLEMENT TRANSMISE au pas d activation (et non le seuil nominal,
    // comme l ecretage `min(sig, ft)` de la loi historique), si bien que
    // int t_m d(delta_m) = G_C independamment du maillage ET du depassement.
    bool tslCamacho_ = false;
    // ---- jointTSLRise — §2.4, « discontinuite temporelle des lois
    // extrinseques » (Papoulia, Sam & Vavasis 2003) ------------------------
    // delta_m^0 / delta_m^f : longueur de la tres courte branche ASCENDANTE
    // que la note prescrit (« delta_m^0 = 1e-3 delta_m^f »), en fraction de
    // delta_m^f. Elle BORNE la raideur de charge/decharge a t_ins/dm0 =
    // t_ins^2/(2 rise G_C) et supprime le Dirac d insertion. Mesure du
    // 2026-09-11 (out_note2026, _dtsafe, _fricoff) : a rise = 0 la secante
    // t_ins/(D dmF) atteint p99,9 = 210 x pj et max = 549 x pj a la trame 9,
    // hors de tout budget CFL (dt/3,49 n a fait que retarder l explosion :
    // bloc de 40 mm, noeuds a 18 m a ~91 us ; eGc de -45 J a -1707 J pour
    // 3,2 J injectes). Defaut 1e-3 sous jointTSL = camacho (cle neuve, opt-in,
    // aucun deck ne depend de rise = 0) ; rise = 0 explicite accepte avec
    // AVERTISSEMENT (reglage de banc) ; refusee hors camacho (inerte).
    double rise_ = 0.0;
    // ---- jointMixLaw = none (defaut) | bk, jointBKEta — §2.4 eq. 18 -----
    // G_C = G_Ic + (G_IIc - G_Ic) m^eta, ratio de mode GELE a l insertion.
    // Benzeggagh-Kenane n existait nulle part dans le depot (zero occurrence
    // avant ce lot).
    bool bkOn_ = false;
    double bkEta_ = 2.0;
    // ---- jointFrictionMobilised = off (defaut) | damage — §2.4 ----------
    // `damage` : le terme frottant mu <-t_n> est multiplie par D. Sans lui
    // il vaut PLEINE VALEUR des D = 0 et le cisaillement est compte deux
    // fois, une fois par la viscoplasticite de la matrice et une fois par
    // le joint — ce que la note interdit (« un mecanisme par physique »).
    bool fricMob_ = false;
    // ---- jointEtaN / jointEtaS, jointViscousInCriterion — §2.5 eq. 20 ---
    // Option B : amortisseur APRES insertion, en Pa.s/m, ALTERNATIVE au DIF
    // (l exclusion est verifiee a la lecture des cles). Le terme TANGENTIEL
    // n existait pas du tout en g0. `jointViscousInCriterion = off` evalue
    // le cap de Coulomb sur la traction ELASTIQUE et non sur la traction
    // totale : sinon l amortisseur recree un effet de vitesse sur le seuil,
    // exactement ce que l option B veut eviter.
    double etaN_ = 0.0, etaS_ = 0.0;
    bool viscInCrit_ = true;
    // ---- gbCombine = mean (defaut) | min — §2.6 eq. 21 -------------------
    // t_n0^F = chi_t min(t_n0^p, t_n0^q) sur les facettes INTER-granulaires.
    // g0 en fait une MOYENNE : un grain fort protege son voisin faible,
    // alors que la note veut le maillon faible.
    bool gbMin_ = false;
    // ---- jointWeibullXu / jointWeibullScale — §2.6 eq. 22 ----------------
    // x_u : seuil du Weibull a TROIS parametres, en FRACTION du seuil de
    // facette. x_0 est recalcule pour que la moyenne d ensemble reste 1 :
    // sans cela toute calibration existante — qui lit la valeur calibree
    // comme la MOYENNE de la population — serait silencieusement decalee.
    double wbXu_ = 0.0;
    // `area` : le facteur d echelle porte sur l AIRE de la facette et non
    // sur le VOLUME des deux voisins. C est ce que demande la note (le
    // maillon faible d une facette echantillonne une SURFACE), et c est la
    // convention des essais de rupture fragile.
    bool wbArea_ = false;
    // ---- compBandLength = solver (defaut) | tetEdge — §1.3 eq. 6 --------
    bool lcCompTet_ = false;
    // ---- energyBreakdown = off (defaut) | on — §3.2 eq. 26 --------------
    bool eBreak_ = false;
    // Pas de temps que le MEME deck aurait sous jointTSL = penalty (budget
    // kPara pj sur toutes les facettes). Comparaison HONNETE avec dt_ sous
    // camacho, quel qu en soit le signe : depuis le 2026-09-11 le budget
    // camacho porte max(k0, kPara pj) sur TOUTES les facettes, liees comprises
    // (le « facteur 3,49 sur le cout du run » du premier lot etait achete a
    // credit — il excluait les facettes liees, et le run a explose).
    // Ecrit une fois par computeStableDt.
    double dtBonded_ = -1.0;
    // §3.2, controle propre a l insertion adaptative : nombre de facettes
    // INSEREES mais jamais ouvertes au-dela de 0,05 delta_m^f. Si cette part
    // n est pas marginale, n_h ou le seuil sont trop bas.
    long nDormant_ = 0;
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
    // ---- B8 : diagnostic de l ENDOMMAGEMENT SOUS-CRITIQUE ----------------
    // Le « diffuse ratcheting » de la penalite intrinseque : la part des
    // joints qui portent deja un D > 0,01 alors qu AUCUN n a rompu. C est la
    // signature de la penalite intrinseque — chaque joint cede un peu, la
    // structure perd de la raideur, et la resistance apparente baisse sans
    // qu une seule fissure soit apparue.
    //
    // POURQUOI IL ARRIVE ICI. Le diagnostic existait depuis 2026-08-25, mais
    // en 2D SEULEMENT et enferme dans `if (scen_ == BRAZILIAN)` — donc sur un
    // essai de CALIBRATION, et pas du tout sur la percussion, c est-a-dire pas
    // sur le cas compare a Imperial. Le contre-audit du 2026-08-30 l a releve :
    // on mesurait l injection d energie de l insertion intrinseque sans pouvoir
    // dire quelle fraction des joints ratchetait. C est la jauge qui manquait.
    //
    // QUAND LE RELEVE EST FIGE. Percussion : au PIC d effort d outil, l instant
    // ou la structure est la plus chargee — analogue exact du pic brésilien.
    // Traction / cisaillement (pas d outil) : au dernier pas ou nBroken_ == 0,
    // soit juste avant la premiere fissure, comme en 2D.
    // Les valeurs de FIN de run sont imprimees en plus : en percussion le pic
    // precede le gros de la fissuration, et l ecart entre les deux se lit.
    //
    // Cout mesure : voir chantier_imperial_2026-08-29/B08_*.md. Le releve ne
    // touche AUCUNE force ; il lit J.D et ecrit trois scalaires.
    //
    // ============ AVERTISSEMENT DE LECTURE, MESURE LE 2026-08-30 ============
    // CE CHIFFRE NE SE COMPARE PAS D UNE PENALITE A L AUTRE.
    // D est un deplacement NORMALISE : dnE = ft h / (pf E) retrecit quand pf
    // monte (Fdem3dSolver.cpp:1588 et :1813), si bien que la MEME ouverture
    // physique se lit comme PLUS d endommagement. Mesure sur
    // configs/fdem3d_percussion.cfg, T = 2e-5, insertion intrinseque :
    //
    //   pf =  20  ->  0,1571 % des joints > D=0,01 ; D moyen 8,099e-4
    //   pf =  60  ->  0,1929 %                     ; D moyen 1,187e-3
    //   pf = 180  ->  0,2100 %                     ; D moyen 1,334e-3
    //   (OMP_NUM_THREADS = 1, convention de la suite ; a 4 fils la tendance
    //    est la meme, +36 % au lieu de +34 %)
    //
    // La jauge MONTE de 34 % (et le D moyen de 65 %) pour un pf multiplie par
    // NEUF, alors que le pic d effort d outil reste PLAT a 1,2 % pres
    // (8729 / 8557 / 8629 N). J avais fait l hypothese inverse — « raidir la
    // penalite doit reduire le ratcheting » — et la mesure la REFUTE.
    //
    // CE QUI, LUI, SE COMPARE : deux runs a MEME penalite, schema d insertion
    // different. Meme deck, pf = 20 :
    //   intrinseque -> 0,1571 % , D moyen 8,099e-4 , pic 8729 N
    //   adaptative  -> 0,0486 % , D moyen 2,498e-4 , pic 9495 N
    // soit x3,24 de joints ratchetant ET un pic INFERIEUR de 8,1 %, sans
    // qu AUCUN joint n ait rompu (0/70000 des deux cotes). Les deux grandeurs
    // bougent ENSEMBLE d un schema a l autre, et SEPAREMENT d une penalite a
    // l autre : c est ce qui distingue l effet de schema de l artefact de
    // normalisation.
    //
    // ET LA JAUGE EST PLUS ROBUSTE QUE LE PIC. Le meme run a 1 et 4 fils donne
    // 8729 contre 8566 N de pic (1,9 % d ecart, REPRODUCTIBLE : peakF_ est un
    // max sur un signal oscillant, la moindre divergence de sommation change
    // quel maximum local est attrape), mais seulement 0,157143 contre
    // 0,155714 % de jauge (0,9 %) et 0,23 % sur le D moyen. Des deux
    // observables, c est la jauge qu il faut citer.
    //
    // CONSEQUENCE POUR B1 (poser jointPenaltyFactor ~ 9,6 au deck de replique) :
    // comparer gDfrac_ avant/apres ce changement N A AUCUN SENS. Le pic
    // d effort et le partage d energie, eux, se comparent.
    // =======================================================================
    double gDfrac_ = 0.0, gDmean_ = 0.0;   // au pic (ou avant 1re rupture)
    double gDpeakRef_ = -1.0;              // effort d outil du dernier releve
    bool   gDlatched_ = false;             // un releve a-t-il eu lieu ?
    long   gDscans_ = 0;                   // combien de balayages (cout)
    void scanSubCriticalDamage();
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
    // lawPhase = <nom> (2026-09-11 nuit) : la loi de volume ne s applique
    // qu aux elements de CETTE phase (la roche), les autres restent
    // elastiques (acier, carbure du train de frappe). -1 = toutes (historique,
    // qui exige une seule phase). Motif : l impact adaptatif laisse le
    // continuum sous l insert porter 1,4 GPa sans rien qui le fasse ceder ;
    // la loi de volume est la correction, et le deck a trois phases.
    int lawPhase_ = -1;
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
    // ---- DIF en schema INTRINSEQUE : le gel a l AMORCAGE (2026-08-25) -----
    // Vrai quand strainRateDIF est arme ET insertion = intrinsic. Le meme
    // balayage d enveloppe que l insertion adaptative est alors execute, mais
    // au franchissement il ne CREE pas le joint (il existe deja) : il se
    // contente d y STAMPER le DIF. Les deux schemas partagent donc le critere
    // exact, et ne different que par sa consequence — c est le temoin
    // recherche pour la comparaison a Yang et al. (joints intrinseques AVEC
    // DIF). Combinaison auparavant REFUSEE par une exception : aucune config
    // valide ne change de comportement (principe I).
    // ---- bulkModel : la loi de VOLUME (2026-08-25) ------------------------
    // `corotational` (defaut, historique) : decomposition polaire, deformation
    // de Biot, sigma = lambda tr(eps) I + 2 mu eps, assemblage P = R sigma.
    // Exact en grandes ROTATIONS, valable en petites DEFORMATIONS seulement.
    // `neohookean` : la loi de Guo (these Imperial 2014, eq. 2.6),
    //   T = (mu/J)(B - I) + (lambda/J) ln(J) I,
    // celle du code Solidity de Yang et al., avec l assemblage EXACT
    // P = J T F^-T. Elle redonne l elasticite lineaire au premier ordre et
    // n en diverge qu aux grandes deformations — sous l insert det F tombe a
    // 0,5-0,7. Le terme (lambda/J) ln J diverge quand J -> 0 : le materiau se
    // raidit sans borne a l ecrasement et l element ne peut plus s inverser,
    // ce que la loi lineaire ne fait pas (d ou le crushCap, garde-fou qui n
    // existe chez personne d autre).
    // ---- Loi de joint : les deux dernieres conventions de Guo -------------
    // jointElastic = linear (defaut) | parabolic : la branche elastique du
    // joint. Guo eq. 2.31 la pose PARABOLIQUE, sigma = ft(2r - r^2) avec
    // r = dn/dnE, ce qui donne une tangente NULLE au pic (transition douce
    // vers l adoucissement) et une pente initiale 2 pj, valable des DEUX
    // cotes de dn = 0 (loi C1 a l origine). rockim posait une droite de
    // pente pj, avec un coude au pic.
    // ---- jointQuadrature : les points d integration du joint -------------
    // `vertex` (defaut) : aux NOEUDS. C est la regle de Newton-Cotes nodale,
    // celle des elements cohesifs d Abaqus, retenue contre les oscillations
    // parasites du champ de traction a forte penalite (Schellekens & de
    // Borst). `midedge` : aux MILIEUX D ARETES, poids 1/3 — la regle de Guo
    // (Table 2.2), exacte a l ordre 2 au lieu de 1.
    // Les deux coincident EXACTEMENT en chargement uniforme ; elles different
    // de 50 a 200 % la ou l ouverture a un gradient a travers la facette,
    // c est-a-dire au front de fissure. Ce n est donc PAS un raffinement de
    // second ordre : c est un effet nul la ou on le mesure d habitude et fort
    // la ou la fissure se decide.
    bool midEdge_ = false;
    bool paraElastic_ = false;
    // jointDeltaC = exact (defaut) | guo : l ouverture critique. Guo
    // eq. 2.30 pose delta_c = 3 Gf/f mesure depuis ZERO, en approchant
    // l integrale de la z-curve par 1/3 la ou elle vaut 0,386307. Son
    // modele dissipe donc 3/0,386307 = 1,159 fois son Gf nominal. C est
    // SA convention, et ses Gf publies ont ete calibres avec : il faut la
    // reproduire pour retrouver ses chiffres.
    bool guoDeltaC_ = false;
    // jointDeltaC = solidity : ce que fait REELLEMENT leur code, qui n est
    // pas ce qu ecrit la these. Y3Dfd.c lignes 1098-1099 et 1125-1126 :
    //     op = 2 el ft / p0                   (l ouverture au pic)
    //     ot = MAXIM(2 op, 3 Gf/ft)           (la PLAGE d adoucissement)
    // la rupture est donc a op + ot, pas a 3 Gf/ft : la convention `guo`
    // oublie et l offset op et le plancher 2 op. Le plancher mord quand le
    // maillage est fin devant Gf/ft — exactement le regime d un impact.
    // Leur propre commentaire dit `/*need further investigation*/`, deux fois.
    bool solidityDeltaC_ = false;
    bool neoHooke_ = false;
    bool difIntrinsic_ = false;
    // strainRateDIFArm = continuous : le DIF recalcule a CHAQUE pas depuis
    // les proprietes de base (ftB, cohB, GfB, GfIIB), comme leur dpeftdif.
    // Exclusif de `insertion` et `envelope`, qui figent le facteur une fois.
    bool difContinuous_ = false;
    // jointFailRule = majority : le joint ne meurt qu une fois PLUS D UN
    // point d integration a D >= 1 (leur `nfail>1`) — 2 sur 3 en 3D, 2 sur 2
    // en 2D. Implique un endommagement PAR POINT (Joint::Dk), car la regle
    // n a pas de sens sur un scalaire partage qui est deja le max des points.
    bool majorityFail_ = false;
    // ---- gcBirth : la naissance d un contact sur un joint qui vient de mourir
    // `ramp` (defaut, historique) : la force part de ZERO et monte sur
    // gcBirthTau (on retranche un volume de reference qui decroit). Sur
    // l impact, ou les joints meurent EN COMPRESSION sous forte charge, c est
    // une perte de portance momentanee.
    // `penalty` : leur solution (Y3Did.c l. 915-964). Au pas exact de la
    // naissance ils lisent la force du joint mourant et RE-ECHELONNENT la
    // penalite de la paire, d1pepe = penalty*fn_joint/fn_contact, bornee a
    // [0,01 ; 3], de sorte que la force est CONTINUE. Le facteur persiste
    // ensuite pour la paire, et la raideur tangentielle le suit.
    bool birthPenalty_ = false;
    // gcBirth = relay (conseil du 12/09, D4/N7) : continuite de FORCE pour les
    // paires nees d un joint MORT (le calage de penalite de Solidity, comme
    // `penalty`), rampe de naissance (vRef, gcBirthTau) pour les paires SANS
    // joint mort (premier contact de deux corps) — `penalty` y injectait
    // ½ k delta0² (1 J a la naissance piston/bit, abort 2,9 us).
    bool birthRelay_ = false;
    // potStiffnessByPhase = max | min (defaut max = bit-identique) — audit B #10.
    // potP_ et potKt_ sont batis sur phases_.maxE() (600 GPa, carbure) pour
    // TOUTES les paires : le contact roche/roche est 10x trop raide. `min` :
    // la paire porte potPenaltyFactor * min(E_A, E_B) (et k_t au prorata).
    bool potByPhase_ = false;
    double potPF_ = 1.0;
    // bulkDamagePhase = <nom> (defaut : toutes) — audit D : l endommagement de
    // Yang n avait aucune garde de phase (acier et carbure endommageables).
    int bdPhase_ = -1;
    double birthPenMin_ = 0.01, birthPenMax_ = 3.0;   // leurs deux bornes
    long nBirthScaled_ = 0;                           // mesure : paires calees
    double birthScaleSum_ = 0.0;                      // ... et facteur moyen
    // ---- strainRateFilter : le taux qui alimente le DIF ------------------
    // `exponential` (defaut, historique) : passe-bas du premier ordre de
    // constante strainRateTau, parce que le taux brut par element est bruite.
    // `none` : le taux BRUT, ce que fait leur code — dpeftdif est calcule sur
    // le taux de l element tel quel, sans lissage (Y3Dfd.c l. 1448).
    bool srFilterOff_ = false;
    long nDifStamped_ = 0;                 // joints ayant recu le gel
    long nDifLate_ = 0;                    // dont DEJA endommages au gel
    // ---- jointDeath : QUAND le joint passe la main au contact -------------
    // `separation` (defaut, historique) : le joint ne meurt qu une fois
    // FRANCHEMENT ouvert (dnMax > 3 dnF). En compression il ne meurt donc
    // JAMAIS, et l algorithme de contact — qui porte le glissement contactMu,
    // le 0,6 de leur Table 4 — ne prend jamais le relais sous l insert.
    // `damage` : le joint meurt des que D >= 1, quel que soit le signe de
    // l ouverture. C est la regle de Guo (these Imperial 2014, §2.3.3) : « the
    // stress-displacement relation is not applied to this failed joint element
    // anymore ; instead, the interaction between the fracture walls will be
    // counted as contact forces that are calculated by the contact algorithm ».
    bool deathOnDamage_ = false;

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
        // ---- gcBirth = penalty (Solidity Y3Did.c l. 915-964) --------------
        // Facteur de penalite PROPRE A CETTE PAIRE, fige au pas de naissance
        // pour que la force du contact naissant egale celle du joint mourant.
        // C est leur d1pepe[icoup]. < 0 = pas encore ne. Inerte sous `ramp`.
        double penScale = -1.0;
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
    // groupContinuum.<corps> = true (2026-09-11) : corps sans joints, ses
    // facettes interieures (et celles de son interface avec un autre corps
    // continu lie par groupBond) sont Joint::perm. nPerm_ > 0 arme la
    // liaison de noeuds par groupes dans integrate(), comme l adaptatif.
    // Vide / 0 partout = aucun changement (bit-identique).
    std::vector<char> grpCont_;
    long nPerm_ = 0;
    bool anyGroupVel_ = false;             // un groupVel.<corps> est pose
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
    long nanEvery_ = 256;                  // nanCheckEvery (C4, w20)
    double work_ = 0.0, peakF_ = 0.0;
    long nBroken_ = 0;
    int nFrag_ = 1;
    double detachedVol_ = 0.0;
    Eigen::Vector3d gripF_{0, 0, 0};
    double sigmaPeak_ = 0.0;
};

} // namespace rockim
