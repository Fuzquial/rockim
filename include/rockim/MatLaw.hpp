#pragma once
// ---------------------------------------------------------------------------
// MatLaw: pluggable constitutive laws for the 3D continuum FEM (mode fem3d).
// One integration point per tet: the solver hands the law the CO-ROTATED
// total (Biot) strain and the element size, the law returns the Cauchy
// stress in the same frame and updates its internal state. Sign convention:
// tension positive, p = tr(sigma)/3 (tension positive).
//
// Laws (config key `law = ...`):
//
//  * elastic  — linear elasticity (lambda, mu from E, nu). The control case.
//
//  * dpr      — rate-INdependent elastoplastic damage:
//               Drucker-Prager cone F = sqrt(J2) + alpha I1 - k (alpha, k
//               matched to Mohr-Coulomb in TRIAXIAL COMPRESSION from
//               cohesion/frictionDeg), deviatoric radial return (psi = 0);
//               Rankine tensile damage: when the max principal EFFECTIVE
//               stress exceeds ft, a scalar damage D grows with exponential
//               softening regularized by the crack band (Gf smeared over the
//               element size lc), nominal stress = (1-D) sigma_eff.
//               The 3D sibling of the 2D fem module's DP/Rankine + damage.
//
//  * saksala2011 — FAITHFUL PORT of the thesis' vumat_saksala_2011.f90
//               (Saksala, IJNAMG 35 (2011) 1483-1505): DP cone with
//               non-associated potential (skBetaDP), modified Rankine
//               (associated), parabolic cap (associated) joined to the
//               cone, viscoplastic CONSISTENCY law with two viscosities
//               (skSdp compression, skSmr tension), confinement-dependent
//               cohesion degradation (skNd, confinement frozen at the
//               elastic trial state), logarithmic volumetric cap hardening
//               (skDcap, skWcap), exponential tensile damage driven by the
//               positive norm of the viscoplastic strain increment (skAt,
//               skBetaT), unilateral crack closure, Koiter corner rule.
//               Deliberately NO erosion, NO fracture-energy regularization,
//               NO compressive damage — as in the published law and the
//               VUMAT. Verified by trace superposition against the Fortran
//               reference on the three loading paths of its test harness.
//
//  * dpdfh    — FAITHFUL PORT of the thesis' VUMATS/dfh/vumat_kstdfh.f
//               (DP-DFH Bohus, Shariati/Saadati/Hild arXiv:2201.01870;
//               obscuration: Forquin & Hild 2010): Drucker-Prager
//               plasticity on the EFFECTIVE stress (f = q - p tan(beta) -
//               d, non-associated psi, radial return + apex), gated to
//               COMPRESSION only (OPTION-3: in tension the principal
//               stress rises elastically to the DFH initiation threshold);
//               ANISOTROPIC local DFH damage: three directional damages
//               D_i in a frame FROZEN at first initiation (Euler ZYX),
//               initiation thresholds = three sorted Weibull draws
//               sig_k = sigw (Zeff/V_el)^(1/m), V_el = lc^3, drawn
//               DETERMINISTICALLY from a 64-bit spatial hash of the
//               element's initial centroid (micrometre resolution,
//               xorshift64 — bit-compatible with the VUMAT's kst_seed);
//               growth by the cube-root obscuration integrator
//               D_i = 1 - exp(-x^3), dx = (S lam)^(1/3) k c dt with
//               lam = (sigma_i/sigw)^m / Zeff floored at 1/V_el and
//               c = sqrt(E/rho); nominal stress: tensile principal
//               components scaled by (1-D_i), shears by min(f_i, f_j)
//               (closed crack = full shear transfer). NO regularization,
//               deletion OFF by default (dfhDeld = 1e9), exactly like the
//               VUMAT (percussion practice: damaged elements stay as the
//               rubble bed). ftScale multiplies the Weibull scale (the
//               solver's FIELD mechanism — sandbox extension, neutral at
//               ftScale = 1). Verified by trace superposition against the
//               ifx-compiled Fortran reference on four loading paths.
//
//  * dfhplus  — DFH+ ETAPE 1 (2026-09-06) : MEME perimetre physique que
//               dpdfh, CADRE different. La loi est dans un fichier SEPARE,
//               src/MatLawDfhPlus.cpp (fabrique : include/rockim/
//               MatLawDfhPlus.hpp) — DpDfhLaw ci-dessus n'est pas touchee.
//               Elle est batie sur une ENERGIE LIBRE POSTULEE unique, avec
//               DECOMPOSITION SPECTRALE de la deformation elastique (seule la
//               partie POSITIVE est degradee, unilateralite naturelle) :
//                 rho psi = (lambda/2)[g_v <tr e>+^2 + <tr e>-^2]
//                         + G A:(e+)^2 + G ||e-||^2,
//               A = somme_i (1-D_i) n_i (x) n_i (repere FIGE), g_v = moyenne
//               HARMONIQUE des trois integrites. D'ou, ANALYTIQUEMENT,
//               sigma = d(rho psi)/d e (projection spectrale de Daleckii-
//               Krein) et Y_i = -d(rho psi)/dD_i >= 0 INCONDITIONNELLEMENT,
//               qui pilote l'obscuration par sigma^eq = sqrt(2 E Y_i / c_nu).
//               Plasticite DP non associee sur la contrainte effective
//               d(rho psi)/d e a D = 0 = C : e^e — meme critere, meme retour,
//               meme apex que dpdfh (compression et triaxiaux identiques a
//               2,5e-10). MEMES CLES dfh* : les cartes sont interchangeables.
//               Cles propres : dfhpPsiClamp (true = ecretage d'admissibilite
//               de la dilatance), dfhpVolInteg (harmonic | min | none).
//               Expose son energie libre au banc (hasFreeEnergy). Resultat :
//               `rockim thermobench dfhplus` PASS, 0 violation, la ou dpdfh en
//               a 7 172. Bancs : `rockim selftest-dfhplus`. Compte rendu :
//               docs/DFHPLUS_etape1.md.
//
//  * saksala  — rate-DEPENDENT damage-viscoplasticity in the spirit of
//               Saksala's model for percussive drilling: the SAME DP cone
//               but with a PERZYNA viscoplastic return (linear overstress,
//               viscosity saksalaEta [Pa s]: consistency F = eta dlam/dt,
//               closed form dlam = F_tr / (G + eta/dt)), a compression CAP
//               on the pressure (crush pressure capP0, linear hardening
//               capH: volumetric plastic return with pc += H |deps_v^p|)
//               and the same crack-band Rankine tensile damage.
//               SIMPLIFICATIONS stated plainly: Perzyna overstress instead
//               of Saksala's consistency formulation, scalar tensile damage
//               only (no separate compressive damage variable), no
//               viscosity on the damage side, zero dilatancy. The rate
//               effect on strength is exact and verified: steady uniaxial
//               overstress = sqrt(3) eta epdot / (1/sqrt(3) - alpha).
//
//  * cdp      — CONCRETE DAMAGED PLASTICITY (Lubliner 1989, Lee & Fenves
//               1998), la loi native d'Abaqus, portee le 2026-09-04 pour la
//               calibration triaxiale Red Bohus (CONTINUUM/calib_bohus_triax/
//               cdp_rockim). Surface de Lubliner en espace EFFECTIF
//               F = [q - 3 alpha p + beta <s_max> - gamma <-s_max>]/(1-alpha)
//               - sigma_bar_c (p = -tr/3 COMPRESSION POSITIVE), alpha de
//               fb0/fc0, gamma de Kc, beta de sigma_bar_c/sigma_bar_t ;
//               potentiel hyperbolique G = sqrt((ecc ft tan psi)^2 + q^2)
//               - p tan psi (non associe) ; retour spectral pleinement
//               implicite (crochet + regula falsi Illinois, jamais de Newton
//               pur : F n'est que C0) ; ecrouissages eps_t_pl / eps_c_pl de
//               Lee-Fenves (h = [r, 0, -(1-r)], clampes) ; endommagement
//               scalaire d = 1 - (1 - s_t d_c)(1 - s_c d_t) recalcule a
//               chaque appel (d depend de r, pas une variable d'etat),
//               sigma_nom = (1-d) sigma_bar. Tables Abaqus : compression en
//               deformation inelastique (abscisses fusionnees, converties
//               eps_c_pl = eps_in - d_c/(1-d_c) sigma_c/E0 a l'init) ;
//               traction en ouverture u_ck (GFI par defaut : sigma_t =
//               ft (1 - u/u_t0), u_t0 = 2 Gf/ft, ou cdpTension = table =
//               type DISPLACEMENT ; type STRAIN NON couvert), convertie A
//               LA VOLEE par element eps_t_pl = u/lc - d_t/(1-d_t) sigma_t/E0
//               (« lc partout » = Abaqus avec l0 = lc), interpolation
//               LINEAIRE EN eps_pl entre noeuds convertis (lecture B = doc
//               Abaqus ; l'aire tractive x lc vaut Gf + lc wDamT, +0,36 % a
//               1 mm) ; plancher sigma_t >= 1e-3 ft (nominal). Gardes a
//               l'init : G1 lc <= 2 E Gf/ft^2, G2 noeuds convertis
//               croissants, G3 snap-back EFFECTIF (pente locale de
//               sigma_bar_t, lc < 10,7 mm sur la carte Bohus). Etat opt-in
//               MatState::Cdp {epsTpl, epsCpl, d} ; champs partages
//               s.D = d_t (monotone), s.Dc = d_c, s.epvEq = eps_c_pl,
//               s.kappa = eps_t_pl. Erosion par les cles existantes
//               (erodeD sur d_t en traction nette — avertissement si
//               erodeD > d_t,max : canal inerte ; erodeWfrac inerte en CDP,
//               la dissipation tractive est plastique ; erodeEpv sur
//               eps_c_pl ; erodeDc sur d_c/d_c,max, garde compDamage levee).
//               POINT AVEUGLE de parite : l'adoucissement compressif est en
//               deformation, NON regularise en maillage — la « chute
//               verticale » Abaqus post-pic est reproduite a l'identique.
//               Cles : cdpDilationDeg (35), cdpEcc (0,1), cdpFbFc (1,16),
//               cdpKc (0,667), cdpWt (0), cdpWc (1), cdpHardening
//               (« sigma:eps_in ... »), cdpCompDamage (« d:eps_in ... »),
//               cdpTension gfi | table, cdpTensionTable (« sigma:u_ck »),
//               cdpTensionDamage (« d:u_ck ... ») ; defauts = carte
//               historique Red Bohus des decks. Revue du 2026-09-04 :
//               (1) ftScale (Weibull du solveur, matWeibullM) est HONORE
//               comme par dpr (E1) : ft, a = ecc ft tan psi et la table de
//               traction (sigma_t x s, u_ck x kappa avec kappa = 1/s si
//               weibullScope = strength — Gf conserve — ou 1 si strengthGf)
//               sont mis a l'echelle par element ; la table de COMPRESSION
//               n'est pas touchee (comme la cohesion de dpr) ; les gardes
//               G1-G3 sont re-verifiees au premier appel de tout element a
//               ftScale != 1 (G3 se durcit en s^2 : lc < 10,7 mm / s^2 sur la
//               carte Bohus). (2) cdpViscosity (0 = rate-independant, bit-
//               identique) : regularisation viscoplastique de Duvaut-Lions
//               de la CDP d'Abaqus/STANDARD, eps_pl,v += dt/(mu+dt) (eps_pl
//               - eps_pl,v), d_v idem, sigma = (1 - d_v) D0 (eps - eps_pl,v).
//               Les decks de reference (*Dynamic, Explicit) portent mu =
//               5e-5 s mais Abaqus/EXPLICIT IGNORE ce parametre (mot-cle
//               *CONCRETE DAMAGED PLASTICITY, 5e donnee : « ignored in
//               Abaqus/Explicit ») : la parite avec les runs de reference est
//               le defaut rate-independant, PAS cdpViscosity = 5e-5. Si on
//               l'active tout de meme (Standard), a 0,75/s la surcontrainte
//               triaxiale n'est PAS E0 mu epsdot ~ 3 MPa mais ~6 fois plus
//               (la laterale visqueuse doit etre compensee par un confinement
//               inviscide que la pente m amplifie) : chiffre dans §5.19.
//               (3) cles cdp* inconnues REFUSEES (faute de frappe = exception).
//               (4) PLANCHER sigma_t >= 1e-3 ft : au-dela de la fin de table
//               le potentiel hyperbolique coule dans les TROIS directions
//               (d eps_lat = lambda (-rho/2 + t/3) > 0 des que rho < 2t/3) :
//               un element fissure en traction GONFLE lateralement et sa
//               deformation volumique plastique croit sans borne — parite
//               Abaqus (meme potentiel, meme residu) ; le garde-fou est
//               erodeD <= d_t,max (retrait de l'element fissure).
//               (5) CRACK BAND EN COMPRESSION, cle cdpCompLength (m ; 0 =
//               absent = defaut = bit-identique, tables lues telles quelles)
//               — 2026-09-04 soir, decision de Fernando : la chute post-pic
//               d'un triaxial est STRUCTURELLE (bande localisee), la table CDP
//               en deformation n'est pas regularisee, l'energie dissipee dans
//               la bande dependait de la maille. Si L_ref > 0, les abscisses
//               eps_in des tables de compression (sigma_c et d_c) sont lues
//               comme definies a L_ref et SEULE LA BRANCHE POST-PIC est remise
//               a l'echelle de l'element : eps_in,loc = eps_in,pic + (eps_in -
//               eps_in,pic) L_ref/lc (pic = premier noeud ou sigma_c nominal
//               est maximal), branche pre-pic intacte, puis conversion
//               d'Abaqus eps_c_pl = eps_in,loc - d_c/(1-d_c) sigma_c/E0 A LA
//               VOLEE par element (compState(epl, lc)). Le deplacement
//               inelastique post-pic u_in = (eps_in - eps_in,pic) L_ref est
//               invariant, donc G_c = int (sigma_c - sigma_res) du_in aussi
//               (Bazant-Oh en compression, logique de compGIIc). Gardes a
//               l'init a lcMax, seulement si la cle est > 0 : (G2c) noeuds
//               convertis croissants, (G3) pente effective <= limC, (G4)
//               snap-back au point materiel -d sigma_bar_c/d eps_c_pl <= E0
//               localement (<=> d eps_total/d eps_pl >= 0 ; en corde
//               -Delta sigma_c/Delta eps_in,loc <= E0). NB : le critere
//               nominal « -d sigma_c/d eps_pl <= E0 » n'est PAS le snap-back
//               (80,5 GPa a lc 4 mm / L_ref 2 mm sur la carte historique,
//               corde 41,6 GPa) : G4 implemente le critere exact. Banc (m) :
//               aire adoucie post-pic x lc invariante a lc = 1 / 2 / 4 mm et
//               egale a G_c(L_ref), branche pre-pic identique a 1e-9, sans
//               la cle rapport 4 entre 1 et 4 mm.
//               (6) CAP VOLUMIQUE DE COMPACTION, cle cdpCap (false = absent
//               = bit-identique), cdpCapP0 (Pa, obligatoire si cdpCap),
//               cdpCapH (Pa, defaut K) — 2026-09-04 soir, constat de Fernando
//               sur la percussion quart de bloc : la CDP n'a pas de cap et son
//               meridien est une droite, sous l'insert (~8 GPa de contact) la
//               roche reste elastique, rien n'est broye. Meme brique que le
//               cap de dpr/saksala (PlasticDamageLaw::stress), placee ENTRE le
//               predicteur elastique et le retour de Lubliner, sur la pression
//               EFFECTIVE p_bar = -tr(sigma_bar_tr)/3 : si p_bar > pc (pc
//               initialisee a cdpCapP0 au premier appel, etat MatState::pc),
//               dev = (p_bar - pc)/(K + H), eps_pl -= dev/3 I (compaction),
//               pc += H dev, sigma_bar_tr += K dev I ; puis Lubliner sur le
//               predicteur corrige. Le cap ne touche NI eps_t_pl/eps_c_pl NI
//               d_t/d_c (comme dpr : ecrouissage propre, la compaction
//               n'alimente pas l'endommagement de la table) ; compte dans
//               wPlas. Cles orphelines (cdpCapP0 sans cdpCap) refusees.
//               REVUE (nuit) : le cap est impose au PREDICTEUR seulement ;
//               le retour de Lubliner qui suit est DILATANT (p = p_tr + K tan
//               psi lambda), donc quand cap et cone agissent dans le meme pas
//               p_bar de fin de pas = pc + K (d tr eps_pl + d pc/H) > pc
//               (depassement borne, recappe au predicteur suivant) — ce n'est
//               PAS le cap exact en fin de pas de dpr (retour deviatorique).
//               Le seuil NOMINAL est (1 - d) pc (0,08 pc0 a d_c = 0,92) :
//               cdpCapP0 se calibre en pression EFFECTIVE. La compaction
//               seule n'erode jamais (erodeEpv lit eps_c_pl, que le cap ne
//               nourrit pas) : en percussion garder erodeDetMin. Champs VTU
//               capPc / epsVpl via vtkCap = true (Fem3dSolver, opt-in).
//               Bancs (n) hydrostatique (p = K eps_v puis pente K H/(K+H),
//               cdpCapH absent = K) et (o) cap + cone : o1 20 MPa/cap 440
//               jamais actif (p_bar max 256 MPa) et q identique a 1e-9 ;
//               o2 20 MPa/cap 200 : identite avant activation puis DOIT
//               diverger ; o3 100 MPa/cap 100 : cap et cone ensemble, q_pic
//               invariant, decalage eps_v^pl/3, borne du depassement pas a pas.
//               Bancs :
//               `rockim selftest-cdp` (15 groupes de controles chiffres, code
//               retour 1 si un pilotage lateral ne converge pas) ; pilote
//               generique `rockim matpoint <cfg> [csv]` (toutes les lois ;
//               cles mp* ; chemins triax | tension | biaxial | uniaxial |
//               hydro (deformation isotrope imposee, sortie p_nom, eps_v, pc,
//               eps_v_pl) ; colonnes `residu` et `eps_in_c` ; code retour 2
//               si un pas n'a pas converge).
//
// The crack band imposes the classical limit lc <= E Gf / ft^2 (snapback):
// the law throws at init if the mesh is too coarse for (Gf, ft).
//
// ---- BRIQUES opt-in du noyau dpr/saksala (etude « briques constitutives »,
//      2026-09-03 ; defaut = bit-identique cles absentes) ---------------------
//  * meridian = linear (defaut) | power : en compression triaxiale le cone
//    lineaire donne q = UCS + (N_phi - 1) sigma3 ; `power` le remplace, pour
//    sigma3 >= 0, par la loi puissance de Saksala-Hokka / Hoek-Brown
//        q = fc0 + merB * sigma3^merN      (merB en MPa^(1-n), sigma3 en MPa)
//    avec fc0 = UCS du cone (sigmaCdp), lue dans le plan (p, q) via la
//    pseudo-sigma3 = -p - q/3 (exacte en triaxial de compression), Newton 1D
//    sur q a p^trial fixe, retour radial deviatorique inchange. Cote tractif
//    (pseudo-sigma3 < 0) et forme de Lode (cercle) IDENTIQUES au cone lineaire :
//    l'ablation isole la courbure. Continu en sigma3 = 0 (q = fc0 des deux cotes).
//  * compDamage = none (defaut) | crackband : endommagement COMPRESSIF scalaire
//    omega_c = compAc (1 - exp(-b_c epvEq)), b_c = fc0 h_e / compGIIc (h_e = lc),
//    applique a la partie spectrale NEGATIVE de la contrainte effective
//    (omega_t / D reste sur la partie positive) — la brique omega_c de la loi
//    MH de Saksala 2018 et le d_c de CDP. Etat expose : MatState::Dc.
//  * capP0 / capH sont desormais lus AUSSI par law = dpr (defaut capP0 = 0 =
//    sans cap, bit-identique) : l'ablation du cap change UNE cle sans changer
//    de loi ni de viscosite.
//  * erodeDc (0 = off) : troisieme canal d'erosion, sur omega_c / compAc
//    NORMALISE (un seuil 0,99 est atteignable ; un seuil brut >= compAc ne le
//    serait jamais — c'est le dmax du VUMAT MH).
//  * Compteurs de dissipation cumules (densites, J/m^3, toujours calcules,
//    sans effet sur la contrainte) : wPlas = int sigma_nom : d eps_p,
//    wDamT = int Y_t dD, wDamC = int Y_c d omega_c, Y = 1/2 sigma^+- : C^-1 sigma^+-.
//    Sur la barre E1 : sum(wDamT V0) = Gf * A (theoreme du banc 4).
//  * rockim selftest-triax [out.csv] : empreinte triaxiale des briques a la
//    carte Red Bohus (cibles analytiques + donnees 404,8 / 599 / 704 / 799,3).
//  * tensionDamage = scalar (defaut) | fixed (2026-09-05, w19) : endommagement
//    de traction A DIRECTION FIGEE (fixed crack model, Rashid 1968 ; Rots &
//    Blaauwendraad 1989) pour dpr/saksala. Au premier depassement du seuil
//    de Rankine (sigma1_eff > ft, ou eps1 > k0 en pilotage deformation) le
//    repere principal est FIGE (Euler ZYX, MatState::Fcm::eul, comme
//    Dfh::eul) ; jusqu'a trois directions s'amorcent successivement, chacune
//    dans le complement orthogonal des precedentes (rotation de (n2, n3)
//    autour de n1 vers la direction principale du complement). Chaque
//    direction i porte son d_i avec la MEME cinetique de bande de fissuration
//    que le scalaire (kappa_i = max de la composante normale n_i.T.n_i du
//    tenseur pilote, d_i = 1 - k0/kappa_i exp(-(kappa_i - k0)/kf), kf =
//    Gf/(lc ft) - k0/2 : ouverture u_i = kappa_i lc). Contrainte nominale
//    dans le repere fige : S_ii -> (1 - d_i) S_ii si S_ii > 0 (unilateral :
//    composante normale COMPRESSIVE intacte), S_ij -> (1 - beta max(d_i^ouv,
//    d_j^ouv)) S_ij avec beta = tensionShearRetention (defaut 1 = pas de
//    retention, la convention min(f_i, f_j) de la VUMAT DP-DFH ; beta < 1 =
//    retention de cisaillement de Rots, d^ouv = d si la fissure est ouverte,
//    0 sinon). Raideur perdue normalement au plan de fissure, conservee
//    parallelement. s.D = max_i d_i (history, VTU « damage », spall),
//    s.kappa = max_i kappa_i ; champs VTU dFix1..3 et colonnes probes
//    dFix1..3 avec la cle. wDamT += 1/2 <S_ii>^2/E dd_i (trois directions).
//    Plasticite DP, cap, omega_c (sur la partie spectrale negative de sigma_eff,
//    soustraite apres l'operateur fige), erosion : inchanges. Refuse pour
//    saksala2011 (port fidele, endommagement pilote par la deformation
//    viscoplastique, pas de bande (ft, Gf, lc)), cdp (ses propres tables) et
//    dpdfh (deja directionnel). rockim selftest-fixed [out.csv] : bancs (a)-(e).
// ---------------------------------------------------------------------------
#include <memory>
#include <string>

#include <Eigen/Dense>

#include "rockim/Config.hpp"
#include "rockim/Material.hpp"

namespace rockim {

struct MatState {
    Eigen::Matrix3d epsP = Eigen::Matrix3d::Zero();    // plastic strain
    double D = 0.0;                                    // tensile damage
    double kappa = 0.0;                                // damage driving strain
    double epvEq = 0.0;                                // equiv. viscopl. strain
    double pc = 0.0;                                   // current cap pressure
    bool eroded = false;

    // briques opt-in (2026-09-03) : endommagement compressif omega_c
    // (compDamage = crackband) et compteurs de dissipation cumules [J/m^3]
    double Dc = 0.0;
    double wPlas = 0.0, wDamT = 0.0, wDamC = 0.0;
    // canal qui a erode l'element (0 vivant, 1 spall D >= erodeD en traction
    // nette, 2 broyage epvEq >= erodeEpv, 3 omega_c / Ac >= erodeDc,
    // 4 dfhDeld) — diagnostic, sans effet sur la contrainte
    int eroCode = 0;

    // per-element strength factor set by the SOLVER before the first law
    // call (Weibull heterogeneity — the sandbox version of the VUMAT's
    // FIELD mechanism); 1 = homogeneous
    double ftScale = 1.0;

    // initial element centroid, set by the solver after meshing (dpdfh
    // seeds its Weibull draws from a spatial hash of these coordinates,
    // exactly like the VUMAT hashes coordMp)
    Eigen::Vector3d x0 = Eigen::Vector3d::Zero();

    // dpdfh sub-state (mirrors the VUMAT's 16-SDV layout; Voigt-6 order
    // 11 22 33 12 23 13, tensorial shears)
    struct Dfh {
        bool seeded = false;
        double snom[6] = {0, 0, 0, 0, 0, 0};   // last NOMINAL stress
        double epsPrev[6] = {0, 0, 0, 0, 0, 0};
        double Dv[3] = {0, 0, 0};              // SDV 4-6
        double eul[3] = {0, 0, 0};             // SDV 7-9 (Euler ZYX, rad)
        double sc[3] = {0, 0, 0};              // SDV 10-12 (sorted draws)
        double ti[3] = {0, 0, 0};              // SDV 13-15 (init times)
        double peeq = 0.0;                     // SDV 3
        double smaxh = 0.0;                    // SDV 16
        double t = 0.0;                        // accumulated total time
        bool dead = false;                     // SDV 1 (STATUS)
    } dfh;

    // saksala2011 sub-state (Voigt-6, order 11 22 33 12 23 13, tensorial
    // shears — mirrors the VUMAT's SDV layout, incl. the local strengths
    // SDV15/16)
    struct Sk11 {
        bool init = false;
        double sbar[6] = {0, 0, 0, 0, 0, 0};           // effective stress
        double epsPrev[6] = {0, 0, 0, 0, 0, 0};        // last total strain
        double kapDP = 0.0, kapMR = 0.0, eqvt = 0.0, epsv = 0.0, pp = 0.0;
        double ftLoc = 0.0, c0Loc = 0.0;               // local strengths
        int active = 0, failed = 0;
    } sk;

    // cdp sub-state (2026-09-04, opt-in : zero et sans effet pour les autres
    // lois) : variables d'ecrouissage de Lee-Fenves et d total (diagnostic ;
    // d_t, d_c se derivent des tables, ne sont pas stockes).
    // Revue du 2026-09-04 : `guarded` = gardes G1-G3 deja verifiees pour CET
    // element (elles dependent de lc et de ftScale ; a ftScale = 1 la
    // verification de l'init a lcMax suffit) ; `epsPv`, `dv` = deformation
    // plastique et endommagement VISQUEUX de Duvaut-Lions (cle cdpViscosity,
    // la regularisation viscoplastique de la CDP d'Abaqus/Standard, ignoree
    // par Abaqus/Explicit ; inertes a 0).
    struct Cdp {
        double epsTpl = 0.0, epsCpl = 0.0, d = 0.0;
        bool guarded = false;
        double dv = 0.0;
        Eigen::Matrix3d epsPv = Eigen::Matrix3d::Zero();
    } cdp;

    // endommagement de traction a direction figee (tensionDamage = fixed,
    // 2026-09-05, opt-in : zero et sans effet pour le chemin scalaire) :
    // nAct directions amorcees (0..3), d[i] / kap[i] = endommagement et
    // deformation pilote maximale de la direction i, eul = repere fige
    // (Euler ZYX, colonnes n1 n2 n3 de R = reul(eul), comme Dfh::eul)
    struct Fcm {
        int nAct = 0;
        double d[3] = {0, 0, 0};
        double kap[3] = {0, 0, 0};
        double eul[3] = {0, 0, 0};
    } fcm;

    // ---- law = dfhplus (2026-09-06, etape 1 du chantier DFH+) -------------
    // Sous-etat OPT-IN, nul et sans effet pour toutes les autres lois (croissance
    // par addition, principe VIII). La loi dfhplus est ecrite sur une ENERGIE
    // LIBRE POSTULEE (cf. src/MatLawDfhPlus.cpp) : sa seule variable d'etat
    // reversible est la deformation elastique, portee par le champ PARTAGE
    // MatState::epsP (eps^e = eps - eps^p) ; ici ne restent que les variables
    // IRREVERSIBLES et le repere fige.
    struct Dfhp {
        bool seeded = false;      // tirages de Weibull faits
        bool frozen = false;      // repere d'endommagement fige
        bool dead = false;        // erosion dfhDeld
        double Dv[3] = {0, 0, 0};        // endommagements directionnels D_i
        double eul[3] = {0, 0, 0};       // repere fige (Euler ZYX, comme Dfh)
        double sc[3] = {0, 0, 0};        // seuils de Weibull tries
        double ti[3] = {0, 0, 0};        // instants d'amorcage
        double peeq = 0.0;               // multiplicateur plastique cumule
        double smaxh = 0.0;              // max de la contrainte equivalente [Pa]
        double t = 0.0;                  // temps total cumule
        double psi = 0.0;                // rho psi au dernier appel [J/m^3]
        double wDam = 0.0;               // int sum Y_i dD_i [J/m^3]
        double clamp = 0.0;              // tan(Psi) effectif du dernier retour
        int nClamp = 0;                  // nb de pas ou l'ecretage a agi
    } dfhp;
};

class MatLaw {
public:
    virtual ~MatLaw() = default;
    // eps: co-rotated Biot strain; lc: element size (V^(1/3)) for the crack
    // band; dt: time step (rate laws). Returns nominal Cauchy stress.
    virtual Eigen::Matrix3d stress(const Eigen::Matrix3d& eps, MatState& s,
                                   double dt, double lc) const = 0;
    virtual std::string name() const = 0;
    double cP() const { return mat_.cP(); }

    // ---- ENERGIE LIBRE EXPOSEE (2026-09-06, AJOUT du chantier DFH+) -------
    // Une loi construite sur une energie libre POSTULEE peut la rendre : le
    // banc thermodynamique (`rockim thermobench`) cesse alors d'estimer rho psi
    // par une sonde de decharge (estimateur biaise, increments contamines et
    // non relaches exclus) et lit la valeur EXACTE. Defaut : « non exposee »,
    // donc STRICTEMENT sans effet sur les lois existantes (dpdfh comprise) —
    // croissance par addition, principe VIII.
    virtual bool hasFreeEnergy() const { return false; }
    // rho psi(eps, etat) [J/m^3] a variables internes FIGEES ; n'a de sens que
    // si hasFreeEnergy() est vrai. `eps` est la deformation TOTALE, comme pour
    // stress() ; l'etat porte la part plastique et l'endommagement.
    virtual double freeEnergy(const Eigen::Matrix3d&, const MatState&) const {
        return 0.0;
    }

    // Forces motrices ANALYTIQUES Y_k = -d(rho psi)/d D_k [J/m^3], remplies
    // dans `Y` (au plus 3) ; rend le nombre ecrit, 0 = « non exposees ».
    // Meme logique additive que freeEnergy : le banc peut alors confronter la
    // formule analytique aux differences finies sur l'energie libre exposee.
    virtual int damageForces(const Eigen::Matrix3d&, const MatState&,
                             double*) const {
        return 0;
    }

    // uniaxial compressive yield of the DP cone (analytic, for verification)
    double sigmaCdp() const { return kdp_ / (1.0 / std::sqrt(3.0) - adp_); }
    // steady uniaxial viscous overstress at axial strain rate epdot
    // (0 for rate-independent laws)
    virtual double viscousOverstress(double) const { return 0.0; }

    static std::unique_ptr<MatLaw> make(const std::string& kind,
                                        const Material& m, const Config& c,
                                        double lcMax);

protected:
    explicit MatLaw(const Material& m) : mat_(m) {
        lam_ = m.E * m.nu / ((1.0 + m.nu) * (1.0 - 2.0 * m.nu));
        G_ = m.G();
        K_ = m.K();
        double sf = std::sin(m.phiDeg * M_PI / 180.0);
        double cf = std::cos(m.phiDeg * M_PI / 180.0);
        adp_ = 2.0 * sf / (std::sqrt(3.0) * (3.0 - sf));
        kdp_ = 6.0 * mat_.cohesion * cf / (std::sqrt(3.0) * (3.0 - sf));
    }
    Eigen::Matrix3d elastic(const Eigen::Matrix3d& eps) const {
        return lam_ * eps.trace() * Eigen::Matrix3d::Identity()
               + 2.0 * G_ * eps;
    }
    Material mat_;
    double lam_ = 0.0, G_ = 0.0, K_ = 0.0;
    double adp_ = 0.0, kdp_ = 0.0;                     // DP cone (MC triax)

    // ---- PORTEE de l'heterogeneite de Weibull (cle weibullScope) ----------
    // Le tirage par element (ftScale) multiplie la RESISTANCE. Reste une
    // question que le depot tranchait jusqu'ici differemment a trois endroits
    // sans que ce soit un choix : la TENACITE suit-elle ?
    //
    //   strength   (defaut) — seule ft est mise a l'echelle. C'est la
    //     convention d'une population de DEFAUTS : le materiau a partout la
    //     meme energie de rupture, seuls les seuils d'amorcage sont
    //     distribues. Consequence : la longueur cohesive locale
    //     E Gf / ft^2 varie en 1/s^2 et l'ouverture critique en 1/s — un
    //     element faible a une zone cohesive PLUS GRANDE. C'est aussi ce que
    //     fait deja la statistique de JOINT (J.ft et J.coh seuls).
    //
    //   strengthGf — ft ET Gf suivent le meme facteur. C'est la convention
    //     du materiau AUTO-SEMBLABLE : l'ouverture critique dnF ~ Gf/ft est
    //     invariante, la longueur cohesive varie en 1/s. C'est exactement ce
    //     que fait deja le DIF de Yang (leurs eq. 2-3 amplifient ft et Gf du
    //     meme facteur, precisement pour que dnF ne bouge pas et que le
    //     compteur d'endommagement ne soit pas corrompu).
    //
    // Aucune des deux n'est « la bonne » : elles decrivent deux materiaux
    // differents. Le defaut reproduit le comportement anterieur du depot.
    bool wScaleGf_ = false;
};

// Replays the three material-point loading paths of the Fortran reference
// harness (VUMATS/saksala/test_saksala_2011.f90) through the C++ port and
// writes the same CSV trace, for quantitative superposition against the
// ifx-compiled reference. Returns 0.
int saksala2011Selftest(const std::string& csvPath);

// Autotest point-materiel de la loi mc (Mohr-Coulomb de Ye et al. 2025) :
// compare les plateaux plastiques aux formules exactes de Mohr-Coulomb en
// traction, compression uniaxiale et compression triaxiale. Retourne 0 si
// l'ecart max est < 1 %.
int mcSelftest(const std::string& csvPath);

// Same idea for the DP-DFH port: replays the four material-point paths of
// VUMATS/dfh/test_kstdfh.f90 (two tension rates, deviatoric compression,
// oblique tension with shears) and writes the same CSV (stresses in MPa).
int dpdfhSelftest(const std::string& csvPath);

// Empreinte TRIAXIALE des briques du noyau (2026-09-03) : compression
// triaxiale pilotee en deformation a sigma3 = 0/20/50/75/100/300/500 MPa,
// carte Red Bohus (E 77,66 GPa, nu 0,29, c 22,77 MPa, phi 50,4 deg = corde
// UCS -> 100 MPa), pour dpr lineaire, dpr puissance (B 56,59, n 0,538), mc,
// et le controle QUI DOIT ECHOUER (puissance a n = 1, B = 6,727 = corde).
// Puis compression uniaxiale avec compDamage = crackband : l'aire adoucie
// nominale x h_e doit valoir compAc * compGIIc (crack band en compression).
// Retourne 0 si tous les criteres passent.
int triaxSelftest(const std::string& csvPath);
// Banc point-materiel de la loi cdp (2026-09-04) : pics triaxiaux et
// identite de translation, pic equibiaxial, traction (pic, energie, lc = 1 et
// 2 mm), decharges, dilatance hyperbolique, controle qui DOIT rater les
// donnees, carte inverse, rapport des meridiens = Kc, proprietes du retour,
// garde de snap-back. Retourne 0 si tout passe.
int cdpSelftest(const std::string& csvPath);

// Bancs falsifiants de l'endommagement a direction figee (tensionDamage =
// fixed, 2026-09-05) au point materiel, carte Red Bohus dpr (dpApex,
// dpTension = off, rankineDrive = stress) : (a) traction x jusqu'a d = 0,9,
// decharge, traction y : raideur E intacte (fixed) contre (1 - 0,9) E
// (scalar, REJETE) ; (b) traction x puis compression x : E dans les deux ;
// (c) traction a 45 deg : repere fige a 45 deg + objectivite ; (d) rotation
// des directions principales apres amorcage : cisaillement parasite pour
// beta = 1 et 0,1 (donnee) ; (e) energie dissipee = Gf/lc a 2 % (les deux).
// Retourne 0 si (a), (b), (c), (e) et la bit-identite scalar/cle passent.
int fixedCrackSelftest(const std::string& csvPath);

// Bancs de la loi dfhplus (2026-09-06, etape 1 du chantier DFH+), au point
// materiel : (1) sigma = d(rho psi)/d eps verifie contre les differences finies
// centrees ; (2) Y_i = -d(rho psi)/d D_i idem ; (3) reduction elastique exacte a
// D = 0 ; (4) objectivite sous rotation rigide ; (5) exposant de vitesse de
// l'obscuration 3/(m+3) pour m = 6, 12, 24 (dfhplus ET dpdfh, cote a cote) ;
// (6) comparaison dfhplus / dpdfh en traction, compression, triaxial et sur un
// chemin NON COAXIAL ; (7) controle qui DOIT echouer : sans l'ecretage de
// dilatance (dfhpPsiClamp = false) la dissipation plastique devient negative.
// Retourne 0 si (1)-(5) passent.
int dfhPlusSelftest(const std::string& csvPath);

// rockim matpoint <cfg> [out.csv] : pilote point-materiel generique (toutes
// les lois), chemins triax | tension | biaxial | uniaxial, cles mp* — le
// moteur de la calibration (quelques ms par confinement). Retourne 0.
int matpointDrive(const Config& cfg, const std::string& csvPath);


} // namespace rockim
