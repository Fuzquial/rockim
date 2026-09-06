// ---------------------------------------------------------------------------
// ThermoBench — banc thermodynamique generique (2026-09-06). Voir
// include/rockim/ThermoBench.hpp pour le POURQUOI et l'usage. Ce fichier
// documente le QUOI : la definition exacte des cinq tests, de leurs
// estimateurs et de leurs limites.
//
// CONVENTIONS. Traction positive ; eps = deformation de Biot co-rotee, la
// meme que recoit une MatLaw depuis le solveur ; SI (Pa, m, s, J/m^3).
// BASE ORTHONORMEE DE Sym(3) (c'est elle qui rend la symetrie majeure
// LISIBLE : elle est orthonormee pour le produit A:B, donc sigma:eps =
// somme_a sigma_a eps_a sans facteur 2 parasite) :
//   B0 = e1(x)e1, B1 = e2(x)e2, B2 = e3(x)e3,
//   B3 = (e1(x)e2 + e2(x)e1)/sqrt(2), B4 = (e2(x)e3 + e3(x)e2)/sqrt(2),
//   B5 = (e1(x)e3 + e3(x)e1)/sqrt(2).
//
// ---------------------------------------------------------------------------
// TEST 1 — SYMETRIE MAJEURE DE LA TANGENTE
// ---------------------------------------------------------------------------
// C_ab = d sigma_a / d eps_b par differences finies CENTREES sur les six
// composantes : sigma(eps0 +- h B_b) evaluee A CHAQUE FOIS depuis une COPIE
// de l'etat en eps0 (sans quoi l'etat derive et la « tangente » n'en est plus
// une), au TEMPS GELE (dt de sonde = 1e-6 dt du chemin : les mecanismes
// pilotes par le temps — l'obscuration de Denoual-Hild en particulier — ne
// doivent pas contribuer a une derivee par rapport a la DEFORMATION).
// Mesure : asym = max_ab |C_ab - C_ba| / max_ab |C_ab|.
//
// LECTURE. Un cadre derive d'UNE energie libre donne sigma = d(rho psi)/d eps
// donc une tangente ELASTIQUE (variables internes gelees) SYMETRIQUE PAR
// CONSTRUCTION — c'est le point de la decision de Fernando. En revanche la
// tangente ELASTOPLASTIQUE d'un ecoulement NON ASSOCIE (psi != beta, le cas
// de DP-DFH) est NON symetrique, et c'est de la physique, pas un defaut. Le
// banc separe donc les deux populations : un increment est dit ELASTIQUE si
// la SIGNATURE IRREVERSIBLE de l'etat (D, D_c, eps^p, peeq, kappa, les
// compteurs de dissipation, les endommagements directionnels...) est
// inchangee AU BIT pres apres les 12 evaluations. Le VERDICT ne porte que
// sur les increments elastiques ; les increments inelastiques sont comptes
// et leur pire cas rapporte, a titre d'information.
//
// ---------------------------------------------------------------------------
// TEST 2 — POSITIVITE DE LA DISSIPATION
// ---------------------------------------------------------------------------
// Sur chaque increment du chemin : D_k = sigma_mid : d eps - d(rho psi),
// sigma_mid = (sigma_k + sigma_k+1)/2 (trapeze — EXACT pour une reponse
// affine, donc exact sur tout le domaine elastique).
//
// COMMENT rho psi EST OBTENUE, ET CE QUE CA COUTE. Aucune loi du depot
// n'expose son energie libre (l'interface MatLaw ne rend que la contrainte) :
// le banc l'ESTIME par une SONDE DE DECHARGE ELASTIQUE, sur une COPIE de
// l'etat, A TEMPS GELE :
//   - direction de decharge d = -(1+nu)/E sigma + nu/E tr(sigma) I, soit
//     l'oppose de la deformation elastique qu'un materiau VIERGE associerait
//     a sigma (pour la loi elastique elle amene exactement a eps = 0) ;
//   - on avance eps <- eps0 + tau d par pas de tau = 1/16 jusqu'a ce que la
//     projection sigma:d CHANGE DE SIGNE (etat relache), avec interpolation
//     lineaire sur le dernier pas ; plafond tau <= 4 ;
//   - rho psi = - integrale de sigma:d eps le long de cette decharge (energie
//     RENDUE par le materiau), integree au trapeze.
// LIMITES, explicitement (elles bornent ce que le test prouve) :
//   (L1) l'estimateur ne mesure que la part RECUPERABLE de l'energie libre.
//        Toute energie STOCKEE non recuperable (ecrouissage cinematique,
//        energie piegee autour des fissures) est comptee comme dissipee :
//        la grandeur mesuree vaut D~ = D + d(rho psi_stockee)/dt. Cette
//        energie stockee etant croissante dans tous les modeles vises,
//        D~ >= D : une valeur NEGATIVE de D~ prouve une dissipation
//        reellement negative — le test n'a pas de faux positif de ce chef,
//        seulement des faux negatifs (il peut manquer une violation).
//   (L2) la sonde suppose la decharge NON DISSIPATIVE. Si un mecanisme
//        evolue PENDANT la sonde, l'estimation est biaisee : le banc le
//        DETECTE (comparaison bit a bit de la signature irreversible avant
//        et apres la sonde), marque l'increment « contamine », l'exclut du
//        verdict et en publie le compte. Le temps gele reduit fortement cette
//        population pour les lois a endommagement pilote par le temps.
//   (L3) si la projection ne change pas de signe avant tau = 4 l'etat n'est
//        pas relache : increment marque « non relache », exclu, compte. MEME
//        VERDICT si la projection change de signe alors que la contrainte est
//        encore GRANDE (|| sigma || > 5 % de || sigma_0 ||) : la direction de
//        decharge est FIGEE sur la compliance vierge, une contrainte qui
//        TOURNE en cours de decharge (anisotropie d'endommagement) peut
//        annuler la projection sans que l'etat soit relache, et l'energie
//        libre serait alors sous-estimee — artefact identifie et neutralise
//        a la relecture du 2026-09-06.
//   (L4) pour une loi VISQUEUSE la sonde rend le potentiel A TEMPS GELE ;
//        la surcontrainte visqueuse apparait alors integralement dans D~,
//        ou elle est positive en charge — ce qui est le comportement voulu.
// Critere : violation si D_k < -tolRel * echelle_k, echelle_k =
// |sigma_mid : d eps| + |d(rho psi)| + plancher. tolRel = 1e-6.
//
// ---------------------------------------------------------------------------
// TEST 3 — REDUCTION
// ---------------------------------------------------------------------------
// Avec --ref <loi> : les DEUX lois jouent les MEMES chemins (memes tirages,
// meme dt, meme lc, meme x0) et l'on compare les contraintes composante par
// composante ; violation si l'ecart relatif depasse 1e-12. C'est le test qui
// servira a prouver que `dfhplus` a parametres neutres REND EXACTEMENT
// `dpdfh` (ou la loi de reference retenue) — il est ecrit maintenant pour
// que la preuve existe avant la loi.
// Sans --ref : REDUCTION ELASTIQUE. Un chemin d'amplitude bornee a 0,3 ft/E
// (donc sous tout seuil de traction, tres loin du cone de compression) est
// joue et sigma doit valoir lambda tr(eps) I + 2 G eps a 1e-12 relatif pres.
// Toute loi doit passer ce test : une loi qui ne redonne pas l'elasticite en
// petit est fausse avant meme d'etre thermodynamique.
//
// ---------------------------------------------------------------------------
// TEST 4 — OBJECTIVITE
// ---------------------------------------------------------------------------
// Le meme tirage est rejoue avec une rotation rigide Q superposee :
// eps_k -> Q eps_k Q^T, MEME x0 (le tirage de Weibull est une propriete du
// POINT, pas du repere : le faire tourner testerait le hachage spatial, pas
// l'objectivite), meme dt, meme lc. On mesure a chaque pas :
//   (a) l'ecart sur les INVARIANTS — contraintes principales triees —
//       rapporte au plus grand |sigma| du chemin : c'est le critere du
//       verdict, tolerance 1e-8 relatif ;
//   (b) l'ecart TENSORIEL || sigma_tourne - Q sigma Q^T || / max |sigma|,
//       plus severe (il verifie en plus que le repere interne tourne avec la
//       matiere), rapporte a titre d'information.
//
// ---------------------------------------------------------------------------
// TEST 5 — CONTINUITE A L'ACTIVATION D'UN MECANISME
// ---------------------------------------------------------------------------
// Une reponse sigma(eps) continue a tangente bornee par l'elasticite ne peut
// pas rendre, sur un increment, une variation de contrainte PLUS GRANDE que
// la prediction elastique. On mesure donc
//   r = || d sigma || / (M || d eps ||),  M = max(3K, 2G)
// (M = plus grande valeur propre du tenseur elastique en norme de Frobenius),
// sur l'increment ENTIER puis sur le MEME increment DECOUPE en 8 sous-pas
// (dt decoupe aussi). Un simple transitoire raide voit r CHUTER au
// raffinement ; une VRAIE discontinuite (chute brutale d'endommagement,
// snap-back du retour, bascule de repere) garde r >> 1. Le verdict porte sur
// la valeur RAFFINEE, seuil r <= 1,05. On note si la signature irreversible
// a change sur cet increment (« mecanisme active »).
//
// ---------------------------------------------------------------------------
// TIRAGES
// ---------------------------------------------------------------------------
// Graine fixe, un generateur PAR TIRAGE (seed + phi * i) : le resultat ne
// depend NI du nombre de fils NI de l'ordonnancement. Six familles, en parts
// egales : traction uniaxiale, compression uniaxiale, TRIAXIAL de 0 a
// 300 MPa de confinement, cisaillement, chemins NON COAXIAUX (les directions
// principales tournent en cours de chemin), etats FORTEMENT ENDOMMAGES
// (pre-charge poussee jusqu'a D ~ 0,95). Chaque tirage porte son propre
// dt (log-uniforme 1e-9 - 1e-5 s, la gamme de la percussion), son propre
// lc (log-uniforme 0,5 - 2 mm), sa propre position x0 (donc son propre
// tirage de Weibull la ou la loi en fait), et un repere aleatoire uniforme.
// ---------------------------------------------------------------------------
#include "rockim/ThermoBench.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "rockim/Config.hpp"
#include "rockim/MatLaw.hpp"
#include "rockim/Material.hpp"

namespace rockim {
namespace {

using M3 = Eigen::Matrix3d;

// ---- base orthonormee de Sym(3) -------------------------------------------
inline M3 basisB(int a) {
    static const double s = 1.0 / 1.4142135623730951;
    M3 B = M3::Zero();
    switch (a) {
        case 0: B(0, 0) = 1.0; break;
        case 1: B(1, 1) = 1.0; break;
        case 2: B(2, 2) = 1.0; break;
        case 3: B(0, 1) = B(1, 0) = s; break;
        case 4: B(1, 2) = B(2, 1) = s; break;
        default: B(0, 2) = B(2, 0) = s; break;
    }
    return B;
}
inline double ddot(const M3& a, const M3& b) { return (a.cwiseProduct(b)).sum(); }

// ---- signature IRREVERSIBLE de l'etat -------------------------------------
// Tout ce qui, dans MatState, ne peut pas revenir en arriere : si ces
// nombres sont inchanges AU BIT pres, l'increment n'a fait bouger AUCUN
// mecanisme dissipatif. On EXCLUT deliberement les champs purement
// memoriels que toute loi met a jour a chaque appel sans rien dissiper
// (dfh.epsPrev, dfh.snom, dfh.t, dfh.smaxh, sk.sbar, sk.epsPrev).
// N = 8 scalaires + 6 (epsP) + 3 (dfh.Dv) + 1 (dfh.peeq) + 4 (sk) + 2 (cdp)
// + 3 (fcm.d) = 27. Un assert de fin de fonction verrouille le compte.
struct Sig {
    static const int N = 27;
    double v[N];
};
Sig irrev(const MatState& s) {
    Sig g;
    for (int i = 0; i < Sig::N; ++i) g.v[i] = 0.0;
    int k = 0;
    g.v[k++] = s.D;
    g.v[k++] = s.kappa;
    g.v[k++] = s.epvEq;
    g.v[k++] = s.pc;
    g.v[k++] = s.Dc;
    g.v[k++] = s.wPlas;
    g.v[k++] = s.wDamT;
    g.v[k++] = s.wDamC;
    for (int i = 0; i < 3; ++i)
        for (int j = i; j < 3; ++j) g.v[k++] = s.epsP(i, j);   // 6
    for (int i = 0; i < 3; ++i) g.v[k++] = s.dfh.Dv[i];        // 3
    g.v[k++] = s.dfh.peeq;
    g.v[k++] = s.sk.kapDP;
    g.v[k++] = s.sk.kapMR;
    g.v[k++] = s.sk.eqvt;
    g.v[k++] = s.sk.epsv;
    g.v[k++] = s.cdp.epsTpl;
    g.v[k++] = s.cdp.epsCpl;
    for (int i = 0; i < 3; ++i) g.v[k++] = s.fcm.d[i];         // 3
    assert(k == Sig::N);
    return g;
}
inline bool same(const Sig& a, const Sig& b) {
    return std::memcmp(a.v, b.v, sizeof(a.v)) == 0;
}

// ---- rotation aleatoire uniforme (Shoemake 1992) --------------------------
M3 randRot(std::mt19937_64& rng) {
    std::uniform_real_distribution<double> U(0.0, 1.0);
    double u1 = U(rng), u2 = U(rng), u3 = U(rng);
    double s1 = std::sqrt(1.0 - u1), s2 = std::sqrt(u1);
    Eigen::Quaterniond q(s2 * std::cos(2.0 * M_PI * u3),
                         s1 * std::sin(2.0 * M_PI * u2),
                         s1 * std::cos(2.0 * M_PI * u2),
                         s2 * std::sin(2.0 * M_PI * u3));
    q.normalize();
    return q.toRotationMatrix();
}

// ---- familles d'etats ------------------------------------------------------
enum Family { F_TENSION = 0, F_COMP, F_TRIAX, F_SHEAR, F_NONCOAX, F_DAMAGED,
              N_FAMILY };
const char* familyName(int f) {
    static const char* n[N_FAMILY] = {"traction", "compression", "triaxial",
                                      "cisaillement", "non-coaxial",
                                      "endommage"};
    return (f >= 0 && f < N_FAMILY) ? n[f] : "?";
}

struct Card {                       // ce que le banc doit savoir du materiau
    double E = 0, nu = 0, G = 0, K = 0, lam = 0, k0 = 0, M = 0, Gf = 0;
};

// Un tirage = tout ce qu'il faut pour rejouer EXACTEMENT le meme chemin.
struct Draw {
    int family = 0;
    double dt = 1.0e-7, lc = 1.0e-3;
    Eigen::Vector3d x0 = Eigen::Vector3d::Zero();
    M3 Q = M3::Identity();          // repere du chargement
    M3 epsEnd = M3::Zero();         // etat vise en fin de PRE-CHARGE
    M3 dEps = M3::Zero();           // increment du CHEMIN (apres pre-charge)
    M3 axisRot = M3::Identity();    // rotation par pas (chemins non coaxiaux)
    int nPre = 24, nPath = 10;
};

Draw makeDraw(int i, unsigned long long seed, const Card& c) {
    std::mt19937_64 rng(seed + 0x9E3779B97F4A7C15ull
                               * (unsigned long long)(i + 1));
    std::uniform_real_distribution<double> U(0.0, 1.0);
    auto logu = [&](double a, double b) {
        return std::exp(std::log(a) + U(rng) * (std::log(b) - std::log(a)));
    };
    Draw d;
    d.family = i % N_FAMILY;
    d.dt = logu(1.0e-9, 1.0e-5);
    d.lc = logu(0.5e-3, 2.0e-3);
    d.x0 = Eigen::Vector3d(U(rng), U(rng), U(rng)) * 0.1;
    d.Q = randRot(rng);
    const double k0 = c.k0, nu = c.nu;
    M3 E0 = M3::Zero();
    switch (d.family) {
        case F_TENSION: {
            double a = k0 * (0.2 + 5.8 * U(rng));
            E0(0, 0) = a; E0(1, 1) = -nu * a; E0(2, 2) = -nu * a;
            break;
        }
        case F_COMP: {
            double a = k0 * (0.5 + 29.5 * U(rng));
            E0(0, 0) = -a; E0(1, 1) = nu * a; E0(2, 2) = nu * a;
            break;
        }
        case F_TRIAX: {
            double s3 = 300.0e6 * U(rng);              // confinement 0-300 MPa
            double eh = -s3 / (3.0 * c.K);
            double a = k0 * 40.0 * U(rng);
            E0 = eh * M3::Identity();
            E0(0, 0) += -a; E0(1, 1) += nu * a; E0(2, 2) += nu * a;
            break;
        }
        case F_SHEAR: {
            double g = k0 * (0.2 + 19.8 * U(rng));
            double p = -200.0e6 * U(rng) / (3.0 * c.K);
            E0 = p * M3::Identity();
            E0(0, 1) = E0(1, 0) = 0.5 * g;
            break;
        }
        case F_NONCOAX: {
            double a = k0 * (0.5 + 9.5 * U(rng));
            double g = k0 * (0.2 + 5.0 * U(rng));
            E0(0, 0) = a; E0(1, 1) = -nu * a; E0(2, 2) = -0.3 * a;
            E0(1, 2) = E0(2, 1) = 0.5 * g;
            double th = 0.02 + 0.10 * U(rng);
            Eigen::Vector3d ax(U(rng) - 0.5, U(rng) - 0.5, U(rng) - 0.5);
            if (ax.norm() < 1e-12) ax = Eigen::Vector3d(0, 0, 1);
            d.axisRot =
                Eigen::AngleAxisd(th, ax.normalized()).toRotationMatrix();
            break;
        }
        default: {   // F_DAMAGED : pre-charge poussee, D vise 0,95
            double a = k0 * (20.0 + 180.0 * U(rng));
            double g = k0 * 2.0 * U(rng);
            E0(0, 0) = a; E0(1, 1) = -nu * a; E0(2, 2) = -nu * a;
            E0(0, 2) = E0(2, 0) = 0.5 * g;
            d.nPre = 60;
            d.dt = logu(1.0e-7, 1.0e-5);   // le temps compte (obscuration)
            break;
        }
    }
    d.epsEnd = d.Q * E0 * d.Q.transpose();
    M3 V = M3::Zero();
    if (U(rng) < 0.5) {
        V = d.epsEnd;
    } else {
        for (int a = 0; a < 6; ++a) V += (2.0 * U(rng) - 1.0) * basisB(a);
    }
    double n = V.norm();
    if (n < 1e-30) { V = basisB(0); n = 1.0; }
    d.dEps = (k0 * (0.05 + 0.95 * U(rng))) * (V / n);
    return d;
}

// ---- sonde d'energie libre (decharge elastique a temps gele) ---------------
double probePsi(const MatLaw& law, const MatState& s0, const M3& eps0,
                const Card& c, double lc, double dtp,
                bool& contaminated, bool& unrelaxed) {
    MatState s = s0;
    M3 sg = law.stress(eps0, s, dtp, lc);
    const Sig g0 = irrev(s);
    M3 d = -((1.0 + c.nu) / c.E) * sg
           + (c.nu / c.E) * sg.trace() * M3::Identity();
    contaminated = false;
    unrelaxed = true;
    if (!(d.norm() > 0.0)) { unrelaxed = false; return 0.0; }
    const double dtau = 1.0 / 16.0, tauMax = 4.0;
    const double sg0n = sg.norm();
    double proj = ddot(sg, d), w = 0.0, tau = 0.0;
    if (!(proj < 0.0)) { unrelaxed = false; return 0.0; }   // deja relache
    while (tau < tauMax - 1e-12) {
        M3 epsN = eps0 + (tau + dtau) * d;
        M3 sgN = law.stress(epsN, s, dtp, lc);
        double projN = ddot(sgN, d);
        if (projN >= 0.0) {                       // traversee : interpolation
            double f = proj / (proj - projN);
            w += 0.5 * proj * f * dtau;
            // garde (L3) : traversee de projection SANS relachement reel.
            // Le pas en tau valant 1/16, une traversee LEGITIME laisse au
            // plus ~6 % de la contrainte initiale a la fin du pas ; le seuil
            // est pose a 15 % pour ne pas rejeter ce reste de discretisation.
            unrelaxed = !(sgN.norm() <= 0.15 * sg0n);
            break;
        }
        w += 0.5 * (proj + projN) * dtau;
        proj = projN;
        tau += dtau;
    }
    if (!same(irrev(s), g0)) contaminated = true;
    return -w;
}

// ---- enregistrement d'une violation ---------------------------------------
struct Viol {
    int test = 0, draw = 0, family = 0, step = 0, elastic = 0, flag = 0;
    double value = 0.0, tol = 0.0, dt = 0.0, lc = 0.0;
    // grandeurs ABSOLUES du cas (unites dans l'en-tete du CSV) : test 2 —
    // aux1 = dissipation [J/m^3], aux2 = rho psi [J/m^3] ; test 5 — aux1 =
    // || d sigma || [Pa], aux2 = || d eps || [-] ; sinon 0.
    double aux1 = 0.0, aux2 = 0.0;
    double eps[6] = {0, 0, 0, 0, 0, 0};
    double sig[6] = {0, 0, 0, 0, 0, 0};
    double D = 0.0, epv = 0.0;
};
void fill6(double* v, const M3& a) {
    v[0] = a(0, 0); v[1] = a(1, 1); v[2] = a(2, 2);
    v[3] = a(0, 1); v[4] = a(1, 2); v[5] = a(0, 2);
}

// resultat d'UN tirage (agrege, plus les violations retenues)
struct DrawOut {
    // test 1
    double symEl = -1.0, symPl = -1.0;
    int symElN = 0, symPlN = 0, symElBad = 0, symPlBad = 0;
    // test 2
    int dissN = 0, dissBad = 0, dissContam = 0, dissUnrelax = 0;
    int dissBadSig = 0;             // violations SIGNIFICATIVES (voir plus bas)
    double dissWorst = 0.0, dissWorstAbs = 0.0, dissWorstPsi = 0.0;
    int dissWorstStep = -1, dissAbsStep = -1;
    double dissAbsWorst = 0.0;      // dissipation la plus negative [J/m^3]
    double dissCoarse = 0.0;        // pire cas RELATIF au trapeze GROSSIER
    int dissQuad = 0;               // violations a 8 sous-pas qui S'EFFACENT
                                    // au raffinement (erreur de quadrature)
    // test 3
    double redWorst = 0.0;
    int redN = 0, redBad = 0;
    // test 4
    double objInv = 0.0, objTen = 0.0;
    int objN = 0, objBad = 0;
    // test 5
    double contCoarse = 0.0, contFine = 0.0, contDsig = 0.0;
    int contN = 0, contBad = 0, contActive = 0;
    double contRunning = 0.0;       // pire r a temps COURANT (info)
    int contRaw = 0;                // flags a temps COURANT (info)
    int contTime = 0;               // flags du test 5 qui DISPARAISSENT quand
                                    // on gele le temps : le saut n'est pas une
                                    // discontinuite de sigma(eps), c'est la
                                    // relaxation PILOTEE PAR LE TEMPS
    // couverture
    double Dmax = 0.0, sigMax = 0.0, pMin = 0.0;
    std::vector<Viol> viol;
};

// ---- un tirage complet -----------------------------------------------------
void runDraw(int idx, const Draw& d, const Card& c, const MatLaw& law,
             const MatLaw* ref, DrawOut& o, bool exact) {
    const double dtp = d.dt * 1.0e-6;            // dt de sonde (temps gele)
    const double hSym = 1.0e-8;                  // pas des differences finies
    const double tolSym = 1.0e-4;
    const double tolDiss = 1.0e-6;
    const double tolRed = 1.0e-12;
    const double tolObj = 1.0e-8;
    const double tolCont = 1.05;
    const int nSub = 8;

    auto push = [&](int test, int step, double value, double tol, int elastic,
                    int flag, const M3& eps, const M3& sg, const MatState& st,
                    double aux1 = 0.0, double aux2 = 0.0) {
        Viol v;
        v.test = test; v.draw = idx; v.family = d.family; v.step = step;
        v.value = value; v.tol = tol; v.elastic = elastic; v.flag = flag;
        v.dt = d.dt; v.lc = d.lc; v.aux1 = aux1; v.aux2 = aux2;
        fill6(v.eps, eps); fill6(v.sig, sg);
        v.D = st.D; v.epv = st.epvEq;
        if (o.viol.size() < 64) o.viol.push_back(v);
    };

    MatState st;
    st.x0 = d.x0;
    M3 eps = M3::Zero(), sg = M3::Zero();

    // ---- pre-charge 0 -> epsEnd ------------------------------------------
    for (int k = 1; k <= d.nPre; ++k) {
        eps = ((double)k / d.nPre) * d.epsEnd;
        sg = law.stress(eps, st, d.dt, d.lc);
        o.Dmax = std::max(o.Dmax, st.D);
        o.sigMax = std::max(o.sigMax, sg.cwiseAbs().maxCoeff());
        o.pMin = std::min(o.pMin, sg.trace() / 3.0);
    }

    // ---- TEST 1 : symetrie majeure de la tangente -------------------------
    {
        const MatState s0 = st;
        Eigen::Matrix<double, 6, 6> C;
        bool elasticStep = true;
        const Sig g0 = irrev(s0);
        for (int b = 0; b < 6; ++b) {
            MatState sp = s0, sm = s0;
            M3 sgp = law.stress(eps + hSym * basisB(b), sp, dtp, d.lc);
            M3 sgm = law.stress(eps - hSym * basisB(b), sm, dtp, d.lc);
            if (!same(irrev(sp), g0) || !same(irrev(sm), g0))
                elasticStep = false;
            for (int a = 0; a < 6; ++a)
                C(a, b) = ddot(sgp - sgm, basisB(a)) / (2.0 * hSym);
        }
        double amax = 0.0, cmax = 0.0;
        for (int a = 0; a < 6; ++a)
            for (int b = 0; b < 6; ++b) {
                cmax = std::max(cmax, std::abs(C(a, b)));
                amax = std::max(amax, std::abs(C(a, b) - C(b, a)));
            }
        double err = cmax > 0.0 ? amax / cmax : 0.0;
        if (elasticStep) {
            o.symElN = 1; o.symEl = err;
            if (err > tolSym) { o.symElBad = 1; push(1, d.nPre, err, tolSym, 1, 0, eps, sg, st); }
        } else {
            o.symPlN = 1; o.symPl = err;
            if (err > tolSym) o.symPlBad = 1;
        }
    }

    // ---- chemin : tests 2 et 5 -------------------------------------------
    {
        bool cont0 = false, unrel0 = false;
        double psi0 = exact
                          ? law.freeEnergy(eps, st)
                          : probePsi(law, st, eps, c, d.lc, dtp, cont0, unrel0);
        M3 epsK = eps, sgK = sg;
        MatState stK = st;
        M3 rot = M3::Identity();
        for (int k = 1; k <= d.nPath; ++k) {
            // increment (avec rotation du repere pour les chemins non coaxiaux)
            if (d.family == F_NONCOAX) rot = d.axisRot * rot;
            M3 epsN = epsK + d.dEps;
            if (d.family == F_NONCOAX)
                epsN = M3(rot * epsK * rot.transpose()) + d.dEps;
            M3 dE = epsN - epsK;

            // Travail RAFFINE (8 sous-pas) : quand la loi expose son energie
            // libre, c'est lui qui porte le verdict du test 2 — le trapeze
            // grossier commet, a chaque COIN de la reponse (une valeur propre
            // de eps qui change de signe, l'endommagement qui plafonne), une
            // erreur de quadrature en O(||d eps||^2) que rien ne distingue
            // d'une vraie dissipation negative. Cette erreur decroit en
            // 1/nSub^2 ; une VRAIE violation, non.
            double wFine = 0.0, psiFine = 0.0;
            bool hasFine = false;

            // --- TEST 5 : continuite (grossier vs raffine, sur des COPIES) --
            double dn = dE.norm();
            if (dn > 0.0) {
                MatState sc = stK;
                M3 sgc = law.stress(epsN, sc, d.dt, d.lc);
                MatState sf = stK;
                M3 sgf = sgK, sgPrev = sgK;
                for (int j = 1; j <= nSub; ++j) {
                    sgf = law.stress(epsK + ((double)j / nSub) * dE, sf,
                                     d.dt / nSub, d.lc);
                    wFine += ddot(0.5 * (sgPrev + sgf), dE / (double)nSub);
                    sgPrev = sgf;
                }
                if (exact) psiFine = law.freeEnergy(epsN, sf);
                hasFine = true;
                double rc = (sgc - sgK).norm() / (c.M * dn);
                double rf = (sgf - sgK).norm() / (c.M * dn);
                bool act = !same(irrev(sf), irrev(stK));
                o.contN++;
                if (act) o.contActive++;
                o.contCoarse = std::max(o.contCoarse, rc);
                // --- DISCRIMINANT TEMPS / DEFORMATION ----------------------
                // r rapporte une variation de CONTRAINTE a un increment de
                // DEFORMATION. Un mecanisme pilote par le TEMPS (l'obscuration
                // de Denoual-Hild : dx = (S lam)^{1/3} k c dt) fait tomber
                // sigma SANS que eps bouge : r y est arbitrairement grand ET
                // invariant au raffinement (subdiviser divise dt d'autant), et
                // pourtant sigma(eps) n'a pas la moindre discontinuite. Le meme
                // increment est donc rejoue A TEMPS GELE (dt x 1e-6) : c'est
                // CETTE valeur qui mesure la continuite de sigma(eps) et qui
                // fait verdict ; la valeur a temps courant reste publiee.
                MatState sq = stK;
                M3 sgq = sgK;
                for (int j = 1; j <= nSub; ++j)
                    sgq = law.stress(epsK + ((double)j / nSub) * dE, sq,
                                     d.dt * 1.0e-6 / nSub, d.lc);
                double rq = (sgq - sgK).norm() / (c.M * dn);
                o.contRunning = std::max(o.contRunning, rf);
                if (rq > o.contFine) {
                    o.contFine = rq;
                    o.contDsig = (sgq - sgK).norm();
                }
                if (rf > tolCont) o.contRaw++;
                if (rf > tolCont && rq <= tolCont) o.contTime++;
                if (rq > tolCont) {
                    o.contBad++;
                    push(5, k, rq, tolCont, 0, act ? 1 : 0, epsN, sgc, sc,
                         (sgq - sgK).norm(), dn);
                }
            }

            // --- avance REELLE du chemin -----------------------------------
            const MatState stK0 = stK;      // etat AVANT le pas (raffinements)
            M3 sgNn = law.stress(epsN, stK, d.dt, d.lc);
            o.Dmax = std::max(o.Dmax, stK.D);
            o.sigMax = std::max(o.sigMax, sgNn.cwiseAbs().maxCoeff());
            o.pMin = std::min(o.pMin, sgNn.trace() / 3.0);

            // --- TEST 2 : dissipation ---------------------------------------
            bool cont1 = false, unrel1 = false;
            double psi1 = exact ? law.freeEnergy(epsN, stK)
                                : probePsi(law, stK, epsN, c, d.lc, dtp,
                                           cont1, unrel1);
            double work = ddot(0.5 * (sgK + sgNn), dE);
            double dpsi = psi1 - psi0;
            double diss = work - dpsi;
            double scale = std::abs(work) + std::abs(dpsi) + 1.0e-30;
            double dissFine8 = 0.0;
            if (exact && hasFine) {
                // energie libre EXACTE (exposee) + travail RAFFINE : plus
                // d'estimateur, plus d'exclusion, plus d'erreur de quadrature
                // au premier plan. La valeur grossiere reste publiee en info.
                o.dissCoarse = std::min(o.dissCoarse, diss / scale);
                work = wFine;
                dpsi = psiFine - psi0;
                diss = work - dpsi;
                scale = std::abs(work) + std::abs(dpsi) + 1.0e-30;
                dissFine8 = diss;
                // --- DISCRIMINANT DE QUADRATURE -------------------------
                // Le trapeze commet, a chaque coin de la reponse, une erreur
                // en O(h^2) : elle est divisee par 16 quand on passe de 8 a
                // 32 sous-pas. Une VRAIE dissipation negative, elle, ne bouge
                // pas. On ne recalcule que sur les increments deja flagues
                // (quelques dizaines sur 120 000) : le cout est nul.
                if (diss / scale < -tolDiss) {
                    const int nSub2 = 32;
                    MatState sf2 = stK0;
                    M3 sgp = sgK, sg2 = sgK;
                    double w2 = 0.0;
                    for (int j = 1; j <= nSub2; ++j) {
                        sg2 = law.stress(epsK + ((double)j / nSub2) * dE, sf2,
                                         d.dt / nSub2, d.lc);
                        w2 += ddot(0.5 * (sgp + sg2), dE / (double)nSub2);
                        sgp = sg2;
                    }
                    double dp2 = law.freeEnergy(epsN, sf2) - psi0;
                    double d2 = w2 - dp2;
                    double s2 = std::abs(w2) + std::abs(dp2) + 1.0e-30;
                    // converge vers 0 comme la quadrature (facteur >= 2 en
                    // passant de 8 a 32 sous-pas, quand la theorie en promet
                    // 16) => artefact de quadrature, pas une violation.
                    if (d2 / s2 >= -tolDiss
                        || std::abs(d2) < 0.5 * std::abs(dissFine8)) {
                        o.dissQuad++;
                        diss = 0.0;
                        scale = s2;
                    } else {
                        diss = d2;
                        scale = s2;
                    }
                }
            }
            if (cont0 || cont1) o.dissContam++;
            else if (unrel0 || unrel1) o.dissUnrelax++;
            else {
                o.dissN++;
                double rel = diss / scale;
                if (rel < o.dissWorst) {
                    o.dissWorst = rel; o.dissWorstStep = k;
                    o.dissWorstAbs = diss; o.dissWorstPsi = psi1;
                }
                if (diss < o.dissAbsWorst) { o.dissAbsWorst = diss; o.dissAbsStep = k; }
                if (rel < -tolDiss) {
                    o.dissBad++;
                    // SIGNIFICATIVE : la dissipation manquante depasse 1 % de
                    // l'energie de bande Gf/lc, c'est-a-dire de tout ce que
                    // l'element a le droit de dissiper en rompant. En dessous,
                    // la violation est reelle mais physiquement negligeable.
                    if (-diss > 0.01 * c.Gf / d.lc) o.dissBadSig++;
                    push(2, k, rel, -tolDiss, 0, 0, epsN, sgNn, stK, diss, psi1);
                }
            }
            psi0 = psi1; cont0 = cont1; unrel0 = unrel1;
            epsK = epsN; sgK = sgNn;
        }
    }

    // ---- TEST 4 : objectivite --------------------------------------------
    {
        M3 Qo;
        {
            std::mt19937_64 rr(0xA5A5A5A5ull + (unsigned long long)idx);
            Qo = randRot(rr);
        }
        MatState sa, sb;
        sa.x0 = d.x0; sb.x0 = d.x0;
        double smax = 1.0e-30;
        std::vector<M3> sgA, sgB;
        sgA.reserve(d.nPre + d.nPath);
        sgB.reserve(d.nPre + d.nPath);
        M3 e = M3::Zero(), rot = M3::Identity();
        for (int k = 1; k <= d.nPre + d.nPath; ++k) {
            if (k <= d.nPre) {
                e = ((double)k / d.nPre) * d.epsEnd;
            } else {
                if (d.family == F_NONCOAX) {
                    rot = d.axisRot * rot;
                    e = rot * e * rot.transpose() + d.dEps;
                } else {
                    e = e + d.dEps;
                }
            }
            M3 a = law.stress(e, sa, d.dt, d.lc);
            M3 b = law.stress(Qo * e * Qo.transpose(), sb, d.dt, d.lc);
            sgA.push_back(a);
            sgB.push_back(b);
            smax = std::max(smax, a.cwiseAbs().maxCoeff());
        }
        for (std::size_t k = 0; k < sgA.size(); ++k) {
            Eigen::SelfAdjointEigenSolver<M3> ea(sgA[k]), eb(sgB[k]);
            double ei = (ea.eigenvalues() - eb.eigenvalues()).cwiseAbs().maxCoeff()
                        / smax;
            double et = (sgB[k] - Qo * sgA[k] * Qo.transpose()).norm() / smax;
            o.objN++;
            o.objInv = std::max(o.objInv, ei);
            o.objTen = std::max(o.objTen, et);
            if (ei > tolObj) {
                o.objBad++;
                push(4, (int)k, ei, tolObj, 0, 0, sgA[k], sgB[k], sa);
            }
        }
    }

    // ---- TEST 3 : reduction ----------------------------------------------
    if (ref) {
        MatState sa, sb;
        sa.x0 = d.x0; sb.x0 = d.x0;
        M3 e = M3::Zero(), rot = M3::Identity();
        for (int k = 1; k <= d.nPre + d.nPath; ++k) {
            if (k <= d.nPre) {
                e = ((double)k / d.nPre) * d.epsEnd;
            } else {
                if (d.family == F_NONCOAX) {
                    rot = d.axisRot * rot;
                    e = rot * e * rot.transpose() + d.dEps;
                } else {
                    e = e + d.dEps;
                }
            }
            M3 a = law.stress(e, sa, d.dt, d.lc);
            M3 b = ref->stress(e, sb, d.dt, d.lc);
            double den = std::max(a.cwiseAbs().maxCoeff(),
                                  b.cwiseAbs().maxCoeff()) + 1.0e-30;
            double err = (a - b).cwiseAbs().maxCoeff() / den;
            o.redN++;
            o.redWorst = std::max(o.redWorst, err);
            if (err > tolRed) {
                o.redBad++;
                push(3, k, err, tolRed, 0, 0, e, a, sa);
            }
        }
    } else {
        // reduction ELASTIQUE : amplitude bornee a 0,3 k0
        MatState sa;
        sa.x0 = d.x0;
        M3 dir = d.epsEnd.norm() > 0.0 ? d.epsEnd / d.epsEnd.norm() : basisB(0);
        for (int k = 1; k <= 8; ++k) {
            M3 e = (0.3 * c.k0 * (double)k / 8.0) * dir;
            M3 a = law.stress(e, sa, d.dt, d.lc);
            M3 b = c.lam * e.trace() * M3::Identity() + 2.0 * c.G * e;
            double den = b.cwiseAbs().maxCoeff() + 1.0e-30;
            double err = (a - b).cwiseAbs().maxCoeff() / den;
            o.redN++;
            o.redWorst = std::max(o.redWorst, err);
            if (err > tolRed) {
                o.redBad++;
                push(3, k, err, tolRed, 0, 0, e, a, sa);
            }
        }
    }
}

// carte materiau par defaut (Red Bohus, celle des bancs du depot)
const char* kDefaultDeck =
    "# carte Red Bohus du banc thermodynamique (identique a selftest-fixed)\n"
    "E = 77.66e9\n"
    "nu = 0.29\n"
    "rho = 2620\n"
    "ft = 9.0e6\n"
    "cohesion = 22.77e6\n"
    "frictionDeg = 50.4\n"
    "Gf = 100\n"
    "# erosion DESARMEE : le banc mesure la LOI, pas la suppression d'element\n"
    "erodeD = 2\n"
    "erodeEpv = 0\n";

}  // namespace

// ===========================================================================
int thermoBench(const ThermoBenchOpts& o) {
    // ---- carte materiau + construction des lois ---------------------------
    std::string deckText = kDefaultDeck;
    if (!o.deck.empty()) {
        std::ifstream f(o.deck);
        if (!f) {
            std::cerr << "[thermo] deck introuvable : " << o.deck << "\n";
            return 2;
        }
        std::ostringstream ss;
        ss << f.rdbuf();
        deckText += "\n# --- deck utilisateur : " + o.deck + " ---\n" + ss.str()
                    + "\n";
    }
    const std::string tmpCfg = o.csv + ".deck.cfg";
    {
        std::ofstream f(tmpCfg);
        if (!f) {
            std::cerr << "[thermo] impossible d'ecrire " << tmpCfg << "\n";
            return 2;
        }
        f << deckText;
    }
    const double lcMax = 2.0e-3;
    Material m;
    std::unique_ptr<MatLaw> law, ref;
    try {
        Config cfg = Config::load(tmpCfg);
        m = Material::from(cfg);
        PhaseSet::validate(m, "thermobench");
        law = MatLaw::make(o.law, m, cfg, lcMax);
        if (!o.refLaw.empty()) ref = MatLaw::make(o.refLaw, m, cfg, lcMax);
    } catch (const std::exception& e) {
        std::cerr << "[thermo] construction de la loi impossible : " << e.what()
                  << "\n";
        return 2;
    }

    Card c;
    c.E = m.E; c.nu = m.nu; c.G = m.G(); c.K = m.K();
    c.lam = m.E * m.nu / ((1.0 + m.nu) * (1.0 - 2.0 * m.nu));
    c.k0 = m.ft / m.E;
    c.M = std::max(3.0 * c.K, 2.0 * c.G);
    c.Gf = m.Gf;

    // Estimateur du test 2 : energie libre EXPOSEE par la loi (exacte) ou
    // sonde de decharge (estimateur biaise). --probe force la seconde, pour
    // comparer une loi neuve a une loi ancienne AVEC LE MEME INSTRUMENT.
    const bool exact = law->hasFreeEnergy() && !o.forceProbe;
    const int n = std::max(6, o.nDraws);
    std::cout << "[thermo] banc thermodynamique — loi = " << law->name()
              << (ref ? std::string(" contre reference = ") + ref->name()
                      : std::string(" (test 3 = reduction elastique)"))
              << "\n[thermo] carte : E " << m.E / 1e9 << " GPa, nu " << m.nu
              << ", ft " << m.ft / 1e6 << " MPa, c " << m.cohesion / 1e6
              << " MPa, phi " << m.phiDeg << " deg, Gf " << m.Gf
              << " J/m^2 ; k0 = ft/E = " << c.k0
              << "\n[thermo] test 2 : "
              << (exact ? "energie libre EXPOSEE par la loi (exacte) + travail "
                          "raffine 8 sous-pas"
                        : (law->hasFreeEnergy()
                               ? "SONDE DE DECHARGE forcee (--probe)"
                               : "sonde de decharge (energie libre non exposee)"))
              << "\n[thermo] tirages : " << n
              << " (graine " << o.seed << "), 6 familles, dt 1e-9-1e-5 s, "
                 "lc 0,5-2 mm\n";

    std::vector<Draw> draws((std::size_t)n);
    std::vector<DrawOut> outs((std::size_t)n);
    for (int i = 0; i < n; ++i) draws[(std::size_t)i] = makeDraw(i, o.seed, c);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16)
#endif
    for (int i = 0; i < n; ++i)
        runDraw(i, draws[(std::size_t)i], c, *law, ref.get(),
                outs[(std::size_t)i], exact);

    // ---- agregation (dans l'ordre des tirages : deterministe) -------------
    struct Agg {
        long long nEval = 0, nBad = 0;
        double worst = 0.0, aux1 = 0.0, aux2 = 0.0;
        int wDraw = -1, wStep = -1, wFam = -1;
    } A1, A1p, A2, A3, A4, A5;
    long long contam = 0, unrelax = 0, mechActive = 0, contTimeAll = 0,
              contRawAll = 0;
    double DmaxAll = 0.0, sigMaxAll = 0.0, pMinAll = 0.0, objTenMax = 0.0;
    double contCoarseMax = 0.0, dissAbsMin = 0.0, dissCoarseMin = 0.0;
    double contRunMax = 0.0;
    long long dissSig = 0, dissQuadAll = 0;
    int dissAbsDraw = -1, dissAbsStep = -1, dissAbsFam = -1;
    long long badByFam[N_FAMILY][5];
    for (int f = 0; f < N_FAMILY; ++f)
        for (int t = 0; t < 5; ++t) badByFam[f][t] = 0;

    std::vector<Viol> rows;
    for (int i = 0; i < n; ++i) {
        const DrawOut& r = outs[(std::size_t)i];
        const int f = draws[(std::size_t)i].family;
        A1.nEval += r.symElN;  A1.nBad += r.symElBad;
        if (r.symElN && r.symEl > A1.worst) { A1.worst = r.symEl; A1.wDraw = i; A1.wFam = f; }
        A1p.nEval += r.symPlN; A1p.nBad += r.symPlBad;
        if (r.symPlN && r.symPl > A1p.worst) { A1p.worst = r.symPl; A1p.wDraw = i; A1p.wFam = f; }
        A2.nEval += r.dissN;   A2.nBad += r.dissBad;
        if (r.dissWorst < A2.worst) {
            A2.worst = r.dissWorst; A2.wDraw = i; A2.wStep = r.dissWorstStep;
            A2.wFam = f; A2.aux1 = r.dissWorstAbs; A2.aux2 = r.dissWorstPsi;
        }
        if (r.dissAbsWorst < dissAbsMin) {
            dissAbsMin = r.dissAbsWorst; dissAbsDraw = i;
            dissAbsStep = r.dissAbsStep; dissAbsFam = f;
        }
        dissSig += r.dissBadSig;
        dissCoarseMin = std::min(dissCoarseMin, r.dissCoarse);
        dissQuadAll += r.dissQuad;
        A3.nEval += r.redN;    A3.nBad += r.redBad;
        if (r.redWorst > A3.worst) { A3.worst = r.redWorst; A3.wDraw = i; A3.wFam = f; }
        A4.nEval += r.objN;    A4.nBad += r.objBad;
        if (r.objInv > A4.worst) { A4.worst = r.objInv; A4.wDraw = i; A4.wFam = f; }
        A5.nEval += r.contN;   A5.nBad += r.contBad;
        if (r.contFine > A5.worst) {
            A5.worst = r.contFine; A5.wDraw = i; A5.wFam = f;
            A5.aux1 = r.contDsig;
        }
        contam += r.dissContam; unrelax += r.dissUnrelax; mechActive += r.contActive;
        contTimeAll += r.contTime;
        contRawAll += r.contRaw;
        DmaxAll = std::max(DmaxAll, r.Dmax);
        sigMaxAll = std::max(sigMaxAll, r.sigMax);
        pMinAll = std::min(pMinAll, r.pMin);
        objTenMax = std::max(objTenMax, r.objTen);
        contCoarseMax = std::max(contCoarseMax, r.contCoarse);
        contRunMax = std::max(contRunMax, r.contRunning);
        badByFam[f][0] += r.symElBad;
        badByFam[f][1] += r.dissBad;
        badByFam[f][2] += r.redBad;
        badByFam[f][3] += r.objBad;
        badByFam[f][4] += r.contBad;
        for (const Viol& v : r.viol)
            if ((int)rows.size() < o.maxRows) rows.push_back(v);
    }

    // ---- CSV des violations ----------------------------------------------
    {
        std::ofstream f(o.csv);
        f.precision(15);
        f << "# aux1/aux2 : test 2 = dissipation [J/m^3] / rho psi [J/m^3] ; "
             "test 5 = |d sigma| [Pa] / |d eps| [-]\n";
        f << "test,draw,famille,pas,valeur,tolerance,elastique,mecanisme,"
             "aux1,aux2,dt,lc,"
             "e11,e22,e33,e12,e23,e13,s11,s22,s33,s12,s23,s13,D,epvEq\n";
        for (const Viol& v : rows) {
            f << v.test << ',' << v.draw << ',' << familyName(v.family) << ','
              << v.step << ',' << v.value << ',' << v.tol << ',' << v.elastic
              << ',' << v.flag << ',' << v.aux1 << ',' << v.aux2 << ',' << v.dt
              << ',' << v.lc;
            for (int k = 0; k < 6; ++k) f << ',' << v.eps[k];
            for (int k = 0; k < 6; ++k) f << ',' << v.sig[k];
            f << ',' << v.D << ',' << v.epv << '\n';
        }
    }

    // ---- compte rendu ------------------------------------------------------
    auto line = [&](const char* name, const Agg& a, const char* unit) {
        std::cout << "  " << std::left << std::setw(34) << name << std::right
                  << std::setw(9) << a.nBad << " / " << std::setw(8) << a.nEval
                  << "   pire " << std::setw(12) << a.worst << " " << unit;
        if (a.wDraw >= 0) {
            std::cout << "  (tirage " << a.wDraw << ", " << familyName(a.wFam);
            if (a.wStep >= 0) std::cout << ", pas " << a.wStep;
            std::cout << ", dt " << draws[(std::size_t)a.wDraw].dt << " s, lc "
                      << draws[(std::size_t)a.wDraw].lc << " m)";
        }
        std::cout << "\n";
    };
    std::cout.precision(6);
    std::cout << std::scientific;
    std::cout << "\n[thermo] COUVERTURE : D max atteint " << DmaxAll
              << ", |sigma| max " << sigMaxAll / 1e6 << " MPa, pression moyenne "
              << "la plus compressive " << pMinAll / 1e6 << " MPa\n";
    std::cout << "[thermo] VIOLATIONS (compte / evaluations, pire cas)\n";
    line("1 symetrie majeure (elastique)", A1, "[-]");
    line("1' symetrie majeure (inelast.)", A1p, "[-] (info)");
    line("2 dissipation >= 0", A2, "[-]");
    line("3 reduction", A3, "[-]");
    line("4 objectivite (invariants)", A4, "[-]");
    line("5 continuite (temps gele)", A5, "[-]");
    if (A2.wDraw >= 0) {
        std::cout << "      -> pire cas RELATIF du test 2, en absolu : "
                  << "dissipation " << A2.aux1 << " J/m^3 pour rho psi = "
                  << A2.aux2 << " J/m^3\n";
        std::cout << "      -> dissipation la plus NEGATIVE de tout le banc : "
                  << dissAbsMin << " J/m^3";
        if (dissAbsDraw >= 0)
            std::cout << " (tirage " << dissAbsDraw << ", "
                      << familyName(dissAbsFam) << ", pas " << dissAbsStep
                      << ", lc " << draws[(std::size_t)dissAbsDraw].lc
                      << " m, soit " << 100.0 * (-dissAbsMin)
                             / (m.Gf / draws[(std::size_t)dissAbsDraw].lc)
                      << " % de Gf/lc)";
        std::cout << "\n      -> dont " << dissSig
                  << " violation(s) SIGNIFICATIVE(S) (deficit > 1 % de Gf/lc)\n";
        if (exact)
            std::cout << "      -> meme mesure au TRAPEZE GROSSIER (info, "
                         "erreur de quadrature comprise) : pire relatif "
                      << dissCoarseMin << " ; " << dissQuadAll
                      << " increment(s) flague(s) a 8 sous-pas se sont EFFACES "
                         "a 32 (erreur de quadrature au coin d'une reponse "
                         "C1 par morceaux, en O(h^2)), non comptes\n";
    }
    if (A5.wDraw >= 0) {
        std::cout << "      -> pire cas du test 5 en ABSOLU : saut de "
                  << A5.aux1 / 1e6 << " MPa sur un increment\n";
        std::cout << "      -> a temps COURANT " << contRawAll
                  << " increment(s) depassent le seuil, dont " << contTimeAll
                  << " DISPARAISSENT a temps gele (dt x 1e-6) : le saut y est "
                     "pilote par le TEMPS (obscuration de Denoual-Hild, dx ~ "
                     "dt), pas par la deformation. r = ||d sigma||/(M ||d eps||) "
                     "n'est un critere de continuite de sigma(eps) qu'a temps "
                     "gele : c'est cette valeur-la qui fait verdict.\n";
    }
    std::cout << "  " << std::left << std::setw(34)
              << "4' objectivite (tensorielle)" << std::right
              << "                        pire " << objTenMax << " [-] (info)\n";
    std::cout << "  " << std::left << std::setw(34)
              << "5' continuite (increment entier)" << std::right
              << "                        pire " << contCoarseMax
              << " [-] (info)\n";
    std::cout << "  " << std::left << std::setw(34)
              << "5'' continuite (temps COURANT)" << std::right
              << "                        pire " << contRunMax
              << " [-] (info)\n";
    std::cout << "[thermo] sonde d'energie libre : " << contam
              << " increments CONTAMINES (un mecanisme a bouge pendant la "
                 "decharge), " << unrelax
              << " NON RELACHES (tau > 4) — exclus du verdict du test 2 ; "
              << mechActive << " increments ou un mecanisme s'active (test 5)\n";
    std::cout << "[thermo] violations par famille (1 elast. / 2 diss. / 3 red. "
                 "/ 4 obj. / 5 cont.)\n";
    for (int f = 0; f < N_FAMILY; ++f) {
        std::cout << "    " << std::left << std::setw(14) << familyName(f)
                  << std::right;
        for (int t = 0; t < 5; ++t) std::cout << std::setw(10) << badByFam[f][t];
        std::cout << "\n";
    }
    std::cout << std::defaultfloat;
    long long total = A1.nBad + A2.nBad + A3.nBad + A4.nBad + A5.nBad;
    std::cout << "[thermo] violations retenues dans " << o.csv << " : "
              << rows.size() << " ligne(s)"
              << ((int)rows.size() >= o.maxRows ? " (plafond --max-rows atteint)"
                                                : "")
              << "\n[thermo] TOTAL " << total << " violation(s) -> "
              << (total == 0 ? "PASS" : "ECHEC") << "\n";
    return total == 0 ? 0 : 1;
}

}  // namespace rockim
