// ---------------------------------------------------------------------------
// MatLawDfhPlus — la loi `law = dfhplus` (2026-09-06, ETAPE 1 du chantier
// DFH+). Voir docs/DFHPLUS_etape1.md pour le compte rendu chiffre.
//
// CE QUE CE FICHIER EST, ET CE QU'IL N'EST PAS
// -------------------------------------------
// `dfhplus` reprend le PERIMETRE PHYSIQUE de `dpdfh` (plasticite de
// Drucker-Prager non associee, endommagement de traction anisotrope a trois
// directions figees a l'amorcage, cinetique d'obscuration de Denoual-Hild avec
// tirages de Weibull par element, unilateralite) mais elle est construite sur
// un CADRE DIFFERENT. Elle n'est PAS un portage de la VUMAT : `dpdfh` l'est et
// le reste, intacte. Les mecanismes (B) ouverture axiale directionnelle et
// (C) pulverisation de `CONTINUUM/loi_dfh_plus/loi_DFH_plus.tex` ne sont PAS
// codes ici : l'etape 1 prouve que le cadre tient sur le perimetre connu.
//
// POURQUOI UN CADRE NEUF (relecture adverse de loi_DFH_plus.pdf)
// --------------------------------------------------------------
// Le cadre de DP-DFH est SUR-DETERMINE : contrainte effective a la Lemaitre,
// encadrement en racine de l'equivalence en energie, et energie libre ecrite
// independamment — trois enonces pour deux libertes. Le banc thermodynamique
// (`rockim thermobench dpdfh`) le mesure : 12,1 % des increments juges y ont
// une dissipation NEGATIVE, les cinq pires deficits sont tous sur des chemins
// NON COAXIAUX, et la reponse y est discontinue a saturation de
// l'endommagement. Decision de Fernando (2026-09-06) : postuler UNE energie
// libre et tout en deriver.
//
// ===========================================================================
// 1. L'ENERGIE LIBRE POSTULEE
// ===========================================================================
// Notation de loi_DFH_plus.tex : \epse = deformation elastique, D_i^t = les
// trois endommagements de traction portes par le repere FIGE {n_1, n_2, n_3},
// \lambda et G les coefficients de Lame, \Cel le tenseur elastique sain.
//
// DECOMPOSITION SPECTRALE (Miehe et al. 2010) de la deformation elastique :
//
//     \epse = \sum_a \varepsilon_a  m_a (x) m_a ,
//     \epse^{+/-} = \sum_a \Mac{\varepsilon_a}_{+/-}  m_a (x) m_a ,
//
// {m_a} = repere PRINCIPAL courant de \epse (il tourne avec le chargement ; il
// n'a pas a coincider avec le repere fige {n_i}, et c'est precisement la ou
// DP-DFH perdait la positivite).
//
// TENSEUR D'INTEGRITE porte par le repere fige :
//
//     \A = \sum_i (1 - D_i^t)  n_i (x) n_i ,
//     g(m) = m . \A m = \sum_i (1 - D_i^t) (n_i . m)^2   dans [0, 1]
//
// (g(n_i) = 1 - D_i^t exactement ; g est l'integrite VUE par une facette de
// normale m, c'est-a-dire la fraction d'aire portante qui reste quand trois
// familles de fissures orthogonales sont ouvertes.)
//
// INTEGRITE VOLUMIQUE, moyenne HARMONIQUE des trois integrites :
//
//     g_v(D_1, D_2, D_3) = 3 / [ 1/(1-D_1) + 1/(1-D_2) + 1/(1-D_3) ]
//
// C'est le couplage en SERIE : une dilatation isotrope doit etre transmise a
// travers les TROIS familles de facettes, les compliances s'ajoutent, et le
// ressort casse des qu'un maillon casse (g_v -> 0 des qu'un D_i -> 1). Aucune
// autre moyenne ne rend a la fois g_v = 1-D pour D_1=D_2=D_3=D (reduction au
// cas isotrope) et g_v -> 0 pour un seul D_i -> 1 (contrainte NULLE en
// traction uniaxiale saturee, comme la VUMAT).
//
//     +-----------------------------------------------------------------+
//     | rho psi(\epse, D_1, D_2, D_3, {n_i}) =                          |
//     |     (lambda/2) [ g_v <tr \epse>_+^2 + <tr \epse>_-^2 ]          |
//     |   + G \A : (\epse^+)^2                                          |
//     |   + G || \epse^- ||^2                                           |
//     +-----------------------------------------------------------------+
//
// La deuxieme ligne se lit aussi G \sum_a g(m_a) <\varepsilon_a>_+^2 :
// l'identite \A : (\epse^+)^2 = \sum_a g(m_a) <\varepsilon_a>_+^2 vient de
// (\epse^+)^2 = \sum_a <\varepsilon_a>_+^2 m_a (x) m_a. Ecrite sous la forme
// \A : (\epse^+)^2 elle est manifestement objective et se derive sans jamais
// deriver les vecteurs propres.
//
// TROIS PROPRIETES PAR CONSTRUCTION :
//  (a) D = 0 : g_v = 1, \A = I, rho psi = (lambda/2)(tr \epse)^2 + G \epse:\epse
//      = l'elasticite lineaire EXACTEMENT (les deux moities spectrales se
//      recollent). La reduction elastique est exacte, pas approchee.
//  (b) UNILATERALITE naturelle : seule la partie POSITIVE est degradee. Une
//      fissure fermee (\varepsilon_a < 0) transmet tout. Aucun interrupteur.
//  (c) SYMETRIE MAJEURE de la tangente elastique : sigma est un gradient, la
//      tangente est une hessienne. Ce n'est pas verifie a posteriori, c'est
//      impossible autrement.
//
// ===========================================================================
// 2. LA CONTRAINTE = d(rho psi) / d \epse
// ===========================================================================
//     sigma = lambda [ g_v <tr \epse>_+ + <tr \epse>_- ] I
//           + \Pplus : [ G (\A \epse^+ + \epse^+ \A) ]
//           + 2 G \epse^-
//
// ou \Pplus = d\epse^+/d\epse est la projection spectrale positive (Daleckii-
// Krein) : dans le repere propre de \epse, (\Pplus : T)_{aa} = H(\varepsilon_a)
// T_{aa} et (\Pplus : T)_{ab} = theta_{ab} T_{ab} avec
// theta_{ab} = (<\varepsilon_a>_+ - <\varepsilon_b>_+)/(\varepsilon_a -
// \varepsilon_b) dans [0,1], prolonge par H a la degenerescence. \Pplus a la
// symetrie majeure (c'est la hessienne de \epse -> (1/2)||\epse^+||^2), ce qui
// autorise l'ecriture ci-dessus.
// Verification : \A = I donne \Pplus : (2G \epse^+) = 2G \epse^+, donc
// sigma = lambda tr(\epse) I + 2G \epse. Le banc le mesure a 1e-12.
//
// ===========================================================================
// 3. LES FORCES MOTRICES Y_i = - d(rho psi) / d D_i
// ===========================================================================
//     Y_i = (lambda/2) [ g_v^2 / (3 (1-D_i)^2) ] <tr \epse>_+^2
//         + G || \epse^+ n_i ||^2                              >= 0
//
// (le premier terme vient de dg_v/dD_i = - g_v^2 / (3 (1-D_i)^2) ; le second
// de d\A/dD_i = - n_i (x) n_i contracte avec (\epse^+)^2.)
// POSITIVITE INCONDITIONNELLE : les deux termes sont des carres multiplies par
// lambda > 0 et G > 0. Comme la cinetique d'obscuration est monotone
// (dD_i >= 0), la dissipation d'endommagement \sum_i Y_i dD_i est positive
// SANS AUCUNE CONDITION — c'est le point que DP-DFH ne pouvait pas tenir.
// BORNE : g_v <= 3(1-D_i) toujours, donc g_v^2/(1-D_i)^2 <= 9 : Y_i reste
// borne quand D_i -> 1 (pas de singularite a saturation, contrairement au
// Y_i^o en (1-D)^{-1/2} de la variante R de loi_DFH_plus.tex).
//
// CONTRAINTE EQUIVALENTE D'ENERGIE. La cinetique de Denoual-Hild est ecrite
// en CONTRAINTE ; on la branche sur Y_i par l'equivalence en energie
//
//     sigma_i^eq = sqrt( 2 E Y_i / c_nu ),
//     c_nu = (2/(1+nu)) [ nu(1-2nu)/6 + 1/2 ]        (0,8067 a nu = 0,29)
//
// c_nu est choisi pour que, en TRACTION UNIAXIALE libre laterale et a D = 0,
// sigma_i^eq soit EXACTEMENT la contrainte effective \bar\sigma_{ii}. La
// cinetique et ses tirages de Weibull sont alors inchanges mot pour mot, et
// l'exposant de vitesse 3/(m+3) est conserve (banc (5) du selftest).
//
// ===========================================================================
// 4. LA PLASTICITE, ET SUR QUELLE CONTRAINTE
// ===========================================================================
// La contrainte EFFECTIVE du cadre est celle qui rend le squelette SAIN :
//
//     \bar\sigma := d(rho psi)/d\epse |_{D = 0} = \Cel : \epse .
//
// C'est le SEUL choix compatible avec le cadre : (i) c'est une derivee de la
// meme energie libre, prise a endommagement nul ; (ii) c'est la seule mesure
// de contrainte INDEPENDANTE de D, donc la seule qui donne une surface de
// charge dont la position ne bouge pas quand la roche se fissure — le
// frottement se joue sur les contacts intacts, pas sur l'aire perdue ;
// (iii) elle est LINEAIRE en \epse, donc le critere f = q - p tan(beta) - d
// et le retour radial de `dpdfh` s'y transposent SANS AUCUNE MODIFICATION.
// C'est aussi ce que fait la VUMAT (retour sur la contrainte effective) : la
// reponse plastique de dfhplus et celle de dpdfh coincident tant qu'aucune
// direction n'est endommagee.
//
// ECRETAGE DE DILATANCE (cle `dfhpPsiClamp`, defaut ACTIF). Avec une energie
// libre unique, la dissipation vaut D = sigma^nom : d\epsp + \sum_i Y_i dD_i :
// c'est la contrainte NOMINALE, pas l'effective, qui est appariee a
// d\epsp. Avec l'ecoulement DP non associe
// d\epsp = dgamma [ (3/(2 \bar q)) \bar s + (tan Psi / 3) I ] on obtient
//
//     D_A / dgamma = q_proj - tan(Psi) p^nom,
//     q_proj := (3/(2 \bar q)) \bar s : s^nom,   p^nom := -tr(sigma^nom)/3,
//
// qui peut etre NEGATIF des que l'endommagement est fortement anisotrope. On
// impose donc tan(Psi) <= max(q_proj, 0) / p^nom quand p^nom > 0 : c'est
// exactement l'ecretage d'admissibilite eq:T:clamp / eq:D1:psibound de
// loi_DFH_plus.tex, transpose au cone. Sur la carte Red Bohus il est INACTIF
// partout ou l'endommagement est faible (sur la surface \bar q/\bar p > tan
// beta = 1,266 alors que tan Psi = 0,268) : il ne touche ni les triaxiaux de
// calibration ni la compression. `dfhpPsiClamp = false` le desarme — c'est la
// VARIANTE QUI DOIT ECHOUER du selftest.
//
// ===========================================================================
// 5. CE QUI CHANGE PAR RAPPORT A dpdfh, ET POURQUOI
// ===========================================================================
//  * Traction uniaxiale saturee : dpdfh rend sigma_11 = (1-D) \bar\sigma_11 ;
//    dfhplus rend sigma_11/(E eps) = [nu g_v + (1-D)] / (1+nu), plus RAIDE de
//    nu g_v/(1+nu) (nul a D = 0 et a D = 1, maximal a mi-chemin). C'est le
//    prix du couplage volumique : le terme lambda ne peut pas etre degrade
//    direction par direction sans casser la symetrie majeure.
//  * Contrainte laterale residuelle en traction : le decoupage spectral laisse
//    sigma_22 = (E nu/(1+nu)) (g_v - 1) eps < 0 (artefact connu du split de
//    Miehe, nul a D = 0).
//  * Cisaillement : la VUMAT postule un facteur min(f_i, f_j) ; ici le facteur
//    tombe du potentiel et vaut la moyenne ARITHMETIQUE des integrites des
//    deux directions. Ce n'est pas un choix, c'est une consequence.
//  * Compression et triaxiaux : IDENTIQUES a dpdfh tant qu'aucune direction ne
//    s'amorce (elasticite exacte + meme retour DP sur la meme effective).
// ---------------------------------------------------------------------------
#include "rockim/MatLawDfhPlus.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>

#include <Eigen/Dense>

namespace rockim {
namespace dfhp {

using M3 = Eigen::Matrix3d;
using V3 = Eigen::Vector3d;

// Plafond d'endommagement, celui de la VUMAT (DCAP) : garde le potentiel fini.
constexpr double DCAP = 0.9999;

struct Props {
    double E = 77.66e9, nu = 0.29, rho = 2620.0;
    double betaDeg = 51.7, dcoh = 153.3e6, psiDeg = 15.0;
    double m = 24.0, sigw = 120.0e6, zeff = 1.0e-9;
    double k = 0.38, S = 4.18879;
    double deld = 1.0e9;
    bool psiVar = false;
    double psi0 = 160.345, kpsi = 0.213793, psiMax = 51.7;
    bool psiClamp = true;          // ecretage d'admissibilite de la dilatance
    int volInteg = 0;              // 0 = harmonique, 1 = min, 2 = aucun (g_v = 1)
};

// --- tirage de Weibull : MEME hachage spatial que dpdfh (kst_seed) ----------
// Transliteration a l'identique (et non un appel : `dpdfh` ne doit etre ni
// modifie ni rendu dependant de ce fichier). Les deux lois tirent donc
// EXACTEMENT les memes seuils pour le meme x0, ce qui rend la comparaison
// dfhplus / dpdfh honnete.
inline void seed(double x, double y, double z, double sigk, double m,
                 double sc[3]) {
    auto q = [](double v) {
        return (uint64_t)(int64_t)std::llround(v * 1.0e6);
    };
    uint64_t h = (q(x) * 73856093ull) ^ (q(y) * 19349663ull);
    h ^= q(z) * 83492791ull;
    h ^= 1234567891234567891ull;
    if (h == 0ull) h = 88172645463325252ull;
    for (int i = 0; i < 3; ++i) {
        h ^= h << 13;
        h ^= h >> 7;
        h ^= h << 17;
        double u = ((double)(h >> 11) + 0.5) / 9007199254740992.0;
        if (u < 1.0e-12) u = 1.0e-12;
        if (u > 1.0 - 1.0e-12) u = 1.0 - 1.0e-12;
        sc[i] = sigk * std::pow(-std::log(1.0 - u), 1.0 / m);
    }
    for (int i = 0; i < 2; ++i)
        for (int j = i + 1; j < 3; ++j)
            if (sc[j] < sc[i]) std::swap(sc[i], sc[j]);
}

// --- repere fige <-> Euler ZYX (memes formules que dpdfh : kst_eulr/kst_reul)
inline void eulrOf(const M3& R, double eul[3]) {
    double cb = std::sqrt(R(2, 1) * R(2, 1) + R(2, 2) * R(2, 2));
    eul[1] = std::atan2(-R(2, 0), cb);
    if (cb > 1.0e-9) {
        eul[0] = std::atan2(R(1, 0), R(0, 0));
        eul[2] = std::atan2(R(2, 1), R(2, 2));
    } else {
        eul[0] = std::atan2(-R(0, 1), R(1, 1));
        eul[2] = 0.0;
    }
}
inline M3 reulOf(const double eul[3]) {
    double ca = std::cos(eul[0]), sa = std::sin(eul[0]);
    double cb = std::cos(eul[1]), sb = std::sin(eul[1]);
    double cg = std::cos(eul[2]), sg = std::sin(eul[2]);
    M3 R;
    R << ca * cb, ca * sb * sg - sa * cg, ca * sb * cg + sa * sg,
         sa * cb, sa * sb * sg + ca * cg, sa * sb * cg - ca * sg,
         -sb,     cb * sg,                cb * cg;
    return R;
}

inline double pos(double x) { return x > 0.0 ? x : 0.0; }
inline double neg(double x) { return x < 0.0 ? x : 0.0; }

// --- decomposition spectrale + projection \Pplus ----------------------------
struct Spec {
    V3 ev;            // valeurs propres de \epse (ordre croissant, Eigen)
    M3 Q;             // vecteurs propres en COLONNES
    M3 epsPos, epsNeg;
};

inline Spec spectral(const M3& e) {
    Spec s;
    Eigen::SelfAdjointEigenSolver<M3> es(e);
    s.ev = es.eigenvalues();
    s.Q = es.eigenvectors();
    V3 p(pos(s.ev(0)), pos(s.ev(1)), pos(s.ev(2)));
    M3 ep = s.Q * p.asDiagonal() * s.Q.transpose();
    s.epsPos = 0.5 * (ep + ep.transpose());
    s.epsNeg = e - s.epsPos;
    return s;
}

// (\Pplus : T) — T symetrique. Formule de Daleckii-Krein pour la derivee de la
// fonction tensorielle isotrope \epse -> \epse^+.
inline M3 projPlus(const Spec& s, const M3& T) {
    M3 Tt = s.Q.transpose() * T * s.Q;
    M3 Rt = M3::Zero();
    const double scale = s.ev.cwiseAbs().maxCoeff() + 1.0e-30;
    for (int a = 0; a < 3; ++a)
        Rt(a, a) = (s.ev(a) > 0.0 ? 1.0 : 0.0) * Tt(a, a);
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) {
            if (a == b) continue;
            double da = s.ev(a) - s.ev(b);
            double th;
            if (std::abs(da) > 1.0e-12 * scale)
                th = (pos(s.ev(a)) - pos(s.ev(b))) / da;
            else
                th = (0.5 * (s.ev(a) + s.ev(b)) > 0.0) ? 1.0 : 0.0;
            Rt(a, b) = th * Tt(a, b);
        }
    M3 R = s.Q * Rt * s.Q.transpose();
    return 0.5 * (R + R.transpose());
}

// --- integrite volumique et sa derivee --------------------------------------
inline double gVol(const double D[3], int mode) {
    if (mode == 2) return 1.0;
    double a[3];
    for (int i = 0; i < 3; ++i) a[i] = 1.0 - std::min(D[i], DCAP);
    if (mode == 1) return std::min(a[0], std::min(a[1], a[2]));
    return 3.0 / (1.0 / a[0] + 1.0 / a[1] + 1.0 / a[2]);
}
// -dg_v/dD_i  (>= 0)
inline double dGVol(const double D[3], int i, int mode) {
    if (mode == 2) return 0.0;
    double a[3];
    for (int j = 0; j < 3; ++j) a[j] = 1.0 - std::min(D[j], DCAP);
    if (mode == 1) {
        double mn = std::min(a[0], std::min(a[1], a[2]));
        return (a[i] <= mn + 1.0e-15) ? 1.0 : 0.0;
    }
    double gv = 3.0 / (1.0 / a[0] + 1.0 / a[1] + 1.0 / a[2]);
    return gv * gv / (3.0 * a[i] * a[i]);
}

} // namespace dfhp

namespace {

class DfhPlusLaw : public MatLaw {
public:
    DfhPlusLaw(const Material& m, const dfhp::Props& p) : MatLaw(m), P_(p) {
        cnu_ = (2.0 / (1.0 + P_.nu))
               * (P_.nu * (1.0 - 2.0 * P_.nu) / 6.0 + 0.5);
        cwav_ = std::sqrt(P_.E / P_.rho);
        tanb_ = std::tan(P_.betaDeg * M_PI / 180.0);
    }

    std::string name() const override { return "dfhplus"; }
    bool hasFreeEnergy() const override { return true; }

    // rho psi(eps, etat) — variables internes FIGEES. Exacte, pas estimee.
    double freeEnergy(const Eigen::Matrix3d& eps,
                      const MatState& s) const override {
        if (s.dfhp.dead) return 0.0;
        return freeEnergyRaw(eps - s.epsP, s.dfhp);
    }

    int damageForces(const Eigen::Matrix3d& eps, const MatState& s,
                     double* Y) const override {
        drivingForces(eps - s.epsP, s.dfhp, Y);
        return 3;
    }

    // Y_i = -d(rho psi)/dD_i, force motrice de la direction i du repere fige.
    void drivingForces(const Eigen::Matrix3d& epsE, const MatState::Dfhp& st,
                       double Y[3]) const {
        dfhp::Spec sp = dfhp::spectral(epsE);
        const double tr = epsE.trace();
        const double pp = dfhp::pos(tr);
        Eigen::Matrix3d R = st.frozen ? dfhp::reulOf(st.eul)
                                      : Eigen::Matrix3d::Identity();
        for (int i = 0; i < 3; ++i) {
            Eigen::Vector3d ni = R.col(i);
            double v = (sp.epsPos * ni).squaredNorm();
            Y[i] = 0.5 * lam_ * dfhp::dGVol(st.Dv, i, P_.volInteg) * pp * pp
                   + G_ * v;
        }
    }

    Eigen::Matrix3d stress(const Eigen::Matrix3d& eps, MatState& s, double dt,
                           double lc) const override {
        auto& st = s.dfhp;
        const double dtv = std::max(dt, 1.0e-30);
        double Vel = lc * lc * lc;
        if (Vel < 1.0e-12) Vel = 1.0;
        if (!st.seeded) {
            st.seeded = true;
            double sigk = P_.sigw * std::pow(P_.zeff / Vel, 1.0 / P_.m)
                          * s.ftScale;
            dfhp::seed(s.x0.x(), s.x0.y(), s.x0.z(), sigk, P_.m, st.sc);
        }
        if (st.dead) return Eigen::Matrix3d::Zero();
        st.t += dtv;

        // ---- 1. predicteur elastique sur la contrainte EFFECTIVE ----------
        Eigen::Matrix3d eeTr = eps - s.epsP;
        Eigen::Matrix3d sbar = elastic(eeTr);          // \Cel : \epse^trial

        // ---- 2. retour de Drucker-Prager (effective, compression seule) ----
        plasticReturn(sbar, s, eeTr);

        // ---- 3. amorcage / croissance de l'endommagement -------------------
        Eigen::Matrix3d ee = eps - s.epsP;
        damageUpdate(ee, st, Vel, s.ftScale, dtv);

        // ---- 4. contrainte NOMINALE = d(rho psi)/d\epse --------------------
        Eigen::Matrix3d sig = nominal(ee, st);

        // ---- 5. sorties partagees ------------------------------------------
        s.D = std::max(st.Dv[0], std::max(st.Dv[1], st.Dv[2]));
        s.epvEq = st.peeq;
        s.kappa = st.smaxh;
        s.wDamT = st.wDam;
        st.psi = freeEnergyRaw(ee, st);
        if (s.D >= P_.deld) {
            st.dead = true;
            s.eroded = true;
            s.eroCode = 4;
            return Eigen::Matrix3d::Zero();
        }
        return sig;
    }

private:
    Eigen::Matrix3d integrity(const MatState::Dfhp& st) const {
        if (!st.frozen) return Eigen::Matrix3d::Identity();
        Eigen::Matrix3d R = dfhp::reulOf(st.eul);
        Eigen::Vector3d a(1.0 - std::min(st.Dv[0], dfhp::DCAP),
                          1.0 - std::min(st.Dv[1], dfhp::DCAP),
                          1.0 - std::min(st.Dv[2], dfhp::DCAP));
        return R * a.asDiagonal() * R.transpose();
    }

    double freeEnergyRaw(const Eigen::Matrix3d& ee,
                         const MatState::Dfhp& st) const {
        dfhp::Spec sp = dfhp::spectral(ee);
        Eigen::Matrix3d A = integrity(st);
        const double tr = ee.trace();
        const double gv = dfhp::gVol(st.Dv, P_.volInteg);
        const double pp = dfhp::pos(tr), pn = dfhp::neg(tr);
        return 0.5 * lam_ * (gv * pp * pp + pn * pn)
               + G_ * (A * sp.epsPos * sp.epsPos).trace()
               + G_ * sp.epsNeg.squaredNorm();
    }

    // sigma = d(rho psi)/d\epse, analytique.
    Eigen::Matrix3d nominal(const Eigen::Matrix3d& ee,
                            const MatState::Dfhp& st) const {
        dfhp::Spec sp = dfhp::spectral(ee);
        Eigen::Matrix3d A = integrity(st);
        const double tr = ee.trace();
        const double gv = dfhp::gVol(st.Dv, P_.volInteg);
        Eigen::Matrix3d T = G_ * (A * sp.epsPos + sp.epsPos * A);
        T = 0.5 * (T + T.transpose());
        Eigen::Matrix3d sig =
            lam_ * (gv * dfhp::pos(tr) + dfhp::neg(tr))
                * Eigen::Matrix3d::Identity()
            + dfhp::projPlus(sp, T) + 2.0 * G_ * sp.epsNeg;
        return 0.5 * (sig + sig.transpose());
    }

    // Retour DP non associe sur \bar\sigma, avec ecretage d'admissibilite de
    // la dilatance. `sbar` entre en predicteur et sort retournee ; s.epsP est
    // mis a jour par \epsp += \Cel^{-1} : (\bar\sigma^tr - \bar\sigma).
    void plasticReturn(Eigen::Matrix3d& sbar, MatState& s,
                       const Eigen::Matrix3d& eeTr) const {
        auto& st = s.dfhp;
        const double I1 = sbar.trace();
        const double pbar = -I1 / 3.0;
        Eigen::Matrix3d dev = sbar + pbar * Eigen::Matrix3d::Identity();
        const double J2 = 0.5 * dev.squaredNorm();
        const double q = std::sqrt(3.0 * std::max(J2, 0.0));
        const double f = q - pbar * tanb_ - P_.dcoh;
        if (!(f > 0.0 && pbar > 0.0)) return;

        double tanp = std::tan(P_.psiDeg * M_PI / 180.0);
        if (P_.psiVar) {
            double ps = P_.psi0 - P_.kpsi * (pbar * 1.0e-6);
            ps = std::min(std::max(ps, 0.0), P_.psiMax);
            tanp = std::tan(ps * M_PI / 180.0);
        }
        const double xK = K_;
        // --- ecretage : tan Psi <= max(q_proj, 0) / p^nom quand p^nom > 0 ---
        // La dissipation qui doit rester positive est celle du PAS ENTIER,
        // ecrite avec la contrainte nominale MOYENNE (regle du trapeze, celle
        // du banc thermodynamique) : un ecretage pose sur la seule contrainte
        // de DEBUT de pas annule D_A au premier ordre et laisse le second
        // ordre le rendre negatif. On fait donc un PREDICTEUR-CORRECTEUR :
        // premier retour avec tan Psi brut -> \epsp de fin de pas -> nominale
        // de fin de pas -> moyenne -> ecretage -> retour definitif.
        if (P_.psiClamp && st.frozen && q > 1.0e-12) {
            Eigen::Matrix3d sn0 = nominal(eeTr, st);
            for (int pass = 0; pass < 2; ++pass) {
                double dl = f / (3.0 * G_ + xK * tanb_ * tanp);
                double qq = q - 3.0 * G_ * dl;
                Eigen::Matrix3d sr;
                if (qq > 0.0) {
                    sr = (qq / q) * dev
                         - (pbar + xK * tanp * dl)
                               * Eigen::Matrix3d::Identity();
                } else {
                    sr = (P_.dcoh / tanb_) * Eigen::Matrix3d::Identity();
                    dl = q / (3.0 * G_);
                }
                Eigen::Matrix3d dSp = sbar - sr;
                double trp = dSp.trace();
                Eigen::Matrix3d dep =
                    (dSp - (trp / 3.0) * Eigen::Matrix3d::Identity())
                        / (2.0 * G_)
                    + (trp / (9.0 * K_)) * Eigen::Matrix3d::Identity();
                Eigen::Matrix3d sn1 = nominal(eeTr - dep, st);
                Eigen::Matrix3d sm = 0.5 * (sn0 + sn1);
                const double pnom = -sm.trace() / 3.0;
                if (!(pnom > 0.0)) break;
                Eigen::Matrix3d sdev = sm + pnom * Eigen::Matrix3d::Identity();
                const double qproj = 1.5 * (dev.cwiseProduct(sdev)).sum() / q;
                const double tmax = std::max(qproj, 0.0) / pnom;
                if (tanp > tmax) {
                    tanp = tmax;
                    st.clamp = tmax;
                    if (pass == 0) ++st.nClamp;
                } else {
                    break;
                }
            }
        }
        double dlam = f / (3.0 * G_ + xK * tanb_ * tanp);
        const double qn = q - 3.0 * G_ * dlam;
        double facd, pnew;
        if (qn > 0.0 && q > 1.0e-12) {
            facd = qn / q;
            pnew = pbar + xK * tanp * dlam;
        } else {
            facd = 0.0;
            pnew = -P_.dcoh / tanb_;
            dlam = q / (3.0 * G_);
        }
        Eigen::Matrix3d sret =
            facd * dev - pnew * Eigen::Matrix3d::Identity();
        Eigen::Matrix3d dS = sbar - sret;
        const double trD = dS.trace();
        Eigen::Matrix3d dSdev =
            dS - (trD / 3.0) * Eigen::Matrix3d::Identity();
        s.epsP += dSdev / (2.0 * G_)
                  + (trD / (9.0 * K_)) * Eigen::Matrix3d::Identity();
        st.peeq += dlam;
        sbar = sret;
    }

    // Amorcage (repere fige au premier depassement) + croissance par la
    // cinetique d'obscuration de Denoual-Hild, pilotee par sigma_i^eq(Y_i).
    void damageUpdate(const Eigen::Matrix3d& ee, MatState::Dfhp& st,
                      double Vel, double ftScale, double dt) const {
        const double THIRD = 1.0 / 3.0;
        // --- gel du repere : directions principales de \epse (= celles de
        // \bar\sigma tant que D = 0, elasticite isotrope), ordre max/mid/min
        if (!st.frozen) {
            dfhp::Spec sp = dfhp::spectral(ee);
            const double tr = ee.trace(), pp = dfhp::pos(tr);
            const double volTerm =
                0.5 * lam_ * dfhp::dGVol(st.Dv, 0, P_.volInteg) * pp * pp;
            double seqMax = 0.0;
            for (int a = 0; a < 3; ++a) {
                const double v = dfhp::pos(sp.ev(a));
                const double Ya = volTerm + G_ * v * v;
                const double sq =
                    std::sqrt(std::max(2.0 * P_.E * Ya / cnu_, 0.0));
                if (sq > seqMax) seqMax = sq;
            }
            if (seqMax > st.smaxh) st.smaxh = seqMax;
            if (seqMax >= st.sc[0]) {
                Eigen::Matrix3d R;
                R.col(0) = sp.Q.col(2);         // max (Eigen trie croissant)
                R.col(1) = sp.Q.col(1);         // mid
                R.col(2) = sp.Q.col(0);         // min
                if (R.determinant() < 0.0) R.col(2) = -R.col(2);
                dfhp::eulrOf(R, st.eul);
                st.frozen = true;
                st.ti[0] = std::max(st.t, 1.0e-30);
            }
        }
        if (!st.frozen) return;
        double Y[3];
        drivingForces(ee, st, Y);
        for (int i = 0; i < 3; ++i) {
            const double seq =
                std::sqrt(std::max(2.0 * P_.E * Y[i] / cnu_, 0.0));
            if (seq > st.smaxh) st.smaxh = seq;
            if (st.ti[i] <= 0.0 && seq >= st.sc[i])
                st.ti[i] = std::max(st.t, 1.0e-30);
            if (st.ti[i] > 0.0 && seq > 0.0 && st.Dv[i] < dfhp::DCAP) {
                const double sigwLoc = P_.sigw * ftScale;
                double xlam = std::pow(seq / sigwLoc, P_.m) / P_.zeff;
                if (xlam * Vel < 1.0) xlam = 1.0 / Vel;
                double xx = std::pow(-std::log(1.0 - st.Dv[i]), THIRD);
                xx += std::pow(P_.S * xlam, THIRD) * P_.k * cwav_ * dt;
                const double Dnew =
                    std::min(1.0 - std::exp(-xx * xx * xx), dfhp::DCAP);
                st.wDam += Y[i] * std::max(Dnew - st.Dv[i], 0.0);
                st.Dv[i] = Dnew;
            }
        }
    }

    dfhp::Props P_;
    double cnu_ = 1.0, cwav_ = 0.0, tanb_ = 0.0;
};

} // namespace

std::unique_ptr<MatLaw> makeDfhPlusLaw(const Material& m, const Config& c) {
    dfhp::Props p;
    p.E = m.E;
    p.nu = m.nu;
    p.rho = m.rho;
    // MEMES CLES que dpdfh : les cartes sont interchangeables.
    p.betaDeg = c.getd("dfhBetaDeg", 51.7);
    p.dcoh = c.getd("dfhDCoh", 153.3e6);
    p.psiDeg = c.getd("dfhPsiDeg", 15.0);
    p.m = c.getd("dfhWeibullM", 24.0);
    p.sigw = c.getd("dfhSigW", 120.0e6);
    p.zeff = c.getd("dfhZeff", 1.0e-9);
    p.k = c.getd("dfhK", 0.38);
    p.S = c.getd("dfhS", 4.18879);
    p.deld = c.getd("dfhDeld", 1.0e9);
    p.psiVar = c.getb("dfhPsiVar", false);
    p.psi0 = c.getd("dfhPsi0", 160.345);
    p.kpsi = c.getd("dfhKPsi", 0.213793);
    p.psiMax = c.getd("dfhPsiMax", 51.7);
    // ---- cles PROPRES au cadre neuf ------------------------------------
    // dfhpPsiClamp : ecretage d'admissibilite de la dilatance (cf. §4 de
    //   l'en-tete). Defaut TRUE. false = la variante qui DOIT echouer.
    // dfhpVolInteg : forme de l'integrite volumique g_v — harmonic (defaut,
    //   couplage en serie), min (borne inferieure, non lisse) ou none
    //   (g_v = 1 : le terme lambda n'est pas degrade — ablation qui isole le
    //   couplage volumique, celui que le contre-exemple de loi_DFH_plus.tex
    //   §T:B accuse).
    p.psiClamp = c.getb("dfhpPsiClamp", true);
    {
        std::string vi = c.gets("dfhpVolInteg", "harmonic");
        if (vi == "harmonic") p.volInteg = 0;
        else if (vi == "min") p.volInteg = 1;
        else if (vi == "none") p.volInteg = 2;
        else throw std::runtime_error("dfhpVolInteg must be harmonic | min | "
                                      "none (got '" + vi + "')");
    }
    if (!(p.m > 1.0) || !(p.sigw > 0.0) || !(p.zeff > 0.0) || !(p.k > 0.0)
        || !(p.S > 0.0) || !(p.dcoh > 0.0)
        || !(p.betaDeg > 0.0 && p.betaDeg < 89.0))
        throw std::runtime_error("dfhplus: requires dfhWeibullM > 1, "
                                 "dfhSigW/dfhZeff/dfhK/dfhS/dfhDCoh > 0, "
                                 "dfhBetaDeg in (0, 89)");
    if (!(m.nu > 0.0 && m.nu < 0.5))
        throw std::runtime_error("dfhplus: nu must be in (0, 0.5) — le cadre "
                                 "suppose lambda > 0 (positivite de Y_i)");
    return std::unique_ptr<MatLaw>(new DfhPlusLaw(m, p));
}

} // namespace rockim
