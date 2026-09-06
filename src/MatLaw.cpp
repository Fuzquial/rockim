// ---------------------------------------------------------------------------
// MatLaw — the fem3d constitutive laws. See the header for the model cards.
// ---------------------------------------------------------------------------
#include "rockim/MatLaw.hpp"
#include "rockim/MatLawDfhPlus.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace rockim {

namespace {

// max principal value of a symmetric 3x3
double maxPrincipal(const Eigen::Matrix3d& S) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(S);
    return es.eigenvalues().maxCoeff();
}

// Euler ZYX <-> matrice de rotation (colonnes = n1 n2 n3), memes formules
// que dfhk::eulr / reul (repere fige de la DP-DFH) — pour tensionDamage = fixed
void fcmEuler(const Eigen::Matrix3d& R, double eul[3]) {
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
Eigen::Matrix3d fcmFrame(const double eul[3]) {
    double ca = std::cos(eul[0]), sa = std::sin(eul[0]);
    double cb = std::cos(eul[1]), sb = std::sin(eul[1]);
    double cg = std::cos(eul[2]), sg = std::sin(eul[2]);
    Eigen::Matrix3d R;
    R(0, 0) = ca * cb;  R(0, 1) = ca * sb * sg - sa * cg;  R(0, 2) = ca * sb * cg + sa * sg;
    R(1, 0) = sa * cb;  R(1, 1) = sa * sb * sg + ca * cg;  R(1, 2) = sa * sb * cg - ca * sg;
    R(2, 0) = -sb;      R(2, 1) = cb * sg;                 R(2, 2) = cb * cg;
    return R;
}

// ---------------------------------------------------------------------------
class ElasticLaw : public MatLaw {
public:
    explicit ElasticLaw(const Material& m) : MatLaw(m) {}
    Eigen::Matrix3d stress(const Eigen::Matrix3d& eps, MatState&,
                           double, double) const override {
        return elastic(eps);
    }
    std::string name() const override { return "elastic"; }
};

// ---------------------------------------------------------------------------
// MOHR-COULOMB elasto-plastique — la loi de Ye, Zhang, Chen & Li (IJRMMS 194,
// 2025, 106233) « MC-FDEM » : elements solides elasto-plastiques Mohr-Coulomb
// collaborant avec les joints cohesifs, ou toute la fissuration reste dans les
// joints et toute la dissipation PLASTIQUE dans le bulk (leur §2.2.1.1,
// eq. 1-3). Cinq parametres : E, nu (elastique) puis c, phi, psi (plastique).
//
// Surface, en convention TRACTION POSITIVE et s1 >= s2 >= s3 :
//     f = (s1 - s3) + (s1 + s3) sin(phi) - 2 c cos(phi)
// (compression uniaxiale : sc = 2 c cos/(1-sin) ; traction uniaxiale :
//  st = 2 c cos/(1+sin) — le cut-off de traction est celui du critere).
// Potentiel NON ASSOCIE (dilatance psi) : g = (s1 - s3) + (s1 + s3) sin(psi).
//
// Retour en contraintes principales facon Clausen & Krabbenhoft : plan
// principal, puis les deux ARETES (compression / extension) quand l'ordre
// s1 >= s2 >= s3 est viole, puis l'APEX quand le point depasse le sommet du
// cone. C'est le vrai Mohr-Coulomb a aretes, pas l'approximation lisse de
// Drucker-Prager (law = dpr).
// ---------------------------------------------------------------------------
class MohrCoulombLaw : public MatLaw {
public:
    MohrCoulombLaw(const Material& m, double coh, double phiDeg, double psiDeg)
        : MatLaw(m), c_(coh) {
        sphi_ = std::sin(phiDeg * M_PI / 180.0);
        cphi_ = std::cos(phiDeg * M_PI / 180.0);
        spsi_ = std::sin(psiDeg * M_PI / 180.0);
        // a^T D b du retour sur le plan principal (D isotrope) :
        //   4 [(lam + G) sin(phi) sin(psi) + G]
        denom_ = 4.0 * ((lam_ + G_) * sphi_ * spsi_ + G_);
        apex_ = (sphi_ > 1e-12) ? c_ * cphi_ / sphi_ : 0.0;   // c cot(phi)
    }

    Eigen::Matrix3d stress(const Eigen::Matrix3d& eps, MatState& s,
                           double, double) const override {
        Eigen::Matrix3d sigTr = elastic(eps - s.epsP);
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(sigTr);
        Eigen::Vector3d ev = es.eigenvalues();               // croissant
        Eigen::Matrix3d V = es.eigenvectors();
        // reordonner en s1 >= s2 >= s3 (colonnes des vecteurs suivent)
        Eigen::Vector3d p(ev(2), ev(1), ev(0));
        Eigen::Matrix3d Q;
        Q.col(0) = V.col(2); Q.col(1) = V.col(1); Q.col(2) = V.col(0);

        double f = yieldF(p);
        if (f <= 0.0) return sigTr;                          // elastique

        Eigen::Vector3d pRet = returnMap(p);
        Eigen::Matrix3d sig = Q * pRet.asDiagonal() * Q.transpose();
        // increment plastique : ce que l'elasticite ne porte plus
        s.epsP = eps - complianceApply(sig);
        return sig;
    }

    std::string name() const override { return "mc"; }

private:
    double yieldF(const Eigen::Vector3d& p) const {
        return (p(0) - p(2)) + (p(0) + p(2)) * sphi_ - 2.0 * c_ * cphi_;
    }

    // eps = D^-1 : sig pour un isotrope
    Eigen::Matrix3d complianceApply(const Eigen::Matrix3d& sig) const {
        double tr = sig.trace();
        return (sig - (lam_ / (2.0 * G_ + 3.0 * lam_)) * tr
                      * Eigen::Matrix3d::Identity()) / (2.0 * G_);
    }

    // retour sur le plan principal : sig = sigTr - dlam D b, b = dg/dsig
    Eigen::Vector3d planeReturn(const Eigen::Vector3d& p) const {
        double dlam = yieldF(p) / denom_;
        double b1 = 1.0 + spsi_, b3 = -(1.0 - spsi_);
        Eigen::Vector3d Db;                                  // D b (isotrope)
        Db(0) = (lam_ + 2.0 * G_) * b1 + lam_ * b3;
        Db(1) = lam_ * (b1 + b3);
        Db(2) = lam_ * b1 + (lam_ + 2.0 * G_) * b3;
        return p - dlam * Db;
    }

    // Retour sur une ARETE : deux surfaces actives. l = 1 : arete de
    // COMPRESSION (s1 = s2), l = 2 : arete d'EXTENSION (s2 = s3). Systeme
    // 2x2 sur les deux multiplicateurs, ecrit avec les memes gradients.
    Eigen::Vector3d edgeReturn(const Eigen::Vector3d& p, int l) const {
        Eigen::Vector3d a1(1.0 + sphi_, 0.0, -(1.0 - sphi_));
        Eigen::Vector3d b1(1.0 + spsi_, 0.0, -(1.0 - spsi_));
        Eigen::Vector3d a2, b2;
        if (l == 1) {                    // f(s2, s3) active en plus
            a2 = Eigen::Vector3d(0.0, 1.0 + sphi_, -(1.0 - sphi_));
            b2 = Eigen::Vector3d(0.0, 1.0 + spsi_, -(1.0 - spsi_));
        } else {                         // f(s1, s2) active en plus
            a2 = Eigen::Vector3d(1.0 + sphi_, -(1.0 - sphi_), 0.0);
            b2 = Eigen::Vector3d(1.0 + spsi_, -(1.0 - spsi_), 0.0);
        }
        auto D = [&](const Eigen::Vector3d& v) {
            Eigen::Vector3d r;
            double tr = v.sum();
            for (int i = 0; i < 3; ++i) r(i) = lam_ * tr + 2.0 * G_ * v(i);
            return r;
        };
        Eigen::Vector3d Db1 = D(b1), Db2 = D(b2);
        Eigen::Matrix2d A;
        A << a1.dot(Db1), a1.dot(Db2), a2.dot(Db1), a2.dot(Db2);
        Eigen::Vector2d rhs(a1.dot(p) - 2.0 * c_ * cphi_,
                            a2.dot(p) - 2.0 * c_ * cphi_);
        Eigen::Vector2d dl = A.fullPivLu().solve(rhs);
        dl(0) = std::max(dl(0), 0.0);
        dl(1) = std::max(dl(1), 0.0);
        return p - dl(0) * Db1 - dl(1) * Db2;
    }

    Eigen::Vector3d returnMap(const Eigen::Vector3d& p) const {
        Eigen::Vector3d r = planeReturn(p);
        if (r(0) >= r(1) && r(1) >= r(2)) {                  // ordre conserve
            return apexClamp(r);
        }
        // l'ordre est viole : le point appartient a une arete
        int l = (r(0) < r(1)) ? 1 : 2;
        Eigen::Vector3d re = edgeReturn(p, l);
        // tri de securite : le retour d'arete rend deux valeurs egales
        double v[3] = {re(0), re(1), re(2)};
        std::sort(v, v + 3, std::greater<double>());
        return apexClamp(Eigen::Vector3d(v[0], v[1], v[2]));
    }

    // APEX : si le retour a franchi le sommet du cone (traction hydrostatique
    // au-dela de c cot(phi)), le point s'y projette — sinon le materiau
    // « tirerait » plus que sa cohesion ne le permet.
    Eigen::Vector3d apexClamp(const Eigen::Vector3d& r) const {
        if (sphi_ <= 1e-12) return r;
        double m = r.sum() / 3.0;
        if (m <= apex_) return r;
        return Eigen::Vector3d::Constant(apex_);
    }

    double c_, sphi_, cphi_, spsi_, denom_, apex_;
};

// ---------------------------------------------------------------------------
// Shared plastic-damage kernel: DP cone return (viscous if eta > 0),
// optional pressure cap, Rankine crack-band tensile damage. dpr and saksala
// are two parameterizations of this kernel.
// ---------------------------------------------------------------------------
// Briques opt-in du noyau (voir MatLaw.hpp) ; tout a zero/false = le noyau
// historique, bit-identique.
struct BrickOpts {
    bool power = false;              // meridian = power
    double fc0 = 0.0;                // UCS du cone (sigmaCdp), continuite en s3 = 0
    double B = 56.59, n = 0.538;     // q = fc0 + 1e6 B (s3/1e6)^n  [MPa^(1-n)]
    bool compDam = false;            // compDamage = crackband
    double Ac = 0.98, GIIc = 1.0e4;  // omega_c = Ac (1 - exp(-b_c epvEq)), b_c = fc0 h/GIIc
    double erodeDc = 0.0;            // erosion sur omega_c / Ac (0 = off)
    // spall sur l'ENERGIE dissipee : wDamT >= erodeWfrac * Gf / lc (0 = off).
    // Le seuil erodeD = 0,98 sur D supprime l'element a kappa ~ 50 k0, ou
    // il porte encore 78 % de ft (D = 1 - k0/kappa mesure la perte de
    // raideur, pas la chute de contrainte) : 80 % de Gf jamais dissipes
    // (barre E1, 2026-09-04). Ce critere-ci retire l'element quand la bande
    // a rendu la fraction voulue de son energie de rupture.
    double erodeWfrac = 0.0;
    // dpApex : au-dela de l'apex du cone (k_eff = k - 3 alpha p <= 0, traction
    // hydrostatique > k/(3 alpha) = 18,8 MPa sur la carte Bohus) le retour
    // radial deviatorique a p fixe donne dlam > sqrtJ2/G : le deviateur CHANGE
    // DE SIGNE, et la part compressive ainsi fabriquee echappe au split
    // unilateral — injection d'energie (mini-bloc B3 : 1e155 J ; A2b : travail
    // plastique -2e4 J). Le seuil erodeD = 0,98 masquait le defaut en
    // supprimant tot les elements tendus. Avec dpApex, aucun retour DP
    // au-dela de l'apex : la traction est l'affaire du cut-off de Rankine
    // (l'OPTION-3 « DP en compression seule » de la VUMAT DP-DFH).
    bool apex = false;
    // dpTension = off : le cone DP n'agit qu'en COMPRESSION (p <= 0) — la
    // vraie OPTION-3 de la VUMAT DP-DFH (pbar > 0). Entre p = 0 et l'apex, la
    // jambe tractive du cone (23,1 MPa effectifs en uniaxial) coulait dans
    // l'espace EFFECTIF pendant que Rankine endommageait : sur la barre E1
    // cette « plasticite tractive » dissipait 3,3 Gf A (revue du 2026-09-04).
    // Avec dpTension = off la traction est l'affaire du cut-off seul.
    bool tensionOff = false;
    // rankineDrive = stress : le cut-off de Rankine est pilote par la
    // CONTRAINTE effective principale maximale (kappa = sigma1_eff / E) et non
    // par la deformation principale maximale (noyau d'origine — malgre le
    // commentaire de l'en-tete). En deformation, la dilatation de Poisson
    // d'un element COMPRIME (2 nu P / E = 7,5e-4 sous 100 MPa lateraux, ou
    // nu sigma / E = 4,7e-4 sous l'UCS) depasse k0 = ft / E = 1,2e-4 et met
    // D = 1 - k0/kappa ~ 0,85 partout, sans aucune traction (matrice C_T1_R_P100 :
    // 100 % du bloc a D >= 0,5 avant l'impact). Les deux pilotages coincident
    // en traction uniaxiale (E1 identique).
    bool rankineStress = false;
    // tensionDamage = fixed : endommagement de traction a direction figee
    // (fixed crack model, Rashid 1968 ; Rots & Blaauwendraad 1989) — voir
    // MatLaw.hpp et fixedCrackUpdate / nominal ci-dessous. shearRet = beta
    // (tensionShearRetention, 1 = pas de retention).
    bool fixedCrack = false;
    double shearRet = 1.0;
};

class PlasticDamageLaw : public MatLaw {
public:
    PlasticDamageLaw(const Material& m, double eta, double capP0, double capH,
                     double erodeD, double erodeEpv,
                     const BrickOpts& br = BrickOpts())
        : MatLaw(m), eta_(eta), capP0_(capP0), capH_(capH), erodeD_(erodeD),
          erodeEpv_(erodeEpv), br_(br) {}

    // q de la surface puissance a p (trial) fixe : q = fc0 + B s3^n avec la
    // pseudo-sigma3 s3 = -p - q/3 (exacte en compression triaxiale). Rend < 0
    // si la racine tombe cote tractif (s3 < 0) : le cone lineaire s'applique.
    // g(q) = q - fc0 - B s3(q)^n est croissante et convexe : Newton depuis
    // q0 >= racine converge de facon monotone.
    double yieldQPower(double p) const {
        double s3max = -p - br_.fc0 / 3.0;              // s3 a q = fc0
        if (!(s3max > 0.0)) return -1.0;
        // On resout en s3 (pas en q) : h(s3) = fc0 + B s3^n + 3 p + 3 s3 = 0,
        // croissante et concave sur [0, s3max], h(0) = -3 s3max < 0,
        // h(s3max) = B s3max^n > 0. Newton depuis s3max (un pas a gauche de la
        // racine, puis convergence monotone par en dessous), borne dans le
        // crochet avec repli bissection — la pente infinie de s3^n en 0
        // (n < 1) faisait caler le Newton en q.
        double lo = 0.0, hi = s3max, s3 = s3max;
        for (int it = 0; it < 60; ++it) {
            double h = br_.fc0 + Bpow(s3) + 3.0 * p + 3.0 * s3;
            if (h > 0.0) hi = s3; else lo = s3;
            if (std::abs(h) < 1.0e-9 * br_.fc0) break;
            double dh = dBpow(s3) + 3.0;
            double sn = s3 - h / dh;
            if (!(sn > lo && sn < hi)) sn = 0.5 * (lo + hi);   // bissection
            if (std::abs(sn - s3) < 1.0e-12 * std::max(s3max, 1.0)) { s3 = sn; break; }
            s3 = sn;
        }
        return -3.0 * p - 3.0 * s3;                     // q = 3 (-p - s3)
    }

    Eigen::Matrix3d stress(const Eigen::Matrix3d& eps, MatState& s,
                           double dt, double lc) const override {
        if (s.eroded) return Eigen::Matrix3d::Zero();
        if (capP0_ > 0.0 && s.pc == 0.0) s.pc = capP0_;   // lazy init

        // ---- elastic predictor on the effective (undamaged) skeleton ----
        Eigen::Matrix3d sig = elastic(eps - s.epsP);

        // ---- pressure cap (compression crush) ---------------------------
        double p = sig.trace() / 3.0;
        double capDv = 0.0;                            // compaction (compteur)
        if (capP0_ > 0.0 && p < -s.pc) {
            double over = -(p + s.pc);                 // > 0
            double dev_ = over / (K_ + capH_);         // |volumetric return|
            s.epsP -= (dev_ / 3.0) * Eigen::Matrix3d::Identity();
            s.pc += capH_ * dev_;
            p += K_ * dev_;
            sig += K_ * dev_ * Eigen::Matrix3d::Identity();
            capDv = dev_;
        }

        // ---- DP cone, deviatoric (visco)plastic return ------------------
        Eigen::Matrix3d dev = sig - p * Eigen::Matrix3d::Identity();
        double sj2 = std::sqrt(0.5) * dev.norm();      // sqrt(J2)
        double F = sj2 + adp_ * 3.0 * p - kdp_;        // I1 = 3 p
        if (br_.power) {
            // meridien puissance : sqrt(J2) <= q_y(p) / sqrt(3), meme retour
            double qy = yieldQPower(p);
            if (qy > 0.0) F = sj2 - qy / std::sqrt(3.0);
        }
        // dpApex : pas de retour au-dela de l'apex (k_eff <= 0) ; dpTension =
        // off : pas de retour en traction (p > 0) — voir BrickOpts
        if ((br_.apex && kdp_ - adp_ * 3.0 * p <= 0.0)
            || (br_.tensionOff && p > 0.0)) F = -1.0;
        Eigen::Matrix3d dEpsP = Eigen::Matrix3d::Zero();   // increment (compteur)
        if (capDv > 0.0) dEpsP -= (capDv / 3.0) * Eigen::Matrix3d::Identity();
        if (F > 0.0 && sj2 > 1e-12) {
            // linear Perzyna: F_{n+1} = eta dlam/dt with radial deviatoric
            // return sqrtJ2 -> sqrtJ2 - G dlam  =>  closed form
            double dlam = F / (G_ + eta_ / std::max(dt, 1e-30));
            Eigen::Matrix3d nfl = dev / (2.0 * sj2);   // flow direction
            s.epsP += dlam * nfl;
            dev *= (1.0 - G_ * dlam / sj2);
            sig = dev + p * Eigen::Matrix3d::Identity();
            s.epvEq += dlam / std::sqrt(3.0);
            dEpsP += dlam * nfl;
        }

        // ---- Rankine crack-band tensile damage --------------------------
        // E1 (2026-08-19) : ftScale est desormais HONORE. Avant ce correctif,
        // le facteur de Weibull par element etait tire, ECRIT DANS LES VTU,
        // et sans le moindre effet sous les lois dpr et saksala — seules
        // dpdfh et saksala2011 le lisaient. On croyait donc avoir une
        // heterogeneite, on n'en avait pas. Audit du 19/08 : aucune config du
        // depot ne combine dpr/saksala et matWeibullM, le correctif ne change
        // donc aucun resultat existant.
        // Gf N'EST PAS mis a l'echelle : le facteur de Weibull porte sur la
        // RESISTANCE (mecanisme FIELD des VUMAT), pas sur la tenacite. La
        // longueur cohesive locale E Gf / ft^2 varie donc en 1/ftScale^2, ce
        // qui est le comportement attendu d'un materiau dont seuls les defauts
        // sont distribues.
        const double ftLoc = mat_.ft * s.ftScale;
        const double GfLoc = wScaleGf_ ? mat_.Gf * s.ftScale : mat_.Gf;
        double k0 = ftLoc / mat_.E;
        double kf = GfLoc / (lc * ftLoc) - 0.5 * k0;
        if (kf <= 0.05 * k0)
            throw std::runtime_error("MatLaw: element size " + std::to_string(lc)
                + " m exceeds the crack-band limit E Gf / ft^2 — refine the "
                  "mesh or raise Gf");
        const double Dold = s.D;
        double dOldF[3] = {0.0, 0.0, 0.0};             // tensionDamage = fixed
        if (!br_.fixedCrack) {
        double e1 = maxPrincipal(eps - s.epsP);        // driving strain
        if (br_.rankineStress)                         // pilotage en contrainte
            e1 = maxPrincipal(sig) / mat_.E;           // sigma1_eff / E
        if (e1 > k0 && e1 > s.kappa) s.kappa = e1;
        if (s.kappa > k0) {
            double D = 1.0 - (k0 / s.kappa)
                             * std::exp(-(s.kappa - k0) / kf);
            if (D > s.D) s.D = std::min(1.0, D);
        }
        } else {
            // direction figee : amorcage / croissance par direction, meme
            // cinetique (k0, kf) ; s.D = max d_i, s.kappa = max kappa_i
            fixedCrackUpdate(eps, sig, s, k0, kf, dOldF);
        }

        // ---- endommagement COMPRESSIF crack-band (compDamage = crackband) --
        // omega_c = Ac (1 - exp(-b_c epvEq)), b_c = fc0 h_e / G_IIc, monotone ;
        // porte sur la partie spectrale negative (brique omega_c de MH 2018)
        const double DcOld = s.Dc;
        if (br_.compDam) {
            double bc = br_.fc0 * lc / br_.GIIc;
            double wc = br_.Ac * (1.0 - std::exp(-bc * s.epvEq));
            if (wc > s.Dc) s.Dc = wc;
        }

        // UNILATERAL damage: degrade the TENSILE principal components only.
        // A scalar (1-D) on the full stress also kills the COMPRESSIVE
        // bearing capacity: under an indenter the damaged surface elements
        // then collapse and the tool tunnels through the block at constant
        // velocity (measured before the fix). Cracked rock still bears
        // compression — the same reason percussion VUMATs keep damaged
        // elements alive instead of deleting them.
        Eigen::Matrix3d nominal = sig;
        if (br_.fixedCrack) {
            // ---- tensionDamage = fixed : operateur d'endommagement dans le
            // repere FIGE (jamais atteint a cle absente). S = R^T sigma R ;
            // S_ii (1 - d_i) si S_ii > 0 (fissure ouverte ; la compression
            // normale passe intacte : unilateral), S_ij (1 - beta max(d_i^ouv,
            // d_j^ouv)) ; retour au repere co-rotationnel. omega_c (compDamage)
            // reste sur la partie spectrale negative de sigma_eff, soustraite
            // apres l'operateur fige (sans d_i : (1 - Dc) sigC + sigT).
            if (s.fcm.nAct > 0) {
                const Eigen::Matrix3d R = fcmFrame(s.fcm.eul);
                const Eigen::Matrix3d S = R.transpose() * sig * R;
                double dop[3];
                for (int i = 0; i < 3; ++i)
                    dop[i] = S(i, i) > 0.0 ? s.fcm.d[i] : 0.0;
                Eigen::Matrix3d Sn = S;
                for (int i = 0; i < 3; ++i) {
                    Sn(i, i) = S(i, i) * (1.0 - dop[i]);
                    for (int j = i + 1; j < 3; ++j) {
                        double fac = 1.0 - br_.shearRet * std::max(dop[i], dop[j]);
                        Sn(i, j) = S(i, j) * fac;
                        Sn(j, i) = Sn(i, j);
                    }
                }
                nominal = R * Sn * R.transpose();
                // compteur : Y_i = 1/2 <S_ii>^2 / E (energie normale effective
                // de la direction i ; le cisaillement retenu n'est pas compte)
                for (int i = 0; i < 3; ++i)
                    if (s.fcm.d[i] > dOldF[i] && S(i, i) > 0.0)
                        s.wDamT += 0.5 * S(i, i) * S(i, i) / mat_.E
                                   * (s.fcm.d[i] - dOldF[i]);
            }
            if (s.Dc > 0.0) {
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(sig);
                Eigen::Matrix3d sigC = Eigen::Matrix3d::Zero();
                for (int q = 0; q < 3; ++q) {
                    double lq = es.eigenvalues()(q);
                    if (lq < 0.0)
                        sigC += lq * es.eigenvectors().col(q)
                                   * es.eigenvectors().col(q).transpose();
                }
                nominal -= s.Dc * sigC;
                if (s.Dc > DcOld) s.wDamC += damageForce(sigC) * (s.Dc - DcOld);
            }
        } else if (s.Dc > 0.0) {
            // branche compDamage (jamais atteinte a cles absentes : Dc = 0)
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(sig);
            Eigen::Matrix3d sigT = Eigen::Matrix3d::Zero();
            for (int q = 0; q < 3; ++q) {
                double lq = es.eigenvalues()(q);
                if (lq > 0.0)
                    sigT += lq * es.eigenvectors().col(q)
                               * es.eigenvectors().col(q).transpose();
            }
            Eigen::Matrix3d sigC = sig - sigT;
            nominal = (1.0 - s.Dc) * sigC + (1.0 - s.D) * sigT;
            if (s.D > Dold)   s.wDamT += damageForce(sigT) * (s.D - Dold);
            if (s.Dc > DcOld) s.wDamC += damageForce(sigC) * (s.Dc - DcOld);
        } else if (s.D > 0.0) {
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(sig);
            Eigen::Matrix3d sigT = Eigen::Matrix3d::Zero();
            for (int q = 0; q < 3; ++q) {
                double lq = es.eigenvalues()(q);
                if (lq > 0.0)
                    sigT += lq * es.eigenvectors().col(q)
                               * es.eigenvectors().col(q).transpose();
            }
            nominal = (sig - sigT) + (1.0 - s.D) * sigT;
            if (s.D > Dold) s.wDamT += damageForce(sigT) * (s.D - Dold);
        }
        // travail plastique en contrainte NOMINALE (compteur, sans effet) :
        // sigma_nom : d eps_p — la pondération (1 - omega_c) seule
        // sur-comptait les etats mixtes (revue du 2026-09-04)
        s.wPlas += nominal.cwiseProduct(dEpsP).sum();

        // erosion = material REMOVAL, only where it is physical: fully
        // damaged elements in net TENSION (spalled chips leaving the
        // surface) or over-crushed elements (comminuted material squeezed
        // out, erodeEpv on the equivalent viscoplastic strain). A damaged
        // element under compression stays: it is the rubble bed.
        // Troisieme canal opt-in (erodeDc) : omega_c / Ac normalise.
        bool spall = s.D >= erodeD_ && nominal.trace() >= 0.0;
        // spall sur l'energie de la bande (erodeWfrac, opt-in)
        if (!spall && br_.erodeWfrac > 0.0 && nominal.trace() >= 0.0
            && s.wDamT >= br_.erodeWfrac * GfLoc / lc)
            spall = true;
        bool crush = erodeEpv_ > 0.0 && s.epvEq >= erodeEpv_;
        bool crushDc = br_.erodeDc > 0.0 && br_.compDam
                       && s.Dc / br_.Ac >= br_.erodeDc;
        if (spall || crush || crushDc) {
            s.eroded = true;
            s.eroCode = spall ? 1 : (crush ? 2 : 3);
            return Eigen::Matrix3d::Zero();
        }
        return nominal;
    }

    std::string name() const override {
        return eta_ > 0.0 ? "saksala" : "dpr";
    }

    // dsig = sqrt(3) eta epdot / (1/sqrt(3) - alpha): lambda-dot = sqrt(3)
    // epdot at steady uniaxial flow, F = eta lambda-dot, dF/dsigma = 1/sqrt3 - a
    double viscousOverstress(double epdot) const override {
        return std::sqrt(3.0) * eta_ * epdot / (1.0 / std::sqrt(3.0) - adp_);
    }

    const BrickOpts& bricks() const { return br_; }

private:
    // loi puissance en unites MPa : 1e6 B (s3/1e6)^n et sa derivee en s3
    double Bpow(double s3) const {
        return 1.0e6 * br_.B * std::pow(std::max(s3, 0.0) / 1.0e6, br_.n);
    }
    double dBpow(double s3) const {
        double x = std::max(s3, 1.0) / 1.0e6;          // garde s3 -> 0 (n < 1)
        return br_.B * br_.n * std::pow(x, br_.n - 1.0);
    }
    // force thermodynamique d'endommagement Y = 1/2 sigma : C^-1 : sigma
    // (isotrope) pour une partie spectrale de la contrainte effective
    double damageForce(const Eigen::Matrix3d& sg) const {
        double tr = sg.trace();
        return 0.5 * ((1.0 + mat_.nu) / mat_.E * sg.squaredNorm()
                      - mat_.nu / mat_.E * tr * tr);
    }

    // ---- tensionDamage = fixed : amorcage et croissance par direction -----
    // Tenseur pilote T = meme grandeur que le cut-off scalaire (sigma_eff/E
    // en rankineDrive = stress, eps - eps_p en deformation). Amorcage : la
    // plus grande valeur propre de T restreinte au COMPLEMENT ORTHOGONAL des
    // directions deja figees depasse k0 = ft/E -> la direction est figee
    // (1re : repere principal complet, n1 = vecteur propre max, main droite ;
    // 2e : rotation de (n2, n3) autour de n1 vers le vecteur propre max du
    // 2x2 restreint ; 3e : n3 seul). Croissance : kappa_i = max(kappa_i,
    // n_i.T.n_i), d_i = 1 - k0/kappa_i exp(-(kappa_i - k0)/kf) monotone,
    // identique a la formule scalaire (u_i = kappa_i lc). dOld recoit les
    // d_i d'entree (compteur wDamT).
    void fixedCrackUpdate(const Eigen::Matrix3d& eps, const Eigen::Matrix3d& sig,
                          MatState& s, double k0, double kf, double dOld[3]) const {
        auto& f = s.fcm;
        for (int i = 0; i < 3; ++i) dOld[i] = f.d[i];
        const Eigen::Matrix3d T = br_.rankineStress ? Eigen::Matrix3d(sig / mat_.E)
                                                    : Eigen::Matrix3d(eps - s.epsP);
        Eigen::Matrix3d R = f.nAct > 0 ? fcmFrame(f.eul) : Eigen::Matrix3d::Identity();
        if (f.nAct == 0) {
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(T);   // croissant
            const double e = es.eigenvalues()(2);
            if (e > k0) {
                R.col(0) = es.eigenvectors().col(2);
                R.col(1) = es.eigenvectors().col(1);
                R.col(2) = es.eigenvectors().col(0);
                if (R.determinant() < 0.0) R.col(2) = -R.col(2);
                fcmEuler(R, f.eul);
                R = fcmFrame(f.eul);           // le repere STOCKE fait foi
                f.nAct = 1;
                f.kap[0] = e;
            }
        } else if (f.nAct == 1) {
            const Eigen::Vector3d n2 = R.col(1), n3 = R.col(2);
            const double a = n2.dot(T * n2), b = n3.dot(T * n3), c = n2.dot(T * n3);
            const double e = 0.5 * (a + b) + std::sqrt(0.25 * (a - b) * (a - b) + c * c);
            if (e > k0) {
                const double th = 0.5 * std::atan2(2.0 * c, a - b);
                const double ct = std::cos(th), st = std::sin(th);
                R.col(1) = ct * n2 + st * n3;
                R.col(2) = -st * n2 + ct * n3;
                fcmEuler(R, f.eul);
                R = fcmFrame(f.eul);
                f.nAct = 2;
                f.kap[1] = e;
            }
        } else if (f.nAct == 2) {
            const Eigen::Vector3d n3 = R.col(2);
            const double e = n3.dot(T * n3);
            if (e > k0) { f.nAct = 3; f.kap[2] = e; }
        }
        for (int i = 0; i < f.nAct; ++i) {
            const Eigen::Vector3d ni = R.col(i);
            const double ei = ni.dot(T * ni);
            if (ei > k0 && ei > f.kap[i]) f.kap[i] = ei;
            if (f.kap[i] > k0) {
                double d = 1.0 - (k0 / f.kap[i]) * std::exp(-(f.kap[i] - k0) / kf);
                if (d > f.d[i]) f.d[i] = std::min(1.0, d);
            }
        }
        s.D = std::max({f.d[0], f.d[1], f.d[2]});
        s.kappa = std::max({f.kap[0], f.kap[1], f.kap[2]});
    }

    double eta_, capP0_, capH_, erodeD_, erodeEpv_;
    BrickOpts br_;
};

} // namespace

// ===========================================================================
// saksala2011 — line-by-line port of the thesis' vumat_saksala_2011.f90
// (module saksala_2011_model). Voigt-6 order 11 22 33 12 23 13, tensorial
// shears, tension positive. Names and structure kept close to the Fortran
// for auditability; the principal decomposition uses Eigen instead of the
// hand-rolled Jacobi (equivalent to fp accuracy, checked by the trace
// superposition test `rockim selftest-saksala2011`).
// ===========================================================================
namespace sk11 {

using V6 = std::array<double, 6>;
constexpr double tiny_num = 1.0e-14;

inline void invariants(const V6& s, double& i1, double& q) {
    i1 = s[0] + s[1] + s[2];
    double p = i1 / 3.0;
    double d1 = s[0] - p, d2 = s[1] - p, d3 = s[2] - p;
    double j2 = 0.5 * (d1 * d1 + d2 * d2 + d3 * d3)
                + s[3] * s[3] + s[4] * s[4] + s[5] * s[5];
    q = std::sqrt(std::max(j2, 0.0));
}

inline double tensorDot(const V6& a, const V6& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
           + 2.0 * (a[3] * b[3] + a[4] * b[4] + a[5] * b[5]);
}

inline void principal(const V6& s, Eigen::Vector3d& eig, Eigen::Matrix3d& vec) {
    Eigen::Matrix3d A;
    A << s[0], s[3], s[5],
         s[3], s[1], s[4],
         s[5], s[4], s[2];
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(A);
    eig = es.eigenvalues();
    vec = es.eigenvectors();
}

inline void rebuild(const Eigen::Vector3d& eig, const Eigen::Matrix3d& vec,
                    V6& s) {
    Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
    for (int k = 0; k < 3; ++k)
        A += eig(k) * vec.col(k) * vec.col(k).transpose();
    s = {A(0, 0), A(1, 1), A(2, 2), A(0, 1), A(1, 2), A(0, 2)};
}

inline void elasticApply(const V6& e, double bulk, double shear, V6& s) {
    double lame = bulk - 2.0 * shear / 3.0;
    double tr = e[0] + e[1] + e[2];
    for (int i = 0; i < 3; ++i) s[i] = lame * tr + 2.0 * shear * e[i];
    for (int i = 3; i < 6; ++i) s[i] = 2.0 * shear * e[i];
}

inline void dpGradients(const V6& s, double alpha, double beta,
                        V6& nf, V6& ng) {
    double i1, q;
    invariants(s, i1, q);
    double p = i1 / 3.0;
    V6 dev = s;
    for (int i = 0; i < 3; ++i) dev[i] -= p;
    double qtol = 1.0e-12 * std::max(1.0, std::abs(i1));
    if (q > qtol)
        for (int i = 0; i < 6; ++i) { nf[i] = dev[i] / (2.0 * q); ng[i] = nf[i]; }
    else
        for (int i = 0; i < 6; ++i) { nf[i] = 0.0; ng[i] = 0.0; }
    for (int i = 0; i < 3; ++i) { nf[i] += alpha; ng[i] += beta; }
}

inline void mrGradient(const V6& s, V6& nmr, double& qpos) {
    Eigen::Vector3d eig;
    Eigen::Matrix3d vec;
    principal(s, eig, vec);
    qpos = 0.0;
    for (int i = 0; i < 3; ++i)
        qpos += std::max(eig(i), 0.0) * std::max(eig(i), 0.0);
    qpos = std::sqrt(qpos);
    Eigen::Vector3d nv = Eigen::Vector3d::Zero();
    if (qpos > tiny_num)
        for (int i = 0; i < 3; ++i) nv(i) = std::max(eig(i), 0.0) / qpos;
    rebuild(nv, vec, nmr);
}

inline void capGradient(const V6& s, double c1, double c2, V6& ncap) {
    double i1, q;
    invariants(s, i1, q);
    double p = i1 / 3.0;
    V6 dev = s;
    for (int i = 0; i < 3; ++i) dev[i] -= p;
    double qtol = 1.0e-12 * std::max(1.0, std::abs(i1));
    if (q > qtol)
        for (int i = 0; i < 6; ++i) ncap[i] = dev[i] / (2.0 * q);
    else
        for (int i = 0; i < 6; ++i) ncap[i] = 0.0;
    double hydro = -(2.0 * c1 * i1 + c2);
    for (int i = 0; i < 3; ++i) ncap[i] += hydro;
}

inline void capCoefficients(double alpha, double beta, double intercept,
                            double ptr0, double pp0, double pp,
                            double& c1, double& c2, double& c3, double& ptr) {
    ptr = ptr0 * pp / pp0;                 // transition tracks the closure
    double itr = -ptr, ic = -pp;
    double delta = itr - ic;
    if (std::abs(delta) <= tiny_num) {
        c1 = 0.0;
        c2 = -beta;
        c3 = beta * ic;
        return;
    }
    c1 = (alpha * itr - intercept + beta * (ic - itr)) / (delta * delta);
    c2 = -beta - 2.0 * c1 * itr;
    c3 = -c1 * ic * ic - c2 * ic;
}

inline void positiveTensorNorm(const V6& t, double& value) {
    Eigen::Vector3d eig;
    Eigen::Matrix3d vec;
    principal(t, eig, vec);
    value = 0.0;
    for (int i = 0; i < 3; ++i)
        value += std::max(eig(i), 0.0) * std::max(eig(i), 0.0);
    value = std::sqrt(value);
}

inline void nominalStress(const V6& eff, double omega, V6& nom) {
    Eigen::Vector3d eig;
    Eigen::Matrix3d vec;
    principal(eff, eig, vec);
    for (int i = 0; i < 3; ++i)
        eig(i) = (1.0 - omega) * std::max(eig(i), 0.0) + std::min(eig(i), 0.0);
    rebuild(eig, vec, nom);
}

struct Props {
    double young, nu, phiDeg, betaDP, ft0, c0, cres0;
    double hdp0, sdp, smr, at, betaT, pp0, ptr0, dcap, wcap, nd;
};

// the full update: transcription of saksala_update_point
inline void updatePoint(double dt, const V6& deps, const Props& P,
                        MatState::Sk11& st, V6& stressNew) {
    double phi = P.phiDeg * M_PI / 180.0;
    double shear = P.young / (2.0 * (1.0 + P.nu));
    double bulk = P.young / (3.0 * (1.0 - 2.0 * P.nu));
    double lame = bulk - 2.0 * shear / 3.0;
    double alpha_dp = 2.0 * std::sin(phi) / (std::sqrt(3.0) * (3.0 - std::sin(phi)));
    double k_dp = 6.0 * std::cos(phi) / (std::sqrt(3.0) * (3.0 - std::sin(phi)));
    double kc = std::sqrt((1.0 + 6.0 * P.betaDP * P.betaDP) / 3.0);

    // local strengths from the state (SDV15/16 in the VUMAT: the FIELD
    // heterogeneity mechanism), falling back to the homogeneous props
    double ft0 = st.ftLoc > 0.0 ? st.ftLoc : P.ft0;
    double c0 = st.c0Loc > 0.0 ? st.c0Loc : P.c0;
    // cres scales with the local cohesion, as in the VUMAT's
    // cres0 = cres_ref * c0 / c_ref
    double cres0 = P.cres0 * c0 / std::max(P.c0, 1e-30);
    double kap_dp_old = std::max(st.kapDP, 0.0);
    double kap_mr_old = std::max(st.kapMR, 0.0);
    double eqvt_old = std::max(st.eqvt, 0.0);
    double epsv = std::max(st.epsv, 0.0);
    double pp = std::max(st.pp, P.pp0);

    V6 sbar, depvp = {0, 0, 0, 0, 0, 0};
    {
        double tr = deps[0] + deps[1] + deps[2];
        for (int i = 0; i < 3; ++i)
            sbar[i] = st.sbar[i] + lame * tr + 2.0 * shear * deps[i];
        for (int i = 3; i < 6; ++i)
            sbar[i] = st.sbar[i] + 2.0 * shear * deps[i];
    }

    // Eq. (17): confinement from the elastic trial state, frozen during
    // the local return mapping
    double sigma_conf = 0.0;
    {
        Eigen::Vector3d eig;
        Eigen::Matrix3d vec;
        principal(sbar, eig, vec);
        std::array<double, 3> e = {eig(0), eig(1), eig(2)};
        std::sort(e.begin(), e.end(), std::greater<double>());
        if (e[0] < 0.0 && e[1] < 0.0 && e[2] < 0.0)
            sigma_conf = -0.5 * (e[0] + e[1]);
    }
    double scale = P.nd > 0.0 ? std::exp(-P.nd * sigma_conf) : 1.0;
    double hdp = P.hdp0 * scale;
    double cres = std::max(cres0, (1.0 - scale) * c0);

    double lam_dp = 0.0, lam_mr = 0.0;
    double tol = 1.0e-9 * std::max({1.0, std::abs(c0), std::abs(ft0),
                                    std::abs(P.pp0)});
    int active = 0, niter = 0, failed = 0;

    auto strengths = [&](double lamDP, double lamMR, double& c_dyn,
                         double& c_static, double& ft_static, double& ft_dyn) {
        c_static = std::max(cres, c0 + hdp * (kap_dp_old + kc * lamDP));
        c_dyn = std::max(cres, c0 + hdp * (kap_dp_old + kc * lamDP)
                               + P.sdp * kc * lamDP / std::max(dt, tiny_num));
        ft_static = ft0 * c_static / std::max(c0, tiny_num);
        ft_dyn = ft_static + P.smr * lamMR / std::max(dt, tiny_num);
    };

    for (int outer = 1; outer <= 30; ++outer) {
        bool did_work = false;
        double c_dyn, c_static, ft_static, ft_dyn;
        strengths(lam_dp, lam_mr, c_dyn, c_static, ft_static, ft_dyn);

        double i1, q, c1, c2, c3, ptr;
        invariants(sbar, i1, q);
        capCoefficients(alpha_dp, P.betaDP, k_dp * c_dyn, P.ptr0, P.pp0, pp,
                        c1, c2, c3, ptr);
        double fcap = q - (c1 * i1 * i1 + c2 * i1 + c3);
        double fn = alpha_dp > tiny_num
                        ? q - i1 / alpha_dp
                              - (ptr * (alpha_dp + 1.0 / alpha_dp)
                                 + k_dp * c_dyn)
                        : -1.0;

        if (fn > 0.0) {
            if (fcap > tol) {
                // ---- return_cap (generalized cutting plane) -------------
                int it;
                for (it = 1; it <= 80; ++it) {
                    pp = P.pp0 + std::log(1.0 + epsv / P.wcap) / P.dcap;
                    capCoefficients(alpha_dp, P.betaDP, k_dp * c_dyn, P.ptr0,
                                    P.pp0, pp, c1, c2, c3, ptr);
                    invariants(sbar, i1, q);
                    double f = q - (c1 * i1 * i1 + c2 * i1 + c3);
                    if (f <= 10.0 * tol) break;
                    V6 ncap, cncap;
                    capGradient(sbar, c1, c2, ncap);
                    elasticApply(ncap, bulk, shear, cncap);
                    double trn = ncap[0] + ncap[1] + ncap[2];
                    double hprobe = std::max(
                        1.0e-9, 1.0e-6 * std::max(1.0, std::abs(f))
                                    / std::max(tensorDot(ncap, cncap), 1.0));
                    V6 sp;
                    for (int i = 0; i < 6; ++i)
                        sp[i] = sbar[i] - hprobe * cncap[i];
                    double epsv_p = std::max(epsv - hprobe * trn, 0.0);
                    double pp_p = P.pp0
                                  + std::log(1.0 + epsv_p / P.wcap) / P.dcap;
                    double c1p, c2p, c3p, ptrp, i1p, qp;
                    capCoefficients(alpha_dp, P.betaDP, k_dp * c_dyn, P.ptr0,
                                    P.pp0, pp_p, c1p, c2p, c3p, ptrp);
                    invariants(sp, i1p, qp);
                    double fp = qp - (c1p * i1p * i1p + c2p * i1p + c3p);
                    double denom = -(fp - f) / hprobe;
                    if (denom <= tiny_num) {
                        if (f > 100.0 * tol) failed = 1;
                        break;
                    }
                    double dl = f / denom;
                    if (dl <= 0.0) { failed = 1; break; }
                    for (int i = 0; i < 6; ++i) {
                        sbar[i] -= dl * cncap[i];
                        depvp[i] += dl * ncap[i];
                    }
                    epsv = std::max(epsv - dl * trn, 0.0);
                }
                pp = P.pp0 + std::log(1.0 + epsv / P.wcap) / P.dcap;
                niter += std::min(it, 80);
                if (it > 80) failed = 1;
                active = 4;
                did_work = true;
            }
        } else {
            // yield_values
            double fdp = q + alpha_dp * i1 - k_dp * c_dyn;
            double qpos;
            V6 nmr;
            mrGradient(sbar, nmr, qpos);
            double fmr = qpos - ft_dyn;

            if (fdp > tol && fmr > tol) {
                // ---- return_corner (Koiter, coupled Newton) -------------
                int it;
                for (it = 1; it <= 80; ++it) {
                    double cdraw = c0 + hdp * (kap_dp_old + kc * lam_dp)
                                   + P.sdp * kc * lam_dp / std::max(dt, tiny_num);
                    double cdyn = std::max(cres, cdraw);
                    double cstat = std::max(
                        cres, c0 + hdp * (kap_dp_old + kc * lam_dp));
                    double ftstat = ft0 * cstat / std::max(c0, tiny_num);
                    double ftdyn = ftstat
                                   + P.smr * lam_mr / std::max(dt, tiny_num);
                    invariants(sbar, i1, q);
                    mrGradient(sbar, nmr, qpos);
                    double f1 = q + alpha_dp * i1 - k_dp * cdyn;
                    double f2 = qpos - ftdyn;
                    if (std::max(f1, f2) <= tol) break;
                    V6 nf, ng, cng, cnmr;
                    dpGradients(sbar, alpha_dp, P.betaDP, nf, ng);
                    elasticApply(ng, bulk, shear, cng);
                    elasticApply(nmr, bulk, shear, cnmr);
                    double dh = cdraw > cres
                                    ? hdp + P.sdp / std::max(dt, tiny_num)
                                    : 0.0;
                    double a11 = tensorDot(nf, cng) + k_dp * kc * dh;
                    double a12 = tensorDot(nf, cnmr);
                    double a21 = tensorDot(nmr, cng);
                    if (cstat > cres + tiny_num)
                        a21 += (ft0 / c0) * hdp * kc;
                    double a22 = tensorDot(nmr, cnmr)
                                 + P.smr / std::max(dt, tiny_num);
                    double det = a11 * a22 - a12 * a21;
                    if (std::abs(det) <= tiny_num) { failed = 1; break; }
                    double dl1 = (f1 * a22 - f2 * a12) / det;
                    double dl2 = (a11 * f2 - a21 * f1) / det;
                    if (dl1 < 0.0 && dl2 >= 0.0) {
                        dl1 = 0.0;
                        dl2 = std::max(f2 / std::max(a22, tiny_num), 0.0);
                    } else if (dl2 < 0.0 && dl1 >= 0.0) {
                        dl2 = 0.0;
                        dl1 = std::max(f1 / std::max(a11, tiny_num), 0.0);
                    } else if (dl1 < 0.0 && dl2 < 0.0) {
                        failed = 1;
                        break;
                    }
                    if (dl1 + dl2 <= tiny_num) { failed = 1; break; }
                    for (int i = 0; i < 6; ++i) {
                        sbar[i] -= dl1 * cng[i] + dl2 * cnmr[i];
                        depvp[i] += dl1 * ng[i] + dl2 * nmr[i];
                    }
                    lam_dp += dl1;
                    lam_mr += dl2;
                }
                niter += std::min(it, 80);
                if (it > 80) failed = 1;
                active = 3;
                did_work = true;
            } else if (fdp > tol) {
                // ---- return_dp ------------------------------------------
                int it;
                for (it = 1; it <= 60; ++it) {
                    double draw = c0 + hdp * (kap_dp_old + kc * lam_dp)
                                  + P.sdp * kc * lam_dp / std::max(dt, tiny_num);
                    double ctrial = std::max(cres, draw);
                    invariants(sbar, i1, q);
                    double f = q + alpha_dp * i1 - k_dp * ctrial;
                    if (f <= tol) break;
                    V6 nf, ng, cng;
                    dpGradients(sbar, alpha_dp, P.betaDP, nf, ng);
                    elasticApply(ng, bulk, shear, cng);
                    double deriv_h = draw > cres
                                         ? hdp + P.sdp / std::max(dt, tiny_num)
                                         : 0.0;
                    double denom = tensorDot(nf, cng) + k_dp * kc * deriv_h;
                    if (denom <= tiny_num) { failed = 1; break; }
                    double dl = f / denom;
                    if (dl <= 0.0) { failed = 1; break; }
                    for (int i = 0; i < 6; ++i) {
                        sbar[i] -= dl * cng[i];
                        depvp[i] += dl * ng[i];
                    }
                    lam_dp += dl;
                }
                niter += std::min(it, 60);
                if (it > 60) failed = 1;
                active = 2;
                did_work = true;
            } else if (fmr > tol) {
                // ---- return_mr ------------------------------------------
                double cstat = std::max(
                    cres, c0 + hdp * (kap_dp_old + kc * lam_dp));
                double ftstat = ft0 * cstat / std::max(c0, tiny_num);
                int it;
                for (it = 1; it <= 60; ++it) {
                    mrGradient(sbar, nmr, qpos);
                    double ftdyn = ftstat
                                   + P.smr * lam_mr / std::max(dt, tiny_num);
                    double f = qpos - ftdyn;
                    if (f <= tol) break;
                    if (qpos <= tiny_num) { failed = 1; break; }
                    V6 cnmr;
                    elasticApply(nmr, bulk, shear, cnmr);
                    double denom = tensorDot(nmr, cnmr)
                                   + P.smr / std::max(dt, tiny_num);
                    if (denom <= tiny_num) { failed = 1; break; }
                    double dl = f / denom;
                    if (dl <= 0.0) { failed = 1; break; }
                    for (int i = 0; i < 6; ++i) {
                        sbar[i] -= dl * cnmr[i];
                        depvp[i] += dl * nmr[i];
                    }
                    lam_mr += dl;
                }
                niter += std::min(it, 60);
                if (it > 60) failed = 1;
                active = 1;
                did_work = true;
            }
        }

        if (!did_work || failed != 0) break;
    }

    // Eq. (9)-(10): damage from the positive part of the vp increment,
    // unilateral nominal stress
    double dep_eq;
    positiveTensorNorm(depvp, dep_eq);
    double eqvt = eqvt_old + dep_eq;
    double omega = 1.0 - (1.0 - P.at + P.at * std::exp(-P.betaT * eqvt));
    omega = std::min(std::max(omega, 0.0), P.at);
    nominalStress(sbar, omega, stressNew);

    st.kapDP = kap_dp_old + kc * lam_dp;
    st.kapMR = kap_mr_old + lam_mr;
    st.eqvt = eqvt;
    st.epsv = epsv;
    st.pp = pp;
    for (int i = 0; i < 6; ++i) st.sbar[i] = sbar[i];
    st.active = active;
    st.failed = failed;
    (void)niter;
}

} // namespace sk11

// ---------------------------------------------------------------------------
// dfhk — line-by-line port of VUMATS/dfh/vumat_kstdfh.f (DP-DFH Bohus).
// The helpers (Jacobi eigensolver, rotations, Euler ZYX, xorshift64 seed)
// are transliterated VERBATIM from the kst_* Fortran subroutines so the
// selftest traces superpose to floating-point accuracy: Eigen's
// eigensolver would give the same frame only up to roundoff, and the
// frozen-frame freeze is a one-shot decision that amplifies any difference.
// Voigt-6 order 11 22 33 12 23 13, tensorial shears, tension positive.
// ---------------------------------------------------------------------------
namespace dfhk {

using V6 = std::array<double, 6>;
using M3 = std::array<std::array<double, 3>, 3>;

struct Props {
    double E = 52.0e9, nu = 0.25, rho = 2620.0;
    double betaDeg = 51.7, dcoh = 153.3e6, psiDeg = 15.0;
    double m = 24.0, sigw = 120.0e6, zeff = 1.0e-9;   // SI: Pa, m^3
    double k = 0.38, S = 4.18879;                     // thesis card 4pi/3
    double deld = 1.0e9;                              // deletion OFF
    // ---- dilatance VARIABLE (2026-08-24, OPT-IN) ---------------------
    // Porte de la branche insertion-pointe (commit c2c2d42).
    // vumat_hole.f / vumat_kstdfh_psivar_phicap.f l.336-339 :
    //   psi = clamp(PSI0 - KPSI * pbar[MPa], 0, PSIMAX)
    // La constante #5 de la carte (15 deg) est alors MORTE. Avec les
    // valeurs de la these, psi ne descend sous PSIMAX qu au-dela de
    // pbar = (PSI0 - PSIMAX)/KPSI = 509 MPa : en dessous l ecoulement est
    // ASSOCIE (psi = beta). psiVar = false laisse psiDeg fixe, chemin
    // d origine bit-identique.
    bool psiVar = false;
    double psi0 = 160.345, kpsi = 0.213793, psiMax = 51.7;   // deg, deg/MPa
};

constexpr double DCAP = 0.9999;

// --- kst_seed: 64-bit spatial hash -> 3 sorted Weibull draws --------------
// Fortran hashes anint(x_mm * 1e3) (micrometre integers); in SI the same
// integers come from llround(x_m * 1e6). uint64 arithmetic reproduces the
// integer*8 wraparound bit for bit (ishft(h,-n) is a logical shift).
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

// --- kst_rot6: rotate a stress 6-vector ------------------------------------
// idir >= 0: components in the frozen frame (b = R^T a R);
// idir <  0: back to the co-rotational frame (b = R a R^T).
inline void rot6(const V6& s, const M3& R, int idir, V6& so) {
    double a[3][3] = {{s[0], s[3], s[5]},
                      {s[3], s[1], s[4]},
                      {s[5], s[4], s[2]}};
    double t[3][3], b[3][3];
    if (idir >= 0) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                t[i][j] = a[i][0] * R[0][j] + a[i][1] * R[1][j]
                        + a[i][2] * R[2][j];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                b[i][j] = R[0][i] * t[0][j] + R[1][i] * t[1][j]
                        + R[2][i] * t[2][j];
    } else {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                t[i][j] = a[i][0] * R[j][0] + a[i][1] * R[j][1]
                        + a[i][2] * R[j][2];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                b[i][j] = R[i][0] * t[0][j] + R[i][1] * t[1][j]
                        + R[i][2] * t[2][j];
    }
    so[0] = b[0][0];
    so[1] = b[1][1];
    so[2] = b[2][2];
    so[3] = 0.5 * (b[0][1] + b[1][0]);
    so[4] = 0.5 * (b[1][2] + b[2][1]);
    so[5] = 0.5 * (b[0][2] + b[2][0]);
}

// --- kst_eulr / kst_reul: rotation matrix <-> Euler ZYX --------------------
inline void eulr(const M3& R, double eul[3]) {
    double cb = std::sqrt(R[2][1] * R[2][1] + R[2][2] * R[2][2]);
    eul[1] = std::atan2(-R[2][0], cb);
    if (cb > 1.0e-9) {
        eul[0] = std::atan2(R[1][0], R[0][0]);
        eul[2] = std::atan2(R[2][1], R[2][2]);
    } else {
        eul[0] = std::atan2(-R[0][1], R[1][1]);
        eul[2] = 0.0;
    }
}

inline void reul(const double eul[3], M3& R) {
    double ca = std::cos(eul[0]), sa = std::sin(eul[0]);
    double cb = std::cos(eul[1]), sb = std::sin(eul[1]);
    double cg = std::cos(eul[2]), sg = std::sin(eul[2]);
    R[0][0] = ca * cb;
    R[0][1] = ca * sb * sg - sa * cg;
    R[0][2] = ca * sb * cg + sa * sg;
    R[1][0] = sa * cb;
    R[1][1] = sa * sb * sg + ca * cg;
    R[1][2] = sa * sb * cg - ca * sg;
    R[2][0] = -sb;
    R[2][1] = cb * sg;
    R[2][2] = cb * cg;
}

// --- kst_eig3: cyclic Jacobi for a symmetric 3x3 ---------------------------
inline void eig3(const V6& s, double pv[3], M3& vec) {
    double a[3][3] = {{s[0], s[3], s[5]},
                      {s[3], s[1], s[4]},
                      {s[5], s[4], s[2]}};
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k) vec[i][k] = (i == k) ? 1.0 : 0.0;
    for (int isweep = 0; isweep < 50; ++isweep) {
        double off = std::abs(a[0][1]) + std::abs(a[0][2])
                   + std::abs(a[1][2]);
        if (off < 1.0e-20) break;
        for (int p = 0; p < 2; ++p)
            for (int q = p + 1; q < 3; ++q) {
                if (std::abs(a[p][q]) <= 1.0e-20) continue;
                double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                double t = (theta >= 0.0 ? 1.0 : -1.0)
                           / (std::abs(theta)
                              + std::sqrt(theta * theta + 1.0));
                double c = 1.0 / std::sqrt(t * t + 1.0);
                double sn = t * c;
                double tau = sn / (1.0 + c);
                double h = t * a[p][q];
                a[p][p] -= h;
                a[q][q] += h;
                a[p][q] = 0.0;
                a[q][p] = 0.0;
                for (int i = 0; i < 3; ++i) {
                    if (i == p || i == q) continue;
                    double aip = a[i][p], aiq = a[i][q];
                    a[i][p] = aip - sn * (aiq + tau * aip);
                    a[p][i] = a[i][p];
                    a[i][q] = aiq + sn * (aip - tau * aiq);
                    a[q][i] = a[i][q];
                }
                for (int i = 0; i < 3; ++i) {
                    double aip = vec[i][p], aiq = vec[i][q];
                    vec[i][p] = aip - sn * (aiq + tau * aip);
                    vec[i][q] = aiq + sn * (aip - tau * aiq);
                }
            }
    }
    pv[0] = a[0][0];
    pv[1] = a[1][1];
    pv[2] = a[2][2];
}

// --- the point update (vumat_dfh main loop body, one point) ----------------
// st.t must already hold the END-of-increment total time (caller adds dt
// before the call — the driver passes the same convention as totalTime).
inline void updatePoint(double dt, const V6& deps, const Props& P,
                        double Vel, double ftScale, MatState::Dfh& st,
                        V6& snomOut) {
    const double THIRD = 1.0 / 3.0;
    double G = P.E / (2.0 * (1.0 + P.nu));
    double xK = P.E / (3.0 * (1.0 - 2.0 * P.nu));
    double alam = P.E * P.nu / ((1.0 + P.nu) * (1.0 - 2.0 * P.nu));
    double tanb = std::tan(P.betaDeg * M_PI / 180.0);
    double tanp = std::tan(P.psiDeg * M_PI / 180.0);   // ecrase si psiVar

    if (st.dead) {
        for (int j = 0; j < 6; ++j) snomOut[j] = 0.0;
        return;
    }
    double cwav = std::sqrt(P.E / P.rho);

    bool lfroz = st.ti[0] > 0.0;
    M3 R;
    double fi[3];

    // 1. effective stress at increment start (exact inverse of D)
    V6 sbar;
    for (int j = 0; j < 6; ++j) sbar[j] = st.snom[j];
    if (lfroz) {
        reul(st.eul, R);
        V6 sfr;
        rot6(sbar, R, 1, sfr);
        for (int i = 0; i < 3; ++i) {
            fi[i] = 1.0;
            if (sfr[i] > 0.0) fi[i] = 1.0 - std::min(st.Dv[i], DCAP);
        }
        sfr[0] /= fi[0];
        sfr[1] /= fi[1];
        sfr[2] /= fi[2];
        sfr[3] /= std::min(fi[0], fi[1]);
        sfr[4] /= std::min(fi[1], fi[2]);
        sfr[5] /= std::min(fi[0], fi[2]);
        rot6(sfr, R, -1, sbar);
    }

    // 2. elastic predictor
    double tr = deps[0] + deps[1] + deps[2];
    sbar[0] += alam * tr + 2.0 * G * deps[0];
    sbar[1] += alam * tr + 2.0 * G * deps[1];
    sbar[2] += alam * tr + 2.0 * G * deps[2];
    sbar[3] += 2.0 * G * deps[3];
    sbar[4] += 2.0 * G * deps[4];
    sbar[5] += 2.0 * G * deps[5];

    // 3. Drucker-Prager return (effective), COMPRESSION only (OPTION-3)
    double xI1 = sbar[0] + sbar[1] + sbar[2];
    double pbar = -xI1 * THIRD;
    double dev1 = sbar[0] + pbar;
    double dev2 = sbar[1] + pbar;
    double dev3 = sbar[2] + pbar;
    double sJ2 = 0.5 * (dev1 * dev1 + dev2 * dev2 + dev3 * dev3)
               + sbar[3] * sbar[3] + sbar[4] * sbar[4] + sbar[5] * sbar[5];
    double q = std::sqrt(3.0 * std::max(sJ2, 0.0));
    double f = q - pbar * tanb - P.dcoh;
    if (P.psiVar) {                       // psi(pbar), cf. Props::psiVar
        double ps = P.psi0 - P.kpsi * (pbar * 1.0e-6);     // pbar en MPa
        if (ps > P.psiMax) ps = P.psiMax;
        if (ps < 0.0) ps = 0.0;
        tanp = std::tan(ps * M_PI / 180.0);
    }
    if (f > 0.0 && pbar > 0.0) {
        double dlam = f / (3.0 * G + xK * tanb * tanp);
        double qn = q - 3.0 * G * dlam;
        double facd, pnew;
        if (qn > 0.0 && q > 1.0e-12) {
            facd = qn / q;
            pnew = pbar + xK * tanp * dlam;
        } else {
            facd = 0.0;
            pnew = -P.dcoh / tanb;
            dlam = q / (3.0 * G);
        }
        sbar[0] = facd * dev1 - pnew;
        sbar[1] = facd * dev2 - pnew;
        sbar[2] = facd * dev3 - pnew;
        sbar[3] *= facd;
        sbar[4] *= facd;
        sbar[5] *= facd;
        st.peeq += dlam;
    }

    // 4a. freeze the frame at first initiation
    if (!lfroz) {
        double pv[3];
        M3 vec;
        eig3(sbar, pv, vec);
        int imax = 0;
        if (pv[1] > pv[imax]) imax = 1;
        if (pv[2] > pv[imax]) imax = 2;
        if (pv[imax] > st.smaxh) st.smaxh = pv[imax];
        if (pv[imax] >= st.sc[0]) {
            int imin = imax;
            for (int i = 0; i < 3; ++i) {
                if (i == imax) continue;
                if (imin == imax) imin = i;
                if (pv[i] < pv[imin]) imin = i;
            }
            int imid = 3 - imax - imin;
            for (int i = 0; i < 3; ++i) {
                R[i][0] = vec[i][imax];
                R[i][1] = vec[i][imid];
                R[i][2] = vec[i][imin];
            }
            double cx = R[1][0] * R[2][1] - R[2][0] * R[1][1];
            double cy = R[2][0] * R[0][1] - R[0][0] * R[2][1];
            double cz = R[0][0] * R[1][1] - R[1][0] * R[0][1];
            if (cx * R[0][2] + cy * R[1][2] + cz * R[2][2] < 0.0) {
                R[0][2] = -R[0][2];
                R[1][2] = -R[1][2];
                R[2][2] = -R[2][2];
            }
            eulr(R, st.eul);
            st.ti[0] = std::max(st.t, 1.0e-30);
            lfroz = true;
        }
    }

    // 4b. per-direction initiation/growth + 5. nominal stress
    for (int j = 0; j < 6; ++j) snomOut[j] = sbar[j];
    if (lfroz) {
        reul(st.eul, R);
        V6 sfr;
        rot6(sbar, R, 1, sfr);
        for (int i = 0; i < 3; ++i) {
            double sn = sfr[i];
            if (sn > st.smaxh) st.smaxh = sn;
            if (st.ti[i] <= 0.0 && sn >= st.sc[i])
                st.ti[i] = std::max(st.t, 1.0e-30);
            if (st.ti[i] > 0.0 && sn > 0.0 && st.Dv[i] < DCAP) {
                double sigwLoc = P.sigw * ftScale;
                double xlam = std::pow(sn / sigwLoc, P.m) / P.zeff;
                if (xlam * Vel < 1.0) xlam = 1.0 / Vel;
                double xx = std::pow(-std::log(1.0 - st.Dv[i]), THIRD);
                xx += std::pow(P.S * xlam, THIRD) * P.k * cwav * dt;
                st.Dv[i] = 1.0 - std::exp(-xx * xx * xx);
                if (st.Dv[i] > DCAP) st.Dv[i] = DCAP;
            }
        }
        for (int i = 0; i < 3; ++i) {
            fi[i] = 1.0;
            if (sfr[i] > 0.0) fi[i] = 1.0 - st.Dv[i];
        }
        sfr[0] *= fi[0];
        sfr[1] *= fi[1];
        sfr[2] *= fi[2];
        sfr[3] *= std::min(fi[0], fi[1]);
        sfr[4] *= std::min(fi[1], fi[2]);
        sfr[5] *= std::min(fi[0], fi[2]);
        rot6(sfr, R, -1, snomOut);
    }

    // optional deletion (dfhDeld; OFF at the VUMAT default 1e9)
    double dmax = std::max({st.Dv[0], st.Dv[1], st.Dv[2]});
    if (dmax >= P.deld) {
        st.dead = true;
        for (int j = 0; j < 6; ++j) snomOut[j] = 0.0;
    }
    for (int j = 0; j < 6; ++j) st.snom[j] = snomOut[j];
}

} // namespace dfhk

namespace {

class Saksala2011Law : public MatLaw {
public:
    Saksala2011Law(const Material& m, const sk11::Props& p)
        : MatLaw(m), P_(p) {}

    Eigen::Matrix3d stress(const Eigen::Matrix3d& eps, MatState& s,
                           double dt, double) const override {
        auto& st = s.sk;
        if (!st.init) {
            st.init = true;
            st.pp = P_.pp0;
            // FIELD-like initialization: the solver's per-element factor
            // scales BOTH local strengths (one common draw; the paper uses
            // two independent shifted-Weibull draws — one factor is the
            // simpler sandbox variant)
            st.ftLoc = P_.ft0 * s.ftScale;
            st.c0Loc = P_.c0 * s.ftScale;
        }
        sk11::V6 eps6 = {eps(0, 0), eps(1, 1), eps(2, 2),
                         eps(0, 1), eps(1, 2), eps(0, 2)};
        sk11::V6 deps;
        for (int i = 0; i < 6; ++i) {
            deps[i] = eps6[i] - st.epsPrev[i];
            st.epsPrev[i] = eps6[i];
        }
        sk11::V6 sig;
        sk11::updatePoint(std::max(dt, 1e-30), deps, P_, st, sig);
        // expose the standard outputs
        double omega = 1.0 - (1.0 - P_.at
                              + P_.at * std::exp(-P_.betaT * st.eqvt));
        s.D = std::min(std::max(omega, 0.0), P_.at);
        s.epvEq = st.eqvt;
        s.pc = st.pp;
        s.kappa = st.kapDP;    // the compression-side field of the 2011
                               // model (it has no separate omega_c)
        Eigen::Matrix3d out;
        out << sig[0], sig[3], sig[5],
               sig[3], sig[1], sig[4],
               sig[5], sig[4], sig[2];
        return out;
    }

    std::string name() const override { return "saksala2011"; }

private:
    sk11::Props P_;
};

// ---------------------------------------------------------------------------
class DpDfhLaw : public MatLaw {
public:
    DpDfhLaw(const Material& m, const dfhk::Props& p) : MatLaw(m), P_(p) {}

    Eigen::Matrix3d stress(const Eigen::Matrix3d& eps, MatState& s,
                           double dt, double lc) const override {
        auto& st = s.dfh;
        double Vel = lc * lc * lc;
        if (Vel < 1.0e-12) Vel = 1.0;
        if (!st.seeded) {
            st.seeded = true;
            double sigk = P_.sigw * std::pow(P_.zeff / Vel, 1.0 / P_.m)
                          * s.ftScale;
            dfhk::seed(s.x0.x(), s.x0.y(), s.x0.z(), sigk, P_.m, st.sc);
        }
        dfhk::V6 eps6 = {eps(0, 0), eps(1, 1), eps(2, 2),
                         eps(0, 1), eps(1, 2), eps(0, 2)};
        dfhk::V6 deps;
        for (int i = 0; i < 6; ++i) {
            deps[i] = eps6[i] - st.epsPrev[i];
            st.epsPrev[i] = eps6[i];
        }
        st.t += std::max(dt, 1.0e-30);
        dfhk::V6 snom;
        dfhk::updatePoint(std::max(dt, 1.0e-30), deps, P_, Vel, s.ftScale,
                          st, snom);
        // standard outputs: damage slot = DMAX (SDV2), epvEq = PEEQ (SDV3),
        // kappa slot = SMAXH (SDV16 diagnostic)
        s.D = std::max({st.Dv[0], st.Dv[1], st.Dv[2]});
        s.epvEq = st.peeq;
        s.kappa = st.smaxh;
        if (st.dead) { s.eroded = true; s.eroCode = 4; }
        Eigen::Matrix3d out;
        out << snom[0], snom[3], snom[5],
               snom[3], snom[1], snom[4],
               snom[5], snom[4], snom[2];
        return out;
    }

    std::string name() const override { return "dpdfh"; }

private:
    dfhk::Props P_;
};

// ===========================================================================
// cdp — Concrete Damaged Plasticity (Lubliner et al. 1989 ; Lee & Fenves
// 1998), la loi NATIVE d'Abaqus (*CONCRETE DAMAGED PLASTICITY), portee le
// 2026-09-04 pour la calibration triaxiale Red Bohus
// (CONTINUUM/calib_bohus_triax/cdp_pente, cdp_inverse, cdp_rockim). Parite
// Abaqus visee : memes cartes, memes conversions de tables, meme
// interpolation (lineaire dans la variable plastique INTERNE — la « lecture
// B » de la doc « Defining tension stiffening / compressive behavior »).
//
// CONVENTIONS. Tenseurs traction positive ; p_bar = -tr(sigma_bar)/3 (donc
// COMPRESSION POSITIVE, comme pbar du port dpdfh) ; q_bar = sqrt(3/2 s:s) ;
// mu = G_ (cisaillement), K = K_ ; SI. Eigen rend les valeurs propres
// CROISSANTES : indice 2 = sigma_hat_max, indice 0 = sigma_hat_min ; les
// d eps_hat suivent le meme ordre (coaxialite, fonction croissante de s_i).
//
// SURFACE (espace effectif) :
//   F = [q - 3 alpha p + beta <s_max> - gamma <-s_max>]/(1 - alpha) - sigma_bar_c
//   alpha = (fb0/fc0 - 1)/(2 fb0/fc0 - 1),  gamma = 3 (1 - Kc)/(2 Kc - 1),
//   beta  = sigma_bar_c/sigma_bar_t (1 - alpha) - (1 + alpha)
// En compression triaxiale (s_max <= 0, seul gamma agit) la surface est la
// DROITE q_pic = fc0_pic + m sigma3, m = (3 alpha + gamma)/(1 - alpha) ; en
// nominal (r = 0, d = d_c, confinement pilote nominal) toute la courbe est
// la table sigma_c(eps_c_pl) translatee de m sigma3 — EXACTEMENT (banc (a)).
// POTENTIEL (non associe, hyperbolique — pas Mohr-Coulomb : psi = 35 deg
// donne -d eps_v/d eps_1 = 0,913, soit psi_MC equivalent 18,3 deg) :
//   G = sqrt(a^2 + q^2) - p tan(psi),  a = ecc ft tan(psi)
//   dG/dsigma = 1,5 s/sqrt(a^2 + q^2) + tan(psi)/3 I     (tr = tan psi > 0)
// RETOUR en espace spectral (coaxial : D0 et G isotropes), pleinement
// implicite en lambda (r, beta, sigma_bar_c evalues a n+1, Lee & Fenves 2001) :
//   p = p_tr + K tan(psi) lambda            <- SIGNE + : la dilatance FAIT
//                                              MONTER la pression effective
//   q (1 + 3 mu lambda/sqrt(a^2 + q^2)) = q_tr   (g croissante concave :
//       Newton depuis q0 = max(q_tr - 3 mu lambda, 0), monotone par la gauche)
//   d eps_hat_i = lambda [1,5 sd_i/sqrt(a^2 + q^2) + tan(psi)/3]
//   d eps_t_pl = max(0, r d eps_hat_max), d eps_c_pl = max(0, -(1-r) d eps_hat_min)
//   r = sum <s_i>/sum |s_i|  (Lee-Fenves : h = [r, 0, -(1-r)] sur les
//   principales EFFECTIVES ; clamp car d eps_hat_min > 0 pres de l'axe
//   hydrostatique et d eps_hat_max >= lambda t/3 > 0 toujours)
// F(lambda) n'est que C0 (noeuds de tables, s_max = 0, r) : crochet par
// doublement de hi puis regula falsi Illinois, JAMAIS de Newton pur. F -> -inf
// quand lambda -> inf (p -> +inf, q -> 0) : le crochet existe toujours.
// ENDOMMAGEMENT, a chaque appel (d depend de r, donc de la contrainte finale :
// ce n'est PAS une variable d'etat, seuls eps_t_pl / eps_c_pl le sont) :
//   d = 1 - (1 - s_t d_c)(1 - s_c d_t), s_t = 1 - w_t r, s_c = 1 - w_c (1-r)
//   sigma_nom = (1 - d) sigma_bar   (r = 0 : d = d_c, fissures refermees ;
//   r = 1 : 1-d = (1-d_c)(1-d_t)). Plafond d, d_t, d_c <= 0,99 (max
//   degradation d'Abaqus/Explicit).
// TABLES. Compression (*CONCRETE COMPRESSION HARDENING / DAMAGE) en
// deformation INELASTIQUE eps_in, abscisses FUSIONNEES (Abaqus accepte des
// abscisses differentes), converties a l'init :
//   eps_c_pl = eps_in - d_c/(1 - d_c) sigma_c/E0       (doc Abaqus)
// puis sigma_c, d_c lineaires en eps_c_pl, constantes apres le dernier
// noeud ; sigma_bar_c = sigma_c/(1 - d_c) (rapport de deux affines).
// Traction en OUVERTURE u_ck (*CONCRETE TENSION STIFFENING type=GFI :
// sigma_t = ft (1 - u/u_t0), u_t0 = 2 Gf/ft ; ou type=DISPLACEMENT =
// cdpTension = table ; le type STRAIN n'est PAS couvert) et *CONCRETE
// TENSION DAMAGE type=DISPLACEMENT ; noeuds fusionnes {0} U {u sigma_t} U
// {u d_t} U {u_t0}, convertis A LA VOLEE par element :
//   eps_t_pl,k = u_k/lc - c_k,  c_k = d_t,k/(1 - d_t,k) sigma_t,k/E0
// (convention « lc partout » = Abaqus avec l0 = lc ; identique aux decks mm
// a < 1 %, et le seul choix sense en SI ou l0 = 1 m vaudrait 1000 lc).
// Plancher sigma_t >= 1e-3 ft applique au NOMINAL avant division par
// (1 - d_t) (sigma_bar_t,res = 170 kPa, beta ~ 650 : terme raide dans F).
// CHOIX lc : celui du solveur, lc = V0^(1/3) (Abaqus : longueur
// caracteristique du tetraedre, V^(1/3) de memoire du manuel).
// ENERGIE (lecture B, carte GFI + d_t a 2 noeuds) : en traction uniaxiale
// W_tot x lc - lc x wDamT = Gf (identite exacte a O(c_k), banc (c)) ; l'aire
// totale x lc vaut Gf + lc wDamT (+0,36 % a 1 mm, +1,1 % a 3 mm). La
// dissipation tractive est PLASTIQUE (eps_t_pl ~ u_ck/lc), wDamT n'en porte
// que 0,36 J/m^2 : le canal erodeWfrac est donc INERTE en CDP (documente).
// POINT AVEUGLE (parite, pas une erreur) : l'adoucissement COMPRESSIF est en
// deformation, non regularise en maillage ; remettre la carte a l'echelle en
// contrainte sans toucher aux abscisses raidit la chute nominale post-pic
// (-0,37 E historique, -2,0 E carte inverse) : a l'echelle d'une eprouvette
// c'est un snap-back structurel, la « chute verticale » Abaqus (391,8 ->
// 19,2 MPa en un cadre) est attendue A L'IDENTIQUE dans rockim. Une
// regularisation (cdpCompBand) serait une cle opt-in ULTERIEURE.
// GARDES a l'init (lcMax, exception nommant le point) : (G1) lcMax <=
// 2 E0 Gf/ft^2 [gfi] ou max |d sigma_t/du| <= E0/lcMax [table] ; (G2)
// noeuds convertis strictement croissants ; (G3) snap-back EFFECTIF :
// pente locale -d sigma_bar_t/d eps_t_pl (aux deux bouts de chaque segment,
// sigma_bar_t = sigma_t/(1 - d_t) n'etant pas affine) <= 0,9 (2 mu + K t)/
// (1 + t/3), et -d sigma_bar_c/d eps_c_pl <= 0,9 (3 mu + 3 alpha K t)/
// ((1 - alpha)(1 - t/3)). Bornes demontrees pour les chemins uniaxiaux ;
// facteur 0,9 de prudence hors uniaxial. Carte Bohus : lc < 10,7 mm.
// EROSION (cles existantes) : spall = d_t >= erodeD en traction nette
// (avertissement si erodeD > d_t,max : canal inerte, poser erodeD <= 0,95),
// erodeWfrac (garde, inerte), erodeEpv sur eps_c_pl, erodeDc sur
// d_c/d_c,max NORMALISE (garde compDamage levee pour cdp) ; eroCode 1/2/3.
// CHAMPS PARTAGES : s.D = d_t (monotone : VTU « damage », V_D09, canal
// spall), s.Dc = d_c, s.epvEq = eps_c_pl, s.kappa = eps_t_pl (VTU kapDP),
// s.cdp.d = d total (diagnostic).
// ---- REVUE du 2026-09-04 (constats traites) --------------------------------
// HETEROGENEITE (ftScale du solveur, matWeibullM) : HONOREE comme par dpr
// (E1). Par element, s = s.ftScale : ft_loc = s ft, a_loc = ecc ft_loc tan psi,
// table de traction sigma_t x s et u_ck x kappa avec kappa = Gf_loc/(Gf s) :
// 1/s en portee `strength` (Gf conserve, u_t0 en 1/s) ou 1 en `strengthGf`
// (ft et Gf suivent s, u_t0 invariant) ; d_t(u) suit la meme dilatation des
// abscisses (courbe auto-semblable en u). La table de COMPRESSION n'est pas
// touchee (comme la cohesion de dpr : le facteur de Weibull porte sur la
// resistance tractive). Les gardes G1-G3 dependent de (lc, s) : verifiees a
// l'init a (lcMax, 1) et RE-VERIFIEES au premier appel de tout element a
// s != 1 (exception nommant lc et ftScale) ; G3 se durcit en s^2 (l'element
// se comporte comme un element de taille s^2 lc) : carte Bohus lc < 10,7 mm
// / s^2 — baisser d_t,max (0,95 -> 0,9 double la marge) ou raffiner.
// A s = 1 l'arithmetique est IDENTIQUE (x 1,0 exact) : selftest-cdp inchange.
// VISCOSITE (cdpViscosity, s ; 0 = rate-independant, defaut, bit-identique) :
// regularisation viscoplastique de Duvaut-Lions de la CDP d'Abaqus/STANDARD
// (Theory 4.5.2) : eps_pl,v += dt/(mu + dt) (eps_pl - eps_pl,v), d_v idem,
// sigma_nom = (1 - d_v) D0 (eps - eps_pl,v) ; le retour inviscide (eps_pl,
// d, F) est inchange. Euler implicite : le retard stationnaire vaut
// EXACTEMENT mu x vitesse (banc (k)). ATTENTION (revue, 2e passe) : les decks
// de reference (*Dynamic, Explicit) portent mu = 5e-5 s sur la ligne
// *CONCRETE DAMAGED PLASTICITY, mais Abaqus/EXPLICIT IGNORE ce parametre
// (Keywords Reference, 5e donnee : « used ... in Abaqus/Standard analyses.
// This parameter is ignored in Abaqus/Explicit ») : la parite avec les runs
// de reference est le DEFAUT rate-independant, pas cdpViscosity = 5e-5. Si
// on l'active tout de meme (Standard), a 0,75/s au pic la surcontrainte
// n'est PAS E0 mu epsdot (~3 MPa) mais ~6 fois plus — la contrainte laterale
// visqueuse mu (D0 epsdot_pl)_lat doit etre compensee par un confinement
// inviscide c que la pente m amplifie : Delta q = (m+1) c + mu lamdot
// |(D0 n)_ax|, c = mu lamdot (D0 n)_lat, n = direction d'ecoulement (banc (k),
// formule exacte a l'etat stationnaire ; 9,9 MPa = +4,9 % sur le pic
// triaxial a 20 MPa). dt doit etre > 0 quand mu > 0 (matpoint : mpDt =
// Delta eps / epsdot).
// CLES : toute cle commencant par « cdp » qui n'est pas dans la liste connue
// est REFUSEE (exception nommant la cle) — Config ignore les inconnues en
// silence, un outil de calibration ne peut pas se le permettre.
// PLANCHER sigma_t >= 1e-3 ft, ARTEFACT documente : une fois le dernier noeud
// de traction atteint, rho = q/sqrt(a^2+q^2) ~ 0,3 et le potentiel
// hyperbolique coule dans les TROIS directions, d eps_hat_lat = lambda
// (-rho/2 + t/3) > 0 des que rho < 2 t/3 = 0,467 : l'element fissure GONFLE
// lateralement (eps_lat remonte de -4,5e-3 a +1,9e-3 entre eps_ax 0,021 et
// 0,06 a lc = 1 mm) et sa deformation volumique plastique croit sans borne
// (lambda t par pas), wPlas depasse Gf/lc. C'est la parite Abaqus (meme
// potentiel, meme residu de table) ; en fem3d un element fissure pousse
// donc ses voisins : le garde-fou est erodeD <= d_t,max (retrait de
// l'element en traction nette), a poser explicitement.
// ---- CRACK BAND EN COMPRESSION (cdpCompLength, m ; 0 = absent = defaut =
//      bit-identique, tables lues telles quelles) — 2026-09-04, soir --------
// Decision de Fernando : la chute post-pic d'un triaxial est STRUCTURELLE
// (bande localisee), pas une propriete du point materiel ; la table CDP de
// compression etant en deformation et non regularisee, l'energie dissipee
// dans la bande depend de la taille d'element (le point aveugle ci-dessus).
// Avec cdpCompLength = L_ref > 0, les abscisses eps_in des tables de
// compression (sigma_c ET d_c, noeuds fusionnes) sont lues comme definies a
// la longueur L_ref et SEULE LA BRANCHE POST-PIC est remise a l'echelle de
// l'element (pic = premier noeud fusionne ou sigma_c nominal est maximal) :
//   eps_in,loc = eps_in,pic + (eps_in - eps_in,pic) L_ref/lc   (eps_in > pic)
// la branche pre-pic (ecrouissage homogene) n'est PAS touchee, puis la
// conversion d'Abaqus s'applique au noeud remis a l'echelle :
//   eps_c_pl,k(lc) = eps_in,loc,k - c_k,  c_k = d_c,k/(1-d_c,k) sigma_c,k/E0
// (c_k ne depend pas de lc : la conversion est faite A LA VOLEE par element,
// comme la traction — compNode / compState(epl, lc) ; a cle absente,
// compState lit les noeuds convertis a l'init, arithmetique IDENTIQUE).
// Consequence : le deplacement inelastique post-pic u_in = (eps_in,loc -
// eps_in,pic) lc = (eps_in - eps_in,pic) L_ref est INVARIANT, donc l'energie
// d'adoucissement par unite d'aire de bande G_c = int (sigma_c - sigma_res)
// du_in est la meme quelle que soit la maille — crack band de Bazant-Oh en
// compression, meme logique que compGIIc de la brique omega_c (§5.18). En
// compression uniaxiale le point materiel restitue EXACTEMENT la table
// (sigma_ax = sigma_c(eps_c_pl), eps_ax = eps_in,loc + sigma/E0), d'ou la
// mesure du banc (m) : [int (sigma - sigma_res) d eps + (sigma_pic -
// sigma_res)^2/(2 E0)] x lc = G_c. Gardes a lcMax (la plus grande maille =
// segments les plus courts), levees a l'init seulement si la cle est > 0 :
//   (G2c) noeuds convertis strictement croissants ;
//   (G3)  pente effective -d sigma_bar_c/d eps_c_pl <= limC (retour) ;
//   (G4)  SNAP-BACK au point materiel : la deformation totale eps = eps_c_pl
//         + sigma_bar_c/E0 = eps_in,loc + sigma_c/E0 doit croitre avec
//         eps_c_pl, soit -d sigma_bar_c/d eps_c_pl <= E0 localement (les
//         deux bouts de chaque segment, sigma_bar = sigma/(1-d) n'etant pas
//         affine) ; en corde c'est -Delta sigma_c/Delta eps_in,loc <= E0,
//         la condition d'Abaqus sur sa table. NOTE : le libelle « -d sigma_c/
//         d eps_pl <= E0 » (nominal sur la plastique) n'est PAS ce critere :
//         il vaut 80,5 GPa > E0 a lc = 4 mm / L_ref 2 mm sur la carte
//         historique alors que la corde en eps_in vaut 41,6 GPa (aucun
//         snap-back) — la garde implemente le critere exact.
// Le seuil limC du retour (141,6 GPa) est plus lache que E0 (77,7) : entre
// les deux le retour converge mais la reponse en deformation totale SAUTE.
// ===========================================================================

// « v:x v:x ... » -> paires (v, x) dans l'ordre ecrit, lecture stricte par
// jeton (std::stod avec consommation complete : une virgule decimale ou un
// jeton sans ':' est signale en nommant la cle ET le jeton). Config::gets
// conserve les espaces internes, aucun ';' n'est necessaire.
std::vector<std::pair<double, double>> cdpParsePairs(const std::string& key,
                                                     const std::string& val) {
    std::vector<std::pair<double, double>> out;
    std::istringstream ss(val);
    std::string tok;
    auto num = [&](const std::string& s) {
        std::size_t used = 0;
        double d = 0.0;
        try { d = std::stod(s, &used); }
        catch (const std::exception&) { used = 0; }
        if (used == 0 || used != s.size())
            throw std::runtime_error("cdp: key '" + key + "' — token '" + tok
                                     + "' is not a number (decimal COMMA "
                                       "instead of point?)");
        return d;
    };
    while (ss >> tok) {
        auto c = tok.find(':');
        if (c == std::string::npos || c == 0 || c + 1 >= tok.size())
            throw std::runtime_error("cdp: key '" + key + "' — token '" + tok
                                     + "' must read value:abscissa");
        double v = num(tok.substr(0, c));
        double x = num(tok.substr(c + 1));
        out.emplace_back(v, x);
    }
    if (out.empty())
        throw std::runtime_error("cdp: key '" + key + "' is empty");
    return out;
}

// interpolation lineaire y(x) sur une table a abscisses croissantes,
// constante au-dela des deux bouts (regle Abaqus pour TOUTE table)
double cdpInterp(const std::vector<double>& x, const std::vector<double>& y,
                 double xq) {
    if (xq <= x.front()) return y.front();
    if (xq >= x.back()) return y.back();
    std::size_t k = 1;
    while (k + 1 < x.size() && x[k] < xq) ++k;
    double w = (xq - x[k - 1]) / (x[k] - x[k - 1]);
    return y[k - 1] + w * (y[k] - y[k - 1]);
}

class CdpLaw : public MatLaw {
public:
    CdpLaw(const Material& m, const Config& c, double lcMax, double erodeD,
           double erodeEpv, double erodeDc, double erodeWfrac)
        : MatLaw(m), erodeD_(erodeD), erodeEpv_(erodeEpv), erodeDc_(erodeDc),
          erodeWfrac_(erodeWfrac) {
        // ---- 1. cles (defauts = carte historique Red Bohus des decks) ----
        double psiDeg = c.getd("cdpDilationDeg", 35.0);
        ecc_ = c.getd("cdpEcc", 0.1);
        fbfc_ = c.getd("cdpFbFc", 1.16);
        Kc_ = c.getd("cdpKc", 0.667);
        wt_ = c.getd("cdpWt", 0.0);
        wc_ = c.getd("cdpWc", 1.0);
        std::string hard = c.gets("cdpHardening",
            "50.6e6:0 126.6e6:0.0008 60e6:0.004 10e6:0.012");
        std::string cdam = c.gets("cdpCompDamage",
            "0:0 0:0.0008 0.5:0.004 0.92:0.012");
        std::string tmode = c.gets("cdpTension", "gfi");
        std::string ttab = c.gets("cdpTensionTable", "");
        std::string tdam = c.gets("cdpTensionDamage", "0:0 0.95:2.35e-5");
        muV_ = c.getd("cdpViscosity", 0.0);
        // crack band en compression (2026-09-04 soir) : 0 = absent = tables
        // lues telles quelles (bit-identique)
        compLen_ = c.getd("cdpCompLength", 0.0);
        ft_ = m.ft; Gf_ = m.Gf; E0_ = m.E;
        // cap volumique de compaction (2026-09-04, percussion quart de bloc) :
        // cdpCap = false = absent = bit-identique ; la meme brique que le cap
        // de dpr/saksala (PlasticDamageLaw::stress, « pressure cap »), sur la
        // pression EFFECTIVE du predicteur, AVANT le retour de Lubliner
        cap_ = c.getb("cdpCap", false);
        capP0_ = c.getd("cdpCapP0", 0.0);
        capH_ = c.getd("cdpCapH", K_);                  // defaut K = E/(3(1-2nu))
        // revue 2026-09-04 : une cle cdp* mal orthographiee (cdpKC,
        // cdpHardenning...) rendrait la carte historique en silence
        {
            static const char* known[] = {
                "cdpDilationDeg", "cdpEcc", "cdpFbFc", "cdpKc", "cdpWt", "cdpWc",
                "cdpHardening", "cdpCompDamage", "cdpTension", "cdpTensionTable",
                "cdpTensionDamage", "cdpViscosity", "cdpCompLength",
                "cdpCap", "cdpCapP0", "cdpCapH"};
            for (const std::string& k : c.keysWithPrefix("cdp")) {
                bool okKey = false;
                for (const char* kk : known) if (k == kk) okKey = true;
                if (!okKey) {
                    std::string list;
                    for (const char* kk : known) list += std::string(" ") + kk;
                    throw std::runtime_error("cdp: unknown key '" + k
                        + "' (typo?) — known cdp keys :" + list);
                }
            }
        }
        if (!(muV_ >= 0.0))
            throw std::runtime_error("cdp: cdpViscosity >= 0 [s] required (0 = "
                                     "rate-independent)");
        if (cap_) {
            if (!c.has("cdpCapP0") || !(capP0_ > 0.0) || !std::isfinite(capP0_))
                throw std::runtime_error("cdp: cdpCap = true requires cdpCapP0 > 0 "
                                         "[Pa] (initial crush pressure of the "
                                         "volumetric cap, effective pressure)");
            if (!(capH_ >= 0.0) || !std::isfinite(capH_))
                throw std::runtime_error("cdp: cdpCapH >= 0 [Pa] required (cap "
                                         "hardening modulus, default K)");
        } else if (c.has("cdpCapP0") || c.has("cdpCapH")) {
            throw std::runtime_error("cdp: cdpCapP0 / cdpCapH given without "
                                     "cdpCap = true (set cdpCap = true or remove "
                                     "the keys — nothing is applied in silence)");
        }
        if (!(compLen_ >= 0.0) || !std::isfinite(compLen_))
            throw std::runtime_error("cdp: cdpCompLength >= 0 [m] required (0 = "
                                     "absent : tables de compression lues "
                                     "telles quelles, non regularisees)");

        // ---- 2. validation ----------------------------------------------
        if (!(psiDeg > 0.0 && psiDeg < 56.0))
            throw std::runtime_error("cdp: 0 < cdpDilationDeg < 56 required "
                "(psi > 0 : retour hydrostatique ; psi >= 56,3 deg degenere "
                "l'ecrouissage equibiaxial)");
        if (!(ecc_ > 0.0))
            throw std::runtime_error("cdp: cdpEcc > 0 required (sommet lisse, "
                                     "aucun cas apex)");
        if (!(fbfc_ > 1.0))
            throw std::runtime_error("cdp: cdpFbFc (fb0/fc0) > 1 required");
        if (!(Kc_ > 0.5 && Kc_ <= 1.0))
            throw std::runtime_error("cdp: 0.5 < cdpKc <= 1 required");
        if (!(wt_ >= 0.0 && wt_ <= 1.0 && wc_ >= 0.0 && wc_ <= 1.0))
            throw std::runtime_error("cdp: cdpWt, cdpWc in [0, 1]");
        if (tmode != "gfi" && tmode != "table")
            throw std::runtime_error("cdp: cdpTension must be gfi | table");
        if (tmode == "table" && ttab.empty())
            throw std::runtime_error("cdp: cdpTension = table requires "
                                     "cdpTensionTable = \"sigma:u_ck ...\"");
        gfi_ = tmode == "gfi";

        // ---- 3. constantes ------------------------------------------------
        t_ = std::tan(psiDeg * M_PI / 180.0);
        alpha_ = (fbfc_ - 1.0) / (2.0 * fbfc_ - 1.0);
        gamma_ = 3.0 * (1.0 - Kc_) / (2.0 * Kc_ - 1.0);
        a_ = ecc_ * ft_ * t_;
        mSlope_ = (3.0 * alpha_ + gamma_) / (1.0 - alpha_);

        // ---- 4. table de compression (abscisses fusionnees, converties) ----
        auto H = cdpParsePairs("cdpHardening", hard);      // (sigma_c, eps_in)
        auto Dc = cdpParsePairs("cdpCompDamage", cdam);    // (d_c, eps_in)
        std::vector<double> hx, hs, dx, dd;
        for (std::size_t k = 0; k < H.size(); ++k) {
            if (!(H[k].first > 0.0))
                throw std::runtime_error("cdp: cdpHardening point "
                    + std::to_string(k) + " : sigma_c must be > 0");
            if (k > 0 && !(H[k].second > H[k - 1].second))
                throw std::runtime_error("cdp: cdpHardening point "
                    + std::to_string(k) + " : eps_in must increase strictly");
            hs.push_back(H[k].first); hx.push_back(H[k].second);
        }
        for (std::size_t k = 0; k < Dc.size(); ++k) {
            if (!(Dc[k].first >= 0.0 && Dc[k].first <= 0.99))
                throw std::runtime_error("cdp: cdpCompDamage point "
                    + std::to_string(k) + " : 0 <= d_c <= 0.99 required");
            if (k > 0 && !(Dc[k].second > Dc[k - 1].second))
                throw std::runtime_error("cdp: cdpCompDamage point "
                    + std::to_string(k) + " : eps_in must increase strictly");
            dd.push_back(Dc[k].first); dx.push_back(Dc[k].second);
        }
        if (hx.front() != 0.0 || dx.front() != 0.0)
            throw std::runtime_error("cdp: the first abscissa of cdpHardening "
                                     "and cdpCompDamage must be 0 (Abaqus)");
        std::vector<double> ein = hx;
        ein.insert(ein.end(), dx.begin(), dx.end());
        std::sort(ein.begin(), ein.end());
        ein.erase(std::unique(ein.begin(), ein.end(),
                              [](double a, double b) {
                                  return std::abs(a - b) <= 1e-12 * std::max(a, b);
                              }), ein.end());
        fc0pic_ = 0.0; dcMax_ = 0.0;
        for (std::size_t k = 0; k < ein.size(); ++k) {
            double sc = cdpInterp(hx, hs, ein[k]);
            double dc = cdpInterp(dx, dd, ein[k]);
            double epl = ein[k] - dc / (1.0 - dc) * sc / E0_;
            if (k == 0 && std::abs(epl) > 1e-15)
                throw std::runtime_error("cdp: converted eps_c_pl(0) != 0 "
                                         "(d_c(0) must be 0)");
            if (k > 0 && !(epl > cEpl_.back()))
                throw std::runtime_error("cdp: converted compression node "
                    + std::to_string(k) + " (eps_in = " + std::to_string(ein[k])
                    + ") gives eps_c_pl = " + std::to_string(epl)
                    + " <= previous " + std::to_string(cEpl_.back())
                    + " — plastic strain values decreasing (Abaqus refuses "
                      "it too)");
            cEpl_.push_back(k == 0 ? 0.0 : epl);
            cSig_.push_back(sc);
            cD_.push_back(dc);
            // crack band en compression : eps_in fusionne et c_k = eps_in -
            // eps_c_pl conserves pour la remise a l'echelle par element
            cEin_.push_back(ein[k]);
            cC_.push_back(k == 0 ? 0.0 : ein[k] - epl);
            fc0pic_ = std::max(fc0pic_, sc / (1.0 - dc));
            dcMax_ = std::max(dcMax_, dc);
        }
        // pic NOMINAL de la table (premier noeud fusionne ou sigma_c est
        // maximal) : l'abscisse a partir de laquelle cdpCompLength agit
        kPic_ = 0;
        for (std::size_t k = 1; k < cSig_.size(); ++k)
            if (cSig_[k] > cSig_[kPic_]) kPic_ = k;
        einPic_ = cEin_[kPic_];

        // ---- 5. table de traction (noeuds en u_ck, c_k precalcules) -------
        std::vector<double> su, ss;                        // sigma_t(u)
        if (gfi_) {
            ut0_ = 2.0 * Gf_ / ft_;
            su = {0.0, ut0_}; ss = {ft_, 0.0};
        } else {
            auto T = cdpParsePairs("cdpTensionTable", ttab);   // (sigma_t, u)
            for (std::size_t k = 0; k < T.size(); ++k) {
                if (!(T[k].first > 0.0) && k == 0)
                    throw std::runtime_error("cdp: cdpTensionTable sigma_t(0) "
                                             "must be > 0");
                if (T[k].first < 0.0)
                    throw std::runtime_error("cdp: cdpTensionTable point "
                        + std::to_string(k) + " : sigma_t must be >= 0");
                if (k > 0 && !(T[k].second > T[k - 1].second))
                    throw std::runtime_error("cdp: cdpTensionTable point "
                        + std::to_string(k) + " : u_ck must increase strictly");
                ss.push_back(T[k].first); su.push_back(T[k].second);
            }
            if (su.front() != 0.0)
                throw std::runtime_error("cdp: cdpTensionTable must start at "
                                         "u_ck = 0");
            if (std::abs(ss.front() - ft_) > 1e-6 * ft_)
                throw std::runtime_error("cdp: cdpTensionTable sigma_t(0) = "
                    + std::to_string(ss.front()) + " must equal ft = "
                    + std::to_string(ft_) + " (sigma_t0 of the potential is "
                      "the Material ft)");
            ut0_ = su.back();
        }
        auto Dt = cdpParsePairs("cdpTensionDamage", tdam);     // (d_t, u)
        std::vector<double> du, dtv;
        for (std::size_t k = 0; k < Dt.size(); ++k) {
            if (!(Dt[k].first >= 0.0 && Dt[k].first <= 0.99))
                throw std::runtime_error("cdp: cdpTensionDamage point "
                    + std::to_string(k) + " : 0 <= d_t <= 0.99 required");
            if (k > 0 && !(Dt[k].second > Dt[k - 1].second))
                throw std::runtime_error("cdp: cdpTensionDamage point "
                    + std::to_string(k) + " : u_ck must increase strictly");
            dtv.push_back(Dt[k].first); du.push_back(Dt[k].second);
        }
        if (du.front() != 0.0 || dtv.front() != 0.0)
            throw std::runtime_error("cdp: cdpTensionDamage must start at "
                                     "d_t = 0 for u_ck = 0");
        std::vector<double> U = {0.0};
        U.insert(U.end(), su.begin(), su.end());
        U.insert(U.end(), du.begin(), du.end());
        std::sort(U.begin(), U.end());
        U.erase(std::unique(U.begin(), U.end(),
                            [](double a, double b) {
                                return std::abs(a - b) <= 1e-12 * std::max(a, b);
                            }), U.end());
        if (U.size() > 16)
            throw std::runtime_error("cdp: at most 16 tension nodes (got "
                                     + std::to_string(U.size()) + ")");
        dtMax_ = 0.0;
        for (double u : U) {
            double st = gfi_ ? ft_ * (1.0 - u / ut0_) : cdpInterp(su, ss, u);
            st = std::max(st, 1.0e-3 * ft_);                  // plancher NOMINAL
            double dt = std::min(cdpInterp(du, dtv, u), 0.99);
            tU_.push_back(u);
            tSig_.push_back(st);
            tD_.push_back(dt);
            tC_.push_back(dt / (1.0 - dt) * st / E0_);
            dtMax_ = std::max(dtMax_, dt);
        }

        // ---- 6. gardes a lcMax --------------------------------------------
        // (G3) pente locale de sigma_bar = sigma/(1-d), sigma et d affines
        // en eps sur un segment : sigma_bar' = [B (1-d) + D sigma]/(1-d)^2,
        // extremale a l'un des deux bouts. Seuils demontres en uniaxial,
        // facteur 0,9 de prudence.
        {
            const double mu = G_, K = K_, t = t_;
            limT_ = 0.9 * (2.0 * mu + K * t) / (1.0 + t / 3.0);
            limC_ = 0.9 * (3.0 * mu + 3.0 * alpha_ * K * t)
                    / ((1.0 - alpha_) * (1.0 - t / 3.0));
        }
        if (compLen_ <= 0.0) {
            // la compression ne depend ni de lc ni de ftScale : une fois pour toutes
            for (std::size_t k = 1; k < cEpl_.size(); ++k) {
                double sl = worstSlope(cEpl_[k - 1], cEpl_[k], cSig_[k - 1],
                                       cSig_[k], cD_[k - 1], cD_[k]);
                if (sl > limC_)
                    throw std::runtime_error("cdp (G3): effective snap-back in "
                        "compression on segment " + std::to_string(k)
                        + " : -d sigma_bar_c/d eps_c_pl = " + std::to_string(sl / 1e9)
                        + " GPa > " + std::to_string(limC_ / 1e9)
                        + " GPa — soften the hardening table or lower d_c");
            }
        } else {
            // crack band en compression : les noeuds post-pic dependent de lc
            // (pas de ftScale) ; la plus grande maille a les segments les plus
            // courts — gardes (G2c)(G3)(G4) a lcMax, voir checkCompGuards
            if (!(lcMax > 0.0))
                throw std::runtime_error("cdp: cdpCompLength > 0 requires "
                                         "lcMax > 0");
            checkCompGuards(lcMax);
        }
        // traction : (G1)(G2)(G3) a (lcMax, ftScale = 1) ; re-verifiees par
        // element a ftScale != 1 (revue 2026-09-04, voir checkGuards)
        checkGuards(lcMax, 1.0);
        if (erodeD_ > dtMax_ + 1e-12)
            std::cout << "[cdp] WARNING: erodeD = " << erodeD_ << " > d_t,max = "
                      << dtMax_ << " : canal spall INERTE (poser erodeD <= "
                      << dtMax_ << " pour qu'il agisse)\n";
        std::cout << "[cdp] psi " << psiDeg << " deg, ecc " << ecc_ << ", fb0/fc0 "
                  << fbfc_ << ", Kc " << Kc_ << " -> alpha " << alpha_
                  << ", gamma " << gamma_ << ", pente triaxiale m = " << mSlope_
                  << " ; fc0_pic " << fc0pic_ / 1e6 << " MPa, fb0 "
                  << fbfc_ * fc0pic_ / 1e6 << " MPa, ft " << ft_ / 1e6
                  << " MPa, u_t0 " << ut0_ << " m (" << (gfi_ ? "gfi" : "table")
                  << "), " << cEpl_.size() << " noeuds de compression, "
                  << tU_.size() << " de traction ; lcMax " << lcMax
                  << " m (gardes G1-G3 OK, seuils " << limT_ / 1e9 << " / "
                  << limC_ / 1e9 << " GPa)"
                  << (muV_ > 0.0 ? " ; viscosite Duvaut-Lions mu = "
                                   + std::to_string(muV_) + " s" : std::string())
                  << (compLen_ > 0.0 ? " ; crack band en compression L_ref = "
                                       + std::to_string(compLen_) + " m (pic a eps_in "
                                       + std::to_string(einPic_) + ", gardes G2c/G3/G4 OK a lcMax)"
                                     : std::string())
                  << (cap_ ? " ; cap volumique p_bar <= pc, pc0 = "
                             + std::to_string(capP0_ / 1e6) + " MPa, H = "
                             + std::to_string(capH_ / 1e9) + " GPa (pente K H/(K+H) = "
                             + std::to_string(K_ * capH_ / (K_ + capH_) / 1e9)
                             + " GPa sur p(eps_v), sans effet sur d_t/d_c)"
                           : std::string())
                  << "\n";
    }

    // (sigma_c, d_c) a eps_c_pl : lineaires entre noeuds convertis, a la
    // longueur de reference (= tables telles quelles ; c'est la lecture
    // d'origine, et celle de tout element quand cdpCompLength est absent)
    void compState(double epl, double& sig, double& d) const {
        sig = cdpInterp(cEpl_, cSig_, epl);
        d = cdpInterp(cEpl_, cD_, epl);
    }
    // noeud k de compression converti pour un element de taille lc (crack
    // band, cdpCompLength > 0) : pre-pic = noeud de l'init EXACTEMENT,
    // post-pic = eps_in,pic + (eps_in,k - eps_in,pic) L_ref/lc - c_k
    double compNode(std::size_t k, double lc) const {
        if (k <= kPic_) return cEpl_[k];
        return einPic_ + (cEin_[k] - einPic_) * (compLen_ / lc) - cC_[k];
    }
    // (sigma_c, d_c) a eps_c_pl pour un element de taille lc : a cle absente,
    // delegue a la lecture d'origine (arithmetique identique) ; sinon noeuds
    // post-pic remis a l'echelle a la volee (interpolation lineaire en
    // eps_c_pl, constante au-dela des deux bouts — regle Abaqus)
    void compState(double epl, double lc, double& sig, double& d) const {
        if (compLen_ <= 0.0) { compState(epl, sig, d); return; }
        const std::size_t n = cEpl_.size();
        if (epl <= cEpl_[0]) { sig = cSig_[0]; d = cD_[0]; return; }
        if (epl >= compNode(n - 1, lc)) { sig = cSig_[n - 1]; d = cD_[n - 1]; return; }
        std::size_t k = 1;
        double ek = compNode(1, lc);
        while (k + 1 < n && ek < epl) { ++k; ek = compNode(k, lc); }
        const double ekm = compNode(k - 1, lc);
        const double w = (epl - ekm) / (ek - ekm);
        sig = cSig_[k - 1] + w * (cSig_[k] - cSig_[k - 1]);
        d = cD_[k - 1] + w * (cD_[k] - cD_[k - 1]);
    }
    // Gardes de compression pour un element de taille lc (cdpCompLength > 0
    // seulement ; appelees a l'init a lcMax, la maille la plus defavorable) :
    // (G2c) noeuds convertis strictement croissants ; (G3) pente effective
    // <= limC (existence du retour) ; (G4) snap-back au point materiel :
    // -d sigma_bar_c/d eps_c_pl <= E0 localement (<=> d eps_total/d eps_pl
    // >= 0). Les segments pre-pic ne dependent pas de lc mais sont verifies
    // aussi (G3, G4 : un ecrouissage est toujours conforme).
    void checkCompGuards(double lc) const {
        const std::size_t n = cEpl_.size();
        const std::string where = " (cdpCompLength = " + std::to_string(compLen_)
            + " m, lc = " + std::to_string(lc) + " m)";
        for (std::size_t k = 1; k < n; ++k) {
            const double e0 = compNode(k - 1, lc), e1 = compNode(k, lc);
            if (!(e1 > e0))
                throw std::runtime_error("cdp (G2c): converted compression node "
                    + std::to_string(k) + " (eps_in = " + std::to_string(cEin_[k])
                    + " at L_ref) is not increasing in eps_c_pl once rescaled "
                      "to the element : eps_c_pl = " + std::to_string(e1)
                    + " <= previous " + std::to_string(e0)
                    + " — raise cdpCompLength, refine, or lower d_c" + where);
            const double sl = worstSlope(e0, e1, cSig_[k - 1], cSig_[k],
                                         cD_[k - 1], cD_[k]);
            // G4 (E0) est plus stricte que G3 (limC ~ 1,8 E0 sur la carte
            // Bohus) : testee d'abord, G3 reste en garde-fou
            if (sl > E0_) {
                const double einLoc = (k > kPic_)
                    ? (cEin_[k] - cEin_[k - 1]) * (compLen_ / lc)
                    : cEin_[k] - cEin_[k - 1];
                throw std::runtime_error("cdp (G4): snap-back of the compressive "
                    "softening at the material point on rescaled segment "
                    + std::to_string(k) + " : -d sigma_bar_c/d eps_c_pl = "
                    + std::to_string(sl / 1e9) + " GPa > E0 = "
                    + std::to_string(E0_ / 1e9) + " GPa (chord -Delta sigma_c/"
                      "Delta eps_in,loc = "
                    + std::to_string((cSig_[k - 1] - cSig_[k]) / einLoc / 1e9)
                    + " GPa) : the total strain would decrease while eps_c_pl "
                      "grows — refine the mesh (roughly lc < "
                    + std::to_string(lc * E0_ / sl)
                    + " m), raise cdpCompLength, or soften the table" + where);
            }
            if (sl > limC_)
                throw std::runtime_error("cdp (G3): effective snap-back in "
                    "compression on rescaled segment " + std::to_string(k)
                    + " : -d sigma_bar_c/d eps_c_pl = " + std::to_string(sl / 1e9)
                    + " GPa > " + std::to_string(limC_ / 1e9)
                    + " GPa — soften the hardening table, lower d_c, refine "
                      "or raise cdpCompLength" + where);
        }
    }
    // eps_c_pl du pic nominal (independant de lc) et du dernier noeud a lc
    double compEplPic() const { return cEpl_[kPic_]; }
    double compEplLast(double lc) const {
        return compLen_ > 0.0 ? compNode(cEpl_.size() - 1, lc) : cEpl_.back();
    }
    double compLength() const { return compLen_; }
    // noeuds de traction convertis pour (lc, s = ftScale) : e_k = kappa u_k/lc
    // - s c_k, sigma_k -> s sigma_k, d_k inchange, kappa = Gf_loc/(Gf s) (1/s en
    // portee strength, 1 en strengthGf). A s = 1 : arithmetique IDENTIQUE a la
    // version d'origine (x 1,0 et / 1,0 exacts en IEEE).
    void tenNodes(double lc, double sF, double* e, double* sg) const {
        const double kap = wScaleGf_ ? 1.0 : 1.0 / sF;
        for (std::size_t k = 0; k < tU_.size(); ++k) {
            e[k] = (tU_[k] * kap) / lc - sF * tC_[k];
            sg[k] = sF * tSig_[k];
        }
    }
    // (sigma_t, d_t) a eps_t_pl : noeuds convertis a la volee pour (lc, s)
    void tenState(double epl, double lc, double sF, double& sig, double& d) const {
        const std::size_t n = tU_.size();
        double e[16], sg[16];
        tenNodes(lc, sF, e, sg);
        if (epl <= e[0]) { sig = sg[0]; d = tD_[0]; return; }
        if (epl >= e[n - 1]) { sig = sg[n - 1]; d = tD_[n - 1]; return; }
        std::size_t k = 1;
        while (k + 1 < n && e[k] < epl) ++k;
        double w = (epl - e[k - 1]) / (e[k] - e[k - 1]);
        sig = sg[k - 1] + w * (sg[k] - sg[k - 1]);
        d = tD_[k - 1] + w * (tD_[k] - tD_[k - 1]);
    }
    void tenState(double epl, double lc, double& sig, double& d) const {
        tenState(epl, lc, 1.0, sig, d);
    }
    // pente adoucissante max de sigma/(1-d) sur un segment (voir 6.)
    static double worstSlope(double e0, double e1, double s0, double s1,
                             double d0, double d1) {
        double B = (s1 - s0) / (e1 - e0), D = (d1 - d0) / (e1 - e0);
        double g0 = (B * (1.0 - d0) + D * s0) / ((1.0 - d0) * (1.0 - d0));
        double g1 = (B * (1.0 - d1) + D * s1) / ((1.0 - d1) * (1.0 - d1));
        return std::max(-g0, -g1);
    }
    // Gardes de traction (G1)(G2)(G3) pour un element de taille lc et de
    // facteur de Weibull sF (revue 2026-09-04) : a l'init (lcMax, 1), puis au
    // premier appel de chaque element a sF != 1 — G1 et G3 se durcissent en
    // s^2 en portee strength (u_t0 en 1/s : l'element se comporte comme un
    // element de taille s^2 lc).
    void checkGuards(double lc, double sF) const {
        const double kap = wScaleGf_ ? 1.0 : 1.0 / sF;
        const double ftL = ft_ * sF, GfL = Gf_ * (wScaleGf_ ? sF : 1.0);
        const std::string where = (sF != 1.0)
            ? " (element : lc = " + std::to_string(lc) + " m, ftScale = "
                  + std::to_string(sF) + ")"
            : std::string();
        if (gfi_) {
            double lcSB = 2.0 * E0_ * GfL / (ftL * ftL);
            if (lc > lcSB)
                throw std::runtime_error("cdp (G1): largest element "
                    + std::to_string(lc) + " m exceeds the nominal snap-back "
                      "limit 2 E Gf/ft^2 = " + std::to_string(lcSB)
                    + " m — refine the mesh or raise Gf" + where);
        } else {
            for (std::size_t k = 1; k < tU_.size(); ++k) {
                double slope = -(tSig_[k] - tSig_[k - 1]) * sF
                               / ((tU_[k] - tU_[k - 1]) * kap);
                if (slope * lc > E0_)
                    throw std::runtime_error("cdp (G1): cdpTensionTable segment "
                        + std::to_string(k) + " : |d sigma_t/du| x lc > E0 "
                          "(nominal snap-back) — refine or soften the table"
                        + where);
            }
        }
        for (std::size_t k = 1; k < tU_.size(); ++k)
            if (!((tU_[k] - tU_[k - 1]) * kap > lc * sF * (tC_[k] - tC_[k - 1])))
                throw std::runtime_error("cdp (G2): converted tension node "
                    + std::to_string(k) + " (u_ck = " + std::to_string(tU_[k])
                    + ") is not increasing in eps_t_pl at lc = "
                    + std::to_string(lc) + " m" + where);
        double e[16], sg[16];
        tenNodes(lc, sF, e, sg);
        for (std::size_t k = 1; k < tU_.size(); ++k) {
            double sl = worstSlope(e[k - 1], e[k], sg[k - 1], sg[k],
                                   tD_[k - 1], tD_[k]);
            if (sl > limT_)
                throw std::runtime_error("cdp (G3): effective snap-back in "
                    "tension on segment " + std::to_string(k) + " at lc = "
                    + std::to_string(lc) + " m : -d sigma_bar_t/d eps_t_pl = "
                    + std::to_string(sl / 1e9) + " GPa > "
                    + std::to_string(limT_ / 1e9) + " GPa — raise Gf, refine, "
                      "or lower the last d_t" + where);
        }
    }
    double alpha() const { return alpha_; }
    double gamma() const { return gamma_; }
    double slopeM() const { return mSlope_; }
    double fc0pic() const { return fc0pic_; }
    double fb0() const { return fbfc_ * fc0pic_; }
    double tanPsi() const { return t_; }
    double aPot() const { return a_; }
    double ut0() const { return ut0_; }
    double Kc() const { return Kc_; }
    bool capOn() const { return cap_; }
    double capP0() const { return capP0_; }
    double capH() const { return capH_; }

    Eigen::Matrix3d stress(const Eigen::Matrix3d& eps, MatState& s,
                           double dt, double lc) const override {
        if (s.eroded) return Eigen::Matrix3d::Zero();
        const double mu = G_, K = K_, t = t_;
        // ---- A'. heterogeneite par element (revue 2026-09-04, E1 de dpr) ----
        // sF = 1 (defaut, homogene) : toutes les expressions ci-dessous se
        // reduisent EXACTEMENT a la version d'origine (x 1,0 exact)
        const double sF = s.ftScale;
        if (sF != 1.0 && !s.cdp.guarded) {
            if (!(sF > 0.0))
                throw std::runtime_error("cdp: ftScale must be > 0");
            checkGuards(lc, sF);
            s.cdp.guarded = true;
        }
        const double aLoc = a_ * sF;                       // ecc ft_loc tan psi
        const double GfLoc = wScaleGf_ ? Gf_ * sF : Gf_;
        if (muV_ > 0.0 && !(dt > 0.0))
            throw std::runtime_error("cdp: cdpViscosity > 0 requires dt > 0 "
                                     "(matpoint : mpDt = strain step / rate)");

        // ---- B. predicteur elastique, decomposition spectrale ------------
        Eigen::Matrix3d sigTr = elastic(eps - s.epsP);
        // ---- B'. cap volumique de compaction (cdpCap = true seulement) -----
        // Meme brique que le « pressure cap » de PlasticDamageLaw::stress :
        // sur la pression EFFECTIVE du predicteur p_bar = -tr(sigma_bar_tr)/3
        // (compression positive), si p_bar > pc (pc initialise a cdpCapP0 au
        // premier appel, etat partage MatState::pc) : retour volumique
        // dev = (p_bar - pc)/(K + H) > 0, eps_pl -= dev/3 I (la trace de
        // eps_pl DIMINUE : compaction), pc += H dev, sigma_bar_tr += K dev I
        // (p_bar redescend a pc + H dev). Le retour de Lubliner s'applique
        // ensuite sur le predicteur corrige. CHOIX documente : le cap ne
        // touche NI eps_t_pl / eps_c_pl NI d_t / d_c — comme dans dpr, la
        // compaction a son propre ecrouissage (pc) et n'alimente pas
        // l'endommagement de la table (l'ecrasement sous l'insert est une
        // densification, pas une fissuration de Lee-Fenves) ; elle est
        // comptee dans wPlas (sigma_nom : d eps_pl) comme tout increment
        // plastique. Sans la cle : aucune instruction executee, sigTr intact.
        // REVUE (2026-09-04 nuit) : le cap est impose au PREDICTEUR seulement.
        // Contrairement a dpr (retour du cone purement deviatorique, p intact,
        // cap exact en fin de pas), le retour de Lubliner est DILATANT
        // (p = p_tr + K tan(psi) lambda, signe +) : quand cap ET cone sont
        // actifs dans le meme pas, p_bar en fin de pas = pc + K tan(psi) lambda
        // = pc + K (d tr eps_pl + d pc/H) > pc. Depassement borne, recappe au
        // predicteur suivant (lag d'un pas), verifie pas a pas par le banc (o3)
        // (sigma3 = 100 MPa, cap 100 : cap et cone actifs ensemble). Le seuil
        // NOMINAL du cap est (1 - d) pc : a d_c = 0,92, 0,08 pc0.
        Eigen::Matrix3d dEpsCap = Eigen::Matrix3d::Zero();
        double capDv = 0.0;
        if (cap_) {
            if (s.pc == 0.0) s.pc = capP0_;                  // lazy init
            const double pBarTr = -sigTr.trace() / 3.0;
            if (pBarTr > s.pc) {
                const double dev = (pBarTr - s.pc) / (K + capH_);
                dEpsCap = -(dev / 3.0) * Eigen::Matrix3d::Identity();
                s.epsP += dEpsCap;
                s.pc += capH_ * dev;
                sigTr += K * dev * Eigen::Matrix3d::Identity();
                capDv = dev;
            }
        }
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(sigTr);
        const Eigen::Vector3d shTr = es.eigenvalues();          // croissant
        const Eigen::Matrix3d N = es.eigenvectors();
        const double pTr = -(shTr(0) + shTr(1) + shTr(2)) / 3.0;
        double sdTr[3];
        for (int i = 0; i < 3; ++i) sdTr[i] = shTr(i) + pTr;
        const double qTr = std::sqrt(1.5 * (sdTr[0] * sdTr[0] + sdTr[1] * sdTr[1]
                                            + sdTr[2] * sdTr[2]));
        const double qTol = 1e-12 * std::max(fc0pic_, std::abs(3.0 * pTr));
        const double tolF = 1e-13 * fc0pic_;

        // ---- C. etat n ------------------------------------------------------
        double sigCn, dCn, sigTn, dTn;
        compState(s.cdp.epsCpl, lc, sigCn, dCn);
        tenState(s.cdp.epsTpl, lc, sF, sigTn, dTn);
        const double sbcN = sigCn / (1.0 - dCn), sbtN = sigTn / (1.0 - dTn);
        const double betaN = sbcN / sbtN * (1.0 - alpha_) - (1.0 + alpha_);
        double sh[3] = {shTr(0), shTr(1), shTr(2)};
        const double Ftr = phi(sh, betaN, sbcN);

        Eigen::Matrix3d sigBar = sigTr;
        Eigen::Matrix3d dEpsP = Eigen::Matrix3d::Zero();
        if (Ftr > tolF) {
            // ---- E. R(lambda) : tout est fonction du multiplicateur --------
            struct Ret { double sh[3], dEh[3], epsT, epsC, F, r; };
            auto evalR = [&](double lam, Ret& o) {
                double p = pTr + K * t * lam;                    // SIGNE +
                double q = 0.0, Rq = aLoc, sd[3] = {0.0, 0.0, 0.0};
                if (qTr > qTol) {
                    q = solveQ(qTr, lam, aLoc);
                    Rq = std::sqrt(aLoc * aLoc + q * q);
                    for (int i = 0; i < 3; ++i) sd[i] = sdTr[i] * (q / qTr);
                }
                for (int i = 0; i < 3; ++i) {
                    o.sh[i] = sd[i] - p;
                    o.dEh[i] = lam * (1.5 * sd[i] / Rq + t / 3.0);
                }
                o.r = rFactor(o.sh);
                o.epsT = s.cdp.epsTpl + std::max(0.0, o.r * o.dEh[2]);
                o.epsC = s.cdp.epsCpl + std::max(0.0, -(1.0 - o.r) * o.dEh[0]);
                double sc, dc, st, dtt;
                compState(o.epsC, lc, sc, dc);
                tenState(o.epsT, lc, sF, st, dtt);
                double sbc = sc / (1.0 - dc), sbt = st / (1.0 - dtt);
                double beta = sbc / sbt * (1.0 - alpha_) - (1.0 + alpha_);
                o.F = phi(o.sh, beta, sbc);
            };
            // ---- F. crochet par doublement, puis regula falsi Illinois -----
            double lo = 0.0, Flo = Ftr;
            double hi = (1.0 - alpha_) * Ftr / (3.0 * mu);
            Ret Rh;
            evalR(hi, Rh);
            int nd = 0;
            while (Rh.F > 0.0) {
                lo = hi; Flo = Rh.F; hi *= 2.0;
                evalR(hi, Rh);
                if (++nd > 60)
                    throw std::runtime_error("cdp: no bracket for the return "
                        "mapping (F_tr = " + std::to_string(Ftr) + " Pa, p_tr = "
                        + std::to_string(pTr) + ", q_tr = " + std::to_string(qTr)
                        + ", lambda = " + std::to_string(hi) + ")");
            }
            double Fhi = Rh.F;
            Ret cur = Rh;
            int side = 0, it = 0;
            for (;; ++it) {
                if (std::abs(cur.F) <= tolF
                    || hi - lo <= 1e-15 * std::max(hi, 1e-300)) break;
                if (it >= 200)
                    throw std::runtime_error("cdp: return mapping failed after "
                        "200 iterations (F/fc0 = " + std::to_string(cur.F / fc0pic_)
                        + ", lambda in [" + std::to_string(lo) + ", "
                        + std::to_string(hi) + "], p_tr = " + std::to_string(pTr)
                        + ", q_tr = " + std::to_string(qTr) + ")");
                double lam = (it % 4 == 3) ? 0.5 * (lo + hi)   // bissection de garde
                                           : (lo * Fhi - hi * Flo) / (Fhi - Flo);
                if (!(lam > lo && lam < hi)) lam = 0.5 * (lo + hi);
                evalR(lam, cur);
                if (cur.F > 0.0) {
                    lo = lam; Flo = cur.F;
                    if (side == -1) Fhi *= 0.5;           // Illinois
                    side = -1;
                } else {
                    hi = lam; Fhi = cur.F;
                    if (side == +1) Flo *= 0.5;
                    side = +1;
                }
            }
            // ---- G. reconstruction dans la base propre du predicteur -------
            sigBar.setZero();
            for (int i = 0; i < 3; ++i) {
                sh[i] = cur.sh[i];
                Eigen::Matrix3d NN = N.col(i) * N.col(i).transpose();
                sigBar += cur.sh[i] * NN;
                dEpsP += cur.dEh[i] * NN;
            }
            s.epsP += dEpsP;
            s.cdp.epsTpl = cur.epsT;
            s.cdp.epsCpl = cur.epsC;
        }

        // ---- H. endommagement (fonction de r a l'etat final) ---------------
        double sigC, dC, sigT, dT;
        compState(s.cdp.epsCpl, lc, sigC, dC);
        tenState(s.cdp.epsTpl, lc, sF, sigT, dT);
        const double r = rFactor(sh);
        const double stf = 1.0 - wt_ * r, scf = 1.0 - wc_ * (1.0 - r);
        double d = 1.0 - (1.0 - stf * dC) * (1.0 - scf * dT);
        d = std::min(std::max(d, 0.0), 0.99);
        Eigen::Matrix3d nominal = (1.0 - d) * sigBar;
        // ---- H'. viscosite de Duvaut-Lions (cdpViscosity > 0 seulement) ----
        // eps_pl,v et d_v relaxent vers les valeurs inviscides avec le temps
        // mu (Euler implicite : x_v += dt/(mu+dt) (x - x_v)) ; la contrainte
        // nominale est celle de l'etat visqueux (CDP d'Abaqus/Standard, Theory
        // 4.5.2 ; Abaqus/Explicit ignore mu : les decks de reference sont inviscides)
        if (muV_ > 0.0) {
            const double f = dt / (muV_ + dt);
            s.cdp.epsPv += f * (s.epsP - s.cdp.epsPv);
            s.cdp.dv += f * (d - s.cdp.dv);
            nominal = (1.0 - s.cdp.dv) * elastic(eps - s.cdp.epsPv);
        }

        // ---- I. compteurs (sans effet sur la contrainte) -------------------
        Eigen::Matrix3d sigTp = Eigen::Matrix3d::Zero();
        for (int i = 0; i < 3; ++i)
            if (sh[i] > 0.0) sigTp += sh[i] * N.col(i) * N.col(i).transpose();
        Eigen::Matrix3d sigCp = sigBar - sigTp;
        s.wDamT += damageForce(sigTp) * std::max(0.0, dT - dTn);
        s.wDamC += damageForce(sigCp) * std::max(0.0, dC - dCn);
        s.wPlas += nominal.cwiseProduct(dEpsP).sum();
        // compaction du cap (terme separe : la ligne ci-dessus reste bit-
        // identique a cle absente) — -p_nom tr(d eps_pl) > 0 en compression
        if (capDv > 0.0) s.wPlas += nominal.cwiseProduct(dEpsCap).sum();

        // ---- J. champs partages -------------------------------------------
        s.D = dT;
        s.Dc = dC;
        s.epvEq = s.cdp.epsCpl;
        s.kappa = s.cdp.epsTpl;
        s.cdp.d = d;

        // ---- K. erosion (memes canaux que dpr) ----------------------------
        bool spall = s.D >= erodeD_ && nominal.trace() >= 0.0;
        if (!spall && erodeWfrac_ > 0.0 && nominal.trace() >= 0.0
            && s.wDamT >= erodeWfrac_ * GfLoc / lc)
            spall = true;                               // inerte en CDP (doc)
        bool crush = erodeEpv_ > 0.0 && s.cdp.epsCpl >= erodeEpv_;
        bool crushDc = erodeDc_ > 0.0 && dcMax_ > 0.0
                       && dC / dcMax_ >= erodeDc_;
        if (spall || crush || crushDc) {
            s.eroded = true;
            s.eroCode = spall ? 1 : (crush ? 2 : 3);
            return Eigen::Matrix3d::Zero();
        }
        return nominal;
    }

    std::string name() const override { return "cdp"; }

private:
    // Phi(sh, beta, sigma_bar_c) — sh croissant, sh[2] = s_max
    double phi(const double sh[3], double beta, double sbc) const {
        double p = -(sh[0] + sh[1] + sh[2]) / 3.0;
        double q = std::sqrt(1.5 * ((sh[0] + p) * (sh[0] + p) + (sh[1] + p) * (sh[1] + p)
                                    + (sh[2] + p) * (sh[2] + p)));
        double smax = sh[2];
        return (q - 3.0 * alpha_ * p + beta * std::max(smax, 0.0)
                - gamma_ * std::max(-smax, 0.0)) / (1.0 - alpha_) - sbc;
    }
    // r = sum <s_i> / sum |s_i| (0 en compression pure, 1 en traction pure)
    double rFactor(const double sh[3]) const {
        double num = 0.0, den = 0.0;
        for (int i = 0; i < 3; ++i) { num += std::max(sh[i], 0.0); den += std::abs(sh[i]); }
        return den > 1e-12 * fc0pic_ ? num / den : 0.0;
    }
    // q(lambda) : g(q) = q (1 + 3 mu lambda/sqrt(a^2 + q^2)) - q_tr, croissante
    // et concave ; depuis q0 = max(q_tr - 3 mu lambda, 0) (g(q0) <= 0, borne
    // inferieure exacte) le Newton converge de facon monotone par la gauche ;
    // crochet [q0, q_tr] maintenu par le signe de g, bissection si sortie
    double solveQ(double qTr, double lam, double a) const {
        const double c3 = 3.0 * G_ * lam;
        double lo = std::max(qTr - c3, 0.0), hi = qTr, q = lo;
        for (int it = 0; it < 50; ++it) {
            double Rq = std::sqrt(a * a + q * q);
            double g = q * (1.0 + c3 / Rq) - qTr;
            if (std::abs(g) <= 1e-13 * qTr) break;
            if (g > 0.0) hi = q; else lo = q;
            double gp = 1.0 + c3 * a * a / (Rq * Rq * Rq);
            double qn = q - g / gp;
            if (!(qn > lo && qn < hi)) qn = 0.5 * (lo + hi);
            if (std::abs(qn - q) <= 1e-16 * qTr) { q = qn; break; }
            q = qn;
        }
        return q;
    }
    // force thermodynamique Y = 1/2 sigma : C^-1 : sigma (isotrope) d'une
    // partie spectrale — meme formule que dpr (proxy herite, pas une identite)
    double damageForce(const Eigen::Matrix3d& sg) const {
        double tr = sg.trace();
        return 0.5 * ((1.0 + mat_.nu) / mat_.E * sg.squaredNorm()
                      - mat_.nu / mat_.E * tr * tr);
    }

    double erodeD_, erodeEpv_, erodeDc_, erodeWfrac_;
    double ecc_ = 0.1, fbfc_ = 1.16, Kc_ = 0.667, wt_ = 0.0, wc_ = 1.0;
    double ft_ = 0.0, Gf_ = 0.0, E0_ = 0.0;
    double t_ = 0.0, alpha_ = 0.0, gamma_ = 0.0, a_ = 0.0, mSlope_ = 0.0;
    double fc0pic_ = 0.0, dcMax_ = 0.0, dtMax_ = 0.0, ut0_ = 0.0;
    double limT_ = 0.0, limC_ = 0.0;
    double muV_ = 0.0;                                // cdpViscosity [s]
    double compLen_ = 0.0;                            // cdpCompLength [m], 0 = off
    bool cap_ = false;                                // cdpCap (cap volumique)
    double capP0_ = 0.0, capH_ = 0.0;                 // cdpCapP0 [Pa], cdpCapH [Pa]
    std::size_t kPic_ = 0;                            // noeud du pic nominal
    double einPic_ = 0.0;                             // eps_in du pic
    bool gfi_ = true;
    std::vector<double> cEpl_, cSig_, cD_;            // compression convertie
    std::vector<double> cEin_, cC_;                   // eps_in fusionne, c_k (crack band)
    std::vector<double> tU_, tSig_, tD_, tC_;         // traction (u_ck natif)
};

} // namespace

// ---------------------------------------------------------------------------
// Autotest POINT MATERIEL de la loi mc : on charge un point en compression
// uniaxiale, en traction uniaxiale et en compression triaxiale a plusieurs
// confinements, et on compare le plateau plastique aux formules exactes de
// Mohr-Coulomb. Aucun ajustement n'est possible — soit le retour est juste,
// soit il ne l'est pas.
//   sc = 2 c cos/(1 - sin)   st = 2 c cos/(1 + sin)
//   s1 = s3 (1 + sin)/(1 - sin) + sc      (compression triaxiale)
// ---------------------------------------------------------------------------
int mcSelftest(const std::string& csvPath) {
    Material m;
    m.E = 52e9; m.nu = 0.25; m.rho = 2620.0;
    m.cohesion = 25e6; m.phiDeg = 40.0; m.ft = 10e6; m.Gf = 70.0;
    Config cfg;
    auto law = MatLaw::make("mc", m, cfg, 1e-3);

    const double sphi = std::sin(m.phiDeg * M_PI / 180.0);
    const double cphi = std::cos(m.phiDeg * M_PI / 180.0);
    const double sc = 2.0 * m.cohesion * cphi / (1.0 - sphi);
    const double st = 2.0 * m.cohesion * cphi / (1.0 + sphi);
    const double N = (1.0 + sphi) / (1.0 - sphi);
    const double lam = m.E * m.nu / ((1.0 + m.nu) * (1.0 - 2.0 * m.nu));
    const double G = m.G();

    std::ofstream out(csvPath);
    out << "case,sigma3_MPa,plateau_MPa,exact_MPa,err_pct\n";
    double worst = 0.0;

    auto run = [&](const char* tag, double s3, double sgn) {
        // pilotage en deformation axiale, confinement lateral maintenu par
        // correction elastique iterative (etat homogene : un seul point)
        MatState st_;
        Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
        double e11 = 0.0, exact = 0.0;
        for (int k = 0; k < 4000; ++k) {
            e11 += sgn * 2.5e-6;
            eps(0, 0) = e11;
            // deformations laterales telles que s22 = s33 = s3 (elastique
            // isotrope) : e22 = e33 = (s3 (1 - nu) ... ) — resolu par
            // iterations de point fixe sur l'etat courant
            for (int it = 0; it < 40; ++it) {
                sig = law->stress(eps, st_, 1.0, 1e-3);
                double err = 0.5 * (sig(1, 1) + sig(2, 2)) - s3;
                if (std::abs(err) < 1e-3) break;
                double dlat = -err / (2.0 * (lam + G));
                eps(1, 1) += dlat; eps(2, 2) += dlat;
            }
        }
        sig = law->stress(eps, st_, 1.0, 1e-3);
        double plateau = sig(0, 0);
        // traction : s1 = st ; compression (s3 <= 0) : s1 = N s3 - sc
        exact = (sgn > 0) ? st : (N * s3 - sc);
        double err = 100.0 * (plateau - exact) / std::abs(exact);
        worst = std::max(worst, std::abs(err));
        out << tag << "," << s3 / 1e6 << "," << plateau / 1e6 << ","
            << exact / 1e6 << "," << err << "\n";
        std::cout << "[mc] " << tag << " s3 = " << s3 / 1e6
                  << " MPa : plateau " << plateau / 1e6 << " MPa, exact "
                  << exact / 1e6 << " MPa, ecart " << err << " %\n";
    };

    run("traction", 0.0, +1.0);
    run("compression", 0.0, -1.0);
    run("triaxial", -10e6, -1.0);
    run("triaxial", -25e6, -1.0);
    run("triaxial", -50e6, -1.0);
    std::cout << "[mc] ecart max = " << worst << " % ["
              << (worst < 1.0 ? "OK" : "FAIL") << "]\n";
    return worst < 1.0 ? 0 : 1;
}

// ---------------------------------------------------------------------------
// Empreinte TRIAXIALE des briques du noyau (etude « briques constitutives »,
// 2026-09-03). Meme pilotage mixte que mcSelftest (deformation axiale imposee,
// contrainte laterale tenue par point fixe elastique), carte Red Bohus :
//   E 77,66 GPa, nu 0,29, c 22,77 MPa, phi 50,4 deg (corde UCS 126,6 -> 799,3 a
//   100 MPa), ft et Gf infinis (la traction ne s'amorce jamais).
// Cibles ANALYTIQUES (le port est juste ou ne l'est pas) :
//   lineaire : q = UCS + (N - 1) sigma3,  N = (1 + sin)/(1 - sin)
//   puissance: q = UCS + B sigma3^n  (MPa)
//   mc       : idem lineaire (aretes de MC, meme corde)
// Cibles EXPERIMENTALES (information, pas critere) : 404,8 / 599,0 / 704,0 /
// 799,3 MPa a 20 / 50 / 75 / 100 MPa ; et deux controles QUI DOIVENT rater
// les donnees (lineaire a 50 MPa ; pente 3,82 de CDP a 20 MPa).
// Puis compDamage = crackband en compression uniaxiale : l'aire adoucie
// nominale x h_e doit valoir compAc * compGIIc a 2 %.
// ---------------------------------------------------------------------------
int triaxSelftest(const std::string& csvPath) {
    Material m;
    m.E = 77.66e9; m.nu = 0.29; m.rho = 2620.0;
    m.cohesion = 22.77e6; m.phiDeg = 50.4;
    m.ft = 1.0e12; m.Gf = 1.0e12;                     // traction inerte
    const double sphi = std::sin(m.phiDeg * M_PI / 180.0);
    const double cphi = std::cos(m.phiDeg * M_PI / 180.0);
    const double ucs = 2.0 * m.cohesion * cphi / (1.0 - sphi);
    const double N = (1.0 + sphi) / (1.0 - sphi);
    const double lam = m.E * m.nu / ((1.0 + m.nu) * (1.0 - 2.0 * m.nu));
    const double G = m.G();
    const double lc = 1.0e-3;

    struct Case { const char* tag; std::string cfgText; double B, n; bool power; };
    std::vector<Case> cases = {
        {"dpr_lineaire", "", 0.0, 1.0, false},
        {"dpr_puissance", "meridian = power\nmerB = 56.59\nmerN = 0.538\n", 56.59, 0.538, true},
        {"dpr_puissance_n1_corde", "meridian = power\nmerB = 6.727\nmerN = 1\n", 6.727, 1.0, true},
        {"dpr_puissance_pente382", "meridian = power\nmerB = 3.82\nmerN = 1\n", 3.82, 1.0, true},
        {"mc", "", 0.0, 1.0, false},
    };
    const double s3list[] = {0.0, 20e6, 50e6, 75e6, 100e6, 300e6, 500e6};
    const double qexp[] = {126.6e6, 404.8e6, 599.0e6, 704.0e6, 799.3e6, 0.0, 0.0};

    std::ofstream out(csvPath);
    out << "law,sigma3_MPa,q_plateau_MPa,q_formula_MPa,err_formula_pct,"
           "q_exp_MPa,err_exp_pct\n";
    bool ok = true;
    double errLin50 = 0.0, err382at20 = 0.0;

    auto writeCfg = [&](const std::string& text) {
        std::string p = csvPath + ".tmp.cfg";
        std::ofstream f(p);
        f << text;
        f.close();
        return Config::load(p);
    };

    for (const auto& cs : cases) {
        Config cfg = writeCfg(cs.cfgText);
        std::string kind = std::string(cs.tag) == "mc" ? "mc" : "dpr";
        auto law = MatLaw::make(kind, m, cfg, lc);
        for (int k = 0; k < 7; ++k) {
            double s3 = s3list[k];
            MatState st_;
            Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
            double e11 = 0.0;
            // 6 % de deformation axiale : la corde a 500 MPa plafonne a 3,5 GPa
            for (int it = 0; it < 24000; ++it) {
                e11 -= 2.5e-6;
                eps(0, 0) = e11;
                for (int j = 0; j < 40; ++j) {
                    sig = law->stress(eps, st_, 1.0, lc);
                    double err = 0.5 * (sig(1, 1) + sig(2, 2)) + s3;
                    if (std::abs(err) < 1e-3) break;
                    double dlat = -err / (2.0 * (lam + G));
                    eps(1, 1) += dlat; eps(2, 2) += dlat;
                }
            }
            sig = law->stress(eps, st_, 1.0, lc);
            double q = -sig(0, 0) - s3;
            double qf = cs.power ? ucs + 1.0e6 * cs.B * std::pow(s3 / 1.0e6, cs.n)
                                 : ucs + (N - 1.0) * s3;
            double ef = 100.0 * (q - qf) / qf;
            double ee = qexp[k] > 0 ? 100.0 * (q - qexp[k]) / qexp[k] : 0.0;
            out << cs.tag << "," << s3 / 1e6 << "," << q / 1e6 << ","
                << qf / 1e6 << "," << ef << "," << qexp[k] / 1e6 << ","
                << ee << "\n";
            std::cout << "[triax] " << cs.tag << " s3 = " << s3 / 1e6
                      << " MPa : q = " << q / 1e6 << " MPa, formule "
                      << qf / 1e6 << " (" << ef << " %)";
            if (qexp[k] > 0) std::cout << ", exp " << qexp[k] / 1e6
                                       << " (" << ee << " %)";
            std::cout << "\n";
            // la loi puissance a une tangente VERTICALE en sigma3 = 0 : le
            // residu lateral du pilotage (kPa) y vaut ~1,5 % sur q — tolerance
            // 2 % a ce seul point, 1 % ailleurs
            double tol = cs.power && cs.n < 1.0 ? (k == 0 ? 2.0 : 1.0) : 0.1;
            if (std::abs(ef) > tol) ok = false;
            if (std::string(cs.tag) == "dpr_lineaire" && k == 2) errLin50 = ee;
            if (std::string(cs.tag) == "dpr_puissance_pente382" && k == 1)
                err382at20 = ee;
        }
    }
    std::cout << "[triax] controles qui DOIVENT rater les donnees : lineaire a "
                 "50 MPa " << errLin50 << " % (attendu ~ -23) ["
              << (errLin50 < -15.0 ? "OK" : "FAIL") << "], pente 3,82 a 20 MPa "
              << err382at20 << " % (attendu ~ -50) ["
              << (err382at20 < -40.0 ? "OK" : "FAIL") << "]\n";
    if (!(errLin50 < -15.0 && err382at20 < -40.0)) ok = false;

    // ---- compDamage = crackband : aire adoucie x h_e = Ac GIIc --------------
    {
        double Ac = 0.98, GIIc = 1.0e4;
        Config cfg = writeCfg("compDamage = crackband\ncompAc = 0.98\n"
                              "compGIIc = 1e4\n");
        auto law = MatLaw::make("dpr", m, cfg, lc);
        MatState st_;
        Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
        double e11 = 0.0, W = 0.0, sPrev = 0.0;
        for (int it = 0; it < 30000; ++it) {           // 60 % axial
            e11 -= 2.0e-5;
            eps(0, 0) = e11;
            for (int j = 0; j < 60; ++j) {
                sig = law->stress(eps, st_, 1.0, lc);
                double err = 0.5 * (sig(1, 1) + sig(2, 2));
                if (std::abs(err) < 1e-3) break;
                double dlat = -err / (2.0 * (lam + G));
                eps(1, 1) += dlat; eps(2, 2) += dlat;
            }
            double sN = -sig(0, 0);
            W += 0.5 * (sN + sPrev) * 2.0e-5;          // aire nominale
            sPrev = sN;
            if (it % 5000 == 0 || it == 100 || it == 1000)
                std::cout << "[triax]   compDamage trace it " << it
                          << " : eps " << e11 << ", sig_nom " << sN / 1e6
                          << " MPa, omega_c " << st_.Dc << ", epvEq "
                          << st_.epvEq << ", eroded " << st_.eroded << "\n";
        }
        double epsP = -st_.epsP(0, 0);
        double sEnd = -sig(0, 0);
        double Wsoft = W - (1.0 - Ac) * ucs * epsP - 0.5 * sEnd * sEnd / m.E;
        double bc = ucs * lc / GIIc;
        double expect = Ac * GIIc / lc * (1.0 - std::exp(-bc * epsP));
        double err = 100.0 * (Wsoft - expect) / expect;
        double wPl = st_.wPlas - (1.0 - Ac) * ucs * epsP;
        std::cout << "[triax] compDamage : omega_c fin " << st_.Dc
                  << ", eps_p " << epsP << ", aire adoucie x h = "
                  << Wsoft * lc << " J/m^2 (attendu Ac GIIc = " << expect * lc
                  << ", ecart " << err << " %) [" << (std::abs(err) < 2.0 ? "OK" : "FAIL")
                  << "] ; compteur wPlas adouci x h = " << wPl * lc
                  << ", wDamC x h = " << st_.wDamC * lc << "\n";
        out << "compDamage,0," << Wsoft * lc << "," << expect * lc << ","
            << err << ",0,0\n";
        if (std::abs(err) > 2.0) ok = false;
    }
    std::cout << "[triax] " << (ok ? "OK" : "FAIL") << "\n";
    return ok ? 0 : 1;
}

// ---------------------------------------------------------------------------
// Pilotage MIXTE 1D pour les bancs point-materiel (cdpSelftest, matpoint) :
// x = la deformation libre (une composante, ou la paire laterale), residu =
// contrainte a tenir - cible. Premier pas par la compliance ENDOMMAGEE
// (stiff x max(1 - d, 0,02) : la compliance elastique sous-relaxe quand
// (1 - d) est petit, le point fixe de triaxSelftest converge alors comme
// d^k — 0,92^40 = 0,035, jamais a 1e-3 Pa), puis SECANTE sur les deux
// derniers essais, repli compliance si la secante s'egare (residu non
// monotone aux coins elastique/plastique). L'etat est RESTAURE a chaque
// essai : l'iteration porte sur le pas incremental, pas sur une suite
// d'appels qui accumuleraient de la plasticite. triaxSelftest reste
// INTOUCHE (croissance par addition).
// ---------------------------------------------------------------------------
namespace {

// Revue 2026-09-04 : la compliance endommagee + secante DIVERGEAIENT quand la
// raideur nominale laterale saute au passage sigma_lat = 0 (r bascule, (1-d)
// passe de 0,05 a 1) — pas larges en traction post-pic, lc tres petit, psi
// pres de la borne : lignes non convergees ecrites sans marquage. Une premiere
// correction (crochet par le signe + bissection « apres deux essais sans
// progres ») ne changeait RIEN aux trois cas (3 / 653 / 2947 pas non
// converges, memes residus) : une secante qui progresse d'un iota remet le
// compteur a zero, et la recherche du signe ne doublait que dans la direction
// de la compliance — fausse quand le residu n'est PAS monotone en x. Sonde
// (scratch, balayage de 2000 points a etat restaure) : dans les trois cas le
// residu a UNE racine franche (saut < 1e-8 Pa au croisement) mais n'est pas
// monotone en traction post-pic (1150 montees / 850 descentes sur le
// balayage) ; il l'est en equibiaxial. D'ou l'algorithme actuel, en deux
// phases et sans hypothese de monotonie :
//   A. CROCHET : premier pas par la compliance endommagee, puis recherche
//      SYMETRIQUE en doublement (x0 +- |dx| 2^k, direction de la compliance
//      d'abord) jusqu'au changement de signe (|dx| 2^k <= 2, ~30 niveaux) ;
//   B. REGULA FALSI d'Illinois dans le crochet, bissection de garde une
//      iteration sur trois, arret a |residu| < tol ou crochet reduit a
//      1e-16 relatif (discontinuite : echec honnete).
// Le point rendu est celui de plus petit |residu| rencontre, et l'etat est
// celui de CE point (dernier appel a la loi). `last` = residu final du pas
// (colonne `residu` de matpoint).
struct Hold1D {
    double tol = 1e-3;              // Pa
    int maxIt = 400;
    double worst = 0.0;             // pire residu final rencontre
    double last = 0.0;              // residu final du dernier solve
    long calls = 0, fails = 0;
    template <class SetX, class Res>
    bool solve(const MatLaw& law, MatState& s, Eigen::Matrix3d& eps, double lc,
               double dt, double& x, double stiff, SetX setX, Res res,
               Eigen::Matrix3d& sig) {
        const MatState s0 = s;
        int it = 0;
        auto eval = [&](double xx) {
            s = s0;
            setX(eps, xx);
            sig = law.stress(eps, s, dt, lc);
            ++calls; ++it;
            return res(sig);
        };
        double xBest = x, rBest = eval(x);
        if (std::abs(rBest) < tol) {
            last = std::abs(rBest); worst = std::max(worst, last);
            return true;
        }
        // ---- A. crochet [xA, xB] avec residus de signes opposes ----------
        double xA = x, rA = rBest;
        double fac = std::max(1.0 - std::max(s.cdp.d, s.Dc), 0.02);
        double dx0 = -rA / (stiff * fac);
        if (!(std::abs(dx0) > 0.0) || !std::isfinite(dx0))
            dx0 = (rA > 0.0 ? -1.0 : 1.0) * 1e-12;
        double xB = xA, rB = rA;
        bool bracketed = false;
        for (int k = 0; k < 64 && !bracketed && it < maxIt; ++k) {
            double step = dx0 * std::ldexp(1.0, k);
            if (std::abs(step) > 2.0) break;                 // 200 % de deformation
            for (int side = 0; side < 2 && !bracketed; ++side) {
                double xn = x + (side == 0 ? step : -step);
                double rn = eval(xn);
                if (std::abs(rn) < std::abs(rBest)) { xBest = xn; rBest = rn; }
                if (rn == 0.0) { bracketed = true; xA = xB = xn; rA = rB = 0.0; break; }
                if ((rn > 0.0) != (rA > 0.0)) { xB = xn; rB = rn; bracketed = true; }
            }
        }
        // ---- B. Illinois + bissection de garde dans le crochet ------------
        if (bracketed && rBest != 0.0) {
            double lo = std::min(xA, xB), hi = std::max(xA, xB);
            double flo = (xA < xB) ? rA : rB, fhi = (xA < xB) ? rB : rA;
            int side = 0, n = 0;
            while (it < maxIt && std::abs(rBest) >= tol) {
                if (hi - lo <= 1e-16 * std::max({std::abs(lo), std::abs(hi), 1e-30}))
                    break;                                   // discontinuite
                double xm = (n % 3 == 2) ? 0.5 * (lo + hi)
                                         : (lo * fhi - hi * flo) / (fhi - flo);
                if (!(xm > lo && xm < hi)) xm = 0.5 * (lo + hi);
                double rm = eval(xm);
                ++n;
                if (std::abs(rm) < std::abs(rBest)) { xBest = xm; rBest = rm; }
                if (rm == 0.0) break;
                if ((rm > 0.0) == (flo > 0.0)) {
                    lo = xm; flo = rm;
                    if (side == -1) fhi *= 0.5;
                    side = -1;
                } else {
                    hi = xm; fhi = rm;
                    if (side == +1) flo *= 0.5;
                    side = +1;
                }
            }
        }
        // ---- etat du point rendu -------------------------------------------
        x = xBest;
        double r = eval(x);
        last = std::abs(r);
        worst = std::max(worst, last);
        if (!(last < tol)) ++fails;
        return last < tol;
    }
};

// setters / residus usuels (axe pilote = z, indice 2)
inline void setLat(Eigen::Matrix3d& e, double x) { e(0, 0) = x; e(1, 1) = x; }
inline void setZ(Eigen::Matrix3d& e, double x) { e(2, 2) = x; }

} // namespace

// ---------------------------------------------------------------------------
// Banc point-materiel de la loi cdp (2026-09-04), carte historique Red Bohus
// (E 77,66 GPa, nu 0,29, ft 8,5 MPa, Gf 100 J/m^2, tables par defaut), lc =
// 1 mm. Chaque controle a un critere chiffre, verdict [OK]/[FAIL] ; code de
// retour 0 seulement si tout passe. Trace CSV : loi, chemin, sigma3, eps_ax,
// eps_lat, eps_vol, q, d_t, d_c, eps_c_pl, eps_t_pl.
//   (a) pics triaxiaux 0/20/50/75/100 MPa = fc0_pic + m sigma3 (0,2 %) et
//       identite de translation q_nom - m sigma3 = sigma_c(eps_c_pl) sur
//       TOUTE la branche (0,5 % fc0_pic) ;
//   (b) equibiaxial : pic = fb0 = fbfc fc0_pic (0,5 %) ;
//   (c) traction uniaxiale : pic = ft (0,2 %), aire totale x lc = Gf (1 %)
//       et identite exacte W lc - lc wDamT = Gf (0,1 %) a lc = 1 et 2 mm ;
//       eps_t_pl = eps_pl axiale (1e-12) ;
//   (d) decharge : pente (1-d_c) E0 en compression apres d_c > 0, (1-d_t) E0
//       en traction, E0 apres refermeture (1 %) ;
//   (e) dilatance uniaxiale -d eps_v/d eps_1 = tan psi/(rho - tan psi/3),
//       psi 35 (0,913), psi 50 (1,977), ecc 0,5 a q = 20 MPa (0,927) (1 %) ;
//   (f) DOIT RATER : carte historique -49,9 % a 20 MPa (entre -55 et -45),
//       -36,4 % a 100 (< -30) ;
//   (g) carte inverse (Kc 0,6061, fc0 316,8) : q(20) = 451,8 MPa (0,2 %) ;
//   (h) rapport des meridiens q_TE/q_TC = Kc (0,5 %) a p_bar fixe, Kc = 1 -> 1 ;
//   (i) retour : compression hydrostatique jamais plastique, traction
//       hydrostatique -> P = (1-alpha) sigma_bar_c/(3 alpha + beta) (1e-6),
//       p_bar monte de K tan psi lambda pendant l'ecoulement ;
//   (j) garde G3 : lcMax = 20 mm doit lever l'exception, 1 mm passer.
//   (k) cdpViscosity, (l) ftScale + cle inconnue, (m) crack band en compression
//       (revue et soir du 2026-09-04, voir les blocs) ;
//   (n) cap volumique cdpCap (percussion, 2026-09-04) : compression
//       hydrostatique en deformation isotrope, p = K eps_v jusqu'a 440 MPa
//       (0,01 %) puis pente K H/(K+H) = K/2 (0,5 %), pc = pc0 + H eps_v^pl,
//       d = 0 ; sans la cle K eps_v partout ; cdpCap sans cdpCapP0 = exception ;
//       cdpCapH absent = K exactement (trace bit-identique a H = K explicite) ;
//   (o) cap et cone (revue de la nuit) : o1 triaxial 20 MPa, cap 440 : p_bar
//       max mesure 256 MPa, cap JAMAIS actif, q identique a 1e-9 (invariance) ;
//       o2 20 MPa, cap 200 < 256 : identite a 1e-9 avant l'activation (post-
//       pic, d_c ~ 0,8) puis DOIT diverger (|dq| > 0,1 MPa, wPlas plus grand,
//       compaction > 0) ; o3 100 MPa, cap 100 : cap actif des q >= 0 (pas
//       ~540) puis avec le cone, q_pic invariant (0,1 %), eps_ax decale de eps_v^pl/3
//       (10 %), depassement p_bar - pc = K (d tr eps_pl + d pc/H) pas a pas.
// ---------------------------------------------------------------------------
int cdpSelftest(const std::string& csvPath) {
    Material m;
    m.E = 77.66e9; m.nu = 0.29; m.rho = 2620.0;
    m.ft = 8.5e6; m.Gf = 100.0;
    m.cohesion = 22.77e6; m.phiDeg = 50.4;             // inutilises par cdp
    const double lam = m.E * m.nu / ((1.0 + m.nu) * (1.0 - 2.0 * m.nu));
    const double G = m.G(), K = m.K(), E0 = m.E;
    const double lc = 1.0e-3;
    const double kLat = 2.0 * (lam + G), kAx = lam + 2.0 * G;

    std::ofstream out(csvPath);
    out << "loi,chemin,sigma3,eps_ax,eps_lat,eps_vol,q,d_t,d_c,eps_c_pl,eps_t_pl\n";
    out.precision(10);
    bool ok = true;
    auto check = [&](const std::string& what, double meas, double ref, double tolPct) {
        double err = 100.0 * (meas - ref) / std::max(std::abs(ref), 1e-300);
        bool p = std::abs(err) <= tolPct;
        if (!p) ok = false;
        std::cout << "[cdp] " << what << " : " << meas << " (attendu " << ref
                  << ", ecart " << err << " %, tol " << tolPct << " %) ["
                  << (p ? "OK" : "FAIL") << "]\n";
        return p;
    };
    auto verdict = [&](const std::string& what, bool p) {
        if (!p) ok = false;
        std::cout << "[cdp] " << what << " [" << (p ? "OK" : "FAIL") << "]\n";
    };
    auto mk = [&](const std::string& text, double lcMax) {
        std::string p = csvPath + ".tmp.cfg";
        std::ofstream f(p);
        f << "law = cdp\n" << text;
        f.close();
        Config cfg = Config::load(p);
        return MatLaw::make("cdp", m, cfg, lcMax);
    };
    auto row = [&](const char* tag, const char* path, double s3,
                   const Eigen::Matrix3d& eps, double q, const MatState& st) {
        out << tag << "," << path << "," << s3 << "," << eps(2, 2) << ","
            << eps(0, 0) << "," << eps.trace() << "," << q << "," << st.D << ","
            << st.Dc << "," << st.cdp.epsCpl << "," << st.cdp.epsTpl << "\n";
    };
    Hold1D hold;

    // ---- (a) triaxial, carte historique -------------------------------------
    const double s3list[5] = {0.0, 20e6, 50e6, 75e6, 100e6};
    const double qexp[5] = {126.6e6, 404.8e6, 599.0e6, 704.0e6, 799.3e6};
    double qTC[5] = {0, 0, 0, 0, 0};
    double worstTrans = 0.0;
    auto triax = [&](const MatLaw& law, const char* tag, double s3, int nSteps,
                     double de, double& qmax, double* trans) {
        const CdpLaw* cdp = dynamic_cast<const CdpLaw*>(&law);
        MatState st;
        Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
        double e = 0.0, x = 0.0;
        qmax = 0.0;
        for (int k = 0; k < nSteps; ++k) {
            const double eplPrev = st.epvEq;
            e -= de;
            hold.solve(law, st, eps, lc, 1.0, x, kLat,
                       [&](Eigen::Matrix3d& ee, double xx) { setZ(ee, e); setLat(ee, xx); },
                       [&](const Eigen::Matrix3d& sg) { return 0.5 * (sg(0, 0) + sg(1, 1)) + s3; },
                       sig);
            double q = -sig(2, 2) - s3;
            qmax = std::max(qmax, q);
            if (trans && cdp && st.epvEq > eplPrev) {   // sur la surface seulement
                double sc, dc;
                cdp->compState(st.epvEq, sc, dc);
                *trans = std::max(*trans, std::abs(q - cdp->slopeM() * s3 - sc)
                                              / cdp->fc0pic());
            }
            if (k % 25 == 0) row(tag, "triax", s3, eps, q, st);
        }
    };
    {
        auto law = mk("", lc);
        const CdpLaw* cdp = dynamic_cast<const CdpLaw*>(law.get());
        std::cout << "[cdp] constantes : alpha " << cdp->alpha() << ", gamma "
                  << cdp->gamma() << ", m " << cdp->slopeM() << ", fc0_pic "
                  << cdp->fc0pic() / 1e6 << " MPa, fb0 " << cdp->fb0() / 1e6
                  << " MPa, a " << cdp->aPot() / 1e6 << " MPa, u_t0 "
                  << cdp->ut0() << " m\n";
        check("(a) alpha", cdp->alpha(), 0.121212, 0.01);
        check("(a) gamma(Kc 0,667)", cdp->gamma(), 2.99102, 0.01);
        check("(a) m", cdp->slopeM(), 3.81737, 0.01);
        for (int k = 0; k < 5; ++k) {
            triax(*law, "hist", s3list[k], 15000, 1.0e-6, qTC[k], &worstTrans);
            check("(a) pic triaxial s3 = " + std::to_string(s3list[k] / 1e6) + " MPa",
                  qTC[k] / 1e6, (cdp->fc0pic() + cdp->slopeM() * s3list[k]) / 1e6, 0.2);
        }
        std::cout << "[cdp] (a) translation |q_nom - m s3 - sigma_c(eps_c_pl)| max = "
                  << 100.0 * worstTrans << " % fc0_pic\n";
        verdict("(a) identite de translation sur toute la branche (< 0,5 % fc0_pic)",
                worstTrans < 5e-3);
    }

    // ---- (b) equibiaxial ----------------------------------------------------
    {
        auto law = mk("", lc);
        const CdpLaw* cdp = dynamic_cast<const CdpLaw*>(law.get());
        MatState st;
        Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
        double e = 0.0, x = 0.0, smax = 0.0;
        for (int k = 0; k < 6000; ++k) {
            e -= 1.0e-6;
            hold.solve(*law, st, eps, lc, 1.0, x, kAx,
                       [&](Eigen::Matrix3d& ee, double xx) { setLat(ee, e); setZ(ee, xx); },
                       [&](const Eigen::Matrix3d& sg) { return sg(2, 2); }, sig);
            smax = std::max(smax, -sig(0, 0));
            if (k % 25 == 0) row("hist", "biaxial", 0.0, eps, -sig(0, 0), st);
        }
        check("(b) pic equibiaxial", smax / 1e6, cdp->fb0() / 1e6, 0.5);
    }

    // ---- (c) traction uniaxiale, lc = 1 et 2 mm ------------------------------
    for (double lcT : {1.0e-3, 2.0e-3}) {
        auto law = mk("", lcT);
        const CdpLaw* cdp = dynamic_cast<const CdpLaw*>(law.get());
        MatState st;
        Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
        double e = 0.0, x = 0.0, smax = 0.0, W = 0.0, sPrev = 0.0, dEpl = 0.0;
        const double de = 2.0e-6;
        int n = (int)std::ceil((cdp->ut0() / lcT + 5.0e-4) / de);
        for (int k = 0; k < n; ++k) {
            e += de;
            hold.solve(*law, st, eps, lcT, 1.0, x, kLat,
                       [&](Eigen::Matrix3d& ee, double xx) { setZ(ee, e); setLat(ee, xx); },
                       [&](const Eigen::Matrix3d& sg) { return 0.5 * (sg(0, 0) + sg(1, 1)); },
                       sig);
            double sN = sig(2, 2);
            W += 0.5 * (sN + sPrev) * de;
            sPrev = sN;
            smax = std::max(smax, sN);
            dEpl = std::max(dEpl, std::abs(st.cdp.epsTpl - st.epsP(2, 2))
                                      / std::max(st.cdp.epsTpl, 1e-30));
            if (k % 25 == 0) row("hist", "tension", 0.0, eps, sN, st);
        }
        std::string tg = " (lc = " + std::to_string(lcT * 1e3) + " mm)";
        check("(c) pic traction" + tg, smax / 1e6, m.ft / 1e6, 0.2);
        check("(c) aire totale x lc = Gf" + tg, W * lcT, m.Gf, 1.0);
        check("(c) identite W lc - lc wDamT = Gf" + tg, W * lcT - lcT * st.wDamT,
              m.Gf, 0.1);
        std::cout << "[cdp]   wDamT " << st.wDamT << " J/m^3 (x lc = " << st.wDamT * lcT
                  << " J/m^2), wPlas x lc " << st.wPlas * lcT << ", d_t fin " << st.D
                  << ", sigma fin " << sig(2, 2) << " Pa, eps fin " << e << "\n";
        // r = 1 - O(residu lateral/ft) ~ 1e-10 avec le pilotage a 1e-3 Pa :
        // l'identite est exacte a cet ordre (1e-12 supposerait sigma_lat = 0 exact)
        std::cout << "[cdp]   |eps_t_pl - eps_pl,zz| / eps_t_pl max = " << dEpl << "\n";
        verdict("(c) eps_t_pl = eps_pl axiale (relatif < 1e-8)" + tg, dEpl < 1e-8);
    }

    // ---- (d) decharges ------------------------------------------------------
    {
        auto law = mk("", lc);
        MatState st;
        Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
        double e = 0.0, x = 0.0;
        auto stepC = [&](double de) {
            e += de;
            hold.solve(*law, st, eps, lc, 1.0, x, kLat,
                       [&](Eigen::Matrix3d& ee, double xx) { setZ(ee, e); setLat(ee, xx); },
                       [&](const Eigen::Matrix3d& sg) { return 0.5 * (sg(0, 0) + sg(1, 1)); },
                       sig);
        };
        for (int k = 0; k < 20000 && st.epvEq < 3.227e-3; ++k) stepC(-1.0e-6);
        double dc = st.Dc, s0 = sig(2, 2), e0 = e;
        for (int k = 0; k < 10; ++k) stepC(+1.0e-6);
        check("(d) decharge compression, pente / E0 (d_c = " + std::to_string(dc) + ")",
              (sig(2, 2) - s0) / (e - e0) / E0, 1.0 - dc, 1.0);
        row("hist", "unload_c", 0.0, eps, -sig(2, 2), st);
    }
    {
        auto law = mk("", lc);
        const CdpLaw* cdp = dynamic_cast<const CdpLaw*>(law.get());
        MatState st;
        Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
        double e = 0.0, x = 0.0;
        auto stepT = [&](double de) {
            e += de;
            hold.solve(*law, st, eps, lc, 1.0, x, kLat,
                       [&](Eigen::Matrix3d& ee, double xx) { setZ(ee, e); setLat(ee, xx); },
                       [&](const Eigen::Matrix3d& sg) { return 0.5 * (sg(0, 0) + sg(1, 1)); },
                       sig);
        };
        double target = 0.5 * cdp->ut0() / lc;
        for (int k = 0; k < 40000 && st.cdp.epsTpl < target; ++k) stepT(+1.0e-6);
        double dt = st.D, s0 = sig(2, 2), e0 = e;
        for (int k = 0; k < 10; ++k) stepT(-1.0e-6);
        check("(d) decharge traction, pente / E0 (d_t = " + std::to_string(dt) + ")",
              (sig(2, 2) - s0) / (e - e0) / E0, 1.0 - dt, 1.0);
        for (int k = 0; k < 40000 && sig(2, 2) > -0.5e6; ++k) stepT(-1.0e-6);
        s0 = sig(2, 2); e0 = e;
        for (int k = 0; k < 10; ++k) stepT(-1.0e-6);
        check("(d) refermeture (sigma < 0 apres traction), pente / E0 (d = "
              + std::to_string(st.cdp.d) + ")", (sig(2, 2) - s0) / (e - e0) / E0, 1.0, 1.0);
        row("hist", "unload_t", 0.0, eps, sig(2, 2), st);
    }

    // ---- (e) dilatance uniaxiale ---------------------------------------------
    {
        struct V { const char* tag; std::string cfg; };
        std::vector<V> vars = {
            {"psi 35, ecc 0,1", ""},
            {"psi 50", "cdpDilationDeg = 50\n"},
            {"ecc 0,5 a q = 20 MPa", "cdpEcc = 0.5\ncdpHardening = 20e6:0 20e6:0.01\n"},
        };
        for (const V& v : vars) {
            auto law = mk(v.cfg, lc);
            const CdpLaw* cdp = dynamic_cast<const CdpLaw*>(law.get());
            MatState st;
            Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
            double e = 0.0, x = 0.0;
            int nPl = 0;
            double rv = 0.0, rl = 0.0, qb = 0.0;
            for (int k = 0; k < 4000 && nPl < 150; ++k) {
                Eigen::Matrix3d epP = st.epsP;
                e -= 1.0e-6;
                hold.solve(*law, st, eps, lc, 1.0, x, kLat,
                           [&](Eigen::Matrix3d& ee, double xx) { setZ(ee, e); setLat(ee, xx); },
                           [&](const Eigen::Matrix3d& sg) { return 0.5 * (sg(0, 0) + sg(1, 1)); },
                           sig);
                Eigen::Matrix3d dP = st.epsP - epP;
                if (dP(2, 2) < -1e-9) {
                    ++nPl;
                    rv = -dP.trace() / dP(2, 2);
                    rl = dP(0, 0) / dP(2, 2);
                    qb = -sig(2, 2) / (1.0 - st.cdp.d);
                }
            }
            double t = cdp->tanPsi(), a = cdp->aPot();
            double rho = qb / std::sqrt(a * a + qb * qb);
            double expV = t / (rho - t / 3.0);
            double expL = -(rho / 2.0 + t / 3.0) / (rho - t / 3.0);
            check(std::string("(e) dilatance -d eps_v/d eps_1, ") + v.tag + " (q_bar "
                  + std::to_string(qb / 1e6) + " MPa, rho " + std::to_string(rho) + ")",
                  rv, expV, 1.0);
            check(std::string("(e) lat/ax, ") + v.tag, rl, expL, 1.0);
            verdict(std::string("(e) identite 1 + 2 lat/ax = -(-d eps_v/d eps_1), ") + v.tag,
                    std::abs(1.0 + 2.0 * rl + rv) < 1e-3 * std::abs(rv) + 1e-12);
            if (v.cfg.empty()) {
                double sMC = rv / (2.0 + rv);
                std::cout << "[cdp]   psi_MC equivalent = " << std::asin(sMC) * 180.0 / M_PI
                          << " deg (Mohr-Coulomb a psi 35 donnerait 2 sin/(1-sin) = "
                          << 2.0 * std::sin(35 * M_PI / 180) / (1 - std::sin(35 * M_PI / 180))
                          << ")\n";
            }
        }
    }

    // ---- (f) controle qui DOIT rater les donnees --------------------------------
    {
        double e20 = 100.0 * (qTC[1] - qexp[1]) / qexp[1];
        double e100 = 100.0 * (qTC[4] - qexp[4]) / qexp[4];
        std::cout << "[cdp] (f) carte historique vs mesures : 20 MPa " << qTC[1] / 1e6
                  << " vs " << qexp[1] / 1e6 << " (" << e20 << " %), 100 MPa "
                  << qTC[4] / 1e6 << " vs " << qexp[4] / 1e6 << " (" << e100 << " %)\n";
        verdict("(f) DOIT RATER : -55 % < ecart(20) < -45 %", e20 > -55.0 && e20 < -45.0);
        verdict("(f) DOIT RATER : ecart(100) < -30 %", e100 < -30.0);
        for (int k = 0; k < 5; ++k)
            std::cout << "[cdp]   s3 " << s3list[k] / 1e6 << " : q " << qTC[k] / 1e6
                      << " MPa, exp " << qexp[k] / 1e6 << " ("
                      << 100.0 * (qTC[k] - qexp[k]) / qexp[k] << " %)\n";
    }

    // ---- (g) carte inverse de Fernando ---------------------------------------
    {
        auto law = mk("cdpKc = 0.6061\ncdpHardening = 126.6e6:0 316.8e6:0.0008 "
                      "150.2e6:0.004 25e6:0.012\n", lc);
        const CdpLaw* cdp = dynamic_cast<const CdpLaw*>(law.get());
        std::cout << "[cdp] (g) carte inverse : gamma " << cdp->gamma() << ", m "
                  << cdp->slopeM() << ", fc0_pic " << cdp->fc0pic() / 1e6 << "\n";
        double q0, q20;
        triax(*law, "inverse", 0.0, 15000, 1.0e-6, q0, nullptr);
        triax(*law, "inverse", 20e6, 15000, 1.0e-6, q20, nullptr);
        check("(g) carte inverse q(0)", q0 / 1e6, 316.8, 0.2);
        check("(g) carte inverse q(20) = 316,8 + 6,751 x 20", q20 / 1e6,
              316.8 + 6.751 * 20.0, 0.2);
        std::cout << "[cdp]   (404,8 mesure : " << 100.0 * (q20 - 404.8e6) / 404.8e6
                  << " % au point materiel ; 404,8 n'est atteint que par le "
                     "rabattement structurel r = 0,896)\n";
    }

    // ---- (h) rapport des meridiens = Kc ---------------------------------------
    auto meridian = [&](const std::string& cfg, const char* tag, double& ratio, double& Kc) {
        auto law = mk(cfg, lc);
        const CdpLaw* cdp = dynamic_cast<const CdpLaw*>(law.get());
        Kc = cdp->Kc();
        double qtc;
        triax(*law, tag, 50e6, 15000, 1.0e-6, qtc, nullptr);
        double pbar = (qtc + 3.0 * 50e6) / 3.0;
        double sa = pbar - 2.0 * Kc * qtc / 3.0;           // axiale du TE
        MatState st;
        Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
        double e = 0.0, x = 0.0, qte = 0.0;
        for (int k = 0; k < 15000; ++k) {
            e -= 1.0e-6;
            hold.solve(*law, st, eps, lc, 1.0, x, kAx,
                       [&](Eigen::Matrix3d& ee, double xx) { setLat(ee, e); setZ(ee, xx); },
                       [&](const Eigen::Matrix3d& sg) { return sg(2, 2) + sa; }, sig);
            qte = std::max(qte, -sig(0, 0) - sa);
            if (k % 25 == 0) row(tag, "TE", sa, eps, -sig(0, 0) - sa, st);
        }
        ratio = qte / qtc;
        std::cout << "[cdp] (h) " << tag << " : q_TC(50) " << qtc / 1e6 << " MPa, p_bar "
                  << pbar / 1e6 << ", sigma_a(TE) " << sa / 1e6 << " MPa ("
                  << (sa > 0 ? "compressive, regime gamma" : "TRACTIVE : hors validite")
                  << "), q_TE " << qte / 1e6 << " MPa, formule fb0 + sigma_a (3a+g)/(1-2a) = "
                  << (cdp->fb0() + sa * (3.0 * cdp->alpha() + cdp->gamma())
                      / (1.0 - 2.0 * cdp->alpha())) / 1e6 << "\n";
    };
    {
        double ratio, Kc;
        meridian("", "hist", ratio, Kc);
        check("(h) q_TE/q_TC = Kc", ratio, Kc, 0.5);
        meridian("cdpKc = 1\n", "Kc1", ratio, Kc);
        check("(h) falsifiant Kc = 1 -> rapport 1", ratio, 1.0, 0.5);
    }

    // ---- (i) proprietes du retour ------------------------------------------------
    {
        auto law = mk("", lc);
        const CdpLaw* cdp = dynamic_cast<const CdpLaw*>(law.get());
        // compression hydrostatique : jamais plastique (pas de cap)
        MatState st;
        Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
        for (int k = 0; k < 1000; ++k) {
            eps = -(k + 1) * 3.0e-6 * Eigen::Matrix3d::Identity();
            sig = law->stress(eps, st, 1.0, lc);
        }
        verdict("(i) compression hydrostatique a " + std::to_string(-sig(0, 0) / 1e6)
                + " MPa : eps_p = 0, d = 0", st.epsP.norm() == 0.0 && st.cdp.d == 0.0);
        // traction hydrostatique (q_tr = 0) : retour purement volumique,
        // sigma_hat = P = (1-alpha) sigma_bar_c/(3 alpha + beta), nominal (1-d) P
        MatState st2;
        double worstIso = 0.0, worstP = 0.0, Pn = 0.0;
        int nPl = 0;
        for (int k = 0; k < 4000; ++k) {
            eps = (k + 1) * 1.0e-7 * Eigen::Matrix3d::Identity();
            sig = law->stress(eps, st2, 1.0, lc);
            double off = std::max({std::abs(sig(0, 1)), std::abs(sig(0, 2)), std::abs(sig(1, 2)),
                                   std::abs(sig(0, 0) - sig(2, 2)), std::abs(sig(1, 1) - sig(2, 2))});
            worstIso = std::max(worstIso, off / cdp->fc0pic());
            if (st2.epsP.norm() > 0.0) {
                ++nPl;
                double sc, dc, stt, dtt;
                cdp->compState(st2.cdp.epsCpl, sc, dc);
                cdp->tenState(st2.cdp.epsTpl, lc, stt, dtt);
                double sbc = sc / (1.0 - dc), sbt = stt / (1.0 - dtt);
                double beta = sbc / sbt * (1.0 - cdp->alpha()) - (1.0 + cdp->alpha());
                double P = (1.0 - cdp->alpha()) * sbc / (3.0 * cdp->alpha() + beta);
                Pn = (1.0 - st2.cdp.d) * P;
                worstP = std::max(worstP, std::abs(sig(0, 0) - Pn) / Pn);
                if (st2.D >= 0.9) break;
            }
        }
        std::cout << "[cdp] (i) traction hydrostatique : " << nPl << " pas plastiques, "
                  << "anisotropie max " << worstIso << " fc0, ecart max a (1-d) P "
                  << worstP << " (dernier P nominal " << Pn / 1e6 << " MPa, d_t " << st2.D << ")\n";
        verdict("(i) traction hydrostatique isotrope (< 1e-9 fc0) et plastique", worstIso < 1e-9 && nPl > 0);
        verdict("(i) sigma = (1-d)(1-alpha) sigma_bar_c/(3 alpha + beta) (< 1e-6)", worstP < 1e-6);
        // signe de p_bar : compression uniaxiale, p_bar(n+1) - p_bar(tr) = K tr(d eps_p) > 0
        MatState st3;
        eps.setZero();
        double e = 0.0, x = 0.0, worstSign = 0.0;
        int nP = 0;
        bool signOk = true;
        for (int k = 0; k < 3000; ++k) {
            Eigen::Matrix3d epP = st3.epsP;
            e -= 1.0e-6;
            hold.solve(*law, st3, eps, lc, 1.0, x, kLat,
                       [&](Eigen::Matrix3d& ee, double xx) { setZ(ee, e); setLat(ee, xx); },
                       [&](const Eigen::Matrix3d& sg) { return 0.5 * (sg(0, 0) + sg(1, 1)); },
                       sig);
            Eigen::Matrix3d dP = st3.epsP - epP;
            if (dP.norm() > 0.0) {
                ++nP;
                Eigen::Matrix3d ee = eps - epP;
                double pTr = -(lam * ee.trace() * 3.0 + 2.0 * G * ee.trace()) / 3.0;
                double pNew = -sig.trace() / (3.0 * (1.0 - st3.cdp.d));
                double dp = pNew - pTr, ref = K * dP.trace();
                if (!(dp > 0.0)) signOk = false;
                worstSign = std::max(worstSign, std::abs(dp - ref) / std::abs(ref));
            }
        }
        std::cout << "[cdp] (i) compression uniaxiale : " << nP << " pas plastiques, "
                  << "p_bar(n+1) - p_bar(tr) = K tr(d eps_p) a " << worstSign << " pres\n";
        verdict("(i) p_bar MONTE pendant l'ecoulement (signe + de K tan psi lambda)",
                nP > 0 && signOk && worstSign < 1e-6);
    }

    // ---- (j) gardes G1 / G2 / G3 : le MESSAGE est inspecte (revue 2026-09-04) ---
    auto tryMk = [&](const std::string& cfg, double lcMax, std::string& msg) {
        msg.clear();
        try { mk(cfg, lcMax); } catch (const std::exception& ex) { msg = ex.what(); }
        return !msg.empty();
    };
    {
        std::string msg;
        bool threw = tryMk("", 0.02, msg);
        std::cout << "[cdp] (j) lcMax = 20 mm : " << (threw ? msg : "aucune exception") << "\n";
        verdict("(j) lcMax = 20 mm leve l'exception G3 (message '(G3)')",
                threw && msg.find("(G3)") != std::string::npos);
        threw = tryMk("", 0.3, msg);
        std::cout << "[cdp] (j) lcMax = 300 mm : " << (threw ? msg : "aucune exception") << "\n";
        verdict("(j) lcMax = 300 mm > 2 E Gf/ft^2 = 215 mm leve G1 (message '(G1)')",
                threw && msg.find("(G1)") != std::string::npos);
        // d_t = 0,9 des u = 1 um : c_1 = 9 x 0,957 ft/E0 = 9,4e-4, et
        // u_1 - u_0 = 1e-6 < lcMax c_1 = 1,9e-6 a lcMax = 2 mm : noeud converti
        // DECROISSANT (G1 passe : 2 mm << 215 mm ; G2 est testee avant G3)
        threw = tryMk("cdpTensionDamage = 0:0 0.9:1e-6 0.95:2.35e-5\n", 2.0e-3, msg);
        std::cout << "[cdp] (j) d_t 0,9 a u = 1 um, lcMax = 2 mm : "
                  << (threw ? msg : "aucune exception") << "\n";
        verdict("(j) noeud converti decroissant leve G2 (message '(G2)')",
                threw && msg.find("(G2)") != std::string::npos);
        threw = tryMk("", 1.0e-3, msg);
        verdict("(j) lcMax = 1 mm passe", !threw);
    }

    // ---- (k) viscosite de Duvaut-Lions (cdpViscosity, revue 2026-09-04) --------
    // Carte PLATE (sigma_c = 20 MPa, d = 0) en compression uniaxiale pilotee a
    // la vitesse epsdot : a l'etat stationnaire eps_pl - eps_pl,v = mu epsdot_pl
    // EXACTEMENT (Euler implicite), et la contrainte totale est
    // sigma = sigma_inv + mu D0 epsdot_pl avec sigma_inv SUR la surface. La
    // laterale totale est tenue a 0, donc la laterale inviscide vaut
    // -c, c = mu lamdot (D0 n)_lat (confinement inviscide) : l'axiale inviscide
    // est -(fc0 + (m+1) c) et Delta q = (m+1) c + mu lamdot |(D0 n)_ax|, avec
    // n = diag(rho/2 + t/3, rho/2 + t/3, -(rho - t/3)) et lamdot (rho - t/3) =
    // epsdot. Ce n'est PAS E0 mu epsdot : la pente m amplifie la laterale
    // visqueuse ~6 fois.
    {
        const std::string flat = "cdpHardening = 20e6:0 20e6:0.01\n"
                                 "cdpCompDamage = 0:0 0:0.01\n";
        auto uni = [&](const std::string& cfg, double edot, int n, double de,
                       std::vector<double>& qs, const CdpLaw** cdpOut) {
            auto law = mk(cfg, lc);
            if (cdpOut) *cdpOut = dynamic_cast<const CdpLaw*>(law.get());
            MatState st;
            Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
            double e = 0.0, x = 0.0;
            const double dtv = de / edot;
            qs.clear();
            for (int k = 0; k < n; ++k) {
                e -= de;
                hold.solve(*law, st, eps, lc, dtv, x, kLat,
                           [&](Eigen::Matrix3d& ee, double xx) { setZ(ee, e); setLat(ee, xx); },
                           [&](const Eigen::Matrix3d& sg) { return 0.5 * (sg(0, 0) + sg(1, 1)); },
                           sig);
                qs.push_back(-sig(2, 2));
                if (k % 25 == 0) row("visc", "uniaxial", edot, eps, -sig(2, 2), st);
            }
            return law;
        };
        std::vector<double> q0, q0b, q1, q2, q3;
        const double edot = 0.75, de = 1.0e-5;
        const int n = 600;
        const CdpLaw* cdp = nullptr;
        auto law0 = uni(flat, edot, n, de, q0, &cdp);
        uni(flat + "cdpViscosity = 0\n", edot, n, de, q0b, nullptr);
        bool same = q0.size() == q0b.size();
        for (std::size_t k = 0; same && k < q0.size(); ++k) same = q0[k] == q0b[k];
        verdict("(k) cdpViscosity = 0 explicite = cle absente (trace BIT-IDENTIQUE)", same);
        check("(k) carte plate inviscide : q = 20 MPa", q0.back() / 1e6, 20.0, 0.1);
        uni(flat + "cdpViscosity = 5e-5\n", edot, n, de, q1, nullptr);
        uni(flat + "cdpViscosity = 5e-5\n", 2.0 * edot, n, de, q2, nullptr);
        // cible exacte : point fixe sur c (rho depend de q_inv = fc0 + m c)
        const double mu = 5.0e-5, t = cdp->tanPsi(), a = cdp->aPot(), mS = cdp->slopeM();
        double c = 0.0, dq = 0.0;
        for (int itp = 0; itp < 20; ++itp) {
            double qi = 20e6 + mS * c;
            double rho = qi / std::sqrt(a * a + qi * qi);
            double nl = rho / 2.0 + t / 3.0, na = -(rho - t / 3.0);
            double ld = edot / (rho - t / 3.0);
            double D0l = lam * (2.0 * nl + na) + 2.0 * G * nl;
            double D0a = lam * (2.0 * nl + na) + 2.0 * G * na;
            c = mu * ld * D0l;
            dq = (mS + 1.0) * c - mu * ld * D0a;
        }
        double dq1 = q1.back() - q0.back(), dq2 = q2.back() - q0.back();
        std::cout << "[cdp] (k) mu 5e-5 s, epsdot 0,75/s : Delta q " << dq1 / 1e6
                  << " MPa (cible " << dq / 1e6 << ", E0 mu epsdot = "
                  << E0 * mu * edot / 1e6 << " MPa, confinement inviscide c = "
                  << c / 1e6 << " MPa) ; a 1,5/s : " << dq2 / 1e6 << " MPa\n";
        check("(k) surcontrainte stationnaire = (m+1) c + mu lamdot |(D0 n)_ax|",
              dq1 / 1e6, dq / 1e6, 1.0);
        check("(k) lineaire en vitesse : Delta q(2 epsdot)/Delta q(epsdot)", dq2 / dq1, 2.0, 1.0);
        verdict("(k) DOIT RATER l'estimation E0 mu epsdot : Delta q > 3 E0 mu epsdot",
                dq1 > 3.0 * E0 * mu * edot);
        uni(flat + "cdpViscosity = 1e3\n", edot, n, de, q3, nullptr);
        check("(k) falsifiant mu = 1e3 s >> T : reponse elastique q = E0 eps",
              q3.back() / 1e6, E0 * n * de / 1e6, 0.5);
        // ce que mu = 5e-5 DONNERAIT si Abaqus l'appliquait (Standard) : carte
        // historique, triaxial 20 MPa a 0,748/s — les decks de reference sont
        // Explicit et ce parametre y est ignore (le defaut inviscide est la parite)
        {
            auto run = [&](const std::string& cfg, double& qpk) {
                auto law = mk(cfg, lc);
                MatState st;
                Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
                double e = 0.0, x = 0.0;
                const double s3 = 20e6, ed = 0.748, dev = 1.0e-5, dtv = dev / ed;
                qpk = 0.0;
                for (int k = 0; k < 800; ++k) {
                    e -= dev;
                    hold.solve(*law, st, eps, lc, dtv, x, kLat,
                               [&](Eigen::Matrix3d& ee, double xx) { setZ(ee, e); setLat(ee, xx); },
                               [&](const Eigen::Matrix3d& sg) { return 0.5 * (sg(0, 0) + sg(1, 1)) + s3; },
                               sig);
                    qpk = std::max(qpk, -sig(2, 2) - s3);
                }
            };
            double qInv, qVis;
            run("", qInv);
            run("cdpViscosity = 5e-5\n", qVis);
            std::cout << "[cdp] (k) si mu 5e-5 etait actif (Standard ; Explicit l'ignore) : carte historique, s3 = 20 MPa, 0,748/s : "
                      << "pic inviscide " << qInv / 1e6 << " MPa, visqueux " << qVis / 1e6
                      << " MPa, surcontrainte " << (qVis - qInv) / 1e6 << " MPa ("
                      << 100.0 * (qVis - qInv) / qInv << " % ; sur 404,8 mesure : "
                      << 100.0 * (qVis - qInv) / 404.8e6 << " %)\n";
            verdict("(k) mu 5e-5 (Standard) : surcontrainte visqueuse > 0 au pic triaxial", qVis > qInv);
        }
    }

    // ---- (l) heterogeneite ftScale (E1) et cles inconnues (revue 2026-09-04) ----
    for (int scope = 0; scope < 2; ++scope) {
        const bool gfS = scope == 1;
        auto law = mk(gfS ? "weibullScope = strengthGf\n" : "", lc);
        const CdpLaw* cdp = dynamic_cast<const CdpLaw*>(law.get());
        const double sF = 1.3;
        const double GfLoc = gfS ? sF * m.Gf : m.Gf;
        MatState st;
        st.ftScale = sF;
        Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
        double e = 0.0, x = 0.0, smax = 0.0, W = 0.0, sPrev = 0.0;
        const double de = 2.0e-6, kap = gfS ? 1.0 : 1.0 / sF;
        int n = (int)std::ceil((cdp->ut0() * kap / lc + 5.0e-4) / de);
        for (int k = 0; k < n; ++k) {
            e += de;
            hold.solve(*law, st, eps, lc, 1.0, x, kLat,
                       [&](Eigen::Matrix3d& ee, double xx) { setZ(ee, e); setLat(ee, xx); },
                       [&](const Eigen::Matrix3d& sg) { return 0.5 * (sg(0, 0) + sg(1, 1)); },
                       sig);
            W += 0.5 * (sig(2, 2) + sPrev) * de;
            sPrev = sig(2, 2);
            smax = std::max(smax, sig(2, 2));
            if (k % 25 == 0) row(gfS ? "ftScale_strengthGf" : "ftScale_strength", "tension", sF, eps, sig(2, 2), st);
        }
        std::string tg = gfS ? " (weibullScope strengthGf)" : " (weibullScope strength)";
        check("(l) ftScale 1,3 : pic traction = 1,3 ft" + tg, smax / 1e6, sF * m.ft / 1e6, 0.2);
        check("(l) ftScale 1,3 : W lc - lc wDamT = Gf_loc" + tg, W * lc - lc * st.wDamT, GfLoc, 0.1);
        // la table de compression n'est PAS mise a l'echelle (comme la cohesion de dpr)
        MatState st2;
        st2.ftScale = sF;
        eps.setZero(); e = 0.0; x = 0.0;
        double qmax = 0.0;
        for (int k = 0; k < 4000; ++k) {
            e -= 1.0e-6;
            hold.solve(*law, st2, eps, lc, 1.0, x, kLat,
                       [&](Eigen::Matrix3d& ee, double xx) { setZ(ee, e); setLat(ee, xx); },
                       [&](const Eigen::Matrix3d& sg) { return 0.5 * (sg(0, 0) + sg(1, 1)); },
                       sig);
            qmax = std::max(qmax, -sig(2, 2));
        }
        check("(l) ftScale 1,3 : pic de compression uniaxiale = fc0_pic (table non touchee)" + tg,
              qmax / 1e6, cdp->fc0pic() / 1e6, 0.2);
    }
    {
        // garde par element : ftScale 4 (s^2 = 16 : lc < 10,7/16 = 0,67 mm) sur
        // un element de 1 mm doit lever G3 au PREMIER appel, en nommant ftScale
        auto law = mk("", lc);
        MatState st;
        st.ftScale = 4.0;
        std::string msg;
        try { law->stress(Eigen::Matrix3d::Zero(), st, 1.0, lc); }
        catch (const std::exception& ex) { msg = ex.what(); }
        std::cout << "[cdp] (l) ftScale 4 sur lc = 1 mm : " << (msg.empty() ? "aucune exception" : msg) << "\n";
        verdict("(l) garde G3 par element a ftScale 4 (message '(G3)' et 'ftScale')",
                msg.find("(G3)") != std::string::npos && msg.find("ftScale") != std::string::npos);
        bool threw = tryMk("cdpKC = 0.6\n", lc, msg);
        std::cout << "[cdp] (l) cle cdpKC : " << (threw ? msg : "aucune exception") << "\n";
        verdict("(l) cle cdp* inconnue (cdpKC) refusee en nommant la cle",
                threw && msg.find("cdpKC") != std::string::npos);
    }

    // ---- (m) crack band en compression (cdpCompLength, 2026-09-04 soir) ------
    // Compression uniaxiale, carte historique, L_ref = 2 mm, lc = 1 / 2 / 4 mm.
    // Le point materiel restitue EXACTEMENT la table (sigma_ax = sigma_c(eps_c_pl),
    // eps_ax = eps_in,loc + sigma/E0, eps_in,loc = eps_c_pl + c), donc l'aire
    // adoucie post-pic [int (sigma - sigma_res) d eps + ((sigma_pic - sigma_res)^2
    // - (sigma_end - sigma_res)^2)/(2 E0)] x lc = lc int (sigma - sigma_res)
    // d eps_in,loc = L_ref int (sigma_c - sigma_res) d eps_in = G_c : invariante
    // en lc (1 %) et egale a la quadrature de la lecture B a L_ref (1 % ; le
    // trapeze sur la table brute, 933,1 J/m^2, en differe de -0,05 %). La branche
    // PRE-pic (eps_c_pl < eps_pl,pic) est identique entre les trois lc (< 1e-9).
    // DOIT ECHOUER : sans la cle, l'aire x lc est proportionnelle a lc (rapport 4
    // entre 1 et 4 mm, a 2 %).
    {
        const double Lref = 2.0e-3, sRes = 10.0e6;
        const double lcs[3] = {1.0e-3, 2.0e-3, 4.0e-3};
        struct Trace { std::vector<double> sig, epl; double Gc, sPic, sEnd, eplPic; int nPre; };
        auto uniC = [&](const std::string& cfg, double lcE, const char* tag, Trace& T) {
            auto law = mk(cfg, lcE);
            const CdpLaw* cdp = dynamic_cast<const CdpLaw*>(law.get());
            MatState st;
            Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
            double e = 0.0, x = 0.0;
            const double de = 1.0e-6;
            const double eplEnd = cdp->compEplLast(lcE) + 2.0e-4;
            T.eplPic = cdp->compEplPic();
            T.sig.clear(); T.epl.clear();
            for (int k = 0; k < 400000 && st.cdp.epsCpl < eplEnd; ++k) {
                e -= de;
                hold.solve(*law, st, eps, lcE, 1.0, x, kLat,
                           [&](Eigen::Matrix3d& ee, double xx) { setZ(ee, e); setLat(ee, xx); },
                           [&](const Eigen::Matrix3d& sg) { return 0.5 * (sg(0, 0) + sg(1, 1)); },
                           sig);
                T.sig.push_back(-sig(2, 2));
                T.epl.push_back(st.cdp.epsCpl);
                if (k % 50 == 0) row(tag, "uniaxial", lcE, eps, -sig(2, 2), st);
            }
            std::size_t iPic = 0;
            for (std::size_t i = 1; i < T.sig.size(); ++i) if (T.sig[i] > T.sig[iPic]) iPic = i;
            T.sPic = T.sig[iPic]; T.sEnd = T.sig.back();
            double W = 0.0;
            for (std::size_t i = iPic + 1; i < T.sig.size(); ++i)
                W += 0.5 * ((T.sig[i] - sRes) + (T.sig[i - 1] - sRes)) * de;
            W += ((T.sPic - sRes) * (T.sPic - sRes) - (T.sEnd - sRes) * (T.sEnd - sRes))
                 / (2.0 * E0);
            T.Gc = W * lcE;
            T.nPre = 0;
            while (T.nPre < (int)T.epl.size() && T.epl[T.nPre] < T.eplPic) ++T.nPre;
            std::cout << "[cdp] (m) " << tag << " lc = " << lcE * 1e3 << " mm : " << T.sig.size()
                      << " pas, pic " << T.sPic / 1e6 << " MPa, sigma fin " << T.sEnd / 1e6
                      << " MPa, eps_c_pl fin " << st.cdp.epsCpl << " (dernier noeud "
                      << cdp->compEplLast(lcE) << "), " << T.nPre << " pas pre-pic, aire adoucie x lc = "
                      << T.Gc << " J/m^2\n";
        };
        // reference : quadrature de la lecture B a L_ref (dε_in = dε_pl + dc) et
        // trapeze sur la table brute (segments (116,6+50)/2 x 3,2e-3 + 50/2 x 8e-3)
        double GcB = 0.0;
        {
            auto law = mk("cdpCompLength = 0.002\n", Lref);
            const CdpLaw* cdp = dynamic_cast<const CdpLaw*>(law.get());
            const double e0 = cdp->compEplPic(), e1 = cdp->compEplLast(Lref);
            const int nq = 200000;
            double fPrev = 0.0, einPrev = 0.0;
            for (int i = 0; i <= nq; ++i) {
                double epl = e0 + (e1 - e0) * i / nq, sc, dc;
                cdp->compState(epl, Lref, sc, dc);
                double ein = epl + dc / (1.0 - dc) * sc / E0, f = sc - sRes;
                if (i > 0) GcB += 0.5 * (f + fPrev) * (ein - einPrev);
                fPrev = f; einPrev = ein;
            }
            GcB *= Lref;
        }
        const double GcRaw = Lref * (0.5 * (116.6e6 + 50.0e6) * 3.2e-3 + 0.5 * 50.0e6 * 8.0e-3);
        std::cout << "[cdp] (m) reference G_c a L_ref = 2 mm : quadrature lecture B "
                  << GcB << " J/m^2, trapeze table brute " << GcRaw << " J/m^2 ("
                  << 100.0 * (GcB - GcRaw) / GcRaw << " %)\n";
        Trace T[3];
        for (int i = 0; i < 3; ++i) uniC("cdpCompLength = 0.002\n", lcs[i], "compLen", T[i]);
        for (int i = 0; i < 3; ++i) {
            std::string tg = " (lc = " + std::to_string(lcs[i] * 1e3) + " mm)";
            check("(m) cdpCompLength 2 mm : aire adoucie post-pic x lc = G_c(lecture B a L_ref)" + tg,
                  T[i].Gc, GcB, 1.0);
            check("(m) cdpCompLength 2 mm : aire adoucie post-pic x lc = trapeze table brute" + tg,
                  T[i].Gc, GcRaw, 1.0);
            verdict("(m) sigma fin = sigma_res (table epuisee)" + tg,
                    std::abs(T[i].sEnd - sRes) < 1e-6 * sRes);
        }
        check("(m) invariance : G_c(2 mm)/G_c(1 mm)", T[1].Gc / T[0].Gc, 1.0, 1.0);
        check("(m) invariance : G_c(4 mm)/G_c(1 mm)", T[2].Gc / T[0].Gc, 1.0, 1.0);
        // branche pre-pic strictement identique (les noeuds pre-pic ne dependent pas de lc)
        double dPre = 0.0;
        int nPre = std::min({T[0].nPre, T[1].nPre, T[2].nPre});
        for (int i = 0; i < nPre; ++i)
            for (int j = 1; j < 3; ++j)
                dPre = std::max(dPre, std::abs(T[j].sig[i] - T[0].sig[i]) / T[0].sig[i]);
        std::cout << "[cdp] (m) branche pre-pic (" << nPre << " pas, eps_c_pl < " << T[0].eplPic
                  << ") : ecart relatif max entre lc = " << dPre << " ; pas pre-pic "
                  << T[0].nPre << " / " << T[1].nPre << " / " << T[2].nPre << "\n";
        verdict("(m) branche pre-pic identique entre 1 / 2 / 4 mm (< 1e-9)",
                dPre < 1e-9 && T[0].nPre == T[1].nPre && T[1].nPre == T[2].nPre && nPre > 0);
        // DOIT ECHOUER : cle absente, aire x lc proportionnelle a lc
        Trace N1, N4;
        uniC("", 1.0e-3, "noCompLen", N1);
        uniC("", 4.0e-3, "noCompLen", N4);
        std::cout << "[cdp] (m) sans la cle : aire x lc = " << N1.Gc << " (1 mm) / " << N4.Gc
                  << " J/m^2 (4 mm), rapport " << N4.Gc / N1.Gc << "\n";
        check("(m) DOIT ECHOUER sans la cle : G_c(4 mm)/G_c(1 mm) = 4 (aire x lc ~ lc)",
              N4.Gc / N1.Gc, 4.0, 2.0);
        verdict("(m) sans la cle a 1 mm : aire x lc = G_c(L_ref)/2 (table lue a lc)",
                std::abs(N1.Gc - 0.5 * GcB) < 0.01 * GcB);
        // cle absente = 0 explicite (bit-identique) et gardes : G4 (snap-back au
        // point materiel) sur une carte chutant de 126,6 a 10 MPa en 4e-4 d'eps_in
        // (corde 291 GPa) a L_ref = lc ; G2c quand la remise a l'echelle rend un
        // noeud converti decroissant (carte historique, lc = 8 mm, L_ref 1 mm)
        {
            Trace Z;
            uniC("cdpCompLength = 0\n", 1.0e-3, "compLen0", Z);
            bool same = Z.sig.size() == N1.sig.size();
            for (std::size_t i = 0; same && i < Z.sig.size(); ++i) same = Z.sig[i] == N1.sig[i];
            verdict("(m) cdpCompLength = 0 explicite = cle absente (trace BIT-IDENTIQUE)", same);
            std::string msg;
            bool threw = tryMk("cdpCompLength = 0.002\ncdpHardening = 50.6e6:0 126.6e6:0.0008 "
                               "10e6:0.0012\ncdpCompDamage = 0:0 0:0.0008 0:0.0012\n", 2.0e-3, msg);
            std::cout << "[cdp] (m) carte 126,6 -> 10 MPa en 4e-4 (corde 291 GPa), L_ref = lc = 2 mm : "
                      << (threw ? msg : "aucune exception") << "\n";
            verdict("(m) garde G4 : snap-back au point materiel leve l'exception (message '(G4)')",
                    threw && msg.find("(G4)") != std::string::npos);
            threw = tryMk("cdpCompLength = 0.002\ncdpHardening = 50.6e6:0 126.6e6:0.0008 "
                          "10e6:0.0012\ncdpCompDamage = 0:0 0:0.0008 0:0.0012\n", 0.5e-3, msg);
            verdict("(m) la meme carte passe a lc = 0,5 mm (segment etire x4 : 72,9 GPa < E0)", !threw);
            threw = tryMk("cdpCompLength = 0.001\n", 8.0e-3, msg);
            std::cout << "[cdp] (m) carte historique, L_ref 1 mm, lc = 8 mm : "
                      << (threw ? msg : "aucune exception") << "\n";
            verdict("(m) garde G2c : noeud converti decroissant apres remise a l'echelle (message '(G2c)')",
                    threw && msg.find("(G2c)") != std::string::npos);
            threw = tryMk("cdpCompLength = -1\n", lc, msg);
            verdict("(m) cdpCompLength < 0 refuse", threw);
        }
    }

    // ---- (n) cap volumique de compaction (cdpCap, 2026-09-04, percussion) -----
    // Compression hydrostatique pilotee en DEFORMATION isotrope, eps = -e I, e
    // croissant jusqu'a 3 % (eps_v = -3e, p_nom = K |eps_v| = 5,5 GPa en fin de
    // course sans cap). Avec cdpCap = true, cdpCapP0 = 440 MPa, cdpCapH = K :
    // p_nom = K |eps_v| tant que p < 440 MPa (0,01 %), puis p = (K pc0 + K H
    // |eps_v|)/(K + H) — pente K H/(K + H) = K/2 (0,5 %) ; pc = pc0 + H eps_v^pl ;
    // wPlas = int p d eps_v^pl ; d_t = d_c = 0, eps_t_pl = eps_c_pl = 0 (le cap
    // n'alimente pas l'endommagement de la table). Sans la cle : p = K |eps_v|
    // partout (0,01 %) et eps_pl = 0 ; cdpCap = false explicite = trace bit-
    // identique ; cdpCap = true sans cdpCapP0 DOIT lever une exception, et
    // cdpCapP0 sans cdpCap aussi (rien d'applique en silence).
    auto capCfgStr = [](double Kv) {
        std::ostringstream kk;
        kk.precision(17);
        kk << Kv;
        return "cdpCap = true\ncdpCapP0 = 440e6\ncdpCapH = " + kk.str() + "\n";
    };
    {
        const double pc0 = 440.0e6;
        const std::string capCfg = capCfgStr(K);
        struct HTrace { std::vector<double> ev, p, pc, evPl, w; MatState st; };
        auto hydro = [&](const std::string& cfgText, const char* tag, HTrace& T) {
            auto law = mk(cfgText, lc);
            MatState st;
            Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
            const int n = 3000;
            const double de = 0.03 / n;
            T.ev.clear(); T.p.clear(); T.pc.clear(); T.evPl.clear(); T.w.clear();
            for (int k = 0; k < n; ++k) {
                eps = -(k + 1) * de * Eigen::Matrix3d::Identity();
                sig = law->stress(eps, st, 1.0, lc);
                const double p = -sig.trace() / 3.0;
                T.ev.push_back(-eps.trace());              // |eps_v|
                T.p.push_back(p);
                T.pc.push_back(st.pc);
                T.evPl.push_back(-st.epsP.trace());       // compaction plastique > 0
                T.w.push_back(st.wPlas);
                if (k % 25 == 0) row(tag, "hydro", 0.0, eps, p, st);
            }
            T.st = st;
        };
        HTrace C, N0, F0;
        hydro(capCfg, "cap440", C);
        hydro("", "nocap", N0);
        hydro("cdpCap = false\n", "capfalse", F0);
        const CdpLaw* cdpC = nullptr;
        {
            auto lawC = mk(capCfg, lc);
            cdpC = dynamic_cast<const CdpLaw*>(lawC.get());
            std::cout << "[cdp] (n) cap : capOn " << cdpC->capOn() << ", pc0 "
                      << cdpC->capP0() / 1e6 << " MPa, H " << cdpC->capH() / 1e9
                      << " GPa (K = " << K / 1e9 << " GPa, pente attendue K H/(K+H) = "
                      << K * cdpC->capH() / (K + cdpC->capH()) / 1e9 << " GPa)\n";
            const double H = cdpC->capH();
            // branche elastique (p < pc0) et branche cappee (formule fermee)
            double errPre = 0.0, errPost = 0.0, errPc = 0.0;
            int nPre = 0, nPost = 0, kFirst = -1;
            for (std::size_t k = 0; k < C.p.size(); ++k) {
                const double pEl = K * C.ev[k];
                if (pEl <= pc0) {
                    errPre = std::max(errPre, std::abs(C.p[k] - pEl) / pEl);
                    ++nPre;
                } else {
                    if (kFirst < 0) kFirst = (int)k;
                    const double pF = (K * pc0 + H * K * C.ev[k]) / (K + H);
                    errPost = std::max(errPost, std::abs(C.p[k] - pF) / pF);
                    errPc = std::max(errPc, std::abs(C.pc[k] - (pc0 + H * C.evPl[k])) / C.pc[k]);
                    ++nPost;
                }
            }
            std::cout << "[cdp] (n) cap 440 : " << nPre << " pas elastiques (ecart max a K eps_v "
                      << 100.0 * errPre << " %), " << nPost << " pas cappes depuis |eps_v| = "
                      << (kFirst >= 0 ? C.ev[kFirst] : 0.0) << " (ecart max a la formule fermee "
                      << 100.0 * errPost << " %, pc = pc0 + H eps_v^pl a " << 100.0 * errPc
                      << " %) ; fin : p_nom " << C.p.back() / 1e6 << " MPa, pc "
                      << C.pc.back() / 1e6 << " MPa, eps_v^pl " << C.evPl.back()
                      << ", wPlas " << C.w.back() << " J/m^3, d " << C.st.cdp.d
                      << ", eps_c_pl " << C.st.cdp.epsCpl << ", eps_t_pl " << C.st.cdp.epsTpl << "\n";
            verdict("(n) cap 440 : p_nom = K eps_v avant le cap (0,01 %)", nPre > 0 && errPre < 1e-4);
            verdict("(n) cap 440 : p_nom = (K pc0 + K H eps_v)/(K+H) apres (1e-6)",
                    nPost > 0 && errPost < 1e-6);
            // pente par difference finie entre le premier pas cappe + 10 et la fin
            const int k0 = kFirst + 10, k1 = (int)C.p.size() - 1;
            const double slope = (C.p[k1] - C.p[k0]) / (C.ev[k1] - C.ev[k0]);
            check("(n) cap 440 : pente dp/d|eps_v| apres le cap / (K/2)", slope / (0.5 * K), 1.0, 0.5);
            check("(n) cap 440 : pente / (K H/(K+H))", slope / (K * H / (K + H)), 1.0, 0.5);
            verdict("(n) cap 440 : pc fin = pc0 + H eps_v^pl (1e-9)",
                    std::abs(C.pc.back() - (pc0 + H * C.evPl.back())) < 1e-9 * C.pc.back());
            verdict("(n) cap 440 : le cap n'alimente pas la table (d = 0, eps_c_pl = eps_t_pl = 0)",
                    C.st.cdp.d == 0.0 && C.st.cdp.epsCpl == 0.0 && C.st.cdp.epsTpl == 0.0
                    && C.st.D == 0.0 && C.st.Dc == 0.0);
            // wPlas = int p_nom d eps_v^pl (trapeze sur la trace)
            double Wq = 0.0;
            for (std::size_t k = 1; k < C.p.size(); ++k)
                Wq += 0.5 * (C.p[k] + C.p[k - 1]) * (C.evPl[k] - C.evPl[k - 1]);
            check("(n) cap 440 : wPlas = int p d eps_v^pl", C.w.back(), Wq, 1.0);
            // deviateur nul : la reponse reste isotrope (retour de Lubliner inactif)
            verdict("(n) cap 440 : eps_pl isotrope (hors-diagonale et deviateur nuls)",
                    std::abs(C.st.epsP(0, 1)) == 0.0 && std::abs(C.st.epsP(0, 0) - C.st.epsP(2, 2)) == 0.0);
        }
        {
            double errN = 0.0;
            for (std::size_t k = 0; k < N0.p.size(); ++k)
                errN = std::max(errN, std::abs(N0.p[k] - K * N0.ev[k]) / (K * N0.ev[k]));
            std::cout << "[cdp] (n) sans la cle : ecart max de p_nom a K eps_v " << 100.0 * errN
                      << " %, p fin " << N0.p.back() / 1e6 << " MPa, eps_pl " << N0.st.epsP.norm()
                      << ", pc " << N0.pc.back() << "\n";
            verdict("(n) sans la cle : p_nom = K eps_v partout (0,01 %), eps_pl = 0, pc = 0",
                    errN < 1e-4 && N0.st.epsP.norm() == 0.0 && N0.pc.back() == 0.0);
            bool same = F0.p.size() == N0.p.size();
            for (std::size_t k = 0; same && k < N0.p.size(); ++k)
                same = F0.p[k] == N0.p[k] && F0.w[k] == N0.w[k];
            verdict("(n) cdpCap = false explicite = cle absente (trace BIT-IDENTIQUE)", same);
            std::string msg;
            bool threw = tryMk("cdpCap = true\n", lc, msg);
            std::cout << "[cdp] (n) cdpCap = true sans cdpCapP0 : " << (threw ? msg : "aucune exception") << "\n";
            verdict("(n) DOIT ECHOUER : cdpCap = true sans cdpCapP0 leve une exception nommant cdpCapP0",
                    threw && msg.find("cdpCapP0") != std::string::npos);
            threw = tryMk("cdpCapP0 = 440e6\n", lc, msg);
            std::cout << "[cdp] (n) cdpCapP0 sans cdpCap : " << (threw ? msg : "aucune exception") << "\n";
            verdict("(n) cdpCapP0 sans cdpCap = true refuse (rien d'applique en silence)", threw);
            threw = tryMk("cdpCap = true\ncdpCapP0 = -1\n", lc, msg);
            verdict("(n) cdpCapP0 <= 0 refuse", threw);
        }
        {
            // revue : la valeur par DEFAUT de cdpCapH (cle absente -> K_) n'etait
            // verifiee par aucun verdict (capCfgStr ecrit toujours cdpCapH = K)
            const std::string defCfg = "cdpCap = true\ncdpCapP0 = 440e6\n";
            auto lawD = mk(defCfg, lc);
            const CdpLaw* cd = dynamic_cast<const CdpLaw*>(lawD.get());
            std::cout << "[cdp] (n) cdpCapH absent : H = " << cd->capH() / 1e9
                      << " GPa, K = E/(3(1-2nu)) = " << K / 1e9 << " GPa, difference "
                      << cd->capH() - K << " Pa\n";
            verdict("(n) cdpCapH absent = K = E/(3(1-2nu)) exactement", cd->capH() == K);
            HTrace D0;
            hydro(defCfg, "capHdef", D0);
            bool same = D0.p.size() == C.p.size();
            for (std::size_t k = 0; same && k < C.p.size(); ++k)
                same = D0.p[k] == C.p[k] && D0.pc[k] == C.pc[k] && D0.w[k] == C.w[k]
                       && D0.evPl[k] == C.evPl[k];
            verdict("(n) cdpCapH absent : trace hydro BIT-IDENTIQUE a cdpCapH = K explicite", same);
        }
    }

    // ---- (o) cap ET cone de Lubliner : invariance, divergence, borne ---------
    // Compression triaxiale (pilote du banc (a), 15 000 pas de 1e-6), carte
    // historique, sans cap et avec cdpCap (H = K). REVUE de la nuit : la version
    // initiale (cap 440 seul) etait vide — p_bar max MESURE 256 MPa a sigma3 =
    // 20 MPa (sigma3 effectif 20/(1-d) = 250 MPa en fin de branche), le cap ne
    // s'activait jamais et l'identite etait triviale. Trois variantes :
    //   o1  sigma3 = 20 MPa, cap 440 : activation JAMAIS (mesuree : pc reste a
    //       pc0), q identique a 1e-9 sur les 15 000 pas — invariance stricte
    //       d'une cle posee mais inactive ;
    //   o2  sigma3 = 20 MPa, cap 200 < 256 : identite a 1e-9 AVANT le premier
    //       pas ou pc > pc0 (post-pic, d_c ~ 0,8 : le seuil est sur la pression
    //       EFFECTIVE p_nom/(1-d), 57 MPa nominaux), puis les courbes DOIVENT
    //       diverger (|dq| > 0,1 MPa, wPlas plus grand, compaction > 0), q_pic
    //       intact (activation apres le pic) — la variante qui falsifie ;
    //   o3  sigma3 = 100 MPa, cap 100 : cap actif des que p_bar >= 100 MPa
    //       (pas ~540, ou q repasse par 0 : le pilote part de sigma = 0) puis
    //       cone actif en meme temps sur toute la branche plastique : q_pic
    //       invariant (0,1 %), eps_ax au pic decale de +eps_v^pl(cap)/3 (10 %),
    //       et la borne documentee du depassement (cap impose au predicteur,
    //       Lubliner dilatant ensuite) p_bar_fin - pc = K tan(psi) d lambda =
    //       K (d tr eps_pl + d pc/H), verifiee pas a pas a 1e-6 fc0.
    {
        struct TTrace {
            std::vector<double> eAx, q, pBar, pc, trEp, w, dC, epsC;
            MatState st;
        };
        auto tri = [&](const std::string& cfgText, const char* tag, double s3, TTrace& T) {
            auto law = mk(cfgText, lc);
            MatState st;
            Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
            double e = 0.0, x = 0.0;
            T = TTrace{};
            for (int k = 0; k < 15000; ++k) {
                e -= 1.0e-6;
                hold.solve(*law, st, eps, lc, 1.0, x, kLat,
                           [&](Eigen::Matrix3d& ee, double xx) { setZ(ee, e); setLat(ee, xx); },
                           [&](const Eigen::Matrix3d& sg) { return 0.5 * (sg(0, 0) + sg(1, 1)) + s3; },
                           sig);
                const double q = -sig(2, 2) - s3;
                T.eAx.push_back(e);
                T.q.push_back(q);
                T.pBar.push_back(-sig.trace() / 3.0 / (1.0 - st.cdp.d));   // effective
                T.pc.push_back(st.pc);
                T.trEp.push_back(st.epsP.trace());
                T.w.push_back(st.wPlas);
                T.dC.push_back(st.Dc);
                T.epsC.push_back(st.cdp.epsCpl);
                if (k % 25 == 0) row(tag, "triax", s3, eps, q, st);
            }
            T.st = st;
        };
        auto capStr = [&](double pc0) {
            std::ostringstream kk;
            kk.precision(17);
            kk << "cdpCap = true\ncdpCapP0 = " << pc0 << "\ncdpCapH = " << K << "\n";
            return kk.str();
        };
        auto argmax = [](const std::vector<double>& v) {
            return (int)(std::max_element(v.begin(), v.end()) - v.begin());
        };
        // premier pas ou le cap a agi (pc > pc0 : lazy init a pc0, puis += H dev)
        auto firstAct = [](const TTrace& B, double pc0) {
            for (std::size_t k = 0; k < B.pc.size(); ++k)
                if (B.pc[k] > pc0) return (int)k;
            return -1;
        };
        auto compare = [&](const TTrace& A, const TTrace& B, int kAct, double& dqBefore,
                           int& nBefore, double& dqAfter) {
            dqBefore = 0.0; nBefore = 0; dqAfter = 0.0;
            for (std::size_t k = 0; k < A.q.size(); ++k) {
                if (kAct < 0 || (int)k < kAct) {
                    dqBefore = std::max(dqBefore, std::abs(B.q[k] - A.q[k])
                                                      / std::max(std::abs(A.q[k]), 1.0));
                    ++nBefore;
                } else {
                    dqAfter = std::max(dqAfter, std::abs(B.q[k] - A.q[k]));
                }
            }
        };
        const double fc0 = 316.8e6;                         // fc0_pic carte historique
        // ---- o1 : 20 MPa, cap 440 — invariance stricte (cap jamais actif) -----
        {
            TTrace A, B;
            tri("", "o1_nocap", 20e6, A);
            tri(capStr(440e6), "o1_cap440", 20e6, B);
            const int kAct = firstAct(B, 440e6);
            double dqB, dqA, pMax = 0.0;
            int nB;
            compare(A, B, kAct, dqB, nB, dqA);
            for (double p : A.pBar) pMax = std::max(pMax, p);
            std::cout << "[cdp] (o1) triaxial 20 MPa, cap 440 : p_bar max sans cap " << pMax / 1e6
                      << " MPa (sigma3 effectif fin " << 20.0 / (1.0 - A.st.cdp.d) << " MPa, d_c "
                      << A.st.Dc << "), premier pas cappe (pc > pc0) : " << kAct
                      << (kAct < 0 ? " (jamais)" : "") << ", " << nB << " pas compares, ecart relatif max sur q "
                      << dqB << ", pc fin " << B.pc.back() / 1e6 << " MPa ; q pic sans/avec "
                      << A.q[argmax(A.q)] / 1e6 << " / " << B.q[argmax(B.q)] / 1e6 << " MPa\n";
            verdict("(o1) cap 440 : jamais actif (pc reste a pc0) et p_bar max < 440 MPa (MESURE)",
                    kAct < 0 && pMax < 440e6);
            verdict("(o1) cap 440 inactif : q(eps) identique a 1e-9 sur les 15 000 pas (invariance stricte)",
                    nB == 15000 && dqB < 1e-9);
        }
        // ---- o2 : 20 MPa, cap 200 < p_bar max — DOIT diverger apres activation --
        {
            TTrace A, B;
            tri("", "o2_nocap", 20e6, A);
            tri(capStr(200e6), "o2_cap200", 20e6, B);
            const int kAct = firstAct(B, 200e6);
            double dqB, dqA;
            int nB;
            compare(A, B, kAct, dqB, nB, dqA);
            const int kPk = argmax(A.q);
            const double compB = -B.st.epsP.trace(), compA = -A.st.epsP.trace();
            std::cout << "[cdp] (o2) triaxial 20 MPa, cap 200 : premier pas cappe " << kAct
                      << (kAct >= 0 ? " (eps_ax " + std::to_string(B.eAx[kAct]) + ", p_nom "
                                          + std::to_string(-A.pBar[kAct] * (1.0 - A.dC[kAct]) / 1e6)
                                          + " MPa, d_c " + std::to_string(A.dC[kAct]) + ", p_bar "
                                          + std::to_string(A.pBar[kAct] / 1e6) + " MPa)"
                                    : std::string(" (jamais)"))
                      << " ; avant : " << nB << " pas, ecart relatif max sur q " << dqB
                      << " ; apres : ecart max |dq| " << dqA / 1e6 << " MPa, wPlas sans/avec "
                      << A.w.back() << " / " << B.w.back() << " J/m^3, -tr(eps_pl) fin sans/avec "
                      << compA << " / " << compB << ", pc fin " << B.pc.back() / 1e6
                      << " MPa, d_c fin sans/avec " << A.st.Dc << " / " << B.st.Dc
                      << " ; q pic sans/avec " << A.q[kPk] / 1e6 << " / " << B.q[argmax(B.q)] / 1e6
                      << " MPa (pic au pas " << kPk << ")\n";
            verdict("(o2) cap 200 : DOIT s'activer (pc > pc0) sur la branche a 20 MPa", kAct >= 0);
            verdict("(o2) cap 200 : activation apres le pic (seuil sur p_bar effectif = p_nom/(1-d))",
                    kAct > kPk);
            verdict("(o2) cap 200 : q(eps) identique a 1e-9 AVANT le premier pas cappe (" + std::to_string(nB) + " pas)",
                    nB > 1000 && dqB < 1e-9);
            verdict("(o2) DOIT DIVERGER : |dq| > 0,1 MPa apres l'activation (mesure " + std::to_string(dqA / 1e6) + " MPa)",
                    dqA > 0.1e6);
            verdict("(o2) DOIT DIVERGER : wPlas plus grand et compaction supplementaire > 0 avec le cap",
                    B.w.back() > A.w.back() && compB - compA > 0.0);
            check("(o2) q_pic avec cap / sans cap (activation post-pic)", B.q[argmax(B.q)], A.q[kPk], 0.1);
        }
        // ---- o3 : 100 MPa, cap 100 — cap et cone actifs ensemble ----------------
        {
            const double s3 = 100e6, pc0 = 100e6, H = K;
            TTrace A, B;
            tri("", "o3_nocap", s3, A);
            tri(capStr(pc0), "o3_cap100", s3, B);
            const int kAct = firstAct(B, pc0);
            const int kA = argmax(A.q), kB = argmax(B.q);
            // depassement p_bar - pc en fin de pas sur les pas ou le cap a agi
            double errBound = 0.0, overMax = 0.0, overMaxRel = 0.0;
            int nCap = 0, nBoth = 0;
            for (std::size_t k = 0; k < B.q.size(); ++k) {
                const double pcPrev = k ? B.pc[k - 1] : pc0;
                const double trPrev = k ? B.trEp[k - 1] : 0.0;
                const double ecPrev = k ? B.epsC[k - 1] : 0.0;
                if (B.pc[k] > pcPrev) {
                    ++nCap;
                    if (B.epsC[k] > ecPrev) ++nBoth;
                    const double over = B.pBar[k] - B.pc[k];
                    const double pred = K * ((B.trEp[k] - trPrev) + (B.pc[k] - pcPrev) / H);
                    errBound = std::max(errBound, std::abs(over - pred));
                    overMax = std::max(overMax, over);
                    overMaxRel = std::max(overMaxRel, over / B.pc[k]);
                }
            }
            const double shiftMeas = B.eAx[kB] - A.eAx[kA];              // > 0 : decale vers la droite (|eps| plus grand)... signe : eps_ax < 0
            const double shiftPred = -(B.pc[kB] - pc0) / H / 3.0;         // compaction axiale = eps_v^pl(cap)/3, eps_ax negatif
            std::cout << "[cdp] (o3) triaxial 100 MPa, cap 100 : premier pas cappe " << kAct
                      << ", pas cappes " << nCap << " dont " << nBoth << " avec le cone actif (d eps_c_pl > 0) ; "
                      << "q pic sans/avec " << A.q[kA] / 1e6 << " / " << B.q[kB] / 1e6
                      << " MPa a eps_ax " << A.eAx[kA] << " / " << B.eAx[kB] << " (decalage mesure "
                      << shiftMeas << ", attendu -eps_v^pl(cap)/3 = " << shiftPred << ", pc au pic "
                      << B.pc[kB] / 1e6 << " MPa) ; depassement p_bar - pc max " << overMax / 1e6
                      << " MPa (" << 100.0 * overMaxRel << " % de pc), ecart max a K (d tr eps_pl + d pc/H) "
                      << errBound << " Pa ; fin : pc " << B.pc.back() / 1e6 << " MPa, -tr(eps_pl) sans/avec "
                      << -A.st.epsP.trace() << " / " << -B.st.epsP.trace() << ", wPlas sans/avec "
                      << A.w.back() << " / " << B.w.back() << " J/m^3, d_c fin " << A.st.Dc << " / " << B.st.Dc << "\n";
            // le pilote du banc (a) part de sigma = 0 : sigma_lat est tenue a
            // -sigma3 des le pas 1 mais sigma_ax = lambda tr(eps) ~ -58 MPa, p_bar
            // ~ 86 MPa ; l'axial rattrape sigma3 a E de = 77,7 kPa/pas, p_bar
            // atteint pc0 = sigma3 au pas ou q repasse par 0 (~540), bien avant
            // le pic (~7900). En matpoint (confinement d'abord) ce serait le pas 0.
            int kQ0 = -1;
            for (std::size_t k = 0; k < A.q.size() && kQ0 < 0; ++k)
                if (A.q[k] >= 0.0) kQ0 = (int)k;
            std::cout << "[cdp] (o3) premier pas a q >= 0 (sigma_ax = -sigma3, p_bar = pc0) : " << kQ0
                      << ", premier pas cappe : " << kAct << ", pic au pas " << kA << "\n";
            verdict("(o3) cap 100 a sigma3 = 100 : actif des que p_bar >= pc0 (pas ou q >= 0, +-1) et bien avant le pic",
                    kAct >= 0 && kQ0 >= 0 && std::abs(kAct - kQ0) <= 1 && kAct < kA / 5);
            verdict("(o3) cap et cone actifs dans un meme pas (couplage exerce, > 1000 pas)", nBoth > 1000);
            check("(o3) q_pic invariant sous le cap (avec / sans)", B.q[kB], A.q[kA], 0.1);
            check("(o3) decalage de eps_ax au pic = -eps_v^pl(cap)/3", shiftMeas, shiftPred, 10.0);
            verdict("(o3) borne du depassement : p_bar_fin - pc = K (d tr eps_pl + d pc/H) a 1e-6 fc0 ("
                        + std::to_string(errBound) + " Pa)", nCap > 0 && errBound < 1e-6 * fc0);
            verdict("(o3) depassement > 0 (Lubliner dilatant apres le cap) et < 5 % de pc",
                    overMax > 0.0 && overMaxRel < 0.05);
            verdict("(o3) compaction supplementaire > 0 et wPlas plus grand avec le cap",
                    -B.st.epsP.trace() > -A.st.epsP.trace() && B.w.back() > A.w.back());
        }
    }

    std::cout << "[cdp] pilotage lateral : " << hold.calls << " appels a la loi, pire residu "
              << hold.worst << " Pa, pas non converges a 1e-3 Pa : " << hold.fails << "\n";
    verdict("pilotage lateral : 0 pas non converge", hold.fails == 0);
    std::cout << "[cdp] " << (ok ? "OK" : "FAIL") << "\n";
    return ok ? 0 : 1;
}

// ---------------------------------------------------------------------------
// rockim matpoint <cfg> [out.csv] — pilote point-materiel GENERIQUE (toutes
// les lois du noyau : elastic, dpr, mc, saksala, saksala2011, dpdfh, cdp),
// le moteur de la calibration. Cles (toutes optionnelles) :
//   mpPath = triax (defaut) | tension | biaxial | uniaxial
//   mpSigma3 = liste de pressions [Pa] separees par des espaces (defaut 0)
//   mpStrainMax (0,03), mpSteps (3000), mpLc (1e-3 m), mpDt (1 s, lois
//   visqueuses), mpConfineFirst (true : consolidation isotrope a sigma3 puis
//   phase axiale ; eps comptes depuis la fin de la consolidation)
// Chemins (axe pilote z, indice 2) : triax/uniaxial = eps_zz decroissante,
// sigma_xx = sigma_yy = -sigma3 tenues ; tension = eps_zz croissante, memes
// laterales ; biaxial = eps_xx = eps_yy decroissantes, sigma_zz = -sigma3
// tenue. Sortie CSV : sigma3, eps_ax, eps_lat, eps_vol, q, sig_ax, sig_lat,
// d_t, d_c, eps_c_pl, eps_t_pl, wPlas, wDamT, wDamC, residu, eps_in_c avec
// q = -sig_ax - sigma3 (compression positive ; en tension q = sig_ax -
// sig_lat, traction positive), d_t = MatState::D, d_c = Dc, eps_c_pl =
// epvEq, eps_t_pl = kappa (pour cdp : exactement d_t, d_c, eps_c_pl,
// eps_t_pl ; pour dpr : D, omega_c, eps_vp equivalente, kappa de Rankine),
// residu = |residu final du pilotage lateral| [Pa] de la ligne (revue
// 2026-09-04 ; > 1e-3 = pas non converge), eps_in_c = eps_c_pl + d_c/(1-d_c)
// sigma_c/E0 (deformation inelastique d'Abaqus, table lue a mpLc : avec
// cdpCompLength, u_in = (eps_in_c - eps_in,pic) x mpLc est invariant en mpLc
// — 2026-09-04 soir). Code retour 0, ou 2 si au moins un pas n'a pas converge.
// Avec cdpViscosity > 0, mpDt = Delta eps / epsdot (pas de deformation /
// vitesse de deformation) : le defaut 1 s n'a de sens que pour les lois
// rate-independantes. En traction post-pic les pas doivent rester <= ~1e-5.
// ---------------------------------------------------------------------------
int matpointDrive(const Config& cfg, const std::string& csvPath) {
    Material m = Material::from(cfg);
    PhaseSet::validate(m, "global");
    const double lc = cfg.getd("mpLc", 1.0e-3);
    const double dt = cfg.getd("mpDt", 1.0);
    std::string kind = cfg.gets("law", "dpr");
    auto law = MatLaw::make(kind, m, cfg, lc);
    std::string path = cfg.gets("mpPath", "triax");
    if (path != "triax" && path != "tension" && path != "biaxial" && path != "uniaxial"
        && path != "hydro" && path != "oedo")
        throw std::runtime_error("mpPath must be triax | tension | biaxial | uniaxial | hydro | oedo");
    std::vector<double> s3s;
    {
        std::istringstream ss(cfg.gets("mpSigma3", "0"));
        std::string tok;
        while (ss >> tok) {
            std::size_t used = 0;
            double v = 0.0;
            try { v = std::stod(tok, &used); } catch (const std::exception&) { used = 0; }
            if (used != tok.size())
                throw std::runtime_error("mpSigma3: token '" + tok + "' is not a number");
            if (v < 0.0) throw std::runtime_error("mpSigma3: pressures are >= 0 [Pa]");
            s3s.push_back(v);
        }
        if (s3s.empty()) s3s.push_back(0.0);
    }
    if (path == "uniaxial") s3s = {0.0};
    const double emax = cfg.getd("mpStrainMax", 0.03);
    const int nSteps = cfg.geti("mpSteps", 3000);
    const bool confFirst = cfg.getb("mpConfineFirst", true);
    if (!(emax > 0.0) || nSteps < 1 || !(lc > 0.0) || !(dt > 0.0))
        throw std::runtime_error("matpoint: mpStrainMax > 0, mpSteps >= 1, mpLc > 0, mpDt > 0");
    const double lam = m.E * m.nu / ((1.0 + m.nu) * (1.0 - 2.0 * m.nu));
    const double G = m.G(), K = m.K();
    const double kLat = 2.0 * (lam + G), kAx = lam + 2.0 * G;

    // ---- mpPath = hydro (2026-09-04, cap volumique) : deformation ISOTROPE
    // imposee eps = -e I, e de 0 a mpStrainMax en mpSteps pas (aucun pilotage
    // lateral, mpSigma3 et mpConfineFirst ignores). Sortie propre : eps_iso,
    // eps_v (= -3e), p_nom (= -tr sigma/3, compression positive), pc (etat
    // MatState::pc : cap de cdp/dpr/saksala, 0 sans cap), eps_v_pl (= -tr
    // eps_pl, compaction positive), d (cdp : d total ; autres : D), wPlas,
    // eroded. Les autres chemins ne sont pas touches.
    if (path == "hydro") {
        std::ofstream outH(csvPath);
        outH << "eps_iso,eps_v,p_nom,pc,eps_v_pl,d,wPlas,eroded\n";
        outH.precision(10);
        std::cout << "[matpoint] law = " << law->name() << ", path hydro (deformation "
                     "isotrope imposee, mpSigma3 / mpConfineFirst ignores), lc " << lc
                  << " m, " << nSteps << " pas a " << emax << "\n";
        auto t0 = std::chrono::steady_clock::now();
        MatState st;
        Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
        const double de = emax / nSteps;
        double pMax = 0.0, evAtCap = 0.0;
        bool capped = false;
        for (int k = 0; k < nSteps; ++k) {
            eps = -(k + 1) * de * Eigen::Matrix3d::Identity();
            sig = law->stress(eps, st, dt, lc);
            const double p = -sig.trace() / 3.0, ev = eps.trace();
            const double evPl = -st.epsP.trace() + 0.0;      // + 0.0 : pas de « -0 »
            const double d = law->name() == "cdp" ? st.cdp.d : st.D;
            pMax = std::max(pMax, p);
            if (!capped && evPl > 0.0) { capped = true; evAtCap = ev; }
            outH << (k + 1) * de << "," << ev << "," << p << "," << st.pc << "," << evPl
                 << "," << d << "," << st.wPlas << "," << (st.eroded ? 1 : 0) << "\n";
            if (st.eroded) break;
        }
        auto t1 = std::chrono::steady_clock::now();
        std::cout << "[matpoint] hydro : p_nom max " << pMax / 1e6 << " MPa, pc fin "
                  << st.pc / 1e6 << " MPa, eps_v^pl fin " << -st.epsP.trace()
                  << (capped ? " (compaction depuis eps_v = " + std::to_string(evAtCap) + ")"
                             : " (aucune compaction)")
                  << ", wPlas " << st.wPlas << " J/m^3 ("
                  << std::chrono::duration<double, std::milli>(t1 - t0).count() << " ms)\n";
        return 0;
    }

    // ---- mpPath = oedo (2026-09-05, element unique K0) : deformation UNIAXIALE
    // imposee eps = diag(0, 0, -e), faces laterales bloquees (aucun pilotage,
    // mpSigma3 / mpConfineFirst ignores) — le trajet du pole de l'insert.
    // Memes colonnes que les chemins pilotes (eps_lat = 0, sig_lat = sigma_xx).
    if (path == "oedo") {
        std::ofstream outO(csvPath);
        outO << "sigma3,eps_ax,eps_lat,eps_vol,q,sig_ax,sig_lat,d_t,d_c,eps_c_pl,eps_t_pl,"
                "wPlas,wDamT,wDamC,residu,eps_in_c\n";
        outO.precision(10);
        std::cout << "[matpoint] law = " << law->name() << ", path oedo (deformation "
                     "uniaxiale imposee, faces laterales bloquees), lc " << lc << " m, "
                  << nSteps << " pas a " << emax << "\n";
        MatState st;
        Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), sig;
        const double de = emax / nSteps;
        double sMax = 0.0;
        for (int k = 0; k < nSteps; ++k) {
            eps(2, 2) = -(k + 1) * de;
            sig = law->stress(eps, st, dt, lc);
            const double sAx = sig(2, 2), sLat = 0.5 * (sig(0, 0) + sig(1, 1));
            sMax = std::max(sMax, -sAx);
            outO << 0.0 << "," << eps(2, 2) << "," << 0.0 << "," << eps.trace() << ","
                 << (sLat - sAx) << "," << sAx << "," << sLat << "," << st.D << "," << st.Dc << ","
                 << (law->name() == "cdp" ? st.cdp.epsCpl : st.epvEq) << ","
                 << (law->name() == "cdp" ? st.cdp.epsTpl : st.kappa) << ","
                 << st.wPlas << "," << st.wDamT << "," << st.wDamC << ",0,0\n";
            if (st.eroded) break;
        }
        std::cout << "[matpoint] oedo : sig_ax max " << sMax / 1e6 << " MPa, sig_ax fin "
                  << -sig(2, 2) / 1e6 << " MPa, sig_lat fin " << -0.5 * (sig(0, 0) + sig(1, 1)) / 1e6
                  << " MPa, d_c fin " << st.Dc << ", eps_c_pl fin "
                  << (law->name() == "cdp" ? st.cdp.epsCpl : st.epvEq) << "\n";
        return 0;
    }

    std::ofstream out(csvPath);
    // eps_in_c (2026-09-04 soir, en FIN de ligne) : deformation inelastique de
    // compression au sens d'Abaqus, eps_in = eps_c_pl + d_c/(1-d_c) sigma_c/E0,
    // pour reconstruire u_in = (eps_in - eps_in,pic) x mpLc avec cdpCompLength
    // (cdp : sigma_c, d_c lus dans la table a mpLc ; autres lois : epvEq +
    // Dc/(1-Dc) |sig_ax|/E0, proxy uniaxial)
    out << "sigma3,eps_ax,eps_lat,eps_vol,q,sig_ax,sig_lat,d_t,d_c,eps_c_pl,eps_t_pl,"
           "wPlas,wDamT,wDamC,residu,eps_in_c\n";
    out.precision(10);
    const CdpLaw* cdpLaw = dynamic_cast<const CdpLaw*>(law.get());
    std::cout << "[matpoint] law = " << law->name() << ", path " << path << ", lc " << lc
              << " m, " << nSteps << " pas a " << emax << ", confinement d'abord : "
              << (confFirst ? "oui" : "non")
              << (cdpLaw && cdpLaw->compLength() > 0.0
                      ? " ; crack band en compression L_ref = " + std::to_string(cdpLaw->compLength())
                            + " m lu a mpLc"
                      : std::string())
              << "\n";
    Hold1D hold;
    for (double s3 : s3s) {
        auto t0 = std::chrono::steady_clock::now();
        MatState st;
        Eigen::Matrix3d eps = Eigen::Matrix3d::Zero(), eps0 = Eigen::Matrix3d::Zero(), sig;
        sig.setZero();
        if (confFirst && s3 > 0.0) {
            // consolidation isotrope : eps = e I, p = -sigma3 tenue par le
            // pilote (rigidite 3K), 100 pas de rampe puis un ajustement final
            double eIso = 0.0;
            for (int k = 1; k <= 100; ++k) {
                double target = -s3 * k / 100.0;
                hold.solve(*law, st, eps, lc, dt, eIso, 3.0 * K,
                           [&](Eigen::Matrix3d& ee, double xx) { ee = xx * Eigen::Matrix3d::Identity(); },
                           [&](const Eigen::Matrix3d& sg) { return sg.trace() / 3.0 - target; }, sig);
            }
            eps0 = eps;
        }
        const bool ten = path == "tension", bia = path == "biaxial";
        const double de = (ten ? 1.0 : -1.0) * emax / nSteps;
        double e = bia ? eps0(0, 0) : eps0(2, 2);
        double x = bia ? eps0(2, 2) : eps0(0, 0);
        double qmax = -1e300, eAtMax = 0.0;
        for (int k = 0; k < nSteps; ++k) {
            e += de;
            if (bia)
                hold.solve(*law, st, eps, lc, dt, x, kAx,
                           [&](Eigen::Matrix3d& ee, double xx) { setLat(ee, e); setZ(ee, xx); },
                           [&](const Eigen::Matrix3d& sg) { return sg(2, 2) + s3; }, sig);
            else
                hold.solve(*law, st, eps, lc, dt, x, kLat,
                           [&](Eigen::Matrix3d& ee, double xx) { setZ(ee, e); setLat(ee, xx); },
                           [&](const Eigen::Matrix3d& sg) { return 0.5 * (sg(0, 0) + sg(1, 1)) + s3; },
                           sig);
            double sAx = bia ? sig(0, 0) : sig(2, 2);
            double sLat = bia ? sig(2, 2) : 0.5 * (sig(0, 0) + sig(1, 1));
            double eAx = (bia ? eps(0, 0) - eps0(0, 0) : eps(2, 2) - eps0(2, 2));
            double eLat = (bia ? eps(2, 2) - eps0(2, 2) : eps(0, 0) - eps0(0, 0));
            double eVol = (eps - eps0).trace();
            double q = ten ? sAx - sLat : -sAx - s3;
            if (q > qmax) { qmax = q; eAtMax = eAx; }
            double einC;
            if (cdpLaw) {
                double sc, dc;
                cdpLaw->compState(st.cdp.epsCpl, lc, sc, dc);
                einC = st.cdp.epsCpl + dc / (1.0 - dc) * sc / m.E;
            } else {
                einC = st.epvEq + (st.Dc < 1.0 ? st.Dc / (1.0 - st.Dc) * std::max(-sAx, 0.0) / m.E : 0.0);
            }
            out << s3 << "," << eAx << "," << eLat << "," << eVol << "," << q << ","
                << sAx << "," << sLat << "," << st.D << "," << st.Dc << ","
                << (law->name() == "cdp" ? st.cdp.epsCpl : st.epvEq) << ","
                << (law->name() == "cdp" ? st.cdp.epsTpl : st.kappa) << ","
                << st.wPlas << "," << st.wDamT << "," << st.wDamC << ","
                << hold.last << "," << einC << "\n";
            if (st.eroded) break;
        }
        auto t1 = std::chrono::steady_clock::now();
        std::cout << "[matpoint] sigma3 " << s3 / 1e6 << " MPa : q max " << qmax / 1e6
                  << " MPa a eps_ax " << eAtMax << " ("
                  << std::chrono::duration<double, std::milli>(t1 - t0).count() << " ms)\n";
    }
    std::cout << "[matpoint] pilotage : " << hold.calls << " appels, pire residu "
              << hold.worst << " Pa, pas non converges : " << hold.fails << "\n";
    // revue 2026-09-04 : un pilotage rompu est FALSIFIANT au niveau du code de
    // retour (2), pour qu'une boucle de calibration Python le voie ; la
    // colonne `residu` du CSV marque chaque ligne
    if (hold.fails > 0) {
        std::cout << "[matpoint] WARN : " << hold.fails << " pas non converges (residu > "
                  << hold.tol << " Pa) — lignes marquees dans la colonne residu ; "
                     "code retour 2\n";
        return 2;
    }
    return 0;
}

// ===========================================================================
// tensionDamage = fixed — bancs falsifiants au point materiel (2026-09-05).
// Carte Red Bohus dpr (E 77,66 GPa, nu 0,29, ft 9 MPa, Gf 100 J/m^2, lc 1 mm ;
// dpApex, dpTension = off, rankineDrive = stress : la traction est l'affaire
// du cut-off seul, le cone ne coule pas). Pilotage en DEFORMATION : l'etat
// « contrainte uniaxiale E e selon n » est impose par eps = e n n^T
// - nu e (I - n n^T) (exact tant que l'espace effectif reste elastique).
//  (a) traction x jusqu'a d = 0,9, decharge a zero, traction y (< ft) : la
//      pente sigma_yy / eps_yy vaut E pour fixed (a 1 %) et (1 - 0,9) E pour
//      scalar — la reponse scalaire est REJETEE par le meme critere ; puis
//      traction y au-dela de ft : 2e direction figee = y, d1 inchange.
//  (b) traction x jusqu'a d = 0,9, decharge, compression x : E dans les deux
//      modeles (unilateral).
//  (c) traction a 45 deg dans xy : eul = (45, 0, 0) deg, n1 = (1, 1, 0)/sqrt2,
//      et sigma_45 = Rz sigma_x Rz^T (objectivite) a 1e-8.
//  (d) rotation des directions principales apres amorcage (d = 0,5 en x) :
//      (d1) tel qu'ecrit — eps_xx maintenu, eps_yy croissant, sigma_zz = 0 :
//      sigma_xy = 0 par symetrie (axes = repere), y s'amorce a S_22 = ft ;
//      (d2) traction croissante selon 45 deg avec x maintenu : cisaillement
//      parasite S_12 nominal / effectif = 1 - beta max(d_i^ouv), desalignement
//      des directions principales nominales et effectives — pour beta = 1 et
//      0,1, DONNEES sans critere.
//  (e) traction x jusqu'a sigma < 1e-4 ft : int sigma deps x lc = Gf et
//      wDamT x lc = Gf a 2 %, les deux modeles.
//  (0) bit-identite au point : `tensionDamage = scalar` explicite = sans cle
//      (memcmp des contraintes) ; cles invalides refusees.
// ===========================================================================
int fixedCrackSelftest(const std::string& csvPath) {
    Material m;
    m.E = 77.66e9; m.nu = 0.29; m.rho = 2620.0;
    m.cohesion = 22.77e6; m.phiDeg = 50.4; m.ft = 9.0e6; m.Gf = 100.0;
    const double lc = 1.0e-3, dt = 1.0e-6, nu = m.nu;
    const double k0 = m.ft / m.E, kf = m.Gf / (lc * m.ft) - 0.5 * k0;
    const std::string base = "law = dpr\ndpApex = true\ndpTension = off\n"
                             "rankineDrive = stress\nerodeD = 2\nerodeEpv = 0\n";
    auto writeCfg = [&](const std::string& text) {
        std::string p = csvPath + ".tmp.cfg";
        std::ofstream f(p);
        f << text;
        f.close();
        return Config::load(p);
    };
    auto mk = [&](const std::string& extra) {
        Config cfg = writeCfg(base + extra);
        return MatLaw::make("dpr", m, cfg, lc);
    };
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    auto uni = [&](double e, const Eigen::Vector3d& n) -> Eigen::Matrix3d {
        return e * n * n.transpose() - nu * e * (I - n * n.transpose());
    };
    std::ofstream out(csvPath);
    out.precision(15);
    out << "test,model,step,exx,eyy,ezz,exy,sxx,syy,szz,sxy,D,d1,d2,d3,nAct,"
           "eul1deg,eul2deg,eul3deg,wDamT\n";
    auto rec = [&](const char* test, const char* model, int step,
                   const Eigen::Matrix3d& e, const Eigen::Matrix3d& sg,
                   const MatState& st) {
        const double r2d = 180.0 / M_PI;
        out << test << "," << model << "," << step << "," << e(0, 0) << "," << e(1, 1)
            << "," << e(2, 2) << "," << e(0, 1) << "," << sg(0, 0) << "," << sg(1, 1)
            << "," << sg(2, 2) << "," << sg(0, 1) << "," << st.D << "," << st.fcm.d[0]
            << "," << st.fcm.d[1] << "," << st.fcm.d[2] << "," << st.fcm.nAct
            << "," << st.fcm.eul[0] * r2d << "," << st.fcm.eul[1] * r2d
            << "," << st.fcm.eul[2] * r2d << "," << st.wDamT << "\n";
    };
    int nPass = 0, nFail = 0;
    auto check = [&](const std::string& what, double got, double want, double tol,
                     bool rel) {
        double err = std::abs(got - want);
        bool ok = rel ? err <= tol * std::abs(want) : err <= tol;
        std::cout << "[fixed] " << what << " : " << got << " (attendu " << want
                  << ", " << (rel ? "rel " : "abs ") << tol << ") -> "
                  << (ok ? "PASS" : "FAIL") << "\n";
        (ok ? nPass : nFail)++;
        return ok;
    };
    auto note = [&](const std::string& what) { std::cout << "[fixed] " << what << "\n"; };
    const Eigen::Vector3d ex(1, 0, 0), ey(0, 1, 0), ez(0, 0, 1);
    const Eigen::Vector3d n45 = (ex + ey).normalized();
    const double de = k0 / 20.0;
    struct Trace { std::vector<double> e; std::vector<Eigen::Matrix3d> sig; };

    // traction uniaxiale selon n jusqu'a s.D >= dTarget ; rend e atteint
    auto pull = [&](MatLaw& law, MatState& st, const Eigen::Vector3d& n, double dTarget,
                    const char* test, const char* model, Trace* tr) {
        double e = 0.0;
        int step = 0;
        while (st.D < dTarget && e < 0.5) {
            e += de;
            Eigen::Matrix3d eps = uni(e, n);
            Eigen::Matrix3d sg = law.stress(eps, st, dt, lc);
            rec(test, model, step++, eps, sg, st);
            if (tr) { tr->e.push_back(e); tr->sig.push_back(sg); }
        }
        return e;
    };
    // rampe lineaire de epsA a epsB en n pas ; cb(k, eps, sigma, etat)
    auto ramp = [&](MatLaw& law, MatState& st, const Eigen::Matrix3d& epsA,
                    const Eigen::Matrix3d& epsB, int n, const char* test, const char* model,
                    const std::function<void(int, const Eigen::Matrix3d&,
                                             const Eigen::Matrix3d&, const MatState&)>& cb) {
        for (int k = 1; k <= n; ++k) {
            double t = (double)k / n;
            Eigen::Matrix3d eps = (1.0 - t) * epsA + t * epsB;
            Eigen::Matrix3d sg = law.stress(eps, st, dt, lc);
            rec(test, model, k, eps, sg, st);
            if (cb) cb(k, eps, sg, st);
        }
    };
    // pente par moindres carres sigma = a eps + b sur des paires (eps, sigma)
    auto slope = [](const std::vector<double>& x, const std::vector<double>& y) {
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        std::size_t n = x.size();
        for (std::size_t i = 0; i < n; ++i) { sx += x[i]; sy += y[i]; sxx += x[i] * x[i]; sxy += x[i] * y[i]; }
        return (n * sxy - sx * sy) / (n * sxx - sx * sx);
    };

    std::cout << "[fixed] carte : E " << m.E << " Pa, nu " << nu << ", ft " << m.ft
              << " Pa, Gf " << m.Gf << " J/m^2, lc " << lc << " m ; k0 = ft/E = " << k0
              << ", kf = Gf/(lc ft) - k0/2 = " << kf << " ; Gf/lc = " << m.Gf / lc << " J/m^3\n";

    // ---- (0) cles : validation, bit-identite scalar / cle explicite ----------
    {
        auto expectThrow = [&](const std::string& extra, const std::string& kind,
                               const std::string& tag) {
            bool thrown = false;
            std::string msg;
            try {
                Config cfg = writeCfg(base + extra);
                auto l = MatLaw::make(kind, m, cfg, lc);
                (void)l;
            } catch (const std::exception& ex) { thrown = true; msg = ex.what(); }
            std::cout << "[fixed] (0) cle refusee " << tag << " : "
                      << (thrown ? "PASS (" + msg.substr(0, 60) + ")" : "FAIL (acceptee)") << "\n";
            (thrown ? nPass : nFail)++;
        };
        expectThrow("tensionDamage = bogus\n", "dpr", "tensionDamage = bogus");
        expectThrow("tensionShearRetention = 0.5\n", "dpr", "tensionShearRetention sans fixed");
        expectThrow("tensionDamage = fixed\ntensionShearRetention = 1.5\n", "dpr", "beta = 1,5");
        expectThrow("tensionDamage = fixed\n", "cdp", "fixed + law = cdp");
        expectThrow("tensionDamage = fixed\n", "saksala2011", "fixed + law = saksala2011");
        auto A = mk(""), B = mk("tensionDamage = scalar\n");
        MatState sa, sb;
        double maxDiff = 0.0;
        int nNonId = 0, step = 0;
        auto cmp3 = [&](const Eigen::Matrix3d& x, const Eigen::Matrix3d& y) {
            for (int i = 0; i < 9; ++i) {
                double a = x.data()[i], b = y.data()[i];
                if (std::memcmp(&a, &b, sizeof(double)) != 0) ++nNonId;
                maxDiff = std::max(maxDiff, std::abs(a - b));
            }
        };
        Eigen::Matrix3d epsPrev = Eigen::Matrix3d::Zero();
        double e = 0.0;
        while (sa.D < 0.9 && e < 0.5) {
            e += de;
            Eigen::Matrix3d eps = uni(e, ex);
            cmp3(A->stress(eps, sa, dt, lc), B->stress(eps, sb, dt, lc));
            epsPrev = eps; ++step;
        }
        for (int k = 1; k <= 50; ++k) {
            Eigen::Matrix3d eps = (1.0 - k / 50.0) * epsPrev;
            cmp3(A->stress(eps, sa, dt, lc), B->stress(eps, sb, dt, lc)); ++step;
        }
        for (int k = 1; k <= 40; ++k) {
            Eigen::Matrix3d eps = uni(0.8 * k0 * k / 40.0, ey);
            cmp3(A->stress(eps, sa, dt, lc), B->stress(eps, sb, dt, lc)); ++step;
        }
        for (int k = 1; k <= 30; ++k) {
            Eigen::Matrix3d eps = uni(-3.0 * k0 * k / 30.0, ex);
            cmp3(A->stress(eps, sa, dt, lc), B->stress(eps, sb, dt, lc)); ++step;
        }
        std::cout << "[fixed] (0) bit-identite sans cle / tensionDamage = scalar : " << step
                  << " pas, composantes non identiques " << nNonId << ", ecart max " << maxDiff
                  << " Pa -> " << (nNonId == 0 ? "PASS" : "FAIL") << "\n";
        (nNonId == 0 ? nPass : nFail)++;
    }

    // ---- (a) traction x (d = 0,9), decharge, traction y ------------------
    Trace trX;                                   // trace x du modele fige (pour (c))
    {
        auto S = mk(""), F = mk("tensionDamage = fixed\n");
        MatState ss, sf;
        Trace tS;
        double eS = pull(*S, ss, ex, 0.9, "a_pull_x", "scalar", &tS);
        double eF = pull(*F, sf, ex, 0.9, "a_pull_x", "fixed", &trX);
        double dsig = 0.0;
        for (std::size_t i = 0; i < std::min(tS.sig.size(), trX.sig.size()); ++i)
            dsig = std::max(dsig, (tS.sig[i] - trX.sig[i]).norm());
        note("(a) traction x : d = 0,9 atteint a eps_xx = " + std::to_string(eS) + " (scalar, D "
             + std::to_string(ss.D) + ") / " + std::to_string(eF) + " (fixed, d1 "
             + std::to_string(sf.fcm.d[0]) + ", nAct " + std::to_string(sf.fcm.nAct)
             + ", eul deg " + std::to_string(sf.fcm.eul[0] * 180 / M_PI) + " "
             + std::to_string(sf.fcm.eul[1] * 180 / M_PI) + " "
             + std::to_string(sf.fcm.eul[2] * 180 / M_PI) + ")");
        check("(a) meme cinetique en traction uniaxiale x : max |sigma_scalar - sigma_fixed| / ft",
              dsig / m.ft, 0.0, 1.0e-9, false);
        check("(a) meme cinetique : d1(fixed) - D(scalar) a d = 0,9", sf.fcm.d[0] - ss.D, 0.0, 1.0e-12, false);
        const double d1Frozen = sf.fcm.d[0];
        // decharge a zero (50 pas)
        ramp(*S, ss, uni(eS, ex), Eigen::Matrix3d::Zero(), 50, "a_unload", "scalar", nullptr);
        ramp(*F, sf, uni(eF, ex), Eigen::Matrix3d::Zero(), 50, "a_unload", "fixed", nullptr);
        // traction y sous ft : eps = uni(e, y), e -> 0,8 k0
        std::vector<double> xS, yS, xF, yF;
        double sxxMaxF = 0.0;
        ramp(*S, ss, Eigen::Matrix3d::Zero(), uni(0.8 * k0, ey), 40, "a_pull_y", "scalar",
             [&](int, const Eigen::Matrix3d& e, const Eigen::Matrix3d& sg, const MatState&) {
                 xS.push_back(e(1, 1)); yS.push_back(sg(1, 1)); });
        ramp(*F, sf, Eigen::Matrix3d::Zero(), uni(0.8 * k0, ey), 40, "a_pull_y", "fixed",
             [&](int, const Eigen::Matrix3d& e, const Eigen::Matrix3d& sg, const MatState&) {
                 xF.push_back(e(1, 1)); yF.push_back(sg(1, 1));
                 sxxMaxF = std::max(sxxMaxF, std::abs(sg(0, 0))); });
        double EyS = slope(xS, yS), EyF = slope(xF, yF);
        check("(a) FIXED : pente sigma_yy/eps_yy apres fissure x (d1 = 0,9) / E", EyF / m.E, 1.0, 0.01, true);
        bool rejected = std::abs(EyS / m.E - 1.0) > 0.01;
        std::cout << "[fixed] (a) SCALAR : pente sigma_yy/eps_yy / E = " << EyS / m.E
                  << " (attendu (1 - 0,9) = 0,1 ; le critere E a 1 % la REJETTE) -> "
                  << (rejected ? "PASS (rejetee)" : "FAIL (acceptee)") << "\n";
        (rejected ? nPass : nFail)++;
        check("(a) SCALAR : pente / E = 1 - D", EyS / m.E, 1.0 - ss.D, 1.0e-6, true);
        check("(a) FIXED : |sigma_xx| nominal pendant la traction y / ft", sxxMaxF / m.ft, 0.0, 1.0e-9, false);
        check("(a) FIXED : d1 inchange par la traction y", sf.fcm.d[0], d1Frozen, 0.0, false);
        // traction y au-dela de ft : 2e direction = y
        double eInit = -1.0;
        int nActPrev = sf.fcm.nAct;
        ramp(*F, sf, uni(0.8 * k0, ey), uni(4.0 * k0, ey), 160, "a_pull_y_ft", "fixed",
             [&](int, const Eigen::Matrix3d& e, const Eigen::Matrix3d&, const MatState& st) {
                 if (st.fcm.nAct > nActPrev && eInit < 0.0) eInit = e(1, 1);
                 nActPrev = st.fcm.nAct; });
        Eigen::Matrix3d R = fcmFrame(sf.fcm.eul);
        note("(a) traction y au-dela de ft : nAct " + std::to_string(sf.fcm.nAct) + ", amorcage a eps_yy = "
             + std::to_string(eInit) + " (k0 = " + std::to_string(k0) + "), d2 fin " + std::to_string(sf.fcm.d[1])
             + ", d1 fin " + std::to_string(sf.fcm.d[0]));
        check("(a) FIXED : 2e direction amorcee a eps_yy (E eps_yy = ft) / k0", eInit / k0, 1.0, de / k0 + 1e-12, true);
        check("(a) FIXED : |n2 . y|", std::abs(R.col(1).dot(ey)), 1.0, 1.0e-9, false);
        check("(a) FIXED : |n1 . x| conserve", std::abs(R.col(0).dot(ex)), 1.0, 1.0e-9, false);
        check("(a) FIXED : d1 inchange par la fissuration en y", sf.fcm.d[0], d1Frozen, 0.0, false);
    }

    // ---- (b) traction x (d = 0,9), decharge, compression x ---------------
    {
        auto S = mk(""), F = mk("tensionDamage = fixed\n");
        MatState ss, sf;
        double eS = pull(*S, ss, ex, 0.9, "b_pull_x", "scalar", nullptr);
        double eF = pull(*F, sf, ex, 0.9, "b_pull_x", "fixed", nullptr);
        ramp(*S, ss, uni(eS, ex), Eigen::Matrix3d::Zero(), 50, "b_unload", "scalar", nullptr);
        ramp(*F, sf, uni(eF, ex), Eigen::Matrix3d::Zero(), 50, "b_unload", "fixed", nullptr);
        std::vector<double> xS, yS, xF, yF;
        ramp(*S, ss, Eigen::Matrix3d::Zero(), uni(-3.0 * k0, ex), 30, "b_comp_x", "scalar",
             [&](int, const Eigen::Matrix3d& e, const Eigen::Matrix3d& sg, const MatState&) {
                 xS.push_back(e(0, 0)); yS.push_back(sg(0, 0)); });
        ramp(*F, sf, Eigen::Matrix3d::Zero(), uni(-3.0 * k0, ex), 30, "b_comp_x", "fixed",
             [&](int, const Eigen::Matrix3d& e, const Eigen::Matrix3d& sg, const MatState&) {
                 xF.push_back(e(0, 0)); yF.push_back(sg(0, 0)); });
        check("(b) SCALAR : pente compression x apres fissure x / E", slope(xS, yS) / m.E, 1.0, 0.01, true);
        check("(b) FIXED : pente compression x apres fissure x / E", slope(xF, yF) / m.E, 1.0, 0.01, true);
        note("(b) sigma_xx fin : scalar " + std::to_string(yS.back()) + " Pa, fixed " + std::to_string(yF.back())
             + " Pa (E eps = " + std::to_string(-3.0 * k0 * m.E) + ") ; D/d1 conserves "
             + std::to_string(ss.D) + " / " + std::to_string(sf.fcm.d[0]));
    }

    // ---- (c) traction a 45 deg dans xy ----------------------------------
    {
        auto F = mk("tensionDamage = fixed\n");
        MatState sf;
        Trace tr45;
        pull(*F, sf, n45, 0.9, "c_pull_45", "fixed", &tr45);
        Eigen::Matrix3d Rz;
        const double c = std::cos(M_PI / 4), s = std::sin(M_PI / 4);
        Rz << c, -s, 0, s, c, 0, 0, 0, 1;
        double dsig = 0.0, smax = 0.0;
        for (std::size_t i = 0; i < std::min(tr45.sig.size(), trX.sig.size()); ++i) {
            dsig = std::max(dsig, (tr45.sig[i] - Rz * trX.sig[i] * Rz.transpose()).norm());
            smax = std::max(smax, trX.sig[i].norm());
        }
        Eigen::Matrix3d R = fcmFrame(sf.fcm.eul);
        const double r2d = 180.0 / M_PI;
        note("(c) 45 deg : eul = (" + std::to_string(sf.fcm.eul[0] * r2d) + ", " + std::to_string(sf.fcm.eul[1] * r2d)
             + ", " + std::to_string(sf.fcm.eul[2] * r2d) + ") deg, n1 = (" + std::to_string(R(0, 0)) + ", "
             + std::to_string(R(1, 0)) + ", " + std::to_string(R(2, 0)) + "), d1 " + std::to_string(sf.fcm.d[0])
             + ", " + std::to_string(tr45.sig.size()) + " pas");
        check("(c) |n1 . (1,1,0)/sqrt2|", std::abs(R.col(0).dot(n45)), 1.0, 1.0e-9, false);
        double a0 = std::fmod(std::abs(sf.fcm.eul[0] * r2d), 180.0);
        check("(c) |eul1| mod 180 deg (45 ou 135)", std::min(std::abs(a0 - 45.0), std::abs(a0 - 135.0)), 0.0, 1.0e-6, false);
        // n2, n3 = base quelconque du plan propre degenere (z, (1,-1,0)) : eul2, eul3
        // ne sont pas contraints ; seul n1 (et eul1) le sont
        check("(c) |n1_z| (fissure dans le plan xy)", std::abs(R(2, 0)), 0.0, 1.0e-9, false);
        check("(c) angle de n1 dans xy = atan2(n1_y, n1_x) mod 180 deg (45 ou -135 -> 45)",
              std::fmod(std::abs(std::atan2(R(1, 0), R(0, 0)) * r2d) + 180.0, 180.0), 45.0, 1.0e-6, false);
        check("(c) objectivite : max |sigma_45 - Rz sigma_x Rz^T| / max |sigma_x|", dsig / smax, 0.0, 1.0e-8, false);
    }

    // ---- (d) rotation des directions principales apres amorcage -----------
    {
        // (d1) tel qu'ecrit : eps_xx maintenu, eps_yy croissant, sigma_zz = 0
        auto S = mk(""), F = mk("tensionDamage = fixed\n");
        MatState ss, sf;
        double eS = pull(*S, ss, ex, 0.5, "d1_pull_x", "scalar", nullptr);
        double eF = pull(*F, sf, ex, 0.5, "d1_pull_x", "fixed", nullptr);
        auto epsD1 = [&](double exx, double eyy) {
            Eigen::Matrix3d e = Eigen::Matrix3d::Zero();
            e(0, 0) = exx; e(1, 1) = eyy; e(2, 2) = -nu / (1.0 - nu) * (exx + eyy);
            return e;
        };
        double sxyMaxS = 0.0, sxyMaxF = 0.0, eyyInit = -1.0, S22Init = 0.0;
        int nPrev = sf.fcm.nAct;
        Eigen::Matrix3d lastS, lastF;
        const int nD1 = 400;
        const double tolD1 = (m.E * nu / ((1 + nu) * (1 - 2 * nu)) + 2.0 * m.G()) * (3.0 + nu) * eF / nD1 / m.ft;
        ramp(*S, ss, epsD1(eS, -nu * eS), epsD1(eS, 3.0 * eS), nD1, "d1_pull_y", "scalar",
             [&](int, const Eigen::Matrix3d&, const Eigen::Matrix3d& sg, const MatState&) {
                 sxyMaxS = std::max(sxyMaxS, std::abs(sg(0, 1))); lastS = sg; });
        ramp(*F, sf, epsD1(eF, -nu * eF), epsD1(eF, 3.0 * eF), nD1, "d1_pull_y", "fixed",
             [&](int, const Eigen::Matrix3d& e, const Eigen::Matrix3d& sg, const MatState& st) {
                 sxyMaxF = std::max(sxyMaxF, std::abs(sg(0, 1))); lastF = sg;
                 if (st.fcm.nAct > nPrev && eyyInit < 0.0) {
                     eyyInit = e(1, 1);
                     // S_22 effectif = lambda tr + 2G eps_yy
                     double lam = m.E * nu / ((1 + nu) * (1 - 2 * nu)), G = m.G();
                     S22Init = lam * e.trace() + 2.0 * G * e(1, 1);
                 }
                 nPrev = st.fcm.nAct; });
        Eigen::Matrix3d R = fcmFrame(sf.fcm.eul);
        std::cout << "[fixed] (d1) DONNEE tel qu'ecrit (x maintenu a eps_xx = " << eF << ", eps_yy -> " << 3.0 * eF
                  << ", sigma_zz_eff = 0) : |sigma_xy| max scalar " << sxyMaxS << " Pa, fixed " << sxyMaxF
                  << " Pa (axes = repere fige : nul par symetrie) ; fixed : y s'amorce a eps_yy = " << eyyInit
                  << " (S_22 effectif " << S22Init << " Pa, ft " << m.ft << "), nAct fin " << sf.fcm.nAct
                  << ", |n2.y| = " << std::abs(R.col(1).dot(ey)) << ", d1 " << sf.fcm.d[0] << " -> d2 " << sf.fcm.d[1]
                  << " ; sigma_yy nominal fin fixed " << lastF(1, 1) << " Pa / scalar " << lastS(1, 1)
                  << " Pa (D scalar fin " << ss.D << ")\n";
        check("(d1) FIXED : |sigma_xy| max / ft (symetrie)", sxyMaxF / m.ft, 0.0, 1.0e-9, false);
        check("(d1) FIXED : amorcage de y a S_22 = ft (a un pas pres)", S22Init / m.ft, 1.0, tolD1 + 1e-12, true);
        check("(d1) FIXED : |n2 . y|", std::abs(R.col(1).dot(ey)), 1.0, 1.0e-9, false);

        // (d2) rotation : x maintenu (d = 0,5), traction croissante selon 45 deg
        for (double beta : {1.0, 0.1}) {
            std::string extra = "tensionDamage = fixed\n";
            if (beta != 1.0) extra += "tensionShearRetention = " + std::to_string(beta) + "\n";
            auto Fb = mk(extra);
            MatState sb;
            double eb = pull(*Fb, sb, ex, 0.5, "d2_pull_x", beta == 1.0 ? "fixed_b1" : "fixed_b01", nullptr);
            const Eigen::Matrix3d e0 = uni(eb, ex);
            double lam = m.E * nu / ((1 + nu) * (1 - 2 * nu)), G = m.G();
            Eigen::Matrix3d last, eLast;
            std::vector<double> ratio;
            ramp(*Fb, sb, e0, e0 + uni(2.0 * eb, n45), 100, "d2_pull_45", beta == 1.0 ? "fixed_b1" : "fixed_b01",
                 [&](int, const Eigen::Matrix3d& e, const Eigen::Matrix3d& sg, const MatState&) {
                     last = sg; eLast = e;
                     double sxyEff = 2.0 * G * e(0, 1);
                     if (std::abs(sxyEff) > 1e3) ratio.push_back(sg(0, 1) / sxyEff); });
            Eigen::Matrix3d sigEff = lam * eLast.trace() * I + 2.0 * G * eLast;
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> en(last), ee(sigEff);
            Eigen::Vector3d pn = en.eigenvectors().col(2), pe = ee.eigenvectors().col(2);
            double mis = std::acos(std::min(1.0, std::abs(pn.dot(pe)))) * 180.0 / M_PI;
            double angN = std::atan2(pn(1), pn(0)) * 180.0 / M_PI, angE = std::atan2(pe(1), pe(0)) * 180.0 / M_PI;
            Eigen::Matrix3d Rb = fcmFrame(sb.fcm.eul);
            std::cout << "[fixed] (d2) DONNEE beta = " << beta << " : x maintenu a eps_xx = " << eb
                      << " (d1 = 0,5), + traction 45 deg jusqu'a e' = " << 2.0 * eb
                      << " : sigma_xy nominal fin " << last(0, 1) << " Pa / effectif " << 2.0 * G * eLast(0, 1)
                      << " Pa = " << last(0, 1) / (2.0 * G * eLast(0, 1)) << " (1 - beta max d^ouv = "
                      << 1.0 - beta * std::max(sb.fcm.d[0], sb.fcm.d[1]) << ") ; ratio min/max sur la rampe "
                      << *std::min_element(ratio.begin(), ratio.end()) << " / " << *std::max_element(ratio.begin(), ratio.end())
                      << " ; direction principale nominale " << angN << " deg vs effective " << angE
                      << " deg (desalignement " << mis << " deg) ; d1 fin " << sb.fcm.d[0] << ", d2 " << sb.fcm.d[1]
                      << " (nAct " << sb.fcm.nAct << ", |n2.y| " << std::abs(Rb.col(1).dot(ey)) << ") ; sigma fin xx/yy/xy "
                      << last(0, 0) << " / " << last(1, 1) << " / " << last(0, 1) << " Pa\n";
        }
        // scalaire sur le meme trajet (d2), pour comparaison
        {
            auto S2 = mk("");
            MatState s2;
            double eb = pull(*S2, s2, ex, 0.5, "d2_pull_x", "scalar", nullptr);
            const Eigen::Matrix3d e0 = uni(eb, ex);
            Eigen::Matrix3d last, eLast;
            ramp(*S2, s2, e0, e0 + uni(2.0 * eb, n45), 100, "d2_pull_45", "scalar",
                 [&](int, const Eigen::Matrix3d& e, const Eigen::Matrix3d& sg, const MatState&) { last = sg; eLast = e; });
            double G = m.G();
            std::cout << "[fixed] (d2) DONNEE scalar : sigma_xy nominal fin " << last(0, 1) << " Pa / effectif "
                      << 2.0 * G * eLast(0, 1) << " Pa = " << last(0, 1) / (2.0 * G * eLast(0, 1)) << " (= 1 - D sur la partie tendue, D fin "
                      << s2.D << ") ; sigma fin xx/yy " << last(0, 0) << " / " << last(1, 1) << " Pa\n";
        }
    }

    // ---- (e) energie dissipee en traction uniaxiale ----------------------
    for (int mod = 0; mod < 2; ++mod) {
        auto L = mk(mod == 0 ? "" : "tensionDamage = fixed\n");
        const char* tag = mod == 0 ? "scalar" : "fixed";
        MatState st;
        double e = 0.0, W = 0.0, sPrev = 0.0, sMax = 0.0;
        int step = 0;
        while (e < 0.5) {
            e += de;
            Eigen::Matrix3d eps = uni(e, ex);
            Eigen::Matrix3d sg = L->stress(eps, st, dt, lc);
            W += 0.5 * (sg(0, 0) + sPrev) * de;
            sPrev = sg(0, 0);
            sMax = std::max(sMax, sPrev);
            if ((step++ % 50) == 0) rec("e_energy", tag, step, eps, sg, st);
            if (st.D > 0.5 && sPrev < 1.0e-4 * m.ft) break;
        }
        std::cout << "[fixed] (e) " << tag << " : pic " << sMax << " Pa (ft " << m.ft << "), arret a eps_xx = " << e
                  << " (sigma " << sPrev << " Pa, D " << st.D << ", " << step << " pas)\n";
        check(std::string("(e) ") + tag + " : int sigma_xx deps_xx x lc / Gf", W * lc / m.Gf, 1.0, 0.02, true);
        check(std::string("(e) ") + tag + " : wDamT x lc / Gf", st.wDamT * lc / m.Gf, 1.0, 0.02, true);
        check(std::string("(e) ") + tag + " : pic / ft", sMax / m.ft, 1.0, de / k0 + 1e-12, true);
    }

    std::cout << "[fixed] bilan : " << nPass << " PASS, " << nFail << " FAIL\n";
    return nFail == 0 ? 0 : 1;
}

std::unique_ptr<MatLaw> MatLaw::make(const std::string& kind,
                                     const Material& m, const Config& c,
                                     double lcMax) {
    std::unique_ptr<MatLaw> law;
    double erodeD = c.getd("erodeD", 0.98);
    double erodeEpv = c.getd("erodeEpv", 1.5);
    // portee de l'heterogeneite (voir MatLaw.hpp) — validee ici pour que la
    // faute de frappe soit signalee, et non ignoree en silence
    std::string wsc = c.gets("weibullScope", "strength");
    if (wsc != "strength" && wsc != "strengthGf" && wsc != "lcz")
        throw std::runtime_error("weibullScope must be strength | strengthGf | lcz");
    // ---- briques opt-in du noyau dpr/saksala (2026-09-03) ----------------
    // Toutes les cles absentes => BrickOpts par defaut => noyau historique.
    BrickOpts br;
    {
        std::string mer = c.gets("meridian", "linear");
        if (mer != "linear" && mer != "power")
            throw std::runtime_error("meridian must be linear | power (got '"
                                     + mer + "')");
        br.power = mer == "power";
        br.B = c.getd("merB", 56.59);
        br.n = c.getd("merN", 0.538);
        if (br.power && (!(br.B > 0.0) || !(br.n > 0.0 && br.n <= 1.0)))
            throw std::runtime_error("meridian = power requires merB > 0 "
                                     "[MPa^(1-n)] and 0 < merN <= 1");
        // UCS du cone lineaire = 2 c cos(phi) / (1 - sin(phi)) : continuite
        // en sigma3 = 0 ; merFc0 permet de la surcharger explicitement
        double sf = std::sin(m.phiDeg * M_PI / 180.0);
        double cf = std::cos(m.phiDeg * M_PI / 180.0);
        br.fc0 = c.getd("merFc0", 2.0 * m.cohesion * cf / (1.0 - sf));
        std::string cd = c.gets("compDamage", "none");
        if (cd != "none" && cd != "crackband")
            throw std::runtime_error("compDamage must be none | crackband "
                                     "(got '" + cd + "')");
        br.compDam = cd == "crackband";
        br.Ac = c.getd("compAc", 0.98);
        br.GIIc = c.getd("compGIIc", 1.0e4);           // [J/m^2] = 10 N/mm
        if (br.compDam && (!(br.Ac > 0.0 && br.Ac <= 1.0) || !(br.GIIc > 0.0)))
            throw std::runtime_error("compDamage = crackband requires "
                                     "0 < compAc <= 1 and compGIIc > 0 [J/m^2]");
        br.erodeDc = c.getd("erodeDc", 0.0);
        if (br.erodeDc < 0.0 || br.erodeDc > 1.0)
            throw std::runtime_error("erodeDc is a NORMALISED threshold on "
                                     "omega_c / compAc: 0 (off) < erodeDc <= 1");
        if (br.erodeDc > 0.0 && !br.compDam && kind != "cdp")
            throw std::runtime_error("erodeDc > 0 requires compDamage = "
                                     "crackband (no omega_c otherwise)");
        br.erodeWfrac = c.getd("erodeWfrac", 0.0);
        if (br.erodeWfrac < 0.0 || br.erodeWfrac > 1.0)
            throw std::runtime_error("erodeWfrac is a fraction of the band "
                                     "energy Gf/lc: 0 (off) < erodeWfrac <= 1");
        br.apex = c.getb("dpApex", false);
        {
            std::string dt = c.gets("dpTension", "on");
            if (dt != "on" && dt != "off")
                throw std::runtime_error("dpTension must be on | off");
            br.tensionOff = dt == "off";
            std::string rd = c.gets("rankineDrive", "strain");
            if (rd != "strain" && rd != "stress")
                throw std::runtime_error("rankineDrive must be strain | stress");
            br.rankineStress = rd == "stress";
        }
        // ---- tensionDamage = scalar (defaut) | fixed (2026-09-05, w19) ----
        {
            std::string td = c.gets("tensionDamage", "scalar");
            if (td != "scalar" && td != "fixed")
                throw std::runtime_error("tensionDamage must be scalar | fixed "
                                         "(got '" + td + "')");
            br.fixedCrack = td == "fixed";
            if (br.fixedCrack && kind != "dpr" && kind != "saksala")
                throw std::runtime_error("tensionDamage = fixed is implemented for "
                                         "law = dpr | saksala only (saksala2011 : port "
                                         "fidele, endommagement pilote par la deformation "
                                         "viscoplastique, pas de bande (ft, Gf, lc) ; cdp : "
                                         "tables propres ; dpdfh : deja directionnel)");
            br.shearRet = c.getd("tensionShearRetention", 1.0);
            if (c.has("tensionShearRetention") && !br.fixedCrack)
                throw std::runtime_error("tensionShearRetention requires "
                                         "tensionDamage = fixed");
            if (!(br.shearRet >= 0.0 && br.shearRet <= 1.0))
                throw std::runtime_error("tensionShearRetention (beta) must be in "
                                         "[0, 1] (1 = no shear retention)");
        }
    }

    if (kind == "elastic") {
        law = std::make_unique<ElasticLaw>(m);
    } else if (kind == "dpr") {
        // dprCap = true : capP0 / capH sont lus AUSSI par dpr (la compaction
        // volumique de saksala sans sa viscosite). Derriere une cle
        // d'activation parce que des decks du depot (fem3d_percussion_dpr.cfg)
        // portent capP0 en heritage du deck saksala : ils doivent rester
        // bit-identiques. Defaut false = dpr sans cap, comme avant.
        double p0 = 0.0, H = 0.0;
        if (c.getb("dprCap", false)) {
            p0 = c.getd("capP0", 0.0);
            H = c.getd("capH", m.K());
            if (p0 <= 0.0)
                throw std::runtime_error("dprCap = true requires capP0 > 0");
        }
        law = std::make_unique<PlasticDamageLaw>(m, 0.0, p0, H, erodeD,
                                                 erodeEpv, br);
    } else if (kind == "mc") {
        // Mohr-Coulomb elasto-plastique de Ye et al. (IJRMMS 194, 2025) :
        // c, phi du bloc materiau par defaut, dilatance psi = 0 (non associe)
        double coh = c.getd("mcCohesion", m.cohesion);
        double phi = c.getd("mcFrictionDeg", m.phiDeg);
        double psi = c.getd("mcDilationDeg", 0.0);
        if (!(coh > 0.0))
            throw std::runtime_error("mcCohesion must be > 0");
        if (!(phi >= 0.0 && phi < 89.0))
            throw std::runtime_error("mcFrictionDeg must be in [0, 89)");
        if (!(psi >= 0.0 && psi <= phi))
            throw std::runtime_error("mcDilationDeg must be in [0, "
                                     "mcFrictionDeg] (associe si psi = phi)");
        law = std::make_unique<MohrCoulombLaw>(m, coh, phi, psi);
    } else if (kind == "saksala") {
        double eta = c.getd("saksalaEta", 0.05e6);     // [Pa s]
        if (!(eta > 0.0))
            throw std::runtime_error("saksalaEta must be > 0 (use law = dpr "
                                     "for the rate-independent limit)");
        double p0 = c.getd("capP0", 8.0 * m.cohesion);
        double H = c.getd("capH", m.K());
        law = std::make_unique<PlasticDamageLaw>(m, eta, p0, H, erodeD,
                                                 erodeEpv, br);
    } else if (kind == "saksala2011") {
        // defaults = Table I of Saksala (2011), converted from MPa to Pa;
        // E, nu, phi, ft, c come from the shared Material block
        sk11::Props p;
        p.young = m.E;
        p.nu = m.nu;
        p.phiDeg = m.phiDeg;
        p.ft0 = m.ft;
        p.c0 = m.cohesion;
        p.betaDP = c.getd("skBetaDP", 0.0346);
        p.cres0 = c.getd("skCres", 2.89e6);
        p.hdp0 = c.getd("skHdp", -10.0e9);
        p.sdp = c.getd("skSdp", 1.0e4);
        p.smr = c.getd("skSmr", 1.0e4);
        p.at = c.getd("skAt", 0.98);
        p.betaT = c.getd("skBetaT", 5000.0);
        p.pp0 = c.getd("skPp0", 1040.0e6);
        p.ptr0 = c.getd("skPtr0", 377.0e6);
        p.dcap = c.getd("skDcap", 1.0e-9);
        p.wcap = c.getd("skWcap", 0.0433);
        p.nd = c.getd("skNd", 7.5e-8);
        // the VUMAT's failed = 3 input validation, as a hard error
        if (!(p.pp0 > p.ptr0) || !(p.ptr0 > 0.0) || !(p.dcap > 0.0)
            || !(p.wcap > 0.0) || !(p.at >= 0.0 && p.at <= 1.0))
            throw std::runtime_error("saksala2011: requires skPp0 > skPtr0 "
                                     "> 0, skDcap > 0, skWcap > 0, skAt in "
                                     "[0, 1]");
        law = std::make_unique<Saksala2011Law>(m, p);
    } else if (kind == "dpdfh") {
        // defaults = the thesis' Red Bohus DP-DFH card (phd/CONTINUUM.md §2,
        // constants=10, converted to SI); E, nu, rho come from the shared
        // Material block. dfhS default = the CARD's 4.18879 (= 4 pi / 3);
        // the VUMAT's internal <=0 fallback is 3.74 — set dfhS explicitly
        // to reproduce a run that relied on that fallback.
        dfhk::Props p;
        p.E = m.E;
        p.nu = m.nu;
        p.rho = m.rho;
        p.betaDeg = c.getd("dfhBetaDeg", 51.7);
        p.dcoh = c.getd("dfhDCoh", 153.3e6);
        p.psiDeg = c.getd("dfhPsiDeg", 15.0);
        p.m = c.getd("dfhWeibullM", 24.0);
        p.sigw = c.getd("dfhSigW", 120.0e6);
        p.zeff = c.getd("dfhZeff", 1.0e-9);            // 1 mm^3 in SI
        p.k = c.getd("dfhK", 0.38);
        p.S = c.getd("dfhS", 4.18879);
        p.deld = c.getd("dfhDeld", 1.0e9);             // deletion OFF
        // dilatance variable psi(pbar) — la forme de vumat_hole.f. OPT-IN :
        // absente => psiDeg fixe, trajectoires inchangees au bit pres.
        p.psiVar = c.getb("dfhPsiVar", false);
        p.psi0 = c.getd("dfhPsi0", 160.345);
        p.kpsi = c.getd("dfhKPsi", 0.213793);
        p.psiMax = c.getd("dfhPsiMax", 51.7);
        if (!(p.m > 1.0) || !(p.sigw > 0.0) || !(p.zeff > 0.0)
            || !(p.k > 0.0) || !(p.S > 0.0) || !(p.dcoh > 0.0)
            || !(p.betaDeg > 0.0 && p.betaDeg < 89.0))
            throw std::runtime_error("dpdfh: requires dfhWeibullM > 1, "
                                     "dfhSigW/dfhZeff/dfhK/dfhS/dfhDCoh > 0, "
                                     "dfhBetaDeg in (0, 89)");
        law = std::make_unique<DpDfhLaw>(m, p);
    } else if (kind == "dfhplus") {
        // DFH+ etape 1 (2026-09-06) : MEME perimetre physique que dpdfh mais
        // batie sur une ENERGIE LIBRE POSTULEE (decomposition spectrale de la
        // deformation elastique, sigma = d(rho psi)/d eps, Y_i = -d(rho psi)/
        // dD_i). Loi SEPAREE, dans src/MatLawDfhPlus.cpp : DpDfhLaw ci-dessus
        // n'est pas touchee. Memes cles dfh* (cartes interchangeables) + les
        // cles propres au cadre dfhpPsiClamp / dfhpVolInteg.
        law = makeDfhPlusLaw(m, c);
    } else if (kind == "cdp") {
        // Concrete Damaged Plasticity (2026-09-04) : cles cdp*, erosion par
        // les cles existantes (erodeD, erodeEpv, erodeDc, erodeWfrac) ;
        // gardes de snap-back propres (G1-G3) dans le constructeur — la garde
        // crack band de dpr/saksala ci-dessous ne s'applique pas
        law = std::make_unique<CdpLaw>(m, c, lcMax, erodeD, erodeEpv,
                                       br.erodeDc, br.erodeWfrac);
    } else {
        throw std::runtime_error("law must be elastic | dpr | mc | saksala | "
                                 "saksala2011 | dpdfh | dfhplus | cdp (got '"
                                 + kind + "')");
    }
    // crack-band feasibility check at the coarsest element (dpr/saksala
    // only: saksala2011 deliberately has NO fracture-energy regularization,
    // like the published law — its results are mesh-sensitive by design)
    if (kind == "dpr" || kind == "saksala") {
        double k0 = m.ft / m.E;
        double kf = m.Gf / (lcMax * m.ft) - 0.5 * k0;
        if (kf <= 0.05 * k0)
            throw std::runtime_error(
                "crack band: largest element (" + std::to_string(lcMax)
                + " m) exceeds E Gf / ft^2 = "
                + std::to_string(m.E * m.Gf / (m.ft * m.ft))
                + " m — refine the mesh or raise Gf");
    }
    if (law) law->wScaleGf_ = (wsc == "strengthGf");
    return law;
}

int saksala2011Selftest(const std::string& csvPath) {
    // Table I of Saksala (2011) in Pa; identical strain paths and dt to the
    // Fortran driver (stresses reported in MPa in the CSV for direct diff)
    sk11::Props P;
    P.young = 60.0e9;  P.nu = 0.2;      P.phiDeg = 30.0;  P.betaDP = 0.0346;
    P.ft0 = 13.0e6;    P.c0 = 37.5e6;   P.cres0 = 2.89e6; P.hdp0 = -10.0e9;
    P.sdp = 1.0e4;     P.smr = 1.0e4;   P.at = 0.98;      P.betaT = 5000.0;
    P.pp0 = 1040.0e6;  P.ptr0 = 377.0e6;
    P.dcap = 1.0e-9;   P.wcap = 0.0433; P.nd = 7.5e-8;
    double dt = 1.0e-6;

    std::ofstream out(csvPath);
    out << "test,step,s11,s22,s33,s12,s23,s13,omega,kapDP,kapMR,epsv,pp,"
           "active\n";
    out.precision(15);
    out << std::scientific;

    struct Path { int id, steps; std::array<double, 3> d; };
    const Path paths[3] = {
        {1, 1500, {2.0e-6, -0.2 * 2.0e-6, -0.2 * 2.0e-6}},
        {2, 1500, {-2.0e-6, 0.2 * 2.0e-6, 0.2 * 2.0e-6}},
        {3, 2500, {-2.0e-5, -2.0e-5, -2.0e-5}},
    };
    for (const Path& pth : paths) {
        MatState::Sk11 st;
        st.init = true;
        st.pp = P.pp0;
        sk11::V6 deps = {pth.d[0], pth.d[1], pth.d[2], 0.0, 0.0, 0.0};
        sk11::V6 sig;
        for (int i = 1; i <= pth.steps; ++i) {
            sk11::updatePoint(dt, deps, P, st, sig);
            if (st.failed != 0) {
                out.close();
                throw std::runtime_error("saksala2011 selftest: local return "
                                         "mapping failed on path "
                                         + std::to_string(pth.id));
            }
            if (i % 10 == 0) {
                double omega = 1.0 - (1.0 - P.at
                                      + P.at * std::exp(-P.betaT * st.eqvt));
                omega = std::min(std::max(omega, 0.0), P.at);
                out << pth.id << ',' << i;
                for (int j = 0; j < 6; ++j) out << ',' << sig[j] / 1.0e6;
                out << ',' << omega << ',' << st.kapDP << ',' << st.kapMR
                    << ',' << st.epsv << ',' << st.pp / 1.0e6 << ','
                    << st.active << '\n';
            }
        }
    }
    return 0;
}

int dpdfhSelftest(const std::string& csvPath) {
    // Red Bohus card in SI; identical strain paths, dt, element size and
    // material-point coordinates to VUMATS/dfh/test_kstdfh.f90 (mm-MPa):
    // the physics is unit-homogeneous, the spatial hash sees the same
    // micrometre integers, so the traces must superpose (stresses reported
    // in MPa in the CSV for direct diff against the Fortran reference).
    dfhk::Props P;                       // defaults = the thesis card
    const double dt = 1.0e-7;
    const double lc = 2.0e-3;            // charLength 2 mm
    const double Vel = lc * lc * lc;
    const double x0[3] = {1.234e-3, 2.345e-3, 3.456e-3};

    std::ofstream out(csvPath);
    out << "path,step,s11,s22,s33,s12,s23,s13,D1,D2,D3,peeq,"
           "sc1,sc2,sc3,smaxh,eul1,eul2,eul3\n";
    out.precision(15);
    out << std::scientific;

    struct Path { int id, steps; std::array<double, 6> d; };
    const Path paths[4] = {
        {1, 2000, {2.0e-6, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {2, 400, {2.0e-5, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {3, 1500, {-2.4e-5, 0.9e-5, 0.9e-5, 0.0, 0.0, 0.0}},
        {4, 2500, {1.6e-6, -4.0e-7, 2.0e-7, 6.0e-7, 3.0e-7, -2.0e-7}},
    };
    for (const Path& pth : paths) {
        MatState::Dfh st;
        double sigk = P.sigw * std::pow(P.zeff / Vel, 1.0 / P.m);
        dfhk::seed(x0[0], x0[1], x0[2], sigk, P.m, st.sc);
        st.seeded = true;
        dfhk::V6 deps;
        for (int j = 0; j < 6; ++j) deps[j] = pth.d[j];
        dfhk::V6 snom;
        for (int i = 1; i <= pth.steps; ++i) {
            st.t += dt;
            dfhk::updatePoint(dt, deps, P, Vel, 1.0, st, snom);
            if (i % 10 == 0) {
                out << pth.id << ',' << i;
                for (int j = 0; j < 6; ++j) out << ',' << snom[j] / 1.0e6;
                for (int j = 0; j < 3; ++j) out << ',' << st.Dv[j];
                out << ',' << st.peeq;
                for (int j = 0; j < 3; ++j) out << ',' << st.sc[j] / 1.0e6;
                out << ',' << st.smaxh / 1.0e6;
                for (int j = 0; j < 3; ++j) out << ',' << st.eul[j];
                out << '\n';
            }
        }
    }
    return 0;
}

} // namespace rockim
