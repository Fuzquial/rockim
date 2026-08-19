// ---------------------------------------------------------------------------
// MatLaw — the fem3d constitutive laws. See the header for the model cards.
// ---------------------------------------------------------------------------
#include "rockim/MatLaw.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>

namespace rockim {

namespace {

// max principal value of a symmetric 3x3
double maxPrincipal(const Eigen::Matrix3d& S) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(S);
    return es.eigenvalues().maxCoeff();
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
class PlasticDamageLaw : public MatLaw {
public:
    PlasticDamageLaw(const Material& m, double eta, double capP0, double capH,
                     double erodeD, double erodeEpv)
        : MatLaw(m), eta_(eta), capP0_(capP0), capH_(capH), erodeD_(erodeD),
          erodeEpv_(erodeEpv) {}

    Eigen::Matrix3d stress(const Eigen::Matrix3d& eps, MatState& s,
                           double dt, double lc) const override {
        if (s.eroded) return Eigen::Matrix3d::Zero();
        if (capP0_ > 0.0 && s.pc == 0.0) s.pc = capP0_;   // lazy init

        // ---- elastic predictor on the effective (undamaged) skeleton ----
        Eigen::Matrix3d sig = elastic(eps - s.epsP);

        // ---- pressure cap (compression crush) ---------------------------
        double p = sig.trace() / 3.0;
        if (capP0_ > 0.0 && p < -s.pc) {
            double over = -(p + s.pc);                 // > 0
            double dev_ = over / (K_ + capH_);         // |volumetric return|
            s.epsP -= (dev_ / 3.0) * Eigen::Matrix3d::Identity();
            s.pc += capH_ * dev_;
            p += K_ * dev_;
            sig += K_ * dev_ * Eigen::Matrix3d::Identity();
        }

        // ---- DP cone, deviatoric (visco)plastic return ------------------
        Eigen::Matrix3d dev = sig - p * Eigen::Matrix3d::Identity();
        double sj2 = std::sqrt(0.5) * dev.norm();      // sqrt(J2)
        double F = sj2 + adp_ * 3.0 * p - kdp_;        // I1 = 3 p
        if (F > 0.0 && sj2 > 1e-12) {
            // linear Perzyna: F_{n+1} = eta dlam/dt with radial deviatoric
            // return sqrtJ2 -> sqrtJ2 - G dlam  =>  closed form
            double dlam = F / (G_ + eta_ / std::max(dt, 1e-30));
            Eigen::Matrix3d nfl = dev / (2.0 * sj2);   // flow direction
            s.epsP += dlam * nfl;
            dev *= (1.0 - G_ * dlam / sj2);
            sig = dev + p * Eigen::Matrix3d::Identity();
            s.epvEq += dlam / std::sqrt(3.0);
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
        double e1 = maxPrincipal(eps - s.epsP);        // driving strain
        if (e1 > k0 && e1 > s.kappa) s.kappa = e1;
        if (s.kappa > k0) {
            double D = 1.0 - (k0 / s.kappa)
                             * std::exp(-(s.kappa - k0) / kf);
            if (D > s.D) s.D = std::min(1.0, D);
        }

        // UNILATERAL damage: degrade the TENSILE principal components only.
        // A scalar (1-D) on the full stress also kills the COMPRESSIVE
        // bearing capacity: under an indenter the damaged surface elements
        // then collapse and the tool tunnels through the block at constant
        // velocity (measured before the fix). Cracked rock still bears
        // compression — the same reason percussion VUMATs keep damaged
        // elements alive instead of deleting them.
        Eigen::Matrix3d nominal = sig;
        if (s.D > 0.0) {
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(sig);
            Eigen::Matrix3d sigT = Eigen::Matrix3d::Zero();
            for (int q = 0; q < 3; ++q) {
                double lq = es.eigenvalues()(q);
                if (lq > 0.0)
                    sigT += lq * es.eigenvectors().col(q)
                               * es.eigenvectors().col(q).transpose();
            }
            nominal = (sig - sigT) + (1.0 - s.D) * sigT;
        }

        // erosion = material REMOVAL, only where it is physical: fully
        // damaged elements in net TENSION (spalled chips leaving the
        // surface) or over-crushed elements (comminuted material squeezed
        // out, erodeEpv on the equivalent viscoplastic strain). A damaged
        // element under compression stays: it is the rubble bed.
        bool spall = s.D >= erodeD_ && nominal.trace() >= 0.0;
        bool crush = erodeEpv_ > 0.0 && s.epvEq >= erodeEpv_;
        if (spall || crush) { s.eroded = true; return Eigen::Matrix3d::Zero(); }
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

private:
    double eta_, capP0_, capH_, erodeD_, erodeEpv_;
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
    double tanp = std::tan(P.psiDeg * M_PI / 180.0);
    double oxm = 1.0 / P.m;

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
        if (st.dead) s.eroded = true;
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

std::unique_ptr<MatLaw> MatLaw::make(const std::string& kind,
                                     const Material& m, const Config& c,
                                     double lcMax) {
    std::unique_ptr<MatLaw> law;
    double erodeD = c.getd("erodeD", 0.98);
    double erodeEpv = c.getd("erodeEpv", 1.5);
    // portee de l'heterogeneite (voir MatLaw.hpp) — validee ici pour que la
    // faute de frappe soit signalee, et non ignoree en silence
    std::string wsc = c.gets("weibullScope", "strength");
    if (wsc != "strength" && wsc != "strengthGf")
        throw std::runtime_error("weibullScope must be strength | strengthGf");
    if (kind == "elastic") {
        law = std::make_unique<ElasticLaw>(m);
    } else if (kind == "dpr") {
        law = std::make_unique<PlasticDamageLaw>(m, 0.0, 0.0, 0.0, erodeD,
                                                 erodeEpv);
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
                                                 erodeEpv);
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
        if (!(p.m > 1.0) || !(p.sigw > 0.0) || !(p.zeff > 0.0)
            || !(p.k > 0.0) || !(p.S > 0.0) || !(p.dcoh > 0.0)
            || !(p.betaDeg > 0.0 && p.betaDeg < 89.0))
            throw std::runtime_error("dpdfh: requires dfhWeibullM > 1, "
                                     "dfhSigW/dfhZeff/dfhK/dfhS/dfhDCoh > 0, "
                                     "dfhBetaDeg in (0, 89)");
        law = std::make_unique<DpDfhLaw>(m, p);
    } else {
        throw std::runtime_error("law must be elastic | dpr | saksala | "
                                 "saksala2011 | dpdfh (got '" + kind + "')");
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
