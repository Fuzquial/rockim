// ---------------------------------------------------------------------------
// bench_energie_matrice.cpp — banc FALSIFIANT d'energie au point materiel.
//
// Question : la loi de la note de septembre 2026 (law = saksala + meridian =
// power + mhForm = principal + dpFlowForm = mc + compDamage = crackband +
// bulkTensionDamage = off + saksalaEta = 1,5e6 + capP0 = 0 + psi = 5 deg)
// CREE-T-ELLE de l'energie sur un cycle de deformation FERME ?
//
// Critere : W_cycle = oint sigma : d eps sur un cycle ferme. Une loi
// dissipative rend W_cycle >= 0 sur 3 cycles identiques successifs ;
// W_cycle < -1e-6 * max(Psi_el) = energie CREEE -> drapeau CREE (*).
//
// Pilote autonome : n'appelle que MatLaw::make et MatLaw::stress (meme recette
// que obj_LOTA/banc_lotA.cpp). Ne modifie RIEN dans src/ ni include/.
//
//   build : tests_f2\build_bench_energie.cmd
//   run   : obj_ENER\bench_energie.exe [rateBase=1e3] [nSeeds=20] [dir=obj_ENER]
// ---------------------------------------------------------------------------
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "rockim/Config.hpp"
#include "rockim/MatLaw.hpp"
#include "rockim/Material.hpp"

using namespace rockim;
using M3 = Eigen::Matrix3d;

// ---- carte du deck configs/loi_note_2026.cfg (SI) ---------------------------
static const double E_ = 60e9, NU_ = 0.2;
static const double LC_ = 1.5e-3;     // m, taille d'element du run
static const double DT_ = 1.6e-9;     // s, pas de temps du run

static Material card() {
    Material m;
    m.rho = 2630.0; m.E = E_; m.nu = NU_;
    m.ft = 11.4e6; m.cohesion = 30e6; m.phiDeg = 35.0; m.Gf = 100.0;
    m.gfShearFactor = 10.0;
    return m;
}

static const char* MATBLOCK =
    "E = 60e9\nnu = 0.2\nrho = 2630\nft = 11.4e6\ncohesion = 30e6\n"
    "frictionDeg = 35\nGf = 100\n";

// La loi complete de la note, telle que le deck la pose (cles de MatLaw seules)
static const char* FULL =
    "law = saksala\nmeridian = power\nmerB = 30.1\nmerN = 0.676\n"
    "merFc0 = 235e6\nsaksalaEta = 1.5e6\ncapP0 = 0\nmhForm = principal\n"
    "dpFlowForm = mc\ndpDilationDeg = 5\ncompDamage = crackband\n"
    "compAc = 0.98\ncompGIIc = 1.0e4\nbulkTensionDamage = off\n";

struct Variant { std::string tag; std::string body; };

static std::string rep(std::string s, const std::string& from, const std::string& to) {
    size_t p = s.find(from);
    if (p != std::string::npos) s.replace(p, from.size(), to);
    return s;
}

static std::vector<Variant> variants() {
    std::vector<Variant> v;
    v.push_back({"FULL", FULL});
    v.push_back({"-mhPrinc(dp)", rep(FULL, "mhForm = principal", "mhForm = dp")});
    v.push_back({"-flowMC(dp)", rep(FULL, "dpFlowForm = mc", "dpFlowForm = dp")});
    v.push_back({"-compDam", rep(FULL, "compDamage = crackband", "compDamage = none")});
    v.push_back({"+Rankine(on)", rep(FULL, "bulkTensionDamage = off", "bulkTensionDamage = on")});
    v.push_back({"-eta(dpr)", rep(rep(FULL, "law = saksala", "law = dpr"), "saksalaEta = 1.5e6\n", "")});
    v.push_back({"-psi(0)", rep(FULL, "dpDilationDeg = 5", "dpDilationDeg = 0")});
    v.push_back({"principal SEUL", rep(rep(FULL, "dpFlowForm = mc", "dpFlowForm = dp"),
                                        "compDamage = crackband", "compDamage = none")});
    v.push_back({"flowMC SEUL", rep(rep(FULL, "mhForm = principal", "mhForm = dp"),
                                     "compDamage = crackband", "compDamage = none")});
    v.push_back({"TEMOIN elastic", "law = elastic\n"});
    v.push_back({"TEMOIN g0", "law = dpr\nmeridian = power\nmerB = 30.1\nmerN = 0.676\n"
                              "merFc0 = 235e6\ncompDamage = crackband\ncompAc = 0.98\n"
                              "compGIIc = 1.0e4\n"});
    return v;
}

static std::unique_ptr<MatLaw> mk(const Variant& v, const std::string& dir) {
    static int n = 0;
    const std::string path = dir + "/_cfg_" + std::to_string(n++) + ".cfg";
    {
        std::ofstream o(path);
        // pas d'erosion : on juge la LOI seule (erodeD = 0,98 et erodeEpv =
        // 1,5 sont armes par defaut dans MatLaw::make, meme en FDEM)
        o << MATBLOCK << v.body << "erodeEpv = 0\nerodeD = 2\n";
    }
    Config c = Config::load(path);
    return MatLaw::make(c.gets("law", "dpr"), card(), c, LC_);
}

// ---- trajets : eps(u), u dans [0,1[, eps(0) = eps(1) --------------------------
struct Traj {
    std::string tag;
    std::function<M3(double)> eps;   // deformation au parametre u
    double rate;                     // /s, norme de eps_point (Frobenius)
};

static M3 sym(double e11, double e22, double e33, double e12, double e13, double e23) {
    M3 e; e << e11, e12, e13, e12, e22, e23, e13, e23, e33; return e;
}
static double tri(double u) { return u < 0.5 ? 2.0 * u : 2.0 * (1.0 - u); }

// Deformation elastique de compression UNIAXIALE EN CONTRAINTE (lateral libre)
static M3 uniaxStress(double a) { return sym(-a, NU_ * a, NU_ * a, 0, 0, 0); }
// Deviatorique tournant : a [cos th (e1e1 - e2e2) + sin th (e1e2 + e2e1)]
static M3 rotDev(double a, double th) {
    return sym(a * std::cos(th), -a * std::cos(th), 0.0, a * std::sin(th), 0, 0);
}

struct RandPath {
    std::vector<M3> knots;           // ferme : knots.front() == knots.back() == 0
    M3 at(double u) const {
        const double x = u * (knots.size() - 1);
        size_t i = (size_t)std::floor(x);
        if (i >= knots.size() - 1) return knots.back();
        const double t = x - i;
        return (1.0 - t) * knots[i] + t * knots[i + 1];
    }
};
static RandPath randPath(unsigned seed, double amp, int nStep) {
    std::mt19937 g(seed);
    std::normal_distribution<double> N(0.0, amp / std::sqrt((double)nStep));
    RandPath p;
    M3 e = M3::Zero();
    p.knots.push_back(e);
    for (int k = 0; k < nStep; ++k) {
        e += sym(N(g), N(g), N(g), 0.7 * N(g), 0.7 * N(g), 0.7 * N(g));
        p.knots.push_back(e);
    }
    const M3 e0 = e;                 // retour LINEAIRE au depart
    for (int k = 1; k <= nStep / 2; ++k)
        p.knots.push_back(e0 * (1.0 - (double)k / (nStep / 2)));
    p.knots.back() = M3::Zero();
    return p;
}

static std::string pct(double a) {
    char b[32]; std::snprintf(b, sizeof b, "%g%%", a * 100.0); return b;
}

static std::vector<Traj> trajectories(double rb, int nSeeds) {
    std::vector<Traj> T;
    // (i) compression uniaxiale puis retour (uniaxial-contrainte), 1, 3, 10 %
    for (double a : {0.01, 0.03, 0.10})
        T.push_back({"(i) uniax-sig " + pct(a),
                     [a](double u) -> M3 { return uniaxStress(a * tri(u)); }, rb});
    // (i-b) compression uniaxiale EN DEFORMATION (lateral bloque), 1 et 3 %
    for (double a : {0.01, 0.03})
        T.push_back({"(i-b) uniax-eps " + pct(a),
                     [a](double u) -> M3 { return sym(-a * tri(u), 0, 0, 0, 0, 0); }, rb});
    // (ii) compression -> traction -> compression (inversion), 1 %
    T.push_back({"(ii) comp/trac 1%",
                 [](double u) -> M3 { return uniaxStress(0.01 * std::sin(2 * M_PI * u)); }, rb});
    // (iii) cisaillement pur alterne, eps12 = 1 % (gamma = 2 %)
    T.push_back({"(iii) cisaill. alt 1%",
                 [](double u) -> M3 { return sym(0, 0, 0, 0.01 * std::sin(2 * M_PI * u), 0, 0); }, rb});
    // (iv) rotation des axes principaux, deviatorique pur, 1, 2, 5 %
    for (double a : {0.01, 0.02, 0.05})
        T.push_back({"(iv) rotDev " + pct(a),
                     [a](double u) -> M3 { return rotDev(a, 2 * M_PI * u); }, rb});
    // (iv-b) rotation SOUS CONFINEMENT modere : -0,2 % hydro + deviatorique 1 %
    T.push_back({"(iv-b) rotDev 1% hyd-.2%",
                 [](double u) -> M3 { return rotDev(0.01, 2 * M_PI * u) - 0.002 * M3::Identity(); }, rb});
    // (v) meme (iv) a 1e4 et 1e5 /s
    for (double r : {1.0e4, 1.0e5})
        T.push_back({"(v) rotDev 1% @1e" + std::to_string((int)std::log10(r)),
                     [](double u) -> M3 { return rotDev(0.01, 2 * M_PI * u); }, r});
    T.push_back({"(v) rotDev 5% @1e5",
                 [](double u) -> M3 { return rotDev(0.05, 2 * M_PI * u); }, 1.0e5});
    T.push_back({"(i) uniax-sig 3% @1e5",
                 [](double u) -> M3 { return uniaxStress(0.03 * tri(u)); }, 1.0e5});
    // (vi) marche aleatoire fermee, nSeeds graines, amplitude ~1,5 %
    for (int s = 0; s < nSeeds; ++s) {
        auto p = std::make_shared<RandPath>(randPath(1000u + s, 0.015, 60));
        T.push_back({"(vi) alea s" + std::to_string(s),
                     [p](double u) -> M3 { return p->at(u); }, rb});
    }
    // (vi-b) marche aleatoire GRANDE amplitude ~5 %, 5 graines
    for (int s = 0; s < std::min(nSeeds, 5); ++s) {
        auto p = std::make_shared<RandPath>(randPath(2000u + s, 0.05, 60));
        T.push_back({"(vi-b) alea 5% s" + std::to_string(s),
                     [p](double u) -> M3 { return p->at(u); }, rb});
    }
    return T;
}

// ---- integration du travail sur 3 cycles ------------------------------------
struct CycleRes {
    double W[3] = {0, 0, 0}, dPsi[3] = {0, 0, 0}, dDiss[3] = {0, 0, 0}, psiMax[3] = {0, 0, 0};
    double resid = 0.0;              // max ||sig - split(C(eps-epsP))|| / ||.||
    double epvEnd = 0, DcEnd = 0, DEnd = 0, DcStart = 0;
    bool cree = false, nanflag = false;
};

static double psiEff(const M3& ee) {
    const double lam = E_ * NU_ / ((1 + NU_) * (1 - 2 * NU_)), G = E_ / (2 * (1 + NU_));
    return 0.5 * lam * ee.trace() * ee.trace() + G * ee.squaredNorm();
}
static M3 elasticC(const M3& ee) {
    const double lam = E_ * NU_ / ((1 + NU_) * (1 - 2 * NU_)), G = E_ / (2 * (1 + NU_));
    return lam * ee.trace() * M3::Identity() + 2 * G * ee;
}
// la contrainte nominale ATTENDUE d'apres l'etat : (1-D) sig+ + (1-Dc) sig-
static M3 expectedNominal(const M3& eps, const MatState& s) {
    const M3 se = elasticC(eps - s.epsP);
    Eigen::SelfAdjointEigenSolver<M3> es(se);
    M3 out = M3::Zero();
    for (int q = 0; q < 3; ++q) {
        const double l = es.eigenvalues()(q);
        const M3 P = es.eigenvectors().col(q) * es.eigenvectors().col(q).transpose();
        out += (l > 0.0 ? (1.0 - s.D) : (1.0 - s.Dc)) * l * P;
    }
    return out;
}

// trace de la partie spectrale NEGATIVE de la contrainte effective C(eps-epsP)
static double trSigMinus(const M3& eps, const MatState& s) {
    Eigen::SelfAdjointEigenSolver<M3> es(elasticC(eps - s.epsP));
    double t = 0.0;
    for (int q = 0; q < 3; ++q) t += std::min(es.eigenvalues()(q), 0.0);
    return t;
}

struct Driver {
    const MatLaw& L;
    MatState st;
    M3 eps = M3::Zero(), sig = M3::Zero();
    double resid = 0.0;
    // S = somme tr(sig-) d(tr eps) : pour omega_c FIGE et epsP FIGE,
    // W_cycle predit = Dc * lambda/(2G) * S (voir en-tete du tableau 3)
    double S = 0.0;
    bool nanflag = false;
    explicit Driver(const MatLaw& l) : L(l) { st.lcComp = LC_; sig = L.stress(eps, st, DT_, LC_); }
    double step(const M3& eNew) {
        const double tm0 = trSigMinus(eps, st);
        const M3 sNew = L.stress(eNew, st, DT_, LC_);
        if (!sNew.allFinite()) nanflag = true;
        const double dW = 0.5 * (sig + sNew).cwiseProduct(eNew - eps).sum();
        const M3 ex = expectedNominal(eNew, st);
        resid = std::max(resid, (sNew - ex).norm() / std::max(ex.norm(), 1.0));
        S += 0.5 * (tm0 + trSigMinus(eNew, st)) * (eNew.trace() - eps.trace());
        sig = sNew; eps = eNew;
        return dW;
    }
    // segment lineaire eps -> eTarget a la vitesse rate
    double ramp(const M3& eTarget, double rate) {
        const M3 e0 = eps;
        const int Nr = std::max(1, (int)std::ceil((eTarget - e0).norm() / (rate * DT_)));
        double W = 0.0;
        for (int k = 1; k <= Nr; ++k) W += step(e0 + (eTarget - e0) * ((double)k / Nr));
        return W;
    }
};

static int nStepsFor(const std::function<M3(double)>& f, double rate) {
    double len = 0.0;
    M3 prev = f(0.0);
    for (int k = 1; k <= 2000; ++k) { M3 e = f((double)k / 2000); len += (e - prev).norm(); prev = e; }
    return std::max(20, (int)std::ceil(len / (rate * DT_)));
}

static CycleRes cycles(Driver& d, const std::function<M3(double)>& f, double rate) {
    CycleRes R;
    R.DcStart = d.st.Dc;
    const int N = nStepsFor(f, rate);
    d.ramp(f(0.0), rate);                         // rampe vers eps(0), hors cycle
    for (int c = 0; c < 3; ++c) {
        const double psi0 = psiEff(d.eps - d.st.epsP);
        const double d0 = d.st.wPlas + d.st.wDamC + d.st.wDamT;
        double W = 0.0, pm = psi0;
        for (int k = 1; k <= N; ++k) {
            W += d.step(f(k == N ? 0.0 : (double)k / N));
            pm = std::max(pm, psiEff(d.eps - d.st.epsP));
        }
        R.W[c] = W;
        R.dPsi[c] = psiEff(d.eps - d.st.epsP) - psi0;
        R.dDiss[c] = d.st.wPlas + d.st.wDamC + d.st.wDamT - d0;
        R.psiMax[c] = pm;
    }
    R.epvEnd = d.st.epvEq; R.DcEnd = d.st.Dc; R.DEnd = d.st.D;
    R.resid = d.resid; R.nanflag = d.nanflag;
    for (int c = 0; c < 3; ++c)
        if (R.W[c] < -1.0e-6 * std::max(R.psiMax[c], 1.0)) R.cree = true;
    return R;
}

// (vii) element PRE-BROYE (omega_c -> Ac par compression uniaxiale 3 %), puis
// decharge a un etat ELASTIQUE (sigma_eff = -100 MPa uniaxial, sous le seuil),
// puis petits cycles deviatoriques a axes TOURNANTS, dans les deux sens.
// omega_c est alors FIGE : W_cycle mesure la seule reponse elastique-endommagee.
// Une loi hyperelastique a endommagement fige rend W_cycle = 0 ; le signe de
// W_cycle s'inverse avec le sens de rotation si la reponse n'est PAS un
// gradient (raideur secante non symetrique).
static CycleRes precrushed(const MatLaw& L, double a, int dir, double rate, double preRate,
                           double* omegaOut) {
    Driver d(L);
    d.ramp(uniaxStress(0.03), preRate);
    // decharge vers eps_t = epsP + C^-1 (-100 MPa e1e1)
    const M3 eT = d.st.epsP + uniaxStress(100.0e6 / E_);
    d.ramp(eT, preRate);
    if (omegaOut) *omegaOut = d.st.Dc;
    auto f = [eT, a, dir](double u) -> M3 { return eT + rotDev(a, dir * 2 * M_PI * u); };
    return cycles(d, f, rate);
}

// (viii) le meme, mais la rotation deviatorique s'accompagne d'une RESPIRATION
// volumique dephasee : eps = eps_t + rotDev(a, 2 pi u) + b sin(2 pi u + phi0) I.
// Pre-broyage plus long (uniax 5 % a preRate) pour monter omega_c.
// Prediction analytique (omega_c et epsP FIGES, sigma_nom = sig - Dc sig-) :
//   W_cycle = Dc * lambda/(2G) * oint tr(sig-) d(tr eps)      (= Dc lam/2G * S)
// car sig- = [grad g - lambda tr(sig-) I]/(2G), g = 1/2 sum <s_i>-^2, et
// oint grad g : d eps = 0 sur un cycle ferme.
struct BreathRes { CycleRes r; double Wpred[3]; double omega; double startErr, endErr; };
static BreathRes breathing(const MatLaw& L, double a, double b, double phi0, double rate,
                           double preAmp, double preRate) {
    BreathRes B{};
    Driver d(L);
    d.ramp(uniaxStress(preAmp), preRate);
    const M3 eT = d.st.epsP + uniaxStress(100.0e6 / E_);
    d.ramp(eT, preRate);
    B.omega = d.st.Dc;
    auto f = [eT, a, b, phi0](double u) -> M3 {
        return eT + rotDev(a, 2 * M_PI * u) + b * std::sin(2 * M_PI * u + phi0) * M3::Identity();
    };
    // on refait cycles() a la main pour capturer S par cycle
    const double lam = E_ * NU_ / ((1 + NU_) * (1 - 2 * NU_)), G = E_ / (2 * (1 + NU_));
    const int N = nStepsFor(f, rate);
    d.ramp(f(0.0), rate);
    CycleRes& R = B.r;
    R.DcStart = d.st.Dc;
    B.startErr = (d.eps - f(0.0)).norm();
    for (int c = 0; c < 3; ++c) {
        const double psi0 = psiEff(d.eps - d.st.epsP);
        const double d0 = d.st.wPlas + d.st.wDamC + d.st.wDamT;
        const double S0 = d.S;
        const double Dc0 = d.st.Dc;
        double W = 0.0, pm = psi0;
        for (int k = 1; k <= N; ++k) {
            W += d.step(f(k == N ? 0.0 : (double)k / N));
            pm = std::max(pm, psiEff(d.eps - d.st.epsP));
        }
        R.W[c] = W;
        R.dPsi[c] = psiEff(d.eps - d.st.epsP) - psi0;
        R.dDiss[c] = d.st.wPlas + d.st.wDamC + d.st.wDamT - d0;
        R.psiMax[c] = pm;
        B.Wpred[c] = Dc0 * lam / (2 * G) * (d.S - S0);
        if (c == 0) B.endErr = (d.eps - f(0.0)).norm();
    }
    R.epvEnd = d.st.epvEq; R.DcEnd = d.st.Dc; R.DEnd = d.st.D;
    R.resid = d.resid; R.nanflag = d.nanflag;
    for (int c = 0; c < 3; ++c)
        if (R.W[c] < -1.0e-6 * std::max(R.psiMax[c], 1.0)) R.cree = true;
    return B;
}

int main(int argc, char** argv) {
    const double rb = argc > 1 ? std::atof(argv[1]) : 1.0e3;
    const int nSeeds = argc > 2 ? std::atoi(argv[2]) : 20;
    const std::string dir = argc > 3 ? argv[3] : "obj_ENER";
    std::printf("banc energie matrice : rate base %.3g /s, dt %.3g s, lc %.3g m, %d graines\n\n",
                rb, DT_, LC_, nSeeds);
    const auto V = variants();
    const auto T = trajectories(rb, nSeeds);
    std::vector<std::unique_ptr<MatLaw>> laws;
    for (const auto& v : V) laws.push_back(mk(v, dir));
    const bool only3 = argc > 4 && std::string(argv[4]) == "t3";   // tableau 3 seul
    int nCree = 0, nCree2 = 0;

    std::ofstream csv(dir + "/bench_energie.csv");
    csv << "traj,rate,variant,W1,W2,W3,dPsi1,dPsi2,dPsi3,dDiss1,dDiss2,dDiss3,psiMax,resid,epv,Dc,D,cree\n";

    if (!only3) {
    std::printf("=== TABLEAU 1 : trajets x variantes, cellule = min_c W_c / Psi_max (cycles 1-3) ; * = CREE, N = NaN ===\n");
    std::printf("%-26s", "");
    for (const auto& v : V) std::printf(" %14s", v.tag.substr(0, 14).c_str());
    std::printf("\n");
    std::vector<std::vector<CycleRes>> all(T.size());
    for (size_t i = 0; i < T.size(); ++i) {
        std::printf("%-26s", T[i].tag.c_str());
        for (size_t j = 0; j < V.size(); ++j) {
            Driver d(*laws[j]);
            CycleRes r = cycles(d, T[i].eps, T[i].rate);
            all[i].push_back(r);
            double mn = 1e300, pm = 1.0;
            for (int c = 0; c < 3; ++c) { mn = std::min(mn, r.W[c]); pm = std::max(pm, r.psiMax[c]); }
            std::printf(" %13.3e%s", mn / pm, r.cree ? "*" : (r.nanflag ? "N" : " "));
            csv << T[i].tag << ',' << T[i].rate << ',' << V[j].tag;
            for (int c = 0; c < 3; ++c) csv << ',' << r.W[c];
            for (int c = 0; c < 3; ++c) csv << ',' << r.dPsi[c];
            for (int c = 0; c < 3; ++c) csv << ',' << r.dDiss[c];
            csv << ',' << pm << ',' << r.resid << ',' << r.epvEnd << ',' << r.DcEnd
                << ',' << r.DEnd << ',' << (r.cree ? 1 : 0) << '\n';
        }
        std::printf("\n");
        std::fflush(stdout);
    }

    std::printf("\n=== DETAIL des cellules CREE du tableau 1 (J/m3) ===\n");
    nCree = 0;
    for (size_t i = 0; i < T.size(); ++i)
        for (size_t j = 0; j < V.size(); ++j) {
            const CycleRes& r = all[i][j];
            if (!r.cree) continue;
            ++nCree;
            std::printf("%-26s %-16s W = [%.4e %.4e %.4e]  dPsi_eff = [%.3e %.3e %.3e]  diss = [%.3e %.3e %.3e]  PsiMax %.3e  epv %.3g Dc %.3f D %.3f\n",
                        T[i].tag.c_str(), V[j].tag.c_str(), r.W[0], r.W[1], r.W[2],
                        r.dPsi[0], r.dPsi[1], r.dPsi[2], r.dDiss[0], r.dDiss[1], r.dDiss[2],
                        std::max({r.psiMax[0], r.psiMax[1], r.psiMax[2]}), r.epvEnd, r.DcEnd, r.DEnd);
        }
    if (nCree == 0) std::printf("(aucune)\n");

    std::printf("\n=== residu max ||sig_rendue - [(1-D) sig+ + (1-Dc) sig-](C(eps-epsP))|| / ||.|| par variante (doit etre ~1e-12) ===\n");
    for (size_t j = 0; j < V.size(); ++j) {
        double mx = 0.0; size_t imx = 0;
        for (size_t i = 0; i < T.size(); ++i)
            if (all[i][j].resid > mx) { mx = all[i][j].resid; imx = i; }
        std::printf("%-16s %.3e  (max sur %s ; Dc fin %.3f, D fin %.3f)\n", V[j].tag.c_str(), mx,
                    T[imx].tag.c_str(), all[imx][j].DcEnd, all[imx][j].DEnd);
    }

    // ---- TABLEAU 2 : (vii) pre-broye, cycles elastiques a axes tournants -------
    std::printf("\n=== TABLEAU 2 : (vii) element PRE-BROYE (uniax 3 %% -> omega_c), decharge a -100 MPa, puis cycles\n"
                "    deviatoriques tournants d'amplitude a, sens +/-. Cellule = W_2 / Psi_max (cycle 2), * = CREE ===\n");
    std::printf("%-30s", "a, sens, vitesse");
    for (const auto& v : V) std::printf(" %14s", v.tag.substr(0, 14).c_str());
    std::printf("\n");
    nCree2 = 0;
    for (double a : {0.0005, 0.001, 0.002})
        for (double rate : {1.0e3, 1.0e5})
            for (int dir : {+1, -1}) {
                char tag[64];
                std::snprintf(tag, sizeof tag, "a=%.2g%% sens %c @%.0e", a * 100, dir > 0 ? '+' : '-', rate);
                std::printf("%-30s", tag);
                for (size_t j = 0; j < V.size(); ++j) {
                    double om = 0.0;
                    CycleRes r = precrushed(*laws[j], a, dir, rate, 1.0e3, &om);
                    const double pm = std::max({r.psiMax[0], r.psiMax[1], r.psiMax[2], 1.0});
                    std::printf(" %13.3e%s", r.W[1] / pm, r.cree ? "*" : " ");
                    if (r.cree) ++nCree2;
                    csv << "(vii) " << tag << ',' << rate << ',' << V[j].tag;
                    for (int c = 0; c < 3; ++c) csv << ',' << r.W[c];
                    for (int c = 0; c < 3; ++c) csv << ',' << r.dPsi[c];
                    for (int c = 0; c < 3; ++c) csv << ',' << r.dDiss[c];
                    csv << ',' << pm << ',' << r.resid << ',' << r.epvEnd << ',' << r.DcEnd
                        << ',' << r.DEnd << ',' << (r.cree ? 1 : 0) << '\n';
                }
                std::printf("\n");
            }
    // detail complet pour FULL et -compDam a a = 0,1 %, 1e3/s
    std::printf("\n--- detail (vii) a = 0,1 %%, 1e3/s : W_1..3, dPsi_eff, dDiss, omega_c fige, epv ---\n");
    for (size_t j = 0; j < V.size(); ++j)
        for (int dir : {+1, -1}) {
            double om = 0.0;
            CycleRes r = precrushed(*laws[j], 0.001, dir, 1.0e3, 1.0e3, &om);
            std::printf("%-16s sens %c  W = [%+.4e %+.4e %+.4e]  dPsi_eff = [%+.2e %+.2e %+.2e]  dDiss = [%+.2e %+.2e %+.2e]  PsiMax %.3e  omega_c %.4f -> %.4f  D %.3f  epv %.4f -> %.4f\n",
                        V[j].tag.c_str(), dir > 0 ? '+' : '-', r.W[0], r.W[1], r.W[2],
                        r.dPsi[0], r.dPsi[1], r.dPsi[2], r.dDiss[0], r.dDiss[1], r.dDiss[2],
                        std::max({r.psiMax[0], r.psiMax[1], r.psiMax[2]}), om, r.DcEnd, r.DEnd,
                        r.epvEnd, r.epvEnd);
        }
    }   // !only3
    // ---- TABLEAU 3 : (viii) pre-broye, rotation + RESPIRATION volumique -------
    std::printf("\n=== TABLEAU 3 : (viii) element PRE-BROYE (uniax 5 %% a 1e2/s -> omega_c), decharge a -100 MPa, puis\n"
                "    eps = eps_t + rotDev(a) + b sin(2 pi u + phi0) I, cycles ELASTIQUES (omega_c fige).\n"
                "    Cellule = W_2 / Psi_max ; * = CREE. Prediction : W = omega_c lambda/(2G) oint tr(sig-) d(tr eps) ===\n");
    std::printf("%-30s", "a=b, phi0, vitesse");
    for (const auto& v : V) std::printf(" %14s", v.tag.substr(0, 14).c_str());
    std::printf("\n");
    int nCree3 = 0;
    std::vector<std::string> detail;
    for (double a : {0.0003, 0.0005})
        for (double rate : {1.0e3, 1.0e5})
            for (int ip = 0; ip < 4; ++ip) {
                const double phi0 = ip * M_PI / 2;
                char tag[64];
                std::snprintf(tag, sizeof tag, "a=b=%.2g%% phi0=%dpi/2 @%.0e", a * 100, ip, rate);
                std::printf("%-30s", tag);
                for (size_t j = 0; j < V.size(); ++j) {
                    BreathRes B = breathing(*laws[j], a, a, phi0, rate, 0.05, 1.0e2);
                    const CycleRes& r = B.r;
                    const double pm = std::max({r.psiMax[0], r.psiMax[1], r.psiMax[2], 1.0});
                    std::printf(" %13.3e%s", r.W[1] / pm, r.cree ? "*" : " ");
                    if (r.cree) ++nCree3;
                    csv << "(viii) " << tag << ',' << rate << ',' << V[j].tag;
                    for (int c = 0; c < 3; ++c) csv << ',' << r.W[c];
                    for (int c = 0; c < 3; ++c) csv << ',' << r.dPsi[c];
                    for (int c = 0; c < 3; ++c) csv << ',' << r.dDiss[c];
                    csv << ',' << pm << ',' << r.resid << ',' << r.epvEnd << ',' << r.DcEnd
                        << ',' << r.DEnd << ',' << (r.cree ? 1 : 0) << '\n';
                    if (a == 0.0003 && rate == 1.0e3) {
                        char line[512];
                        std::snprintf(line, sizeof line,
                            "%-16s phi0=%dpi/2  W = [%+.4e %+.4e %+.4e]  W_pred = [%+.4e %+.4e %+.4e]  dPsi_eff = [%+.1e %+.1e %+.1e]  dDiss = [%+.1e %+.1e %+.1e]  PsiMax %.3e  omega_c %.3f -> %.3f  D %.3f  epv %.4f  |e0-f0| %.1e |e1-f0| %.1e",
                            V[j].tag.c_str(), ip, r.W[0], r.W[1], r.W[2], B.Wpred[0], B.Wpred[1], B.Wpred[2],
                            r.dPsi[0], r.dPsi[1], r.dPsi[2], r.dDiss[0], r.dDiss[1], r.dDiss[2],
                            std::max({r.psiMax[0], r.psiMax[1], r.psiMax[2]}), B.omega, r.DcEnd, r.DEnd, r.epvEnd,
                            B.startErr, B.endErr);
                        detail.push_back(line);
                    }
                }
                std::printf("\n");
                std::fflush(stdout);
            }
    std::printf("\n--- detail (viii) a = b = 0,03 %%, 1e3/s : mesure contre prediction analytique ---\n");
    for (const auto& l : detail) std::printf("%s\n", l.c_str());

    std::printf("\nCSV : %s/bench_energie.csv   (CREE tableau 1 : %d, tableau 2 : %d, tableau 3 : %d)\n",
                dir.c_str(), nCree, nCree2, nCree3);
    return (nCree + nCree2 + nCree3) > 0 ? 1 : 0;
}
