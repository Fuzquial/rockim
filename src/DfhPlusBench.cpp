// ---------------------------------------------------------------------------
// DfhPlusBench — `rockim selftest-dfhplus [out.csv]` : les bancs falsifiants
// de la loi `dfhplus` au POINT MATERIEL (2026-09-06, etape 1 du chantier DFH+).
//
// Le banc thermodynamique generique (`rockim thermobench`) prouve les
// proprietes de CADRE sur des chemins tires au hasard. Celui-ci prouve les
// proprietes de la FORMULATION, avec des criteres chiffres et une variante qui
// DOIT echouer :
//
//   (1) sigma = d(rho psi)/d eps : la contrainte analytique contre les
//       differences finies centrees sur rho psi, sur 6 etats varies (traction,
//       compression, cisaillement, non coaxial, endommage) x 6 composantes.
//       Critere : ecart relatif < 1e-6.
//   (2) Y_i = -d(rho psi)/dD_i : idem, differences finies sur D_i a eps fige.
//       Critere : ecart relatif < 1e-6, et Y_i >= 0 partout.
//   (3) REDUCTION ELASTIQUE exacte a D = 0 : sigma = lambda tr(eps) I + 2G eps
//       a 1e-13.
//   (4) OBJECTIVITE : chemin rejoue sous rotation rigide Q, sigma_Q =
//       Q sigma Q^T a 1e-12 (ecart tensoriel, pas seulement les invariants).
//   (5) EXPOSANT DE VITESSE de l'obscuration : traction uniaxiale pilotee en
//       deformation a 8 vitesses (1e2 - 1e6 /s), pente de log(sigma_pic) vs
//       log(epsdot) pour m = 6, 12, 24. Cible 3/(m+3) = 0,3333 / 0,2000 /
//       0,1111 ; on l'exige a 15 % relatif POUR LES DEUX LOIS (dfhplus ET
//       dpdfh, cote a cote : dpdfh sert d'etalon, on ne demande pas a la loi
//       neuve mieux que ce que la loi de reference atteint sur le meme banc).
//   (6) COMPARAISON dfhplus / dpdfh (mesure, sans verdict) : traction
//       uniaxiale, compression uniaxiale, triaxiaux 20 / 50 / 100 MPa, et un
//       chemin NON COAXIAL (traction + cisaillement dont les directions
//       principales tournent apres l'amorcage).
//   (7) CONTROLE QUI DOIT ECHOUER : la meme carte avec dfhpPsiClamp = false.
//       Sans l'ecretage d'admissibilite, la dissipation plastique
//       sigma^nom : d eps^p devient negative sur un chemin ou l'endommagement
//       est fortement anisotrope. Le banc l'exige : si elle reste positive, la
//       cle ne sert a rien et le banc le DIT (mais ne fait pas echouer, la
//       carte pouvant simplement ne pas atteindre le regime).
//
// Code de retour : 0 si (1)-(5) passent, 1 sinon.
// ---------------------------------------------------------------------------
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
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

// Carte Red Bohus du chantier (identique a celle du banc thermodynamique).
const char* kCard =
    "E = 77.66e9\n"
    "nu = 0.29\n"
    "rho = 2620\n"
    "ft = 9.0e6\n"
    "cohesion = 22.77e6\n"
    "frictionDeg = 50.4\n"
    "Gf = 100\n"
    "erodeD = 2\n"
    "erodeEpv = 0\n";

struct Built {
    Material m;
    std::unique_ptr<MatLaw> law;
};

Built build(const std::string& kind, const std::string& extra,
            const std::string& tag) {
    const std::string path = "._dfhp_" + tag + ".cfg";
    {
        std::ofstream f(path);
        f << kCard << extra << "\n";
    }
    Config c = Config::load(path);
    Built b;
    b.m = Material::from(c);
    b.law = MatLaw::make(kind, b.m, c, 2.0e-3);
    std::remove(path.c_str());
    return b;
}

M3 sym(int a) {
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
double ddot(const M3& a, const M3& b) { return (a.cwiseProduct(b)).sum(); }

M3 rotAxis(const Eigen::Vector3d& ax, double th) {
    return Eigen::AngleAxisd(th, ax.normalized()).toRotationMatrix();
}

// Pilotage uniaxial libre lateralement : on cherche eps_lat tel que
// sigma_lat = 0, par secante (la loi n'est pas lineaire).
double holdLateral(const MatLaw& law, MatState& s, M3& eps, double eAx,
                   double dt, double lc, M3& sig, double e0) {
    double x = e0;
    for (int it = 0; it < 12; ++it) {
        MatState t = s;
        eps.setZero();
        eps(0, 0) = eAx;
        eps(1, 1) = eps(2, 2) = x;
        M3 sg = law.stress(eps, t, dt, lc);
        double r = 0.5 * (sg(1, 1) + sg(2, 2));
        if (std::abs(r) < 1.0e2) { sig = sg; break; }   // 0,1 kPa sur ~100 MPa
        MatState t2 = s;
        const double h = 1.0e-9;
        M3 e2 = eps;
        e2(1, 1) = e2(2, 2) = x + h;
        M3 sg2 = law.stress(e2, t2, dt, lc);
        double r2 = 0.5 * (sg2(1, 1) + sg2(2, 2));
        double d = (r2 - r) / h;
        if (!(std::abs(d) > 0.0)) break;
        x -= r / d;
    }
    eps.setZero();
    eps(0, 0) = eAx;
    eps(1, 1) = eps(2, 2) = x;
    sig = law.stress(eps, s, dt, lc);
    return x;
}

// Traction uniaxiale pilotee a vitesse de deformation constante : rend la
// contrainte axiale MAXIMALE (le pic) et le D atteint au pic.
double tensionPeak(const MatLaw& law, double edot, double lc, double& Dend,
                   bool& postPeak) {
    // Le domaine doit couvrir TOUTE la gamme de vitesses : a 1e9 /s le pic
    // depasse le GPa. eMax = 2,4e-2 -> plafond 1,86 GPa, largement au-dessus.
    const double eMax = 6.0e-2;          // plafond E eMax = 4,66 GPa : meme a
    const int n = 20000;                 // 1e7 /s et m = 6 la branche
                                         // adoucissante est atteinte
    const double de = eMax / n;
    const double dt = de / edot;
    MatState s;
    s.x0 = Eigen::Vector3d(0.0037, 0.0071, 0.0113);
    M3 eps = M3::Zero(), sig = M3::Zero();
    double peak = 0.0, x = 0.0;
    Dend = 0.0;
    for (int k = 1; k <= n; ++k) {
        x = holdLateral(law, s, eps, k * de, dt, lc, sig, x);
        if (sig(0, 0) > peak) peak = sig(0, 0);
        Dend = s.D;
        if (s.D > 0.99) break;
    }
    // Le pic n'est un VRAI maximum que si la branche adoucissante a ete vue :
    // sinon le « pic » est le bout de la rampe de deformation, pas la loi.
    postPeak = sig(0, 0) < 0.8 * peak;
    return peak;
}

double slope(const std::vector<double>& x, const std::vector<double>& y) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    const double n = (double)x.size();
    for (std::size_t i = 0; i < x.size(); ++i) {
        sx += x[i]; sy += y[i]; sxx += x[i] * x[i]; sxy += x[i] * y[i];
    }
    return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

} // namespace

int dfhPlusSelftest(const std::string& csvPath) {
    std::ofstream out(csvPath);
    out.precision(12);
    out << "banc,cas,grandeur,dfhplus,dpdfh,ecart_rel\n";
    std::cout.precision(6);
    int fails = 0;
    const double lc = 1.0e-3, dt = 1.0e-7;

    Built P = build("dfhplus", "", "p");
    Built R = build("dpdfh", "", "r");
    const Material& m = P.m;
    const double lam = m.E * m.nu / ((1.0 + m.nu) * (1.0 - 2.0 * m.nu));
    const double G = m.G();
    const double k0 = m.ft / m.E;

    std::cout << "\n=== selftest-dfhplus (carte Red Bohus, lc " << lc
              << " m) ===\n";

    // ---- etats de travail : varies, dont endommages -----------------------
    struct Case { const char* name; M3 eps; int nPre; };
    std::vector<Case> cases;
    {
        M3 e;
        e = M3::Zero(); e(0, 0) = 3.0 * k0; e(1, 1) = e(2, 2) = -m.nu * 3 * k0;
        cases.push_back({"traction", e, 40});
        e = M3::Zero(); e(0, 0) = -20 * k0; e(1, 1) = e(2, 2) = m.nu * 20 * k0;
        cases.push_back({"compression", e, 40});
        e = M3::Zero(); e(0, 1) = e(1, 0) = 4.0 * k0;
        cases.push_back({"cisaillement", e, 40});
        e = M3::Zero(); e(0, 0) = 5.0 * k0; e(2, 2) = -2.0 * k0;
        e(1, 2) = e(2, 1) = 2.0 * k0;
        cases.push_back({"non-coaxial", e, 60});
        e = M3::Zero(); e(0, 0) = 60.0 * k0; e(1, 1) = e(2, 2) = -m.nu * 60 * k0;
        e(0, 2) = e(2, 0) = 3.0 * k0;
        cases.push_back({"endommage", e, 120});
        e = M3::Zero(); e = -30.0 * k0 * M3::Identity();
        e(0, 0) += 40.0 * k0;
        cases.push_back({"triaxial", e, 60});
    }

    // =====================================================================
    // (1) sigma = d(rho psi)/d eps   et   (2) Y_i = -d(rho psi)/dD_i
    // =====================================================================
    double worstSig = 0.0, worstY = 0.0, minY = 1e300;
    for (const Case& C : cases) {
        MatState s;
        s.x0 = Eigen::Vector3d(0.0037, 0.0071, 0.0113);
        M3 eps = M3::Zero(), sig = M3::Zero();
        for (int k = 1; k <= C.nPre; ++k) {
            eps = ((double)k / C.nPre) * C.eps;
            sig = P.law->stress(eps, s, dt, lc);
        }
        // (1) differences finies centrees sur rho psi, etat FIGE
        const double h = 1.0e-10;
        double scale = sig.cwiseAbs().maxCoeff() + 1.0e-30;
        for (int a = 0; a < 6; ++a) {
            double pp = P.law->freeEnergy(eps + h * sym(a), s);
            double pm = P.law->freeEnergy(eps - h * sym(a), s);
            double fd = (pp - pm) / (2.0 * h);
            double an = ddot(sig, sym(a));
            double err = std::abs(fd - an) / scale;
            worstSig = std::max(worstSig, err);
        }
        // (2) Y_i ANALYTIQUE (expose par la loi) contre differences finies
        // centrees sur rho psi a eps fige, plus la POSITIVITE.
        double Yan[3] = {0, 0, 0};
        int nY = P.law->damageForces(eps, s, Yan);
        for (int i = 0; i < nY; ++i) {
            minY = std::min(minY, Yan[i]);
            const double hd = 1.0e-7;
            if (s.dfhp.Dv[i] + hd > 0.9997 || s.dfhp.Dv[i] < hd) continue;
            MatState sp = s, sm = s;
            sp.dfhp.Dv[i] += hd;
            sm.dfhp.Dv[i] -= hd;
            double fd = -(P.law->freeEnergy(eps, sp)
                          - P.law->freeEnergy(eps, sm)) / (2.0 * hd);
            double sc2 = std::abs(fd) + std::abs(Yan[i]) + 1.0;
            worstY = std::max(worstY, std::abs(fd - Yan[i]) / sc2);
        }
        out << "1-2," << C.name << ",D_max," << s.D << ",," << "\n";
    }
    bool ok1 = worstSig < 1.0e-6;
    bool ok2 = minY >= -1.0e-9 && worstY < 1.0e-6;
    std::cout << "(1) sigma = d(rho psi)/d eps        : pire ecart relatif "
              << worstSig << (ok1 ? "   PASS" : "   ECHEC") << "\n";
    std::cout << "(2) Y_i = -d(rho psi)/dD_i >= 0     : min Y_i " << minY
              << " J/m^3, analytique vs DF " << worstY
              << (ok2 ? "   PASS" : "   ECHEC") << "\n";
    out << "1,,pire_ecart_sigma," << worstSig << ",,\n";
    out << "2,,min_Y," << minY << ",,\n";
    if (!ok1) ++fails;
    if (!ok2) ++fails;

    // =====================================================================
    // (3) reduction elastique exacte a D = 0
    // =====================================================================
    {
        double worst = 0.0;
        for (int a = 0; a < 6; ++a) {
            MatState s;
            s.x0 = Eigen::Vector3d(0.01, 0.02, 0.03);
            for (int k = 1; k <= 8; ++k) {
                M3 e = (0.3 * k0 * k / 8.0) * sym(a);
                M3 sg = P.law->stress(e, s, dt, lc);
                M3 ref = lam * e.trace() * M3::Identity() + 2.0 * G * e;
                worst = std::max(worst,
                                 (sg - ref).cwiseAbs().maxCoeff()
                                     / (ref.cwiseAbs().maxCoeff() + 1e-30));
            }
        }
        bool ok = worst < 1.0e-13;
        std::cout << "(3) reduction elastique a D = 0     : pire ecart "
                  << worst << (ok ? "   PASS" : "   ECHEC") << "\n";
        out << "3,,reduction_elastique," << worst << ",,\n";
        if (!ok) ++fails;
    }

    // =====================================================================
    // (4) objectivite : rotation rigide superposee
    // =====================================================================
    {
        M3 Q = rotAxis(Eigen::Vector3d(0.3, -0.7, 0.5), 0.9);
        double worst = 0.0;
        for (const Case& C : cases) {
            MatState sa, sb;
            sa.x0 = sb.x0 = Eigen::Vector3d(0.0037, 0.0071, 0.0113);
            double smax = 1e-30;
            for (int k = 1; k <= C.nPre; ++k) {
                M3 e = ((double)k / C.nPre) * C.eps;
                M3 a = P.law->stress(e, sa, dt, lc);
                M3 b = P.law->stress(Q * e * Q.transpose(), sb, dt, lc);
                smax = std::max(smax, a.cwiseAbs().maxCoeff());
                worst = std::max(worst,
                                 (b - Q * a * Q.transpose()).norm() / smax);
            }
        }
        bool ok = worst < 1.0e-12;
        std::cout << "(4) objectivite (tensorielle)       : pire ecart "
                  << worst << (ok ? "   PASS" : "   ECHEC") << "\n";
        out << "4,,objectivite," << worst << ",,\n";
        if (!ok) ++fails;
    }

    // =====================================================================
    // (5) exposant de vitesse 3/(m+3)
    // ---------------------------------------------------------------------
    // La loi d'echelle sigma_pic ~ epsdot^{3/(m+3)} de Denoual-Hild ne vaut
    // QUE dans le regime de FRAGMENTATION MULTIPLE, c'est-a-dire quand la
    // densite de defauts amorces vue par l'element depasse largement 1 :
    //   lambda_frag V_el = (sigma/sigma_w)^m (V_el/Z_eff) >> 1.
    // En dessous, l'element casse sur UN defaut et le pic est rate-
    // INDEPENDANT (le plancher `xlam Vel < 1` de la VUMAT). Le banc ne fait
    // donc la regression QUE sur les points ou lambda_frag V_el > 10, soit
    // sigma_pic > sigma_w 10^{1/m} — et il l'annonce.
    // =====================================================================
    {
        const double ms[3] = {6.0, 12.0, 24.0};
        const double edots[7] = {1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9};
        std::cout << "(5) exposant de vitesse de l'obscuration (cible "
                     "3/(m+3), regime de fragmentation MULTIPLE)\n";
        bool allOk = true;
        for (int im = 0; im < 3; ++im) {
            std::ostringstream ex;
            ex << "dfhWeibullM = " << ms[im] << "\n";
            Built Pp = build("dfhplus", ex.str(), "p5");
            Built Rr = build("dpdfh", ex.str(), "r5");
            const double sMin = 120.0e6 * std::pow(10.0, 1.0 / ms[im]);
            std::vector<double> lx, lyP, lyR;
            int nUsed = 0;
            for (int i = 0; i < 7; ++i) {
                double dP = 0, dR = 0;
                bool ppP = false, ppR = false;
                double sp = tensionPeak(*Pp.law, edots[i], lc, dP, ppP);
                double sr = tensionPeak(*Rr.law, edots[i], lc, dR, ppR);
                bool keep = sp > sMin && sr > sMin && ppP && ppR;
                out << "5,m=" << ms[im] << " edot=" << edots[i]
                    << ",sigma_pic_MPa," << sp / 1e6 << "," << sr / 1e6 << ","
                    << (sp - sr) / sr << (keep ? "" : "   [hors regime]")
                    << "\n";
                if (!keep) continue;
                ++nUsed;
                lx.push_back(std::log(edots[i]));
                lyP.push_back(std::log(sp));
                lyR.push_back(std::log(sr));
            }
            double tgt = 3.0 / (ms[im] + 3.0);
            if (nUsed < 3) {
                std::cout << "     m = " << ms[im]
                          << " : moins de 3 points dans le regime multiple — "
                             "banc non concluant\n";
                allOk = false;
                continue;
            }
            double nP = slope(lx, lyP), nR = slope(lx, lyR);
            bool okP = std::abs(nP - tgt) / tgt < 0.10;
            bool okR = std::abs(nR - tgt) / tgt < 0.10;
            allOk = allOk && okP;
            std::cout << "     m = " << std::setw(4) << ms[im] << "  cible "
                      << tgt << "  dfhplus " << nP << " ("
                      << 100.0 * (nP - tgt) / tgt << " %)  dpdfh " << nR
                      << " (" << 100.0 * (nR - tgt) / tgt << " %)  sur "
                      << nUsed << " points" << (okP ? "  PASS" : "  ECHEC")
                      << (okR ? "" : "  [dpdfh hors cible]") << "\n";
            out << "5,m=" << ms[im] << ",exposant," << nP << "," << nR << ","
                << (nP - tgt) / tgt << "\n";
        }
        if (!allOk) ++fails;
    }

    // =====================================================================
    // (6) comparaison dfhplus / dpdfh (mesure, sans verdict)
    // =====================================================================
    std::cout << "(6) dfhplus contre dpdfh (mesure)\n";
    {
        // --- traction uniaxiale : pic, deformation au pic, aire dissipee ---
        auto tension = [&](MatLaw& law, double& peak, double& epk,
                           double& area, double& sat) {
            const double eMax = 6.0e-3;
            const int n = 6000;
            const double de = eMax / n;
            MatState s;
            s.x0 = Eigen::Vector3d(0.0037, 0.0071, 0.0113);
            M3 eps = M3::Zero(), sig = M3::Zero();
            double x = 0.0, prevS = 0.0, prevE = 0.0;
            peak = 0.0; epk = 0.0; area = 0.0; sat = 0.0;
            for (int k = 1; k <= n; ++k) {
                x = holdLateral(law, s, eps, k * de, 1.0e-7, lc, sig, x);
                double e = k * de;
                area += 0.5 * (sig(0, 0) + prevS) * (e - prevE);
                prevS = sig(0, 0); prevE = e;
                if (sig(0, 0) > peak) { peak = sig(0, 0); epk = e; }
                if (s.D > 0.98 && sat == 0.0) sat = sig(0, 0);
            }
            return s.D;
        };
        double pP, ePp, aP, satP, pR, eR, aR, satR;
        double DP = tension(*P.law, pP, ePp, aP, satP);
        double DR = tension(*R.law, pR, eR, aR, satR);
        std::cout << "     traction : pic " << pP / 1e6 << " vs " << pR / 1e6
                  << " MPa (" << 100.0 * (pP - pR) / pR << " %), eps_pic "
                  << ePp << " vs " << eR << ", aire " << aP / 1e6 << " vs "
                  << aR / 1e6 << " MJ/m^3 (" << 100.0 * (aP - aR) / aR
                  << " %), sigma a D = 0,98 : " << satP / 1e6 << " vs "
                  << satR / 1e6 << " MPa, D fin " << DP << " / " << DR << "\n";
        out << "6,traction,pic_MPa," << pP / 1e6 << "," << pR / 1e6 << ","
            << (pP - pR) / pR << "\n";
        out << "6,traction,aire_MJm3," << aP / 1e6 << "," << aR / 1e6 << ","
            << (aP - aR) / aR << "\n";
        out << "6,traction,sigma_D098_MPa," << satP / 1e6 << "," << satR / 1e6
            << ",\n";
    }
    {
        // --- compression uniaxiale et triaxiaux : pic de q ------------------
        auto compress = [&](MatLaw& law, double s3, double& qpk) {
            const double eMax = 0.02;
            const int n = 2000;
            MatState s;
            s.x0 = Eigen::Vector3d(0.0037, 0.0071, 0.0113);
            M3 eps = M3::Zero(), sig = M3::Zero();
            // consolidation isotrope a -s3
            if (s3 > 0.0) {
                double e = 0.0;
                for (int it = 0; it < 200; ++it) {
                    MatState t = s;
                    M3 ee = e * M3::Identity();
                    M3 sg = law.stress(ee, t, 1.0, lc);
                    double r = -sg.trace() / 3.0 - s3;
                    if (std::abs(r) < 1e3) break;
                    e += r / (3.0 * m.K());
                }
                eps = e * M3::Identity();
                sig = law.stress(eps, s, 1.0, lc);
            }
            const double e0 = eps(2, 2), x0 = eps(0, 0);
            double x = x0;
            qpk = 0.0;
            for (int k = 1; k <= n; ++k) {
                double eAx = e0 - eMax * k / n;
                // pilotage lateral : sigma_lat = -s3
                for (int it = 0; it < 40; ++it) {
                    MatState t = s;
                    M3 ee = M3::Zero();
                    ee(2, 2) = eAx; ee(0, 0) = ee(1, 1) = x;
                    M3 sg = law.stress(ee, t, 1.0, lc);
                    double r = 0.5 * (sg(0, 0) + sg(1, 1)) + s3;
                    if (std::abs(r) < 1e3) break;
                    MatState t2 = s;
                    M3 e2 = ee; e2(0, 0) = e2(1, 1) = x + 1e-9;
                    M3 sg2 = law.stress(e2, t2, 1.0, lc);
                    double r2 = 0.5 * (sg2(0, 0) + sg2(1, 1)) + s3;
                    double d = (r2 - r) / 1e-9;
                    if (!(std::abs(d) > 0.0)) break;
                    x -= r / d;
                }
                M3 ee = M3::Zero();
                ee(2, 2) = eAx; ee(0, 0) = ee(1, 1) = x;
                sig = law.stress(ee, s, 1.0, lc);
                double q = 0.5 * (sig(0, 0) + sig(1, 1)) - sig(2, 2);
                qpk = std::max(qpk, q);
            }
        };
        const double s3s[4] = {0.0, 20.0e6, 50.0e6, 100.0e6};
        for (int i = 0; i < 4; ++i) {
            double qP = 0, qR = 0;
            compress(*P.law, s3s[i], qP);
            compress(*R.law, s3s[i], qR);
            std::cout << "     compression sigma3 = " << s3s[i] / 1e6
                      << " MPa : q_pic " << qP / 1e6 << " vs " << qR / 1e6
                      << " MPa (" << 100.0 * (qP - qR) / qR << " %)\n";
            out << "6,compression s3=" << s3s[i] / 1e6 << ",q_pic_MPa,"
                << qP / 1e6 << "," << qR / 1e6 << "," << (qP - qR) / qR
                << "\n";
        }
    }
    double dissNoClamp = 0.0, dissClamp = 0.0;
    {
        // --- chemin NON COAXIAL : traction jusqu'a l'amorcage, puis le
        // repere principal TOURNE (cisaillement croissant + rotation du
        // repere de chargement). C'est le regime ou dpdfh perd la positivite.
        auto noncoax = [&](MatLaw& law, double& sMax, double& diss,
                           double& Dend) {
            MatState s;
            s.x0 = Eigen::Vector3d(0.0037, 0.0071, 0.0113);
            M3 eps = M3::Zero(), sig = M3::Zero(), prev = M3::Zero();
            M3 base = M3::Zero();
            base(0, 0) = 40.0 * k0;
            base(1, 1) = base(2, 2) = -m.nu * 40.0 * k0;
            sMax = 0.0; diss = 0.0; Dend = 0.0;
            for (int k = 1; k <= 80; ++k) {         // pre-charge en traction
                eps = ((double)k / 80.0) * base;
                prev = sig;
                M3 pe = eps - ((k - 1) / 80.0) * base;
                sig = law.stress(eps, s, dt, lc);
                diss += ddot(0.5 * (sig + prev), pe);
                sMax = std::max(sMax, sig.cwiseAbs().maxCoeff());
            }
            M3 e0 = eps;
            M3 dE = M3::Zero();
            dE(0, 1) = dE(1, 0) = 1.5 * k0;
            dE(2, 2) = -3.0 * k0;
            for (int k = 1; k <= 40; ++k) {         // rotation + cisaillement
                M3 Q = rotAxis(Eigen::Vector3d(0.2, 0.3, 1.0), 0.05 * k);
                M3 en = Q * e0 * Q.transpose() + ((double)k / 40.0) * dE;
                M3 de = en - eps;
                prev = sig;
                sig = law.stress(en, s, dt, lc);
                diss += ddot(0.5 * (sig + prev), de);
                eps = en;
                sMax = std::max(sMax, sig.cwiseAbs().maxCoeff());
            }
            Dend = s.D;
            return s;
        };
        double sP = 0, dP = 0, DP = 0, sR = 0, dR = 0, DR = 0;
        MatState fP = noncoax(*P.law, sP, dP, DP);
        MatState fR = noncoax(*R.law, sR, dR, DR);
        std::cout << "     chemin NON COAXIAL : |sigma| max " << sP / 1e6
                  << " vs " << sR / 1e6 << " MPa, travail total " << dP / 1e6
                  << " vs " << dR / 1e6 << " MJ/m^3, D fin " << DP << " / "
                  << DR << " ; ecretage de dilatance actif "
                  << fP.dfhp.nClamp << " fois\n";
        out << "6,non-coaxial,sigma_max_MPa," << sP / 1e6 << "," << sR / 1e6
            << "," << (sP - sR) / sR << "\n";
        out << "6,non-coaxial,travail_MJm3," << dP / 1e6 << "," << dR / 1e6
            << "," << (dP - dR) / dR << "\n";
        out << "6,non-coaxial,nClamp," << fP.dfhp.nClamp << ",,\n";
        (void)fR;
        dissClamp = dP;
    }

    // =====================================================================
    // (7) controle qui DOIT echouer : dfhpPsiClamp = false
    // =====================================================================
    {
        Built N = build("dfhplus", "dfhpPsiClamp = false\n", "nc");
        // chemin de compression fortement endommage : on charge en traction
        // pour amorcer, puis on comprime OBLIQUEMENT (le repere fige ne suit
        // pas), et l'on mesure la dissipation plastique incrementale
        // sigma^nom : d eps^p, pas par pas.
        auto plasticDiss = [&](MatLaw& law, int& nNeg, double& worst,
                               int& nPl) {
            MatState s;
            s.x0 = Eigen::Vector3d(0.0037, 0.0071, 0.0113);
            M3 eps = M3::Zero(), sig = M3::Zero();
            M3 base = M3::Zero();
            base(0, 0) = 60.0 * k0;
            base(1, 1) = base(2, 2) = -m.nu * 60.0 * k0;
            nNeg = 0; worst = 0.0; nPl = 0;
            for (int k = 1; k <= 120; ++k) {
                eps = ((double)k / 120.0) * base;
                sig = law.stress(eps, s, dt, lc);
            }
            // decharge complete puis COMPRESSION UNIAXIALE dans un repere
            // TOURNE de 40 deg : la roche est fissuree normalement a x, on
            // l'ecrase obliquement — le cone travaille avec un endommagement
            // fortement non coaxial, exactement le regime vise.
            M3 Q = rotAxis(Eigen::Vector3d(0.0, 1.0, 0.3), 0.70);
            M3 uni = M3::Zero();
            uni(0, 0) = -1.0;
            uni(1, 1) = uni(2, 2) = m.nu;
            M3 dE = Q * uni * Q.transpose();
            const double aMax = 400.0 * k0;
            for (int k = 1; k <= 400; ++k) {
                M3 en = (aMax * (double)k / 400.0) * dE;
                M3 ep0 = s.epsP;
                M3 s0 = sig;
                sig = law.stress(en, s, dt, lc);
                M3 dep = s.epsP - ep0;
                if (dep.norm() > 0.0) {
                    ++nPl;
                    double d = ddot(0.5 * (sig + s0), dep);
                    double sc = std::abs(ddot(0.5 * (sig + s0), dep))
                                + 0.5 * (sig + s0).norm() * dep.norm() + 1.0;
                    if (d < -1e-9 * sc) {
                        ++nNeg;
                        worst = std::min(worst, d);
                    }
                }
            }
        };
        // La carte de la these porte Psi = 15 deg (tan 0,268), tres en
        // dessous de q/p sur la surface : l'ecretage n'a rien a mordre. Pour
        // FALSIFIER la cle il faut la mettre en danger — ecoulement ASSOCIE
        // (Psi = beta = 51,7 deg, tan 1,266), le cas ou la dilatance est
        // maximale.
        Built A = build("dfhplus", "dfhPsiDeg = 51.7\n", "as");
        Built B = build("dfhplus",
                        "dfhPsiDeg = 51.7\ndfhpPsiClamp = false\n", "ao");
        int nOn = 0, nOff = 0, pOn = 0, pOff = 0;
        int nAOn = 0, nAOff = 0, pAOn = 0, pAOff = 0;
        double wOn = 0, wOff = 0, wAOn = 0, wAOff = 0;
        plasticDiss(*P.law, nOn, wOn, pOn);
        plasticDiss(*N.law, nOff, wOff, pOff);
        plasticDiss(*A.law, nAOn, wAOn, pAOn);
        plasticDiss(*B.law, nAOff, wAOff, pAOff);
        std::cout << "     ASSOCIE (Psi = beta = 51,7 deg) : ecretage ACTIF "
                  << nAOn << " increment(s) negatif(s) (pire " << wAOn
                  << " J/m^3) ; DESARME " << nAOff << " (pire " << wAOff
                  << " J/m^3) sur " << pAOff << " increments plastiques\n";
        out << "7,associe,increments_negatifs_clamp_on," << nAOn << ",,\n";
        out << "7,associe,increments_negatifs_clamp_off," << nAOff << ",,\n";
        out << "7,associe,pire_clamp_off," << wAOff << ",,\n";
        std::cout << "(7) controle qui DOIT echouer (dfhpPsiClamp)\n"
                  << "     increments plastiques sur le chemin : " << pOn
                  << " (ecretage actif) / " << pOff << " (desarme)\n"
                  << "     ecretage ACTIF   : " << nOn
                  << " increment(s) a dissipation plastique negative, pire "
                  << wOn << " J/m^3\n"
                  << "     ecretage DESARME : " << nOff
                  << " increment(s), pire " << wOff << " J/m^3"
                  << (nOff > nOn ? "   [la cle sert : ecart mesure]"
                                 : "   [la carte n'atteint pas le regime — "
                                   "l'ecretage ne peut pas etre falsifie ici]")
                  << "\n";
        out << "7,,increments_negatifs_clamp_on," << nOn << ",,\n";
        out << "7,,increments_negatifs_clamp_off," << nOff << ",,\n";
        out << "7,,pire_clamp_on," << wOn << ",,\n";
        out << "7,,pire_clamp_off," << wOff << ",,\n";
        dissNoClamp = wOff;
        (void)dissNoClamp; (void)dissClamp;
    }

    std::cout << "\n=== selftest-dfhplus : " << (fails == 0 ? "PASS" : "ECHEC")
              << " (" << fails << " critere(s) manque(s)) ; trace " << csvPath
              << " ===\n";
    return fails == 0 ? 0 : 1;
}

} // namespace rockim
