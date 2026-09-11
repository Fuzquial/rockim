// Banc COURT du curseur de Coulomb de la branche jointTSL = camacho (3D),
// correctif du 2026-09-11 (RETOUR_v3 §1.4 quater ; note v3 §3.4).
// Noyau bancheable sans le solveur : rockim::camachoFrictionSlider
// (include/rockim/Fdem3dSolver.hpp) + jtsl::shearCap a cohesion nulle.
//
// Criteres FALSIFIANTS :
//  (i)   trajet en traction pure (t_n > 0) : tau_fric == 0 sur tout le trajet,
//        et slip suit delta_s (aucune raideur de collage cachee) ;
//  (ii)  cycle de glissement alterne sous compression : dissipation par cycle
//        = 2 f_cap (P - 2 f_cap/pj), P = amplitude crete a crete, >= 0 ;
//        travail cumule jamais < -f_cap^2/(2 pj) ;
//  (iii) continuite en delta_s = 0 : |tau(+eps) - tau(-eps)| = 2 pj eps -> 0.
// Variante qui DOIT echouer : la forme retiree, tau = f_cap ds/|ds|
// (saut de 2 f_cap en (iii)) ; et un ressort cape SANS retour de glissement
// (dissipation nulle en (ii)).
//
// Compile (depuis rockim_g1, vcvars64 charge) :
//   cl /nologo /std:c++17 /EHsc /D_USE_MATH_DEFINES /DNOMINMAX /I include
//      /I ..\rockim\eigen-3.4.0 /Fo:<obj>\ /Fe:<obj>\check_camacho_friction3d.exe
//      tests_f2\check_camacho_friction3d.cpp
#include <cmath>
#include <cstdio>
#include <Eigen/Dense>
#include "rockim/Fdem3dSolver.hpp"
using namespace rockim;

int main() {
    int fail = 0;
    auto chk = [&](bool ok, const char* what, double a, double b) {
        std::printf("  %-58s %s   (%.9g vs %.9g)\n", what, ok ? "PASS" : "FAIL",
                    a, b);
        if (!ok) ++fail;
    };
    const double pj = 3.0e14;      // Pa/m, ~10 E/h
    const double mu = 0.7;
    const double tnC = -10.0e6;    // compression 10 MPa
    const double fCapC = jtsl::shearCap(0.0, mu, tnC, 1.0, true);
    chk(std::fabs(fCapC - mu * 10.0e6) < 1e-3, "f_cap = mu <-t_n> (D = 1)",
        fCapC, mu * 10.0e6);
    chk(jtsl::shearCap(0.0, mu, 5.0e6, 1.0, true) == 0.0,
        "f_cap = 0 en traction", jtsl::shearCap(0.0, mu, 5.0e6, 1.0, true), 0.0);
    chk(jtsl::shearCap(0.0, mu, tnC, 0.0, true) == 0.0,
        "f_cap = 0 a D = 0 sous fricMob", jtsl::shearCap(0.0, mu, tnC, 0.0, true),
        0.0);

    // (i) traction pure : trajet de glissement quelconque, t_n = +5 MPa
    {
        Eigen::Vector3d slip = Eigen::Vector3d::Zero();
        double maxTau = 0.0, maxGap = 0.0;
        const double fCap = jtsl::shearCap(0.0, mu, 5.0e6, 1.0, true);
        for (int i = 0; i <= 2000; ++i) {
            double th = 2.0 * M_PI * i / 500.0;
            Eigen::Vector3d ds(2e-6 * std::sin(th), 1e-6 * std::sin(2 * th), 0);
            Eigen::Vector3d tf = camachoFrictionSlider(ds, pj, fCap, slip);
            maxTau = std::max(maxTau, tf.norm());
            maxGap = std::max(maxGap, (slip - ds).norm());
        }
        chk(maxTau == 0.0, "(i) traction pure : max |tau_fric| == 0", maxTau, 0.0);
        chk(maxGap == 0.0, "(i) traction pure : slip suit delta_s", maxGap, 0.0);
    }

    // (ii) cycle alterne sous compression, amplitude crete A, P = 2A
    auto cycleWork = [&](bool withReturn, double A, int nCyc, double& minCum) {
        Eigen::Vector3d slip = Eigen::Vector3d::Zero();
        Eigen::Vector3d dsPrev = Eigen::Vector3d::Zero();
        Eigen::Vector3d tfPrev = Eigen::Vector3d::Zero();
        double W = 0.0, Wlast = 0.0;
        minCum = 0.0;
        const int N = 4000;
        // rodage : un cycle complet avant la mesure (etat stationnaire)
        for (int c = 0; c < nCyc + 1; ++c) {
            Wlast = W;
            for (int i = 1; i <= N; ++i) {
                double th = 2.0 * M_PI * i / N;
                Eigen::Vector3d ds(A * std::sin(th), 0, 0);
                Eigen::Vector3d tf;
                if (withReturn) {
                    tf = camachoFrictionSlider(ds, pj, fCapC, slip);
                } else {                      // ressort cape SANS retour
                    tf = pj * ds;
                    double tt = tf.norm();
                    if (tt > fCapC) tf *= fCapC / tt;
                }
                // travail de la traction de joint contre le glissement :
                // dissipation quand > 0 (tf s oppose au mouvement relatif).
                // Trapeze : le controle « ressort cape sans retour » derive
                // d un potentiel, son travail sur un cycle ferme doit sortir
                // nul a l erreur de quadrature pres (rectangle : 8e-4 relatif).
                W += 0.5 * (tf + tfPrev).dot(ds - dsPrev);
                dsPrev = ds;
                tfPrev = tf;
                minCum = std::min(minCum, W);
            }
        }
        return W - Wlast;                     // travail du DERNIER cycle
    };
    {
        const double A = 1.0e-6, P = 2.0 * A;
        double minCum = 0.0;
        double Wc = cycleWork(true, A, 3, minCum);
        double Wth = 2.0 * fCapC * (P - 2.0 * fCapC / pj);
        chk(std::fabs(Wc - Wth) < 2e-3 * Wth,
            "(ii) dissipation/cycle = 2 f_cap (P - 2 f_cap/pj)", Wc, Wth);
        chk(Wc > 0.0, "(ii) dissipation/cycle > 0", Wc, 0.0);
        chk(minCum >= -fCapC * fCapC / (2.0 * pj) * (1.0 + 1e-6),
            "(ii) travail cumule >= -f_cap^2/(2 pj)", minCum,
            -fCapC * fCapC / (2.0 * pj));
        // asymptote 2 f_cap P a pj -> inf (ici 2 f_cap/pj = 4,7e-8 << P)
        chk(std::fabs(Wc - 2.0 * fCapC * P) < 0.05 * 2.0 * fCapC * P,
            "(ii) ~ 2 f_cap P (pj grand)", Wc, 2.0 * fCapC * P);
        // variante qui DOIT echouer : sans retour de glissement, 0 dissipe
        double mc2 = 0.0;
        double Wbad = cycleWork(false, A, 3, mc2);
        chk(std::fabs(Wbad) < 1e-4 * Wth,
            "(ii) CONTROLE : ressort cape sans retour dissipe 0", Wbad, 0.0);
    }

    // (iii) continuite en delta_s = 0, etat vierge
    {
        const double eps = 1e-12;
        Eigen::Vector3d s1 = Eigen::Vector3d::Zero(), s2 = Eigen::Vector3d::Zero();
        Eigen::Vector3d tp = camachoFrictionSlider(Eigen::Vector3d(eps, 0, 0),
                                                   pj, fCapC, s1);
        Eigen::Vector3d tm = camachoFrictionSlider(Eigen::Vector3d(-eps, 0, 0),
                                                   pj, fCapC, s2);
        double jump = (tp - tm).norm();
        chk(std::fabs(jump - 2.0 * pj * eps) < 1e-9 * 2.0 * pj * eps,
            "(iii) saut en 0 = 2 pj eps (-> 0)", jump, 2.0 * pj * eps);
        chk(jump < 1e-3 * fCapC, "(iii) saut << f_cap", jump, fCapC);
        // forme RETIREE (traction constante dirigee par le deplacement)
        Eigen::Vector3d dsp(eps, 0, 0), dsm(-eps, 0, 0);
        double jumpOld = ((fCapC / dsp.norm()) * dsp
                          - (fCapC / dsm.norm()) * dsm).norm();
        chk(std::fabs(jumpOld - 2.0 * fCapC) < 1e-9 * fCapC,
            "(iii) CONTROLE : forme retiree saute de 2 f_cap", jumpOld,
            2.0 * fCapC);
    }
    std::printf("%s (%d echec(s))\n", fail ? "ECHEC" : "OK", fail);
    return fail ? 1 : 0;
}
