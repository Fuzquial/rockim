// Banc COURT du frottement en cap de la branche jointTSL = camacho (2D),
// correctif du 2026-09-11 (RETOUR_v3 §1.4 quater, v3 §3.4). Exerce
// rockim::jfric::capReturn (include/rockim/FdemSolver.hpp), la fonction
// PURE que jointForces() appelle — pas une transcription.
//  (i)   trajet en traction pure (f_cap = 0) : tau_fric == 0 a chaque pas,
//        et slip_p suit delta_s (pas de raideur de collage en traction) ;
//  (ii)  cycle de glissement alterne sous compression, amplitude A :
//        travail dissipe par cycle = aire de la boucle d'hysteresis =
//        2 f_cap x (2A - 2 f_cap/pj) (parallelogramme : sur chaque
//        demi-cycle la traction remonte de -f_cap a +f_cap sur 2 f_cap/pj)
//        -> tend vers 2 f_cap x 2A quand pj -> inf ; increments JAMAIS negatifs
//        (tau_fric d(slip_p) >= 0 pas a pas) ;
//  (iii) continuite en delta_s = 0 : |tau(+h) - tau(-h)| -> 2 pj h, et non
//        2 f_cap comme dans la premiere ecriture (fr signe par dtg) ;
//  (iv)  variante qui DOIT echouer : la premiere ecriture (+-f_cap selon le
//        signe de delta_s) est discontinue en 0 — le banc le mesure pour
//        prouver que le critere (iii) discrimine.
#include <cmath>
#include <cstdio>
#include "rockim/FdemSolver.hpp"
using namespace rockim;
int main() {
    int fail = 0;
    auto chk = [&](bool ok, const char* what, double a, double b) {
        std::printf("  %-58s %s   (%.9g vs %.9g)\n", what, ok ? "PASS" : "FAIL", a, b);
        if (!ok) ++fail;
    };
    const double pj = 4.0e12 / 1e-3;     // 4 E/h avec E = 40 GPa... ordre 4e15
    const double mu = std::tan(35.0 * M_PI / 180.0);
    // (i) traction pure : t_n = +5 MPa -> shearCap(0, mu, +5e6, D, ...) = 0
    {
        double slip = 0.0, maxTau = 0.0, maxErr = 0.0;
        for (int i = 0; i <= 1000; ++i) {
            const double ds = 1e-6 * std::sin(2.0 * M_PI * i / 250.0);
            const double fCap = jtsl::shearCap(0.0, mu, 5.0e6, 1.0, true);
            const double tau = jfric::capReturn(pj, ds, fCap, slip);
            maxTau = std::max(maxTau, std::fabs(tau));
            maxErr = std::max(maxErr, std::fabs(slip - ds));
        }
        chk(maxTau == 0.0, "(i) traction pure : tau_fric == 0", maxTau, 0.0);
        chk(maxErr == 0.0, "(i) traction pure : slip_p suit delta_s", maxErr, 0.0);
    }
    // (ii) cycle alterne sous compression t_n = -20 MPa, D = 1
    {
        const double tn = -20.0e6;
        const double fCap = jtsl::shearCap(0.0, mu, tn, 1.0, true);
        chk(std::fabs(fCap - mu * 20.0e6) < 1e-6 * fCap, "(ii) f_cap == mu <-t_n>", fCap, mu * 20.0e6);
        const double A = 1e-5;                     // amplitude (m)
        const int N = 4000;                        // pas par cycle
        double slip = 0.0, dsPrev = 0.0, W = 0.0, Wneg = 0.0, Wcyc = 0.0;
        double slipPrev = 0.0;
        for (int c = 0; c < 4; ++c) {
            Wcyc = 0.0;
            for (int i = 1; i <= N; ++i) {
                const double ds = A * std::sin(2.0 * M_PI * i / N);
                const double tau = jfric::capReturn(pj, ds, fCap, slip);
                // travail du terme frottant = tau x d(slip_p) (le reste est
                // elastique reversible dans le ressort pj) ; on mesure les
                // deux : increment plastique et travail total tau x d(ds)
                const double dWp = tau * (slip - slipPrev);
                if (dWp < 0.0) Wneg += dWp;
                Wcyc += tau * (ds - dsPrev);
                W += dWp;
                dsPrev = ds; slipPrev = slip;
            }
        }
        // cycle stabilise : travail total sur un cycle ferme = travail
        // plastique (ressort revenu a son etat) = aire du parallelogramme
        // 2 f_cap x (2A - 2 f_cap/pj)
        const double Wexp = 2.0 * fCap * (2.0 * A - 2.0 * fCap / pj);
        chk(std::fabs(Wcyc - Wexp) < 2e-3 * Wexp, "(ii) travail par cycle == 2 f_cap x 2A (- retour el.)", Wcyc, Wexp);
        chk(Wcyc > 0.0, "(ii) travail par cycle > 0 (dissipe)", Wcyc, 0.0);
        chk(Wneg == 0.0, "(ii) aucun increment plastique negatif", Wneg, 0.0);
        chk(W > 0.0, "(ii) travail plastique cumule > 0", W, 0.0);
    }
    // (iii) continuite en delta_s = 0 (slip_p = 0, compression)
    {
        const double fCap = jtsl::shearCap(0.0, mu, -20.0e6, 1.0, true);
        const double h = 1e-12;
        double s1 = 0.0, s2 = 0.0;
        const double tp = jfric::capReturn(pj, +h, fCap, s1);
        const double tm = jfric::capReturn(pj, -h, fCap, s2);
        chk(std::fabs((tp - tm) - 2.0 * pj * h) < 1e-9 * 2.0 * pj * h, "(iii) saut en 0 == 2 pj h (continu)", tp - tm, 2.0 * pj * h);
        chk(std::fabs(tp - tm) < 1e-3 * fCap, "(iii) saut en 0 << f_cap", tp - tm, fCap);
        // (iv) variante qui DOIT echouer : la premiere ecriture
        auto old = [&](double ds) { return ds > 0.0 ? fCap : (ds < 0.0 ? -fCap : 0.0); };
        const double jumpOld = old(+h) - old(-h);
        chk(std::fabs(jumpOld - 2.0 * fCap) < 1e-12 * fCap, "(iv) premiere ecriture : saut == 2 f_cap (DISCONTINUE, attendu)", jumpOld, 2.0 * fCap);
    }
    // (v) sous fricMob, D = 0 : cap nul, pas de collage
    {
        double slip = 0.0;
        const double fCap = jtsl::shearCap(0.0, mu, -20.0e6, 0.0, true);
        const double tau = jfric::capReturn(pj, 3e-6, fCap, slip);
        chk(fCap == 0.0 && tau == 0.0 && slip == 3e-6, "(v) D = 0 sous fricMob : tau = 0, slip_p = delta_s", tau, 0.0);
    }
    std::printf("%s (%d echec(s))\n", fail ? "FAIL" : "PASS", fail);
    return fail ? 1 : 0;
}
