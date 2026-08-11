// ---------------------------------------------------------------------------
// yan_point — material-point driver for the cohesive law of Yan, Wang, Jiao
// et al. (IJRMMS 169, 2023, 105439), section 2.5.
//
// It links against the SAME header the solver uses (rockim/YanSoftening.hpp),
// so f(D) and I = int_0^1 f(D) dD are the shipped functions, not a copy. The
// driver reproduces the traction update of FdemSolver::jointForces() with
// jointSoftening = yan for a single joint point, and integrates
//
//     GfI  = int_0^ot sigma(o) do      (eq. 13)
//     GfII = int_0^st tau(s)   ds      (eq. 15)
//
// by the trapezoid rule over a monotonic opening (resp. sliding) path. The
// elastic branch of width dnE = ft/pj is NOT part of the fracture energy: the
// article's o is measured from the peak, so the integral starts at dnE.
//
// Output: TSV on stdout, one section per run, so the python side can plot
// figures 4, 5 and 6 straight from the C++ law.
//   build:  clang++ -std=c++17 -O2 -I../include tools/yan_point.cpp -o yan_point
// ---------------------------------------------------------------------------
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "rockim/YanSoftening.hpp"

using rockim::yan::Params;
using rockim::yan::fD;
using rockim::yan::integralFD;

namespace {

// Table 1 of the article
constexpr double FT   = 1.3e6;      // tensile strength [Pa]
constexpr double COH  = 16.4e6;     // cohesion [Pa]
constexpr double PHI  = 23.0;       // friction angle [deg]
constexpr double GFI  = 3.8;        // mode-I fracture energy [N/m]
constexpr double GFII = 84.0;       // mode-II fracture energy [N/m]

// One joint point, jointSoftening = yan. Mirrors the solver's state
// (D = Dmax, omax, slip) and its traction update.
struct Point {
    Params P;
    double pj;                      // joint penalty [Pa/m]
    double ft, coh, tanPhi;
    double dnE, ot, st;             // elastic width, critical opening/slip
    double D = 0.0;                 // = Dmax (irreversible)
    double omax = 0.0;
    double slip = 0.0;
    bool fricScaled = false;

    Point(double pjIn, double ftIn, double cohIn, double phiDeg,
          double gfI, double gfII, Params Pin = Params())
        : P(Pin), pj(pjIn), ft(ftIn), coh(cohIn),
          tanPhi(std::tan(phiDeg * M_PI / 180.0)) {
        double I = integralFD(P);
        dnE = ft / pj;
        ot  = gfI  / (ft  * I);     // eq. 13
        st  = gfII / (coh * I);     // eq. 15
    }

    // normal traction at opening dn (solver's branch, jointSoftening = yan)
    double sigma(double dn) {
        double rn = (dn > dnE) ? (dn - dnE) / ot : 0.0;
        double rs = std::abs(slip) / st;
        double Dnow = std::sqrt(rn * rn + rs * rs);     // eq. 16
        if (Dnow > D) D = std::fmin(1.0, Dnow);         // irreversible
        if (dn < 0.0) return pj * dn;                   // closed: penalty
        if (dn > omax) omax = dn;
        double sMax = std::fmin(pj * omax, fD(D, P) * ft);
        return (omax > 1e-30) ? sMax * dn / omax : 0.0; // eq. 17
    }

    // tangential traction for a prescribed tangential displacement dtg under a
    // fixed normal state, return mapping exactly as the solver does
    double tau(double dtg, double sigN) {
        double tauTr = pj * (dtg - slip);
        double f = fD(D, P);
        double lim = f * coh
                   + (fricScaled ? f : 1.0) * tanPhi * std::fmax(0.0, -sigN);
        double t = std::fmax(-lim, std::fmin(lim, tauTr));
        if (t != tauTr) {
            slip += (tauTr - t) / pj;
            double rn = 0.0;                            // pure shear path
            double rs = std::abs(slip) / st;
            double Dt = std::sqrt(rn * rn + rs * rs);   // eq. 14
            if (Dt > D) D = std::fmin(1.0, Dt);
        }
        return t;
    }
};

}  // namespace

int main(int argc, char** argv) {
    Params P;
    int N = (argc > 1) ? std::atoi(argv[1]) : 2000000;
    double I = integralFD(P);

    // ---- f(D) curve (article fig. 5) --------------------------------------
    std::printf("# SECTION fD\n#D\tf\n");
    for (int i = 0; i <= 500; ++i) {
        double D = i / 500.0;
        std::printf("%.10f\t%.10f\n", D, fD(D, P));
    }
    std::printf("# integral_fD\t%.15f\n", I);

    // ---- mode I: monotonic opening, GfI by trapezoid (eq. 13) -------------
    // pj taken as the solver would with insertionPenaltyFactor = 4 and
    // E = 15 GPa, h = 2 mm: 4 E / h = 3e13 Pa/m.
    const double PJ = 4.0 * 15e9 / 2e-3;
    {
        Point pt(PJ, FT, COH, PHI, GFI, GFII);
        double dnMax = pt.dnE + pt.ot;
        double GI = 0.0, prev = 0.0;
        std::printf("# SECTION modeI\n#o\tsigma\tD\n");
        int stride = N / 800;
        for (int i = 0; i <= N; ++i) {
            double dn = dnMax * i / N;
            double s = pt.sigma(dn);
            if (i > 0) GI += 0.5 * (s + prev) * (dnMax / N);
            prev = s;
            if (i % stride == 0)
                std::printf("%.12e\t%.10e\t%.8f\n", dn, s, pt.D);
        }
        // the elastic branch (0 -> dnE) is stored energy, not fracture energy
        double Gel = 0.5 * pt.ft * pt.dnE;
        std::printf("# dnE\t%.12e\n# ot\t%.12e\n", pt.dnE, pt.ot);
        std::printf("# GfI_total\t%.10f\n# GfI_elastic\t%.10f\n"
                    "# GfI_fracture\t%.10f\n# GfI_target\t%.10f\n",
                    GI, Gel, GI - Gel, GFI);
    }

    // ---- mode II: monotonic sliding at zero normal stress (eq. 15) --------
    {
        Point pt(PJ, FT, COH, PHI, GFI, GFII);
        double sE = pt.coh / pt.pj;             // elastic tangential branch
        double sMax = sE + pt.st;
        double GII = 0.0, prev = 0.0;       // int tau d(dtg), total path
        double GIIs = 0.0, prevS = 0.0, slipPrev = 0.0;  // int tau ds, eq. 15
        std::printf("# SECTION modeII\n#s\ttau\tD\n");
        int stride = N / 800;
        for (int i = 0; i <= N; ++i) {
            double dtg = sMax * i / N;
            double t = pt.tau(dtg, 0.0);        // sigma_n = 0: pure cohesion
            if (i > 0) {
                GII  += 0.5 * (t + prev) * (sMax / N);
                // eq. 15 is written on the article's sliding amount s, which
                // is the PLASTIC slip of the return mapping (the article's
                // extrinsic element has no elastic tangential branch); the
                // difference with int tau d(dtg) is exactly the recoverable
                // c^2 / (2 pj) stored on the elastic branch
                GIIs += 0.5 * (t + prevS) * (pt.slip - slipPrev);
            }
            prev = t; prevS = t; slipPrev = pt.slip;
            if (i % stride == 0)
                std::printf("%.12e\t%.10e\t%.8f\n", dtg, t, pt.D);
        }
        double Gel = 0.5 * pt.coh * sE;
        std::printf("# sE\t%.12e\n# st\t%.12e\n", sE, pt.st);
        std::printf("# GfII_total\t%.10f\n# GfII_elastic\t%.10f\n"
                    "# GfII_fracture\t%.10f\n# GfII_slip\t%.10f\n"
                    "# GfII_target\t%.10f\n",
                    GII, Gel, GII - Gel, GIIs, GFII);
    }

    // ---- load / unload / reload cycle in mode I (article fig. 6, eq. 17) --
    {
        Point pt(PJ, FT, COH, PHI, GFI, GFII);
        double dnMax = pt.dnE + pt.ot;
        double o1 = pt.dnE + 0.30 * pt.ot;      // first unloading point
        double o2 = pt.dnE + 0.62 * pt.ot;      // second unloading point
        std::vector<double> legs = {o1, 0.0, o1, o2, 0.15 * o2, o2, dnMax};
        double cur = 0.0;
        std::printf("# SECTION cycleI\n#o\tsigma\tD\n");
        for (double tgt : legs) {
            int n = 4000;
            for (int i = 1; i <= n; ++i) {
                double dn = cur + (tgt - cur) * i / n;
                std::printf("%.12e\t%.10e\t%.8f\n", dn, pt.sigma(dn), pt.D);
            }
            cur = tgt;
        }
    }

    // ---- load / unload / reload cycle in mode II -------------------------
    {
        Point pt(PJ, FT, COH, PHI, GFI, GFII);
        double sE = pt.coh / pt.pj;
        double sMax = sE + pt.st;
        double s1 = sE + 0.30 * pt.st, s2 = sE + 0.62 * pt.st;
        std::vector<double> legs = {s1, 0.2 * s1, s1, s2, 0.15 * s2, s2, sMax};
        double cur = 0.0;
        std::printf("# SECTION cycleII\n#s\ttau\tD\n");
        for (double tgt : legs) {
            int n = 4000;
            for (int i = 1; i <= n; ++i) {
                double dtg = cur + (tgt - cur) * i / n;
                std::printf("%.12e\t%.10e\t%.8f\n", dtg, pt.tau(dtg, 0.0), pt.D);
            }
            cur = tgt;
        }
    }
    // ---- literal eq. 18 secant cycle in mode II --------------------------
    // rockim's shear is a RETURN MAPPING: unloading follows the pj slope and
    // the slip is frozen, which is not the origin secant of eq. 18. This
    // section evaluates eq. 18 directly (tau = f(Dmax) c |s| / smax) on the
    // same load path, so the figure can show both and the deviation is
    // explicit rather than hidden.
    {
        Params P2;
        double I = integralFD(P2);
        double st = GFII / (COH * I);
        double smax = 0.0, Dmax = 0.0;
        std::vector<double> legs = {0.30 * st, 0.05 * st, 0.30 * st,
                                    0.62 * st, 0.10 * st, 0.62 * st, st};
        double cur = 0.0;
        std::printf("# SECTION cycleII_eq18\n#s\ttau\tD\n");
        for (double tgt : legs) {
            int n = 4000;
            for (int i = 1; i <= n; ++i) {
                double sl = cur + (tgt - cur) * i / n;
                if (sl > smax) smax = sl;
                double D = smax / st;                     // eq. 14
                if (D > Dmax) Dmax = std::fmin(1.0, D);
                double f = fD(Dmax, P2);
                // eq. 18 with sigma_n = 0 (pure cohesion branch)
                double t = (smax > 1e-30) ? COH * f * sl / smax : 0.0;
                std::printf("%.12e\t%.10e\t%.8f\n", sl, t, Dmax);
            }
            cur = tgt;
        }
    }

    return 0;
}
