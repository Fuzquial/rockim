#pragma once
// ---------------------------------------------------------------------------
// YanSoftening — exponential cohesive softening of Yan, Wang, Jiao et al.,
// "A 2D adaptive finite-discrete element method for simulating fracture and
// fragmentation in geomaterials", Int. J. Rock Mech. Min. Sci. 169 (2023)
// 105439, section 2.5.
//
// The article replaces the linear post-peak branch of the classical Munjiza
// joint by a reduction factor f(D) applied to the strengths (its eq. 9-10):
//
//     sigma = f(D) ft                (o >= 0, tension)
//     tau   = f(D) (c - sigma tanPhi)
//
// with (its eq. 11)
//
//                 /        a + b - 1     D (a + c b)      \
//     f(D) =  | 1 - --------- exp( ------------------- ) | (a(1-D) + b(1-D)^c)
//                 \          a + b     (a + b)(1 - a - b) /
//
// fitted on tensile tests: a = 0.63, b = 1.8, c = 6.0. f(0) = 1 exactly (no
// damage, full strength) and f(1) = 0 exactly (broken element). D follows
// eq. 12 (pure tension, D = o/ot), eq. 14 (pure shear, D = |s|/st) and
// eq. 16 (mixed, D = sqrt((o/ot)^2 + (|s|/st)^2)).
//
// The critical opening ot and slip st are NOT free parameters: they are set
// by the fracture energies (its eq. 13 and 15),
//
//     GfI  = int_0^ot sigma(o) do = ft  ot  int_0^1 f(D) dD
//     GfII = int_0^st tau(s)   ds = c   st  int_0^1 f(D) dD
//
// so ot = GfI / (ft I) and st = GfII / (c I) with I = int_0^1 f(D) dD.
// I is a pure function of (a, b, c); integralFD() evaluates it by composite
// Simpson quadrature once at initialisation. For the article's parameters
// I = 0.386307294744 (converged to 1e-15 by n = 4096).
//
// Header-only and free of any solver state on purpose: the material-point
// driver (tools/yan_point.cpp) exercises the very same functions the solver
// calls, so the fracture-energy verification tests the shipped code and not
// a re-implementation of it.
// ---------------------------------------------------------------------------
#include <cmath>

namespace rockim {
namespace yan {

struct Params {
    double a = 0.63;
    double b = 1.8;
    double c = 6.0;
};

// Reduction factor f(D), eq. 11. Clamped to [0, 1] on D outside [0, 1].
inline double fD(double D, const Params& P) {
    if (D <= 0.0) return 1.0;
    if (D >= 1.0) return 0.0;
    const double s = P.a + P.b;
    // s == 1 is the removable singularity of the exponent; the article's fit
    // has s = 2.43, but guard it so a user-supplied (a, b) cannot produce NaN
    if (std::abs(s) < 1e-12 || std::abs(1.0 - s) < 1e-12)
        return P.a * (1.0 - D) + P.b * std::pow(1.0 - D, P.c);
    const double K = (s - 1.0) / s;
    const double alpha = (P.a + P.c * P.b) / (s * (1.0 - s));
    const double g = 1.0 - K * std::exp(alpha * D);
    const double h = P.a * (1.0 - D) + P.b * std::pow(1.0 - D, P.c);
    const double f = g * h;
    return f < 0.0 ? 0.0 : (f > 1.0 ? 1.0 : f);
}

// I = int_0^1 f(D) dD by composite Simpson (n even). f is C^inf on [0, 1] so
// the O(h^4) rule converges to machine precision well before n = 4096.
inline double integralFD(const Params& P, int n = 4096) {
    if (n < 2) n = 2;
    if (n % 2) ++n;
    const double h = 1.0 / n;
    double s = fD(0.0, P) + fD(1.0, P);
    for (int i = 1; i < n; ++i)
        s += fD(i * h, P) * ((i % 2) ? 4.0 : 2.0);
    return s * h / 3.0;
}

}  // namespace yan
}  // namespace rockim
