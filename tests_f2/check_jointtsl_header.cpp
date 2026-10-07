// Compile-check + banc minimal du noyau JointTsl.hpp apres l'ajout de la
// branche ascendante (2026-09-11). Verifie :
//  (a) rise = 0 => forme initiale (dm0 = 0, offset neutre) ;
//  (b) rise > 0 => traction CONTINUE a l'insertion : split(eq. 19) rend
//      (t_n^ins, t_s^ins) a separation geometrique nulle ;
//  (c) raideur de charge bornee = t_ins / dm0 ;
//  (d) integrale de la branche adoucissante = G_C a rise = 0 ET a rise = 1e-3.
#include <cmath>
#include <cstdio>
#include "rockim/JointTsl.hpp"
#include "rockim/YangDif.hpp"
using namespace rockim;
int main() {
    int fail = 0;
    auto chk = [&](bool ok, const char* what, double a, double b) {
        std::printf("  %-52s %s   (%.9g vs %.9g)\n", what, ok ? "PASS" : "FAIL", a, b);
        if (!ok) ++fail;
    };
    const double tn = 8.0e6, ts = 3.0e6, tn0 = 11.4e6, ts0 = 30e6, GIc = 100, GIIc = 1000, eta = 2.0;
    // (a) rise = 0
    jtsl::Stamp S0 = jtsl::stampInsertion(tn, ts, tn0, ts0, GIc, GIIc, eta, 0.0);
    chk(S0.dm0 == 0.0, "rise=0 : dm0 == 0", S0.dm0, 0.0);
    double dn, ds, dm; jtsl::effOffset(S0, 1e-6, 2e-6, dn, ds, dm);
    chk(std::fabs(dm - jtsl::effSep(1e-6, 2e-6, S0.beta)) < 1e-18, "rise=0 : effOffset == effSep", dm, jtsl::effSep(1e-6, 2e-6, S0.beta));
    // (b) rise = 1e-3 : continuite a l'insertion
    jtsl::Stamp S = jtsl::stampInsertion(tn, ts, tn0, ts0, GIc, GIIc, eta, 1e-3);
    jtsl::effOffset(S, 0.0, 0.0, dn, ds, dm);
    chk(std::fabs(dm - S.dm0) < 1e-15 * S.dm0, "rise>0 : dmEff(0,0) == dm0", dm, S.dm0);
    double tm = jtsl::traction(S, dm, dm);
    chk(std::fabs(tm - S.tmIns) < 1e-9 * S.tmIns, "rise>0 : t_m a l'insertion == t_ins", tm, S.tmIns);
    double tnOut, tsScale; jtsl::split(tm, dm, dn, S.beta, tnOut, tsScale);
    chk(std::fabs(tnOut - tn) < 1e-6 * tn, "rise>0 : eq.19 rend t_n^ins a l'insertion", tnOut, tn);
    chk(std::fabs(tsScale * ds - ts) < 1e-6 * ts, "rise>0 : eq.19 rend t_s^ins a l'insertion", tsScale * ds, ts);
    // (c) raideur bornee
    chk(std::fabs(S.loadingStiffness() - S.tmIns / S.dm0) < 1e-9 * S.loadingStiffness(), "k0 == t_ins/dm0", S.loadingStiffness(), S.tmIns / S.dm0);
    chk(S0.loadingStiffness() == 0.0, "rise=0 : k0 signale NON BORNE (0)", S0.loadingStiffness(), 0.0);
    // (d) integrale sur la branche adoucissante (ouverture geometrique 0 -> dmF)
    for (int r = 0; r < 2; ++r) {
        const jtsl::Stamp& X = r ? S : S0;
        const int N = 200000; double I = 0.0, g0 = 0.0;
        for (int i = 0; i < N; ++i) {
            double g = X.dmF * (i + 0.5) / N;             // ouverture geometrique
            double e = g + X.dm0;                         // effective (charge monotone)
            I += jtsl::traction(X, e, e) * (X.dmF / N);
            (void)g0;
        }
        chk(std::fabs(I - X.Gc) < 1e-4 * X.Gc, r ? "int t dδ == G_C (rise=1e-3)" : "int t dδ == G_C (rise=0)", I, X.Gc);
    }
    // decharge secante bornee : depuis g = 1e-9 m (juste ouvert)
    { double e = 1e-9 + S.dm0; double t = jtsl::traction(S, e, e); double k = t / e;
      chk(k <= S.loadingStiffness() * (1 + 1e-9), "secante depuis 1 nm <= k0", k, S.loadingStiffness());
      // a rise = 0 la secante vaut t_ins/delta : 1/delta, non bornee. Testee a
      // 1 pm (dm0/1e-12 = 2,5e4) — a 1 nm le rapport n'est que dm0/1nm = 24,8
      double e0 = 1e-12; double t0 = jtsl::traction(S0, e0, e0); double k0 = t0 / e0;
      chk(k0 > 1e3 * S.loadingStiffness(), "rise=0 : secante depuis 1 pm >> k0 (non bornee)", k0, S.loadingStiffness()); }
    // enveloppes de shearCap
    chk(jtsl::shearCap(1e6, 0.5, -2e6, 0.3, false) == 1e6 + 0.5 * 2e6, "shearCap yan, off", jtsl::shearCap(1e6, 0.5, -2e6, 0.3, false), 1e6 + 0.5 * 2e6);
    chk(jtsl::shearCap(1e6, 0.5, 5e5, 0.3, false, true, 11e6) == 1e6 + 0.5 * (-std::min(5e5, 11e6)), "shearCap yang == coh + mu*mcFrictionTerm", jtsl::shearCap(1e6, 0.5, 5e5, 0.3, false, true, 11e6), 1e6 + 0.5 * mcFrictionTerm(5e5, 11e6, true));
    // YangDif : surcharge bit-identique
    int yf = 0;
    for (double v : {1e-5, 1e-2, 1.0, 50.0, 1e3, 9e3}) if (difCompressionYang(v) != difCompressionYang(v, 0.07)) ++yf;
    std::printf("  %-52s %s\n", "YangDif : dif(e) == dif(e, 0.07) bit a bit, 6 pts", yf ? "FAIL" : "PASS");
    fail += yf;
    std::printf("%s (%d echec)\n", fail ? "BANC EN ECHEC" : "BANC OK", fail);
    return fail ? 1 : 0;
}
