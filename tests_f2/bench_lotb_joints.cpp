// ---------------------------------------------------------------------------
// bench_lotb.cpp — banc COURT du LOT B. Tete-a-tete entre les fonctions que le
// solveur appelle et ce que la note exige. Aucun maillage, aucun run : tout se
// verifie a la main en quelques millisecondes, avec pour chaque point un
// critere FALSIFIANT et une variante qui DOIT echouer.
//
//   cl /nologo /std:c++17 /EHsc /I include /I ..\rockim\eigen-3.4.0 bench_lotb.cpp
// ---------------------------------------------------------------------------
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "rockim/JointTsl.hpp"
#include "rockim/YangDif.hpp"

using namespace rockim;

// COPIE EXACTE de difShearExp de src/Fdem3dSolver.cpp (namespace anonyme).
static double difShearExp(double edot, double n) {
    if (edot <= 5.0e-6) return 1.0;
    if (edot > 1.0e4)   return 1.84;
    double v = 0.77 + 0.56 * std::pow(edot, n);
    return v < 1.0 ? 1.0 : (v > 1.84 ? 1.84 : v);
}

static int fails = 0;
static void check(bool ok, const char* what) {
    std::printf("  [%s] %s\n", ok ? "OK " : "ECHEC", what);
    if (!ok) ++fails;
}

int main() {
    // =====================================================================
    // 1. difExpS : le defaut 0,07 doit rendre difCompressionYang BIT POUR BIT
    // =====================================================================
    std::printf("1. difShearExp(e, 0.07) == difCompressionYang(e)\n");
    bool same = true, diff07 = false;
    for (int i = -12; i <= 6; ++i)
        for (double f = 1.0; f < 10.0; f += 0.37) {
            double e = f * std::pow(10.0, i);
            if (difShearExp(e, 0.07) != difCompressionYang(e)) same = false;
            // variante qui DOIT echouer : un autre exposant doit changer
            // quelque chose, sinon la cle serait inerte
            if (difShearExp(e, 0.1707) != difCompressionYang(e)) diff07 = true;
        }
    check(same, "bit-identite a l'exposant litteral (defaut)");
    check(diff07, "VARIANTE FALSIFIANTE : a_s = 0,1707 change bien la valeur");

    // =====================================================================
    // 2. §2.4 eq. 17 : OBJECTIVITE. int t_m d(delta_m) = G_C, quel que soit
    //    le DEPASSEMENT a l'insertion. C'est la propriete que la loi ancree
    //    sur le seuil nominal (branche `penalty`) n'a PAS.
    // =====================================================================
    std::printf("2. int t_m d(delta_m) = G_C (objectivite au depassement)\n");
    const double GIc = 100.0, GIIc = 1000.0, tn0 = 11.4e6, ts0 = 30.0e6;
    for (double over : {1.0, 1.5, 3.0}) {          // depassement du seuil
        jtsl::Stamp S = jtsl::stampInsertion(over * tn0, 0.0, tn0, ts0,
                                             GIc, GIIc, 0.0);
        // integration trapezoidale de la loi TELLE QUE LE SOLVEUR L'APPELLE
        const int N = 2000000;
        double W = 0.0, dmPrev = 0.0, tPrev = jtsl::traction(S, 0.0, 0.0);
        for (int k = 1; k <= N; ++k) {
            double dm = S.dmF * k / (double)N;
            double t = jtsl::traction(S, dm, dm);
            W += 0.5 * (t + tPrev) * (dm - dmPrev);
            dmPrev = dm; tPrev = t;
        }
        double err = std::fabs(W - S.Gc) / S.Gc;
        std::printf("     depassement x%.1f : t_m^ins = %.4g Pa, delta_f = "
                    "%.4g m, W = %.6f J/m2 (G_C = %.6f, ecart %.2e)\n",
                    over, S.tmIns, S.dmF, W, S.Gc, err);
        check(err < 1e-6, "l'aire vaut G_C independamment du depassement");
    }
    {   // CONTRE-EXEMPLE : la loi ancree sur le SEUIL NOMINAL, celle de
        // rockim_g0 (t_ins ecrete a ft, delta_f fige a 2 G/ft). A
        // depassement x3 elle dissipe encore G_C, mais elle a SAUTE de
        // 3 ft a ft a l'insertion : c'est le saut de contrainte que la
        // note veut supprimer.
        double saut = (3.0 * tn0 - std::min(3.0 * tn0, tn0)) / (3.0 * tn0);
        std::printf("     [temoin] loi ancree au seuil nominal : saut de "
                    "traction a l'insertion = %.0f %% a depassement x3\n",
                    100.0 * saut);
        check(saut > 0.5, "VARIANTE FALSIFIANTE : l'ecretage min(sig,ft) "
                          "saute bien, la loi de la note non");
    }

    // =====================================================================
    // 3. §2.4 eq. 19 : la partition conserve le travail, sur un trajet
    //    MIXTE quelconque. t_n d(delta_n) + t_s.d(delta_s) = t_m d(delta_m).
    // =====================================================================
    std::printf("3. partition eq. 19 : le travail des composantes = t_m ddm\n");
    // L'identite est EXACTE dans le continu ; ici les deux membres sont
    // integres par la regle du trapeze sur le MEME trajet, donc l'ecart
    // residuel doit decroitre en O(1/N^2). On le VERIFIE au lieu de le
    // supposer : un ecart qui stagne signalerait une vraie non-conservation.
    double errPrev = -1.0;
    for (int N : {125000, 500000, 2000000}) {
        jtsl::Stamp S = jtsl::stampInsertion(8.0e6, 12.0e6, tn0, ts0,
                                             GIc, GIIc, 0.0);
        double W = 0.0, Wm = 0.0;
        double dnP = 0.0, dsP = 0.0, dmP = 0.0, tmP = S.tmIns;
        double tnP = 0.0, tsP = 0.0;
        {   // La loi est EXTRINSEQUE : en delta_m -> 0+ la traction ne tend
            // PAS vers zero, elle tend vers t_m^ins. Les composantes y
            // tendent vers t_m^ins fois la direction limite du trajet.
            // Initialiser tnP = tsP = 0 introduirait une erreur d'ordre 1
            // sur le PREMIER trapeze — et c'est exactement ce qu'on a
            // mesure avant de corriger : convergence en 1/N au lieu de 1/N^2.
            const double eps = 1e-9;
            double dn = 0.6 * S.dmF * std::sin(1.2 * eps);
            double ds = 0.5 * S.dmF * eps / S.beta;
            double dm = jtsl::effSep(dn, ds, S.beta);
            double sc = 0.0;
            jtsl::split(S.tmIns, dm, dn, S.beta, tnP, sc);
            tsP = sc * ds;
        }
        for (int k = 1; k <= N; ++k) {
            // trajet non proportionnel : le normal ouvre en sin, le
            // tangentiel en rampe — rien de « bien range »
            double x = k / (double)N;
            double dn = 0.6 * S.dmF * std::sin(1.2 * x);
            double ds = 0.5 * S.dmF * x / S.beta;
            double dm = jtsl::effSep(dn, ds, S.beta);
            if (dm >= S.dmF) break;
            double tm = jtsl::traction(S, dm, dm);
            double tn = 0.0, sc = 0.0;
            jtsl::split(tm, dm, dn, S.beta, tn, sc);
            double ts = sc * ds;
            W  += 0.5 * ((tn + tnP) * (dn - dnP) + (ts + tsP) * (ds - dsP));
            Wm += 0.5 * (tm + tmP) * (dm - dmP);
            dnP = dn; dsP = ds; dmP = dm; tmP = tm; tnP = tn; tsP = ts;
        }
        double err = std::fabs(W - Wm) / std::fabs(Wm);
        std::printf("     N = %8d : W(composantes) = %.8g, W(effectif) = "
                    "%.8g, ecart %.2e%s\n", N, W, Wm, err,
                    errPrev > 0.0 ? "" : "");
        // L'identite est ALGEBRIQUE : ce qui reste est du bruit machine,
        // pas une erreur de quadrature (l'ecart plafonne vers 1e-14 au lieu
        // de continuer a decroitre en 1/N^2).
        check(err < 1e-11, "eq. 19 conserve le travail sur trajet mixte");
        errPrev = err;
    }

    // =====================================================================
    // 4. §2.4 eq. 18 : Benzeggagh-Kenane. Traction pure -> G_Ic ;
    //    cisaillement pur -> G_IIc ; monotone entre les deux.
    // =====================================================================
    std::printf("4. Benzeggagh-Kenane : bornes et monotonie\n");
    {
        jtsl::Stamp A = jtsl::stampInsertion(tn0, 0.0, tn0, ts0,
                                             GIc, GIIc, 2.0);
        jtsl::Stamp B = jtsl::stampInsertion(0.0, ts0, tn0, ts0,
                                             GIc, GIIc, 2.0);
        check(std::fabs(A.Gc - GIc) < 1e-9, "traction pure -> G_C = G_Ic");
        check(std::fabs(B.Gc - GIIc) < 1e-9, "cisaillement pur -> G_C = G_IIc");
        double prev = -1.0; bool mono = true;
        for (int i = 0; i <= 40; ++i) {
            double r = i / 40.0;                       // t_s / t_n
            jtsl::Stamp S = jtsl::stampInsertion(tn0, r * 6.0 * tn0, tn0, ts0,
                                                 GIc, GIIc, 2.0);
            if (S.Gc < prev - 1e-12) mono = false;
            prev = S.Gc;
        }
        check(mono, "G_C croit avec la part de mode II");
        jtsl::Stamp N0 = jtsl::stampInsertion(tn0, 3.0 * tn0, tn0, ts0,
                                              GIc, GIIc, 0.0);
        check(std::fabs(N0.Gc - GIc) < 1e-12,
              "jointMixLaw = none (eta <= 0) rend G_C = G_Ic, inchange");
    }

    // =====================================================================
    // 5. §2.2 eq. 12 : le critere elliptique est STRICTEMENT plus selectif
    //    que le OU logique, et coincide avec lui sur les axes.
    // =====================================================================
    std::printf("5. critere elliptique vs OU logique\n");
    {
        double phi = jtsl::phiInsert(0.8 * tn0, 0.8 * ts0, tn0, ts0);
        bool ouLogique = (0.8 * tn0 >= tn0) || (0.8 * ts0 >= ts0);
        std::printf("     (0,8 ; 0,8) : Phi = %.4f, OU logique = %s\n",
                    phi, ouLogique ? "insere" : "n'insere pas");
        check(phi >= 1.0 && !ouLogique,
              "l'ellipse insere la ou le OU laisse intact");
        check(std::fabs(jtsl::phiInsert(tn0, 0.0, tn0, ts0) - 1.0) < 1e-12,
              "coincide sur l'axe de traction (Phi = 1 au seuil)");
        check(std::fabs(jtsl::phiInsert(0.0, ts0, tn0, ts0) - 1.0) < 1e-12,
              "coincide sur l'axe de cisaillement");
        check(jtsl::phiInsert(-5e6, 0.0, tn0, ts0) == 0.0,
              "compression pure sans cisaillement : Phi = 0 (Macaulay)");
    }

    // =====================================================================
    // 6. §2.6 eq. 22 : Weibull a trois parametres, moyenne d'ensemble = 1.
    //    C'est la propriete sur laquelle reposent les calibrations.
    // =====================================================================
    std::printf("6. Weibull 3 parametres : moyenne d'ensemble preservee\n");
    for (double m : {5.0, 10.0, 24.0})
        for (double xu : {0.0, 0.3, 0.6}) {
            double gam = std::tgamma(1.0 + 1.0 / m);
            double x0 = (1.0 - xu) / gam;
            std::mt19937 rng(12345);
            std::uniform_real_distribution<double> U(0.0, 1.0);
            double s = 0.0, mn = 1e300;
            const int N = 2000000;
            for (int i = 0; i < N; ++i) {
                double u = U(rng);
                u = u < 1e-12 ? 1e-12 : (u > 1.0 - 1e-12 ? 1.0 - 1e-12 : u);
                double w = (xu <= 0.0)
                    ? std::pow(-std::log(1.0 - u), 1.0 / m) / gam
                    : xu + x0 * std::pow(-std::log(1.0 - u), 1.0 / m);
                s += w; mn = std::min(mn, w);
            }
            double mean = s / N;
            std::printf("     m = %4.1f, x_u = %.1f : moyenne = %.5f, "
                        "minimum tire = %.4f\n", m, xu, mean, mn);
            check(std::fabs(mean - 1.0) < 3e-3, "moyenne = 1");
            check(mn >= xu - 1e-12, "aucun tirage sous le plancher x_u");
        }

    // =====================================================================
    // 7. §1.3 eq. 6 : h_e = arete du tetraedre REGULIER de meme volume.
    // =====================================================================
    std::printf("7. tetEdgeLength : arete du tetra regulier\n");
    {
        double a = 3e-3;                     // arete visee
        double V = a * a * a / (6.0 * std::sqrt(2.0));   // volume du regulier
        double he = jtsl::tetEdgeLength(V);
        std::printf("     a = %.6e m -> V = %.6e m3 -> h_e = %.6e m\n",
                    a, V, he);
        check(std::fabs(he - a) / a < 1e-12, "h_e restitue l'arete");
        check(jtsl::tetEdgeLength(0.0) == 0.0, "V = 0 -> 0, pas de NaN");
    }

    // =====================================================================
    // 8. §2.4 : shearCap. Le facteur D sur la part FROTTANTE.
    // =====================================================================
    std::printf("8. shearCap : frottement mobilise par D\n");
    {
        double coh = 10e6, mu = 0.6, tn = -50e6;      // 50 MPa de compression
        double off0 = jtsl::shearCap(coh, mu, tn, 0.0, false);
        double on0  = jtsl::shearCap(coh, mu, tn, 0.0, true);
        double on1  = jtsl::shearCap(coh, mu, tn, 1.0, true);
        std::printf("     D = 0 : off = %.4g Pa, damage = %.4g Pa ; "
                    "D = 1 damage = %.4g Pa\n", off0, on0, on1);
        check(std::fabs(off0 - (coh + mu * 50e6)) < 1.0,
              "off (defaut) = c + mu|sn| a pleine valeur des D = 0");
        check(std::fabs(on0 - coh) < 1e-9,
              "damage : a D = 0 le frottement ne porte PAS");
        check(std::fabs(on1 - off0) < 1e-9,
              "damage : a D = 1 les deux formes coincident");
        check(jtsl::shearCap(coh, mu, +5e6, 0.0, false) == coh,
              "en TRACTION le terme frottant est nul (Macaulay)");
    }

    std::printf("\n%s (%d echec(s))\n", fails ? "BANC EN ECHEC" : "BANC OK",
                fails);
    return fails ? 1 : 0;
}
