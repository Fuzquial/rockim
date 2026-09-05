// ---------------------------------------------------------------------------
// selftest-pdc3d — BANC EN FORME FERMEE du cutter PDC 3D (ToolPdc3d.hpp).
//
// Sans maillage, sans simulation, < 1 s. Motif T0 : il verifie LE noyau que le
// solveur appelle (pdc3d::query), pas une copie. Chaque famille a sa variante
// qui DOIT echouer — un banc que rien ne peut faire tomber ne prouve rien.
//
//  G1 EXACTITUDE. A l'interieur d'un polygone convexe, la distance au bord est
//     le minimum des distances aux SEGMENTS du bord. On la calcule par la
//     formule point-segment (exacte, sans echantillonnage) et on la compare a
//     la penetration rendue par query() sur 20 000 points tires dans le solide.
//     Critere : ecart <= 1e-12 R.  DOIT ECHOUER : le meme noyau prive de sa
//     contrainte de chanfrein (une geometrie chanfreinee testee comme si elle
//     etait vive) — pres du chanfrein la penetration est alors SURESTIMEE.
//  G2 PROJECTION SUR LA PEAU. p + pen.N doit tomber SUR la frontiere : la
//     distance signee y vaut 0 a 1e-12 R pres. C'est le test conjoint de pen ET
//     de la normale — une normale fausse envoie le point ailleurs que sur la
//     peau. Verifie aussi |N| = 1.
//  G3 LARGEUR DE CONTACT a la surface libre : bissection sur la coordonnee
//     laterale, dans le plan de la face de coupe a la hauteur z = E_z + passe,
//     contre la forme fermee 2 sqrt(R^2 - (d/cos b - R)^2). Cas article
//     (R 6,5 / d 1,016 / b 20 deg) : 7,180 mm ; a b = 0 : 6,979 mm, qui doit
//     coincider avec le 2 sqrt(d(D-d)) du deck 2D.  DOIT ECHOUER : la formule
//     naive 2 sqrt(d(D-d)) appliquee a b = 20 deg (ecart 2,9 %).
//  G4 THEOREME DU PLANCHER. Point le plus bas du solide, par echantillonnage
//     de sa surface : pour b < 0 (orientation article / decks 2D) c'est l'arete
//     E, a 1e-9 pres. Pour b > 0 il descend de t.sin(b) sous l'arete — c'est la
//     bande labouree du deck v3 (0,855 mm), et la raison de `cutterFloor`.
//     La variante b > 0 est le DOIT ECHOUER de « rien sous l'arete ».
//  G5 SIGNE DE LA GARDE. Pour b < 0 le corps est AU-DESSUS du niveau de l'arete
//     (z du centre a mi-epaisseur > E_z) ; pour b > 0 il est en dessous. Le
//     meme signe que le 2D est la meme physique.
// ---------------------------------------------------------------------------
#include "rockim/ToolPdc3d.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

namespace {

using Eigen::Vector3d;
using rockim::pdc3d::Frame;
using rockim::pdc3d::Geom;
using rockim::pdc3d::Query;

struct Seg { double r0, s0, r1, s1; };

// Polygone meridien (rho, s) : face de coupe, chanfrein, tranche, face arriere.
std::vector<Seg> meridian(const Geom& g) {
    const double ht = 0.5 * g.thick;
    std::vector<Seg> e;
    if (g.cham > 0.0) {
        double ga = g.chamDeg * M_PI / 180.0;
        double sc = ht - g.cham * std::tan(ga);
        e.push_back({0.0, ht, g.R - g.cham, ht});
        e.push_back({g.R - g.cham, ht, g.R, sc});
        e.push_back({g.R, sc, g.R, -ht});
    } else {
        e.push_back({0.0, ht, g.R, ht});
        e.push_back({g.R, ht, g.R, -ht});
    }
    e.push_back({g.R, -ht, 0.0, -ht});
    return e;
}

double segDist(double r, double s, const Seg& e) {
    double dx = e.r1 - e.r0, dy = e.s1 - e.s0;
    double L2 = dx * dx + dy * dy;
    double t = L2 > 0 ? ((r - e.r0) * dx + (s - e.s0) * dy) / L2 : 0.0;
    t = std::min(1.0, std::max(0.0, t));
    double px = e.r0 + t * dx, py = e.s0 + t * dy;
    return std::hypot(r - px, s - py);
}

double boundaryDist(const Geom& g, const Vector3d& p, const Frame& f,
                    const Vector3d& edge) {
    const Vector3d C = edge + g.R * f.u;
    const Vector3d M = C - 0.5 * g.thick * f.n;
    const Vector3d d = p - M;
    double s = d.dot(f.n);
    double rho = (d - s * f.n).norm();
    double best = 1e300;
    for (const auto& e : meridian(g)) best = std::min(best, segDist(rho, s, e));
    return best;
}

// Tirage uniforme d'un point INTERIEUR (rejet dans la boite englobante locale).
Vector3d drawInside(const Geom& g, const Frame& f, const Vector3d& edge,
                    std::mt19937_64& rng) {
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    const Vector3d C = edge + g.R * f.u;
    const Vector3d M = C - 0.5 * g.thick * f.n;
    for (;;) {
        Vector3d q = M + (U(rng) * g.R) * f.u + (U(rng) * g.R) * f.w
                     + (U(rng) * 0.5 * g.thick) * f.n;
        if (rockim::pdc3d::query(g, f, edge, q).inside) return q;
    }
}

// Point le plus bas de la SURFACE du solide (echantillonnage dense).
double lowestZ(const Geom& g, const Frame& f, const Vector3d& edge) {
    const Vector3d C = edge + g.R * f.u;
    double zmin = 1e300;
    const int NT = 720, NR = 60, NS = 40;
    for (int it = 0; it < NT; ++it) {
        double th = 2.0 * M_PI * it / NT;
        Vector3d rad = std::cos(th) * f.u + std::sin(th) * f.w;
        for (int is = 0; is <= NS; ++is) {                 // tranche
            double s = -g.thick * is / NS;
            zmin = std::min(zmin, (C + g.R * rad + s * f.n).z());
        }
        for (int ir = 0; ir <= NR; ++ir) {                 // deux faces
            double r = g.R * ir / NR;
            zmin = std::min(zmin, (C + r * rad).z());
            zmin = std::min(zmin, (C + r * rad - g.thick * f.n).z());
        }
    }
    return zmin;
}

} // namespace

int pdc3dSelftest(const std::string& csvPath) {
    using rockim::pdc3d::contactWidth;
    using rockim::pdc3d::frame;
    using rockim::pdc3d::query;

    std::ofstream csv(csvPath);
    csv << "cas,parametre,attendu,obtenu,ecart,verdict\n";
    int fails = 0;
    auto check = [&](const char* cas, const char* what, double exp, double got,
                     double tol, bool mustFail = false) {
        double err = std::abs(got - exp);
        bool ok = err <= tol;
        if (mustFail) ok = !ok;                 // la variante DOIT echouer
        csv << cas << "," << what << "," << exp << "," << got << "," << err
            << "," << (ok ? "PASS" : "FAIL") << "\n";
        std::cout << "[pdc3d] " << (ok ? "[PASS] " : "[FAIL] ") << cas << " / "
                  << what << (mustFail ? " (variante qui DOIT echouer)" : "")
                  << " : attendu " << exp << ", obtenu " << got << ", ecart "
                  << err << "\n";
        if (!ok) ++fails;
        return ok;
    };

    // Geometrie de l'article : R 6,5 mm, t 2,5 mm, chanfrein 0,4 mm a 45 deg,
    // garde arriere 20 deg dans l'orientation des decks 2D (b = -20).
    Geom g;  g.R = 0.0065; g.thick = 0.0025; g.cham = 0.0004; g.chamDeg = 45.0;
    g.rakeDeg = -20.0;
    const Vector3d E(0.0, 0.0, 0.0);
    const Frame f = frame(g.rakeDeg);
    std::mt19937_64 rng(20260903);
    const int N = 20000;

    // ---- G1 : penetration = distance exacte au polygone meridien ----------
    double worst = 0.0, worstNaive = 0.0;
    Geom gv = g; gv.cham = 0.0;               // variante : chanfrein ignore
    for (int i = 0; i < N; ++i) {
        Vector3d p = drawInside(g, f, E, rng);
        double ref = boundaryDist(g, p, f, E);
        Query q = query(g, f, E, p);
        worst = std::max(worst, std::abs(q.pen - ref));
        Query qv = query(gv, f, E, p);          // noyau vif sur un solide chanfreine
        worstNaive = std::max(worstNaive, std::abs(qv.pen - ref));
    }
    check("G1", "max |pen - dist_exacte| / R", 0.0, worst / g.R, 1e-12);
    check("G1", "variante sans chanfrein / R", 0.0, worstNaive / g.R, 1e-12, true);

    // ---- G2 : p + pen N tombe sur la peau, |N| = 1 -------------------------
    double worstSkin = 0.0, worstNorm = 0.0;
    for (int i = 0; i < N; ++i) {
        Vector3d p = drawInside(g, f, E, rng);
        Query q = query(g, f, E, p);
        worstNorm = std::max(worstNorm, std::abs(q.normal.norm() - 1.0));
        Vector3d skin = p + q.pen * q.normal;
        worstSkin = std::max(worstSkin, boundaryDist(g, skin, f, E));
    }
    check("G2", "max dist(p + pen N, peau) / R", 0.0, worstSkin / g.R, 1e-12);
    check("G2", "max | |N| - 1 |", 0.0, worstNorm, 1e-12);

    // ---- G3 : largeur de contact a la surface libre ----------------------
    auto widthNum = [&](double rakeDeg, double depth) {
        Geom gg = g; gg.cham = 0.0; gg.rakeDeg = rakeDeg;
        Frame ff = frame(rakeDeg);
        // un point du plan de la face, juste sous la peau, a la hauteur z = depth
        auto inside = [&](double wv) {
            // p = E + a u + wv w - eps n, avec z(p) = depth : a = (depth + eps n_z)/u_z
            double eps = 1e-9;
            double a = (depth + eps * ff.n.z()) / ff.u.z();
            Vector3d p = E + a * ff.u + wv * ff.w - eps * ff.n;
            return query(gg, ff, E, p).inside;
        };
        double lo = 0.0, hi = gg.R;
        if (!inside(lo)) return 0.0;
        for (int it = 0; it < 200; ++it) {
            double mid = 0.5 * (lo + hi);
            if (inside(mid)) lo = mid; else hi = mid;
        }
        return 2.0 * lo;
    };
    // Tolerance 1e-5 mm = 1e-8 m : le point sonde est place eps = 1e-9 m SOUS
    // la peau de la face, ce qui decale la largeur de ~eps.tan(b) a garde non
    // nulle (1,1e-9 m mesure a 20 deg, 0 a b = 0). Reste 20 000 fois plus
    // serre que l ecart de la variante falsifiante (0,2 mm).
    const double dcut = 0.001016;
    double wA = widthNum(20.0, dcut), wA0 = widthNum(0.0, dcut);
    check("G3", "largeur b=20 deg (mm)", contactWidth(g.R, dcut, 20.0) * 1e3,
          wA * 1e3, 1e-5);
    check("G3", "largeur b=0 (mm)", contactWidth(g.R, dcut, 0.0) * 1e3,
          wA0 * 1e3, 1e-5);
    check("G3", "b=0 = formule 2D 2 sqrt(d(D-d)) (mm)",
          2.0 * std::sqrt(dcut * (2.0 * g.R - dcut)) * 1e3, wA0 * 1e3, 1e-5);
    check("G3", "formule 2D appliquee a b=20 (mm)",
          2.0 * std::sqrt(dcut * (2.0 * g.R - dcut)) * 1e3, wA * 1e3, 1e-5, true);

    // ---- G4 : theoreme du plancher -----------------------------------------
    double zNeg = lowestZ(g, frame(-20.0), E);
    Geom gp = g; gp.rakeDeg = 20.0;
    double zPos = lowestZ(gp, frame(20.0), E);
    check("G4", "b=-20 : point le plus bas = arete (mm)", 0.0, zNeg * 1e3, 1e-6);
    check("G4", "b=+20 : rien sous l arete (mm)", 0.0, zPos * 1e3, 1e-6, true);
    check("G4", "b=+20 : bande labouree t.sin(b) (mm)",
          -g.thick * std::sin(20.0 * M_PI / 180.0) * 1e3, zPos * 1e3, 1e-6);

    // ---- G5 : signe de la garde --------------------------------------------
    double zMidNeg = (E - 0.5 * g.thick * frame(-20.0).n).z();
    double zMidPos = (E - 0.5 * g.thick * frame(20.0).n).z();
    check("G5", "b=-20 : corps au-dessus de l arete (signe)", 1.0,
          zMidNeg > 0 ? 1.0 : -1.0, 0.0);
    check("G5", "b=+20 : corps sous l arete (signe)", -1.0,
          zMidPos > 0 ? 1.0 : -1.0, 0.0);

    std::cout << "[pdc3d] selftest cutter PDC 3D : " << (fails == 0 ? "[PASS]" : "[FAIL]")
              << " (" << fails << " echec(s))\n";
    return fails == 0 ? 0 : 1;
}
