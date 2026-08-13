#pragma once
// ---------------------------------------------------------------------------
// Contact par POTENTIEL de Munjiza — le coeur geometrique 2D, PUR (aucun etat
// de solveur) : c'est la forme des eq. 2-3 de Yan, Zheng & Wang (IJRMMS 169,
// 2023), qui sont elles-memes la force de contact distribuee de Munjiza
// (The Combined Finite-Discrete Element Method, 2004 ; Munjiza & Andrews 2000).
//
// Potentiel d'un triangle : phi = 3 min(l1, l2, l3), les l_i etant les
// coordonnees barycentriques — 1 au centroide, 0 sur le bord, lineaire par
// morceaux sur les trois sous-triangles centroidaux (la fonction "tente" de
// Munjiza). Force totale sur le contacteur A recouvrant la cible B :
//
//     F_A = p [ grad phi_A - grad phi_B ] integre sur S = A inter B
//         = p  ∮_{dS} (phi_A - phi_B) n dG        (Gauss, n sortant de S)
//
// et F_B = -F_A par construction (3e loi de Newton EXACTE). Le champ est
// CONSERVATIF : l'energie de contact est une fonction d'etat du recouvrement,
// un rebond elastique restitue le travail — c'est le test decisif
// (selftest-potential2d), et c'est ce que le contact penalite noeud-arete
// quasi-plastique du solveur ne peut pas faire par construction.
//
// Integration EXACTE : le long d'une arete de dS, phi_A - phi_B est lineaire
// par morceaux, avec des cassures la ou l'arete traverse une MEDIANE de A ou
// de B (le lieu l_i = l_j ou l'argmin change). On subdivise a chaque
// traversee (6 fonctions lineaires a tester, 3 par triangle) et le trapeze
// est exact sur chaque morceau. Subdiviser a un point ou la fonction est en
// fait lineaire est inoffensif — donc aucun test d'appartenance n'est requis.
//
// Repartition nodale CONSISTANTE : la charge lineique lineaire de chaque
// morceau est lumpee en deux forces d'extremite (regle du trapeze consistant,
// resultante et moment exacts), chacune appliquee au MEME point spatial sur
// A (+) et sur B (-) via leurs coordonnees barycentriques : somme des forces
// nodales = 0 machine, somme des moments = 0 machine.
// ---------------------------------------------------------------------------
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>

namespace rockim {
namespace pot {

using V2 = Eigen::Vector2d;

inline double cross2(const V2& a, const V2& b) {
    return a.x() * b.y() - a.y() * b.x();
}

// Coordonnees barycentriques d'un triangle CCW (den = 2 aire > 0 exige).
struct Bary {
    V2 P0, e1, e2;                         // X = P0 + l1 e1 + l2 e2
    double den;                            // cross(e1, e2) = 2 A
    bool ok;
    void set(const V2& A, const V2& B, const V2& C) {
        P0 = A;
        e1 = B - A;
        e2 = C - A;
        den = cross2(e1, e2);
        ok = den > 1e-300;
    }
    void lam(const V2& X, double l[3]) const {
        V2 d = X - P0;
        double l1 = cross2(d, e2) / den;   // coefficient de e1
        double l2 = cross2(e1, d) / den;   // coefficient de e2
        l[0] = 1.0 - l1 - l2;
        l[1] = l1;
        l[2] = l2;
    }
    double phi(const V2& X) const {        // 3 min(l) : 1 au centroide, 0 au bord
        double l[3];
        lam(X, l);
        return 3.0 * std::min({l[0], l[1], l[2]});
    }
};

// Clip de Sutherland-Hodgman du triangle A par les demi-plans du triangle B
// (les deux CCW). out doit pouvoir contenir 8 sommets ; retourne leur nombre.
inline int clipTriTri(const V2 A[3], const V2 B[3], V2* out) {
    V2 buf[8];
    int n = 3;
    for (int k = 0; k < 3; ++k) out[k] = A[k];
    for (int c = 0; c < 3; ++c) {                     // arete de clip B[c]->B[c+1]
        const V2& C0 = B[c];
        V2 ce = B[(c + 1) % 3] - C0;
        int m = 0;
        for (int i = 0; i < n; ++i) {
            const V2& P = out[i];
            const V2& Q = out[(i + 1) % n];
            double dp = cross2(ce, P - C0);           // >= 0 : interieur (gauche)
            double dq = cross2(ce, Q - C0);
            if (dp >= 0.0) {
                buf[m++] = P;
                if (dq < 0.0) buf[m++] = P + (Q - P) * (dp / (dp - dq));
            } else if (dq >= 0.0) {
                buf[m++] = P + (Q - P) * (dp / (dp - dq));
            }
            if (m > 7) break;                          // garde (degenere)
        }
        n = m;
        for (int i = 0; i < n; ++i) out[i] = buf[i];
        if (n < 3) return 0;
    }
    return n;
}

struct PairForce {
    V2 fA[3], fB[3];                       // forces nodales (p inclus)
    V2 F;                                  // resultante sur A (= -resultante B)
    V2 cen;                                // centroide du recouvrement
    double area;                           // aire du recouvrement
};

// Force de contact par potentiel entre deux triangles CCW aux positions
// courantes. p = penalite [N/m de phi.longueur — memes unites que les autres
// forces nodales du solveur si p contient l'epaisseur]. Retourne false si
// pas de recouvrement (ou triangle degenere/inverse).
inline bool pairForce(const V2 posA[3], const V2 posB[3], double p,
                      PairForce& R) {
    Bary bA, bB;
    bA.set(posA[0], posA[1], posA[2]);
    bB.set(posB[0], posB[1], posB[2]);
    if (!bA.ok || !bB.ok) return false;    // inverse ou degenere : on passe

    V2 S[8];
    int n = clipTriTri(posA, posB, S);
    if (n < 3) return false;

    // aire et centroide (formules du polygone)
    double A2 = 0.0;
    V2 cen(0.0, 0.0);
    for (int i = 0; i < n; ++i) {
        const V2& P = S[i];
        const V2& Q = S[(i + 1) % n];
        double w = cross2(P, Q);
        A2 += w;
        cen += (P + Q) * w;
    }
    if (A2 <= 1e-300) return false;        // recouvrement nul ou retourne
    R.area = 0.5 * A2;
    R.cen = cen / (3.0 * A2);

    for (int k = 0; k < 3; ++k) {
        R.fA[k].setZero();
        R.fB[k].setZero();
    }
    R.F.setZero();

    // g(X) = phi_A - phi_B ; cassures la ou l_i - l_j change de signe
    auto g = [&](const V2& X) { return bA.phi(X) - bB.phi(X); };

    for (int i = 0; i < n; ++i) {
        const V2& Q0 = S[i];
        const V2& Q1 = S[(i + 1) % n];
        V2 e = Q1 - Q0;
        double L = e.norm();
        if (L < 1e-300) continue;
        V2 nrm(e.y() / L, -e.x() / L);     // sortante du polygone CCW

        // parametres de subdivision : traversees des 3 medianes de A et de B
        double ss[16];
        int ns = 0;
        ss[ns++] = 0.0;
        ss[ns++] = 1.0;
        double lA0[3], lA1[3], lB0[3], lB1[3];
        bA.lam(Q0, lA0);
        bA.lam(Q1, lA1);
        bB.lam(Q0, lB0);
        bB.lam(Q1, lB1);
        auto addCross = [&](double d0, double d1) {
            if ((d0 > 0.0 && d1 < 0.0) || (d0 < 0.0 && d1 > 0.0)) {
                double s = d0 / (d0 - d1);
                if (s > 1e-12 && s < 1.0 - 1e-12 && ns < 16) ss[ns++] = s;
            }
        };
        for (int a = 0; a < 3; ++a) {
            int b = (a + 1) % 3;
            addCross(lA0[a] - lA0[b], lA1[a] - lA1[b]);
            addCross(lB0[a] - lB0[b], lB1[a] - lB1[b]);
        }
        std::sort(ss, ss + ns);

        for (int q = 0; q + 1 < ns; ++q) {
            double sa = ss[q], sb = ss[q + 1];
            double dL = (sb - sa) * L;
            if (dL < 1e-300) continue;
            V2 Xa = Q0 + sa * e;
            V2 Xb = Q0 + sb * e;
            double ga = g(Xa), gb = g(Xb);
            // charge lineique lineaire p g(s) nrm, lumpee au trapeze
            // consistant : Fa en Xa, Fb en Xb (resultante ET moment exacts)
            V2 Fa = (p * dL * (2.0 * ga + gb) / 6.0) * nrm;
            V2 Fb = (p * dL * (ga + 2.0 * gb) / 6.0) * nrm;
            R.F += Fa + Fb;
            // application aux MEMES points spatiaux sur A (+) et B (-)
            double la[3], lb[3];
            for (int pt = 0; pt < 2; ++pt) {
                const V2& X = pt ? Xb : Xa;
                const V2& Fp = pt ? Fb : Fa;
                bA.lam(X, la);
                bB.lam(X, lb);
                for (int k = 0; k < 3; ++k) {
                    R.fA[k] += la[k] * Fp;
                    R.fB[k] -= lb[k] * Fp;
                }
            }
        }
    }
    return true;
}

} // namespace pot

// selftest-potential2d — LE test decisif du chantier A3 : collision
// elastique sans frottement entre deux corps rigides triangulaires, frontale
// puis oblique (couple non nul). Un champ conservatif doit restituer le
// travail : |gcWork| / KE0 doit etre au niveau de l'erreur d'integration du
// saute-mouton, des ordres de grandeur sous le contact penalite
// quasi-plastique (qui dissipe ~80 % par construction). Implante dans
// FdemSolver.cpp, appele par main.cpp.
int potentialSelftest(const std::string& csvPath);

} // namespace rockim
