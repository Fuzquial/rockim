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

// ---------------------------------------------------------------------------
// pot3 — le meme contact par potentiel en 3D : recouvrement TET-TET.
//
// Potentiel d'un tetraedre : phi = 4 min(l0..l3), 1 au centroide, 0 sur les
// faces, lineaire par morceaux sur les quatre sous-tets centroidaux. Force :
//
//     F_A = p [ grad phi_A - grad phi_B ] integre sur le VOLUME S = A inter B
//         = p  ∮_{dS} (phi_A - phi_B) n dG     (Gauss, n sortant de S)
//
// dS est le bord du POLYEDRE de recouvrement : le tet A coupe par les quatre
// demi-espaces du tet B (clip de polyedre convexe, face de coupe reconstruite
// a chaque plan — le tri angulaire du chapeau est valide par convexite).
// Integration EXACTE : chaque face de dS est fan-triangulee puis subdivisee
// par les 12 plans de cassure (l_i = l_j de A et de B, la ou l'argmin du min
// change) ; sur chaque fragment les deux phi sont lineaires et le lumping
// nodal CONSISTANT (F_i = p n Aire (2 g_i + g_j + g_k) / 12 aux sommets du
// fragment) reproduit resultante ET moment exactement. Chaque force est
// appliquee au MEME point spatial sur A (+) et sur B (-) via les
// barycentriques des deux tets : 3e loi machine, comme en 2D.
// ---------------------------------------------------------------------------
namespace pot3 {

using V3 = Eigen::Vector3d;

// Barycentriques d'un tet a orientation positive (V = det/6 > 0).
struct Bary4 {
    V3 P0;
    Eigen::Matrix3d Minv;                  // (l1,l2,l3) = Minv (X - P0)
    double vol = 0.0;                      // volume du tet (det/6)
    bool ok;
    void set(const V3& A, const V3& B, const V3& C, const V3& D) {
        P0 = A;
        Eigen::Matrix3d M;
        M.col(0) = B - A;
        M.col(1) = C - A;
        M.col(2) = D - A;
        double det = M.determinant();
        vol = det / 6.0;
        ok = det > 1e-300;
        if (ok) Minv = M.inverse();
    }
    void lam(const V3& X, double l[4]) const {
        V3 q = Minv * (X - P0);
        l[1] = q[0];
        l[2] = q[1];
        l[3] = q[2];
        l[0] = 1.0 - q[0] - q[1] - q[2];
    }
    double phi(const V3& X) const {
        double l[4];
        lam(X, l);
        return 4.0 * std::min({l[0], l[1], l[2], l[3]});
    }
};

// Polyedre convexe = liste de faces polygonales (sommets CCW vus de
// l'EXTERIEUR). Capacites fixes largement suffisantes pour tet coupe par
// quatre plans.
struct Poly3 {
    static constexpr int MAXF = 16, MAXV = 24;
    int nF = 0;
    int nV[MAXF];
    V3 v[MAXF][MAXV];
    void clear() { nF = 0; }
};

// Coupe le polyedre par le demi-espace n.(X - O) <= 0 (on garde l'interieur
// de la face de B d'outward n) et referme par la face de coupe (normale +n).
inline void clipHalf(Poly3& P, const V3& n, const V3& O, Poly3& out,
                     double tol) {
    out.clear();
    V3 cut[4 * Poly3::MAXF];
    int nCut = 0;
    for (int f = 0; f < P.nF; ++f) {
        int m = P.nV[f];
        V3* poly = P.v[f];
        int k = 0;
        V3 buf[Poly3::MAXV];
        for (int i = 0; i < m; ++i) {
            const V3& A = poly[i];
            const V3& B = poly[(i + 1) % m];
            double da = n.dot(A - O), db = n.dot(B - O);
            if (da <= 0.0) {
                if (k < Poly3::MAXV) buf[k++] = A;
                // un sommet garde RASANT (|d| < tol) appartient aussi a la
                // face de coupe : sans lui, le chapeau d'un clip tangent au
                // sommet manque un coin et le polyedre ne se referme pas —
                // c'etait la source du 5e-3 de la frontale miroir
                if (std::abs(da) < tol && nCut < 4 * Poly3::MAXF)
                    cut[nCut++] = A;
                if (db > 0.0) {                        // sortie
                    V3 X = A + (B - A) * (da / (da - db));
                    if (k < Poly3::MAXV) buf[k++] = X;
                    if (nCut < 4 * Poly3::MAXF) cut[nCut++] = X;
                }
            } else if (db <= 0.0) {                    // entree
                V3 X = A + (B - A) * (da / (da - db));
                if (k < Poly3::MAXV) buf[k++] = X;
                if (nCut < 4 * Poly3::MAXF) cut[nCut++] = X;
            }
        }
        if (k >= 3 && out.nF < Poly3::MAXF) {
            out.nV[out.nF] = k;
            for (int i = 0; i < k; ++i) out.v[out.nF][i] = buf[i];
            ++out.nF;
        }
    }
    // face de coupe : points d'intersection ordonnes par angle autour de
    // leur centre dans le plan (convexite => tri valide), dedup RELATIF a la
    // taille du chapeau.
    if (nCut >= 3 && out.nF < Poly3::MAXF) {
        V3 c = V3::Zero();
        for (int i = 0; i < nCut; ++i) c += cut[i];
        c /= nCut;
        double rmax = 0.0;
        for (int i = 0; i < nCut; ++i)
            rmax = std::max(rmax, (cut[i] - c).norm());
        V3 u = V3::Zero();
        for (int i = 0; i < nCut; ++i)     // premier point non confondu
            if ((cut[i] - c).norm() > 0.5 * rmax) { u = cut[i] - c; break; }
        double un = u.norm();
        if (un > 1e-300) {
            u /= un;
            V3 w = n.cross(u);
            struct AP { double a; V3 x; };
            AP ap[4 * Poly3::MAXF];
            int na = 0;
            for (int i = 0; i < nCut; ++i) {
                V3 d = cut[i] - c;
                ap[na++] = {std::atan2(d.dot(w), d.dot(u)), cut[i]};
            }
            std::sort(ap, ap + na, [](const AP& x, const AP& y) {
                return x.a < y.a;
            });
            int k = 0;
            V3 buf[Poly3::MAXV];
            double dd = 1e-9 * rmax;
            for (int i = 0; i < na; ++i) {
                if (k > 0 && (ap[i].x - buf[k - 1]).norm() < dd) continue;
                if (k < Poly3::MAXV) buf[k++] = ap[i].x;
            }
            while (k > 1 && (buf[k - 1] - buf[0]).norm() < dd) --k;
            if (k >= 3) {
                out.nV[out.nF] = k;
                for (int i = 0; i < k; ++i) out.v[out.nF][i] = buf[i];
                ++out.nF;
            }
        }
    }
}

struct PairForce3 {
    V3 fA[4], fB[4];                       // forces nodales (p inclus)
    V3 F;                                  // resultante sur A
    V3 cen;                                // centroide du VOLUME de S
    double vol;                            // volume de S
};

inline bool pairForce(const V3 pa[4], const V3 pb[4], double p,
                      PairForce3& R) {
    Bary4 bA, bB;
    bA.set(pa[0], pa[1], pa[2], pa[3]);
    bB.set(pb[0], pb[1], pb[2], pb[3]);
    if (!bA.ok || !bB.ok) return false;    // inverse ou degenere

    // faces sortantes d'un tet positif : (1,2,3) (0,3,2) (0,1,3) (0,2,1)
    static const int TF[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
    Poly3 P, Q;
    P.clear();
    for (int f = 0; f < 4; ++f) {
        P.nV[P.nF] = 3;
        for (int i = 0; i < 3; ++i) P.v[P.nF][i] = pa[TF[f][i]];
        ++P.nF;
    }
    // tolerance de rasance RELATIVE a la taille des tets
    double scale = 0.0;
    for (int k = 1; k < 4; ++k)
        scale = std::max({scale, (pa[k] - pa[0]).norm(),
                          (pb[k] - pb[0]).norm()});
    double tol = 1e-12 * scale;
    for (int f = 0; f < 4; ++f) {          // demi-espaces de B
        const V3& A = pb[TF[f][0]];
        V3 n = (pb[TF[f][1]] - A).cross(pb[TF[f][2]] - A);
        double nn = n.norm();
        if (nn < 1e-300) return false;
        clipHalf(P, n / nn, A, Q, tol);
        P = Q;
        if (P.nF < 3) return false;        // plus de volume
    }

    // volume + centroide par decoupage en tets depuis le barycentre, et
    // controle de FERMETURE : pour un polyedre clos, la somme des normales
    // ponderees par l'aire est nulle. Deux tets exactement TANGENTS (voisins
    // par arete ou sommet d'un maillage qui pave l'espace) produisent des
    // slivers a volume quasi nul mais a GRANDES faces mal refermees — sans
    // ces gardes, le residu de fermeture donnait des kN de force parasite
    // au repos (5 joints casses a charge nulle sur la grille 3D, attrape par
    // le controle zeroload).
    V3 g = V3::Zero();
    int ng = 0;
    for (int f = 0; f < P.nF; ++f)
        for (int i = 0; i < P.nV[f]; ++i) {
            g += P.v[f][i];
            ++ng;
        }
    g /= ng;
    double vol = 0.0;
    V3 cen = V3::Zero();
    V3 closure = V3::Zero();
    double aTot = 0.0;
    for (int f = 0; f < P.nF; ++f) {
        V3 nr = V3::Zero();
        int m = P.nV[f];
        for (int i = 0; i < m; ++i)
            nr += P.v[f][i].cross(P.v[f][(i + 1) % m]);
        closure += nr;
        aTot += nr.norm();
        for (int i = 1; i + 1 < m; ++i) {
            const V3& a = P.v[f][0];
            const V3& b = P.v[f][i];
            const V3& c = P.v[f][i + 1];
            double vt = (a - g).cross(b - g).dot(c - g) / 6.0;
            vol += vt;
            cen += vt * (a + b + c + g) / 4.0;
        }
    }
    // plancher PHYSIQUE : un recouvrement reel a un volume mesurable devant
    // celui des tets ; les contacts de mesure nulle sont rejetes
    if (vol <= 1e-12 * std::min(bA.vol, bB.vol)) return false;
    if (aTot > 0.0 && closure.norm() > 1e-6 * aTot) return false;
    R.vol = vol;
    R.cen = cen / vol;
    for (int k = 0; k < 4; ++k) {
        R.fA[k].setZero();
        R.fB[k].setZero();
    }
    R.F.setZero();

    // ---- integration exacte face par face --------------------------------
    // fragments triangulaires subdivises par les 12 plans de cassure ;
    // scratch fixe (les fragments d'une face de tet coupe restent peu
    // nombreux — garde par saturation, la force reste bornee par p A phi<=1)
    struct Frag { V3 a, b, c; };
    static Frag frag[256], tmp[256];
    auto gval = [&](const V3& X) { return bA.phi(X) - bB.phi(X); };

    for (int f = 0; f < P.nF; ++f) {
        int m = P.nV[f];
        // normale sortante de la face (polygone plan, CCW exterieur)
        V3 nr = V3::Zero();
        for (int i = 0; i < m; ++i)
            nr += P.v[f][i].cross(P.v[f][(i + 1) % m]);
        double a2 = nr.norm();
        if (a2 < 1e-300) continue;
        V3 n = nr / a2;
        int nf = 0;
        for (int i = 1; i + 1 < m; ++i)
            if (nf < 256) frag[nf++] = {P.v[f][0], P.v[f][i], P.v[f][i + 1]};
        // subdivision par les 6 + 6 fonctions l_i - l_j
        for (int t = 0; t < 2; ++t) {
            const Bary4& b4 = t ? bB : bA;
            for (int a = 0; a < 4; ++a)
                for (int b = a + 1; b < 4; ++b) {
                    int nt = 0;
                    for (int q = 0; q < nf; ++q) {
                        const Frag& F0 = frag[q];
                        const V3 vv[3] = {F0.a, F0.b, F0.c};
                        double d[3];
                        for (int i = 0; i < 3; ++i) {
                            double l[4];
                            b4.lam(vv[i], l);
                            d[i] = l[a] - l[b];
                        }
                        // decoupe du triangle par d = 0 (les deux cotes)
                        bool pos = d[0] > 0 || d[1] > 0 || d[2] > 0;
                        bool neg = d[0] < 0 || d[1] < 0 || d[2] < 0;
                        if (!(pos && neg)) {
                            if (nt < 256) tmp[nt++] = F0;
                            continue;
                        }
                        // polygone coupe en deux : clip contre d>=0 et d<=0
                        for (int side = 0; side < 2; ++side) {
                            V3 buf[4];
                            int k = 0;
                            for (int i = 0; i < 3; ++i) {
                                const V3& A = vv[i];
                                const V3& B = vv[(i + 1) % 3];
                                double da = side ? -d[i] : d[i];
                                double db = side ? -d[(i + 1) % 3]
                                                 : d[(i + 1) % 3];
                                if (da >= 0.0) {
                                    if (k < 4) buf[k++] = A;
                                    if (db < 0.0 && k < 4)
                                        buf[k++] = A + (B - A) * (da / (da - db));
                                } else if (db >= 0.0 && k < 4) {
                                    buf[k++] = A + (B - A) * (da / (da - db));
                                }
                            }
                            if (k >= 3) {
                                for (int i = 1; i + 1 < k; ++i)
                                    if (nt < 256)
                                        tmp[nt++] = {buf[0], buf[i], buf[i + 1]};
                            }
                        }
                    }
                    nf = nt;
                    for (int q = 0; q < nf; ++q) frag[q] = tmp[q];
                }
        }
        // ---- lumping nodal consistant par fragment -----------------------
        for (int q = 0; q < nf; ++q) {
            const V3 vv[3] = {frag[q].a, frag[q].b, frag[q].c};
            double A2f = (vv[1] - vv[0]).cross(vv[2] - vv[0]).norm();
            if (A2f < 1e-300) continue;
            double Af = 0.5 * A2f;
            double gv[3] = {gval(vv[0]), gval(vv[1]), gval(vv[2])};
            for (int i = 0; i < 3; ++i) {
                double w = Af * (2.0 * gv[i] + gv[(i + 1) % 3]
                                 + gv[(i + 2) % 3]) / 12.0;
                V3 Fp = (p * w) * n;
                R.F += Fp;
                double la[4], lb[4];
                bA.lam(vv[i], la);
                bB.lam(vv[i], lb);
                for (int k = 0; k < 4; ++k) {
                    R.fA[k] += la[k] * Fp;
                    R.fB[k] -= lb[k] * Fp;
                }
            }
        }
    }
    return true;
}

} // namespace pot3

// selftest-potential2d — LE test decisif du chantier A3 : collision
// elastique sans frottement entre deux corps rigides triangulaires, frontale
// puis oblique (couple non nul). Un champ conservatif doit restituer le
// travail : |gcWork| / KE0 doit etre au niveau de l'erreur d'integration du
// saute-mouton, des ordres de grandeur sous le contact penalite
// quasi-plastique (qui dissipe ~80 % par construction). Implante dans
// FdemSolver.cpp, appele par main.cpp.
int potentialSelftest(const std::string& csvPath);

// selftest-potential3d — le meme test en 3D : deux TETS rigides (6 ddl,
// quaternion implicite via Rodrigues), collision frontale puis oblique.
// Implante dans Fdem3dSolver.cpp.
int potentialSelftest3d(const std::string& csvPath);

} // namespace rockim
