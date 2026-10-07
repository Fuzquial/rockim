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

// ---------------------------------------------------------------------------
// (2026-10-07) ENERGIE EXACTE ET FORCES NODALES GRADIENT — correctif 3 de
// ENQUETE_CONTACT.md (§3.2), cle potForceExact, et socle de gcBirth = offset.
//
// L'energie de la fonctionnelle de Munjiza pour la paire (A, B) est
//
//     E = p  int_S (phi_A + phi_B) dA ,     S = A inter B .
//
// Pour des corps RIGIDES, la resultante et le moment de pairForce derivent
// exactement de E. Pour des triangles DEFORMABLES, la repartition nodale de
// pairForce (charges de bord reparties par les barycentriques lambda_a(X))
// vaut
//
//     f_a(Munjiza) = - dE/dx_a  +  p I_A grad(lambda_a)       (noeuds de A)
//     f_b(Munjiza) = - dE/dx_b  +  p I_B grad(lambda_b)       (noeuds de B)
//
// avec I_A = int_S phi_A dA, I_B = int_S phi_B dA. Le terme en plus est
// auto-equilibre (Sum_a grad lambda_a = 0, moment nul) : une « pression »
// interne fictive dont le travail p I_A d(ln aire_A) ne derive d'aucun
// potentiel (boucle fermee sur ddl nodaux : 2,9e-4 p, loop2.cpp de
// l'enquete). Le retrancher rend f = -grad E EXACTEMENT (5e-10 aux
// differences finies, travail sur boucle 1e-16).
//
// I_A, I_B : integration EXACTE par decoupe de S selon les medianes des
// deux triangles (9 regions ou argmin lambda_A = i et argmin lambda_B = j ;
// phi_A = 3 lambda_i y est lineaire, son integrale vaut 3 lambda_i(centroide)
// x aire). Meme decoupe que potx.hpp de l'enquete.
// ---------------------------------------------------------------------------

// Clip du polygone convexe P (n sommets) par le demi-plan g.X + c >= 0.
inline int clipHalfPlane(const V2* P, int n, const V2& g, double c, V2* out) {
    int m = 0;
    for (int i = 0; i < n; ++i) {
        const V2& a = P[i];
        const V2& b = P[(i + 1) % n];
        double da = g.dot(a) + c, db = g.dot(b) + c;
        if (da >= 0.0) {
            out[m++] = a;
            if (db < 0.0) out[m++] = a + (b - a) * (da / (da - db));
        } else if (db >= 0.0) {
            out[m++] = a + (b - a) * (da / (da - db));
        }
        if (m > 14) break;                             // garde (degenere)
    }
    return m;
}

// lambda_i(X) = g[i].X + c[i] pour un triangle CCW (den = 2 aire > 0).
// set() rend false (et ne divise pas) pour un triangle inverse ou degenere,
// avec le meme seuil que Bary (revue M6) : les fonctions publiques
// ci-dessous ne dependent plus de ce que l appelant passe d abord par
// pairForce.
struct AffBary {
    V2 g[3];
    double c[3];
    bool set(const V2 P[3]) {
        double den = cross2(P[1] - P[0], P[2] - P[0]);
        if (!(den > 1e-300)) return false;
        for (int i = 0; i < 3; ++i) {
            const V2& Pj = P[(i + 1) % 3];
            const V2& Pk = P[(i + 2) % 3];
            V2 e = Pk - Pj;                            // l_i = cross(e, X-Pj)/den
            g[i] = V2(-e.y(), e.x()) / den;
            c[i] = -g[i].dot(Pj);
        }
        return true;
    }
    double l(int i, const V2& X) const { return g[i].dot(X) + c[i]; }
};

// I_A = int_S phi_A, I_B = int_S phi_B et l'aire de S. false si S est vide.
inline bool overlapIntegrals(const V2 posA_[3], const V2 posB_[3], double& IA,
                             double& IB, double& area) {
    IA = IB = area = 0.0;
    // REPERE LOCAL (origine au 1er sommet de A) : en coordonnees absolues
    // (|X| ~ 100 m sur le tunnel), lambda = g.X + c et l aire des sous-
    // polygones perdent tout chiffre significatif pour un element fin
    // (h ~ 1e-6 m) — mesure : I_A = -1,46 pour une aire de 2e-13 m2. Les
    // integrales sont invariantes par translation, les gradients aussi.
    const V2 O = posA_[0];
    V2 posA[3], posB[3];
    for (int k = 0; k < 3; ++k) {
        posA[k] = posA_[k] - O;
        posB[k] = posB_[k] - O;
    }
    AffBary fa, fb;
    if (!fa.set(posA) || !fb.set(posB)) return false;  // inverse / degenere
    V2 S[8];
    int n = clipTriTri(posA, posB, S);
    if (n < 3) return false;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            V2 P1[16], P2[16];
            int m = n;
            for (int k = 0; k < n; ++k) P1[k] = S[k];
            // region argmin_A = i : l_k - l_i >= 0 pour k != i ; idem B, j
            for (int k = 0; k < 3 && m >= 3; ++k) {
                if (k == i) continue;
                m = clipHalfPlane(P1, m, fa.g[k] - fa.g[i], fa.c[k] - fa.c[i], P2);
                for (int t = 0; t < m; ++t) P1[t] = P2[t];
            }
            for (int k = 0; k < 3 && m >= 3; ++k) {
                if (k == j) continue;
                m = clipHalfPlane(P1, m, fb.g[k] - fb.g[j], fb.c[k] - fb.c[j], P2);
                for (int t = 0; t < m; ++t) P1[t] = P2[t];
            }
            if (m < 3) continue;
            double A2 = 0.0;
            V2 cen(0.0, 0.0);
            for (int t = 0; t < m; ++t) {
                const V2& P = P1[t];
                const V2& Q = P1[(t + 1) % m];
                double w = cross2(P, Q);
                A2 += w;
                cen += (P + Q) * w;
            }
            if (A2 <= 0.0) continue;
            cen /= 3.0 * A2;
            double a = 0.5 * A2;
            IA += a * 3.0 * fa.l(i, cen);
            IB += a * 3.0 * fb.l(j, cen);
            area += a;
        }
    return true;
}

// Energie exacte de la paire, E = p (I_A + I_B) (0 sans recouvrement).
inline double pairEnergy(const V2 posA[3], const V2 posB[3], double p) {
    double IA, IB, ar;
    if (!overlapIntegrals(posA, posB, IA, IB, ar)) return 0.0;
    // E >= 0 par definition (phi >= 0 sur S) ; sur un recouvrement RASANT,
    // l arrondi des barycentriques peut donner -1e-20 : on le ramene a 0
    return std::max(0.0, p * (IA + IB));
}

// Forces nodales GRADIENT EXACT, en place sur le resultat de pairForce (meme
// p) : f_a -= p I_A grad lambda_a, f_b -= p I_B grad lambda_b. La resultante
// R.F et le moment ne changent pas (terme auto-equilibre). Rend l'energie
// E = p (I_A + I_B) de la paire (sert aussi a gcBirth = offset).
inline double exactGradientForces(const V2 posA[3], const V2 posB[3], double p,
                                  PairForce& R) {
    double IA, IB, ar;
    if (!overlapIntegrals(posA, posB, IA, IB, ar)) return 0.0;
    AffBary fa, fb;
    if (!fa.set(posA) || !fb.set(posB)) return 0.0;   // garde (deja vue)
    for (int k = 0; k < 3; ++k) {
        R.fA[k] -= (p * IA) * fa.g[k];
        R.fB[k] -= (p * IB) * fb.g[k];
    }
    return std::max(0.0, p * (IA + IB));       // meme garde que pairEnergy
}

// gcBirth = offset (2026-10-07) : naissance d'une paire EN recouvrement sans
// creation d'energie. A la naissance on fige e0 = E(S0), l'energie de
// recouvrement preexistante, qui n'a ete payee par aucun travail. La paire
// porte ensuite le potentiel DECALE
//
//     U(E) = (E - e0)^2 / E   si E > e0,    0 sinon,
//
// fonction de E SEULE, donc -grad U = -U'(E) grad E : le champ reste un
// gradient (exactement si les forces nodales le sont, potForceExact). La
// force est CONTINUE (U'(e0) = 0) et rejoint la force pleine quand E >> e0
// (U' = 1 - (e0/E)^2 -> 1). Cliquet absorbant : e0 <- min(e0, E) quand la
// paire se separe, ce qui n'a lieu que sur la branche U = 0 ; aucune energie
// n'est donc creee a geometrie fixee (contrairement a la rampe gcBirthTau,
// dont la relaxation de aRef a geometrie fixee materialise E(S0)).
// ATTENTION (revue C2) : la naissance a lieu a la PREMIERE evaluation ou
// pairForce rend vrai, donc avec une aire > 0 : e0 = E(S1) > 0 pour TOUTE
// paire, y compris une paire nee par approche. Aucune paire n'est donc en
// Munjiza pur sous offset : une paire en contact persistant garde le facteur
// 1 - (e0/E)^2 < 1 et une penetration decalee d'environ celle du premier pas
// (v_rel dt), jusqu'a ce qu'elle se separe (cliquet) — comportement voisin de
// la rampe, qui part aussi de 0. Seul e0 = 0 (inatteignable a la naissance)
// redonnerait Munjiza a l'identique.
// Rend le facteur U'(E) a appliquer aux forces ; met a jour e0 ; U en sortie.
inline double birthOffsetScale(double E, double& e0, double* U = nullptr) {
    if (E <= e0) {                                     // branche libre
        e0 = std::max(0.0, E);                         // cliquet
        if (U) *U = 0.0;
        return 0.0;
    }
    if (e0 <= 0.0) {
        if (U) *U = E;
        return 1.0;
    }
    double r = e0 / E;
    if (U) *U = (E - e0) * (1.0 - r);
    return 1.0 - r * r;
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
// tagIn/tagOut (2026-10-03, potForce = volume) : etiquette d ORIGINE par
// face (0 = face de A, 1 = face de coupe portee par un plan de B). Facultatif :
// sous nullptr (chemin historique) rien n est lu ni ecrit, memes flottants.
inline void clipHalf(Poly3& P, const V3& n, const V3& O, Poly3& out,
                     double tol, const int* tagIn = nullptr,
                     int* tagOut = nullptr) {
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
            if (tagOut) tagOut[out.nF] = tagIn ? tagIn[f] : 0;
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
                if (tagOut) tagOut[out.nF] = 1;
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

// Test d'AXE SEPARATEUR sur les 8 plans de faces, avec CACHE du plan gagnant
// (hint, a la Baraff) : les voisins tangents d'un maillage qui pave l'espace
// dominent les paires candidates, et leur clip complet — qui echoue toujours
// — dominait le cout de la detection (mesure : -6 % seulement en optimisant
// la grille). En regime etabli, une paire tangente coute UN test de plan
// (~16 produits scalaires). STRICTEMENT conservatif : on ne declare separe
// que si tous les sommets de l'autre tet sont du cote exterieur (d >= 0) —
// le recouvrement est alors de mesure nulle et le clip l'aurait rejete :
// bit-neutre par construction. hint = -1 si aucun plan connu.
inline bool separated(const V3 pa[4], const V3 pb[4], int& hint) {
    static const int TF[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
    static const int TE[6][2] = {{0, 1}, {0, 2}, {0, 3},
                                 {1, 2}, {1, 3}, {2, 3}};
    // axes 0-7 : plans de faces (0-3 = A, 4-7 = B) ; axes 8-43 : produits
    // croises d'aretes (8 + 6 iA + iB) — le jeu COMPLET du SAT pour deux
    // polytopes convexes : toute paire disjointe possede un axe separateur
    // dans cet ensemble, donc le clip complet n'est paye que par les paires
    // REELLEMENT en recouvrement. Mesure sans les axes d'aretes : sur la
    // surface deformee de la phase debris, la majorite des paires disjointes
    // echappait aux 8 plans de faces et payait le clip a chaque pas.
    auto test = [&](int axis) {
        if (axis < 8) {
            const V3* S = axis < 4 ? pa : pb;
            const V3* O = axis < 4 ? pb : pa;
            const int f = axis & 3;
            const V3& A = S[TF[f][0]];
            V3 n = (S[TF[f][1]] - A).cross(S[TF[f][2]] - A);
            for (int k = 0; k < 4; ++k)
                if (n.dot(O[k] - A) < 0.0) return false;
            return true;
        }
        const int ia = (axis - 8) / 6, ib = (axis - 8) % 6;
        V3 ea = pa[TE[ia][1]] - pa[TE[ia][0]];
        V3 eb = pb[TE[ib][1]] - pb[TE[ib][0]];
        V3 n = ea.cross(eb);
        double nn = n.squaredNorm();
        if (nn < 1e-300) return false;         // aretes paralleles : axe nul
        double aLo = 1e300, aHi = -1e300, bLo = 1e300, bHi = -1e300;
        for (int k = 0; k < 4; ++k) {
            double da = n.dot(pa[k]);
            double db = n.dot(pb[k]);
            aLo = std::min(aLo, da); aHi = std::max(aHi, da);
            bLo = std::min(bLo, db); bHi = std::max(bHi, db);
        }
        return aHi <= bLo || bHi <= aLo;       // intervalles disjoints
    };
    if (hint >= 0 && test(hint)) return true;
    for (int axis = 0; axis < 44; ++axis) {
        if (axis == hint) continue;
        if (test(axis)) {
            hint = axis;
            return true;
        }
    }
    return false;                           // recouvrement reel (SAT complet)
}

inline bool pairForce(const V3 pa[4], const V3 pb[4], double p,
                      PairForce3& R) {
    Bary4 bA, bB;
    bA.set(pa[0], pa[1], pa[2], pa[3]);
    bB.set(pb[0], pb[1], pb[2], pb[3]);
    if (!bA.ok || !bB.ok) return false;    // inverse ou degenere

    // faces sortantes d'un tet positif : (1,2,3) (0,3,2) (0,1,3) (0,2,1)
    static const int TF[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
    // PING-PONG src/dst : l'ancienne copie `P = Q` apres chaque plan de coupe
    // deplacait ~9 Ko de Poly3 quatre fois par clip — mesure aux compteurs :
    // 51 M de clips VIDES (contacts rasants de la phase debris) payaient
    // ~5 us chacun, 52 % du run percussion. Echanger deux pointeurs produit
    // exactement les memes flottants : bit-neutre par construction.
    Poly3 PA, PB;
    Poly3* src = &PA;
    Poly3* dst = &PB;
    src->clear();
    for (int f = 0; f < 4; ++f) {
        src->nV[src->nF] = 3;
        for (int i = 0; i < 3; ++i) src->v[src->nF][i] = pa[TF[f][i]];
        ++src->nF;
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
        clipHalf(*src, n / nn, A, *dst, tol);
        std::swap(src, dst);
        if (src->nF < 3) return false;     // plus de volume
    }
    const Poly3& P = *src;

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
    // thread_local (2026-09-11) : la boucle des paires du solveur 3D est
    // parallele en phase geometrique ; un scratch statique partage serait une
    // course. En serie, strictement identique.
    static thread_local Frag frag[256], tmp[256];
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
        // ping-pong cur/nxt (audit C, 13/09) : l ancienne recopie tmp -> frag
        // apres chacun des 12 plans coutait jusqu a 96 x 256 x 72 o par face ;
        // les deux tampons s echangent par pointeur. Memes flottants, meme
        // ordre : bit-identique par construction (deja fait sur Poly3).
        Frag* cur = frag;
        Frag* nxt = tmp;
        // subdivision par les 6 + 6 fonctions l_i - l_j
        for (int t = 0; t < 2; ++t) {
            const Bary4& b4 = t ? bB : bA;
            for (int a = 0; a < 4; ++a)
                for (int b = a + 1; b < 4; ++b) {
                    int nt = 0;
                    for (int q = 0; q < nf; ++q) {
                        const Frag& F0 = cur[q];
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
                            if (nt < 256) nxt[nt++] = F0;
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
                                        nxt[nt++] = {buf[0], buf[i], buf[i + 1]};
                            }
                        }
                    }
                    nf = nt;
                    std::swap(cur, nxt);
                }
        }
        // ---- lumping nodal consistant par fragment -----------------------
        for (int q = 0; q < nf; ++q) {
            const V3 vv[3] = {cur[q].a, cur[q].b, cur[q].c};
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

// ---------------------------------------------------------------------------
// potForce = volume (2026-10-03) : force de contact fondee sur le VOLUME de
// recouvrement, Feng et al. (cadre general), forme de Liu, Ma, Liu, Tang &
// Fish, CMAME 395 (2022) 114981, eq. (1)-(4), (9)-(12), (55)-(56).
// Potentiel Phi(V) = 1/2 kn V^2 / V', V' = 2 VA VB / (VA + VB) ; force sur A
// F_A = - dPhi/dx_A = - (kn V / V') sum_{faces de S sur dA} a_i n_i, ou n_i est
// la normale sortante de A (translation rigide de A : dV = sum a_i n_i . dx).
// F_B = - F_A EXACTEMENT : la somme des aires vectorielles d un polyedre clos
// est nulle, donc la part de dB vaut l oppose. Appliquee au CENTROIDE du
// recouvrement et repartie par les fonctions de forme des deux tets (eq. 56).
// Meme decoupage (clip de A par les 4 demi-espaces de B), memes gardes
// (plancher de volume, fermeture) que pairForce ; ce qui disparait, c est
// l integration exacte du potentiel de Munjiza sur les faces subdivisees par
// les 12 plans de cassure, l essentiel du cout en regime de contact persistant.
// Raideur : Liu et al. section 2.6, kn_volume ~ 3 a 10 kn_Munjiza (5 retenu),
// le facteur est porte par l appelant (potVolumeFactor).
// ---------------------------------------------------------------------------
inline bool pairForceVolume(const V3 pa[4], const V3 pb[4], double kn,
                            double Vref, PairForce3& R) {
    Bary4 bA, bB;
    bA.set(pa[0], pa[1], pa[2], pa[3]);
    bB.set(pb[0], pb[1], pb[2], pb[3]);
    if (!bA.ok || !bB.ok || !(Vref > 0.0)) return false;
    static const int TF[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
    Poly3 PA, PB;
    Poly3* src = &PA;
    Poly3* dst = &PB;
    int tA[Poly3::MAXF], tB[Poly3::MAXF];
    int* tsrc = tA;
    int* tdst = tB;
    src->clear();
    for (int f = 0; f < 4; ++f) {
        src->nV[src->nF] = 3;
        for (int i = 0; i < 3; ++i) src->v[src->nF][i] = pa[TF[f][i]];
        tsrc[src->nF] = 0;
        ++src->nF;
    }
    double scale = 0.0;
    for (int k = 1; k < 4; ++k)
        scale = std::max({scale, (pa[k] - pa[0]).norm(),
                          (pb[k] - pb[0]).norm()});
    const double tol = 1e-12 * scale;
    for (int f = 0; f < 4; ++f) {
        const V3& A = pb[TF[f][0]];
        V3 n = (pb[TF[f][1]] - A).cross(pb[TF[f][2]] - A);
        double nn = n.norm();
        if (nn < 1e-300) return false;
        clipHalf(*src, n / nn, A, *dst, tol, tsrc, tdst);
        std::swap(src, dst);
        std::swap(tsrc, tdst);
        if (src->nF < 3) return false;
    }
    const Poly3& P = *src;
    V3 g = V3::Zero();
    int ng = 0;
    for (int f = 0; f < P.nF; ++f)
        for (int i = 0; i < P.nV[f]; ++i) { g += P.v[f][i]; ++ng; }
    g /= ng;
    double vol = 0.0;
    V3 cen = V3::Zero(), closure = V3::Zero(), gradA = V3::Zero();
    double aTot = 0.0;
    for (int f = 0; f < P.nF; ++f) {
        V3 nr = V3::Zero();
        const int m = P.nV[f];
        for (int i = 0; i < m; ++i)
            nr += P.v[f][i].cross(P.v[f][(i + 1) % m]);
        closure += nr;
        aTot += nr.norm();
        if (tsrc[f] == 0) gradA += 0.5 * nr;   // aire vectorielle, face de A
        for (int i = 1; i + 1 < m; ++i) {
            const V3& a = P.v[f][0];
            const V3& b = P.v[f][i];
            const V3& c = P.v[f][i + 1];
            double vt = (a - g).cross(b - g).dot(c - g) / 6.0;
            vol += vt;
            cen += vt * (a + b + c + g) / 4.0;
        }
    }
    if (vol <= 1e-12 * std::min(bA.vol, bB.vol)) return false;
    if (aTot > 0.0 && closure.norm() > 1e-6 * aTot) return false;
    R.vol = vol;
    R.cen = cen / vol;
    // Application FACE PAR FACE (et non au seul centroide, l approximation
    // x_c ~ x_m de Liu et al. eq. 7-8, mesuree a 0,2 % d energie perdue sur
    // la collision oblique du selftest) : chaque face de S recoit
    // -(kn V/V') a_i n_i en son centre d aire, sur le tet qui la porte. La
    // resultante sur A est F_A (eq. 1-2), le moment est EXACTEMENT celui du
    // potentiel (eq. 5-6), et la somme des moments sur la surface close est
    // nulle : quantite de mouvement et moment cinetique conserves.
    const double s = kn * vol / Vref;
    for (int k = 0; k < 4; ++k) { R.fA[k].setZero(); R.fB[k].setZero(); }
    R.F.setZero();
    for (int f = 0; f < P.nF; ++f) {
        const int m = P.nV[f];
        V3 av = V3::Zero(), cf = V3::Zero();
        double at = 0.0;
        for (int i = 1; i + 1 < m; ++i) {
            const V3& a = P.v[f][0];
            const V3& b = P.v[f][i];
            const V3& c = P.v[f][i + 1];
            V3 tr = 0.5 * (b - a).cross(c - a);
            const double ta = tr.norm();
            av += tr;
            cf += ta * (a + b + c) / 3.0;
            at += ta;
        }
        if (at < 1e-300) continue;
        cf /= at;
        const V3 Ff = -s * av;             // sur le tet qui porte la face
        double l[4];
        if (tsrc[f] == 0) {
            bA.lam(cf, l);
            for (int k = 0; k < 4; ++k) R.fA[k] += l[k] * Ff;
            R.F += Ff;
        } else {
            bB.lam(cf, l);
            for (int k = 0; k < 4; ++k) R.fB[k] += l[k] * Ff;
        }
    }
    (void)gradA;
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

// selftest-potcontact2d (2026-10-07) — triangles DEFORMABLES : forces nodales
// exactes (gradient, boucle fermee) et naissance a sommet commun sans
// creation d energie (gcBirth = offset). Implante dans FdemSolver.cpp.
int potentialExactSelftest(const std::string& csvPath);

// selftest-potential3d — le meme test en 3D : deux TETS rigides (6 ddl,
// quaternion implicite via Rodrigues), collision frontale puis oblique.
// Implante dans Fdem3dSolver.cpp.
int potentialSelftest3d(const std::string& csvPath, bool volumeForce = false);

} // namespace rockim
