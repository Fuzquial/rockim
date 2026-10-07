#pragma once
// ---------------------------------------------------------------------------
// Guards.hpp — gardes des ENTREES (C3) et detecteur de NaN/Inf REEL (C4).
// Plan de robustesse du 2026-09-05, chantier C, binaire w20.
//
// C3 (maillage) — trois erreurs NOMMEES, levees a la construction du maillage,
// avant le premier pas, dans TOUS les solveurs qui importent ou construisent
// un maillage :
//   * noeud orphelin : jamais reference par un element. Sur un import Gmsh il
//     recoit une masse nulle et etait, jusqu'a w19, epingle FIXED en silence
//     (la « broche fantome » du 05/09 : point du champ de taille (24, 24, 32)
//     sous le pole d'impact de tous les maillages T1/Q1, 5,9 kN sur un pic de
//     45,6). Il fausse aussi la boite englobante. Aucun maillage a noeud
//     orphelin n'a de sens : erreur, jamais de reparation silencieuse.
//   * masse nodale nulle ou negative apres lumping : un noeud reference par un
//     element mais sans masse (rho <= 0, volume nul) diviserait l'integrateur.
//   * element degenere : volume (aire) <= 0 apres reparation d'orientation, ou
//     < 1e-6 x la mediane du maillage (sliver qui ecrase la CFL et l'inverse
//     du jacobien). L'ancien « degenerate tet » sans identifiant est remplace.
//
// C4 (NaN) — checkFinite() balaie TOUTES les composantes de u, v (et f) de
// TOUS les noeuds tous les nanCheckEvery pas (defaut 256), nomme le premier
// noeud non fini, sa position de reference, le champ et la composante, et
// l'element voisin, puis leve NanError : main.cpp l'attrape, ecrit
// <outputDir>/ERROR.txt et rend le code 3. Remplace le detecteur AVEUGLE de
// Fem3dSolver (u_[0] : noeud 0 possiblement FIXED donc toujours fini) et
// l'echantillonnage a 1/256 des fdem 2D/3D. Lecture PURE : aucun flottant du
// calcul n'est touche, bit-identite par construction.
// ---------------------------------------------------------------------------
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace rockim {

// Exception dediee : reconnue par main.cpp (ERROR.txt + code de retour 3).
struct NanError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

namespace guards {

// « pas d'identifiants d'origine » / « tous les noeuds » (arguments par defaut)
inline const std::vector<long> kNoIds{};
inline const std::vector<char> kNoMask{};

// (x, y[, z]) d'un vecteur Eigen 2D ou 3D
template <class Vec>
inline std::string coords(const Vec& p) {
    std::ostringstream os;
    os.precision(9);
    os << "(";
    for (int k = 0; k < (int)p.size(); ++k) os << (k ? ", " : "") << p[k];
    os << ")";
    return os.str();
}

inline std::string num(double x) {
    std::ostringstream os;
    os.precision(9);
    os << x;
    return os.str();
}

// ---- C3a : noeuds orphelins ------------------------------------------------
// X : positions (indices internes 0..n-1) ; conn : connectivite (indices
// internes) ; ids : identifiants d'origine (Gmsh) ou vide (indice interne +
// 1 est alors imprime tel quel avec la mention « indice »). Leve une erreur
// nommee des qu'un noeud n'est reference par aucun element.
template <class Vec, std::size_t N>
inline void checkOrphans(const std::vector<Vec>& X,
                         const std::vector<std::array<int, N>>& conn,
                         const std::vector<long>& ids) {
    std::vector<char> used(X.size(), 0);
    for (const auto& c : conn)
        for (int a : c)
            if (a >= 0 && (std::size_t)a < used.size()) used[a] = 1;
    std::size_t nOrph = 0;
    long first = -1;
    for (std::size_t i = 0; i < used.size(); ++i)
        if (!used[i]) {
            if (first < 0) first = (long)i;
            ++nOrph;
        }
    if (nOrph == 0) return;
    std::ostringstream os;
    os << "mesh: " << nOrph << " noeud" << (nOrph > 1 ? "s" : "")
       << " orphelin" << (nOrph > 1 ? "s" : "") << " (jamais reference par un "
          "element), premier : id ";
    if (!ids.empty()) os << ids[(std::size_t)first];
    else os << (first + 1) << " (indice interne " << first << ")";
    os << ", " << coords(X[(std::size_t)first])
       << " ; nettoyez le maillage (meshes/drop_orphans.py)";
    throw std::runtime_error(os.str());
}

// ---- C3b : masse nodale nulle ou negative apres lumping --------------------
// referenced : si non vide, seuls les noeuds marques sont controles (les
// noeuds jamais references sont du ressort de checkOrphans, ou d'une
// geometrie taillee qui les epingle en connaissance de cause).
template <class Vec>
inline void checkMasses(const char* tag, const std::vector<double>& m,
                        const std::vector<Vec>& X,
                        const std::vector<char>& referenced,
                        const std::vector<long>& ids) {
    for (std::size_t i = 0; i < m.size(); ++i) {
        if (!referenced.empty() && !referenced[i]) continue;
        if (m[i] > 0.0 && std::isfinite(m[i])) continue;
        std::ostringstream os;
        os << "mesh: masse nodale nulle ou negative apres lumping (" << tag
           << ") : noeud id ";
        if (!ids.empty()) os << ids[i];
        else os << (i + 1) << " (indice interne " << i << ")";
        os << ", " << coords(X[i]) << ", m = " << num(m[i])
           << " kg ; densite rho <= 0 ou element de volume nul";
        throw std::runtime_error(os.str());
    }
}

// masque « reference par au moins un element »
template <std::size_t N>
inline std::vector<char> referencedMask(std::size_t nNodes,
                                        const std::vector<std::array<int, N>>& conn) {
    std::vector<char> used(nNodes, 0);
    for (const auto& c : conn)
        for (int a : c)
            if (a >= 0 && (std::size_t)a < nNodes) used[a] = 1;
    return used;
}

// ---- C3c : element degenere -------------------------------------------------
// Erreur nommee pour UN element : indice interne, ses noeuds (id d'origine
// Gmsh si ids non vide, sinon indice + 1) et leurs coordonnees, son volume
// (aire) et, si connue (med > 0), la mediane du maillage.
template <class Vec, std::size_t N>
[[noreturn]] inline void degenerateError(const char* tag, long e,
                                         const std::array<int, N>& conn,
                                         const std::vector<Vec>& X,
                                         const std::vector<long>& ids,
                                         double vol, double med = -1.0,
                                         double relTol = 1e-6) {
    const char* what = N == 4 ? "tetraedre" : "triangle";
    const char* meas = N == 4 ? "volume" : "aire";
    std::ostringstream os;
    os << "mesh: " << what << " degenere (" << tag << ") : element " << e
       << " (indice interne), noeuds";
    for (std::size_t a = 0; a < N; ++a) {
        int n = conn[a];
        os << (a ? ", " : " ") << "id ";
        if (!ids.empty() && n >= 0 && (std::size_t)n < ids.size())
            os << ids[(std::size_t)n];
        else os << (n + 1);
        if (n >= 0 && (std::size_t)n < X.size()) os << " " << coords(X[(std::size_t)n]);
    }
    os << " ; " << meas << " = " << num(vol);
    if (med > 0.0)
        os << " (mediane du maillage " << num(med) << ", seuil "
           << num(relTol * med) << ")";
    os << " ; " << (!(vol > 0.0) ? "noeuds coplanaires ou confondus"
                                  : "sliver : remaillez (Gmsh Mesh.Optimize)");
    throw std::runtime_error(os.str());
}

// vol : volume (3D) ou aire (2D) de chaque element APRES reparation de
// l'orientation ; conn : connectivite ; X : positions ; ids : identifiants
// d'origine des noeuds (Gmsh) ou vide. Erreur si vol <= 0 ou
// vol < relTol x mediane (defaut 1e-6).
template <class Vec, std::size_t N>
inline void checkDegenerate(const char* tag, const std::vector<double>& vol,
                            const std::vector<std::array<int, N>>& conn,
                            const std::vector<Vec>& X,
                            const std::vector<long>& ids,
                            double relTol = 1e-6) {
    if (vol.empty()) return;
    std::vector<double> s(vol);
    std::nth_element(s.begin(), s.begin() + (std::ptrdiff_t)(s.size() / 2), s.end());
    const double med = s[s.size() / 2];
    for (std::size_t e = 0; e < vol.size(); ++e) {
        const bool bad = !(vol[e] > 0.0) || !std::isfinite(vol[e])
                         || vol[e] < relTol * med;
        if (bad) degenerateError(tag, (long)e, conn[e], X, ids, vol[e], med, relTol);
    }
}

// ---- C4 : NaN / Inf ---------------------------------------------------------
// n noeuds ; field(i, k) rend le vecteur k (0 = u, 1 = v, 2 = f...) du noeud
// i ; names[k] son nom ; nFields leur nombre ; X0(i) la position de reference
// du noeud i ; elemOf(i) un element contenant le noeud (-1 = inconnu) ;
// scalars : grandeurs globales (travail, energies) a controler aussi.
template <class Field, class Pos, class ElemOf>
inline void checkFinite(const char* tag, long step, double t, std::size_t n,
                        int nFields, const char* const* names, Field field,
                        Pos X0, ElemOf elemOf,
                        const std::vector<std::pair<const char*, double>>& scalars) {
    for (std::size_t i = 0; i < n; ++i)
        for (int k = 0; k < nFields; ++k) {
            const auto& a = field(i, k);
            for (int c = 0; c < (int)a.size(); ++c) {
                const double x = a[c];
                if (std::isfinite(x)) continue;
                std::ostringstream os;
                os.precision(9);
                os << tag << " : NaN/Inf detecte au pas " << step << " (t = " << t
                   << " s) : noeud " << i << " (X0 = " << coords(X0(i))
                   << "), champ " << names[k] << " composante " << "xyz"[c]
                   << " = " << x;
                long e = elemOf(i);
                if (e >= 0) os << ", element voisin " << e;
                os << " ; reduisez dtFactor (un nanCheckEvery plus petit "
                      "localise le pas plus finement)";
                throw NanError(os.str());
            }
        }
    for (const auto& sc : scalars)
        if (!std::isfinite(sc.second)) {
            std::ostringstream os;
            os.precision(9);
            os << tag << " : NaN/Inf detecte au pas " << step << " (t = " << t
               << " s) : grandeur globale " << sc.first << " = " << sc.second
               << " (aucune composante nodale non finie) ; reduisez dtFactor";
            throw NanError(os.str());
        }
}

} // namespace guards
} // namespace rockim
