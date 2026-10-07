#pragma once
#include <cmath>
#include <string>

#include <Eigen/Dense>
// ---------------------------------------------------------------------------
// CUTTER PDC 3D — NOYAU GEOMETRIQUE (2026-09-03).
//
// Reproduction de Heilman, Euser, Frash, Meng, Li, Lei & Rougier, ARMA 24-0238
// (Los Alamos, code HOSS) : un cutter PDC UNIQUE, disque de 13 mm de diametre et
// 2,5 mm d'epaisseur, arete chanfreinee, corps rigide a vitesse imposee, sur un
// granite Utah FORGE 40 x 30 x 20 mm.
//
// CE QUE CE FICHIER EST : la geometrie SEULE — « ce noeud est-il dans le
// cutter, a quelle profondeur, et dans quelle direction le repousser ». Pas
// d'Eigen de solveur, pas de masse, pas de pas de temps : testable en forme
// fermee par `rockim selftest-pdc3d` (motif T0, ToolSignorini.hpp — le banc
// verifie CE noyau, pas une transcription).
//
// CE QUE CE FICHIER N'EST PAS : l'extrusion du coin 2D. En 2D (FdemSolver.cpp
// ~6216-6260) le cutter est un COIN infini en profondeur, teste par trois
// booleens sans rapport (dn, dt2, floorFlat) et dont la normale est GELEE a
// rakeNormal() quelle que soit la facette touchee (l. 6259). Ici le cutter est
// un CYLINDRE FINI CHANFREINE, convexe, et l'appartenance sort d'une DISTANCE
// SIGNEE exacte calculee dans le demi-plan meridien (rho, s) : la normale de
// repoussement est celle de la facette reellement la plus proche (face de
// coupe, face arriere, tranche ou chanfrein) et varie continument a l'exterieur
// du solide. C'est la surface de contact QUI S'ELARGIT avec l'enfoncement —
// l'argument meme de l'article pour expliquer que 20 deg de garde donne plus de
// force que 10 deg (« increased area of contact »).
//
// REPERE, convention de Tool.hpp:67-74 RELEVEE y -> z, pour qu'un meme signe de
// backRakeDeg soit la meme physique en 2D et en 3D :
//   e1 = direction de coupe (+x)  ;  e3 = normale a la surface libre (+z)
//   n = (cos b, 0, sin b)   normale SORTANTE de la face de coupe (vers la roche)
//   u = (-sin b, 0, cos b)  de l'arete vers le sommet de la face
//   w = (0, 1, 0)           lateral, dans le plan de la face
//   E = tool_.x             l'ARETE DE COUPE (comme en 2D), a la profondeur de passe
//   C = E + R u             centre du disque de coupe
//   M = C - (t/2) n         centre a mi-epaisseur
// E est le point le plus BAS du cercle de coupe : z(theta) = C_z + R cos(theta)
// cos b, minimal a theta = pi, ou il vaut E_z + R cos b - R cos b = E_z.
// Demonstration, pas hypothese (le selftest G3 la verifie numeriquement).
//
// LE SIGNE DE LA GARDE ARRIERE, A CONNAITRE. Pour b < 0 le sommet de la face
// est EN AVANT de l'arete et le corps (derriere la face, direction -n) monte
// en arriere : RIEN sous l'arete, le point le plus bas du solide est E. C'est
// l'orientation des decks 2D valides (cut2d_v3ter_rake.cfg, backRakeDeg = -20,
// « le corps monte en arriere, rien sous l arete ») et celle d'un vrai PDC de
// forage. Pour b > 0 le dos du cutter DESCEND de t.sin(b) sous la ligne
// d'arete et laboure son propre plancher (0,855 mm pour t = 2,5 mm a 20 deg) :
// c'est le defaut historique du deck v3 que `cutterFloor` corrige en 2D
// (FdemSolver.cpp:6235-6244). Le solveur 3D avertit quand b > 0.
//
// COORDONNEES MERIDIENNES. Pour un point p :
//   s   = (p - M).n              cote le long de l'axe du disque
//   rv  = (p - M) - s n          composante radiale, rho = |rv|, rh = rv/rho
// Le solide est l'intersection de demi-plans du demi-plan (rho >= 0, s), TOUS a
// gradient unitaire, donc chaque d_i est la vraie distance a sa facette :
//   d1 = s - t/2                      face de coupe     normale (0, +1)
//   d2 = -s - t/2                     face arriere      normale (0, -1)
//   d3 = rho - R                      tranche           normale (+1, 0)
//   d4 = (rho - (R - c)) sin g + (s - t/2) cos g
//                                     CHANFREIN         normale (sin g, cos g)
//        c = longueur du chanfrein mesuree sur la face de coupe depuis la
//        tranche, g = son angle depuis le plan de la face. c = 0 -> d4 est
//        redondante avec d1 et d3 et n'est pas evaluee : arete vive.
// Le polygone meridien est CONVEXE, donc a l'INTERIEUR la distance au bord vaut
// exactement -max(d_i), et la normale de repoussement est celle de la facette
// active : N = m_rho rh + m_s n. Exact, sans boucle ni tolerance.
// Degenerescence rho -> 0 (noeud sur l'axe) : rh est indefini, mais seules d1
// et d2 peuvent alors etre actives (d3 = -R, d4 << 0 des que t/2 < R), et la
// normale vaut +-n ; on pose rh = 0 et on refuse toute normale nulle.
//
// LA NORMALE N'EST PAS CONTINUE A L'INTERIEUR : sur l'axe median du polygone
// (ou deux d_i s'egalent) elle saute de 90 deg pour une arete vive, deux fois
// 45 deg avec un chanfrein a 45 deg. Le chanfrein DIVISE le pire saut par deux,
// il ne l'abolit pas. Consequence de conception, chiffree par le panel du
// 2026-09-03 : le noeud doit rester LOIN de l'axe median, penetration << min
// (t/2, R, c). Sous Signorini l'approche par pas vaut v.dt = 10 m/s x 4,4e-9 s
// = 0,044 um ; sous penalite l'ecretage geometrique vaut 0,6 h = 150 um sur la
// maille 0,25 mm de l'article, 3 400 fois plus, et atteint l'axe median dans
// l'anneau meme de l'arete de coupe. Le cutter PDC 3D tourne donc en
// `toolContact = signorini` ; le solveur avertit fort en penalite.
//
// LARGEUR DE CONTACT a la surface libre, en forme fermee (critere G2) : le
// cercle de coupe coupe le plan z = E_z + d (d = passe) en cos(theta) =
// d/(R cos b) - 1, d'ou une largeur 2 R sin(theta) =
//   2 sqrt(R^2 - (d / cos b - R)^2)
// qui redevient 2 sqrt(d (D - d)) a b = 0 — le 6,98 mm du deck 2D. Pour
// l'article (R = 6,5 mm, d = 1,016 mm, b = 20 deg) : 7,18 mm.
// ---------------------------------------------------------------------------
namespace rockim {
namespace pdc3d {

struct Geom {
    double R = 0.0065;      // rayon du disque [m]           (cutterDia / 2)
    double thick = 0.0025;  // epaisseur [m]                  (cutterThick)
    double cham = 0.0;      // longueur du chanfrein [m]      (chamferLen, 0 = vif)
    double chamDeg = 45.0;  // angle du chanfrein depuis la face [deg]
    double rakeDeg = -20.0; // garde arriere, convention Tool.hpp relevee [deg]
};

struct Frame {
    Eigen::Vector3d n, u, w;
};

inline Frame frame(double rakeDeg) {
    double b = rakeDeg * M_PI / 180.0;
    Frame f;
    f.n = {std::cos(b), 0.0, std::sin(b)};
    f.u = {-std::sin(b), 0.0, std::cos(b)};
    f.w = {0.0, 1.0, 0.0};
    return f;
}

struct Query {
    bool inside = false;
    double pen = 0.0;                          // >= 0 si inside
    Eigen::Vector3d normal{0.0, 0.0, 0.0};     // repousse le noeud HORS du cutter
    int facet = 0;                             // 1 face, 2 dos, 3 tranche, 4 chanfrein
};

// Distance signee au cutter dont l'arete de coupe est en `edge`.
inline Query query(const Geom& g, const Frame& f, const Eigen::Vector3d& edge,
                   const Eigen::Vector3d& p) {
    Query q;
    const double ht = 0.5 * g.thick;
    const Eigen::Vector3d C = edge + g.R * f.u;
    const Eigen::Vector3d M = C - ht * f.n;
    const Eigen::Vector3d d = p - M;
    const double s = d.dot(f.n);
    const Eigen::Vector3d rv = d - s * f.n;
    const double rho = rv.norm();
    const Eigen::Vector3d rh = rho > 1e-14 ? Eigen::Vector3d(rv / rho)
                                           : Eigen::Vector3d::Zero();
    double di[4];
    double mr[4], ms[4];
    di[0] = s - ht;        mr[0] = 0.0;  ms[0] = 1.0;
    di[1] = -s - ht;       mr[1] = 0.0;  ms[1] = -1.0;
    di[2] = rho - g.R;     mr[2] = 1.0;  ms[2] = 0.0;
    int nc = 3;
    if (g.cham > 0.0) {
        double ga = g.chamDeg * M_PI / 180.0;
        double sg = std::sin(ga), cg = std::cos(ga);
        di[3] = (rho - (g.R - g.cham)) * sg + (s - ht) * cg;
        mr[3] = sg; ms[3] = cg;
        nc = 4;
    }
    int k = 0;
    for (int i = 1; i < nc; ++i) if (di[i] > di[k]) k = i;
    if (di[k] >= 0.0) return q;                 // dehors ou sur la peau
    q.pen = -di[k];
    q.normal = mr[k] * rh + ms[k] * f.n;
    double nn = q.normal.norm();
    if (nn < 1e-12) return q;                   // axe + facette radiale : impossible, refuse
    q.normal /= nn;
    q.inside = true;
    q.facet = k + 1;
    return q;
}

// Largeur de contact a la surface libre, en forme fermee (voir en-tete).
inline double contactWidth(double R, double depth, double rakeDeg) {
    double cb = std::cos(rakeDeg * M_PI / 180.0);
    double a = depth / cb - R;
    double v = R * R - a * a;
    return v > 0.0 ? 2.0 * std::sqrt(v) : 0.0;
}

} // namespace pdc3d
} // namespace rockim

// Banc en forme fermee (`rockim selftest-pdc3d [csv]`), implemente dans
// src/ToolPdc3d.cpp. 0 = PASS, 1 = FAIL. Voir l'en-tete du .cpp.
int pdc3dSelftest(const std::string& csvPath);
