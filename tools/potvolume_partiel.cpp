// ---------------------------------------------------------------------------
// potvolume_partiel — force de contact de potForce = volume (Liu et al. 2022)
// contre le potentiel de Munjiza, sur deux tetras face contre face a
// enfoncement fixe, l'un decale pour ne couvrir qu'une fraction de la face.
//
// Lie au MEME en-tete que le solveur (rockim/PotentialContact.hpp). Facteur
// 8/3 = equivalence exacte face contre face entre tetras egaux. Resultat
// (2026-10-03) : le rapport volume/Munjiza suit la fraction de face couverte
// (1,01 a 97 %, 0,82 a 79 %, 0,49 a 47 %, 0,25 a 24 %), et 0,77 pour un tetra
// moitie plus petit entierement en contact : la force volume croit comme
// l'aire^2 de contact, celle de Munjiza comme l'aire. Voir DOCUMENTATION,
// ligne potForce.
//   build : g++ -std=c++17 -O2 -Iinclude -I<eigen> tools/potvolume_partiel.cpp -o potvolume_partiel
// ---------------------------------------------------------------------------
#include "rockim/PotentialContact.hpp"
#include <cmath>
#include <cstdio>
#include <initializer_list>
#include <utility>
using namespace rockim::pot3;

static double vol(const V3 p[4]) {
    return std::abs((p[1] - p[0]).dot((p[2] - p[0]).cross(p[3] - p[0]))) / 6;
}

int main() {
    const double h = 1, d = 0.02, p = 1;   // hauteur, enfoncement, penalite
    // A : base dans z = 0, sommet en +z ; B : base dans z = d, sommet en -z,
    // de taille relative sB, decale de s en x.
    for (double sB : {1.0, 0.5}) {
        std::printf("taille B / A = %.1f\n", sB);
        for (double s : {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6}) {
            V3 a[4] = {V3(0, 0, 0), V3(1, 0, 0), V3(0.5, 0.866, 0), V3(0.5, 0.289, h)};
            V3 b[4] = {V3(s, 0, d), V3(s + sB, 0, d), V3(s + 0.5 * sB, 0.866 * sB, d),
                       V3(s + 0.5 * sB, 0.289 * sB, d - h * sB)};
            if ((a[1] - a[0]).dot((a[2] - a[0]).cross(a[3] - a[0])) < 0) std::swap(a[1], a[2]);
            if ((b[1] - b[0]).dot((b[2] - b[0]).cross(b[3] - b[0])) < 0) std::swap(b[1], b[2]);
            PairForce3 Rm{}, Rv{};
            const bool om = pairForce(a, b, p, Rm);
            const double VA = vol(a), VB = vol(b), Vr = 2 * VA * VB / (VA + VB);
            const bool ov = pairForceVolume(a, b, (8.0 / 3) * p, Vr, Rv);
            std::printf("  decalage %.1f  recouvrement/face A %.3f  Fz Munjiza %+.5e"
                        "  Fz volume %+.5e  volume/Munjiza %.3f\n",
                        s, Rv.vol / (d * 0.433), om ? Rm.F.z() : 0, ov ? Rv.F.z() : 0,
                        (om && ov) ? Rv.F.z() / Rm.F.z() : 0);
        }
    }
}
