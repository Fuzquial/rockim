#pragma once
// ---------------------------------------------------------------------------
// YangDif.hpp — facteurs d amplification dynamique de Yang, Xiang, Naderi,
// Wang, Aising, Ugarte & Latham (IJRMMS 191, 2025, 106125), leurs eq. 2 et 3,
// et mesure du taux principal d un tenseur taux 3x3. PARTAGES entre
// FdemSolver (2D) et Fdem3dSolver (3D).
//
// Ces deux fonctions vivaient dans le namespace anonyme de FdemSolver.cpp,
// donc invisibles du 3D. Elles sont promues ici SANS qu une seule expression
// change : il ne doit exister qu UNE transcription de l article dans le
// depot. Une divergence 2D/3D sur les bornes en dur (5e-6, 1e2, 1e4, 1,85,
// 1,84) serait MUETTE et fausserait toute comparaison dimensionnelle — le
// mode de panne exact que la garde de parite du 2026-08-18 cherchait a eviter.
//
// L exposant de la TRACTION est parametre parce que le 0,07 imprime est tres
// probablement une coquille recopiee de l eq. 2 :
//   * avec 0,07 la loi saute de 1 a 1,1245 en 5e-6 /s (+12,4 %) et de 1,5160
//     a 1,85 en 1e2 /s (+22 %) — elle ne se raccorde a aucune de ses bornes ;
//   * la COURBE TRACEE dans leur propre figure 2(b) suit un exposant voisin
//     de 0,17, et 0,1707 raccorde EXACTEMENT les deux bornes (1,0010 en bas,
//     1,8500 en haut). Deux raccords simultanes avec un seul parametre : ce
//     n est pas une coincidence.
//
// CORRECTION 2026-08-30 (contre-audit). Cette ligne portait « 1,0031 en bas ».
// C etait FAUX : 0,95 + 0,41*(5e-6)^0,1707 = 1,001039, soit 1,0010 — valeur
// que le lot 2b (biblio_insertion/2026-08-29_lot2b_...md l. 217) donnait deja
// juste. Le raccord HAUT, lui, est exact : 0,95 + 0,41*(1e2)^0,1707 = 1,84989.
// L argument est inchange, mais il ne faut pas citer 1,0031 : c est le seul
// chiffre errone de cet en-tete, et il etait justement produit en PREUVE.
//
// ET LE MOT « PREDICTION » N EST PAS SOUTENABLE. L article confirmateur est
// IJRMMS 206 (2026) 106660, soit la MEME ANNEE que cette derivation
// (2026-08-18) ; rien ici n etablit qu il n etait pas deja paru. A presenter
// comme une INFERENCE INDEPENDANTE — ce qu elle est, et ce qui suffit :
// l exposant litteral 0,07 ne raccorde AUCUNE de ses deux bornes (1,1245 au
// lieu de 1 ; 1,5160 au lieu de 1,85, soit 22 % de discontinuite), 0,1707
// raccorde les deux, et l article de 2026 imprime 0,17.
// L eq. 2 (compression) est elle CONFIRMEE a 0,07 par sa figure 2(a).
//
// MESURE : le saut de 22 % en 1e2 /s est un ATTRACTEUR dans un schema
// d insertion extrinseque — un joint qui franchit le seuil voit sa resistance
// bondir et cesse de s inserer, si bien que la population inseree s empile
// juste sous 1e2. Mesure du 2026-08-18 sur le banc de traction, meme config a
// l exposant pres : mediane 99,36 /s (max 99,9988) avec l exposant litteral,
// contre 40,22 /s avec 0,1707. Verrouille par dif_yang_litteral_2d et
// dif_yang_fig2_2d dans verify_suite.py.
//
// L ecretage a [1, plateau] est un no-op sur la transcription litterale
// (toutes ses valeurs y sont deja) : il ne sert qu a garantir la monotonie et
// la borne pour la variante figure.
// ---------------------------------------------------------------------------
#include <algorithm>
#include <cmath>

#include <Eigen/Dense>

namespace rockim {

inline double difTensionYang(double edot, double n) {
    if (edot <= 5.0e-6) return 1.0;
    if (edot > 1.0e2)   return 1.85;
    double v = 0.95 + 0.41 * std::pow(edot, n);
    return v < 1.0 ? 1.0 : (v > 1.85 ? 1.85 : v);
}

inline double difCompressionYang(double edot) {
    if (edot <= 5.0e-6) return 1.0;
    if (edot > 1.0e4)   return 1.84;
    double v = 0.77 + 0.56 * std::pow(edot, 0.07);
    return v < 1.0 ? 1.0 : (v > 1.84 ? 1.84 : v);
}

// ---------------------------------------------------------------------------
// max |lambda| d une matrice SYMETRIQUE 3x3, forme fermee trigonometrique
// (Smith 1961). C est le remplacant 3D de la formule de Mohr 2x2 du solveur
// 2D : celle-ci est un cercle de Mohr ecrit a la main et n a AUCUNE
// transposition 3x3.
//
// POURQUOI PAS Eigen ICI. Le chemin VTU du solveur 3D utilise bien
// SelfAdjointEigenSolver, mais il ne le paie qu aux ecritures de frame ; ici
// on est appele PAR ELEMENT ET PAR PAS, soit ~1e11 fois sur un banc de
// percussion. La forme fermee coute environ 4x moins que computeDirect et
// 4,5x moins que le solveur par defaut, sur un coeur d elementForces qui
// tourne autour de 300 ns.
//
// Les valeurs propres d une symetrique sont reelles et se trient
// e1 >= e2 >= e3, donc max|lambda| = max(|e1|, |e3|) : la mediane ne peut
// pas les depasser.
// ---------------------------------------------------------------------------
inline double maxAbsEigSym3(const Eigen::Matrix3d& A) {
    const double p1 = A(0, 1) * A(0, 1) + A(0, 2) * A(0, 2)
                    + A(1, 2) * A(1, 2);
    const double q = A.trace() / 3.0;
    if (!(p1 > 0.0))                       // deja diagonale (inclut A = 0)
        return std::max(std::abs(A(0, 0)),
                        std::max(std::abs(A(1, 1)), std::abs(A(2, 2))));
    const double d0 = A(0, 0) - q, d1 = A(1, 1) - q, d2 = A(2, 2) - q;
    const double p2 = d0 * d0 + d1 * d1 + d2 * d2 + 2.0 * p1;
    const double p = std::sqrt(p2 / 6.0);
    if (!(p > 0.0)) return std::abs(q);    // spectre degenere : triple q
    const Eigen::Matrix3d B = (A - q * Eigen::Matrix3d::Identity()) / p;
    double r = 0.5 * B.determinant();
    r = r < -1.0 ? -1.0 : (r > 1.0 ? 1.0 : r);   // acos hors [-1,1] = NaN
    const double phi = std::acos(r) / 3.0;
    const double e1 = q + 2.0 * p * std::cos(phi);
    const double e3 = q + 2.0 * p * std::cos(phi + 2.0943951023931953);
    return std::max(std::abs(e1), std::abs(e3));
}


// ---------------------------------------------------------------------------
// Enveloppe de Mohr-Coulomb : la part FROTTEMENT, a multiplier par tan(phi).
//
//   YAN (son eq. 8)  : max(0, -sig)      -- le frottement disparait des que le
//                      joint est en traction, la cohesion seule tient.
//   YANG (son eq. 1) : -min(sig, ft)     -- le frottement DECROIT lineairement
//                      en traction jusqu au cut-off, puis reste constant :
//                        fs = c - sig tanphi   si sig < ft
//                        fs = c - ft  tanphi   si sig > ft
//
// Convention rockim : sig > 0 en TRACTION. Les deux formes coincident donc
// exactement en compression ; elles ne different qu en traction, ou celle de
// Yang AFFAIBLIT le cisaillement. Chiffre sur le banc de percussion
// (c = 217 MPa, ft = 87 MPa, phi = 40 deg) : au cut-off, 144 MPa contre 217,
// soit -34 %. Ce n est pas un detail — c est ce qui gouverne le partage entre
// rupture en traction et rupture en cisaillement dans les zones tendues, donc
// le facies radial.
//
// Le cut-off est pris sur le ft NON ENDOMMAGE, par coherence avec le
// glissement au pic s_p du solveur, lui aussi evalue sur l enveloppe non
// endommagee (regle maison : ce qui fixe une resistance ne depend pas de D).
// Yang ne tranche pas ce point.
// ---------------------------------------------------------------------------
inline double mcFrictionTerm(double sig, double ft, bool yang) {
    return yang ? -std::min(sig, ft) : std::max(0.0, -sig);
}

} // namespace rockim
