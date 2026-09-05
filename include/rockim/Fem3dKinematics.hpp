#pragma once
// ---------------------------------------------------------------------------
// Fem3dKinematics — cinematique de HENCKY (deformation logarithmique) pour
// fem3d, cle opt-in `kinematics = hencky` (2026-09-05, rockim_f2w16.exe).
//
// Chemin historique (`kinematics = biot`, defaut) : F = R U par iteration de
// Newton sur R, deformation de BIOT eps = sym(R^T F) - I, la loi rend une
// contrainte sigma_c dans le repere co-rotationnel, forces nodales
// f_a = -V0 (R sigma_c) dN_a/dX : le couple conjugue est (contrainte de Biot,
// deformation de Biot), volume et gradients de REFERENCE.
//
// Chemin hencky : decomposition SPECTRALE de C = F^T F = sum lambda_i^2 n_i n_i^T
// (Eigen::SelfAdjointEigenSolver 3x3, exacte au conditionnement pres) :
//   U       = sum lambda_i n_i n_i^T          (etirement droit)
//   U^-1    = sum (1/lambda_i) n_i n_i^T
//   R       = F U^-1                           (polaire exacte, pas d'iteration)
//   eps     = ln U = sum ln(lambda_i) n_i n_i^T (deformation logarithmique)
// La loi recoit eps = ln U et rend sigma_c interpretee comme contrainte de
// CAUCHY co-rotationnelle ; premier Piola-Kirchhoff P = J sigma F^-T avec
// sigma = R sigma_c R^T et F^-T = R U^-1, soit P = J R sigma_c U^-1, et
// f_a = -V0 P dN_a/dX (equivalent a -V (R sigma_c R^T) dN_a/dx sur la
// configuration COURANTE, V = J V0, dN/dx = F^-T dN/dX ; on garde les
// gradients de reference deja stockes : une decomposition + trois produits
// 3x3 par element, rien de plus).
// Garde-fou : det F <= 1e-9 ou lambda_min^2 <= 1e-18 (element degenere, ou
// echec du solveur propre) -> retour false, l'appelant reprend le chemin
// historique (R = I, eps de Biot) exactement comme aujourd'hui.
// Ce fichier n'est inclus que par Fem3dSolver.cpp (branche hencky) et par le
// test unitaire etude_lois_fem/bitid_w16/selftest_hencky/test_kinematics.cpp.
// ---------------------------------------------------------------------------
#include <cmath>

#include <Eigen/Dense>

namespace rockim {
namespace fem3dkin {

// Decomposition spectrale de C = F^T F ; rend R (polaire exacte), eps = ln U
// et Uinv = U^-1. false si l'element est degenere (l'appelant garde son
// chemin historique). Aucun etat partage : sur en OpenMP.
inline bool hencky(const Eigen::Matrix3d& F, double det, Eigen::Matrix3d& R,
                   Eigen::Matrix3d& eps, Eigen::Matrix3d& Uinv) {
    if (!(det > 1e-9)) return false;
    const Eigen::Matrix3d C = F.transpose() * F;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(C);
    if (es.info() != Eigen::Success) return false;
    const Eigen::Vector3d l2 = es.eigenvalues();        // croissantes, = lambda_i^2
    if (!(l2(0) > 1e-18)) return false;                  // lambda_min <= 1e-9 : degenere
    const Eigen::Matrix3d& N = es.eigenvectors();
    Eigen::Vector3d lam, lnl, inv;
    for (int i = 0; i < 3; ++i) {
        lam(i) = std::sqrt(l2(i));
        lnl(i) = std::log(lam(i));
        inv(i) = 1.0 / lam(i);
    }
    eps  = N * lnl.asDiagonal() * N.transpose();
    Uinv = N * inv.asDiagonal() * N.transpose();
    R    = F * Uinv;
    return true;
}

} // namespace fem3dkin
} // namespace rockim
