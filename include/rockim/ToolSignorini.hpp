#pragma once
#include <string>
// ---------------------------------------------------------------------------
// A1/T0 — NOYAU D'IMPULSION DU CONTACT OUTIL DE SIGNORINI (CD-Lagrange).
//
// EXTRACTION PURE de la branche `toolContact = signorini` de
// FdemSolver::toolContact(). Motif : le banc T0 doit verifier CE noyau, et
// non une transcription — la lecon du 2026-08-18 (« un reglage que l'on peut
// ecrire, que le solveur accepte, et qui ne fait rien, est pire que son
// absence ») vaut aussi pour un test qui verifierait sa propre copie.
//
// CE QUI RESTE DANS LE SOLVEUR : la GEOMETRIE (normale, direction tangente,
// penetration selon la forme PDC / FLAT / DISC). Elle y est dupliquee A
// DESSEIN entre les voies penalite et Signorini, pour que la voie penalite
// reste bit-identique ; ce fichier ne la touche pas.
//
// CE QUI VOYAGE ICI : l'algebre impulsion / saut de vitesse, seule, sans
// Eigen ni etat de solveur — donc testable en forme fermee.
//
// SOURCES. Fekak, Brun, Gravouil et al. (2017), schema CD-Lagrange ;
// Dureisseix, Greco, Brun et al., « Explicit dynamics and non-smooth
// interface behaviors », JTCAM (2024), Algorithme 1 ; Ghesquiere-Dierickx,
// Anciaux, Acary & Molinari, arXiv:2606.01355 (Nonsmooth Newmark), qui
// applique la meme idee a la fragmentation a modele cohesif.
//
// LE THEOREME QUE LE BANC T0 VERIFIE. La masse de rockim est diagonale et
// l'outil est un obstacle RIGIDE a vitesse imposee : l'operateur de Delassus
// H = L M^-1 L^T est diagonal et vaut 1/m_i, donc l'impulsion se calcule en
// FORME FERMEE noeud par noeud, sans systeme a resoudre. La relation
// imposee est le lemme de viabilite de Moreau :
//
//     si g > 0  alors  r = 0  ;  sinon  0 <= v_n  ⊥  r_n >= 0
//
// AVEC relax = 0 (le defaut) la condition est v_n = 0 apres impulsion :
// contact normal PARFAITEMENT INELASTIQUE (restitution e = 0). Consequence
// mesurable, et c'est le critere du banc : un noeud INITIALEMENT AU REPOS
// heurte par l'outil repart exactement a v_outil selon la normale — jamais
// a 2 v_outil (borne du choc elastique), et a fortiori jamais aux 377 m/s
// PAR PAS que la voie penalite autorisait (mesure du 2026-08-18 : penalite
// kp = 4,83e10 N/m, ecretage geometrique a 0,6 h, masse nodale 2,46e-5
// kg/m, dt = 1,27e-9 s).
// ---------------------------------------------------------------------------

namespace rockim {
namespace toolsig {

// Resultat de l'impulsion d'UN noeud contre l'outil, dans le repere local
// (normale sortante de l'outil vers le noeud, tangente de glissement).
struct Impulse {
    double rn = 0.0;        // impulsion normale   [N.s/m en 2D]
    double rt = 0.0;        // impulsion tangentielle
    bool   active = false;  // false = condition de Signorini : aucune impulsion
    bool   sticking = false;// true = le cap de Coulomb n'a pas morde (collage)
};

// pen     : penetration, > 0 quand le noeud est DANS l'outil
// vnFree  : vitesse LIBRE relative a l'outil, projetee sur la normale
//           (< 0 = le noeud s'enfonce). « Libre » = apres application de
//           toutes les autres forces du pas : v* = v + (dt/m) f.
// vtFree  : idem, projetee sur la tangente
// m, dt   : masse nodale lumpee, pas de temps
// mu      : coefficient de frottement (cap de Coulomb sur l'IMPULSION)
// relax   : rattrapage facon Baumgarte de la penetration residuelle,
//           dans [0, 1]. 0 = condition de VITESSE pure (defaut) : on annule
//           l'approche, on ne resorbe pas la penetration deja acquise.
inline Impulse impulse(double pen, double vnFree, double vtFree,
                       double m, double dt, double mu, double relax) {
    Impulse r;
    // (3) le noeud se separe-t-il de lui-meme au pas suivant ? Condition de
    //     Signorini : gap predit g+ = -pen + dt v_n >= 0 -> impulsion NULLE.
    if (-pen + dt * vnFree >= 0.0) return r;
    // (4) impulsion normale. Positive par construction : jamais d'adhesion,
    //     l'outil ne peut que repousser.
    double vTarget = relax * pen / dt;
    double rn = m * (vTarget - vnFree);
    if (rn <= 0.0) return r;
    // (5) frottement : impulsion de COLLAGE, ecretee par le cap de Coulomb.
    //     Pas de regularisation en tanh — le cap porte sur l'impulsion, il
    //     n'a pas besoin d'etre lisse.
    double rt = -m * vtFree;
    double cap = mu * rn;
    r.sticking = (rt <= cap && rt >= -cap);
    if (rt > cap) rt = cap;
    else if (rt < -cap) rt = -cap;
    r.rn = rn;
    r.rt = rt;
    r.active = true;
    return r;
}

}  // namespace toolsig

// T0 : banc du noyau ci-dessus (selftest-toolcontact). Defini dans
// FdemSolver.cpp, aux cotes des autres selftests de la maison.
int toolSignoriniSelftest(const std::string& csvPath);

}  // namespace rockim
