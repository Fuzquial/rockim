#pragma once
// ---------------------------------------------------------------------------
// ThermoBench — banc THERMODYNAMIQUE generique des lois de comportement
// (2026-09-06, etape 1 du chantier DFH+).
//
// POURQUOI. La relecture adverse de `CONTINUUM/loi_dfh_plus/loi_DFH_plus.pdf`
// a etabli que le CADRE de DP-DFH est SUR-DETERMINE (contrainte effective a la
// Lemaitre + equivalence en energie + energie libre ecrite independamment :
// trois enonces pour deux libertes) et que plusieurs dissipations peuvent y
// devenir negatives. La decision de Fernando (2026-09-06) est de reconstruire
// la loi sur une ENERGIE LIBRE POSTULEE UNIQUE, dont on derive la contrainte
// (sigma = d(rho psi)/d eps^e) et les forces motrices (Y_k = -d(rho psi)/d D_k).
// Un tel cadre se PROUVE au point materiel : symetrie majeure de la tangente,
// dissipation positive, reduction exacte au cas connu, objectivite, continuite.
//
// Ce banc est ecrit AVANT la loi, exactement comme la fiche
// « banc court avant run long » l'exige : il doit d'abord ECHOUER sur une loi
// dont on sait qu'elle viole quelque chose (`dpdfh`, plasticite NON ASSOCIEE
// + endommagement pilote par le temps) et PASSER sur `elastic`, dont la
// thermodynamique est exacte. Les deux mesures sont le controle du banc.
//
// USAGE
//   rockim thermobench <loi> [sortie.csv] [options]
//     --ref <loi>     compare la loi testee a une loi de REFERENCE sur les
//                     memes chemins (test 3 « reduction », 1e-12) ; sans
//                     cette option le test 3 verifie la reduction ELASTIQUE
//                     (amplitudes sous tous les seuils => sigma doit valoir
//                     lambda tr(eps) I + 2 G eps a 1e-12 pres)
//     --deck <cfg>    carte materiau / cles de loi supplementaires, ajoutees
//                     APRES la carte par defaut (Red Bohus) : derniere valeur
//                     gagnante, comme tout deck rockim
//     --draws N       nombre de tirages (defaut 12000, minimum utile 1e4)
//     --seed S        graine (defaut 20260906) — reproductibilite exacte
//     --max-rows N    nombre maximum de lignes ecrites dans le CSV (defaut
//                     20000) ; le COMPTE des violations, lui, est complet
//     --probe         force l'estimateur par SONDE DE DECHARGE du test 2 meme
//                     si la loi expose son energie libre (MatLaw::freeEnergy) :
//                     c'est ainsi qu'on compare une loi neuve a `dpdfh` avec le
//                     MEME instrument
//
// Code de retour : 0 si aucun test ne trouve de violation, 1 sinon, 2 si la
// loi ne se construit pas. Le banc N'ECRIT RIEN dans l'etat du depot et ne
// touche a AUCUNE loi existante : il n'utilise que l'interface publique
// `MatLaw::stress(eps, MatState&, dt, lc)`.
// ---------------------------------------------------------------------------
#include <string>

namespace rockim {

struct ThermoBenchOpts {
    std::string law;                       // loi testee (obligatoire)
    std::string refLaw;                    // loi de reference (test 3) ; vide
    std::string csv = "thermobench.csv";   // CSV des violations
    std::string deck;                      // deck materiau optionnel
    unsigned long long seed = 20260906ull;
    int nDraws = 12000;
    int maxRows = 20000;
    // --probe : FORCE l'estimateur par sonde de decharge meme si la loi expose
    // son energie libre. Sert a comparer une loi neuve a une loi ancienne AVEC
    // LE MEME INSTRUMENT (sans quoi les deux colonnes ne sont pas comparables).
    bool forceProbe = false;
};

// Joue le banc et imprime le compte rendu. Voir src/ThermoBench.cpp pour la
// definition PRECISE de chacun des cinq tests et de ses tolerances.
int thermoBench(const ThermoBenchOpts& o);

} // namespace rockim
