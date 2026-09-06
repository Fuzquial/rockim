#pragma once
// ---------------------------------------------------------------------------
// MatLawDfhPlus — fabrique de la loi `law = dfhplus` (2026-09-06, etape 1 du
// chantier DFH+). La loi elle-meme est dans src/MatLawDfhPlus.cpp, dont
// l'en-tete porte l'ENERGIE LIBRE POSTULEE et toutes les derivees.
//
// Pourquoi un fichier separe : la regle absolue du chantier est que
// `DpDfhLaw` (src/MatLaw.cpp) reste INTACTE — code, cles et resultats. Une
// unite de compilation propre rend la separation verifiable a la lecture.
// ---------------------------------------------------------------------------
#include <memory>

#include "rockim/Config.hpp"
#include "rockim/MatLaw.hpp"
#include "rockim/Material.hpp"

namespace rockim {

// Construit la loi dfhplus a partir de la carte materiau partagee et des cles
// du deck (les memes dfh* que dpdfh, plus dfhpPsiClamp et dfhpVolInteg).
// Leve std::runtime_error sur une carte inadmissible.
std::unique_ptr<MatLaw> makeDfhPlusLaw(const Material& m, const Config& c);

} // namespace rockim
