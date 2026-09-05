#pragma once
// ---------------------------------------------------------------------------
// Config: minimal "key = value" file reader ('#' starts a comment).
// Keeps the project dependency-free (no JSON library needed).
//
// C1 / C7 (w21, 2026-09-05, decision de Fernando 20:00) : Config SUIT les cles
// CONSOMMEES. Chaque getter (gets/getd/geti/getb/reqs/reqd/has/keysWithPrefix)
// marque la cle lue et note la valeur par defaut employee quand la cle est
// absente du deck. Les solveurs copient Config par valeur (`Config cfg_`) : le
// stockage est donc PARTAGE entre copies (shared_ptr), pour que les lectures
// faites par le solveur, Material::from, MatLaw... soient vues par le deck
// d'origine que main.cpp audite apres init() (KeyGuard.hpp). keys() (registre
// des modes) ne marque rien.
// w22 (relecture, cout) : apres l'audit, main.cpp appelle seal() : les getters
// ne suivent plus rien (ni verrou, ni chaine du defaut) — lecture directe de la
// table, comme avant w21 (confineGaugeTime est lu a CHAQUE pas dans step()).
// Aucun flottant du calcul n'est touche : bit-identite par construction.
// ---------------------------------------------------------------------------
#include <atomic>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace rockim {

class Config {
public:
    Config();
    static Config load(const std::string& path);

    bool has(const std::string& key) const;

    std::string gets(const std::string& key, const std::string& def) const;
    double      getd(const std::string& key, double def) const;
    int         geti(const std::string& key, int def) const;
    bool        getb(const std::string& key, bool def) const;

    // Required variants (throw with a clear message if missing)
    std::string reqs(const std::string& key) const;
    double      reqd(const std::string& key) const;

    // Cles presentes commencant par `prefix` (ordre lexicographique) — pour
    // qu'une loi refuse une cle mal orthographiee de sa famille au lieu de
    // l'ignorer en silence (revue cdp du 2026-09-04). Marque les cles rendues.
    std::vector<std::string> keysWithPrefix(const std::string& prefix) const;
    // Toutes les cles du deck (ordre lexicographique) — registre des cles par
    // mode (C6, w20). Lecture seule, NE marque PAS.
    std::vector<std::string> keys() const;

    // ---- C1 / C7 : suivi des cles consommees ------------------------------
    // Cles du deck jamais lues par aucun getter (ordre lexicographique).
    std::vector<std::string> unusedKeys() const;
    // Numero de ligne (1-based) de la derniere definition d'une cle du deck ;
    // 0 si absente.
    int lineOf(const std::string& key) const;
    // Cles consommees (deck ou defaut), ordre lexicographique.
    std::vector<std::string> consumedKeys() const;
    // Valeur EFFECTIVE d'une cle consommee : valeur du deck (fromDeck, line)
    // ou defaut du code (defaults : tous les defauts distincts rencontres,
    // dans l'ordre de lecture ; > 1 = deux lectures a defauts differents).
    struct Effective {
        std::string key;
        std::string value;               // deck, ou premier defaut lu
        bool fromDeck = false;
        int line = 0;
        std::vector<std::string> defaults;
    };
    std::vector<Effective> effective() const;   // tri stable par cle
    // Toutes les cles du DECK (ordre lexicographique) avec valeur, ligne et
    // etat consomme — pour config_effective.cfg (C7, w22 : les cles du deck non
    // consommees mais legitimes y restent actives, le fichier est rejouable).
    struct DeckKey {
        std::string key;
        std::string value;
        int line = 0;
        bool consumed = false;
    };
    std::vector<DeckKey> deck() const;
    // Fin du suivi (apres l'audit et config_effective) : les getters lisent la
    // table sans verrou ni marquage. Irreversible ; la table n'est plus modifiee.
    void seal() const;
    bool sealed() const;
    const std::string& path() const;
    std::size_t size() const;

private:
    struct Entry {
        std::string value;
        int line = 0;
        bool consumed = false;
    };
    struct Data {
        std::string path;
        std::map<std::string, Entry> kv;
        // cle -> defauts distincts employes (ordre d'apparition) ; ne contient
        // que des cles ABSENTES du deck lues avec un defaut
        std::map<std::string, std::vector<std::string>> defaults;
        std::mutex mu;                   // getters possiblement appeles hors init
        std::atomic<bool> sealed{false}; // w22 : suivi termine (lecture seule)
    };
    std::shared_ptr<Data> d_;

    // marque la cle ; rend l'entree du deck ou nullptr ; note le defaut
    const Entry* touch(const std::string& key, const std::string* def) const;
    // lecture seule (apres seal) : ni verrou ni marquage
    const Entry* peek(const std::string& key) const;
};

} // namespace rockim
