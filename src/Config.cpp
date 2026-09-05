#include "rockim/Config.hpp"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace rockim {

static std::string trim(const std::string& s) {
    const char* ws = " \t\r\n";
    auto a = s.find_first_not_of(ws);
    if (a == std::string::npos) return "";
    auto b = s.find_last_not_of(ws);
    return s.substr(a, b - a + 1);
}

Config::Config() : d_(std::make_shared<Data>()) {}

Config Config::load(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Config: cannot open '" + path + "'");
    Config c;
    c.d_->path = path;
    std::string line;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        auto h = line.find('#');
        if (h != std::string::npos) line = line.substr(0, h);
        line = trim(line);
        if (line.empty()) continue;
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        // cle repetee = la derniere gagne (overrides de verify_suite) ; la
        // ligne retenue est celle de la derniere definition
        Entry& e = c.d_->kv[trim(line.substr(0, eq))];
        e.value = trim(line.substr(eq + 1));
        e.line = lineNo;
    }
    return c;
}

// Marque `key` consommee. Si la cle est absente du deck et qu'un defaut est
// fourni, le note (C7 : config_effective.cfg). Rend l'entree du deck ou nullptr.
const Config::Entry* Config::touch(const std::string& key, const std::string* def) const {
    std::lock_guard<std::mutex> lock(d_->mu);
    auto it = d_->kv.find(key);
    if (it != d_->kv.end()) {
        it->second.consumed = true;
        return &it->second;
    }
    if (def) {
        auto& v = d_->defaults[key];
        bool seen = false;
        for (const auto& s : v) if (s == *def) { seen = true; break; }
        if (!seen) v.push_back(*def);
    }
    return nullptr;
}

static std::string fmtNum(double d) {
    std::ostringstream os;
    os.precision(15);
    os << d;
    return os.str();
}

// w22 : apres seal(), la table n'est plus modifiee — lecture concurrente sure
// sans verrou (std::map non mute), aucune chaine de defaut construite.
const Config::Entry* Config::peek(const std::string& key) const {
    auto it = d_->kv.find(key);
    return it == d_->kv.end() ? nullptr : &it->second;
}
void Config::seal() const { d_->sealed.store(true); }
bool Config::sealed() const { return d_->sealed.load(); }

bool Config::has(const std::string& key) const {
    if (d_->sealed.load()) return peek(key) != nullptr;
    return touch(key, nullptr) != nullptr;
}

std::string Config::gets(const std::string& key, const std::string& def) const {
    if (d_->sealed.load()) { const Entry* e = peek(key); return e ? e->value : def; }
    const Entry* e = touch(key, &def);
    return e ? e->value : def;
}
// Strict numeric parsing: bare std::stod accepts any parsable PREFIX and
// silently drops the rest, so a French-locale "0,5" reads as 0.0 and a typo
// "5Oe9" as 5.0 — wrong physics with no diagnostic. Require the whole value
// to be consumed and name the offending key in the error.
static double parseNum(const std::string& key, const std::string& val) {
    std::size_t used = 0;
    double d;
    try { d = std::stod(val, &used); }
    catch (const std::exception&) {
        throw std::runtime_error("Config: key '" + key + "' has non-numeric "
                                 "value '" + val + "'");
    }
    while (used < val.size()
           && (val[used] == ' ' || val[used] == '\t')) ++used;
    if (used != val.size())
        throw std::runtime_error("Config: key '" + key + "' has trailing "
                                 "garbage in value '" + val + "' (decimal "
                                 "COMMA instead of point?)");
    return d;
}

double Config::getd(const std::string& key, double def) const {
    if (d_->sealed.load()) { const Entry* e = peek(key); return e ? parseNum(key, e->value) : def; }
    const std::string sdef = fmtNum(def);
    const Entry* e = touch(key, &sdef);
    return e ? parseNum(key, e->value) : def;
}
int Config::geti(const std::string& key, int def) const {
    const Entry* e;
    if (d_->sealed.load()) e = peek(key);
    else { const std::string sdef = std::to_string(def); e = touch(key, &sdef); }
    if (!e) return def;
    double d = parseNum(key, e->value);
    int i = (int)d;
    if ((double)i != d)
        throw std::runtime_error("Config: key '" + key + "' expects an "
                                 "integer, got '" + e->value + "'");
    return i;
}
bool Config::getb(const std::string& key, bool def) const {
    const Entry* e;
    if (d_->sealed.load()) e = peek(key);
    else { const std::string sdef = def ? "true" : "false"; e = touch(key, &sdef); }
    if (!e) return def;
    const std::string& v = e->value;
    return v == "1" || v == "true" || v == "yes" || v == "on";
}
std::string Config::reqs(const std::string& key) const {
    const Entry* e = d_->sealed.load() ? peek(key) : touch(key, nullptr);
    if (!e) throw std::runtime_error("Config: missing required key '" + key + "'");
    return e->value;
}
double Config::reqd(const std::string& key) const {
    return parseNum(key, reqs(key));
}
std::vector<std::string> Config::keys() const {
    std::lock_guard<std::mutex> lock(d_->mu);
    std::vector<std::string> out;
    out.reserve(d_->kv.size());
    for (const auto& kv : d_->kv) out.push_back(kv.first);
    return out;
}
std::vector<std::string> Config::keysWithPrefix(const std::string& prefix) const {
    std::lock_guard<std::mutex> lock(d_->mu);
    std::vector<std::string> out;
    const bool mark = !d_->sealed.load();
    for (auto& kv : d_->kv)
        if (kv.first.compare(0, prefix.size(), prefix) == 0) {
            if (mark) kv.second.consumed = true;
            out.push_back(kv.first);
        }
    return out;
}

// ---- C1 / C7 ---------------------------------------------------------------
std::vector<std::string> Config::unusedKeys() const {
    std::lock_guard<std::mutex> lock(d_->mu);
    std::vector<std::string> out;
    for (const auto& kv : d_->kv)
        if (!kv.second.consumed) out.push_back(kv.first);
    return out;
}
int Config::lineOf(const std::string& key) const {
    std::lock_guard<std::mutex> lock(d_->mu);
    auto it = d_->kv.find(key);
    return it == d_->kv.end() ? 0 : it->second.line;
}
std::vector<std::string> Config::consumedKeys() const {
    std::vector<std::string> out;
    for (const auto& e : effective()) out.push_back(e.key);
    return out;
}
std::vector<Config::Effective> Config::effective() const {
    std::lock_guard<std::mutex> lock(d_->mu);
    std::vector<Effective> out;
    // std::map : deja trie par cle, tri stable par construction
    for (const auto& kv : d_->kv)
        if (kv.second.consumed) {
            Effective e;
            e.key = kv.first;
            e.value = kv.second.value;
            e.fromDeck = true;
            e.line = kv.second.line;
            out.push_back(e);
        }
    for (const auto& kd : d_->defaults) {
        if (d_->kv.count(kd.first)) continue;      // lue du deck : deja listee
        Effective e;
        e.key = kd.first;
        e.defaults = kd.second;
        e.value = kd.second.empty() ? "" : kd.second.front();
        out.push_back(e);
    }
    std::stable_sort(out.begin(), out.end(),
                     [](const Effective& a, const Effective& b) { return a.key < b.key; });
    return out;
}
std::vector<Config::DeckKey> Config::deck() const {
    std::lock_guard<std::mutex> lock(d_->mu);
    std::vector<DeckKey> out;
    out.reserve(d_->kv.size());
    for (const auto& kv : d_->kv) {
        DeckKey k;
        k.key = kv.first;
        k.value = kv.second.value;
        k.line = kv.second.line;
        k.consumed = kv.second.consumed;
        out.push_back(k);
    }
    return out;
}
const std::string& Config::path() const { return d_->path; }
std::size_t Config::size() const {
    std::lock_guard<std::mutex> lock(d_->mu);
    return d_->kv.size();
}

} // namespace rockim
