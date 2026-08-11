#include "rockim/Config.hpp"
#include <fstream>
#include <stdexcept>

namespace rockim {

static std::string trim(const std::string& s) {
    const char* ws = " \t\r\n";
    auto a = s.find_first_not_of(ws);
    if (a == std::string::npos) return "";
    auto b = s.find_last_not_of(ws);
    return s.substr(a, b - a + 1);
}

Config Config::load(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Config: cannot open '" + path + "'");
    Config c;
    std::string line;
    while (std::getline(in, line)) {
        auto h = line.find('#');
        if (h != std::string::npos) line = line.substr(0, h);
        line = trim(line);
        if (line.empty()) continue;
        auto eq = line.find('=');
        if (eq == std::string::npos) continue;
        c.kv_[trim(line.substr(0, eq))] = trim(line.substr(eq + 1));
    }
    return c;
}

std::string Config::gets(const std::string& key, const std::string& def) const {
    auto it = kv_.find(key);
    return it == kv_.end() ? def : it->second;
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
    auto it = kv_.find(key);
    return it == kv_.end() ? def : parseNum(key, it->second);
}
int Config::geti(const std::string& key, int def) const {
    auto it = kv_.find(key);
    if (it == kv_.end()) return def;
    double d = parseNum(key, it->second);
    int i = (int)d;
    if ((double)i != d)
        throw std::runtime_error("Config: key '" + key + "' expects an "
                                 "integer, got '" + it->second + "'");
    return i;
}
bool Config::getb(const std::string& key, bool def) const {
    auto it = kv_.find(key);
    if (it == kv_.end()) return def;
    const std::string& v = it->second;
    return v == "1" || v == "true" || v == "yes" || v == "on";
}
std::string Config::reqs(const std::string& key) const {
    auto it = kv_.find(key);
    if (it == kv_.end()) throw std::runtime_error("Config: missing required key '" + key + "'");
    return it->second;
}
double Config::reqd(const std::string& key) const {
    return parseNum(key, reqs(key));
}

} // namespace rockim
