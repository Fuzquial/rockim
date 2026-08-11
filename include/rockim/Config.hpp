#pragma once
// ---------------------------------------------------------------------------
// Config: minimal "key = value" file reader ('#' starts a comment).
// Keeps the project dependency-free (no JSON library needed).
// ---------------------------------------------------------------------------
#include <map>
#include <string>

namespace rockim {

class Config {
public:
    static Config load(const std::string& path);

    bool has(const std::string& key) const { return kv_.count(key) > 0; }

    std::string gets(const std::string& key, const std::string& def) const;
    double      getd(const std::string& key, double def) const;
    int         geti(const std::string& key, int def) const;
    bool        getb(const std::string& key, bool def) const;

    // Required variants (throw with a clear message if missing)
    std::string reqs(const std::string& key) const;
    double      reqd(const std::string& key) const;

private:
    std::map<std::string, std::string> kv_;
};

} // namespace rockim
