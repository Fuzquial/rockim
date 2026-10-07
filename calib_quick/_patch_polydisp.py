import io, os
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
H="include/rockim/Tessellation.hpp"; C="src/Tessellation.cpp"; F="src/FdemSolver.cpp"
h=io.open(H,encoding="utf-8").read(); c=io.open(C,encoding="utf-8").read(); f=io.open(F,encoding="utf-8").read()
def rep(s,a,b,name):
    assert s.count(a)==1,(name,s.count(a)); return s.replace(a,b,1)
# ---- hpp : signature ----
h=rep(h,"""                              std::mt19937& rng, bool randomSeeds = false,
                              bool useDelaunay = false, double elemSize = 0.0);""",
"""                              std::mt19937& rng, bool randomSeeds = false,
                              bool useDelaunay = false, double elemSize = 0.0,
                              double sizeSpread = 0.0,
                              const std::vector<double>& phaseSize = {});
    // sizeSpread : POLYDISPERSITE (2026-09-02, opt-in). Ecart-type de
    //            ln(taille) : chaque graine de Poisson porte son propre
    //            espacement s_i = s L_i, L_i log-normale de moyenne 1, et deux
    //            graines s acceptent a distance >= 0,35 (s_i + s_j). Moyenne
    //            conservee, cellules petites et grandes melangees. Exige
    //            randomSeeds. 0 = comportement historique, bit-identique.
    // phaseSize : taille cible PAR PHASE [m] (<= 0 = indifferent). Les grains
    //            sont alors affectes par ordre d aire decroissante avec une
    //            affinite log-normale a la taille de la phase, le deficit
    //            restant pondere par l AIRE : les fractions globales sont
    //            conservees par construction, quelle que soit la distribution.""","hpp sig")
# ---- cpp : signature de la definition ----
c=rep(c,"                                 bool useDelaunay, double elemSize) {",
"                                 bool useDelaunay, double elemSize,\n                                 double sizeSpread,\n                                 const std::vector<double>& phaseSize) {","cpp sig")
# ---- cpp : tirage des graines polydisperse ----
c=rep(c,"""        double dmin = 0.7 * s;
        long target = std::lround(W * H / (std::sqrt(3.0) / 2.0 * s * s));
        double cs = dmin;
""","""        double dmin = 0.7 * s;
        long target = std::lround(W * H / (std::sqrt(3.0) / 2.0 * s * s));
        // Polydispersite : espacement PROPRE a chaque graine. Le tirage de la
        // loi normale n est consomme que si sizeSpread > 0 : sans la cle, la
        // sequence aleatoire et l arithmetique sont strictement inchangees.
        std::vector<double> sp;
        std::normal_distribution<double> N01(0.0, 1.0);
        auto drawSp = [&]() {
            if (sizeSpread <= 0.0) return s;
            return s * std::exp(sizeSpread * N01(rng) - 0.5 * sizeSpread * sizeSpread);
        };
        double cs = (sizeSpread > 0.0) ? dmin * std::exp(2.0 * sizeSpread) : dmin;
""","cpp seeds a")
c=rep(c,"""            Vec2 p(Ux(rng), Uy(rng));
            int ci = std::clamp((int)(p.x() / W * gx), 0, gx - 1);
            int cj = std::clamp((int)(p.y() / H * gy), 0, gy - 1);
            bool ok = true;
            for (int dj = -1; dj <= 1 && ok; ++dj)
                for (int di = -1; di <= 1 && ok; ++di) {
                    int cx = ci + di, cy = cj + dj;
                    if (cx < 0 || cy < 0 || cx >= gx || cy >= gy) continue;
                    for (int q : acc[(std::size_t)cy * gx + cx])
                        if ((seeds[q] - p).norm() < dmin) { ok = false; break; }
                }
            if (!ok) continue;
            acc[(std::size_t)cj * gx + ci].push_back((int)seeds.size());
            seeds.push_back(p);
""","""            Vec2 p(Ux(rng), Uy(rng));
            const double spi = drawSp();
            int ci = std::clamp((int)(p.x() / W * gx), 0, gx - 1);
            int cj = std::clamp((int)(p.y() / H * gy), 0, gy - 1);
            bool ok = true;
            for (int dj = -1; dj <= 1 && ok; ++dj)
                for (int di = -1; di <= 1 && ok; ++di) {
                    int cx = ci + di, cy = cj + dj;
                    if (cx < 0 || cy < 0 || cx >= gx || cy >= gy) continue;
                    for (int q : acc[(std::size_t)cy * gx + cx]) {
                        const double lim = (sizeSpread > 0.0) ? 0.35 * (spi + sp[q]) : dmin;
                        if ((seeds[q] - p).norm() < lim) { ok = false; break; }
                    }
                }
            if (!ok) continue;
            acc[(std::size_t)cj * gx + ci].push_back((int)seeds.size());
            seeds.push_back(p);
            sp.push_back(spi);
""","cpp seeds b")
c=rep(c,"""    } else {
        std::uniform_real_distribution<double> U(-1.0, 1.0);
        for (int j = 0;; ++j) {""","""    } else {
        if (sizeSpread > 0.0)
            throw std::runtime_error("Tessellation: grainSizeSpread exige "
                                     "grainSeeding = random (le reseau "
                                     "hexagonal n a qu une taille)");
        std::uniform_real_distribution<double> U(-1.0, 1.0);
        for (int j = 0;; ++j) {""","cpp hex guard")
# ---- cpp : affectation des phases par classe de taille ----
c=rep(c,"""    std::vector<int> order(cells.size());
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), rng);
    std::vector<double> assigned(frac.size(), 0.0);
    T.phaseOfGrain.assign(cells.size(), 0);
    for (int g : order) {
        int best = 0;
        double bestDef = -1e300;
        for (std::size_t p = 0; p < frac.size(); ++p) {
            double def = frac[p] * Atot - assigned[p];
            if (def > bestDef) { bestDef = def; best = (int)p; }
        }
        T.phaseOfGrain[g] = best;
        assigned[best] += T.grainArea[g];
    }
""","""    std::vector<int> order(cells.size());
    std::iota(order.begin(), order.end(), 0);
    std::vector<double> assigned(frac.size(), 0.0);
    T.phaseOfGrain.assign(cells.size(), 0);
    bool anyPhaseSize = false;
    for (double v : phaseSize) if (v > 0.0) anyPhaseSize = true;
    if (!anyPhaseSize) {
        // chemin historique : grains melanges, glouton sur le deficit d aire
        std::shuffle(order.begin(), order.end(), rng);
        for (int g : order) {
            int best = 0;
            double bestDef = -1e300;
            for (std::size_t p = 0; p < frac.size(); ++p) {
                double def = frac[p] * Atot - assigned[p];
                if (def > bestDef) { bestDef = def; best = (int)p; }
            }
            T.phaseOfGrain[g] = best;
            assigned[best] += T.grainArea[g];
        }
    } else {
        // TAILLE PAR PHASE (2026-09-02) : grains pris du plus grand au plus
        // petit ; chaque phase a une affinite log-normale (sigma 0,35) a sa
        // taille cible, multipliee par son deficit RELATIF d aire. Le deficit
        // est en aire, donc les fractions globales sont respectees quelle que
        // soit la distribution des tailles ; l affinite ne fait que choisir
        // QUELS grains vont a quelle phase. Une phase saturee (deficit <= 0)
        // ne recrute plus tant qu une autre a du deficit.
        if (phaseSize.size() != frac.size())
            throw std::runtime_error("Tessellation: phaseSize doit avoir une "
                                     "entree par phase");
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            if (T.grainArea[a] != T.grainArea[b]) return T.grainArea[a] > T.grainArea[b];
            return a < b;
        });
        for (int g : order) {
            const double dg = 2.0 * std::sqrt(T.grainArea[g] / M_PI);
            int best = -1;
            double bestScore = -1e300;
            for (std::size_t p = 0; p < frac.size(); ++p) {
                const double def = frac[p] * Atot - assigned[p];
                if (def <= 0.0) continue;
                double aff = 1.0;
                if (phaseSize[p] > 0.0) {
                    const double z = std::log(dg / phaseSize[p]) / 0.35;
                    aff = std::exp(-0.5 * z * z);
                }
                const double score = aff * def / (frac[p] * Atot);
                if (score > bestScore) { bestScore = score; best = (int)p; }
            }
            if (best < 0) {                        // toutes saturees : max deficit
                double bestDef = -1e300;
                for (std::size_t p = 0; p < frac.size(); ++p) {
                    const double def = frac[p] * Atot - assigned[p];
                    if (def > bestDef) { bestDef = def; best = (int)p; }
                }
            }
            T.phaseOfGrain[g] = best;
            assigned[best] += T.grainArea[g];
        }
    }
""","cpp phases")
# ---- FdemSolver : cles, passage, journal ----
f=rep(f,"""    double gh = cfg_.getd("grainElemSize", 0.0);
    Tessellation T = Tessellation::build(W_, H_, d, jit, lloyd, mf, refine,
                                         phases_.fraction, rng,
                                         seeding == "random",
                                         gm == "delaunay", gh);
""","""    double gh = cfg_.getd("grainElemSize", 0.0);
    // ---- POLYDISPERSITE (2026-09-02, opt-in) ---------------------------------
    // grainSizeSpread = ecart-type de ln(taille) des grains (0 = historique) ;
    // phase.<nom>.grainSize = taille cible par phase (fractions conservees).
    const double gSpread = cfg_.getd("grainSizeSpread", 0.0);
    if (gSpread < 0.0 || gSpread > 1.5)
        throw std::runtime_error("grainSizeSpread doit etre dans [0 ; 1,5] "
                                 "(ecart-type de ln(taille) ; 0,3 = modere, 0,6 = fort)");
    if (gSpread > 0.0 && seeding != "random")
        throw std::runtime_error("grainSizeSpread exige grainSeeding = random");
    std::vector<double> pSize;
    bool anyPSize = false;
    for (const std::string& nm : phases_.name) {
        const double v = cfg_.getd("phase." + nm + ".grainSize", -1.0);
        pSize.push_back(v);
        if (v > 0.0) anyPSize = true;
    }
    if (anyPSize && phases_.n() < 2)
        throw std::runtime_error("phase.<nom>.grainSize n a de sens qu avec "
                                 "plusieurs phases");
    if (!anyPSize) pSize.clear();
    Tessellation T = Tessellation::build(W_, H_, d, jit, lloyd, mf, refine,
                                         phases_.fraction, rng,
                                         seeding == "random",
                                         gm == "delaunay", gh,
                                         gSpread, pSize);
    if (gSpread > 0.0 || anyPSize) {
        // journal : distribution REALISEE (apres Lloyd et contraction), par
        // phase — c est ce chiffre qui compte, pas la cle demandee.
        std::vector<double> aSum(phases_.n(), 0.0), dSum(phases_.n(), 0.0),
                            d2Sum(phases_.n(), 0.0);
        std::vector<int> nG(phases_.n(), 0);
        double aTot = 0.0, lnSum = 0.0, ln2Sum = 0.0;
        for (int g = 0; g < T.nGrains; ++g) {
            const int p = T.phaseOfGrain[g];
            const double a = T.grainArea[g], dg = 2.0 * std::sqrt(a / M_PI);
            aSum[p] += a; dSum[p] += dg; d2Sum[p] += dg * dg; ++nG[p]; aTot += a;
            lnSum += std::log(dg); ln2Sum += std::log(dg) * std::log(dg);
        }
        const double n = (double)T.nGrains;
        const double sdLn = std::sqrt(std::max(0.0, ln2Sum / n - (lnSum / n) * (lnSum / n)));
        std::cout << "[FDEM] POLYDISPERSITE : " << T.nGrains << " grains, ecart-type de "
                     "ln(d_eq) REALISE = " << sdLn << " (demande " << gSpread
                  << ", Lloyd " << lloyd << " iterations homogeneise)" << std::endl;
        for (int p = 0; p < phases_.n(); ++p) {
            const double m = nG[p] ? dSum[p] / nG[p] : 0.0;
            const double sd = nG[p] ? std::sqrt(std::max(0.0, d2Sum[p] / nG[p] - m * m)) : 0.0;
            std::cout << "[FDEM]   " << phases_.name[p] << " : fraction d aire "
                      << 100.0 * aSum[p] / aTot << " % (cible " << 100.0 * phases_.fraction[p]
                      << " %), " << nG[p] << " grains, d_eq = " << 1000.0 * m << " +- "
                      << 1000.0 * sd << " mm" << (pSize.size() && pSize[p] > 0.0
                          ? " (cible " + std::to_string(1000.0 * pSize[p]) + " mm)" : "")
                      << std::endl;
        }
        if (gSpread > 0.0 && lloyd > 1)
            std::cout << "[FDEM] WARNING: lloydIters = " << lloyd << " avec grainSizeSpread : "
                         "chaque iteration de Lloyd rapproche les cellules d une taille "
                         "unique — comparer l ecart-type REALISE a la demande, et "
                         "reduire lloydIters si l ecart est trop grand" << std::endl;
    }
""","fdem call")
io.open(H,"w",encoding="utf-8").write(h); io.open(C,"w",encoding="utf-8").write(c); io.open(F,"w",encoding="utf-8").write(f)
print("polydispersite : 7 greffes appliquees (hpp 1, Tessellation.cpp 5, FdemSolver.cpp 1)")
