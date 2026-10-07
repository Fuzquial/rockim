import io, os
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
H="include/rockim/Tessellation.hpp"; C="src/Tessellation.cpp"
h=io.open(H,encoding="utf-8").read(); c=io.open(C,encoding="utf-8").read()
start = "    if (randomSeeds) {\n        // Poisson-disc dart throwing: uniform candidates, accepted when at\n"
end = "            sp.push_back(spi);\n        }\n    } else {\n"
i0 = c.index(start); i1 = c.index(end, i0) + len(end)
assert c.count(start) == 1
NL = chr(92) + "n"
new = '''    if (randomSeeds && sizeSpread > 0.0) {
        // POLYDISPERSITE (2026-09-02) : addition sequentielle TRIEE. Les N
        // espacements s_i = s L_i (L_i log-normale de moyenne 1, ecart-type
        // de ln = sizeSpread) sont tires d abord, puis places du plus GRAND
        // au plus petit, chaque graine acceptee a distance
        // >= 0,35 (s_i + s_j) de toutes les precedentes. Tirer la taille a
        // chaque essai (Poisson-disc ordinaire) rejette surtout les grandes
        // graines, qui exigent le plus de place : l ecart-type realise
        // s effondrait (0,5 demande -> 0,20 obtenu). N est divise par
        // E[L^2] = exp(sigma^2) pour conserver la taille moyenne.
        long target = std::lround(W * H / (std::sqrt(3.0) / 2.0 * s * s
                                           * std::exp(sizeSpread * sizeSpread)));
        target = std::max(2L, target);
        std::normal_distribution<double> N01(0.0, 1.0);
        std::vector<double> want((std::size_t)target);
        for (double& v : want)
            v = s * std::exp(sizeSpread * N01(rng) - 0.5 * sizeSpread * sizeSpread);
        std::sort(want.begin(), want.end(), std::greater<double>());
        const double cs = 0.7 * want.front();     // >= toute distance d acceptation
        int gx = std::max(1, (int)(W / cs)), gy = std::max(1, (int)(H / cs));
        std::vector<std::vector<int>> acc((std::size_t)gx * gy);
        std::uniform_real_distribution<double> Ux(eps, W - eps);
        std::uniform_real_distribution<double> Uy(eps, H - eps);
        int skipped = 0;
        for (double spi : want) {
            bool placed = false;
            for (int tr = 0; tr < 4000 && !placed; ++tr) {
                Vec2 p(Ux(rng), Uy(rng));
                int ci = std::clamp((int)(p.x() / W * gx), 0, gx - 1);
                int cj = std::clamp((int)(p.y() / H * gy), 0, gy - 1);
                bool ok = true;
                for (int dj = -1; dj <= 1 && ok; ++dj)
                    for (int di = -1; di <= 1 && ok; ++di) {
                        int cx = ci + di, cy = cj + dj;
                        if (cx < 0 || cy < 0 || cx >= gx || cy >= gy) continue;
                        for (int q : acc[(std::size_t)cy * gx + cx])
                            if ((seeds[q] - p).norm() < 0.35 * (spi + sp[q])) { ok = false; break; }
                    }
                if (!ok) continue;
                acc[(std::size_t)cj * gx + ci].push_back((int)seeds.size());
                seeds.push_back(p);
                sp.push_back(spi);
                placed = true;
            }
            if (!placed) ++skipped;
        }
        if (skipped > 0)
            std::printf("[tess] polydispersite: %d graine(s) sur %ld non placee(s) "
                        "(domaine sature) - l ecart-type realise en tient compte''' + NL + '''",
                        skipped, target);
    } else if (randomSeeds) {
        // Poisson-disc dart throwing: uniform candidates, accepted when at
        // least dmin from every accepted seed. Unlike the jittered hex
        // lattice (whose Voronoi keeps three preferred boundary
        // orientations even after Lloyd), the resulting boundary
        // orientation distribution is isotropic.
        double dmin = 0.7 * s;
        long target = std::lround(W * H / (std::sqrt(3.0) / 2.0 * s * s));
        double cs = dmin;
        int gx = std::max(1, (int)(W / cs)), gy = std::max(1, (int)(H / cs));
        std::vector<std::vector<int>> acc((std::size_t)gx * gy);
        std::uniform_real_distribution<double> Ux(eps, W - eps);
        std::uniform_real_distribution<double> Uy(eps, H - eps);
        long attempts = 0, maxAttempts = 400 * std::max(1L, target);
        while ((long)seeds.size() < target && attempts < maxAttempts) {
            ++attempts;
            Vec2 p(Ux(rng), Uy(rng));
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
        }
    } else {
'''
c = c[:i0] + new + c[i1:]
a = "    //            espacement s_i = s L_i, L_i log-normale de moyenne 1, deux\n    //            graines s acceptent a distance >= 0,35 (s_i + s_j), et les\n"
b = "    //            espacement s_i = s L_i, L_i log-normale de moyenne 1 ; les\n    //            graines sont placees du plus GRAND au plus petit (addition\n    //            sequentielle triee), acceptees a >= 0,35 (s_i + s_j), et les\n"
assert h.count(a) == 1; h = h.replace(a, b, 1)
io.open(H,"w",encoding="utf-8").write(h); io.open(C,"w",encoding="utf-8").write(c)
print("tri sequentiel : branche polydisperse separee, Poisson-disc historique restaure verbatim")
