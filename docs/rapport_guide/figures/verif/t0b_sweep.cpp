// Rejeu hors depot du banc T0b (C8/C9 de selftest-toolcontact) : copie
// textuelle de la lambda `bar` de src/FdemSolver.cpp (lignes ~10836-10900),
// appel du VRAI noyau include/rockim/ToolSignorini.hpp, balayage en N.
#include "rockim/ToolSignorini.hpp"
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>
using namespace rockim::toolsig;
int main() {
    const double Ebar = 48.26e9, rhoBar = 2610.0, Lbar = 0.1;
    const double cBar = std::sqrt(Ebar / rhoBar);
    const double vW = 10.0;
    const double tContactTh = 2.0 * Lbar / cBar;
    FILE* hist = std::fopen("t0b_histoire_N100.csv", "w");
    std::fprintf(hist, "t_us,v_moy,v_noeud0,v_noeudN,F_mur_N_m2,Ec_J_m2,Ee_J_m2,W_mur_J_m2\n");
    auto bar = [&](int N, double damping, double& vOut, double& tCon,
                   double& lossRel, bool rec) {
        const double M = rhoBar * Lbar;
        const double m = M / N;
        const double kSp = N * Ebar / Lbar;
        const double dt = 0.2 * 2.0 / (2.0 * std::sqrt(kSp / m));
        std::vector<double> u(N, 0.0), v(N, 0.0), f(N, 0.0);
        double xW = -1e-9;
        double work = 0.0;
        long firstTouch = -1, lastTouch = -1;
        const long nSteps = (long)(6.0 * tContactTh / dt);
        for (long s = 0; s < nSteps; ++s) {
            for (int i = 0; i < N; ++i) f[i] = 0.0;
            for (int i = 0; i + 1 < N; ++i) {
                double fs = kSp * (u[i + 1] - u[i]);
                f[i] += fs; f[i + 1] -= fs;
            }
            if (damping > 0.0)
                for (int i = 0; i < N; ++i)
                    if (v[i] != 0.0)
                        f[i] -= damping * std::abs(f[i]) * (v[i] > 0.0 ? 1.0 : -1.0);
            std::vector<double> vFree(N);
            for (int i = 0; i < N; ++i) vFree[i] = v[i] + dt * f[i] / m;
            double pen = xW - u[0];
            double rnStep = 0.0;
            if (pen > 0.0) {
                Impulse r = impulse(pen, vFree[0] - vW, 0.0, m, dt, 0.0, 0.0);
                if (r.active) {
                    vFree[0] += r.rn / m;
                    work += r.rn * vW;
                    rnStep = r.rn;
                    if (firstTouch < 0) firstTouch = s;
                    lastTouch = s;
                }
            }
            for (int i = 0; i < N; ++i) { v[i] = vFree[i]; u[i] += dt * v[i]; }
            xW += dt * vW;
            if (rec && (s % 20 == 0)) {
                double vm = 0, ke = 0, ue = 0;
                for (int i = 0; i < N; ++i) { vm += v[i]; ke += 0.5 * m * v[i] * v[i]; }
                for (int i = 0; i + 1 < N; ++i) { double e = u[i+1]-u[i]; ue += 0.5*kSp*e*e; }
                std::fprintf(hist, "%.6f,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n", s*dt*1e6, vm/N, v[0], v[N-1],
                             rnStep/dt, ke, ue, work);
            }
            if (lastTouch >= 0 && s > lastTouch + 200 && u[0] > xW) break;
        }
        double vm = 0.0, ke = 0.0, ue = 0.0;
        for (int i = 0; i < N; ++i) { vm += v[i]; ke += 0.5 * m * v[i] * v[i]; }
        for (int i = 0; i + 1 < N; ++i) { double e = u[i + 1] - u[i]; ue += 0.5 * kSp * e * e; }
        vOut = vm / N;
        tCon = (lastTouch - firstTouch + 1) * dt;
        lossRel = (work - ke - ue) / (0.5 * M * vW * vW);
    };
    std::printf("N,amort,v_sortie_m_s,v_sur_2v,t_contact_us,t_sur_2Lc,perte_rel,perte_x_N\n");
    int Ns[] = {10, 25, 50, 100, 200, 400, 800, 1600};
    for (double d : {0.0, 0.05})
      for (int N : Ns) {
        double vo, tc, lr;
        bar(N, d, vo, tc, lr, (N == 100 && d == 0.0));
        std::printf("%d,%.2f,%.6f,%.6f,%.4f,%.6f,%.6e,%.5f\n", N, d, vo, vo/(2*vW), tc*1e6, tc/tContactTh, lr, lr*N);
      }
    std::fclose(hist);
    std::fprintf(stderr, "2L/c = %.4f us, c = %.2f m/s\n", tContactTh*1e6, cBar);
}
