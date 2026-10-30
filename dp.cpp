// Cut-rank DP executor:  Z(w) = sum_x prod_j w_j(x_j) (-1)^{sum_{i<j} A_ij x_i x_j}
// Plan (from plan.py, written by engine.py):  m nsteps maxr, then per step: v r rn rho pivmask nker ker...
// Jobs (stdin-like file):  njobs; per job: nslices; per slice: pre_re pre_im then m lines w0re w0im w1re w1im.
// Output: one line per job, "re im" of sum over slices of pre * Z(w_slice).
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <complex>
#include <vector>
#include <immintrin.h>
typedef std::complex<double> cd;
struct Step { int v, r, rn; uint64_t rho, pivmask; std::vector<uint64_t> kspan; };

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: dp plan.txt jobs.txt\n"); return 1; }
    FILE* f = fopen(argv[1], "r"); int m, ns, maxr;
    if (fscanf(f, "%d %d %d", &m, &ns, &maxr) != 3) return 2;
    std::vector<Step> st(ns);
    for (auto& s : st) {
        unsigned long long rho, pm; int nk;
        if (fscanf(f, "%d %d %d %llu %llu %d", &s.v, &s.r, &s.rn, &rho, &pm, &nk) != 6) return 3;
        s.rho = rho; s.pivmask = pm; s.kspan.assign(1, 0);
        for (int i = 0; i < nk; i++) { unsigned long long k; if (fscanf(f, "%llu", &k) != 1) return 4;
            size_t L = s.kspan.size(); for (size_t j = 0; j < L; j++) s.kspan.push_back(s.kspan[j] ^ k); }
    }
    fclose(f);
    std::vector<cd> a((size_t)1 << maxr), b((size_t)1 << maxr);
    FILE* g = fopen(argv[2], "r"); int nj; if (fscanf(g, "%d", &nj) != 1) return 5;
    std::vector<cd> w0(m), w1(m);
    for (int j = 0; j < nj; j++) {
        int nsl; if (fscanf(g, "%d", &nsl) != 1) return 6; cd tot = 0;
        for (int sl = 0; sl < nsl; sl++) {
            double pr, pi; if (fscanf(g, "%lf %lf", &pr, &pi) != 2) return 7;
            for (int q = 0; q < m; q++) { double x0, y0, x1, y1;
                if (fscanf(g, "%lf %lf %lf %lf", &x0, &y0, &x1, &y1) != 4) return 8; w0[q] = cd(x0, y0); w1[q] = cd(x1, y1); }
            cd* old = a.data(); cd* nw = b.data(); old[0] = 1;
            for (const Step& s : st) {
                const uint64_t lowmask = (s.r == 64) ? ~0ULL : ((1ULL << s.r) - 1), xbit = 1ULL << s.r;
                const cd W0 = w0[s.v], W1 = w1[s.v]; const int nk = s.kspan.size();
                const uint64_t* ks = s.kspan.data(); const uint64_t pm = s.pivmask, rho = s.rho;
                const int64_t N = (int64_t)1 << s.rn;
#pragma omp parallel for schedule(static)
                for (int64_t cp = 0; cp < N; cp++) {
                    uint64_t u0 = _pdep_u64((uint64_t)cp, pm); cd acc = 0;
                    for (int t = 0; t < nk; t++) {
                        uint64_t u = u0 ^ ks[t], c = u & lowmask;
                        if (u & xbit) acc += (__builtin_parityll(c & rho) ? -W1 : W1) * old[c];
                        else acc += W0 * old[c];
                    }
                    nw[cp] = acc;
                }
                std::swap(old, nw);
            }
            tot += cd(pr, pi) * old[0];
        }
        printf("%.17g %.17g\n", tot.real(), tot.imag()); fflush(stdout);
    }
    return 0;
}
