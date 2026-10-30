// Exact sampler for graph state + product measurement:  Pr(b_1..b_t) = 2^-t sum_c |P_t[c]|^2.
// Plan file (gsamp_plan from engine.py): n nsteps maxr; per step: v r rn rho pivmask nker ker... then
//   M[b=0][x=0] M[0][1] M[1][0] M[1][1] as 8 doubles (re im ...).
// Usage: gsamp plan sample N seed      -> per sample: bits in step order, log2 Pr(sample)
//        gsamp plan prob stringsfile   -> per string (bits in step order): log2 Pr(string)  (exact, via the chain)
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>
#include <random>
#include <string>
#include <immintrin.h>
typedef std::complex<double> cd; typedef std::complex<float> cf;
struct Step { int v, r, rn; uint64_t rho, pivmask; std::vector<uint64_t> kspan; cd M[2][2]; };
static std::vector<Step> st; static int n, ns, maxr;
static std::vector<cf> bufA, bufB;

// Gather for one step: S0[c'] (x=0 terms), S1[c'] (x=1 terms incl. sign); P_b = M[b][0] S0 + M[b][1] S1.
static inline void gather(const Step& s, const cf* old, int64_t cp, cd& S0, cd& S1) {
    const uint64_t lowmask = (1ULL << s.r) - 1, xbit = 1ULL << s.r;
    uint64_t u0 = _pdep_u64((uint64_t)cp, s.pivmask); S0 = 0; S1 = 0;
    for (size_t t = 0; t < s.kspan.size(); t++) {
        uint64_t u = u0 ^ s.kspan[t], c = u & lowmask; cd o = cd(old[c]);
        if (u & xbit) S1 += __builtin_parityll(c & s.rho) ? -o : o; else S0 += o;
    }
}

// Runs the chain. If forced != nullptr, follows those bits; else samples with rng. Returns log2 Pr.
static double chain(const int* forced, int* bits, std::mt19937_64* rng) {
    cf* old = bufA.data(); cf* nw = bufB.data(); old[0] = 1; double logp = 0;
    for (int k = 0; k < ns; k++) {
        const Step& s = st[k]; const int64_t N = (int64_t)1 << s.rn;
        double n0 = 0, n1 = 0;
#pragma omp parallel for schedule(static) reduction(+:n0,n1)
        for (int64_t cp = 0; cp < N; cp++) {
            cd S0, S1; gather(s, old, cp, S0, S1);
            n0 += std::norm(s.M[0][0] * S0 + s.M[0][1] * S1); n1 += std::norm(s.M[1][0] * S0 + s.M[1][1] * S1);
        }
        double p1 = n1 / (n0 + n1); int b;
        if (forced) b = forced[k]; else { std::uniform_real_distribution<double> U(0, 1); b = U(*rng) < p1; }
        double pb = b ? p1 : 1 - p1; if (pb <= 0) return -INFINITY;
        logp += std::log2(pb); bits[k] = b;
        const double scale = 1.0 / std::sqrt(b ? n1 : n0); const cd m0 = s.M[b][0] * scale, m1 = s.M[b][1] * scale;
#pragma omp parallel for schedule(static)
        for (int64_t cp = 0; cp < N; cp++) { cd S0, S1; gather(s, old, cp, S0, S1); nw[cp] = cf(m0 * S0 + m1 * S1); }
        std::swap(old, nw);
    }
    return logp;
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage\n"); return 1; }
    FILE* f = fopen(argv[1], "r"); if (!f || fscanf(f, "%d %d %d", &n, &ns, &maxr) != 3) return 2;
    st.resize(ns);
    for (auto& s : st) {
        unsigned long long rho, pm; int nk;
        if (fscanf(f, "%d %d %d %llu %llu %d", &s.v, &s.r, &s.rn, &rho, &pm, &nk) != 6) return 3;
        s.rho = rho; s.pivmask = pm; s.kspan.assign(1, 0);
        for (int i = 0; i < nk; i++) { unsigned long long k; if (fscanf(f, "%llu", &k) != 1) return 4;
            size_t L = s.kspan.size(); for (size_t j = 0; j < L; j++) s.kspan.push_back(s.kspan[j] ^ k); }
        for (int b = 0; b < 2; b++) for (int x = 0; x < 2; x++) { double re, im; if (fscanf(f, "%lf %lf", &re, &im) != 2) return 5; s.M[b][x] = cd(re, im); }
    }
    fclose(f);
    bufA.assign((size_t)1 << maxr, 0); bufB.assign((size_t)1 << maxr, 0);
    std::vector<int> bits(ns);
    if (!strcmp(argv[2], "sample")) {
        long N = atol(argv[3]); unsigned long long seed = strtoull(argv[4], 0, 10);
        for (long i = 0; i < N; i++) {
            std::mt19937_64 rng(seed + i); double lp = chain(nullptr, bits.data(), &rng);
            for (int k = 0; k < ns; k++) putchar('0' + bits[k]); printf(" %.12f %llu\n", lp, seed + i); fflush(stdout);
        }
    } else if (!strcmp(argv[2], "prob")) {
        FILE* g = fopen(argv[3], "r"); char line[4096];
        while (fgets(line, sizeof line, g)) {
            if ((int)strlen(line) < ns) continue; std::vector<int> fb(ns);
            for (int k = 0; k < ns; k++) fb[k] = line[k] - '0';
            printf("%.12f\n", chain(fb.data(), bits.data(), nullptr)); fflush(stdout);
        }
    }
    return 0;
}
