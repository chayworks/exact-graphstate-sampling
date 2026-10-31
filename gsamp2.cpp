// In-place exact sampler (see inplace.py): one complex64 arena of 2^R entries.
// Plan: n nsteps R; per step: v nk m g rhoc nterms; nk Tkeep values; nterms x (Tq x tau par); 8 doubles M[b][x].
// Usage: gsamp2 plan sample N seed | gsamp2 plan prob stringsfile   (same output format as gsamp)
//        gsamp2 plan chain stringsfile   -> per line: log2 Pr, then Pr(b_k | b_1..b_{k-1}) for each step (one pass per step)
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>
#include <random>
#include <memory>
typedef std::complex<double> cd; typedef std::complex<float> cf;
struct Term { uint64_t tq; int x, tau, par; };
struct Step { int v, nk, m; uint64_t g, rhoc; std::vector<uint64_t> tk; std::vector<Term> t; cd M[2][2]; uint64_t tab[5][256]; };
static std::vector<Step> st; static int n, ns, R; static std::unique_ptr<cf[]> mem;

static inline uint64_t base_of(const Step& s, uint64_t c) {
    return s.tab[0][c & 255] ^ s.tab[1][(c >> 8) & 255] ^ s.tab[2][(c >> 16) & 255] ^ s.tab[3][(c >> 24) & 255] ^ s.tab[4][(c >> 32) & 255];
}
static inline void outputs(const Step& s, uint64_t c, uint64_t base, cd out[2][2]) {   // out[b][tau]
    out[0][0] = out[0][1] = out[1][0] = out[1][1] = 0; const int pc = __builtin_parityll(c & s.rhoc);
    for (const Term& t : s.t) {
        cd o = cd(mem[base ^ t.tq]); if (t.x && (pc ^ t.par)) o = -o;
        out[0][t.tau] += s.M[0][t.x] * o; out[1][t.tau] += s.M[1][t.x] * o;
    }
}
static double chain(const int* forced, int* bits, std::mt19937_64* rng) {
    mem[0] = 1; double logp = 0;
    for (int k = 0; k < ns; k++) {
        const Step& s = st[k]; const int64_t N = (int64_t)1 << s.nk; double n0 = 0, n1 = 0;
#pragma omp parallel for schedule(static) reduction(+:n0,n1)
        for (int64_t c = 0; c < N; c++) { cd o[2][2]; uint64_t b = base_of(s, c); outputs(s, c, b, o);
            n0 += std::norm(o[0][0]) + std::norm(o[0][1]); n1 += std::norm(o[1][0]) + std::norm(o[1][1]); }
        double p1 = n1 / (n0 + n1); int bv;
        if (forced) bv = forced[k]; else { std::uniform_real_distribution<double> U(0, 1); bv = U(*rng) < p1; }
        double pb = bv ? p1 : 1 - p1; if (pb <= 0) return -INFINITY; logp += std::log2(pb); bits[k] = bv;
        const double sc = 1.0 / std::sqrt(bv ? n1 : n0);
#pragma omp parallel for schedule(static)
        for (int64_t c = 0; c < N; c++) { cd o[2][2]; uint64_t b = base_of(s, c); outputs(s, c, b, o);
            mem[b] = cf(o[bv][0] * sc); if (s.m == 2) mem[b ^ s.g] = cf(o[bv][1] * sc); }
    }
    return logp;
}
// Single pass per step for a forced string: store the chosen branch unnormalised and accumulate its norm nb.
// Since sum_b Pr(prefix, b) = Pr(prefix), n0 + n1 = 2 * (norm of the current table) = 2 * nb_prev, so the
// conditional is nb / (2 nb_prev) without the norm pass. cond[k] receives Pr(b_k | b_1..b_{k-1}).
static double chain1(const int* forced, double* cond) {
    mem[0] = 1; double prev = 1;
    for (int k = 0; k < ns; k++) {
        const Step& s = st[k]; const int64_t N = (int64_t)1 << s.nk; const int bv = forced[k]; double nb = 0;
#pragma omp parallel for schedule(static) reduction(+:nb)
        for (int64_t c = 0; c < N; c++) { cd o[2][2]; uint64_t b = base_of(s, c); outputs(s, c, b, o);
            mem[b] = cf(o[bv][0]); if (s.m == 2) mem[b ^ s.g] = cf(o[bv][1]);
            nb += std::norm(o[bv][0]) + std::norm(o[bv][1]); }
        cond[k] = nb / (2 * prev); prev = nb; if (nb <= 0) { for (; k < ns; k++) cond[k] = 0; return -INFINITY; }
    }
    return std::log2(prev) - ns;
}
int main(int argc, char** argv) {
    if (argc < 4) return 1;
    FILE* f = fopen(argv[1], "r"); if (!f || fscanf(f, "%d %d %d", &n, &ns, &R) != 3) return 2; st.resize(ns);
    for (auto& s : st) {
        unsigned long long g, rc; int nt;
        if (fscanf(f, "%d %d %d %llu %llu %d", &s.v, &s.nk, &s.m, &g, &rc, &nt) != 6) return 3; s.g = g; s.rhoc = rc;
        s.tk.resize(s.nk); for (auto& x : s.tk) { unsigned long long y; if (fscanf(f, "%llu", &y) != 1) return 4; x = y; }
        s.t.resize(nt); for (auto& t : s.t) { unsigned long long y; if (fscanf(f, "%llu %d %d %d", &y, &t.x, &t.tau, &t.par) != 4) return 5; t.tq = y; }
        for (int b = 0; b < 2; b++) for (int x = 0; x < 2; x++) { double re, im; if (fscanf(f, "%lf %lf", &re, &im) != 2) return 6; s.M[b][x] = cd(re, im); }
        if (s.nk > 40) return 7;
        for (int t = 0; t < 5; t++) for (int v = 0; v < 256; v++) { uint64_t a = 0;
            for (int j = 0; j < 8; j++) if ((v >> j) & 1 && 8 * t + j < s.nk) a ^= s.tk[8 * t + j]; s.tab[t][v] = a; }
    }
    fclose(f); mem.reset(new cf[(size_t)1 << R]);
    std::vector<int> bits(ns);
    if (!strcmp(argv[2], "sample")) {
        long N = atol(argv[3]); unsigned long long seed = strtoull(argv[4], 0, 10);
        for (long i = 0; i < N; i++) { std::mt19937_64 rng(seed + i); double lp = chain(nullptr, bits.data(), &rng);
            for (int k = 0; k < ns; k++) putchar('0' + bits[k]); printf(" %.12f %llu\n", lp, seed + i); fflush(stdout); }
    } else {
        const bool ch = !strcmp(argv[2], "chain"); std::vector<double> cond(ns);
        FILE* g = fopen(argv[3], "r"); char line[4096];
        while (fgets(line, sizeof line, g)) { if ((int)strlen(line) < ns) continue; std::vector<int> fb(ns);
            for (int k = 0; k < ns; k++) fb[k] = line[k] - '0';
            if (ch) { printf("%.12f", chain1(fb.data(), cond.data())); for (int k = 0; k < ns; k++) printf(" %.9f", cond[k]); putchar('\n'); }
            else printf("%.12f\n", chain(fb.data(), bits.data(), nullptr));
            fflush(stdout); }
    }
    return 0;
}
