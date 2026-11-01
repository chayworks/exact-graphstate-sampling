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

// One pass per step for a forced string: write the chosen branch unnormalised and accumulate its norm nb.
// sum_b Pr(prefix, b) = Pr(prefix) gives n0 + n1 = 2 * (norm of the current table) = 2 * nb_prev, so
// Pr(b_k | prefix) = nb / (2 nb_prev) without the norm pass. cond[k] receives that conditional.
static double chain1(const int* forced, double* cond) {
    cf* old = bufA.data(); cf* nw = bufB.data(); old[0] = 1; double prev = 1;
    for (int k = 0; k < ns; k++) {
        const Step& s = st[k]; const int64_t N = (int64_t)1 << s.rn; const int b = forced[k];
        const cd m0 = s.M[b][0], m1 = s.M[b][1]; double nb = 0;
#pragma omp parallel for schedule(static) reduction(+:nb)
        for (int64_t cp = 0; cp < N; cp++) { cd S0, S1; gather(s, old, cp, S0, S1); cd o = m0 * S0 + m1 * S1; nw[cp] = cf(o); nb += std::norm(o); }
        cond[k] = nb / (2 * prev); prev = nb; std::swap(old, nw);
        if (nb <= 0) { for (; k < ns; k++) cond[k] = 0; return -INFINITY; }
    }
    return std::log2(prev) - ns;
}

// Hybrid machinery (as in gsamp2.cpp): steps k0..k1-1 on the current table `old` of norm *prev, stored scaled so that
// true table = stored * 2^*logscale. forced[k] in {0,1}: single pass; forced[k] == -1: norm pass, draw, write pass.
static cf *cur = nullptr, *nxt = nullptr;
static double run_steps(int k0, int k1, const int* forced, int* bits, double* cond, std::mt19937_64* rng, double* prev, double* logscale) {
    double logp = 0;
    for (int k = k0; k < k1; k++) {
        const Step& s = st[k]; const int64_t N = (int64_t)1 << s.rn; int b = forced[k]; const double sc = 1.0 / std::sqrt(*prev);
        if (b < 0) {
            double n0 = 0, n1 = 0;
#pragma omp parallel for schedule(static) reduction(+:n0,n1)
            for (int64_t cp = 0; cp < N; cp++) { cd S0, S1; gather(s, cur, cp, S0, S1);
                n0 += std::norm(s.M[0][0] * S0 + s.M[0][1] * S1); n1 += std::norm(s.M[1][0] * S0 + s.M[1][1] * S1); }
            std::uniform_real_distribution<double> U(0, 1); b = U(*rng) < n1 / (n0 + n1);
        }
        const cd m0 = s.M[b][0] * sc, m1 = s.M[b][1] * sc; double nb = 0;
#pragma omp parallel for schedule(static) reduction(+:nb)
        for (int64_t cp = 0; cp < N; cp++) { cd S0, S1; gather(s, cur, cp, S0, S1); cd o = m0 * S0 + m1 * S1; nxt[cp] = cf(o); nb += std::norm(o); }
        *logscale += 0.5 * std::log2(*prev); cond[k] = nb / 2; bits[k] = b; *prev = nb; std::swap(cur, nxt);
        if (nb <= 0) { for (; k < k1; k++) { cond[k] = 0; bits[k] = forced[k] < 0 ? 0 : forced[k]; } return -INFINITY; }
        logp += std::log2(cond[k]);
    }
    return logp;
}
static double hybrid(const int* forced, int* bits, double* cond, std::mt19937_64* rng) {
    cur = bufA.data(); nxt = bufB.data(); cur[0] = 1; double prev = 1, logscale = 0; return run_steps(0, ns, forced, bits, cond, rng, &prev, &logscale);
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
    } else if (!strcmp(argv[2], "hybrid")) {
        // gsamp plan hybrid stringsfile seed: '?' = sampled from the exact conditional (only in a trailing data suffix)
        unsigned long long seed = strtoull(argv[4], 0, 10); std::vector<double> cond(ns); FILE* g = fopen(argv[3], "r"); char line[4096]; long li = 0;
        while (fgets(line, sizeof line, g)) { if ((int)strlen(line) < ns) continue; std::vector<int> fb(ns);
            for (int k = 0; k < ns; k++) fb[k] = line[k] == '?' ? -1 : line[k] - '0';
            std::mt19937_64 r2(seed + (unsigned long long)li++); double lp = hybrid(fb.data(), bits.data(), cond.data(), &r2);
            for (int k = 0; k < ns; k++) putchar('0' + bits[k]); printf(" %.12f", lp); for (int k = 0; k < ns; k++) printf(" %.9f", cond[k]); putchar('\n'); fflush(stdout); }
    } else if (!strcmp(argv[2], "prefix")) {
        // gsamp plan prefix stringsfile K outprefix: run forced steps 0..K-1, dump the contiguous table read by step K
        // (binary: int32 r, int32 1, double logscale, 2^r complex64) to outprefix_<line>.bin.  K == ns dumps the amplitude.
        int K = atoi(argv[4]); const char* outp = argv[5]; std::vector<double> cond(ns); FILE* g = fopen(argv[3], "r"); char line[4096]; long li = 0;
        while (fgets(line, sizeof line, g)) { if ((int)strlen(line) < K) continue; std::vector<int> fb(ns, 0);
            for (int k = 0; k < K; k++) { if (line[k] == '?') { fprintf(stderr, "'?' before K\n"); return 8; } fb[k] = line[k] - '0'; }
            cur = bufA.data(); nxt = bufB.data(); cur[0] = 1; double prev = 1, logscale = 0; double lp = run_steps(0, K, fb.data(), bits.data(), cond.data(), nullptr, &prev, &logscale);
            int r = K < ns ? st[K].r : 0, one = 1; char fn[1024]; snprintf(fn, sizeof fn, "%s_%ld.bin", outp, li++); FILE* o = fopen(fn, "wb");
            fwrite(&r, 4, 1, o); fwrite(&one, 4, 1, o); fwrite(&logscale, 8, 1, o); fwrite(cur, sizeof(cf), (size_t)1 << r, o); fclose(o);
            printf("%ld %.12f %.6f\n", li - 1, lp, logscale); fflush(stdout); }
    } else if (!strcmp(argv[2], "resume")) {
        // gsamp plan resume stringsfile K seed dumpprefix1 [dumpprefix2 ...]: sum the slices' tables, continue steps K..ns-1
        int K = atoi(argv[4]); unsigned long long seed = strtoull(argv[5], 0, 10); int nsl = argc - 6; std::vector<double> cond(ns); FILE* g = fopen(argv[3], "r"); char line[4096]; long li = 0;
        int r = K < ns ? st[K].r : 0; int64_t Nt = (int64_t)1 << r; std::vector<std::vector<cf>> tabs(nsl, std::vector<cf>(Nt)); std::vector<double> ls(nsl);
        while (fgets(line, sizeof line, g)) { if ((int)strlen(line) < ns) continue; std::vector<int> fb(ns);
            for (int k = 0; k < ns; k++) fb[k] = line[k] == '?' ? -1 : line[k] - '0';
            double L = -1e300;
            for (int s = 0; s < nsl; s++) { char fn[1024]; snprintf(fn, sizeof fn, "%s_%ld.bin", argv[6 + s], li); FILE* o = fopen(fn, "rb"); if (!o) { fprintf(stderr, "missing %s\n", fn); return 9; }
                int r2, one; if (fread(&r2, 4, 1, o) != 1 || r2 != r || fread(&one, 4, 1, o) != 1 || fread(&ls[s], 8, 1, o) != 1 || fread(tabs[s].data(), sizeof(cf), Nt, o) != (size_t)Nt) { fprintf(stderr, "bad dump %s\n", fn); return 10; } fclose(o); L = std::max(L, ls[s]); }
            cur = bufA.data(); nxt = bufB.data(); double norm = 0;
            for (int64_t c = 0; c < Nt; c++) { cd acc = 0; for (int s = 0; s < nsl; s++) acc += cd(tabs[s][c]) * std::exp2(ls[s] - L); cur[c] = cf(acc); norm += std::norm(acc); }
            double lp_prefix = -K + std::log2(norm) + 2 * L; double prev = norm, logscale = L;
            std::mt19937_64 r2(seed + (unsigned long long)li); double lps = run_steps(K, ns, fb.data(), bits.data(), cond.data(), &r2, &prev, &logscale);
            for (int k = 0; k < ns; k++) putchar('0' + (k < K ? fb[k] : bits[k])); printf(" %.12f %.12f", lp_prefix + lps, lp_prefix);
            for (int k = K; k < ns; k++) printf(" %.9f", cond[k]); putchar('\n'); fflush(stdout); li++; }
    } else if (!strcmp(argv[2], "prob") || !strcmp(argv[2], "chain")) {
        const bool ch = !strcmp(argv[2], "chain"); std::vector<double> cond(ns);
        FILE* g = fopen(argv[3], "r"); char line[4096];
        while (fgets(line, sizeof line, g)) {
            if ((int)strlen(line) < ns) continue; std::vector<int> fb(ns);
            for (int k = 0; k < ns; k++) fb[k] = line[k] - '0';
            if (ch) { printf("%.12f", chain1(fb.data(), cond.data())); for (int k = 0; k < ns; k++) printf(" %.9f", cond[k]); putchar('\n'); }
            else printf("%.12f\n", chain(fb.data(), bits.data(), nullptr));
            fflush(stdout);
        }
    }
    return 0;
}
