// In-place exact sampler (see inplace.py): one complex64 arena of 2^R entries.
// Plan: n nsteps R; per step: v nk m g rhoc nterms; nk Tkeep values; nterms x (Tq x tau par); 8 doubles M[b][x].
// Usage: gsamp2 plan sample N seed | gsamp2 plan prob stringsfile   (same output format as gsamp)
//        gsamp2 plan chain stringsfile   -> per line: log2 Pr, then Pr(b_k | b_1..b_{k-1}) for each step (one pass per step)
//        gsamp2 plan hybrid stringsfile seed -> lines may contain '?' (sampled from the exact conditional; only valid in a
//                                               trailing suffix of data vertices); prints bits, log2 Pr, all conditionals
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>
#include <random>
#include <memory>
#include <algorithm>
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
// Hybrid mode: forced[k] in {0,1} fixes step k's bit (single pass, lazily rescaled table); forced[k] == -1 draws the
// bit from its exact conditional Pr(b_k | b_1..b_{k-1}) (two passes). The conditional at a forced step marginalises
// every later vertex, so sampled positions are only meaningful after all projected (ancilla) vertices, i.e. at the
// end of the order with a suffix of data vertices. cond[k] receives Pr(b_k | prefix) for every step.
// Steps k0..k1-1 of the hybrid loop on the current table of norm *prev (stored table = true table * 2^-*logscale).
// Returns sum of log2 cond over the steps run; updates *prev and *logscale; -INFINITY on a zero-probability branch.
static double run_steps(int k0, int k1, const int* forced, int* bits, double* cond, std::mt19937_64* rng, double* prev, double* logscale) {
    double logp = 0;
    for (int k = k0; k < k1; k++) {
        const Step& s = st[k]; const int64_t N = (int64_t)1 << s.nk; int bv = forced[k]; const double sc = 1.0 / std::sqrt(*prev);
        if (bv < 0) {
            double n0 = 0, n1 = 0;
#pragma omp parallel for schedule(static) reduction(+:n0,n1)
            for (int64_t c = 0; c < N; c++) { cd o[2][2]; uint64_t b = base_of(s, c); outputs(s, c, b, o);
                n0 += std::norm(o[0][0]) + std::norm(o[0][1]); n1 += std::norm(o[1][0]) + std::norm(o[1][1]); }
            std::uniform_real_distribution<double> U(0, 1); bv = U(*rng) < n1 / (n0 + n1);
        }
        double nb = 0;
#pragma omp parallel for schedule(static) reduction(+:nb)
        for (int64_t c = 0; c < N; c++) { cd o[2][2]; uint64_t b = base_of(s, c); outputs(s, c, b, o);
            cd v0 = o[bv][0] * sc, v1 = o[bv][1] * sc; mem[b] = cf(v0); if (s.m == 2) mem[b ^ s.g] = cf(v1);
            nb += std::norm(v0) + std::norm(v1); }
        *logscale += 0.5 * std::log2(*prev); cond[k] = nb / 2; bits[k] = bv; *prev = nb;
        if (nb <= 0) { for (; k < k1; k++) { cond[k] = 0; bits[k] = forced[k] < 0 ? 0 : forced[k]; } return -INFINITY; }
        logp += std::log2(cond[k]);
    }
    return logp;
}
static double hybrid(const int* forced, int* bits, double* cond, std::mt19937_64* rng) {
    mem[0] = 1; double prev = 1, logp = 0;
    for (int k = 0; k < ns; k++) {
        const Step& s = st[k]; const int64_t N = (int64_t)1 << s.nk; int bv = forced[k]; const double sc = 1.0 / std::sqrt(prev);
        if (bv < 0) {
            double n0 = 0, n1 = 0;
#pragma omp parallel for schedule(static) reduction(+:n0,n1)
            for (int64_t c = 0; c < N; c++) { cd o[2][2]; uint64_t b = base_of(s, c); outputs(s, c, b, o);
                n0 += std::norm(o[0][0]) + std::norm(o[0][1]); n1 += std::norm(o[1][0]) + std::norm(o[1][1]); }
            std::uniform_real_distribution<double> U(0, 1); bv = U(*rng) < n1 / (n0 + n1);
        }
        double nb = 0;
#pragma omp parallel for schedule(static) reduction(+:nb)
        for (int64_t c = 0; c < N; c++) { cd o[2][2]; uint64_t b = base_of(s, c); outputs(s, c, b, o);
            cd v0 = o[bv][0] * sc, v1 = o[bv][1] * sc; mem[b] = cf(v0); if (s.m == 2) mem[b ^ s.g] = cf(v1);
            nb += std::norm(v0) + std::norm(v1); }
        cond[k] = nb / 2; bits[k] = bv; prev = nb;      // table now has norm nb = 2 * Pr(b_k | prefix) (previous norm scaled to 1)
        if (nb <= 0) { for (; k < ns; k++) { cond[k] = 0; bits[k] = forced[k] < 0 ? 0 : forced[k]; } return -INFINITY; }
        logp += std::log2(cond[k]);
    }
    return logp;
}
// Distinct coset offsets tq of step K's terms: the old table's live addresses are base_of(st[K], c) ^ tq. {0} past the end.
static std::vector<uint64_t> coset_offsets(int K) {
    std::vector<uint64_t> v; if (K >= ns) { v.push_back(0); return v; }
    for (const Term& t : st[K].t) v.push_back(t.tq); std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end()); return v;
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
    } else if (!strcmp(argv[2], "prefix")) {
        // gsamp2 plan prefix stringsfile K outprefix: run the forced steps 0..K-1 of each line and dump the table
        // (binary: int32 nk, double logscale, 2^nk complex64 at logical indices c) to outprefix_<line>.bin.
        // Slicing: each slice plan shares the step structure, so the dumps of the slices add up (see resume).
        int K = atoi(argv[4]); const char* outp = argv[5]; FILE* g = fopen(argv[3], "r"); char line[4096]; long li = 0; std::vector<double> cond(ns);
        while (fgets(line, sizeof line, g)) { if ((int)strlen(line) < K) continue; std::vector<int> fb(ns, 0);
            for (int k = 0; k < K; k++) { if (line[k] == '?') { fprintf(stderr, "'?' before K in prefix mode\n"); return 8; } fb[k] = line[k] - '0'; }
            mem[0] = 1; double prev = 1, logscale = 0; double lp = run_steps(0, K, fb.data(), bits.data(), cond.data(), nullptr, &prev, &logscale);
            // live table = { base_of(st[K], c) ^ tq : c < 2^nk, tq among the coset offsets of step K's terms }
            int nk = K < ns ? st[K].nk : 0; std::vector<uint64_t> tqs = coset_offsets(K); int ntq = (int)tqs.size();
            char fn[1024]; snprintf(fn, sizeof fn, "%s_%ld.bin", outp, li++); FILE* o = fopen(fn, "wb");
            fwrite(&nk, 4, 1, o); fwrite(&ntq, 4, 1, o); fwrite(&logscale, 8, 1, o);
            for (int64_t c = 0; c < ((int64_t)1 << nk); c++) for (uint64_t tq : tqs) { cf v = K < ns ? mem[base_of(st[K], c) ^ tq] : mem[0]; fwrite(&v, sizeof(cf), 1, o); }
            fclose(o); printf("%ld %.12f %.6f\n", li - 1, lp, logscale); fflush(stdout); }
    } else if (!strcmp(argv[2], "resume")) {
        // gsamp2 plan resume stringsfile K seed dumpprefix1 [dumpprefix2 ...]: sum the dumped tables of all slices for each
        // line, then run steps K..ns-1 with '?' sampled. Prints: bits, log2 Pr(full string), log2 Pr(prefix of K steps),
        // then the conditionals of steps K..ns-1.  Pr(prefix) = 2^-K * || summed true table ||^2.
        int K = atoi(argv[4]); unsigned long long seed = strtoull(argv[5], 0, 10); int nsl = argc - 6; FILE* g = fopen(argv[3], "r"); char line[4096]; long li = 0;
        std::vector<double> cond(ns); int nk = K < ns ? st[K].nk : 0; std::vector<uint64_t> tqs = coset_offsets(K); int ntq = (int)tqs.size();
        int64_t Nt = ((int64_t)1 << nk) * ntq; std::vector<std::vector<cf>> tabs(nsl, std::vector<cf>(Nt)); std::vector<double> ls(nsl);
        while (fgets(line, sizeof line, g)) { if ((int)strlen(line) < ns) continue; std::vector<int> fb(ns);
            for (int k = 0; k < ns; k++) fb[k] = line[k] == '?' ? -1 : line[k] - '0';
            double L = -1e300;
            for (int s = 0; s < nsl; s++) { char fn[1024]; snprintf(fn, sizeof fn, "%s_%ld.bin", argv[6 + s], li); FILE* o = fopen(fn, "rb"); if (!o) { fprintf(stderr, "missing %s\n", fn); return 9; }
                int nk2, ntq2; if (fread(&nk2, 4, 1, o) != 1 || nk2 != nk || fread(&ntq2, 4, 1, o) != 1 || ntq2 != ntq || fread(&ls[s], 8, 1, o) != 1 || fread(tabs[s].data(), sizeof(cf), Nt, o) != (size_t)Nt) { fprintf(stderr, "bad dump %s\n", fn); return 10; } fclose(o); L = std::max(L, ls[s]); }
            double norm = 0; int64_t idx = 0;
            for (int64_t c = 0; c < ((int64_t)1 << nk); c++) for (uint64_t tq : tqs) { cd acc = 0; for (int s = 0; s < nsl; s++) acc += cd(tabs[s][idx]) * std::exp2(ls[s] - L);
                if (K < ns) mem[base_of(st[K], c) ^ tq] = cf(acc); else mem[0] = cf(acc); norm += std::norm(acc); idx++; }
            double lp_prefix = -K + std::log2(norm) + 2 * L; double prev = norm, logscale = L;
            std::mt19937_64 r2(seed + (unsigned long long)li); double lps = run_steps(K, ns, fb.data(), bits.data(), cond.data(), &r2, &prev, &logscale);
            for (int k = 0; k < ns; k++) putchar('0' + (k < K ? fb[k] : bits[k])); printf(" %.12f %.12f", lp_prefix + lps, lp_prefix);
            for (int k = K; k < ns; k++) printf(" %.9f", cond[k]); putchar('\n'); fflush(stdout); li++; }
    } else if (!strcmp(argv[2], "hybrid")) {
        // gsamp2 plan hybrid stringsfile seed: lines of ns chars in {0,1,?}; per line: bits, log2 Pr(string), then every conditional
        unsigned long long seed = strtoull(argv[4], 0, 10); std::vector<double> cond(ns); std::mt19937_64 rng(seed);
        FILE* g = fopen(argv[3], "r"); char line[4096]; long li = 0;
        while (fgets(line, sizeof line, g)) { if ((int)strlen(line) < ns) continue; std::vector<int> fb(ns);
            for (int k = 0; k < ns; k++) fb[k] = line[k] == '?' ? -1 : line[k] - '0';
            std::mt19937_64 r2(seed + (unsigned long long)li++); double lp = hybrid(fb.data(), bits.data(), cond.data(), &r2);
            for (int k = 0; k < ns; k++) putchar('0' + bits[k]); printf(" %.12f", lp); for (int k = 0; k < ns; k++) printf(" %.9f", cond[k]); putchar('\n'); fflush(stdout); }
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
