// Case (1) of minweight_exact.py: enumerate c = Av with popcount(c) <= K and find the lightest stabilizer S_v with
// Av = c, where A may be singular: v ranges over the coset G c + ker A (G a generalized inverse, ker A given), and c
// must lie in the column space (annihilated by every kernel vector, A symmetric). weight = popcount(v | c).
// Usage: ./minweight_case1 G.txt ker.txt K      (G: n lines of n chars; ker.txt: one kernel vector per line)
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <omp.h>
typedef unsigned long long u64;
struct B2 { u64 a, b; };
static inline int pc(B2 x) { return __builtin_popcountll(x.a) + __builtin_popcountll(x.b); }
static inline B2 x2(B2 x, B2 y) { return {x.a ^ y.a, x.b ^ y.b}; }
static inline B2 o2(B2 x, B2 y) { return {x.a | y.a, x.b | y.b}; }
int n, nk; std::vector<B2> Gcol, E, coset; std::vector<std::vector<unsigned char>> Kbit;   // Kbit[j][i] = kernel_j[i]
struct Res { int w; std::vector<int> c; B2 v; };
static std::vector<B2> read(const char* p) {
    FILE* f = fopen(p, "r"); char line[1024]; std::vector<std::string> rows; std::vector<B2> out;
    while (fgets(line, sizeof line, f)) { std::string s(line); while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back(); if (!s.empty()) rows.push_back(s); }
    for (auto& s : rows) { B2 b{0, 0}; for (int j = 0; j < (int)s.size(); j++) if (s[j] == '1') { if (j < 64) b.a |= 1ULL << j; else b.b |= 1ULL << (j - 64); } out.push_back(b); }
    return out;
}
static void rec(int k, int start, B2 c, B2 v0, unsigned par, int* idx, int depth, Res& best) {
    if (k == 0) {
        if (par) return;                                   // c not in the column space
        for (const B2& kv : coset) { B2 v = x2(v0, kv); int w = pc(o2(v, c)); if (w < best.w) { best.w = w; best.c.assign(idx, idx + depth); best.v = v; } }
        return;
    }
    for (int j = start; j <= n - k; j++) {
        idx[depth] = j; unsigned p2 = par; for (int t = 0; t < nk; t++) if (Kbit[t][j]) p2 ^= 1u << t;
        rec(k - 1, j + 1, x2(c, E[j]), x2(v0, Gcol[j]), p2, idx, depth + 1, best);
    }
}
int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: minweight_case1 G.txt ker.txt K\n"); return 1; }
    std::vector<B2> Grow = read(argv[1]); std::vector<B2> ker = read(argv[2]); n = (int)Grow.size(); nk = (int)ker.size();
    Gcol.assign(n, {0, 0}); E.assign(n, {0, 0});
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { bool g = j < 64 ? (Grow[i].a >> j) & 1 : (Grow[i].b >> (j - 64)) & 1; if (g) { if (i < 64) Gcol[j].a |= 1ULL << i; else Gcol[j].b |= 1ULL << (i - 64); } }
    for (int i = 0; i < n; i++) { if (i < 64) E[i].a |= 1ULL << i; else E[i].b |= 1ULL << (i - 64); }
    Kbit.assign(nk, std::vector<unsigned char>(n, 0));
    for (int t = 0; t < nk; t++) for (int j = 0; j < n; j++) Kbit[t][j] = j < 64 ? (ker[t].a >> j) & 1 : (ker[t].b >> (j - 64)) & 1;
    for (unsigned m = 0; m < (1u << nk); m++) { B2 s{0, 0}; for (int t = 0; t < nk; t++) if ((m >> t) & 1) s = x2(s, ker[t]); coset.push_back(s); }
    int K = atoi(argv[3]); int global = 1 << 30; Res gbest{1 << 30, {}, {0, 0}};
    printf("n = %d, kernel dimension %d (coset size %zu)\n", n, nk, coset.size());
    for (int k = 1; k <= K; k++) {
        double t0 = omp_get_wtime(); Res best{1 << 30, {}, {0, 0}};
        #pragma omp parallel
        { Res loc{1 << 30, {}, {0, 0}}; int idx[64];
          #pragma omp for schedule(dynamic, 1)
          for (int j = 0; j <= n - k; j++) { idx[0] = j; unsigned p = 0; for (int t = 0; t < nk; t++) if (Kbit[t][j]) p ^= 1u << t; rec(k - 1, j + 1, E[j], Gcol[j], p, idx, 1, loc); }
          #pragma omp critical
          { if (loc.w < best.w) best = loc; } }
        printf("|c| = %d: min weight %d at c =", k, best.w); for (int j : best.c) printf(" %d", j);
        printf("  v ="); for (int j = 0; j < n; j++) if (j < 64 ? (best.v.a >> j) & 1 : (best.v.b >> (j - 64)) & 1) printf(" %d", j);
        printf("   (%.1f s)\n", omp_get_wtime() - t0); fflush(stdout);
        if (best.w < global) { global = best.w; gbest = best; }
    }
    printf("case (1): minimum over 1 <= popcount(Av) <= %d is weight %d\n", K, global);
    return 0;
}
