// Exhaustive minimum stabilizer weight of a graph state, over all v with popcount(v) <= K.
// weight(S_v) = popcount(v | Av) with A the adjacency over GF(2). Since weight >= popcount(v), the global
// minimum d_min equals the minimum over popcount(v) <= K whenever that minimum is <= K + 1.
// Usage: ./minweight A.txt K       (A.txt: n lines of n characters '0'/'1', symmetric, zero diagonal)
// Prints, per popcount k, the minimum weight and one minimiser; parallel over the first chosen index.
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
int n; std::vector<B2> R, E;
struct Res { int w; std::vector<int> v; };
// enumerate combinations of size k from [start, n) extending (v, Av); track best
static void rec(int k, int start, B2 v, B2 Av, int* idx, int depth, Res& best) {
    if (k == 0) { int w = pc(o2(v, Av)); if (w < best.w) { best.w = w; best.v.assign(idx, idx + depth); } return; }
    for (int j = start; j <= n - k; j++) { idx[depth] = j; rec(k - 1, j + 1, x2(v, E[j]), x2(Av, R[j]), idx, depth + 1, best); }
}
int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: minweight A.txt K\n"); return 1; }
    FILE* f = fopen(argv[1], "r"); char line[1024]; std::vector<std::string> rows;
    while (fgets(line, sizeof line, f)) { std::string s(line); while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back(); if (!s.empty()) rows.push_back(s); }
    n = (int)rows.size(); if (n > 128) { fprintf(stderr, "n > 128\n"); return 1; }
    R.assign(n, {0, 0}); E.assign(n, {0, 0});
    for (int i = 0; i < n; i++) { if ((int)rows[i].size() != n) { fprintf(stderr, "bad row %d\n", i); return 1; }
        for (int j = 0; j < n; j++) if (rows[i][j] == '1') { if (j < 64) R[i].a |= 1ULL << j; else R[i].b |= 1ULL << (j - 64); }
        if (i < 64) E[i].a |= 1ULL << i; else E[i].b |= 1ULL << (i - 64); }
    bool nodiag = argc > 3 && !strcmp(argv[3], "nodiag");   // for A^{-1}: symmetric but the diagonal may be nonzero
    for (int i = 0; i < n; i++) { if (!nodiag && rows[i][i] != '0') { fprintf(stderr, "diagonal\n"); return 1; } for (int j = 0; j < n; j++) if (rows[i][j] != rows[j][i]) { fprintf(stderr, "asym\n"); return 1; } }
    int K = atoi(argv[2]); int global = 1 << 30; std::vector<int> gv;
    for (int k = 1; k <= K; k++) {
        double t0 = omp_get_wtime(); Res best{1 << 30, {}};
        #pragma omp parallel
        { Res loc{1 << 30, {}}; int idx[64];
          #pragma omp for schedule(dynamic, 1)
          for (int j = 0; j <= n - k; j++) { idx[0] = j; rec(k - 1, j + 1, E[j], R[j], idx, 1, loc); }
          #pragma omp critical
          { if (loc.w < best.w) best = loc; } }
        printf("|v| = %d: min weight %d at v =", k, best.w); for (int j : best.v) printf(" %d", j);
        printf("   (%.1f s)\n", omp_get_wtime() - t0); fflush(stdout);
        if (best.w < global) { global = best.w; gv = best.v; }
    }
    printf("minimum over popcount(v) <= %d: weight %d, v =", K, global); for (int j : gv) printf(" %d", j); printf("\n");
    if (global <= K + 1) printf("PROVEN: d_min = %d (every v with popcount > %d has weight > %d)\n", global, K, K);
    else printf("NOT PROVEN: d_min <= %d; a lighter stabilizer could still exist with popcount(v) in [%d, %d]\n", global, K + 1, global - 1);
    return 0;
}
