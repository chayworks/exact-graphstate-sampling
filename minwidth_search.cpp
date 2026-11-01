// Vertex-order search for the cut-rank dynamic program on large graphs (n <= 1024): simulated annealing with
// local moves and incremental recomputation of the affected cuts. Objective: log2 sum_k 2^{r_k} (dominated by the
// width); reports the best (width, log2 sum) found. Several independent chains run in parallel (OpenMP), each from
// the given start order with its own seed.
// Usage: ./minwidth_search A.txt start_order.json iters seed_base out_order.json [max_move=24] [T0=0.6]
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <omp.h>
typedef unsigned long long u64;
static int n, W;                                   // vertices, words per row
static std::vector<u64> A;                         // n rows of W words
static inline const u64* row(int i) { return &A[(size_t)i * W]; }
// rank over GF(2) of A[P, F] where P = order[0..k), F = order[k..n): build rows restricted to F
static int cut_rank(const std::vector<int>& ord, int k, std::vector<u64>& scratch, const std::vector<u64>& fmask) {
    int rows = k; scratch.assign((size_t)rows * W, 0);
    for (int i = 0; i < rows; i++) { const u64* r = row(ord[i]); for (int w = 0; w < W; w++) scratch[(size_t)i * W + w] = r[w] & fmask[w]; }
    int rk = 0;
    for (int w = 0; w < W && rk < rows; w++) for (int b = 0; b < 64 && rk < rows; b++) {
        u64 bit = 1ULL << b; int p = -1;
        for (int i = rk; i < rows; i++) if (scratch[(size_t)i * W + w] & bit) { p = i; break; }
        if (p < 0) continue;
        if (p != rk) for (int x = 0; x < W; x++) std::swap(scratch[(size_t)p * W + x], scratch[(size_t)rk * W + x]);
        for (int i = 0; i < rows; i++) if (i != rk && (scratch[(size_t)i * W + w] & bit)) for (int x = w; x < W; x++) scratch[(size_t)i * W + x] ^= scratch[(size_t)rk * W + x];
        rk++;
    }
    return rk;
}
// cheaper: rank(A[P,F]) = rank(A[F,P]) ; use the smaller side as rows
static int cut(const std::vector<int>& ord, int k, std::vector<u64>& scratch) {
    std::vector<u64> fmask(W, 0), pmask(W, 0);
    for (int i = k; i < n; i++) fmask[ord[i] >> 6] |= 1ULL << (ord[i] & 63);
    for (int i = 0; i < k; i++) pmask[ord[i] >> 6] |= 1ULL << (ord[i] & 63);
    if (k <= n - k) return cut_rank(ord, k, scratch, fmask);
    std::vector<int> rev(ord.rbegin(), ord.rend()); return cut_rank(rev, n - k, scratch, pmask);
}
static double objective(const std::vector<int>& r) { double m = 0; for (int x : r) m = std::max(m, (double)x); double s = 0; for (int x : r) s += std::exp2(x - m); return m + std::log2(s); }
int main(int argc, char** argv) {
    if (argc < 6) { fprintf(stderr, "usage: minwidth_search A.txt start_order.json iters seed_base out.json [max_move] [T0]\n"); return 1; }
    { FILE* f = fopen(argv[1], "r"); char line[4096]; std::vector<std::string> rows;
      while (fgets(line, sizeof line, f)) { std::string s(line); while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back(); if (!s.empty()) rows.push_back(s); }
      n = (int)rows.size(); W = (n + 63) / 64; A.assign((size_t)n * W, 0);
      for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (rows[i][j] == '1') A[(size_t)i * W + (j >> 6)] |= 1ULL << (j & 63); fclose(f); }
    std::vector<int> start; { FILE* f = fopen(argv[2], "r"); int v; char c; while (fscanf(f, "%c", &c) == 1) { if (c >= '0' && c <= '9') { ungetc(c, f); if (fscanf(f, "%d", &v) == 1) start.push_back(v); } } fclose(f); }
    if ((int)start.size() != n) { fprintf(stderr, "order has %zu entries, graph has %d\n", start.size(), n); return 1; }
    long iters = atol(argv[3]); int seed0 = atoi(argv[4]); int maxmove = argc > 6 ? atoi(argv[6]) : 24; double T0 = argc > 7 ? atof(argv[7]) : 0.6;
    int nth = omp_get_max_threads(); std::vector<std::vector<int>> best_ord(nth); std::vector<double> best_obj(nth, 1e300); std::vector<int> best_w(nth, 1 << 30);
    #pragma omp parallel
    {
        int tid = omp_get_thread_num(); std::mt19937_64 rng(seed0 + tid); std::vector<u64> scratch;
        std::vector<int> ord = start, r(n, 0);
        for (int k = 1; k < n; k++) r[k] = cut(ord, k, scratch);
        double cur = objective(r); best_obj[tid] = cur; best_ord[tid] = ord; best_w[tid] = *std::max_element(r.begin(), r.end());
        std::uniform_real_distribution<double> U(0, 1);
        for (long it = 0; it < iters; it++) {
            double T = T0 * (1.0 - (double)it / iters) + 0.02;
            int i = (int)(U(rng) * n); int d = (int)(U(rng) * (2 * maxmove + 1)) - maxmove; int j = std::min(n - 1, std::max(0, i + d)); if (i == j) continue;
            std::vector<int> nord = ord; int v = nord[i]; nord.erase(nord.begin() + i); nord.insert(nord.begin() + j, v);
            int lo = std::min(i, j), hi = std::max(i, j); std::vector<int> nr = r;
            for (int k = std::max(1, lo); k <= hi && k < n; k++) nr[k] = cut(nord, k, scratch);
            double nobj = objective(nr);
            if (nobj <= cur || U(rng) < std::exp((cur - nobj) / T)) {
                ord.swap(nord); r.swap(nr); cur = nobj;
                int w = *std::max_element(r.begin(), r.end());
                if (cur < best_obj[tid]) { best_obj[tid] = cur; best_ord[tid] = ord; best_w[tid] = w; }
            }
            if (it % 2000 == 0 && tid == 0) { printf("it %ld: current width %d obj %.2f; best %d / %.2f\n", it, *std::max_element(r.begin(), r.end()), cur, best_w[tid], best_obj[tid]); fflush(stdout); }
        }
    }
    int b = 0; for (int t = 1; t < nth; t++) if (best_obj[t] < best_obj[b]) b = t;
    printf("best over %d chains: width %d, log2 sum 2^r = %.2f\n", nth, best_w[b], best_obj[b]);
    FILE* f = fopen(argv[5], "w"); fprintf(f, "["); for (int i = 0; i < n; i++) fprintf(f, "%s%d", i ? ", " : "", best_ord[b][i]); fprintf(f, "]\n"); fclose(f);
    return 0;
}
