// Enumerate the 2^m stabilizers supported on a K-qubit prefix (K <= 64) and compute
//   E[z^2] = sum over supports U of 3^{-|U|} (sum_{S: supp S = U} v(S))^2,   v(S) = sign(S) prod_{j in U} c_j[S_j],
// with c_j[.] = +-1 the signs of the measurement Bloch components (|n_j| = 1/sqrt 3). Input from capvar.py export:
//   K m / K lines "cX cY cZ" / m lines "xmask zmask phase"   (phase p: S = i^p X^x Z^z, p in 0..3)
// Elements are enumerated by a Gray code over the generators with the phase tracked by the symplectic rule
//   (i^p X^a Z^b)(i^q X^c Z^d) = i^{p+q+2|b&c|} X^{a^c} Z^{b^d};  Hermitian sign = i^{p - |x&z|} = +-1.
// Supports are grouped by hashing; to bound memory the enumeration is repeated for each of 2^PB partitions of the
// support mask's top bits (a support is handled by exactly one pass). Threads split the top generators.
// Output: E[z^2], Var(z), number of distinct supports, largest group, and the "diagonal" sum_{S != I} 3^{-|S|}.
// Usage: capvar <gens.txt> [passbits=4] [threads]
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <thread>
#include <mutex>
#include <algorithm>
using u64 = unsigned long long;
struct Gen { u64 x, z; int p; };
int main(int argc, char** argv) {
    FILE* f = fopen(argv[1], "r"); int K, m; if (fscanf(f, "%d %d", &K, &m) != 2) return 1;
    std::vector<int> cX(K), cY(K), cZ(K);
    for (int j = 0; j < K; j++) if (fscanf(f, "%d %d %d", &cX[j], &cY[j], &cZ[j]) != 3) return 1;
    std::vector<Gen> g(m);
    for (int i = 0; i < m; i++) if (fscanf(f, "%llu %llu %d", &g[i].x, &g[i].z, &g[i].p) != 3) return 1;
    fclose(f);
    int PB = argc > 2 ? atoi(argv[2]) : 4; int nthr = argc > 3 ? atoi(argv[3]) : (int)std::thread::hardware_concurrency();
    int TB = 0; while ((1 << TB) < nthr && TB < m) TB++;            // top TB generators split across threads
    // sign table per qubit per letter: letter code 1 = X (x=1,z=0), 2 = Z (0,1), 3 = Y (1,1)
    auto vsign = [&](u64 x, u64 z, int p) -> int {
        u64 u = x | z; int cnt = 0; int prod = 1;
        for (u64 t = u; t; t &= t - 1) { int j = __builtin_ctzll(t); u64 b = 1ULL << j; int xx = (x & b) != 0, zz = (z & b) != 0;
            prod *= (xx && zz) ? cY[j] : (xx ? cX[j] : cZ[j]); cnt += xx && zz; }
        int q = ((p - cnt) % 4 + 4) % 4; if (q != 0 && q != 2) { fprintf(stderr, "non-Hermitian element\n"); exit(2); }
        return (q == 0 ? 1 : -1) * prod;
    };
    long double Ez2 = 0, diag = 0; u64 nsupp = 0, maxgroup = 0; std::mutex mtx;
    u64 passes = 1ULL << PB;
    for (u64 pass = 0; pass < passes; pass++) {
        // each thread: fixed top bits (thread id in 0..2^TB-1), Gray code over the lower m-TB generators
        std::vector<std::thread> th; std::vector<std::unordered_map<u64, long long>> maps(1 << TB); std::vector<u64> wsum(1 << TB, 0);
        for (int t = 0; t < (1 << TB); t++) th.emplace_back([&, t]() {
            if (t >= (1 << std::min(TB, m))) return;
            u64 x = 0, z = 0; int p = 0;
            for (int i = 0; i < TB && i < m; i++) if ((t >> i) & 1) { int gi = m - 1 - i; p = (p + g[gi].p + 2 * __builtin_popcountll(z & g[gi].x)) % 4; x ^= g[gi].x; z ^= g[gi].z; }
            int low = std::max(m - TB, 0); u64 N = 1ULL << low; auto& mp = maps[t];
            auto visit = [&](u64 xx, u64 zz, int pp) {
                u64 u = xx | zz; if (PB && (u >> (K - PB)) != pass) return;
                int v = vsign(xx, zz, pp); mp[u] += v;
            };
            visit(x, z, p);
            for (u64 i = 1; i < N; i++) { int gi = __builtin_ctzll(i); p = (p + g[gi].p + 2 * __builtin_popcountll(z & g[gi].x)) % 4; x ^= g[gi].x; z ^= g[gi].z; visit(x, z, p); }
        });
        for (auto& t : th) t.join();
        // merge thread maps (same support may appear in several threads)
        std::unordered_map<u64, long long> all;
        for (auto& mp : maps) { for (auto& kv : mp) all[kv.first] += kv.second; std::unordered_map<u64, long long>().swap(mp); }   // merge, then free
        for (auto& kv : all) { int w = __builtin_popcountll(kv.first); long double term = powl(3.0L, -w) * (long double)kv.second * kv.second;
            Ez2 += term; nsupp++; u64 gsz = (u64)llabs(kv.second); if (gsz > maxgroup) maxgroup = gsz; }
        fprintf(stderr, "pass %llu/%llu: %zu supports, running E[z^2] = %.9Lf\n", pass + 1, passes, all.size(), Ez2);
    }
    // diagonal part sum_{S != I} 3^{-|S|} (every element once) for comparison: recompute quickly in one pass without grouping
    {
        u64 x = 0, z = 0; u64 N = 1ULL << m;
        for (u64 i = 1; i < N; i++) { int gi = __builtin_ctzll(i); x ^= g[gi].x; z ^= g[gi].z; diag += powl(3.0L, -__builtin_popcountll(x | z)); }
    }
    printf("K = %d, group size 2^%d = %llu, distinct supports = %llu, largest |group sum| = %llu\n", K, m, 1ULL << m, nsupp, maxgroup);
    printf("E[z^2] = %.12Lf   Var(z) = %.6Le   (diagonal sum_{S!=I} 3^-|S| = %.6Le)\n", Ez2, Ez2 - 1, diag);
    for (double c : {1.01, 1.02, 1.03, 1.05, 1.10}) printf("Chebyshev: Pr_uniform(z > %.2f) <= %.3Le\n", c, (Ez2 - 1) / ((c - 1) * (c - 1)));
    return 0;
}
