// CUDA port of the in-place cut-rank sampler gsamp2.cpp: one complex64 arena of 2^R entries in device memory.
// Same plan format and the same modes (sample, prob, chain, hybrid, prefix, resume); one thread per coset
// representative c, which reads its own coset's slots and then overwrites its own slots (as the CPU version does),
// so the in-place butterfly is race-free across threads. Accumulations are in double.
// Build: nvcc -O3 -arch=native -o gsamp_cuda gsamp_cuda.cu      Usage: as gsamp2 (see that file's header).
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <complex>
#include <vector>
#include <random>
#include <algorithm>
#include <cuda_runtime.h>
typedef std::complex<double> cd;
struct Term { uint64_t tq; int x, tau, par; };
struct StepH { int v, nk, m; uint64_t g, rhoc; std::vector<uint64_t> tk; std::vector<Term> t; cd M[2][2]; };
// device-side step description (plain data)
struct StepD { int nk, m, nt; uint64_t g, rhoc; uint64_t tq[16]; int x[16], tau[16], par[16]; double Mre[2][2], Mim[2][2]; uint64_t tab[5][256]; };
static std::vector<StepH> st; static int n, ns, R; static float2* d_mem = nullptr; static StepD* d_steps = nullptr;
static double *d_acc = nullptr; static const int BLOCK = 256; static int GRID = 2048;
#define CK(x) do { cudaError_t e = (x); if (e != cudaSuccess) { fprintf(stderr, "CUDA error %s at line %d\n", cudaGetErrorString(e), __LINE__); exit(99); } } while (0)

__device__ __forceinline__ uint64_t base_of_d(const StepD& s, uint64_t c) {
    return s.tab[0][c & 255] ^ s.tab[1][(c >> 8) & 255] ^ s.tab[2][(c >> 16) & 255] ^ s.tab[3][(c >> 24) & 255] ^ s.tab[4][(c >> 32) & 255];
}
// o[b][tau] for coset c; reads the coset's slots
__device__ __forceinline__ void outputs_d(const StepD& s, const float2* mem, uint64_t c, uint64_t base, double ore[2][2], double oim[2][2]) {
    ore[0][0] = ore[0][1] = ore[1][0] = ore[1][1] = 0; oim[0][0] = oim[0][1] = oim[1][0] = oim[1][1] = 0;
    const int pc = __popcll(c & s.rhoc) & 1;
    for (int i = 0; i < s.nt; i++) {
        float2 v = mem[base ^ s.tq[i]]; double vr = v.x, vi = v.y;
        if (s.x[i] && (pc ^ s.par[i])) { vr = -vr; vi = -vi; }
        const int xi = s.x[i], ti = s.tau[i];
        // out[b][tau] += M[b][x] * v
        for (int b = 0; b < 2; b++) { double mr = s.Mre[b][xi], mi = s.Mim[b][xi]; ore[b][ti] += mr * vr - mi * vi; oim[b][ti] += mr * vi + mi * vr; }
    }
}
// pass 1: norms of both branches (acc[0] = n0, acc[1] = n1)
__global__ void k_norms(const StepD* sp, const float2* mem, int64_t N, double* acc) {
    const StepD& s = *sp; double n0 = 0, n1 = 0; double ore[2][2], oim[2][2];
    for (int64_t c = blockIdx.x * (int64_t)blockDim.x + threadIdx.x; c < N; c += (int64_t)gridDim.x * blockDim.x) {
        uint64_t b = base_of_d(s, c); outputs_d(s, mem, c, b, ore, oim);
        n0 += ore[0][0] * ore[0][0] + oim[0][0] * oim[0][0] + ore[0][1] * ore[0][1] + oim[0][1] * oim[0][1];
        n1 += ore[1][0] * ore[1][0] + oim[1][0] * oim[1][0] + ore[1][1] * ore[1][1] + oim[1][1] * oim[1][1];
    }
    __shared__ double s0[BLOCK], s1[BLOCK]; s0[threadIdx.x] = n0; s1[threadIdx.x] = n1; __syncthreads();
    for (int o = BLOCK / 2; o > 0; o >>= 1) { if (threadIdx.x < o) { s0[threadIdx.x] += s0[threadIdx.x + o]; s1[threadIdx.x] += s1[threadIdx.x + o]; } __syncthreads(); }
    if (threadIdx.x == 0) { atomicAdd(&acc[0], s0[0]); atomicAdd(&acc[1], s1[0]); }
}
// pass 2: write branch bv scaled by sc into the coset's own slots; acc[2] += norm of what was written
__global__ void k_write(const StepD* sp, float2* mem, int64_t N, int bv, double sc, double* acc) {
    const StepD& s = *sp; double nb = 0; double ore[2][2], oim[2][2];
    for (int64_t c = blockIdx.x * (int64_t)blockDim.x + threadIdx.x; c < N; c += (int64_t)gridDim.x * blockDim.x) {
        uint64_t b = base_of_d(s, c); outputs_d(s, mem, c, b, ore, oim);
        double v0r = ore[bv][0] * sc, v0i = oim[bv][0] * sc, v1r = ore[bv][1] * sc, v1i = oim[bv][1] * sc;
        mem[b] = make_float2((float)v0r, (float)v0i); if (s.m == 2) mem[b ^ s.g] = make_float2((float)v1r, (float)v1i);
        nb += v0r * v0r + v0i * v0i + v1r * v1r + v1i * v1i;
    }
    __shared__ double sb[BLOCK]; sb[threadIdx.x] = nb; __syncthreads();
    for (int o = BLOCK / 2; o > 0; o >>= 1) { if (threadIdx.x < o) sb[threadIdx.x] += sb[threadIdx.x + o]; __syncthreads(); }
    if (threadIdx.x == 0) atomicAdd(&acc[2], sb[0]);
}
// gather / scatter the live table at step K: live = { base_of(st[K], c) ^ tq_j }, index c * ntq + j
__global__ void k_gather(const StepD* sp, const float2* mem, int64_t Nc, const uint64_t* tqs, int ntq, float2* out) {
    const StepD& s = *sp;
    for (int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x; i < Nc * ntq; i += (int64_t)gridDim.x * blockDim.x) { int64_t c = i / ntq; int j = (int)(i % ntq); out[i] = mem[base_of_d(s, c) ^ tqs[j]]; }
}
__global__ void k_scatter(const StepD* sp, float2* mem, int64_t Nc, const uint64_t* tqs, int ntq, const float2* in) {
    const StepD& s = *sp;
    for (int64_t i = blockIdx.x * (int64_t)blockDim.x + threadIdx.x; i < Nc * ntq; i += (int64_t)gridDim.x * blockDim.x) { int64_t c = i / ntq; int j = (int)(i % ntq); mem[base_of_d(s, c) ^ tqs[j]] = in[i]; }
}
static void norms(int k, double& n0, double& n1) {
    double z[3] = {0, 0, 0}; CK(cudaMemcpy(d_acc, z, sizeof z, cudaMemcpyHostToDevice));
    k_norms<<<GRID, BLOCK>>>(d_steps + k, d_mem, (int64_t)1 << st[k].nk, d_acc); CK(cudaGetLastError());
    CK(cudaMemcpy(z, d_acc, sizeof z, cudaMemcpyDeviceToHost)); n0 = z[0]; n1 = z[1];
}
static double write(int k, int bv, double sc) {
    double z[3] = {0, 0, 0}; CK(cudaMemcpy(d_acc, z, sizeof z, cudaMemcpyHostToDevice));
    k_write<<<GRID, BLOCK>>>(d_steps + k, d_mem, (int64_t)1 << st[k].nk, bv, sc, d_acc); CK(cudaGetLastError());
    CK(cudaMemcpy(z, d_acc, sizeof z, cudaMemcpyDeviceToHost)); return z[2];
}
static void set_one() { float2 one = make_float2(1.f, 0.f); CK(cudaMemcpy(d_mem, &one, sizeof one, cudaMemcpyHostToDevice)); }
// --- the same algorithms as gsamp2.cpp, with the two passes replaced by kernels ---
static double chain(const int* forced, int* bits, std::mt19937_64* rng) {     // normalised every step; two passes per step
    set_one(); double logp = 0;
    for (int k = 0; k < ns; k++) {
        double n0, n1; norms(k, n0, n1); double p1 = n1 / (n0 + n1); int bv;
        if (forced) bv = forced[k]; else { std::uniform_real_distribution<double> U(0, 1); bv = U(*rng) < p1; }
        double pb = bv ? p1 : 1 - p1; if (pb <= 0) return -INFINITY; logp += std::log2(pb); bits[k] = bv;
        write(k, bv, 1.0 / std::sqrt(bv ? n1 : n0));
    }
    return logp;
}
static double chain1(const int* forced, double* cond) {                       // single pass per step, unnormalised (as gsamp2)
    set_one(); double prev = 1;
    for (int k = 0; k < ns; k++) {
        double nb = write(k, forced[k], 1.0); cond[k] = nb / (2 * prev); prev = nb;
        if (nb <= 0) { for (; k < ns; k++) cond[k] = 0; return -INFINITY; }
    }
    return std::log2(prev) - ns;
}
static double run_steps(int k0, int k1, const int* forced, int* bits, double* cond, std::mt19937_64* rng, double* prev, double* logscale) {
    double logp = 0;
    for (int k = k0; k < k1; k++) {
        int bv = forced[k]; const double sc = 1.0 / std::sqrt(*prev);
        if (bv < 0) { double n0, n1; norms(k, n0, n1); std::uniform_real_distribution<double> U(0, 1); bv = U(*rng) < n1 / (n0 + n1); }
        double nb = write(k, bv, sc);
        *logscale += 0.5 * std::log2(*prev); cond[k] = nb / 2; bits[k] = bv; *prev = nb;
        if (nb <= 0) { for (; k < k1; k++) { cond[k] = 0; bits[k] = forced[k] < 0 ? 0 : forced[k]; } return -INFINITY; }
        logp += std::log2(cond[k]);
    }
    return logp;
}
static double hybrid(const int* forced, int* bits, double* cond, std::mt19937_64* rng) {
    set_one(); double prev = 1, logscale = 0; return run_steps(0, ns, forced, bits, cond, rng, &prev, &logscale);
}
static std::vector<uint64_t> coset_offsets(int K) {
    std::vector<uint64_t> v; if (K >= ns) { v.push_back(0); return v; }
    for (const Term& t : st[K].t) v.push_back(t.tq); std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end()); return v;
}
int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: gsamp_cuda plan mode ...\n"); return 1; }
    FILE* f = fopen(argv[1], "r"); if (!f || fscanf(f, "%d %d %d", &n, &ns, &R) != 3) return 2; st.resize(ns);
    std::vector<StepD> hd(ns);
    for (int k = 0; k < ns; k++) { StepH& s = st[k]; unsigned long long g, rc; int nt;
        if (fscanf(f, "%d %d %d %llu %llu %d", &s.v, &s.nk, &s.m, &g, &rc, &nt) != 6) return 3; s.g = g; s.rhoc = rc;
        s.tk.resize(s.nk); for (auto& x : s.tk) { unsigned long long y; if (fscanf(f, "%llu", &y) != 1) return 4; x = y; }
        s.t.resize(nt); for (auto& t : s.t) { unsigned long long y; if (fscanf(f, "%llu %d %d %d", &y, &t.x, &t.tau, &t.par) != 4) return 5; t.tq = y; }
        for (int b = 0; b < 2; b++) for (int x = 0; x < 2; x++) { double re, im; if (fscanf(f, "%lf %lf", &re, &im) != 2) return 6; s.M[b][x] = cd(re, im); }
        if (s.nk > 40 || nt > 16) return 7;
        StepD& d = hd[k]; memset(&d, 0, sizeof d); d.nk = s.nk; d.m = s.m; d.nt = nt; d.g = s.g; d.rhoc = s.rhoc;
        for (int i = 0; i < nt; i++) { d.tq[i] = s.t[i].tq; d.x[i] = s.t[i].x; d.tau[i] = s.t[i].tau; d.par[i] = s.t[i].par; }
        for (int b = 0; b < 2; b++) for (int x = 0; x < 2; x++) { d.Mre[b][x] = s.M[b][x].real(); d.Mim[b][x] = s.M[b][x].imag(); }
        for (int t = 0; t < 5; t++) for (int v = 0; v < 256; v++) { uint64_t a = 0; for (int j = 0; j < 8; j++) if ((v >> j) & 1 && 8 * t + j < s.nk) a ^= s.tk[8 * t + j]; d.tab[t][v] = a; }
    }
    fclose(f);
    int dev = 0; cudaDeviceProp prop; CK(cudaGetDeviceProperties(&prop, dev)); GRID = prop.multiProcessorCount * 8;
    CK(cudaMalloc(&d_mem, ((size_t)1 << R) * sizeof(float2))); CK(cudaMemset(d_mem, 0, ((size_t)1 << R) * sizeof(float2)));
    CK(cudaMalloc(&d_steps, ns * sizeof(StepD))); CK(cudaMemcpy(d_steps, hd.data(), ns * sizeof(StepD), cudaMemcpyHostToDevice)); CK(cudaMalloc(&d_acc, 3 * sizeof(double)));
    fprintf(stderr, "gsamp_cuda: %s, arena 2^%d complex64 = %.2f GiB\n", prop.name, R, ((double)(1ull << R)) * 8 / (1 << 30));
    std::vector<int> bits(ns); std::vector<double> cond(ns);
    if (!strcmp(argv[2], "sample")) {
        long N = atol(argv[3]); unsigned long long seed = strtoull(argv[4], 0, 10);
        for (long i = 0; i < N; i++) { std::mt19937_64 rng(seed + i); double lp = chain(nullptr, bits.data(), &rng);
            for (int k = 0; k < ns; k++) putchar('0' + bits[k]); printf(" %.12f %llu\n", lp, seed + i); fflush(stdout); }
    } else if (!strcmp(argv[2], "prefix")) {
        int K = atoi(argv[4]); const char* outp = argv[5]; FILE* g = fopen(argv[3], "r"); char line[4096]; long li = 0;
        int nk = K < ns ? st[K].nk : 0; std::vector<uint64_t> tqs = coset_offsets(K); int ntq = (int)tqs.size(); int64_t Nc = (int64_t)1 << nk;
        uint64_t* d_tqs; CK(cudaMalloc(&d_tqs, ntq * 8)); CK(cudaMemcpy(d_tqs, tqs.data(), ntq * 8, cudaMemcpyHostToDevice));
        float2* d_buf; CK(cudaMalloc(&d_buf, Nc * ntq * sizeof(float2))); std::vector<float2> buf(Nc * ntq);
        while (fgets(line, sizeof line, g)) { if ((int)strlen(line) < K) continue; std::vector<int> fb(ns, 0);
            for (int k = 0; k < K; k++) { if (line[k] == '?') { fprintf(stderr, "'?' before K\n"); return 8; } fb[k] = line[k] - '0'; }
            set_one(); double prev = 1, logscale = 0; double lp = run_steps(0, K, fb.data(), bits.data(), cond.data(), nullptr, &prev, &logscale);
            if (K < ns) { k_gather<<<GRID, BLOCK>>>(d_steps + K, d_mem, Nc, d_tqs, ntq, d_buf); CK(cudaGetLastError()); CK(cudaMemcpy(buf.data(), d_buf, Nc * ntq * sizeof(float2), cudaMemcpyDeviceToHost)); }
            else CK(cudaMemcpy(buf.data(), d_mem, sizeof(float2), cudaMemcpyDeviceToHost));
            char fn[1024]; snprintf(fn, sizeof fn, "%s_%ld.bin", outp, li++); FILE* o = fopen(fn, "wb");
            fwrite(&nk, 4, 1, o); fwrite(&ntq, 4, 1, o); fwrite(&logscale, 8, 1, o); fwrite(buf.data(), sizeof(float2), Nc * ntq, o); fclose(o);
            printf("%ld %.12f %.6f\n", li - 1, lp, logscale); fflush(stdout); }
    } else if (!strcmp(argv[2], "resume")) {
        int K = atoi(argv[4]); unsigned long long seed = strtoull(argv[5], 0, 10); int nsl = argc - 6; FILE* g = fopen(argv[3], "r"); char line[4096]; long li = 0;
        int nk = K < ns ? st[K].nk : 0; std::vector<uint64_t> tqs = coset_offsets(K); int ntq = (int)tqs.size(); int64_t Nc = (int64_t)1 << nk; int64_t Nt = Nc * ntq;
        uint64_t* d_tqs; CK(cudaMalloc(&d_tqs, ntq * 8)); CK(cudaMemcpy(d_tqs, tqs.data(), ntq * 8, cudaMemcpyHostToDevice));
        float2* d_buf; CK(cudaMalloc(&d_buf, Nt * sizeof(float2))); std::vector<std::vector<float2>> tabs(nsl, std::vector<float2>(Nt)); std::vector<double> ls(nsl); std::vector<float2> sum(Nt);
        while (fgets(line, sizeof line, g)) { if ((int)strlen(line) < ns) continue; std::vector<int> fb(ns);
            for (int k = 0; k < ns; k++) fb[k] = line[k] == '?' ? -1 : line[k] - '0';
            double L = -1e300;
            for (int s = 0; s < nsl; s++) { char fn[1024]; snprintf(fn, sizeof fn, "%s_%ld.bin", argv[6 + s], li); FILE* o = fopen(fn, "rb"); if (!o) { fprintf(stderr, "missing %s\n", fn); return 9; }
                int nk2, ntq2; if (fread(&nk2, 4, 1, o) != 1 || nk2 != nk || fread(&ntq2, 4, 1, o) != 1 || ntq2 != ntq || fread(&ls[s], 8, 1, o) != 1 || fread(tabs[s].data(), sizeof(float2), Nt, o) != (size_t)Nt) { fprintf(stderr, "bad dump %s\n", fn); return 10; } fclose(o); L = std::max(L, ls[s]); }
            double norm = 0;
            for (int64_t i = 0; i < Nt; i++) { double ar = 0, ai = 0; for (int s = 0; s < nsl; s++) { double f2 = std::exp2(ls[s] - L); ar += tabs[s][i].x * f2; ai += tabs[s][i].y * f2; } sum[i] = make_float2((float)ar, (float)ai); norm += ar * ar + ai * ai; }
            if (K < ns) { CK(cudaMemcpy(d_buf, sum.data(), Nt * sizeof(float2), cudaMemcpyHostToDevice)); k_scatter<<<GRID, BLOCK>>>(d_steps + K, d_mem, Nc, d_tqs, ntq, d_buf); CK(cudaGetLastError()); }
            else CK(cudaMemcpy(d_mem, sum.data(), sizeof(float2), cudaMemcpyHostToDevice));
            double lp_prefix = -K + std::log2(norm) + 2 * L; double prev = norm, logscale = L;
            std::mt19937_64 r2(seed + (unsigned long long)li); double lps = run_steps(K, ns, fb.data(), bits.data(), cond.data(), &r2, &prev, &logscale);
            for (int k = 0; k < ns; k++) putchar('0' + (k < K ? fb[k] : bits[k])); printf(" %.12f %.12f", lp_prefix + lps, lp_prefix);
            for (int k = K; k < ns; k++) printf(" %.9f", cond[k]); putchar('\n'); fflush(stdout); li++; }
    } else if (!strcmp(argv[2], "hybrid")) {
        unsigned long long seed = strtoull(argv[4], 0, 10); FILE* g = fopen(argv[3], "r"); char line[4096]; long li = 0;
        while (fgets(line, sizeof line, g)) { if ((int)strlen(line) < ns) continue; std::vector<int> fb(ns);
            for (int k = 0; k < ns; k++) fb[k] = line[k] == '?' ? -1 : line[k] - '0';
            std::mt19937_64 r2(seed + (unsigned long long)li++); double lp = hybrid(fb.data(), bits.data(), cond.data(), &r2);
            for (int k = 0; k < ns; k++) putchar('0' + bits[k]); printf(" %.12f", lp); for (int k = 0; k < ns; k++) printf(" %.9f", cond[k]); putchar('\n'); fflush(stdout); }
    } else {
        const bool ch = !strcmp(argv[2], "chain"); FILE* g = fopen(argv[3], "r"); char line[4096];
        while (fgets(line, sizeof line, g)) { if ((int)strlen(line) < ns) continue; std::vector<int> fb(ns);
            for (int k = 0; k < ns; k++) fb[k] = line[k] - '0';
            if (ch) { printf("%.12f", chain1(fb.data(), cond.data())); for (int k = 0; k < ns; k++) printf(" %.9f", cond[k]); putchar('\n'); }
            else printf("%.12f\n", chain(fb.data(), bits.data(), nullptr));
            fflush(stdout); }
    }
    return 0;
}
