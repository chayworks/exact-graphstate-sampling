// Reproduce the sampler's random stream: for seeds base..base+N-1, std::mt19937_64 rng(seed) and
// std::uniform_real_distribution<double>(0,1) drawn 70 times; print the bit (U < 0.5) per draw.
// Mode "dump base N": one line per seed with 70 bits (to compare with the first 22 bits of the samples, whose
// conditionals are exactly 1/2 so the sampled bit is exactly U < 0.5).
// Mode "scan N reps": for reps bases spaced far apart, the chi2 statistic sum_k z_k^2 over the 70 draws
// (z_k = mean of +-1 bits * sqrt(N)); its mean should be 70 if consecutive seeds give independent streams.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>
#include <cmath>
int main(int argc, char** argv) {
    if (!strcmp(argv[1], "dump")) {
        unsigned long long base = strtoull(argv[2], 0, 10); int N = atoi(argv[3]);
        for (int i = 0; i < N; i++) { std::mt19937_64 rng(base + i); std::uniform_real_distribution<double> U(0, 1);
            for (int k = 0; k < 70; k++) putchar(U(rng) < 0.5 ? '1' : '0'); printf(" %llu\n", base + i); }
    } else {
        int N = atoi(argv[2]), reps = atoi(argv[3]); double tot = 0, tot2 = 0;
        for (int r = 0; r < reps; r++) { unsigned long long base = 1000000ULL + 7919ULL * 1000 * r; std::vector<double> s(70, 0);
            for (int i = 0; i < N; i++) { std::mt19937_64 rng(base + i); std::uniform_real_distribution<double> U(0, 1);
                for (int k = 0; k < 70; k++) s[k] += (U(rng) < 0.5) ? -1 : 1; }
            double chi = 0; for (int k = 0; k < 70; k++) { double z = s[k] / sqrt((double)N); chi += z * z; }
            tot += chi; tot2 += chi * chi; }
        double m = tot / reps, sd = sqrt(tot2 / reps - m * m);
        printf("consecutive seeds, N=%d, %d reps: mean chi2 over 70 draws = %.2f (expected 70.0, sd of mean %.2f)\n", N, reps, m, sd / sqrt((double)reps));
    }
    return 0;
}
