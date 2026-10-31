// Unit test of gsamp2's base_of with nk up to 40: the 5 byte-tables must reproduce the direct XOR of T columns.
#define main gsamp2_main
#include "gsamp2.cpp"
#undef main
#include <random>
int main() {
    std::mt19937_64 rng(5); int bad = 0;
    for (int trial = 0; trial < 200; trial++) {
        Step s; s.nk = 1 + rng() % 40; s.tk.resize(s.nk);
        for (auto& x : s.tk) x = rng() & ((1ULL << 40) - 1);
        for (int t = 0; t < 5; t++) for (int v = 0; v < 256; v++) { uint64_t a = 0;
            for (int j = 0; j < 8; j++) if ((v >> j) & 1 && 8 * t + j < s.nk) a ^= s.tk[8 * t + j]; s.tab[t][v] = a; }
        for (int k = 0; k < 50; k++) { uint64_t c = rng() & ((1ULL << s.nk) - 1); uint64_t ref = 0;
            for (int j = 0; j < s.nk; j++) if ((c >> j) & 1) ref ^= s.tk[j];
            if (base_of(s, c) != ref) bad++; }
    }
    // negative control: drop the 5th table and nk = 36 must fail
    Step s; s.nk = 36; s.tk.resize(36); for (auto& x : s.tk) x = rng() & ((1ULL << 40) - 1);
    for (int t = 0; t < 5; t++) for (int v = 0; v < 256; v++) { uint64_t a = 0;
        for (int j = 0; j < 8; j++) if ((v >> j) & 1 && 8 * t + j < s.nk && t < 4) a ^= s.tk[8 * t + j]; s.tab[t][v] = a; }
    int cbad = 0; for (int k = 0; k < 50; k++) { uint64_t c = rng() & ((1ULL << 36) - 1); uint64_t ref = 0;
        for (int j = 0; j < 36; j++) if ((c >> j) & 1) ref ^= s.tk[j]; if (base_of(s, c) != ref) cbad++; }
    printf("base_of nk<=40: %d mismatches of 10000; control (4 tables only, nk=36): %d mismatches of 50\n", bad, cbad);
    return (bad == 0 && cbad > 40) ? 0 : 1;
}
