"""Check that the sampler's output is reproduced bit for bit by its random stream where the theory fixes the law.

For the first 22 steps of the production order the cut rank equals the prefix length, so every conditional is
exactly 1/2 and the sampled bit is exactly (U < 0.5) for the uniform draw U of std::mt19937_64 seeded with
1000000 + i (rngbits.cpp reproduces those draws). This script checks that identity for every sample, reports the
per-step agreement beyond step 22 (where conditionals move away from 1/2), and gives the chi-square of the
single-bit means over step ranges. All single-bit marginals of the ideal distribution are exactly uniform
(minweight.py: no stabilizer of weight < 8), so sum z^2 over k bits is chi-square with k degrees of freedom.
Usage: python -P rng_check.py results/samples_nq70.txt work/rngbits_1M.txt results/laptop_samples_w30.txt work/rngbits_2M.txt"""
import sys, numpy as np
from scipy.stats import chi2
def samples(p): return [l.split() for l in open(p) if len(l.split()) >= 3 and len(l.split()[0]) == 70]
def rng(p): return {int(l.split()[1]): l.split()[0] for l in open(p) if l.strip()}
S, R = samples(sys.argv[1]), rng(sys.argv[2])
m22 = sum(R[int(s[2])][:22] == s[0][:22] for s in S); mall = sum(R[int(s[2])] == s[0] for s in S)
agree = np.mean([[R[int(s[2])][k] == s[0][k] for k in range(70)] for s in S], 0)
print(f"{len(S)} production samples: first 22 bits equal the RNG sign bits in {m22}/{len(S)}; whole 70-bit string equal in {mall}")
print("per-step agreement with (U < 0.5):", " ".join(f"{a:.2f}" for a in agree))
X = 1 - 2 * np.array([[int(c) for c in s[0]] for s in S], np.int8); z = X.mean(0) * np.sqrt(len(X))
for a, b in [(0, 22), (22, 70), (0, 70)]:
    s = float((z[a:b] ** 2).sum()); print(f"single-bit means, steps {a + 1}-{b}: sum z^2 = {s:.1f} on {b - a} dof, p = {chi2.sf(s, b - a):.3f}")
if len(sys.argv) > 4:
    L, R2 = samples(sys.argv[3]), rng(sys.argv[4])
    print(f"laptop samples: first 22 bits equal the RNG sign bits in {sum(R2[int(s[2])][:22] == s[0][:22] for s in L)}/{len(L)}")
# control: a stream with the wrong seed offset must not match
bad = sum(R.get(int(s[2]) + 1, "x" * 70)[:22] == s[0][:22] for s in S)
print(f"control (seed + 1): first 22 bits match in {bad}/{len(S)}")
