"""Engine-independent check of the exact samples. While the cut rank of the measured prefix equals its length,
the reduced state of the graph state on that prefix is maximally mixed, so the prefix marginal is exactly uniform.
Hence the first t0 step-bits of the samples (t0 read from the plan) must show no bias and no pairwise or
three-way parity correlations beyond sampling error. This uses only the plan's rank profile and the sample file,
not the sampler or the amplitude engine. A sign or kernel error that correlated early bits would fail it.
Usage: python check_uniform_prefix.py <plan.txt> <samples.txt>"""
import sys, itertools, numpy as np
rows = [l.split() for l in open(sys.argv[1])][1:]
rn = [int(r[2]) for r in rows]; t0 = next(t for t in range(1, len(rn) + 1) if rn[t - 1] < t) - 1
S = np.array([[int(c) for c in l.split()[0]] for l in open(sys.argv[2]) if len(l.split()) > 1 and len(l.split()[0]) == len(rn)])
N = len(S); X = 1 - 2 * S[:, :t0]
z1 = X.mean(0) * np.sqrt(N)
z2 = np.array([(X[:, i] * X[:, j]).mean() for i, j in itertools.combinations(range(t0), 2)]) * np.sqrt(N)
z3 = np.array([(X[:, i] * X[:, j] * X[:, k]).mean() for i, j, k in itertools.combinations(range(t0), 3)]) * np.sqrt(N)
print(f"uniform prefix length t0 = {t0} (cut rank = prefix length); N = {N} samples; z = sqrt(N) * sample mean")
print(f"bit means      : {t0:5d} tests, max |z| = {abs(z1).max():.2f}, frac |z| > 2 = {(abs(z1) > 2).mean():.3f}")
print(f"pair parities  : {len(z2):5d} tests, max |z| = {abs(z2).max():.2f}, frac |z| > 2 = {(abs(z2) > 2).mean():.3f}")
print(f"triple parities: {len(z3):5d} tests, max |z| = {abs(z3).max():.2f}, frac |z| > 2 = {(abs(z3) > 2).mean():.3f}")
print("expected for uniform: frac |z| > 2 about 0.046; max |z| over 4000 tests about 3.5-4")
# negative control: a deliberately biased copy must fail (flip bit 3 of every sample whose bit 7 is 1)
Y = X.copy(); Y[:, 3] = Y[:, 7]
zc = abs((Y[:, 3] * Y[:, 7]).mean()) * np.sqrt(N); print(f"control (bit 3 replaced by a copy of bit 7): pair |z| for (3,7) = {zc:.1f}, must exceed 4.5")
ok = abs(z1).max() < 4 and abs(z2).max() < 4.5 and abs(z3).max() < 5 and (abs(z3) > 2).mean() < 0.08 and zc > 4.5
print("PASS" if ok else "FAIL")
