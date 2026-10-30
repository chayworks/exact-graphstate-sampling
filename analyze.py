"""Summarise gsamp output for nq70: count, per-sample time, and the distribution of 2^n Pr(sample).

For exact samples, E[2^n p(x)] = 2^n sum_x p(x)^2 (the collision probability times 2^n), so the linear
cross-entropy F = 2^n <p> - 1 measured on our own samples estimates 2^n sum p^2 - 1. This is a
self-consistency statistic, not an independent check: the independent check is recomputing p(x) for
some samples with dp.cpp (xcheck70.py).
"""
import sys, json, math, numpy as np
lines = [l.split() for l in open(sys.argv[1]) if l.strip()]
n = 70; steps = json.load(open(sys.argv[2]))
assert all(len(b) == n for b, _, _ in lines), "truncated line"
lp = np.array([float(x[1]) for x in lines])
bits = [b for b, _, _ in lines]
print(f"samples: {len(lines)}  distinct: {len(set(bits))}  seeds: {lines[0][2]}..{lines[-1][2]}")
print(f"log2 Pr: mean {lp.mean():.3f}  min {lp.min():.3f}  max {lp.max():.3f}  (uniform would be -70)")
xeb = np.mean(np.exp2(lp + n)) - 1
se = np.std(np.exp2(lp + n)) / math.sqrt(len(lp))
print(f"linear XEB on own samples: {xeb:.4f} +/- {se:.4f}")
# qubit-order marginals (bits are in DP step order; map to circuit qubits)
B = np.array([[int(c) for c in b] for b in bits])
ones = np.zeros(n)
for k, q in enumerate(steps): ones[q] = B[:, k].mean()
print("mean bit value per circuit qubit, min/max:", ones.min().round(3), ones.max().round(3))
