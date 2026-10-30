"""Pure-Python executor of a plan, and brute force, for testing."""
import itertools, numpy as np

def brute(A, w):
    n = len(A); tot = 0
    for x in itertools.product((0, 1), repeat=n):
        ph = sum(int(A[i][j]) * x[i] * x[j] for i in range(n) for j in range(i + 1, n)) & 1
        p = (-1) ** ph
        for j in range(n): p *= w[j][x[j]]
        tot += p
    return tot

def run(steps, w):
    old = np.array([1.0 + 0j])
    for st in steps:
        r, rn, rho, piv, ker = st["r"], st["rn"], st["rho"], st["piv"], st["ker"]
        kspan = [0]
        for k in ker: kspan += [s ^ k for s in kspan]
        new = np.zeros(1 << rn, complex)
        for cp in range(1 << rn):
            u0 = 0
            for j in range(rn):
                if (cp >> j) & 1: u0 |= 1 << piv[j]
            acc = 0
            for k in kspan:
                u = u0 ^ k; c = u & ((1 << r) - 1); x = (u >> r) & 1
                sgn = -1 if (x and bin(c & rho).count("1") & 1) else 1
                acc += old[c] * w[st["v"]][x] * sgn
            new[cp] = acc
        old = new
    return old[0]
