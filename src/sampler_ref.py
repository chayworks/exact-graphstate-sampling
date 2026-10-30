"""Exact qubit-by-qubit sampler (reference, pure Python).

Measuring the first t qubits of the DP order leaves (up to 2^{-t/2}) sum_c P_t[c] Z^{B_t c}|G_F>, and the
states Z^v|G_F> are orthonormal, so  Pr(b_1..b_t) = 2^{-t} sum_c |P_t[c]|^2.  No right environment is needed.
"""
import numpy as np

def extend(old, st, wv):
    """One DP step with weights wv = (w(0), w(1)) for the processed vertex."""
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
            acc += old[c] * wv[x] * (-1 if (x and bin(c & rho).count("1") & 1) else 1)
        new[cp] = acc
    return new

def conditionals(am, prefix_bits):
    """Pr(b_v = 0 / 1 | earlier bits) at the next step, plus the two candidate prefix vectors."""
    P = np.array([1.0 + 0j]); t = 0
    for st, bv in zip(am.steps, prefix_bits):
        v = st["v"]; P = extend(P, st, (am.M[v][bv, 0], am.M[v][bv, 1])) / np.sqrt(2)
    st = am.steps[len(prefix_bits)]; v = st["v"]
    cand = [extend(P, st, (am.M[v][b, 0], am.M[v][b, 1])) / np.sqrt(2) for b in (0, 1)]
    nrm = np.array([np.vdot(c, c).real for c in cand])
    return nrm / nrm.sum(), nrm

def chain_distribution(am):
    """Exact distribution implied by the sampler: product of its conditionals, for every string (small n)."""
    n = am.n; out = {}
    def rec(prefix, p):
        if len(prefix) == n:
            b = [0] * n
            for st, bv in zip(am.steps, prefix): b[st["v"]] = bv
            out[tuple(b)] = p; return
        cond, _ = conditionals(am, prefix)
        for bv in (0, 1):
            if cond[bv] > 0: rec(prefix + [bv], p * cond[bv])
    rec([], 1.0); return out

def sample(am, rng):
    P = np.array([1.0 + 0j]); b = [0] * am.n
    for st in am.steps:
        v = st["v"]
        cand = [extend(P, st, (am.M[v][x, 0], am.M[v][x, 1])) / np.sqrt(2) for x in (0, 1)]
        nrm = np.array([np.vdot(c, c).real for c in cand]); p1 = nrm[1] / nrm.sum()
        bv = int(rng.random() < p1); b[v] = bv; P = cand[bv] / np.sqrt(nrm[bv] / nrm.sum())
    return b
