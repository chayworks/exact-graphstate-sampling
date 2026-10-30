"""Search a qubit order minimising the DP cost  sum_k 2^{r_k}  (r_k = GF(2) rank of A[prefix, suffix])."""
import sys, random, math, numpy as np
from tableau import clifford_body, cut_rank
from graphform import graph_form

def rank_bits(rows):
    rows = [r for r in rows if r]; rk = 0
    while rows:
        p = max(rows); rk += 1; hb = p.bit_length() - 1
        rows = [r ^ p if (r >> hb) & 1 else r for r in rows if r != p]
        rows = [r for r in rows if r]
    return rk

def profile(A, order):
    n = len(order); pos = {q: i for i, q in enumerate(order)}; out = []
    for k in range(1, n):
        P = order[:k]; F = order[k:]
        rows = [sum(int(A[i, j]) << t for t, j in enumerate(F)) for i in P]
        out.append(rank_bits(rows))
    return out

def cost(prof): return sum(2.0 ** r for r in prof)

if __name__ == "__main__":
    n, T, tail = clifford_body(sys.argv[1]); A, U = graph_form(T)
    nat = list(range(n)); pn = profile(A, nat)
    # cross-check against the tableau (independent formula, before any local Clifford)
    tab = [cut_rank(T, range(k)) for k in range(1, n)]
    assert pn == tab, "graph-form cut ranks disagree with the tableau"
    print("edges", int(A.sum()) // 2, " natural order: max", max(pn), " log2 cost %.2f" % math.log2(cost(pn)))
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
    cur = nat[:]; cc = cost(pn); best = (cc, cur[:]); Tt = 2.0
    for it in range(int(sys.argv[3]) if len(sys.argv) > 3 else 3000):
        a, b = sorted(rng.sample(range(n), 2)); new = cur[:]
        if rng.random() < 0.5: new[a], new[b] = new[b], new[a]
        else: new[a:b + 1] = new[a:b + 1][::-1]
        nc = cost(profile(A, new))
        if nc <= cc or rng.random() < math.exp((math.log2(cc) - math.log2(nc)) / Tt):
            cur, cc = new, nc
            if cc < best[0]: best = (cc, cur[:])
        Tt = max(0.02, Tt * 0.998)
    pb = profile(A, best[1])
    print("best: max", max(pb), " log2 cost %.2f" % math.log2(best[0]))
    print("profile", pb); print("order", best[1])
    np.save("A_nq70.npy", A)
