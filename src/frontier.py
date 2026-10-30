"""Hardness frontier for IBM-style random graph state sampling (#151 construction): LNN chain, alternating
odd/even CZ brickwork, random sqrt(X) or S.sqrt(X) on every qubit after each layer, until CZ depth n.
For each n: graph form, natural-order width, and the best width from randomized greedy + short local
search (an upper bound on the true linear rank-width). Prints the exact-sampler cost model:
memory = 2 complex64 tables of 2^width; work = 4 * sum_k 2^{r_k} table updates per sample."""
import sys, math, random, json, numpy as np
from tableau import Tableau
from graphform import graph_form
from order import profile
from order2 import greedy, key

def ibm_like(n, seed):
    rng = random.Random(seed); T = Tableau(n)
    for q in range(n): T.h(q)
    for layer in range(n):
        for a in range(layer % 2, n - 1, 2): T.cz(a, a + 1)
        for q in range(n):
            if rng.random() < 0.5: T.s_(q)
            T.sx(q)
    return T

def best_width(A, seed, ngreedy, iters):
    n = len(A); nb = [sum(int(A[i][j]) << j for j in range(n)) for i in range(n)]; rng = random.Random(seed)
    best = min((key(profile(A, o)), o) for o in (greedy(nb, n, rng, 0.9) for _ in range(ngreedy)))
    cur, ck = best[1], best[0]
    for it in range(iters):
        a, b = sorted(rng.sample(range(n), 2)); b = min(b, a + 12); new = cur[:]
        x = new.pop(a); new.insert(b, x); nk = key(profile(A, new))
        if nk <= ck: cur, ck = new, nk
    return min(best[0], ck)

if __name__ == "__main__":
    out = []
    for n in [int(x) for x in sys.argv[1].split(",")]:
        A, _ = graph_form(ibm_like(n, 1000 + n))
        nat = max(profile(A, list(range(n))))
        k = best_width(A, n, 4, int(sys.argv[2]))
        rec = dict(n=n, edges=int(A.sum()) // 2, natural=nat, width=k[0], log2_work=round(math.log2(4) + k[2], 2),
                   mem_GB=round(2 * 8 * 2 ** k[0] / 1e9, 2))
        print(json.dumps(rec), flush=True); out.append(rec)
