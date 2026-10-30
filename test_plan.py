"""Plan + reference executor == brute force on random graphs, weights and orders.
Kills: dropping the rho sign (edges inside the processed set ignored), wrong kernel span,
wrong pivot scatter in R. Negative control: perturbing one edge after planning must disagree."""
import random, numpy as np
from plan import plan
from ref import brute, run
rng = random.Random(3); worst = 0; ctrl = 1e9
for t in range(120):
    n = rng.randint(1, 10); p = rng.random()
    A = np.zeros((n, n), np.uint8)
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < p: A[i, j] = A[j, i] = 1
    w = [(complex(rng.gauss(0, 1), rng.gauss(0, 1)), complex(rng.gauss(0, 1), rng.gauss(0, 1))) for _ in range(n)]
    order = list(range(n)); rng.shuffle(order)
    zb = brute(A, w); zd = run(plan(A, order), w)
    worst = max(worst, abs(zb - zd) / max(1, abs(zb)))
    if n >= 3 and A[0, 1]:
        B = A.copy(); B[0, 1] = B[1, 0] = 0; ctrl = min(ctrl, abs(brute(B, w) - zd) / max(1, abs(zb)))
print(f"120 random cases: max rel err {worst:.2e}; negative control min rel diff {ctrl:.3f}")
assert worst < 1e-10 and ctrl > 1e-3
