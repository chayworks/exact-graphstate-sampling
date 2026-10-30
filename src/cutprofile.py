import sys, random, numpy as np
from tableau import clifford_body, cut_rank
n, T, tail = clifford_body(sys.argv[1])
print("n", n, "tail ops", len(tail), "tail gate set", sorted({t[0] for t in tail}))
prof = [cut_rank(T, range(k)) for k in range(1, n)]
print("prefix cut ranks (natural order), max", max(prof), ":", prof)
rng = random.Random(0)
rs = [cut_rank(T, rng.sample(range(n), n // 2)) for _ in range(300)]
print("random balanced bipartitions: min", min(rs), "mean", np.mean(rs), "max", max(rs))
# local search: swap pairs to lower the balanced cut rank
A = set(rng.sample(range(n), n // 2)); best = cut_rank(T, A)
for it in range(4000):
    a = rng.choice(sorted(A)); b = rng.choice(sorted(set(range(n)) - A))
    A2 = (A - {a}) | {b}; c = cut_rank(T, A2)
    if c <= best: A, best = A2, c
print("local-search balanced cut rank:", best)
