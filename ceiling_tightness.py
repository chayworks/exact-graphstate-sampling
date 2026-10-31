"""Is the floor(n/2) ceiling tight? Best linear cut-rank width found for random graphs G(n, 1/2) and G(n, 0.3),
with the same order search as frontier.py, compared with the #151-style brickwork instances and the ceiling.
Usage: python ceiling_tightness.py 40,50,60,70 <iters> [seeds=3]   -> results/ceiling_tightness.jsonl"""
import sys, json, random, numpy as np
from frontier import best_width, ibm_like
from graphform import graph_form
from order import profile
ns = [int(x) for x in sys.argv[1].split(",")]; iters = int(sys.argv[2]); seeds = int(sys.argv[3]) if len(sys.argv) > 3 else 3
with open("results/ceiling_tightness.jsonl", "a") as f:
    for n in ns:
        for fam, p in [("G(n,1/2)", 0.5), ("G(n,0.3)", 0.3), ("brickwork", None)]:
            for s in range(seeds):
                if p is None:
                    A, _ = graph_form(ibm_like(n, 1000 + n + 17 * s))
                else:
                    rng = random.Random(100 * n + s); A = np.zeros((n, n), np.uint8)
                    for i in range(n):
                        for j in range(i + 1, n):
                            if rng.random() < p: A[i, j] = A[j, i] = 1
                nat = max(profile(A, list(range(n)))); k = best_width(A, n + s, 4, iters)
                rec = dict(n=n, family=fam, seed=s, edges=int(A.sum()) // 2, natural=nat, width=int(k[0]), ceiling=n // 2, gap=n // 2 - int(k[0]))
                print(json.dumps(rec), flush=True); f.write(json.dumps(rec) + "\n")
