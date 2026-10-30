"""Simulated annealing over qubit orders, objective (max cut rank, log2 sum 2^{r_k}), several seeds in parallel.
Writes the best order per seed to work/order_seed<k>.json."""
import sys, json, math, random, multiprocessing as mp, numpy as np
from order import profile

def obj(pr): return max(pr) + math.log2(sum(2.0 ** r for r in pr)) / 64.0

def anneal(args):
    seed, iters, start = args
    A = np.load("A_nq70.npy"); n = len(A); rng = random.Random(seed)
    cur = start[:] if seed % 2 == 0 else rng.sample(range(n), n)
    cc = obj(profile(A, cur)); best = (cc, cur[:]); T = 1.0
    for it in range(iters):
        a, b = sorted(rng.sample(range(n), 2)); new = cur[:]
        u = rng.random()
        if u < 0.4: new[a], new[b] = new[b], new[a]
        elif u < 0.7: new[a:b + 1] = new[a:b + 1][::-1]
        else: x = new.pop(a); new.insert(b, x)
        nc = obj(profile(A, new))
        if nc <= cc or rng.random() < math.exp((cc - nc) / T):
            cur, cc = new, nc
            if cc < best[0]:
                best = (cc, cur[:]); json.dump(best[1], open(f"work/order_seed{seed}.json", "w"))
        T = max(0.01, T * (1 - 4.0 / iters))
    pr = profile(A, best[1])
    return seed, max(pr), math.log2(sum(2.0 ** r for r in pr))

if __name__ == "__main__":
    iters = int(sys.argv[1]); start = json.load(open("order_best.json"))
    with mp.Pool(4) as p:
        for seed, w, c in p.imap_unordered(anneal, [(s, iters, start) for s in range(4)]):
            print(f"seed {seed}: max width {w}, log2 cost {c:.2f}", flush=True)
