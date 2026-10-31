"""Simulated-annealing search for a low-width vertex order (linear rank-width upper bound), multi-seed.

State: a permutation. Cost: (max cut rank, number of steps at the max, log2 sum_k 2^{r_k}) compared
lexicographically, flattened to a scalar for annealing as  max*1e6 + count*1e3 + log2sum. Moves: relocate
one vertex by up to 12 positions, or swap two vertices within 12 positions. Starts from the noisy greedy
order of order2.greedy. Writes the best order per seed to <out_prefix>_s<seed>.json and prints a JSON line
per seed; the overall best goes to <out_prefix>_best.json.

Usage: python order_anneal.py <qasm> <out_prefix> <nseeds> <iters> <nproc> [T0=0.6] [T1=0.02]"""
import sys, json, math, random, multiprocessing as mp
from tableau import clifford_body
from graphform import graph_form
from order import profile
from order2 import greedy, key

def scalar(k): return k[0] * 1e6 + k[1] * 1e3 + k[2]

def anneal(args):
    A, n, seed, iters, T0, T1 = args
    rng = random.Random(seed); nb = [sum(int(A[i][j]) << j for j in range(n)) for i in range(n)]
    cur = min((greedy(nb, n, rng, 0.9) for _ in range(6)), key=lambda o: key(profile(A, o)))
    ck = key(profile(A, cur)); cs = scalar(ck); best, bk = cur[:], ck
    for it in range(iters):
        T = T0 * (T1 / T0) ** (it / max(1, iters - 1))
        new = cur[:]
        if rng.random() < 0.5:
            a, b = sorted(rng.sample(range(n), 2)); b = min(b, a + 12); x = new.pop(a); new.insert(b, x)
        else:
            a = rng.randrange(n); b = min(n - 1, a + rng.randint(1, 12)); new[a], new[b] = new[b], new[a]
        nk = key(profile(A, new)); ns = scalar(nk)
        if ns <= cs or rng.random() < math.exp(-(ns - cs) / (T * 1e3)):
            cur, ck, cs = new, nk, ns
            if ns < scalar(bk): best, bk = new[:], nk
    return seed, best, bk

if __name__ == "__main__":
    F, out, nseeds, iters, nproc = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    T0 = float(sys.argv[6]) if len(sys.argv) > 6 else 0.6; T1 = float(sys.argv[7]) if len(sys.argv) > 7 else 0.02
    n, T, _ = clifford_body(F); A, _ = graph_form(T)
    overall = None
    with mp.Pool(nproc) as p:
        for seed, best, bk in p.imap_unordered(anneal, [(A, n, s, iters, T0, T1) for s in range(nseeds)]):
            json.dump(best, open(f"{out}_s{seed}.json", "w"))
            print(json.dumps(dict(seed=seed, width=bk[0], at_max=bk[1], log2_sum=round(bk[2], 2))), flush=True)
            if overall is None or scalar(bk) < scalar(overall[1]): overall = (best, bk)
    json.dump(overall[0], open(f"{out}_best.json", "w"))
    print(json.dumps(dict(best_width=overall[1][0], at_max=overall[1][1], log2_sum=round(overall[1][2], 2), file=f"{out}_best.json")))
