"""Better vertex orders: randomized greedy construction + local search on (max width, #steps at max, cost).
Writes work/order2_seed<k>.json for any order improving on max width 31."""
import sys, json, math, random, multiprocessing as mp, numpy as np
from order import profile

def rank_rows(rows):
    rk = 0; rows = [r for r in rows if r]
    while rows:
        p = rows.pop(); rk += 1; hb = p.bit_length() - 1
        rows = [r ^ p if (r >> hb) & 1 else r for r in rows]; rows = [r for r in rows if r]
    return rk

def cut(nb, P, Fmask):
    return rank_rows([nb[i] & Fmask for i in P])

def greedy(nb, n, rng, noise):
    P = []; full = (1 << n) - 1; Pm = 0
    for k in range(n):
        best = None
        for v in rng.sample(range(n), n):
            if (Pm >> v) & 1: continue
            Fm = full & ~(Pm | (1 << v)); c = cut(nb, P + [v], Fm) + noise * rng.random()
            if best is None or c < best[0]: best = (c, v)
        P.append(best[1]); Pm |= 1 << best[1]
    return P

def key(pr):
    m = max(pr); return (m, sum(1 for r in pr if r == m), math.log2(sum(2.0 ** r for r in pr)))

def work(seed):
    A = np.load("A_nq70.npy"); n = len(A); nb = [sum(int(A[i][j]) << j for j in range(n)) for i in range(n)]
    rng = random.Random(seed); best = None
    for g in range(6):
        o = greedy(nb, n, rng, 0.9); k = key(profile(A, o))
        if best is None or k < best[0]: best = (k, o)
    cur, ck = best[1], best[0]
    for it in range(25000):
        a, b = sorted(rng.sample(range(n), 2)); new = cur[:]; u = rng.random()
        if b - a > 12 and u < 0.7: b = a + rng.randint(1, 12)
        if u < 0.35: new[a], new[b] = new[b], new[a]
        elif u < 0.65: x = new.pop(a); new.insert(b, x)
        else: new[a:b + 1] = new[a:b + 1][::-1]
        nk = key(profile(A, new))
        if nk <= ck or rng.random() < 0.002:
            cur, ck = new, nk
            if ck < best[0]:
                best = (ck, cur[:])
                if ck[0] <= 30: json.dump(cur, open(f"work/order2_seed{seed}.json", "w"))
    return seed, best[0]

if __name__ == "__main__":
    with mp.Pool(int(sys.argv[2])) as p:
        for seed, k in p.imap_unordered(work, range(int(sys.argv[1]))):
            print(f"seed {seed}: max {k[0]} (x{k[1]} steps), log2 cost {k[2]:.2f}", flush=True)
