"""Minimum stabilizer weight of the instance's graph state, and what it implies for marginals.

In the graph frame every stabilizer is S_v = X^v Z^{Av} (v in GF(2)^n, v != 0) with support v | Av, so its weight
is popcount(v | Av). By the stabilizer-group formula (stabprefix.py), the marginal of the output distribution on a
qubit set F is exactly uniform unless a nontrivial stabilizer is supported inside F; the measured Bloch vectors
here have no zero component, so the converse holds too. Hence every marginal on fewer than d_min qubits is
exactly uniform, where d_min = min_v popcount(v | Av).

Since weight(S_v) >= popcount(v), an exhaustive enumeration of all v with popcount(v) <= w - 1 that finds no
stabilizer lighter than w proves d_min = w once a weight-w stabilizer is known. This script does the search in
numpy for popcount(v) <= 4, then a randomised local search for an upper bound; minweight.cpp does the exhaustive
part for larger popcounts.

Usage: python minweight.py <qasm> [out_prefix]      (writes <out_prefix>_A.txt with the adjacency rows as 0/1 lines)"""
import sys, itertools, json, numpy as np
from tableau import clifford_body
from graphform import graph_form
from amp import tail_matrix
from stabprefix import bloch

def packrows(A):
    n = A.shape[0]; W = (n + 63) // 64
    R = np.zeros((n, W), np.uint64)
    for i in range(n):
        for j in np.nonzero(A[i])[0]: R[i, j // 64] |= np.uint64(1) << np.uint64(j % 64)
    return R

def popcount64(x):
    x = x - ((x >> np.uint64(1)) & np.uint64(0x5555555555555555))
    x = (x & np.uint64(0x3333333333333333)) + ((x >> np.uint64(2)) & np.uint64(0x3333333333333333))
    x = (x + (x >> np.uint64(4))) & np.uint64(0x0f0f0f0f0f0f0f0f)
    return (x * np.uint64(0x0101010101010101)) >> np.uint64(56)

def weight_of(R, E, idx):
    """weight of S_v for v = set of row indices idx (array of shape (m, k)); E = packed unit vectors."""
    v = np.zeros((idx.shape[0], R.shape[1]), np.uint64); Av = np.zeros_like(v)
    for c in range(idx.shape[1]):
        Av ^= R[idx[:, c]]; v ^= E[idx[:, c]]
    return popcount64(v | Av).sum(1)

if __name__ == "__main__":
    qasm = sys.argv[1]; out = sys.argv[2] if len(sys.argv) > 2 else "work/minweight"
    n, T, tail = clifford_body(qasm); A, U = graph_form(T)
    per = [[] for _ in range(n)]
    for name, qs in tail: per[qs[0]].append((name, qs))
    blochs = [bloch(tail_matrix(per[q]) @ U[q].conj().T) for q in range(n)]
    minabs = min(np.abs(b).min() for b in blochs)
    print(f"n = {n}, edges = {A.sum() // 2}, degrees min/mean/max = {A.sum(1).min()}/{A.sum(1).mean():.1f}/{A.sum(1).max()}")
    print(f"measured Bloch vectors in the graph frame: smallest |component| = {minabs:.6f} (all nonzero: {minabs > 1e-9})")
    with open(out + "_A.txt", "w") as f:
        for i in range(n): f.write("".join(str(int(x)) for x in A[i]) + "\n")
    R = packrows(A); E = packrows(np.eye(n, dtype=np.uint8)); best = (10 ** 9, None)
    # self-check of the packed arithmetic against the dense definition on random v
    An = A.astype(np.uint8); rng0 = np.random.default_rng(0)
    for _ in range(50):
        v = (rng0.random(n) < 0.3).astype(np.uint8); idx = np.nonzero(v)[0][None, :]
        if idx.shape[1] == 0: continue
        assert int(weight_of(R, E, idx)[0]) == int(np.count_nonzero(v | (An @ v % 2)))
    for k in range(1, 5):
        combos = np.array(list(itertools.combinations(range(n), k)), np.int64)
        w = weight_of(R, E, combos); i = int(w.argmin())
        print(f"|v| = {k}: {len(combos)} stabilizers, min weight {int(w[i])} at v = {combos[i].tolist()}")
        if w[i] < best[0]: best = (int(w[i]), combos[i].tolist())
    # randomised local search over larger |v|: flip one bit at a time, accept non-worsening moves
    rng = np.random.default_rng(1); An = A.astype(np.uint8)
    def wt(v): return int(np.count_nonzero(v | (An @ v % 2)))
    for trial in range(400):
        v = (rng.random(n) < rng.choice([0.05, 0.1, 0.2, 0.3])).astype(np.uint8)
        if not v.any(): continue
        cur = wt(v)
        for it in range(3000):
            j = rng.integers(n); v[j] ^= 1
            if not v.any(): v[j] ^= 1; continue
            w2 = wt(v)
            if w2 <= cur: cur = w2
            else: v[j] ^= 1
        if cur < best[0]: best = (cur, np.nonzero(v)[0].tolist()); print(f"local search trial {trial}: weight {cur}, |v| = {len(best[1])}")
    print(f"best stabilizer found: weight {best[0]}, v = {best[1]} (|v| = {len(best[1])})")
    json.dump({"n": n, "edges": int(A.sum() // 2), "best_weight": best[0], "v": best[1]}, open(out + "_best.json", "w"))
