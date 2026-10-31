"""Exact minimum stabilizer weight by a two-sided enumeration.

A stabilizer S_v = X^v Z^{Av} has weight w = |v ∪ Av| = |v| + |Av| − |v ∩ Av|. To decide whether any stabilizer has
weight ≤ W once all v with |v| ≤ t are exhausted (minweight.cpp), split the remaining cases:
  (1) |Av| ≤ t: enumerate c = Av with |c| ≤ t and v = A^-1c (minweight.cpp on A^-1; needs A invertible over GF(2));
  (2) |v| ≥ t+1 and |Av| ≥ t+1: then |v \\ Av| ≤ W − (t+1) and |Av \\ v| ≤ W − (t+1), so e = v ⊕ Av = (A+I)v has
      |e| ≤ 2(W − t − 1); enumerate e and solve (A+I)v = e (needs A+I invertible, else a small kernel coset).
With t = 11 and W = 14 case (2) is |e| ≤ 4, 970k systems. If neither case yields weight ≤ 14, d_min = 15 exactly
(a weight-15 stabilizer is known). Usage: python minweight_exact.py <A.txt> [t=11] [W=14]"""
import sys, itertools, numpy as np
A = np.array([[int(c) for c in l.strip()] for l in open(sys.argv[1]) if l.strip()], np.uint8); n = len(A)
t = int(sys.argv[2]) if len(sys.argv) > 2 else 11; W = int(sys.argv[3]) if len(sys.argv) > 3 else 14

def gf2_inverse(M):
    n = len(M); aug = np.concatenate([M.copy() % 2, np.eye(n, dtype=np.uint8)], 1); r = 0; piv = []
    for c in range(n):
        p = [i for i in range(r, n) if aug[i, c]]
        if not p: continue
        aug[[r, p[0]]] = aug[[p[0], r]]
        for i in range(n):
            if i != r and aug[i, c]: aug[i] ^= aug[r]
        piv.append(c); r += 1
    return (r, aug[:, n:] if r == n else None)

def weight(v): return int(np.count_nonzero(v | (A @ v % 2)))

rA, Ainv = gf2_inverse(A); B = (A + np.eye(n, dtype=np.uint8)) % 2; rB, Binv = gf2_inverse(B)
print(f"n = {n}; rank A = {rA}; rank (A+I) = {rB}")
if Ainv is not None:
    assert ((A @ Ainv) % 2 == np.eye(n, dtype=np.uint8)).all()
    with open(sys.argv[1].replace("_A.txt", "_Ainv.txt"), "w") as f:
        for i in range(n): f.write("".join(str(int(x)) for x in Ainv[i]) + "\n")
    print(f"A^-1 written (symmetric: {(Ainv == Ainv.T).all()}, nonzero diagonal entries: {int(Ainv.diagonal().sum())}); "
          f"case (1): run minweight.cpp on it with K = {t}, nodiag")
else:
    # generalized inverse: E A = R (RREF); for c in colspace(A), v0 = G c solves A v0 = c, solutions v0 + ker A;
    # consistency: the rows of E below the pivots (left kernel = ker A, A symmetric) must annihilate c.
    aug = np.concatenate([A.copy(), np.eye(n, dtype=np.uint8)], 1); r = 0; piv = []
    for c in range(n):
        p = [i for i in range(r, n) if aug[i, c]]
        if not p: continue
        aug[[r, p[0]]] = aug[[p[0], r]]
        for i in range(n):
            if i != r and aug[i, c]: aug[i] ^= aug[r]
        piv.append(c); r += 1
    E = aug[:, n:]; G = np.zeros((n, n), np.uint8)
    for i, pc in enumerate(piv): G[pc] = E[i]
    H = E[r:]                                   # left-kernel vectors (n - r of them)
    assert ((H @ A) % 2 == 0).all() and ((A @ H.T) % 2 == 0).all()
    rng = np.random.default_rng(0)
    for _ in range(200):                        # G c solves A v = c for every c in the column space
        v = (rng.random(n) < 0.5).astype(np.uint8); c = (A @ v) % 2; v0 = (G @ c) % 2
        assert ((A @ v0) % 2 == c).all() and ((H @ c) % 2 == 0).all()
    ks = [H[i] for i in range(len(H))]
    allk = [np.bitwise_xor.reduce([ks[i] for i in range(len(ks)) if (m >> i) & 1], 0) for m in range(1, 2 ** len(ks))]
    print(f"A is singular, kernel dimension {len(ks)}: the {len(allk)} nonzero kernel elements are stabilizers X^k of weight "
          f"{[int(k.sum()) for k in allk]} (these are case (1) with c = 0)")
    base = sys.argv[1].replace("_A.txt", "")
    with open(base + "_G.txt", "w") as f:
        for i in range(n): f.write("".join(str(int(x)) for x in G[i]) + "\n")
    with open(base + "_ker.txt", "w") as f:
        for k in ks: f.write("".join(str(int(x)) for x in k) + "\n")
    print(f"generalized inverse G and kernel written ({base}_G.txt, {base}_ker.txt); case (1): minweight_case1 with K = {t}")
if Binv is not None:
    assert ((B @ Binv) % 2 == np.eye(n, dtype=np.uint8)).all()
    emax = 2 * (W - t - 1); best = (10 ** 9, None); count = 0
    cols = Binv.T.copy()                    # v = Binv e = XOR of columns of Binv indexed by e
    for k in range(1, emax + 1):
        for idx in itertools.combinations(range(n), k):
            v = np.bitwise_xor.reduce(cols[list(idx)], 0); count += 1
            if v.any():
                w = weight(v)
                if w < best[0]: best = (w, np.nonzero(v)[0].tolist(), list(idx))
    print(f"case (2): {count} vectors e with 1 <= |e| <= {emax} solved; lightest stabilizer among them: weight {best[0]} "
          f"at v = {best[1]} (e = {best[2]})")
    if best[0] > W: print(f"case (2) excludes every stabilizer of weight <= {W} with |v| >= {t + 1} and |Av| >= {t + 1}")
    # self-check: the known weight-15 stabilizer must be reproduced by its own e
    v15 = np.zeros(n, np.uint8); v15[[0, 5, 27, 42, 48, 49, 52, 58, 61]] = 1; e15 = (B @ v15) % 2
    vb = np.bitwise_xor.reduce(cols[np.nonzero(e15)[0]], 0); print(f"self-check: known weight-15 v recovered from its e (|e| = {int(e15.sum())}): {(vb == v15).all()}")
else:
    print("A+I is singular: case (2) needs the kernel coset; not handled here")
