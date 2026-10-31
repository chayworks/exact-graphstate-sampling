"""Exact minimum stabilizer weight by a three-way split (the Python part; the enumerations run in C++).

A stabilizer S_v = X^v Z^{Av} has weight w = |v ∪ Av| = |v| + |Av| − |v ∩ Av|. To decide whether any stabilizer has
weight ≤ W, split by t:
  (0) |v| ≤ t: enumerate v (minweight.cpp with K = t);
  (1) |Av| ≤ t: enumerate c = Av with |c| ≤ t, v ∈ G c + ker A (minweight_case1.cpp, G a generalised inverse; the
      case c = 0 is the kernel itself, listed here);
  (2) |v| ≥ t+1 and |Av| ≥ t+1: then |v \\ Av| ≤ W − (t+1) and |Av \\ v| ≤ W − (t+1), so e = v ⊕ Av = (A+I)v has
      |e| ≤ 2(W − t − 1); enumerate e, v ∈ G_B e + ker(A+I) (done here).
The minimum over (0), (1), (2) is exact for weights ≤ W; if it exceeds W and a weight-(W+1) stabilizer is known,
d_min = W + 1. Usage: python -P minweight_exact.py <base_A.txt> [t=11] [W=14]; writes <base>_G.txt, <base>_ker.txt."""
import sys, itertools, numpy as np
A = np.array([[int(c) for c in l.strip()] for l in open(sys.argv[1]) if l.strip()], np.uint8); n = len(A)
t = int(sys.argv[2]) if len(sys.argv) > 2 else 11; W = int(sys.argv[3]) if len(sys.argv) > 3 else 14
base = sys.argv[1][:-6] if sys.argv[1].endswith("_A.txt") else sys.argv[1]

def geninv(M):
    """(rank, G, H): E M = RREF; G c solves M v = c whenever c is in the column space, i.e. H c = 0 (H = left kernel)."""
    m = len(M); aug = np.concatenate([M.copy() % 2, np.eye(m, dtype=np.uint8)], 1); r = 0; piv = []
    for c in range(m):
        p = [i for i in range(r, m) if aug[i, c]]
        if not p: continue
        aug[[r, p[0]]] = aug[[p[0], r]]
        for i in range(m):
            if i != r and aug[i, c]: aug[i] ^= aug[r]
        piv.append(c); r += 1
    E = aug[:, m:]; G = np.zeros((m, m), np.uint8)
    for i, pc in enumerate(piv): G[pc] = E[i]
    H = E[r:]
    rng = np.random.default_rng(0)
    for _ in range(100):                      # G c solves M v = c on the column space; H annihilates it
        v = (rng.random(m) < 0.5).astype(np.uint8); c = (M @ v) % 2
        assert ((M @ ((G @ c) % 2)) % 2 == c).all() and ((H @ c) % 2 == 0).all()
    return r, G, H

def weight(v): return int(np.count_nonzero(v | (A @ v % 2)))
def span(H): return [np.bitwise_xor.reduce([H[i] for i in range(len(H)) if (m >> i) & 1], 0) for m in range(1, 2 ** len(H))]

rA, G, H = geninv(A); B = (A + np.eye(n, dtype=np.uint8)) % 2; rB, GB, HB = geninv(B)
print(f"n = {n}; rank A = {rA}; rank (A+I) = {rB}")
kw = [int(k.sum()) for k in span(H)]
print(f"kernel of A: dimension {len(H)}; its {len(kw)} nonzero elements are stabilizers X^k of weight {kw} (case (1) with c = 0)")
with open(base + "_G.txt", "w") as f:
    for i in range(n): f.write("".join(str(int(x)) for x in G[i]) + "\n")
with open(base + "_ker.txt", "w") as f:
    for k in H: f.write("".join(str(int(x)) for x in k) + "\n")
print(f"generalised inverse and kernel written ({base}_G.txt, {base}_ker.txt); case (1): minweight_case1 with K = {t}")
# case (2)
emax = 2 * (W - t - 1); best = (10 ** 9, None, None); count = 0; skipped = 0
cols = GB.T.copy(); cosetB = [np.zeros(n, np.uint8)] + span(HB)
for k in range(1, emax + 1):
    for idx in itertools.combinations(range(n), k):
        e = np.zeros(n, np.uint8); e[list(idx)] = 1
        if len(HB) and ((HB @ e) % 2).any(): skipped += 1; continue   # e not in the column space of A+I
        v0 = np.bitwise_xor.reduce(cols[list(idx)], 0); count += 1
        for kb in cosetB:
            v = v0 ^ kb
            if v.any():
                w = weight(v)
                if w < best[0]: best = (w, np.nonzero(v)[0].tolist(), list(idx))
print(f"case (2): e with 1 <= |e| <= {emax}: {count} in the column space of A+I solved ({skipped} outside it), kernel of A+I "
      f"of dimension {len(HB)}; lightest stabilizer among them: weight {best[0]} at v = {best[1]} (e = {best[2]})")
print(f"case (2) {'excludes' if best[0] > W else 'FINDS'} a stabilizer of weight <= {W} with |v| >= {t + 1} and |Av| >= {t + 1}")
vr = (np.random.default_rng(1).random(n) < 0.3).astype(np.uint8); er = (B @ vr) % 2
vb = np.bitwise_xor.reduce(cols[np.nonzero(er)[0]], 0)
print(f"self-check: a random v recovered from its e up to ker(A+I): {any(((vb ^ kb) == vr).all() for kb in cosetB)}")
