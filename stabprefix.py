"""Exact prefix probabilities from the stabilizer group alone, with no dynamic program.

The state before the single-qubit tails is a stabilizer state rho = 2^-n sum_S S over its 2^n signed
stabilizers. Measuring qubit j after the tail W_j in the Z basis is the POVM E_j(b) = (I + (-1)^b n_j.sigma)/2
with Bloch vector n_j = Tr(sigma W_j^dag Z W_j)/2. For a set F of t qubits,

    Pr(b_F) = 2^-t * sum_{S supported on F} sign(S) * prod_{j in supp S} (-1)^{b_j} n_j[S_j].

The stabilizers supported on F form a group of 2^(t - cutrank(F)) elements, found as a GF(2) kernel; while
the cut rank equals t the sum is the identity alone and the marginal is exactly uniform. With deficit
d = t - cutrank(F) <= ~16 the sum is enumerated exactly. This gives every conditional Pr(b_k | b_1..b_{k-1})
along a vertex order for the first steps independently of the cut-rank DP, and explains why those
conditionals sit so close to 1/2: each nontrivial term is a product of |n| <= 1 components over the support
of a stabilizer, so a stabilizer of weight w contributes at most max|n|^w.

Usage: python stabprefix.py <qasm> <order.json> <strings.txt> <chain_output.txt> [dmax]
Compares the stabilizer-only conditionals with gsamp/gsamp2 chain-mode output for every string."""
import sys, json, numpy as np

SIG = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]

def bloch(W):
    O = W.conj().T @ np.diag([1.0, -1.0]) @ W
    n = np.array([np.trace(s @ O).real / 2 for s in SIG]); assert abs(np.linalg.norm(n) - 1) < 1e-12; return n

def kernel_gf2(M):
    """Basis of {a : M a = 0 mod 2} for M (m x n) uint8; returns list of uint8 vectors."""
    M = M.copy() % 2; m, n = M.shape; piv = []; r = 0
    for c in range(n):
        rows = [i for i in range(r, m) if M[i, c]]
        if not rows: continue
        M[[r, rows[0]]] = M[[rows[0], r]]
        for i in range(m):
            if i != r and M[i, c]: M[i] ^= M[r]
        piv.append(c); r += 1
        if r == m: break
    free = [c for c in range(n) if c not in piv]; basis = []
    for f in free:
        a = np.zeros(n, np.uint8); a[f] = 1
        for i, c in enumerate(piv):
            if M[i, f]: a[c] = 1
        basis.append(a)
    return basis

def _mul(A, B):                       # Paulis as (x, z, p): i^p prod_j X^x_j Z^z_j
    x, z, p = A; x2, z2, p2 = B
    return x ^ x2, z ^ z2, (p + p2 + 2 * int(np.count_nonzero(z & x2))) % 4

def stab_group_on(T, F):
    """All signed stabilizers of the tableau state supported on the qubit set F: list of (sign, x, z)."""
    n = T.n; Fc = [j for j in range(n) if j not in set(F)]
    gens = [(T.xs[i].copy(), T.zs[i].copy(), (2 * int(T.r[i]) + int(np.count_nonzero(T.xs[i] & T.zs[i]))) % 4) for i in range(n)]
    if Fc:
        M = np.concatenate([T.xs[:, Fc].T, T.zs[:, Fc].T], 0).astype(np.uint8)   # constraints x rows: generators
        basis = kernel_gf2(M)
    else:
        basis = [np.eye(n, dtype=np.uint8)[i] for i in range(n)]
    elems = [(np.zeros(n, np.uint8), np.zeros(n, np.uint8), 0)]
    for a in basis:
        B = (np.zeros(n, np.uint8), np.zeros(n, np.uint8), 0)
        for i in np.nonzero(a)[0]: B = _mul(B, gens[i])
        elems += [_mul(e, B) for e in elems]
    out = []
    for x, z, p in elems:
        q = (p - int(np.count_nonzero(x & z))) % 4; assert q in (0, 2), "non-Hermitian product"
        assert not (x[Fc].any() or z[Fc].any())
        out.append((1 if q == 0 else -1, x, z))
    return out

def prefix_prob_factors(T, blochs, F, B):
    """2^t * Pr(b_F) for each row of B (bits of the qubits F, in F order). blochs[j] is n_j."""
    G = stab_group_on(T, F); t = len(F); A = np.zeros(len(G)); Supp = np.zeros((len(G), t), np.uint8)
    for k, (s, x, z) in enumerate(G):
        c = s
        for pos, j in enumerate(F):
            if x[j] or z[j]:
                comp = 0 if (x[j] and not z[j]) else (1 if (x[j] and z[j]) else 2)
                c *= blochs[j][comp]; Supp[k, pos] = 1
        A[k] = c
    par = (B.astype(np.int64) @ Supp.T.astype(np.int64)) % 2
    return (1 - 2 * par).astype(float) @ A, len(G)

def conditionals(T, blochs, order, B, dmax=16):
    """Pr(b_k | b_1..b_{k-1}) for k = 1.. while the deficit <= dmax; B rows are strings in step order.
    Returns (cond array N x K, K)."""
    N = B.shape[0]; prev = np.ones(N); cond = []
    for t in range(1, len(order) + 1):
        F = order[:t]; r = _cut_rank(T, F)
        if t - r > dmax: break
        f, ng = prefix_prob_factors(T, blochs, F, B[:, :t]); assert ng == 2 ** (t - r), (ng, t, r)
        cur = f / 2.0 ** t; cond.append(cur / prev); prev = cur
    return np.array(cond).T, len(cond)

def _cut_rank(T, F):
    from tableau import gf2_rank
    n = T.n; Fc = [j for j in range(n) if j not in set(F)]
    if not Fc: return 0
    M = np.concatenate([T.xs[:, Fc].T, T.zs[:, Fc].T], 0)
    return len(F) - (n - gf2_rank(M))          # kernel dim of the constraint = t - cutrank  =>  cutrank = t - (n - rank M)

if __name__ == "__main__":
    from tableau import clifford_body
    from amp import tail_matrix
    F, OJ, SF, CF = sys.argv[1:5]; dmax = int(sys.argv[5]) if len(sys.argv) > 5 else 16
    n, T, tail = clifford_body(F); order = json.load(open(OJ))
    per = [[] for _ in range(n)]
    for name, qs in tail: per[qs[0]].append((name, qs))
    blochs = [bloch(tail_matrix(per[q])) for q in range(n)]
    rows = [l.split() for l in open(CF) if l.strip()]; rows = [r for r in rows if len(r) == n + 1]
    chain = np.array([[float(x) for x in r[1:]] for r in rows])
    strings = [l.strip() for l in open(SF) if l.strip()][:len(chain)]
    B = np.array([[int(c) for c in s] for s in strings], np.uint8)
    cond, K = conditionals(T, blochs, order, B, dmax)
    print(f"{len(B)} strings; stabilizer-only conditionals for steps 1..{K} (deficit <= {dmax})")
    print(f"{'step':>4} {'deficit':>7} {'max|c-1/2| stab':>16} {'max|c_stab - c_chain|':>22}")
    for k in range(K):
        d = (k + 1) - _cut_rank(T, order[:k + 1])
        print(f"{k + 1:>4} {d:>7} {np.abs(cond[:, k] - 0.5).max():>16.3e} {np.abs(cond[:, k] - chain[:, k]).max():>22.3e}")
    print(f"overall max |c_stab - c_chain| over {K} steps x {len(B)} strings: {np.abs(cond - chain[:, :K]).max():.3e}")
