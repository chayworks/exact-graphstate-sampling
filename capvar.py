"""Exact second moment of the hybrid sampler's prefix marginal, from the stabilizer group alone.

For a prefix P of K qubits measured in product bases with Bloch vectors n_j (all components +-1/sqrt(3) for the
face-state basis), z(b) = 2^K Pr(b_P) has E_uniform[z] = 1 and
    E_uniform[z^2] = 2^K sum_b Pr(b_P)^2 = sum over pairs (S, S') of stabilizers supported in P with supp S = supp S'
                     of sign(S) sign(S') prod_{j in supp} n_j[S_j] n_j[S'_j]
(the sum over b kills every pair whose supports differ: a qubit in one support only contributes sum_b (-1)^b = 0).
With |n_j[.]| = 1/sqrt(3) this is  sum over supports U of 3^{-|U|} (sum_{S: supp S = U} v(S))^2,  v(S) = +-1.
Var(z) = E[z^2] - 1, and Chebyshev gives Pr_uniform(z > c) <= Var(z) / (c - 1)^2: a rigorous bound on the mass
the capped sampler cannot produce, hence on its total-variation distance from the Born law of the prefix.

Usage (from graphsamp/):
  python capvar.py export <qasm> <order.json> <K> <out.txt>   generators of the group supported on the first K
        qubits of the order (x mask, z mask, phase) and the per-qubit sign table, for capvar.cpp
  python capvar.py grouped <qasm> <order.json> <K>             E[z^2] by the grouped formula in Python (K - rank <= ~20)
The group has 2^(K - cutrank) elements; capvar.cpp enumerates them by Gray code with phase tracking.
test_capvar.py checks the grouped formula against brute force over all prefixes on small random circuits."""
import sys, json, numpy as np
from tableau import clifford_body
from amp import tail_matrix
from stabprefix import bloch, kernel_gf2, _mul, stab_group_on

def setup(qasm, order_json, K):
    n, T, tail = clifford_body(qasm); order = json.load(open(order_json)); P = order[:K]
    per = [[] for _ in range(n)]
    for name, qs in tail: per[qs[0]].append((name, qs))
    blochs = [bloch(tail_matrix(per[q])) for q in range(n)]
    return n, T, blochs, P

def generators_on(T, P):
    """Generators (x, z, phase p with S = i^p X^x Z^z) of the stabilizers supported on P, as a GF(2) basis."""
    n = T.n; Pc = [j for j in range(n) if j not in set(P)]
    gens = [(T.xs[i].copy(), T.zs[i].copy(), (2 * int(T.r[i]) + int(np.count_nonzero(T.xs[i] & T.zs[i]))) % 4) for i in range(n)]
    M = np.concatenate([T.xs[:, Pc].T, T.zs[:, Pc].T], 0).astype(np.uint8)
    basis = kernel_gf2(M) if Pc else [np.eye(n, dtype=np.uint8)[i] for i in range(n)]
    out = []
    for a in basis:
        B = (np.zeros(n, np.uint8), np.zeros(n, np.uint8), 0)
        for i in np.nonzero(a)[0]: B = _mul(B, gens[i])
        assert not (B[0][Pc].any() or B[1][Pc].any()); out.append(B)
    return out

def v_of(sign, x, z, blochs, P):
    c = sign
    for j in P:
        if x[j] or z[j]:
            comp = 0 if (x[j] and not z[j]) else (1 if (x[j] and z[j]) else 2)
            c *= int(np.sign(blochs[j][comp]))
    return c

def grouped_Ez2(T, blochs, P, flip=None):
    """E[z^2] by the grouped formula; flip = index of one element whose sign is negated (negative control)."""
    G = stab_group_on(T, P); groups = {}
    for k, (s, x, z) in enumerate(G):
        if flip is not None and k == flip: s = -s
        U = tuple(int(x[j] or z[j]) for j in P); groups[U] = groups.get(U, 0) + v_of(s, x, z, blochs, P)
    return sum(3.0 ** (-sum(U)) * tot ** 2 for U, tot in groups.items()), len(G), groups

def export_group(T, blochs, P, path):
    gens = generators_on(T, P); pos = {q: i for i, q in enumerate(P)}
    with open(path, "w") as f:
        f.write(f"{len(P)} {len(gens)}\n")
        for j in P: f.write(" ".join(str(int(np.sign(blochs[j][c]))) for c in range(3)) + "\n")   # c_j[X], c_j[Y], c_j[Z]
        for x, z, p in gens:
            xm = sum(1 << pos[j] for j in P if x[j]); zm = sum(1 << pos[j] for j in P if z[j]); f.write(f"{xm} {zm} {p}\n")
    return len(gens)

if __name__ == "__main__":
    mode, qasm, oj, K = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    n, T, blochs, P = setup(qasm, oj, K)
    for j in P: assert np.allclose(np.abs(blochs[j]), 1 / np.sqrt(3), atol=1e-9), "not a face-state basis"
    if mode == "export":
        m = export_group(T, blochs, P, sys.argv[5])
        print(f"K = {K}, group of 2^{m} stabilizers supported on the prefix, generators written to {sys.argv[5]}")
    elif mode == "grouped":
        Ez2, ng, groups = grouped_Ez2(T, blochs, P)
        print(f"K = {K}: group size {ng}, distinct supports {len(groups)}, E[z^2] = {Ez2:.12f}, Var(z) = {Ez2 - 1:.6e}")
