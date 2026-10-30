"""Local-Clifford-equivalent graph form of a stabilizer state.

graph_form(T) returns (A, U) with A the symmetric 0/1 adjacency matrix and U a list of per-qubit 2x2
unitaries such that (U_0 x ... x U_{n-1}) |S> = e^{i phi} |G_A>, |G_A> = 2^{-n/2} sum_x (-1)^{sum_{i<j} A_ij x_i x_j} |x>.
"""
import numpy as np
from tableau import Tableau, gf2_rank

MAT = {"h": np.array([[1, 1], [1, -1]]) / np.sqrt(2), "sdg": np.diag([1, -1j]), "z": np.diag([1, -1])}

def _rowmul(T, i, k):
    """Generator i <- generator k * generator i, with the exact sign (Aaronson-Gottesman g function)."""
    x1, z1, x2, z2 = T.xs[k].astype(int), T.zs[k].astype(int), T.xs[i].astype(int), T.zs[i].astype(int)
    g = np.where((x1 == 1) & (z1 == 1), z2 - x2,
        np.where((x1 == 1) & (z1 == 0), z2 * (2 * x2 - 1),
        np.where((x1 == 0) & (z1 == 1), x2 * (1 - 2 * z2), 0)))
    ph = (2 * int(T.r[k]) + 2 * int(T.r[i]) + int(g.sum())) % 4
    assert ph in (0, 2)
    T.r[i] = ph // 2; T.xs[i] ^= T.xs[k]; T.zs[i] ^= T.zs[k]

def graph_form(T0):
    n = T0.n
    T = Tableau(n); T.xs, T.zs, T.r = T0.xs.copy(), T0.zs.copy(), T0.r.copy()
    gates = []
    def do(g, q): gates.append((g, q)); {"h": T.h, "sdg": T.sdg, "z": T.z}[g](q)
    # 1. Hadamards until the X block is invertible (greedy; asserted below)
    for _ in range(3):
        for j in range(n):
            if gf2_rank(T.xs) == n: break
            before = gf2_rank(T.xs); do("h", j)
            if gf2_rank(T.xs) <= before: do("h", j)   # undo
    assert gf2_rank(T.xs) == n, "could not make the X block invertible"
    # 2. Row-reduce X block to identity (sign-exact row products)
    for c in range(n):
        p = next(i for i in range(c, n) if T.xs[i, c])
        if p != c:
            for arr in (T.xs, T.zs): arr[[c, p]] = arr[[p, c]]
            T.r[[c, p]] = T.r[[p, c]]
        for i in range(n):
            if i != c and T.xs[i, c]: _rowmul(T, i, c)
    # 3. Clear Y's on the diagonal, then fix signs
    for i in range(n):
        if T.zs[i, i]: do("sdg", i)
    for i in range(n):
        if T.r[i]: do("z", i)
    assert (T.xs == np.eye(n, dtype=np.uint8)).all() and not T.r.any() and not np.diag(T.zs).any()
    A = T.zs.astype(np.uint8); assert (A == A.T).all()
    U = [np.eye(2, dtype=complex) for _ in range(n)]
    for g, q in gates: U[q] = MAT[g] @ U[q]
    return A, U
