"""Check Tableau against dense statevectors: every generator must stabilize C|0> with the right sign.

Gate matrices are the qelib1/Qiskit definitions. Kills: wrong sign rule in h/s/cz (sign check),
wrong sx conjugation (Z -> -Y), swapped x/z columns. Negative control: flipping one sign must fail.
"""
import numpy as np, random
from tableau import Tableau

I2 = np.eye(2); X = np.array([[0, 1], [1, 0]]); Z = np.diag([1, -1]); Y = 1j * X @ Z
G1 = {"h": np.array([[1, 1], [1, -1]]) / np.sqrt(2), "s": np.diag([1, 1j]), "sdg": np.diag([1, -1j]),
      "sx": np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]]) / 2, "x": X, "z": Z}

def apply1(psi, U, q, n):
    psi = psi.reshape([2] * n); psi = np.moveaxis(np.tensordot(U, psi, axes=([1], [q])), 0, q); return psi.reshape(-1)

def cz(psi, a, b, n):
    psi = psi.reshape([2] * n).copy(); idx = [slice(None)] * n; idx[a] = 1; idx[b] = 1; psi[tuple(idx)] *= -1; return psi.reshape(-1)

def pauli_op(xr, zr, sign, n):
    psi_op = np.array([[1]])
    for q in range(n):  # qubit 0 is the most significant axis, matching reshape([2]*n)
        P = {(0, 0): I2, (1, 0): X, (0, 1): Z, (1, 1): Y}[(int(xr[q]), int(zr[q]))]
        psi_op = np.kron(psi_op, P)
    return (-1) ** int(sign) * psi_op

def check(n, ops, flip=False):
    T = Tableau(n); psi = np.zeros(2 ** n, complex); psi[0] = 1
    for name, qs in ops:
        T.apply(name, qs)
        psi = cz(psi, *qs, n) if name == "cz" else apply1(psi, G1[name], qs[0], n)
    if flip: T.r[0] ^= 1
    return max(np.linalg.norm(pauli_op(T.xs[i], T.zs[i], T.r[i], n) @ psi - psi) for i in range(n))

if __name__ == "__main__":
    rng = random.Random(1)
    worst = 0
    for trial in range(200):
        n = rng.randint(2, 6); ops = []
        for _ in range(rng.randint(1, 40)):
            g = rng.choice(["h", "s", "sdg", "sx", "x", "z", "cz", "cz"])
            ops.append((g, rng.sample(range(n), 2)) if g == "cz" else (g, [rng.randrange(n)]))
        worst = max(worst, check(n, ops))
    ctrl = check(4, [("h", [0]), ("cz", [0, 1]), ("sx", [2]), ("s", [3]), ("h", [3])], flip=True)
    print(f"200 random circuits: max |P psi - psi| = {worst:.2e}   negative control (one sign flipped): {ctrl:.2f}")
    assert worst < 1e-12 and ctrl > 1
