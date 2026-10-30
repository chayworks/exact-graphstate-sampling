"""(U x ...)|S> must equal |G_A> up to a global phase, against dense statevectors of random circuits.
Kills: wrong sign handling in _rowmul (sign errors flip graph signs -> overlap < 1), wrong gate matrices,
wrong order of composition in U. Negative control: flipping one edge of A must drop the overlap."""
import numpy as np, random
from tableau import Tableau
from graphform import graph_form
from test_tableau import G1, apply1, cz

def graph_vec(A):
    n = len(A); v = np.zeros(2 ** n, complex)
    for idx in range(2 ** n):
        x = [(idx >> (n - 1 - q)) & 1 for q in range(n)]
        v[idx] = (-1) ** sum(int(A[i][j]) * x[i] * x[j] for i in range(n) for j in range(i + 1, n))
    return v / np.sqrt(2 ** n)

rng = random.Random(7); worst = 1.0; ctrl = 1.0
for trial in range(150):
    n = rng.randint(2, 7); ops = []
    for _ in range(rng.randint(1, 50)):
        g = rng.choice(["h", "s", "sdg", "sx", "x", "z", "cz", "cz"])
        ops.append((g, rng.sample(range(n), 2)) if g == "cz" else (g, [rng.randrange(n)]))
    T = Tableau(n); psi = np.zeros(2 ** n, complex); psi[0] = 1
    for name, qs in ops:
        T.apply(name, qs); psi = cz(psi, *qs, n) if name == "cz" else apply1(psi, G1[name], qs[0], n)
    A, U = graph_form(T)
    for q in range(n): psi = apply1(psi, U[q], q, n)
    worst = min(worst, abs(np.vdot(graph_vec(A), psi)))
    if n >= 3:
        B = A.copy(); B[0, 1] ^= 1; B[1, 0] ^= 1; ctrl = min(ctrl, abs(np.vdot(graph_vec(B), psi)))
print(f"min |<G_A|U|S>| over 150 circuits = {worst:.15f}; negative control (one edge flipped) max overlap <= {ctrl:.3f}")
assert abs(worst - 1) < 1e-12 and ctrl < 0.9
