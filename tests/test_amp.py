"""Born probabilities |<b|psi>|^2 from graph form + DP vs dense statevector of the same circuit.
Circuits mimic nq70: random Clifford body, then tdg; h; rz(theta); h on every qubit.
Kills: using U instead of U^dagger, wrong tail order (G @ M vs M @ G), wrong rz convention.
Negative control: omitting the conjugate transpose on U must disagree."""
import random, numpy as np
from tableau import Tableau
from amp import Amplitudes, tail_matrix
from ref import run
from test_tableau import G1, apply1, cz

TH = 0.9553166181245092
rng = random.Random(11); worst = 0; ctrl = 0
for trial in range(60):
    n = rng.randint(2, 7); ops = []
    for _ in range(rng.randint(5, 60)):
        g = rng.choice(["h", "s", "sx", "x", "z", "cz", "cz", "cz"])
        ops.append((g, rng.sample(range(n), 2)) if g == "cz" else (g, [rng.randrange(n)]))
    tail = []
    for q in range(n): tail += [("tdg", [q]), ("h", [q]), (f"rz({TH})", [q]), ("h", [q])]
    T = Tableau(n); psi = np.zeros(2 ** n, complex); psi[0] = 1
    for name, qs in ops:
        T.apply(name, qs); psi = cz(psi, *qs, n) if name == "cz" else apply1(psi, G1[name], qs[0], n)
    for name, qs in tail: psi = apply1(psi, tail_matrix([(name, qs)]), qs[0], n)
    order = list(range(n)); rng.shuffle(order)
    am = Amplitudes(n, T, tail, order)
    for idx in range(2 ** n):
        b = [(idx >> (n - 1 - q)) & 1 for q in range(n)]
        p = abs(run(am.steps, am.weights(b))) ** 2 / 2 ** n
        worst = max(worst, abs(p - abs(psi[idx]) ** 2))
    if trial < 10:  # negative control: U instead of U^dagger
        from graphform import graph_form
        _, U = graph_form(T); M = [am.W[q] @ U[q] for q in range(n)]
        bad = max(abs(abs(run(am.steps, [(M[q][b[q], 0], M[q][b[q], 1]) for q in range(n)])) ** 2 / 2 ** n - abs(psi[idx]) ** 2)
                  for idx in range(2 ** n) for b in [[(idx >> (n - 1 - q)) & 1 for q in range(n)]])
        ctrl = max(ctrl, bad)
print(f"60 circuits, all 2^n outcomes: max |p_dp - p_sv| = {worst:.2e}; negative control (U not U^dag) max diff = {ctrl:.3f}")
assert worst < 1e-12 and ctrl > 1e-3
