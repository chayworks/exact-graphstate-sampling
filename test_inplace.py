"""In-place reference: the distribution implied by its chain (all branches) must equal the statevector,
with the arena fixed at 2^R entries (R = max cut rank). The executor asserts that no coset writes
another coset's input slots. Kills: a wrong g (output in the wrong slot), a wrong T_new, wrong term
grouping (tau), dropping the rho sign. Negative control: swapping the two outputs (tau) must fail."""
import random, itertools, numpy as np
from tableau import Tableau
from amp import Amplitudes, tail_matrix
from inplace import plan_inplace, sample_chain_probs
from test_tableau import G1, apply1, cz
TH = 0.9553166181245092
rng = random.Random(31); worst = 0; ctrl = 0; ncase = 0
for trial in range(30):
    n = rng.randint(2, 6); ops = []
    for _ in range(rng.randint(8, 60)):
        g = rng.choice(["h", "s", "sx", "x", "z", "cz", "cz", "cz"])
        ops.append((g, rng.sample(range(n), 2)) if g == "cz" else (g, [rng.randrange(n)]))
    tail = []
    for q in range(n): tail += [("tdg", [q]), ("h", [q]), (f"rz({TH})", [q]), ("h", [q])]
    T = Tableau(n); psi = np.zeros(2 ** n, complex); psi[0] = 1
    for name, qs in ops:
        T.apply(name, qs); psi = cz(psi, *qs, n) if name == "cz" else apply1(psi, G1[name], qs[0], n)
    for name, qs in tail: psi = apply1(psi, tail_matrix([(name, qs)]), qs[0], n)
    order = list(range(n)); rng.shuffle(order); am = Amplitudes(n, T, tail, order)
    R = max(max(s["r"], s["rn"]) for s in am.steps)
    steps = plan_inplace(am.A, order, max(R, 1))
    for path in itertools.product((0, 1), repeat=n):
        p = float(np.prod(sample_chain_probs(steps, max(R, 1), am.M, path)))
        b = [0] * n
        for st, bv in zip(steps, path): b[st["v"]] = bv
        idx = int("".join(map(str, b)), 2); worst = max(worst, abs(p - abs(psi[idx]) ** 2))
    ncase += 1
    if trial < 10 and any(s["m"] == 2 and s["r"] > 0 for s in steps):   # negative control: swap tau labels
        bad = [dict(s, terms=[(qv, x, 1 - t if s["m"] == 2 else t) for qv, x, t in s["terms"]]) for s in steps]
        dif = 0
        for path in itertools.product((0, 1), repeat=n):
            p = float(np.prod(sample_chain_probs(bad, max(R, 1), am.M, path)))
            b = [0] * n
            for st, bv in zip(steps, path): b[st["v"]] = bv
            dif = max(dif, abs(p - abs(psi[int("".join(map(str, b)), 2)]) ** 2))
        ctrl = max(ctrl, dif)
print(f"{ncase} circuits, every string, arena 2^R: max |p - p_sv| = {worst:.2e}; negative control (tau swapped) = {ctrl:.3f}")
assert worst < 1e-12 and ctrl > 1e-3
