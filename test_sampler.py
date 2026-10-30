"""The sampler's implied distribution (product of its exact conditionals, every branch enumerated)
must equal |<b|psi>|^2 from a dense statevector, for random nq70-like circuits and random DP orders.
Also: per-step marginal 2^{-t} sum|P_t|^2 must equal the statevector marginal (checked via the products).
Kills: using |sum_c P[c]|^2 instead of sum_c |P[c]|^2 (negative control below), wrong 1/sqrt(2) scaling
(probabilities would not sum to 1 before normalising -- asserted), wrong weight row M[b, :]."""
import random, numpy as np
from tableau import Tableau
from amp import Amplitudes, tail_matrix
from sampler_ref import chain_distribution, conditionals, extend
from test_tableau import G1, apply1, cz
TH = 0.9553166181245092
rng = random.Random(21); worst = 0; worst_unnorm = 0; ctrl = 0
for trial in range(40):
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
    order = list(range(n)); rng.shuffle(order); am = Amplitudes(n, T, tail, order)
    dist = chain_distribution(am)
    for idx in range(2 ** n):
        b = tuple((idx >> (n - 1 - q)) & 1 for q in range(n))
        worst = max(worst, abs(dist.get(b, 0.0) - abs(psi[idx]) ** 2))
    _, nrm = conditionals(am, []); worst_unnorm = max(worst_unnorm, abs(nrm.sum() - 1))   # first-step marginals sum to 1 unnormalised
    if trial < 15:   # negative control: coherent sum instead of incoherent
        st = am.steps[0]; v = st["v"]
        c = [extend(np.array([1.0 + 0j]), st, (am.M[v][x, 0], am.M[v][x, 1])) / np.sqrt(2) for x in (0, 1)]
        wrong = np.array([abs(cc.sum()) ** 2 for cc in c]); wrong = wrong / wrong.sum()
        st2 = am.steps[1] if n > 1 else None
        if n > 1:
            # second-step conditional is where coherent vs incoherent first differ
            P = c[0] / np.linalg.norm(c[0]); v2 = st2["v"]
            c2 = [extend(P, st2, (am.M[v2][x, 0], am.M[v2][x, 1])) for x in (0, 1)]
            wr = np.array([abs(cc.sum()) ** 2 for cc in c2]); wr = wr / max(wr.sum(), 1e-300)
            rt = np.array([np.vdot(cc, cc).real for cc in c2]); rt = rt / rt.sum()
            ctrl = max(ctrl, abs(wr - rt).max())
print(f"40 circuits, every string: max |p_chain - p_sv| = {worst:.2e}; unnormalised first-step total off by {worst_unnorm:.1e}; "
      f"negative control (coherent sum) max conditional diff = {ctrl:.3f}")
assert worst < 1e-12 and worst_unnorm < 1e-12 and ctrl > 1e-2
