"""stabprefix: prefix marginals from the stabilizer group must equal the statevector marginals for every prefix
of a random order and every string, on random small circuits with the #151 tail. Kills: Pauli phase
bookkeeping (Y = iXZ convention, sign of products), kernel enumeration, Bloch-vector sign, bit order.
Negative control: flipping the sign bit of one generator must disagree."""
import random, numpy as np
from amp import tail_matrix
from stabprefix import bloch, conditionals, stab_group_on, _cut_rank
from test_gsamp import build

rng = random.Random(5); worst = 0; ctrl = 0; nsteps = 0
for i in range(30):
    n = rng.randint(2, 8); T, tail, psi = build(n, rng, rng.randint(5, 60))
    order = list(range(n)); rng.shuffle(order)
    per = [[] for _ in range(n)]
    for name, qs in tail: per[qs[0]].append((name, qs))
    blochs = [bloch(tail_matrix(per[q])) for q in range(n)]
    p = np.abs(psi) ** 2                                              # index bit (n-1-q) is qubit q
    idx = np.arange(2 ** n); bits = np.array([(idx >> (n - 1 - q)) & 1 for q in range(n)]).T   # 2^n x n, column q
    B = bits[:, order].astype(np.uint8)                                # all strings, in step order
    cond, K = conditionals(T, blochs, order, B, dmax=n); assert K == n; nsteps += K
    pref = np.cumprod(cond, 1)                                         # Pr(b_1..b_t) for each string
    for t in range(1, n + 1):
        # statevector marginal of the first t qubits in the order, evaluated at each string's prefix
        marg = np.zeros(2 ** n)
        key = sum(B[:, k].astype(int) << (t - 1 - k) for k in range(t))
        acc = np.bincount(key, weights=p, minlength=2 ** t); marg = acc[key]
        worst = max(worst, np.abs(pref[:, t - 1] - marg).max())
    # group size check and negative control on the first circuit few
    assert len(stab_group_on(T, order)) == 2 ** n
    if i < 5:
        T.r[0] ^= 1; condb, _ = conditionals(T, blochs, order, B, dmax=n); T.r[0] ^= 1
        ctrl = max(ctrl, np.abs(np.cumprod(condb, 1)[:, -1] - p).max())
print(f"stabilizer-only prefix marginals vs statevector: 30 circuits, {nsteps} prefixes, every string: max diff {worst:.2e}; "
      f"negative control (one generator sign flipped) = {ctrl:.3f}")
assert worst < 1e-12 and ctrl > 1e-3
