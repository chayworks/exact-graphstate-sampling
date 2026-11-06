"""Independent check of the exact sampler at full size, where no statevector exists.

IQPopt trains on the classical formula (Recio-Armengol et al., arXiv:2503.02934) for |psi> = exp(i sum theta_j X_gj)|0>:
    <Z_a> = E_{z ~ uniform} cos( sum_j 2 theta_j [g_j . a odd] (-1)^{g_j . z} ),
a Monte Carlo average that shares no code with the transfer-matrix sampler. Step 1 checks the formula (evaluated
exactly over all z) against the dense statevector on small circuits; step 2 compares, for a trained 153-qubit model,
every single-qubit <Z_q> and 300 random two-qubit <Z_a Z_b> between the formula (2^16 Monte Carlo z) and our exact
samples (n_s exact samples), in units of their combined standard error.
Usage: python check153.py <IQP-NISQ results dir> [n_s]"""
import sys, numpy as np
from pathlib import Path
from heavyhex import gate_list
from iqp_exact import IQPModel, statevector_probs

def zexp_formula(gates, thetas, a, zs):
    """E_z cos(sum_j 2 theta_j [|g_j & a| odd] (-1)^{g_j . z}) over the rows of zs (bits by qubit)."""
    arg = np.zeros(len(zs))
    for g, t in zip(gates, thetas):
        if sum(1 for q in g if q in a) % 2 == 1:
            par = np.zeros(len(zs), np.int64)
            for q in g: par ^= zs[:, q]
            arg += 2 * t * (1 - 2 * par)
    v = np.cos(arg); return v.mean(), v.std() / np.sqrt(len(v))

rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else 11)
# 1. formula vs dense statevector, exact over all z
worst = 0
for trial in range(20):
    n = int(rng.integers(3, 10)); gates = [(i,) for i in range(n)] + [(i, (i + 1) % n) for i in range(n - 1)] + [(0, n - 1)]
    th = rng.uniform(-np.pi, np.pi, len(gates)); P = statevector_probs(n, gates, th)
    allz = (np.arange(1 << n)[:, None] >> np.arange(n)) & 1
    for _ in range(4):
        a = set(int(q) for q in rng.choice(n, int(rng.integers(1, 3)), replace=False))
        exact = sum(P[y] * (-1) ** sum((y >> q) & 1 for q in a) for y in range(1 << n))
        worst = max(worst, abs(zexp_formula(gates, th, a, allz)[0] - exact))
print(f"step 1: formula vs statevector, 80 expectation values on up to 9 qubits: max |diff| = {worst:.1e}")
assert worst < 1e-12

# 2. full size
B = Path(sys.argv[1]); ns = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
pp = sorted(B.glob("trained_params_iqp/18N_*/*/*/rank_0/params.npy"))[0]
th = np.load(pp, allow_pickle=False); N = 153; g = gate_list(N); L = len(th) // len(g); g = g * L
m = IQPModel(N, g, th); S = np.array([y for y, _ in m.sample(rng, ns)]); spins = 1 - 2 * S.astype(float)
zs = rng.integers(0, 2, (1 << 16, N))
def pull(f, sf, e):
    # error of a mean of ns +-1 values whose expectation is f: sqrt((1 - f^2) / ns), the null hypothesis's own variance
    return (f - e) / np.sqrt(sf ** 2 + max(1 - f * f, 1e-12) / ns)
pulls = []
for q in range(N):
    f, sf = zexp_formula(g, th, {q}, zs); e = spins[:, q].mean(); pulls.append(pull(f, sf, e))
pairs = [tuple(int(x) for x in rng.choice(N, 2, replace=False)) for _ in range(300)]
for a, b in pairs:
    f, sf = zexp_formula(g, th, {a, b}, zs); e = (spins[:, a] * spins[:, b]).mean(); pulls.append(pull(f, sf, e))
pulls = np.array(pulls)
print(f"step 2: {pp.parent.relative_to(B).as_posix()}, {N} qubits, {ns} exact samples vs formula with 2^16 z:")
print(f"  {len(pulls)} expectation values; pulls mean {pulls.mean():+.3f}, sd {pulls.std():.3f} (1 if they agree), "
      f"max |pull| {np.abs(pulls).max():.2f}, fraction beyond 3: {np.mean(np.abs(pulls) > 3):.4f}")
# control: a sampler whose largest coupler angle is scaled by 1.5 must be caught on <Z_a> of that coupler's end
# (a coupler's own angle drops out of <Z_a Z_b>, whose overlap with it is even)
k = N + int(np.argmax(np.abs(th[N:len(th) // L]))); a, b = gate_list(N)[k]
th2 = th.copy(); th2[k] = 1.5 * th2[k]
sp2 = 1 - 2 * np.array([y for y, _ in IQPModel(N, g, th2).sample(rng, ns)]).astype(float)
for q in (a, b):
    f, sf = zexp_formula(g, th, {q}, zs)
    print(f"  control (coupler {a}-{b} angle x1.5 in the sampler): pull on <Z{q}> = {pull(f, sf, sp2[:, q].mean()):.1f}; "
          f"same qubit, correct sampler: {pull(f, sf, spins[:, q].mean()):.1f}")
