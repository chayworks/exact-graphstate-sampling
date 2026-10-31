"""Simulation-free test of IBM's hardware shots: low-order marginals.

For the ideal output distribution, every marginal on fewer than d_min qubits is exactly uniform (d_min = minimum
stabilizer weight of the graph state, see minweight.py). So every single-bit mean, pair parity and triple parity
of the ideal distribution is exactly 0 (for +-1 bits), and any deviation seen in the hardware shots is noise,
detectable without any simulation. Reports z-scores of all 1-, 2- and 3-bit parities of the shots, the fraction
beyond 2 sigma (4.6% expected if at chance), the largest, and the same statistics for our exact samples (which
must be at chance) and for a control with a biased bit.

Usage: python hw_marginals.py <ibm_strings.txt> <ideal_samples.txt> <step_qubits.json> [n]
Strings are in the plan's step order; step_qubits.json maps step -> qubit for attribution."""
import sys, json, itertools, numpy as np
from scipy.stats import chi2
S, I, Q = sys.argv[1:4]; n = int(sys.argv[4]) if len(sys.argv) > 4 else 70
sq = json.load(open(Q))
def load(p):
    rows = [l.split()[0] for l in open(p) if l.strip()]; rows = [r for r in rows if len(r) == n and set(r) <= {"0", "1"}]
    return 1 - 2 * np.array([[int(c) for c in r] for r in rows], np.int8)   # +1 for 0, -1 for 1
ZS = {}
def report(name, X):
    N = len(X); out = {}; ZS[name] = {}
    for k in (1, 2, 3):
        combos = list(itertools.combinations(range(n), k)); z = np.empty(len(combos))
        for i, c in enumerate(combos):
            p = np.ones(N, np.int8)
            for j in c: p = p * X[:, j]
            z[i] = p.mean() * np.sqrt(N)          # SE of a +-1 mean is 1/sqrt(N) under uniformity
        ZS[name][k] = z.copy()
        frac = float((np.abs(z) > 2).mean()); big = int(np.abs(z).argmax()); chi = float((z ** 2).sum()); pch = float(chi2.sf(chi, len(combos)))
        out[k] = dict(tests=len(combos), frac_beyond_2sigma=frac, max_z=float(np.abs(z[big])), where=[sq[j] for j in combos[big]],
                      mean_abs_z=float(np.abs(z).mean()), expected_mean_abs_z=float(np.sqrt(2 / np.pi)), sum_z2=chi, chi2_p=pch)
        print(f"{name:>22} {k}-bit: {len(combos):>6} tests, {100 * frac:4.1f}% beyond 2 sigma (4.6% expected), "
              f"max |z| {np.abs(z[big]):.2f} at qubits {out[k]['where']}, mean |z| {np.abs(z).mean():.3f} (expected 0.798), "
              f"sum z^2 {chi:.0f} on {len(combos)} dof (p = {pch:.3f})")
        if k == 1:
            order = np.argsort(-np.abs(z))[:5]
            print(f"{'':>22}   largest single-bit biases: " + ", ".join(f"q{sq[j]}: mean(b)={0.5 - z[j] / (2 * np.sqrt(N)):.3f} z={z[j]:+.2f}" for j in order))
    # 4-bit parities, vectorised: products of all triples with max index < j, times bit j
    tri = list(itertools.combinations(range(n), 3)); P3 = np.ones((N, len(tri)), np.int8)
    for i, c in enumerate(tri): P3[:, i] = X[:, c[0]] * X[:, c[1]] * X[:, c[2]]
    tmax = np.array([c[2] for c in tri]); chi = 0.0; tests = 0; big = (0.0, None); z4 = []
    for j in range(3, n):
        m = tmax < j; z = (P3[:, m] * X[:, j:j + 1]).mean(0) * np.sqrt(N); chi += float((z ** 2).sum()); tests += int(m.sum()); z4.append(z)
        a = int(np.abs(z).argmax())
        if abs(z[a]) > big[0]: big = (float(abs(z[a])), [sq[q] for q in tri[np.nonzero(m)[0][a]]] + [sq[j]])
    ZS[name][4] = np.concatenate(z4)
    pch = float(chi2.sf(chi, tests)); out[4] = dict(tests=tests, sum_z2=chi, chi2_p=pch, max_z=big[0], where=big[1])
    print(f"{name:>22} 4-bit: {tests:>6} tests, max |z| {big[0]:.2f} at qubits {big[1]}, sum z^2 {chi:.0f} on {tests} dof (p = {pch:.3f})")
    return out
hw, ideal = load(S), load(I)
print(f"hardware shots N = {len(hw)}, ideal samples N = {len(ideal)}, n = {n}")
res = {"hardware": report("hardware shots", hw), "ideal": report("ideal (exact samples)", ideal)}
ctrl = ideal.copy(); rng = np.random.default_rng(0); flip = rng.random(len(ctrl)) < 0.3; ctrl[flip, 7] = 1   # bias one bit (step 8) by 15%
res["control"] = report("control: 1 bit biased 15%", ctrl)
json.dump(res, open("results/hw_marginals.json", "w"), indent=1)
np.savez_compressed("results/hw_marginals_z.npz", **{f"{'hw' if 'hardware' in nm else 'ideal' if 'ideal' in nm else 'ctrl'}_k{k}": z for nm, d in ZS.items() for k, z in d.items()})
