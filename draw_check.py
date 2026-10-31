"""Does the sampler draw each bit with its own conditional? A direct test of the drawing step, given the conditionals.

chain mode prints, for every sample, Pr(b_k = observed | b_1..b_{k-1}) =: c_k at every step. If the sampler drew b_k
as (U < Pr(b_k = 1 | prefix)) with an unbiased uniform U, the residual b_k - p_k (p_k = Pr(b_k = 1 | prefix) = c_k if
b_k = 1 else 1 - c_k) has mean 0 and variance p_k(1 - p_k) exactly, whatever the prefix. So per step
z_k = sum_i (b_ik - p_ik) / sqrt(sum_i p_ik(1 - p_ik)) ~ N(0, 1), and sum_k z_k^2 over steps with non-degenerate
conditionals is chi-square. This removes the prefix-induced variation that the raw bit means carry, and tests
exactly the mechanism that the raw single-bit test (hw_marginals.py) left in doubt.
Usage: python -P draw_check.py results/chain/ideal_fwd_chain.txt cloud4/ideal_fwd.txt"""
import sys, numpy as np
from scipy.stats import chi2
n = 70
rows = [l.split() for l in open(sys.argv[1]) if len(l.split()) == n + 1]
C = np.array([[float(x) for x in r[1:]] for r in rows])          # column 0 is log2 Pr; strings come from the input file
S = [l.strip() for l in open(sys.argv[2]) if len(l.strip()) == n][:len(C)]; assert len(S) == len(C)
B = np.array([[int(ch) for ch in s] for s in S], np.int8)
lp = np.array([float(r[0]) for r in rows]); assert np.abs(np.log2(C).sum(1) - lp).max() < 1e-5, "conditionals do not multiply to the printed log2 Pr"
P = np.where(B == 1, C, 1 - C)                     # Pr(b_k = 1 | prefix) for each sample and step
var = (P * (1 - P)).sum(0); ok = var > 1e-6        # steps whose conditionals are not all exactly 0/1
z = (B - P).sum(0)[ok] / np.sqrt(var[ok]); k = int(ok.sum())
print(f"{len(B)} samples; steps with non-degenerate conditionals: {k}")
chi = float((z ** 2).sum()); print(f"drawing test: sum z^2 = {chi:.1f} on {k} dof, p = {chi2.sf(chi, k):.3f}; max |z| {np.abs(z).max():.2f}")
for a, b in [(0, 22), (22, 50), (50, 70)]:
    m = ok.copy(); m[:a] = False; m[b:] = False; zz = (B - P).sum(0)[m] / np.sqrt(var[m]); s = float((zz ** 2).sum())
    print(f"  steps {a + 1}-{b}: sum z^2 = {s:.1f} on {int(m.sum())} dof, p = {chi2.sf(s, int(m.sum())):.3f}")
# controls: (a) a sampler that draws with the wrong probability p' = p^1.1 (mild bias); (b) bits shifted by one sample
rng = np.random.default_rng(0); Pb = P ** 1.1; Bb = (rng.random(P.shape) < Pb).astype(np.int8)
zb = (Bb - P).sum(0)[ok] / np.sqrt(var[ok]); s = float((zb ** 2).sum()); print(f"control (bias p -> p^1.1, synthetic draws): sum z^2 = {s:.1f} on {k} dof, p = {chi2.sf(s, k):.2e}")
Bg = (rng.random(P.shape) < P).astype(np.int8); zg = (Bg - P).sum(0)[ok] / np.sqrt(var[ok]); s = float((zg ** 2).sum())
print(f"control (correct synthetic draws from the same conditionals): sum z^2 = {s:.1f} on {k} dof, p = {chi2.sf(s, k):.3f}")