"""Is the output distribution Porter-Thomas, as far as our own samples can tell?
For samples drawn from p itself, z = 2^n p(x) over the samples follows the size-biased density z e^{-z} if the
distribution over strings is Exp(1) (Porter-Thomas): E[z] = 2, E[z^2] = 6, E[log2 p] = -n + (1 - gamma)/ln 2.
Prints the moments and the Kolmogorov-Smirnov distance to z e^{-z} (and, for contrast, to Exp(1), which is the wrong
reference for samples drawn by p). Usage, from graphsamp/:
  python porter_thomas_test.py results/samples_nq70.txt 70   -> results/porter_thomas_nq70.txt
  python porter_thomas_test.py results/nq80/samples_nq80.txt 80 --acc"""
import sys, math, numpy as np
fn, n = sys.argv[1], int(sys.argv[2]); acc = "--acc" in sys.argv
rows = [l.split() for l in open(fn) if l.strip() and not l.startswith("#")]
if acc: rows = [r[1:] for r in rows if r[0] == "ACC"]
lp = np.array([float(r[1]) for r in rows if len(r) > 1 and len(r[0]) == n]); N = len(lp); z = np.exp2(lp + n)
zs = np.sort(z); i = np.arange(1, N + 1)
def ks(cdf):
    F = cdf(zs); return max(np.max(i / N - F), np.max(F - (i - 1) / N))
ks_pt = ks(lambda x: 1 - (1 + x) * np.exp(-x)); ks_exp = ks(lambda x: 1 - np.exp(-x)); crit = 1.36 / math.sqrt(N)
print(f"{N} samples of the {n}-qubit instance from {fn}")
print(f"mean z = {z.mean():.4f} +/- {z.std(ddof=1) / math.sqrt(N):.4f} (Porter-Thomas: 2); mean z^2 = {np.mean(z * z):.3f} +/- {np.std(z * z, ddof=1) / math.sqrt(N):.3f} (Porter-Thomas: 6)")
print(f"mean log2 Pr = {lp.mean():.3f} +/- {lp.std(ddof=1) / math.sqrt(N):.3f} (Porter-Thomas: {-n + (1 - 0.5772156649) / math.log(2):.3f}; uniform: {-n})")
print(f"KS distance to z e^-z: {ks_pt:.4f} (5% level {crit:.4f}); to Exp(1), the wrong reference for samples drawn by p: {ks_exp:.4f}")
