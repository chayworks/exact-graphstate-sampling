"""Fidelity of IBM's #151 hardware run from the exact probabilities of its shots, two ways.

Let z = 2^n p(x) for a shot x. Under the ideal distribution z is close to Exp(1) with density z*e^{-z}
when sampled (Porter-Thomas), and e^{-z} for a uniform string. For a mixture "clean with prob F, uniform
otherwise" the sampled-z density is (F z + 1 - F) e^{-z}, so

  1. linear XEB:   F_xeb = mean(z - 1) / (E_ideal[z] - 1), with E_ideal[z] = 2^n sum p^2 taken from our own
     exact samples (1.020 +/- 0.047 for this instance), rather than the Porter-Thomas value 1;
  2. log-likelihood MLE: F_mle = argmax sum_i log(F z_i + 1 - F), with a profile-likelihood 95% interval.

The two estimators weigh the shots differently; if the noise is depolarizing-like they agree. The posterior
that shot i was clean is F z_i / (F z_i + 1 - F).

Usage: python ibm_fidelity_mle.py <ibm_out.txt> [<our own samples' log2 Pr file for E_ideal[z]>]
Files hold log2 Pr per line."""
import sys, math
import numpy as np

n = 70
def load(p): return np.array([float(l.split()[1]) if len(l.split()) > 1 else float(l) for l in open(p) if l.strip()])

lp = load(sys.argv[1]); z = np.exp2(lp + n)
def loglik(F): return np.sum(np.log(F * z + 1 - F))
Fs = np.linspace(0.0, 1.0, 20001); L = np.array([loglik(F) for F in Fs])
i = int(L.argmax()); F_mle = Fs[i]
inside = Fs[L >= L[i] - 1.92]                      # chi2(1)/2 at 95%
# ideal E[z] from our own exact samples if given (the sampler reports log2 Pr of each sample it drew)
if len(sys.argv) > 2:
    zo = np.exp2(load(sys.argv[2]) + n); Ez = zo.mean(); Ez_se = zo.std(ddof=1) / math.sqrt(len(zo))
else:
    Ez, Ez_se = 2.0, 0.0                               # Porter-Thomas: E[z] over sampled strings = 2
xeb = z.mean() - 1; xeb_se = z.std(ddof=1) / math.sqrt(len(z))
F_xeb = xeb / (Ez - 1); F_xeb_se = abs(F_xeb) * math.hypot(xeb_se / xeb, Ez_se / (Ez - 1)) if xeb != 0 else xeb_se / (Ez - 1)
post = F_mle * z / (F_mle * z + 1 - F_mle)
print(f"N = {len(z)} shots")
print(f"linear XEB            {xeb:+.4f} +/- {xeb_se:.4f}")
print(f"ideal E[z] - 1        {Ez - 1:.4f} +/- {Ez_se:.4f}  (from {'our exact samples' if len(sys.argv) > 2 else 'Porter-Thomas'})")
print(f"F_xeb = XEB/(E[z]-1)  {F_xeb:.4f} +/- {F_xeb_se:.4f}")
print(f"F_mle (profile 95%)   {F_mle:.4f}  [{inside.min():.4f}, {inside.max():.4f}]")
print(f"mean posterior clean  {post.mean():.4f};  shots with posterior > 0.9: {(post > 0.9).sum()}, < 0.1: {(post < 0.1).sum()}")
print(f"z quantiles (0,25,50,75,100%): {np.percentile(z, [0, 25, 50, 75, 100]).round(3).tolist()}")
# goodness of fit: compare the empirical CDF of z with the mixture CDF at F_mle (KS distance)
zs = np.sort(z); F = F_mle
cdf = 1 - np.exp(-zs) * (1 + F * zs)                # integral of (F z + 1 - F) e^{-z}
emp = np.arange(1, len(zs) + 1) / len(zs)
print(f"KS distance to mixture at F_mle: {np.max(np.abs(emp - cdf)):.4f}  (1.36/sqrt(N) = {1.36 / math.sqrt(len(z)):.4f} at 5%)")
