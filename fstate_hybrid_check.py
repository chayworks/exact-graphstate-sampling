"""Full-scale validation of the hybrid rejection sampler on the 70-qubit fstate instance, where the exact sampler is
the ground truth. Input: the driver's sample file (ACC/rej lines with log2 Pr and z) produced with the width-30
in-place plan, 50 bits proposed, 20 bits drawn exactly. Checks:
  1. the accepted samples' linear XEB against their own exact probabilities equals the exact samples' value
     (1.020 +/- 0.047 over 1,012 exact samples; a biased sampler shifts this);
  2. the prefix marginal z (an average of 2^20 output probabilities in uniform units) is concentrated near 1, so
     the cap excludes nothing (no proposal above the cap);
  3. the reported log2 Pr of a few accepted samples is reproduced by prob mode of the same binary (strings re-scored).
Usage: python fstate_hybrid_check.py work/fstate_hybrid_samples.txt results/samples_nq70.txt [prob_rescored.txt]"""
import sys, math, numpy as np
n = 70
rows = [l.split() for l in open(sys.argv[1]) if l.startswith(("ACC", "rej"))]
acc = np.array([float(r[2]) for r in rows if r[0] == "ACC"]); zs = np.array([float(r[3][2:]) for r in rows]); above = sum("Z_ABOVE_CAP" in " ".join(r) for r in rows)
ex = np.array([float(l.split()[1]) for l in open(sys.argv[2]) if len(l.split()) >= 2 and len(l.split()[0]) == n])
def xeb(lp): z = np.exp2(lp + n); return z.mean() - 1, z.std(ddof=1) / math.sqrt(len(z))
xa, sa = xeb(acc); xe, se = xeb(ex)
print(f"proposals {len(rows)}, accepted {len(acc)} ({100 * len(acc) / len(rows):.1f}%); z: mean {zs.mean():.5f}, sd {zs.std():.5f}, min {zs.min():.4f}, max {zs.max():.4f}; above cap: {above}")
print(f"linear XEB of accepted hybrid samples vs own exact probabilities: {xa:+.3f} +/- {sa:.3f}; exact sampler (1,012 samples): {xe:+.3f} +/- {se:.3f}; difference {(xa - xe) / math.hypot(sa, se):+.2f} sigma")
# the mean of log2 Pr is a second distributional statistic (the entropy-like moment)
print(f"mean log2 Pr: hybrid {acc.mean():.3f} +/- {acc.std(ddof=1) / math.sqrt(len(acc)):.3f}; exact {ex.mean():.3f} +/- {ex.std(ddof=1) / math.sqrt(len(ex)):.3f}")
if len(sys.argv) > 3:
    re = [float(l.split()[0]) for l in open(sys.argv[3]) if l.strip()]; k = len(re)
    rep = [float(r[2]) for r in rows if r[0] == "ACC"][:k]
    print(f"re-scored {k} accepted samples with prob mode: max |log2 Pr difference| = {max(abs(a - b) for a, b in zip(re, rep)):.1e}")
