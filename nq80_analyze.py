"""Report on a hybrid-sampler run of an undoped n-qubit instance (cloud run 12: 80 qubits), from the driver's sample file
and the optional cross-check file. Checks: acceptance and the spread of z (the prefix marginal in uniform units, which
averages 2^m output probabilities and must stay within a fraction of a percent of 1 with no proposal above the cap);
linear XEB of the accepted samples against their own exact probabilities (an ideal sampler of a Porter-Thomas-like
output gives about 1; a uniform sampler gives 0); mean log2 Pr; time per proposal; and the agreement of the re-scored
samples (independent slicing) with the sampler's reported log2 Pr.
Usage: python nq80_analyze.py <samples.txt> [<xcheck.txt>] [--n 80]"""
import sys, math, numpy as np
n = 80
if "--n" in sys.argv: i = sys.argv.index("--n"); n = int(sys.argv[i + 1]); del sys.argv[i:i + 2]
rows = [l.split() for l in open(sys.argv[1]) if l.startswith(("ACC", "rej"))]
acc = [r for r in rows if r[0] == "ACC"]; z = np.array([float(r[3][2:]) for r in rows]); above = sum("Z_ABOVE_CAP" in " ".join(r) for r in rows)
lp = np.array([float(r[2]) for r in acc]); secs = np.array([float(r[4].rstrip("s")) for r in rows])
print(f"{len(rows)} proposals, {len(acc)} accepted ({100 * len(acc) / len(rows):.1f}%); above cap: {above}")
print(f"z: mean {z.mean():.5f}, sd {z.std(ddof=1) if len(z) > 1 else 0:.5f}, min {z.min():.4f}, max {z.max():.4f}  (2^-10 = {2 ** -10:.5f} expected spread for m = 20)")
zz = np.exp2(lp + n)
print(f"linear XEB of accepted samples vs own exact probabilities: {zz.mean() - 1:+.3f} +/- {zz.std(ddof=1) / math.sqrt(len(zz)) if len(zz) > 1 else 0:.3f}")
print(f"mean log2 Pr: {lp.mean():.3f} +/- {lp.std(ddof=1) / math.sqrt(len(lp)) if len(lp) > 1 else 0:.3f}; min {lp.min():.3f}, max {lp.max():.3f}; -n = {-n}")
print(f"time per proposal: mean {secs.mean():.1f} s, min {secs.min():.1f}, max {secs.max():.1f}")
if len(sys.argv) > 2:
    xc = {l.split()[0]: float(l.split()[1]) for l in open(sys.argv[2]) if len(l.split()) >= 2}
    rep = {r[1]: float(r[2]) for r in acc}; d = [(s, xc[s], rep[s]) for s in xc if s in rep]
    for s, a, b in d: print(f"cross-check {s[:16]}...: pair slicing {a:.9f}, sampler {b:.9f}, diff {a - b:.1e}")
    if d: print(f"max |diff| over {len(d)} re-scored samples: {max(abs(a - b) for _, a, b in d):.1e}")
