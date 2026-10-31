"""Exact XEB of IBM's #151 hardware samples from our cloud scoring, with controls, and a string-by-string
comparison against #260's independently computed FP64 probabilities.

Usage: python ibm_xeb_analyze.py <ibm_out.txt> <uniform_out.txt> <reversed_out.txt> [<#260 ibm_scored_fp64.jsonl>]
Our files hold log2 Pr per line, in the order of the input strings (IBM index order)."""
import sys, json, math, random
import numpy as np

def load(p): return np.array([float(l) for l in open(p) if l.strip()])
def xeb(lp, n=70):
    z = np.exp2(lp + n); m = z.mean() - 1; se = z.std(ddof=1) / math.sqrt(len(z))
    rng = random.Random(424242); bs = []
    for _ in range(5000):
        s = [z[rng.randrange(len(z))] for _ in range(len(z))]; bs.append(np.mean(s) - 1)
    return m, se, np.percentile(bs, 2.5), np.percentile(bs, 97.5)

ibm, uni, rev = load(sys.argv[1]), load(sys.argv[2]), load(sys.argv[3])
for name, lp in [("IBM hardware samples", ibm), ("uniform random (control)", uni), ("IBM bit-reversed (control)", rev)]:
    m, se, lo, hi = xeb(lp); print(f"{name:28s} N={len(lp):5d}  linear XEB {m:+.4f} +/- {se:.4f}  (bootstrap 95% [{lo:+.4f}, {hi:+.4f}])")
if len(sys.argv) > 4:
    rec = [json.loads(l) for l in open(sys.argv[4])][1:]
    theirs = np.array([math.log2(r["probability"]) for r in sorted(rec, key=lambda r: r["index"])])
    k = min(len(theirs), len(ibm)); d = ibm[:k] - theirs[:k]
    print(f"vs #260 (independent code), first {k} IBM strings: max |log2 Pr difference| = {np.abs(d).max():.2e}, "
          f"median {np.median(np.abs(d)):.2e}; their XEB on the same strings {np.mean(np.exp2(theirs[:k] + 70)) - 1:+.4f}")
