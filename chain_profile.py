"""Position-resolved error profile from exact per-step conditionals (gsamp/gsamp2 chain mode).

For a string scored along the DP order, step k yields c_k = Pr(b_k | b_1..b_{k-1}), exact. For an ideal
sample E[2 c_k - 1] = 2 sum_b Pr(b|prefix)^2 - 1 averaged over prefixes, which is >= 0 and is zero exactly
when the prefix marginal is uniform (early steps, where the cut rank equals the prefix length). For a shot
whose bit k is uncorrelated with the circuit, E[2 c_k - 1] = 0. The ratio

    f_k = mean_hardware(2 c_k - 1) / mean_ideal(2 c_k - 1)

is therefore a per-step fidelity-like quantity in [0, 1], conditional on the hardware prefix being what it
is. It depends on the order; running the reverse order (same width) covers the qubits that sit early in the
forward order. Steps where the ideal mean is below `--min-ideal` carry no information and are skipped.

Usage: python chain_profile.py <hw_chain.txt> <ideal_chain.txt> <step_qubits.json> [--min-ideal 0.05] [--out profile.json]
Chain files: per line "log2Pr c_1 ... c_n". Bootstrap SEs over strings."""
import sys, json, math
import numpy as np

def load(p):
    rows = [[float(x) for x in l.split()] for l in open(p) if l.strip()]
    lp = np.array([r[0] for r in rows]); c = np.array([r[1:] for r in rows]); return lp, c

def main(a):
    mi = float(a[a.index("--min-ideal") + 1]) if "--min-ideal" in a else 0.05
    out = a[a.index("--out") + 1] if "--out" in a else None
    lp_h, ch = load(a[0]); lp_i, ci = load(a[1]); qub = json.load(open(a[2]))
    n = ch.shape[1]; assert ci.shape[1] == n == len(qub)
    xh = 2 * ch - 1; xi = 2 * ci - 1                               # per-string, per-step
    mh, mi_ = xh.mean(0), xi.mean(0)
    rng = np.random.default_rng(7); B = 2000
    bh = np.array([xh[rng.integers(0, len(xh), len(xh))].mean(0) for _ in range(B)])
    bi = np.array([xi[rng.integers(0, len(xi), len(xi))].mean(0) for _ in range(B)])
    print(f"hardware strings {len(xh)}, ideal strings {len(xi)}; XEB check: hw {np.exp2(lp_h + n).mean() - 1:+.4f}, ideal {np.exp2(lp_i + n).mean() - 1:+.4f}")
    print(f"{'step':>4} {'qubit':>5} {'ideal<2c-1>':>12} {'hw<2c-1>':>10} {'f_k':>7} {'se':>6}")
    res = []
    for k in range(n):
        if mi_[k] < mi:
            res.append(dict(step=k + 1, qubit=qub[k], ideal=mi_[k], hw=mh[k], f=None)); continue
        f = mh[k] / mi_[k]; fb = bh[:, k] / bi[:, k]; se = fb.std()
        res.append(dict(step=k + 1, qubit=qub[k], ideal=mi_[k], hw=mh[k], f=f, se=se))
        print(f"{k + 1:>4} {qub[k]:>5} {mi_[k]:>12.4f} {mh[k]:>10.4f} {f:>7.3f} {se:>6.3f}")
    # where the information lives: spread of the conditionals by step, and the share of var(log2 Pr) over strings
    # that the last k steps account for (ideal samples); prefix marginals are uniform while the spread is ~0
    sdi, sdh = xi.std(0), xh.std(0); li = np.log2(np.clip(ci, 1e-300, None)); tot = li.sum(1).var()
    share = [li[:, n - k:].sum(1).var() / tot for k in (4, 8, 12, 16)]
    first = next((k + 1 for k in range(n) if sdi[k] > 1e-3), None)
    Lh = np.log2(ch); print(f"hardware: var(log2 Pr) over strings {lp_h.var(ddof=1):.3f}; variance of the last 4 / 8 conditionals' log2 sum "
                            f"{Lh[:, -4:].sum(1).var(ddof=1):.3f} / {Lh[:, -8:].sum(1).var(ddof=1):.3f}")
    print(f"ideal: first step with std(2c-1) > 1e-3: {first}; share of var(log2 Pr) from the last 4/8/12/16 steps: "
          + "/".join(f"{v:.3f}" for v in share))
    print("step : std(2c-1) ideal | hardware   (steps 1..n)")
    print(" ".join(f"{k + 1}:{sdi[k]:.1e}|{sdh[k]:.1e}" for k in range(n)))
    inf = [r for r in res if r["f"] is not None]
    print(f"{len(inf)} informative steps (ideal mean >= {mi}); first informative step {inf[0]['step'] if inf else None}")
    if inf:
        f = np.array([r["f"] for r in inf]); print(f"mean f over informative steps {f.mean():.3f}; min {f.min():.3f} at step {inf[int(f.argmin())]['step']} (qubit {inf[int(f.argmin())]['qubit']})")
    if out: json.dump(res, open(out, "w"), indent=1)

if __name__ == "__main__": main(sys.argv[1:])
