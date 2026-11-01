"""Analysis of cloud run 7 (doped instance). Inputs are the run's own output files.
  shots_scored.txt : "<70-bit string> <log2 Pr_doped> <seconds>" per line; the first 16 strings are uniform controls,
                     the rest are IBM's published doped shots (Qiskit order, as given).
  samples.txt      : "ACC|rej <70-bit string> <log2 Pr_doped> z=<marginal in uniform units> <seconds>[ Z_ABOVE_CAP]"
Reports: linear XEB of the hardware shots and of the controls with bootstrap intervals, the fidelity estimators
(XEB normalised by the samples' own 2^70 sum p^2, and the mixture MLE), the hybrid sampler's acceptance, max z and
self-consistency XEB, and timing. Usage: python doped_analyze.py <run_dir> [n_controls=16]"""
import sys, os, math, numpy as np
from scipy.stats import chi2
from scipy.optimize import minimize_scalar
d = sys.argv[1]; nctrl = int(sys.argv[2]) if len(sys.argv) > 2 else 16; n = 70
rows = [l.split() for l in open(os.path.join(d, "shots_scored.txt")) if l.strip() and not l.startswith("#")]
lp = np.array([float(r[1]) for r in rows]); sec = np.array([float(r[2]) for r in rows]); z = np.exp2(lp + n)
zc, zh = z[:nctrl], z[nctrl:]
def xeb(zz, B=2000):
    if len(zz) == 0: return float("nan"), float("nan"), (float("nan"), float("nan"))
    rng = np.random.default_rng(0); bs = [rng.choice(zz, len(zz)).mean() - 1 for _ in range(B)]
    return zz.mean() - 1, zz.std(ddof=1) / math.sqrt(len(zz)), (np.percentile(bs, 2.5), np.percentile(bs, 97.5))
print(f"scored strings: {len(rows)} ({nctrl} uniform controls + {len(zh)} IBM doped shots); mean time per pass {sec.mean():.0f} s on this node")
xh, seh, cih = xeb(zh); xc, sec_, cic = xeb(zc)
print(f"IBM doped shots      N={len(zh):5d}  linear XEB {xh:+.4f} +/- {seh:.4f}  (bootstrap 95% [{cih[0]:+.3f}, {cih[1]:+.3f}])")
print(f"uniform controls     N={len(zc):5d}  linear XEB {xc:+.4f} +/- {sec_:.4f}  (bootstrap 95% [{cic[0]:+.3f}, {cic[1]:+.3f}])")
# samples
acc = []; zs = []; nprop = 0; above = 0; sacc = []
sp = os.path.join(d, "samples.txt")
if os.path.exists(sp):
    for l in open(sp):
        if l.startswith(("ACC", "rej")):
            f = l.split(); nprop += 1; zz = float(f[3][2:]); zs.append(zz); above += "Z_ABOVE_CAP" in l
            if f[0] == "ACC": acc.append(float(f[2])); sacc.append(float(f[4].rstrip("s")))
    if nprop:
        za = np.exp2(np.array(acc) + n) if acc else np.array([])
        print(f"hybrid sampler: {len(acc)} accepted of {nprop} proposals ({100 * len(acc) / nprop:.1f}%); prefix marginal z: mean {np.mean(zs):.4f}, "
              f"sd {np.std(zs):.4f}, min {min(zs):.4f}, max {max(zs):.4f}; proposals above the cap: {above}; mean time per proposal {np.mean(sacc) if sacc else float('nan'):.0f} s")
        if len(za) >= 2:
            xs, ses, cis = xeb(za); print(f"accepted samples: self-consistency linear XEB {xs:+.3f} +/- {ses:.3f} (ideal: 2^70 sum p^2 - 1, about +1)")
            Ez = za.mean()
            # fidelity of the hardware shots: XEB ratio and mixture maximum likelihood, as for the fstate instance
            F_xeb = xh / (Ez - 1) if Ez > 1 else float("nan")
            nll = lambda F: -np.sum(np.log(np.clip(F * zh + 1 - F, 1e-300, None)))
            r = minimize_scalar(nll, bounds=(0, 1), method="bounded"); Fm = r.x
            lo = minimize_scalar(lambda F: abs(nll(F) - nll(Fm) - 1.92), bounds=(0, Fm), method="bounded").x if Fm > 0 else 0.0
            hi = minimize_scalar(lambda F: abs(nll(F) - nll(Fm) - 1.92), bounds=(Fm, 1), method="bounded").x if Fm < 1 else 1.0
            print(f"fidelity of IBM's doped shots: F_xeb = XEB / (E[z] - 1) = {F_xeb:.3f} (E[z] - 1 = {Ez - 1:.3f} from {len(za)} samples, large uncertainty); "
                  f"F_mle (clean-or-uniform mixture, profile 95%) = {Fm:.3f} [{lo:.3f}, {hi:.3f}]; IBM reports 0.26-0.32")
