"""Which angle convention reproduces IQPopt? Score the samples IQPopt's own simulator drew from each trained 28-qubit
model under exp(i f theta X) for f in {1, 1/2, 2}: the true law maximises the expected log-likelihood (Gibbs), so the
right f gives the highest mean log2 p, and mean 2^n p should match sum p^2 computed exactly from the model.
Usage: python convention_check.py <IQP-NISQ results dir>"""
import sys, glob, os, numpy as np
from pathlib import Path
from heavyhex import gate_list
from iqp_exact import IQPModel

B = Path(sys.argv[1])
sims = sorted(B.glob("simulation_results_iqp/8N_*/*/*/rank_*/simulated_samples_20260128.npy"))
print(len(sims), "simulated sample files")
for sp in sims[:6]:
    rel = sp.parent.relative_to(B / "simulation_results_iqp")
    pp = B / "trained_params_iqp" / rel / "params.npy"
    if not pp.exists(): print("no params for", rel); continue
    th = np.load(pp, allow_pickle=False); S = np.load(sp, allow_pickle=False); N = S.shape[1]
    g = gate_list(N); L = len(th) // len(g); assert len(th) == L * len(g), (len(th), len(g))
    res = []
    for f in (1, 0.5, 2):
        m = IQPModel(N, g * L, th * f); lp = np.mean([np.log2(max(m.prob(s), 1e-300)) for s in S[:300]])
        res.append((f, round(lp, 3)))
    print(rel.as_posix(), "N", N, "layers", L, "samples", S.shape, S.dtype, "mean log2 p by f:", res)

# bit order and calibration on the first model: IQPopt's samples vs our exact samples from the same parameters
sp = sims[0]; rel = sp.parent.relative_to(B / "simulation_results_iqp")
th = np.load(B / "trained_params_iqp" / rel / "params.npy", allow_pickle=False); S = np.load(sp, allow_pickle=False)
N = S.shape[1]; m = IQPModel(N, gate_list(N), th); rng = np.random.default_rng(3)
ours = m.sample(rng, 2000)
lp_ours = np.array([np.log2(p) for _, p in ours]); lp_sim = np.array([np.log2(m.prob(s)) for s in S])
lp_rev = np.array([np.log2(max(m.prob(s[::-1]), 1e-300)) for s in S])
mean_ours = np.mean([p for _, p in ours]) * 2 ** N; mean_sim = np.mean([m.prob(s) for s in S]) * 2 ** N
print(f"{rel.as_posix()}: mean log2 p  ours {lp_ours.mean():.3f} +/- {lp_ours.std() / np.sqrt(len(lp_ours)):.3f}, "
      f"IQPopt {lp_sim.mean():.3f} +/- {lp_sim.std() / np.sqrt(len(lp_sim)):.3f}, IQPopt reversed bits {lp_rev.mean():.3f}")
print(f"mean 2^n p: ours {mean_ours:.1f}, IQPopt {mean_sim:.1f}; per-qubit means max |diff| "
      f"{np.max(np.abs(np.mean([y for y, _ in ours], 0) - S.mean(0))):.3f}")
