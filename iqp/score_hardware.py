"""Exact scoring of IQP-NISQ's IBM hardware samples (arXiv:2511.05267, Zenodo 18983816) and exact classical sampling
of the same trained models.

For every hardware sample file with trained parameters: build the model (gate list of src/utils/hardware.py, angle
convention exp(i theta X) checked in convention_check.py), draw K exact samples with our sampler (timed), and compute
  - mean log2 p of the hardware shots and of exact samples (the ideal value), and of uniform strings (-N);
  - the normalised linear cross-entropy F = (2^N <p>_hw - 1) / (2^N <p>_ideal - 1), <p>_ideal = sum_y p(y)^2
    estimated from the exact samples (an unbiased estimate of the collision probability).
Usage: python score_hardware.py <IQP-NISQ results dir> [K]"""
import sys, time, numpy as np
from pathlib import Path
from collections import defaultdict
from heavyhex import gate_list
from iqp_exact import IQPModel

B = Path(sys.argv[1]); K = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
files = sorted(B.glob("evaluation_results_iqp/*/*/*/rank_*/samples_*.npy"))
rows = []; rng = np.random.default_rng(20261106)
for f in files:
    rel = f.parent.relative_to(B / "evaluation_results_iqp"); pp = B / "trained_params_iqp" / rel / "params.npy"
    if not pp.exists(): continue
    th = np.load(pp, allow_pickle=False); S = np.load(f, allow_pickle=False); N = S.shape[1]
    g = gate_list(N); L = len(th) // len(g); assert len(th) == L * len(g), (rel, len(th), len(g))
    m = IQPModel(N, g * L, th)
    t0 = time.time(); ex = m.sample(rng, K); tsamp = (time.time() - t0) / K
    pe = np.array([p for _, p in ex]); ph = np.array([m.prob(s) for s in S])
    coll = pe.mean(); scale = 2.0 ** N
    F = (scale * ph.mean() - 1) / (scale * coll - 1)
    # standard error of F from both sample means (delta method)
    sF = np.hypot(scale * ph.std() / np.sqrt(len(ph)), F * scale * pe.std() / np.sqrt(len(pe))) / (scale * coll - 1)
    rows.append(dict(rel=rel.as_posix(), N=N, shots=len(S), width=m.width, lp_hw=np.log2(ph).mean(),
                     lp_ideal=np.log2(pe).mean(), F=F, sF=sF, t=tsamp))
    r = rows[-1]
    print(f"{r['rel']:45s} N={N:3d} shots={len(S):5d} w={m.width} log2p hw {r['lp_hw']:8.2f} ideal {r['lp_ideal']:8.2f} "
          f"uniform {-N:5d}  F_XEB {F:6.3f} +/- {sF:.3f}  exact sampler {1e3 * tsamp:.1f} ms/sample", flush=True)
print("\nby size (mean over models; F_XEB with the spread across models):")
by = defaultdict(list)
for r in rows: by[r["N"]].append(r)
for N in sorted(by):
    R = by[N]; Fs = np.array([r["F"] for r in R])
    print(f"N = {N:3d}: {len(R):2d} models, width {max(r['width'] for r in R)}, F_XEB mean {Fs.mean():.3f} "
          f"(models {Fs.min():.3f} .. {Fs.max():.3f}); log2 p hw {np.mean([r['lp_hw'] for r in R]):.1f} vs ideal "
          f"{np.mean([r['lp_ideal'] for r in R]):.1f} vs uniform {-N}; exact sampling {1e3 * np.mean([r['t'] for r in R]):.1f} ms/sample")
