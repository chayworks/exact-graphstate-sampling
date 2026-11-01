"""What the 2,051 exact probabilities published by tracker entry #242 (Manabe, Gu, Pan; Zenodo 10.5281/zenodo.21912448,
CC BY 4.0) say about IBM's doped run and about the circuit's ideal output distribution. Each batch holds 256 amplitudes:
the IBM shot and its 255 siblings that differ on qubits q0..q7. Writes nothing; prints
  (1) a Porter-Thomas test of the 523,005 sibling strings (z = 2^70 p should be Exp(1) for a string not sampled from
      the circuit), which justifies E_ideal[z] = 2 in the fidelity estimators of ibm_fidelity_mle.py;
  (2) the IBM shots' mean z and their rank among their siblings;
  (3) a local 8-qubit XEB from the conditional distribution of q0..q7 given the other 62 bits (ideal 1, uniform 0, mixture F).
Usage: python doped_tn242_analysis.py <amplitude_batches.npz>"""
import sys, math, numpy as np
z = np.load(sys.argv[1]); rf = float(z["recovery_factor"]); rv = z["raw_vectors"]; sel = z["selected_vector_indices"]
P = np.abs(rv.astype(np.complex128) * rf) ** 2 * 2.0 ** 70
mask = np.ones(P.shape, bool); mask[np.arange(len(sel)), sel] = False
zs = P[mask].ravel(); print(f"#242 batches: {P.shape[0]} x {P.shape[1]} amplitudes; non-shot siblings N = {zs.size}")
print(f"mean z (ideal uniform-string expectation 1): {zs.mean():.5f} +/- {zs.std(ddof=1) / math.sqrt(zs.size):.5f};  mean z^2 (PT: 2): {(zs ** 2).mean():.4f}")
s = np.sort(zs); emp = np.arange(1, s.size + 1) / s.size; ks = np.max(np.abs(emp - (1 - np.exp(-s))))
print(f"KS distance to Exp(1): {ks:.5f}  (1.36/sqrt(N) = {1.36 / math.sqrt(s.size):.5f} at 5%; siblings of one shot are not independent, so this is indicative)")
q = [.01, .1, .5, .9, .99, .999]
print("quantiles z (1,10,50,90,99,99.9%):", np.percentile(zs, [100 * x for x in q]).round(4).tolist(), " Exp(1):", (-np.log(1 - np.array(q))).round(4).tolist())
bm = P[mask].reshape(P.shape[0], P.shape[1] - 1).mean(1)
print(f"per-batch mean z over 255 siblings: mean {bm.mean():.4f}, sd {bm.std(ddof=1):.4f} (iid Exp(1) would give sd {1 / math.sqrt(255):.4f})")
sh = P[np.arange(len(sel)), sel]; rank = (P < sh[:, None]).sum(1) / 255.0
print(f"IBM shots: mean z {sh.mean():.4f} +/- {sh.std(ddof=1) / math.sqrt(len(sh)):.4f}; mean rank of the shot's z among its 256 siblings {rank.mean():.4f} (uniform sampling 0.5, ideal 2/3)")
loc = 256 * sh / P.sum(1) - 1
print(f"local 8-qubit XEB: mean of 256 p(shot)/sum_siblings p - 1 = {loc.mean():.4f} +/- {loc.std(ddof=1) / math.sqrt(len(sh)):.4f}")
