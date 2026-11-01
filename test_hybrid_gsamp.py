"""The two-table sampler's new modes (hybrid, prefix, resume) against the in-place sampler gsamp2 on the hy21 T-gadget
case (files from test_hybrid.py) and against dense statevectors. Both engines use the same vertex order, so hybrid
conditionals and log2 Pr must agree; sliced prefix+resume (vertex 8, K = 10) must reproduce the unsliced result."""
import os, json, numpy as np
from tableau import clifford_body
from amp import Amplitudes
from engine import write_sampler_plan
from plan import plan
from test_gsamp import docker
assert os.path.exists("work/hy21.qasm"), "run test_hybrid.py first"
steps = json.load(open("work/hy21_steps.json")); N = len(steps); n = 5; K = max(i for i, v in enumerate(steps) if v < 3) + 1
nn, T, tl = clifford_body("work/hy21.qasm"); order = steps
am = Amplitudes(N, T, tl, order); write_sampler_plan("work/hy21g.txt", am); A = am.A.copy(); M0 = [m.copy() for m in am.M]
vS = 8; A2 = A.copy(); A2[vS, :] = 0; A2[:, vS] = 0; am.A = A2; am.steps = plan(A2, order)
for xv in (0, 1):
    Ms = [m.copy() for m in M0]; Ms[vS][:, 1 - xv] = 0
    if xv:
        for u in np.nonzero(A[vS])[0]: Ms[u][:, 1] *= -1
    am.M = Ms; write_sampler_plan(f"work/hy21g_sl{xv}.txt", am)
out = docker("./gsamp2 work/hy21.txt hybrid work/hy21_full.txt 5; echo SEP; ./gsamp work/hy21g.txt hybrid work/hy21_full.txt 5; echo SEP; "
             f"./gsamp work/hy21g_sl0.txt prefix work/hy21_full.txt {K} work/hy21g_d0; ./gsamp work/hy21g_sl1.txt prefix work/hy21_full.txt {K} work/hy21g_d1; "
             f"./gsamp work/hy21g_sl0.txt resume work/hy21_full.txt {K} 3 work/hy21g_d0 work/hy21g_d1").split("SEP")
h2 = [l.split() for l in out[0].strip().split("\n")]; h1 = [l.split() for l in out[1].strip().split("\n")]
res = [l.split() for l in out[2].strip().split("\n") if len(l.split()[0]) == N]
d_lp = max(abs(float(a[1]) - float(b[1])) for a, b in zip(h1, h2)); d_c = max(abs(float(x) - float(y)) for a, b in zip(h1, h2) for x, y in zip(a[2:], b[2:]))
d_s = max(abs(float(r[1]) - float(b[1])) for r, b in zip(res, h2)); d_p = max(abs(float(r[2]) - sum(np.log2(float(c)) for c in b[2:2 + K])) for r, b in zip(res, h2))
# suffix sampling through the two-table engine: 20,000 draws of the last 2 bits for a fixed prefix vs the exact conditional
from test_hybrid import p_exact, line, pos, m
prefix = [1, 0, 1, 1, 0]; l = line(prefix, set(range(n - m, n))); open("work/hy21g_suf.txt", "w").write((l + "\n") * 20000)
draws = docker("./gsamp work/hy21g.txt hybrid work/hy21g_suf.txt 7").strip().split("\n"); cnt = np.zeros(2 ** m)
for d in draws:
    bits = d.split()[0]; cnt[int("".join(bits[pos[q]] for q in range(n - m, n)), 2)] += 1
idx = [i for i in range(2 ** n) if all(((i >> (n - 1 - q)) & 1) == prefix[q] for q in range(n - m))]; pc = np.array([p_exact[i] for i in idx]); pc /= pc.sum()
tv = 0.5 * np.abs(cnt / cnt.sum() - pc).sum()
print(f"gsamp vs gsamp2 hybrid: log2 Pr {d_lp:.1e}, conditionals {d_c:.1e}; sliced prefix+resume: log2 Pr {d_s:.1e}, prefix {d_p:.1e}; suffix draws TV {tv:.4f}")
assert d_lp < 1e-6 and d_c < 1e-6 and d_s < 1e-6 and d_p < 1e-6 and tv < 0.01
