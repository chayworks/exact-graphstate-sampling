"""Independent check of full 70-qubit samples: recompute Pr(sample) = |Z(w)|^2 / 2^70 with dp.cpp
(sliced amplitude engine, a different algorithm from gsamp's chain of conditionals) and compare with the
log2 Pr that gsamp reported. A string with one bit flipped is included as a negative control.

Usage: python xcheck70.py <nq70.qasm> <samples_nq70.txt> <k samples> <n sliced vertices>
"""
import sys, json, math, os, subprocess, numpy as np
from tableau import clifford_body
from amp import Amplitudes
from engine import Sliced, write_plan, write_jobs
from order import profile, cost

F, SF, k, s = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
n, T, tail = clifford_body(F); order = json.load(open("order_best2.json")); am = Amplitudes(n, T, tail, order)
steps_v = [st["v"] for st in am.steps]
assert steps_v == json.load(open("cloud/nq70_step_qubits.json")), "step order differs from the cloud plan"
# greedy slicing for this order (same rule as slicesearch.py)
S = []
def c(S):
    keep = [q for q in order if q not in S]; return cost(profile(am.A[np.ix_(keep, keep)], list(range(len(keep)))))
while len(S) < s:
    S.append(min((q for q in order if q not in S), key=lambda v: c(S + [v])))
sl = Sliced(am.A, order, S); maxr = write_plan("work/x70_plan.txt", sl.steps, sl.m)
print(f"sliced {S}: width {maxr}, buffers {2 * 16 * 2 ** maxr / 1e9:.1f} GB, 2^{s} slices", flush=True)
rows = [l.split() for l in open(SF) if l.strip()][:k]
def bvec(bits):
    b = [0] * n
    for kk, q in enumerate(steps_v): b[q] = int(bits[kk])
    return b
ws = [am.weights(bvec(bits)) for bits, _, _ in rows]
flip = list(rows[0][0]); flip[35] = "1" if flip[35] == "0" else "0"; ws.append(am.weights(bvec("".join(flip))))
write_jobs("work/x70_jobs.txt", sl, ws)
out = subprocess.run(["docker", "run", "--rm", "--memory", "9g", "-v", os.path.abspath(".").replace("\\", "/") + ":/w", "-w", "/w",
                      "-e", "OMP_NUM_THREADS=16", "gcc:14", "./dp", "work/x70_plan.txt", "work/x70_jobs.txt"],
                     capture_output=True, text=True, env={**os.environ, "MSYS_NO_PATHCONV": "1"})
assert out.returncode == 0, out.stderr
Z = [complex(*map(float, l.split())) for l in out.stdout.split("\n") if l.strip()]
assert len(Z) == len(ws)
for (bits, lp, seed), z in zip(rows, Z):
    print(f"seed {seed}: gsamp log2 Pr = {float(lp):.7f}   dp.cpp log2 |Z|^2/2^70 = {math.log2(abs(z) ** 2 / 2 ** n):.7f}")
print(f"negative control (sample 1, step-bit 35 flipped): dp.cpp log2 = {math.log2(abs(Z[-1]) ** 2 / 2 ** n):.7f}")
