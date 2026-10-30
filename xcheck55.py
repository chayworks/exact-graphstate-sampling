"""Scale cross-check on real nq70 data (55-vertex induced subgraph of the best order):
gsamp's reported Pr(sample) (product of chain conditionals) vs |Z(w)|^2 / 2^m from dp.cpp (independent
amplitude path, with 6 sliced vertices). Agreement means two different algorithms give the same Born
probability on 55-qubit strings; disagreement on a negative control (one bit flipped in the weights) is required."""
import sys, json, types, subprocess, os, numpy as np
from tableau import clifford_body
from amp import Amplitudes
from engine import Sliced, write_plan, write_jobs
F = sys.argv[1]; m = 55
n, T, tail = clifford_body(F); order = json.load(open("order_best.json")); am = Amplitudes(n, T, tail, order)
keep = order[:m]; Ak = am.A[np.ix_(keep, keep)]; M = [am.M[q] for q in keep]
samples = [l.split() for l in open(sys.argv[2]) if l.strip()]
sl = Sliced(Ak, list(range(m)), [3, 17, 25, 31, 40, 48]); write_plan("work/x55_plan.txt", sl.steps, sl.m)
ws = []
for bits, lp, seed in samples:
    b = [int(c) for c in bits]; ws.append([(M[k][b[k], 0], M[k][b[k], 1]) for k in range(m)])
b = [int(c) for c in samples[0][0]]; b[10] ^= 1; ws.append([(M[k][b[k], 0], M[k][b[k], 1]) for k in range(m)])
write_jobs("work/x55_jobs.txt", sl, ws)
out = subprocess.run(["docker", "run", "--rm", "--memory", "4g", "-v", os.path.abspath(".").replace("\\", "/") + ":/w", "-w", "/w",
                      "-e", "OMP_NUM_THREADS=16", "gcc:14", "./dp", "work/x55_plan.txt", "work/x55_jobs.txt"],
                     capture_output=True, text=True, env={**os.environ, "MSYS_NO_PATHCONV": "1"})
assert out.returncode == 0, out.stderr
Z = [complex(*map(float, l.split())) for l in out.stdout.split("\n") if l.strip()]
for (bits, lp, seed), z in zip(samples, Z):
    print(f"seed {seed}: gsamp log2 Pr = {float(lp):.9f}   dp.cpp log2 |Z|^2/2^m = {np.log2(abs(z) ** 2 / 2 ** m):.9f}")
print(f"negative control (bit 10 of sample 1 flipped): dp.cpp log2 = {np.log2(abs(Z[-1]) ** 2 / 2 ** m):.9f}")
