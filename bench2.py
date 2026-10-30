"""Steady-state throughput: one nq70 amplitude, first `k` slices of an s-vertex slicing, one dp run."""
import sys, json, math, random, itertools, numpy as np
from tableau import clifford_body
from amp import Amplitudes
from engine import Sliced, write_plan
F, s, k = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
S = json.load(open("slice_S_long.json"))[:s]
n, T, tail = clifford_body(F); order = json.load(open("order_best.json"))
am = Amplitudes(n, T, tail, order); sl = Sliced(am.A, order, S)
maxr = write_plan("work/nq70_plan.txt", sl.steps, sl.m)
rng = random.Random(0); b = [rng.randint(0, 1) for _ in range(n)]
sls = list(itertools.islice(sl.slices(am.weights(b)), k))
with open("work/nq70_jobk.txt", "w") as f:
    f.write(f"1\n{len(sls)}\n")
    for pre, ww in sls:
        f.write("%.17g %.17g\n" % (pre.real, pre.imag))
        for a, c in ww: f.write("%.17g %.17g %.17g %.17g\n" % (a.real, a.imag, c.real, c.imag))
v = sum(2.0 ** (st["r"] + 1) for st in sl.steps)
print(f"s={s} maxr={maxr} slices run={len(sls)} of {2**s}; visits per slice {v:.3e}; total visits this run {v*len(sls):.3e}")
