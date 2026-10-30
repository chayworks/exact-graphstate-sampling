"""Time one slice of one nq70 amplitude (real circuit, real weights) for a given slicing depth."""
import sys, json, math, time, numpy as np, random
from tableau import clifford_body
from amp import Amplitudes
from engine import Sliced, write_plan
from order import profile, cost
F = sys.argv[1]; s = int(sys.argv[2]); S = json.load(open("slice_S_long.json"))[:s]
n, T, tail = clifford_body(F); order = json.load(open("order_best.json"))
am = Amplitudes(n, T, tail, order)
assert (am.A == np.load("A_nq70.npy")).all()
sl = Sliced(am.A, order, S); maxr = write_plan("work/nq70_plan.txt", sl.steps, sl.m)
rng = random.Random(0); b = [rng.randint(0, 1) for _ in range(n)]
pre, ww = next(sl.slices(am.weights(b)))
with open("work/nq70_job1.txt", "w") as f:
    f.write("1\n1\n%.17g %.17g\n" % (pre.real, pre.imag))
    for a, c in ww: f.write("%.17g %.17g %.17g %.17g\n" % (a.real, a.imag, c.real, c.imag))
pr = profile(sl.A if hasattr(sl, 'Ak') else am.A[np.ix_(sl.keep, sl.keep)], list(range(sl.m)))
print(f"s={s} maxr={maxr} states visited per slice 2^{math.log2(sum(2.0 ** (st['r'] + 1) for st in sl.steps)):.2f}", flush=True)
