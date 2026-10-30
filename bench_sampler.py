"""Throughput of gsamp on real nq70 data: induced subgraph on the first m vertices of the best order,
with the real measurement matrices. Not a sample of the full instance -- a speed benchmark only."""
import sys, json, math, types, numpy as np
from tableau import clifford_body
from amp import Amplitudes
from plan import plan
from engine import write_sampler_plan
F, m = sys.argv[1], int(sys.argv[2])
n, T, tail = clifford_body(F); order = json.load(open("order_best.json")); am = Amplitudes(n, T, tail, order)
keep = order[:m]; Ak = am.A[np.ix_(keep, keep)]
sub = types.SimpleNamespace(n=m, M=[am.M[q] for q in keep], steps=plan(Ak, list(range(m))))
maxr = write_sampler_plan("work/bench_sampler_plan.txt", sub)
visits = sum(2.0 ** (s["r"] + 2) for s in sub.steps)
print(f"m={m} maxr={maxr} visits per sample (2 passes) {visits:.3e}")
