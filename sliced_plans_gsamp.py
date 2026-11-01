"""Sliced plans for the two-table sampler gsamp.cpp (same slicing transform as sliced_plans.py: the sliced vertices lose
their edges, keep only the column x_v of their weight, and flip M_u[:, 1] of their neighbours when x_v = 1). All
assignments share the plan (plan.py on A'), so their 'prefix' dumps add up in 'resume'. Memory: two tables of
2^maxr complex64.  Usage: python sliced_plans_gsamp.py <gadget.qasm> <order.json> <out_prefix> [v1,v2,...]
Writes <out_prefix>_s<assignment>.txt and <out_prefix>_slices.json (K = first step of the pure-data suffix; for an
instance without ancillas K = number of steps). Checks that every sliced vertex and all its neighbours lie before K."""
import sys, json, itertools, numpy as np
from tableau import clifford_body
from amp import Amplitudes
from plan import plan
from engine import write_sampler_plan

if __name__ == "__main__":
    qasm, oj, out = sys.argv[1:4]; order = json.load(open(oj)); S = [int(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else []
    N, T, tail = clifford_body(qasm); am = Amplitudes(N, T, tail, order); A = am.A.copy(); M0 = [m.copy() for m in am.M]
    A2 = A.copy()
    for v in S: A2[v, :] = 0; A2[:, v] = 0
    am.A = A2; am.steps = plan(A2, order); maxr = max(max(s["r"], s["rn"]) for s in am.steps)
    print(f"N = {N}, width {maxr} with slices {S} (two tables of 2^{maxr} complex64 = {2 * 2 ** maxr * 8 / 1e9:.1f} GB)")
    sv = [s["v"] for s in am.steps]; K = next((i for i in range(len(sv)) if all(v < 70 for v in sv[i:])), len(sv))
    if K == 0: K = len(sv)
    pre = set(sv[:K])
    for v in S: assert v in pre and set(np.nonzero(A[v])[0]) <= pre, f"slice vertex {v} or a neighbour lies in the suffix"
    for bits in itertools.product([0, 1], repeat=len(S)):
        Ms = [m.copy() for m in M0]
        for v, xv in zip(S, bits):
            Ms[v][:, 1 - xv] = 0
            if xv:
                for u in np.nonzero(A[v])[0]: Ms[u][:, 1] *= -1
        am.M = Ms; write_sampler_plan(f"{out}_s{''.join(map(str, bits))}.txt", am)
    json.dump({"slices": S, "maxr": maxr, "K_suffix_start": K, "step_vertices": sv}, open(f"{out}_slices.json", "w"))
    print(f"wrote {2 ** max(len(S), 0)} plan(s); K = {K} of {len(sv)}")
