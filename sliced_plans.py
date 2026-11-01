"""Sliced in-place plans. Slicing vertex v fixes x_v: the Ising sum splits into the two values of x_v, each a sum on the
graph with v's edges removed, with weights w_v(x_v) kept and the phases (-1)^{x_v x_u} folded into M_u[:, 1] for every
neighbour u. All 2^{|S|} slice assignments share the same graph A' and therefore the same plan structure (addresses,
step sizes); only the M matrices differ, so their dumped tables add up (gsamp2 prefix/resume).
Usage: python sliced_plans.py <gadget.qasm> <order.json> <out_prefix> [v1,v2,...]
With no slice list: scores candidate vertices (those whose removal lowers the width most, searched among the vertices
processed around the widest cuts) and prints the best; with a list: writes <out_prefix>_s<assignment>.txt plans and
<out_prefix>_slices.json."""
import sys, json, itertools, numpy as np
from tableau import clifford_body, gf2_rank
from amp import Amplitudes
from inplace import plan_inplace, write_inplace_plan

def profile(A, order, ks=None):
    n = len(A); ks = ks or range(1, n); return {k: gf2_rank(A[np.ix_(order[k:], order[:k])]) for k in ks}

if __name__ == "__main__":
    qasm, oj, out = sys.argv[1:4]; order = json.load(open(oj))
    N, T, tail = clifford_body(qasm); am = Amplitudes(N, T, tail, order); A = am.A.copy(); M = [m.copy() for m in am.M]
    pr = [s["r"] for s in am.steps]; w = max(pr); hi = [k for k, r in enumerate(pr) if r >= w - 1]
    print(f"N = {N}, width {w}, widest cuts at steps {hi[0]}..{hi[-1]} ({len(hi)} steps)")
    if len(sys.argv) < 5:
        cand = sorted(set(order[max(0, hi[0] - 40):min(N, hi[-1] + 40)])); ks = list(range(max(1, hi[0] - 5), min(N, hi[-1] + 6)))
        best = []
        for v in cand:
            A2 = A.copy(); A2[v, :] = 0; A2[:, v] = 0; p2 = profile(A2, order, ks); best.append((max(p2.values()), v))
        best.sort(); print("best single slices (new max over the widest region, vertex):", best[:12])
        # greedy second slice on top of the best first
        v1 = best[0][1]; A2 = A.copy(); A2[v1, :] = 0; A2[:, v1] = 0; best2 = []
        for v in cand:
            if v == v1: continue
            A3 = A2.copy(); A3[v, :] = 0; A3[:, v] = 0; p3 = profile(A3, order, ks); best2.append((max(p3.values()), v))
        best2.sort(); print(f"best second slice with {v1}:", best2[:8])
        raise SystemExit
    S = [int(x) for x in sys.argv[4].split(",")]
    A2 = A.copy()
    for v in S: A2[v, :] = 0; A2[:, v] = 0
    st = plan_inplace(A2, order, max(1, max(profile(A2, order).values()))); R = max(max(s["r"], s["rn"]) for s in st) if "rn" in st[0] else None
    R = max(len(s["keep"]) for s in st) + 1
    full = profile(A2, order); R = max(full.values()); print(f"sliced graph: width {R} (was {w}); plan with R = {R}")
    st = plan_inplace(A2, order, R)
    for bits in itertools.product([0, 1], repeat=len(S)):
        Ms = [m.copy() for m in M]
        for v, xv in zip(S, bits):
            Ms[v][:, 1 - xv] = 0                                   # only x_v = xv contributes
            if xv:
                for u in np.nonzero(A[v])[0]: Ms[u][:, 1] *= -1     # (-1)^{x_v x_u} for every neighbour u
        tag = "".join(map(str, bits)); write_inplace_plan(f"{out}_s{tag}.txt", st, R, Ms, N)
    K = next(i for i, s in enumerate(st) if all(t["v"] < 70 for t in st[i:]))     # first step of the pure-data suffix
    pre = set(s["v"] for s in st[:K])
    for v in S: assert v in pre and set(np.nonzero(A[v])[0]) <= pre, f"slice vertex {v} or a neighbour lies in the suffix; tables could not be summed at K"
    json.dump({"slices": S, "R": R, "width_unsliced": w, "K_suffix_start": K, "step_vertices": [s["v"] for s in st]}, open(f"{out}_slices.json", "w"))
    print(f"wrote {2 ** len(S)} slice plans; suffix (data only) starts at step K = {K} of {len(st)}; table 2^{R} complex64 = {2 ** R * 8 / 1e9:.1f} GB")
