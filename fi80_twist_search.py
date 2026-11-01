"""Slice sets for the 80-qubit hybrid run: K = 60 (20 bits sampled). A slice vertex may have suffix neighbours when the
twist of sliced_plans.py exists (its suffix neighbourhood in the span of the sliced cut basis at K). Prints the widths of
all valid single vertices and of the best pairs. Usage: python fi80_twist_search.py <qasm> <order.json> <K> [pairs]"""
import sys, json, itertools, numpy as np
from tableau import clifford_body, gf2_rank
from amp import Amplitudes
from inplace import plan_inplace
from sliced_plans import gf2_solve, profile
qasm, oj, K = sys.argv[1], sys.argv[2], int(sys.argv[3]); order = json.load(open(oj))
N, T, tail = clifford_body(qasm); am = Amplitudes(N, T, tail, order); A = am.A.copy(); n = len(A)
def width(A2): return max(gf2_rank(A2[np.ix_(order[k:], order[:k])]) for k in range(1, n))
def valid(S):
    A2 = A.copy()
    for v in S: A2[v, :] = 0; A2[:, v] = 0
    w = width(A2); st = plan_inplace(A2, order, max(1, w)); B = st[K]["B"]; suf = set(order[K:])
    ok = all(gf2_solve(B, sum(1 << int(u) for u in np.nonzero(A[v])[0] if u in suf)) is not None for v in S)
    return w, ok
w0 = width(A); print(f"unsliced width {w0}; K = {K}")
single = sorted((valid([v]), v) for v in order[:K]); print("single vertices (width, twist valid), vertex:", [(w, ok, v) for (w, ok), v in single[:12]])
json.dump({"single": [(int(w), bool(ok), int(v)) for (w, ok), v in single]}, open("work/fi80_twist_single.json", "w"))
if "pairs" in sys.argv:
    cand = [v for (w, ok), v in single if ok][:40]; best = []
    for a, b in itertools.combinations(cand, 2):
        (w, ok) = valid([a, b])
        if ok: best.append((w, a, b))
    best.sort(); print("best valid pairs (width, a, b):", best[:15])
    json.dump({"pairs": best[:60]}, open("work/fi80_twist_pairs.json", "w"))
