"""Greedy slice-set search with the twist constraint, for the hybrid sampler on an undoped instance: starting from the
empty set, add the vertex (among the first K steps) that lowers the width most while every slice vertex keeps a valid
twist (its suffix neighbourhood in the span of the sliced cut basis at K). Prints width after 0..maxk slices.
Usage: python twist_greedy.py <circuit.qasm> <order.json> <K> <maxk> [out.json]"""
import sys, json, numpy as np
from tableau import clifford_body, gf2_rank
from amp import Amplitudes
from inplace import plan_inplace
from sliced_plans import gf2_solve
qasm, oj, K, maxk = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]); order = json.load(open(oj))
N, T, tail = clifford_body(qasm); am = Amplitudes(N, T, tail, order); A = am.A.copy(); n = len(A)
def prof(A2): return [gf2_rank(A2[np.ix_(order[k:], order[:k])]) for k in range(1, n)]
def sliced(S):
    A2 = A.copy()
    for v in S: A2[v, :] = 0; A2[:, v] = 0
    return A2
def valid(S, A2, w):
    st = plan_inplace(A2, order, max(1, w)); B = st[K]["B"]; suf = set(order[K:])
    return all(gf2_solve(B, sum(1 << int(u) for u in np.nonzero(A[v])[0] if u in suf)) is not None for v in S)
S = []; p = prof(A); hist = [(0, max(p), [])]; print(f"n = {n}, K = {K}, width {max(p)}, widest steps {[k + 1 for k, r in enumerate(p) if r == max(p)]}")
for k in range(1, maxk + 1):
    best = None
    for v in order[:K]:
        if v in S: continue
        S2 = S + [v]; A2 = sliced(S2); p2 = prof(A2); w2 = max(p2); score = (w2, sum(1 for r in p2 if r == w2))   # width, then how many steps sit at it
        if best is None or score < best[0]: best = (score, v, A2, w2)
    (w2, nw), v, A2, _ = best
    if not valid(S + [v], A2, w2):
        # fall back to the best valid candidate
        cands = []
        for u in order[:K]:
            if u in S or u == v: continue
            A3 = sliced(S + [u]); p3 = prof(A3); w3 = max(p3)
            if valid(S + [u], A3, w3): cands.append(((w3, sum(1 for r in p3 if r == w3)), u))
        if not cands: print("no valid candidate"); break
        (w2, nw), v = min(cands)
    S.append(v); hist.append((k, w2, list(S))); print(f"{k} slices {S}: width {w2} ({nw} steps at that width)", flush=True)
if len(sys.argv) > 5: json.dump({"n": n, "K": K, "history": hist}, open(sys.argv[5], "w"))
