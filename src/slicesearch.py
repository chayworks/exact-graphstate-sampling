"""Greedy choice of sliced vertices: repeatedly remove the vertex that most lowers the DP cost of
(order minus S). Prints, per |S|, max width, log2 per-slice cost and log2 total cost (2^|S| slices)."""
import sys, json, math, numpy as np
from order import profile, cost
A = np.load("A_nq70.npy"); order = json.load(open(sys.argv[1])); smax = int(sys.argv[2])
def sub(S):
    keep = [q for q in order if q not in S]
    Ak = A[np.ix_(keep, keep)]; pr = profile(Ak, list(range(len(keep)))) if len(keep) > 1 else [0]
    return pr
S = []
for s in range(smax + 1):
    pr = sub(S); c = cost(pr)
    print(f"|S|={s:2d} max width {max(pr):2d}  log2 cost/slice {math.log2(c):6.2f}  log2 total {math.log2(c) + s:6.2f}  S={S}", flush=True)
    if s == smax: break
    best = None
    for v in order:
        if v in S: continue
        cc = cost(sub(S + [v]))
        if best is None or cc < best[0]: best = (cc, v)
    S.append(best[1])
json.dump(S, open("slice_S.json", "w"))
