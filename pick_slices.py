"""Greedy choice of k sliced vertices for dp.cpp on any instance: repeatedly remove the vertex that most lowers the
DP cost of the order restricted to the remaining vertices. Usage: python pick_slices.py <qasm> <order.json> <k> <out.json>"""
import sys, json, numpy as np
from tableau import clifford_body
from graphform import graph_form
from order import profile, cost
F, OJ, k, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
n, T, _ = clifford_body(F); A, _ = graph_form(T); order = json.load(open(OJ))
def prof(S):
    keep = [q for q in order if q not in S]; return profile(A[np.ix_(keep, keep)], list(range(len(keep))))
S = []; print(json.dumps(dict(S=[], width=max(prof(S)), log2_cost=round(float(np.log2(cost(prof(S)))), 2))))
for _ in range(k):
    best = min((cost(prof(S + [v])), v) for v in order if v not in S); S.append(best[1])
    pr = prof(S); print(json.dumps(dict(S=S, width=max(pr), log2_cost_per_slice=round(float(np.log2(cost(pr))), 2), log2_total=round(float(np.log2(cost(pr))) + len(S), 2))))
json.dump(S, open(out, "w"))
