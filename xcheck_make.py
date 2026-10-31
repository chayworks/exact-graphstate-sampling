"""Independent-engine job maker for any instance: writes a sliced dp.cpp plan and a jobs file for the given sampled
strings, so dp.cpp (a different algorithm from the samplers) can recompute their exact probabilities.
Usage: python xcheck_make.py <qasm> <order.json> <sliced_vertices.json> <samples.txt> <nstrings> <out_prefix>
Samples file: "<bits in step order> <log2Pr> <seed>" per line (gsamp/gsamp2 sample output).
Writes <out_prefix>_plan.txt and <out_prefix>_jobs.txt and prints the sampler's log2 Pr of each chosen string."""
import sys, json
from tableau import clifford_body
from amp import Amplitudes
from engine import Sliced, write_plan, write_jobs

F, OJ, SJ, SF, k, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6]
n, T, tail = clifford_body(F); order = json.load(open(OJ)); S = json.load(open(SJ))
am = Amplitudes(n, T, tail, order); steps_v = [st["v"] for st in am.steps]
sl = Sliced(am.A, order, S); write_plan(f"{out}_plan.txt", sl.steps, sl.m)
rows = [l.split() for l in open(SF) if len(l.split()) >= 2 and len(l.split()[0]) == n][:k]
ws = []
for r in rows:
    b = [0] * n
    for pos, v in enumerate(steps_v): b[v] = int(r[0][pos])
    ws.append(am.weights(b))
write_jobs(f"{out}_jobs.txt", sl, ws)
print(json.dumps(dict(n=n, sliced=S, dp_width=max(max(s["r"], s["rn"]) for s in sl.steps), strings=[(r[0], float(r[1])) for r in rows])))
