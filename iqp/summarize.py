"""Per-size summary of results/score_hardware.txt: log-likelihood gap (bits per shot by which IBM's shots are less
likely under the trained model than exact samples of it), with the spread across models; median linear XEB."""
import re, sys, numpy as np
from collections import defaultdict
by = defaultdict(list)
for l in open(sys.argv[1]):
    m = re.search(r"N=\s*(\d+).*w=(\d+) log2p hw\s+(\S+) ideal\s+(\S+).*F_XEB\s+(\S+).*sampler (\S+) ms", l)
    if m: N, w, hw, ide, F, t = m.groups(); by[int(N)].append((int(w), float(hw), float(ide), float(F), float(t)))
print("| qubits | models | width | log-likelihood gap, bits/shot: mean [min, max] | median linear XEB | exact sampler ms/sample |")
print("|---|---|---|---|---|---|")
for N in sorted(by):
    a = np.array(by[N]); gap = a[:, 2] - a[:, 1]
    print(f"| {N} | {len(a)} | {int(a[:,0].max())} | {gap.mean():.2f} [{gap.min():.2f}, {gap.max():.2f}] | {np.median(a[:,3]):.2f} | {a[:,4].mean():.1f} |")
