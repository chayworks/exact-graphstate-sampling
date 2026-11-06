"""The paper's own graph metrics (src/analysis_scripts/analyze_nisq_results.py of IQP-NISQ, reimplemented without
networkx or pickles): edge probability, percentage of bipartite graphs, and the TVD between the degree distribution
and Binomial(N_nodes - 1, p_target), for IBM's hardware shots and for exact classical samples of the same models.
A string of n(n-1)/2 bits is the strict upper triangle of the adjacency matrix in np.triu_indices order.
Usage: python graph_metrics.py <IQP-NISQ results dir> <datasets_summary.yml> [K exact samples per model]"""
import sys, math, numpy as np, yaml
from pathlib import Path
from collections import defaultdict
from heavyhex import gate_list
from iqp_exact import IQPModel

def bipartite(adj):
    n = len(adj); col = [-1] * n
    for s in range(n):
        if col[s] >= 0: continue
        col[s] = 0; st = [s]
        while st:
            u = st.pop()
            for v in np.nonzero(adj[u])[0]:
                if col[v] < 0: col[v] = 1 - col[u]; st.append(v)
                elif col[v] == col[u]: return False
    return True

def metrics(S, nodes, p):
    iu = np.triu_indices(nodes, 1); E = S.mean(); bip = 0; deg = np.zeros(nodes, np.int64)
    for s in S:
        A = np.zeros((nodes, nodes), np.uint8); A[iu] = s; A[(iu[1], iu[0])] = s
        bip += bipartite(A); deg += np.bincount(A.sum(1), minlength=nodes)
    emp = deg / deg.sum(); k = np.arange(nodes)
    binom = np.array([math.comb(nodes - 1, int(j)) * p ** j * (1 - p) ** (nodes - 1 - j) for j in k])
    return E, 100 * bip / len(S), 0.5 * np.abs(emp - binom).sum()

B = Path(sys.argv[1]); summ = yaml.safe_load(open(sys.argv[2])); K = int(sys.argv[3]) if len(sys.argv) > 3 else 512
target = {(e["nodes"], e["graph_type"], e["density_category"]): e["average_density"] for e in summ}
rng = np.random.default_rng(5); out = defaultdict(list)
print(f"{'model':42s} {'target p':>8s} | {'hardware: p, bip%, TVD':>24s} | {'exact classical: p, bip%, TVD':>30s}")
for f in sorted(B.glob("evaluation_results_iqp/*/*/*/rank_*/samples_*.npy")):
    rel = f.parent.relative_to(B / "evaluation_results_iqp"); pp = B / "trained_params_iqp" / rel / "params.npy"
    if not pp.exists(): continue
    ds, _, dens, _ = rel.parts; nodes = int(ds.split("N_")[0]); gtype = ds.split("N_")[1]
    p = target.get((nodes, gtype, dens))
    if p is None: print("no target for", rel); continue
    th = np.load(pp, allow_pickle=False); H = np.load(f, allow_pickle=False); N = H.shape[1]
    g = gate_list(N); L = len(th) // len(g); m = IQPModel(N, g * L, th)
    C = np.array([y for y, _ in m.sample(rng, K)])
    mh, mc = metrics(H, nodes, p), metrics(C, nodes, p)
    out[(N, gtype)].append((p, mh, mc))
    print(f"{rel.as_posix():42s} {p:8.3f} | {mh[0]:6.3f} {mh[1]:6.1f} {mh[2]:6.3f}         | {mc[0]:6.3f} {mc[1]:6.1f} {mc[2]:6.3f}", flush=True)
print("\nmeans over models | qubits | type | target p | hardware p, bipartite %, degree TVD | exact classical p, bipartite %, degree TVD")
for (N, gt) in sorted(out):
    a = out[(N, gt)]; P = np.mean([x[0] for x in a]); h = np.mean([x[1] for x in a], 0); c = np.mean([x[2] for x in a], 0)
    print(f"| {N} | {gt} | {P:.3f} | {h[0]:.3f}, {h[1]:.1f}, {h[2]:.3f} | {c[0]:.3f}, {c[1]:.1f}, {c[2]:.3f} |")
