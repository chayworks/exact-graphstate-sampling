"""Design map for graph-state sampling on IBM's heavy-hex chips: the #151 construction (H on every qubit, then CZ layers
with a random sqrt(X) or S sqrt(X) on every qubit between layers, the face-state measurement tail) with the CZ layers
taken from a proper edge colouring of the heavy-hex coupling map of a 156-qubit Heron r2 (ibm_aachen layout, the map
written in IQP-NISQ's src/utils/hardware.py), cycled for CZ depth d, on the first n qubits of the chip.
For each (n, d) prints: the width along the chip's row order (an upper bound by Proposition 2's argument applied to
row cuts), the best width the order search of frontier_instance finds, the ceiling floor(n/2), and the in-place memory
2^width complex64. Usage: python hh_design.py n1,d1 n2,d2 ... [--iters I]"""
import sys, random, numpy as np
sys.path.insert(0, "iqp")
from heavyhex import aachen_edges
from tableau import clifford_body
from graphform import graph_form
from order import profile
from order2 import greedy, key
TH = 0.9553166181245092

def colour_classes(n):
    E = sorted({tuple(sorted(e)) for e in aachen_edges() if e[0] < n and e[1] < n})
    col = {}; used = {q: set() for q in range(n)}
    for (a, b) in E:                                   # greedy edge colouring (heavy-hex has degree <= 3)
        c = 0
        while c in used[a] or c in used[b]: c += 1
        col[(a, b)] = c; used[a].add(c); used[b].add(c)
    k = max(col.values()) + 1
    return [[e for e in E if col[e] == c] for c in range(k)]

def write_qasm(n, d, seed, path):
    rng = random.Random(seed); cls = colour_classes(n)
    L = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{n}];"] + [f"h q[{q}];" for q in range(n)]
    for layer in range(d):
        for a, b in cls[layer % len(cls)]: L.append(f"cz q[{a}],q[{b}];")
        for q in range(n):
            if rng.random() < 0.5: L.append(f"s q[{q}];")
            L.append(f"sx q[{q}];")
    for q in range(n): L += [f"tdg q[{q}];", f"h q[{q}];", f"rz({TH}) q[{q}];", f"h q[{q}];"]
    open(path, "w").write("\n".join(L) + "\n"); return len(cls)

def search(A, start, seed, ngreedy, iters):
    n = len(A); nb = [sum(int(A[i][j]) << j for j in range(n)) for i in range(n)]; rng = random.Random(seed)
    cands = [(key(profile(A, start)), start)] + [(key(profile(A, o)), o) for o in (greedy(nb, n, rng, 0.9) for _ in range(ngreedy))]
    ck, cur = min(cands, key=lambda t: t[0])
    for it in range(iters):
        a, b = sorted(rng.sample(range(n), 2)); b = min(b, a + 12); new = cur[:]
        x = new.pop(a); new.insert(b, x); nk = key(profile(A, new))
        if nk <= ck: cur, ck = new, nk
    return cur, ck

if __name__ == "__main__":
    import tempfile, os
    iters = 2000
    if "--iters" in sys.argv: i = sys.argv.index("--iters"); iters = int(sys.argv[i + 1]); del sys.argv[i:i + 2]
    for arg in sys.argv[1:]:
        n, d = (int(x) for x in arg.split(","))
        path = os.path.join(tempfile.gettempdir(), f"hh_{n}_{d}.qasm"); ncol = write_qasm(n, d, 2000 + 7 * n + d, path)
        N, T, tail = clifford_body(path); A = np.array(graph_form(T)[0], dtype=np.uint8) % 2
        wrow = max(profile(A, list(range(n))))
        order, k = search(A, list(range(n)), n + d, 4, iters)
        w = k[0]; mem = (1 << w) * 8 / 1e9
        print(f"n = {n:3d}, d = {d:3d} ({ncol} CZ colours): row-order width {wrow:3d}, searched width {w:3d}, ceiling {n // 2:3d}, "
              f"in-place table {mem:,.1f} GB", flush=True)
