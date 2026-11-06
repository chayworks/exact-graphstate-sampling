"""Check of the gate-cut bound for the #151 construction with qubit count n and CZ depth d chosen independently.

For a 1D brickwork, the bond between qubits a and a+1 carries the CZs of every layer with layer % 2 == a % 2, i.e.
ceil(d/2) or floor(d/2) of them; single-qubit gates do not change entanglement, and the cut rank of the output graph
state across {0..a} | {a+1..n-1} equals its entanglement in bits, which is at most the number of CZs across the cut.
So along the spatial order the width is at most ceil(d/2), whatever n is; any order is also capped by floor(n/2).
This script builds instances (frontier_instance's write_qasm generalised to depth d), takes the graph form, and
prints the maximum cut rank along the spatial order next to both bounds, plus the count of CZs across the worst bond.
Usage: python gatecut_check.py n1,d1 n2,d2 ..."""
import sys, random, numpy as np
from tableau import clifford_body
from graphform import graph_form
from order import profile
TH = 0.9553166181245092

def write_qasm(n, d, seed, path):
    rng = random.Random(seed); L = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{n}];"]
    for q in range(n): L.append(f"h q[{q}];")
    for layer in range(d):
        for a in range(layer % 2, n - 1, 2): L.append(f"cz q[{a}],q[{a + 1}];")
        for q in range(n):
            if rng.random() < 0.5: L.append(f"s q[{q}];")
            L.append(f"sx q[{q}];")
    for q in range(n): L += [f"tdg q[{q}];", f"h q[{q}];", f"rz({TH}) q[{q}];", f"h q[{q}];"]
    open(path, "w").write("\n".join(L) + "\n")

if __name__ == "__main__":
    import tempfile, os
    for arg in sys.argv[1:]:
        n, d = (int(x) for x in arg.split(","))
        path = os.path.join(tempfile.gettempdir(), f"gatecut_{n}_{d}.qasm"); write_qasm(n, d, 1000 + n + 7 * d, path)
        N, T, tail = clifford_body(path); A = graph_form(T)[0]
        A = np.array(A, dtype=np.uint8) % 2
        prof = profile(A, list(range(n))); w = max(prof.values()) if isinstance(prof, dict) else max(prof)
        print(f"n = {n:3d}, d = {d:3d}: spatial-order width {w:3d}; bound ceil(d/2) = {(d + 1) // 2:3d}; floor(n/2) = {n // 2:3d}", flush=True)
