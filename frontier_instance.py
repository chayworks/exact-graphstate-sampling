"""Build an n-qubit instance of the #151 construction (frontier.ibm_like, seed 1000+n, the same seed as the
frontier table), with the #151 measurement tail on every qubit, search a vertex order, and write:
  work/fi<n>.qasm            the circuit (so every other tool, e.g. the independent engine, can read it)
  work/fi<n>_order.json      the order found (max cut rank = width)
  work/fi<n>_inplace_plan.txt  in-place sampler plan (gsamp2), R = width
  work/fi<n>_step_qubits.json  qubit of each step
Prints width, sum_k 2^{r_k} (work per pass) and the memory of the in-place table.
Usage: python frontier_instance.py <n> <ngreedy> <iters> [seed]"""
import sys, json, math, random
from tableau import Tableau, clifford_body
from graphform import graph_form
from order import profile
from order2 import greedy, key
from amp import Amplitudes
from inplace import plan_inplace, write_inplace_plan
TH = 0.9553166181245092

def write_qasm(n, seed, path):
    rng = random.Random(seed); L = [f"OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{n}];"]
    for q in range(n): L.append(f"h q[{q}];")
    for layer in range(n):
        for a in range(layer % 2, n - 1, 2): L.append(f"cz q[{a}],q[{a + 1}];")
        for q in range(n):
            if rng.random() < 0.5: L.append(f"s q[{q}];")
            L.append(f"sx q[{q}];")
    for q in range(n): L += [f"tdg q[{q}];", f"h q[{q}];", f"rz({TH}) q[{q}];", f"h q[{q}];"]
    open(path, "w").write("\n".join(L) + "\n")

def search(A, seed, ngreedy, iters):
    n = len(A); nb = [sum(int(A[i][j]) << j for j in range(n)) for i in range(n)]; rng = random.Random(seed)
    cands = [(key(profile(A, o)), o) for o in (greedy(nb, n, rng, 0.9) for _ in range(ngreedy))]
    ck, cur = min(cands, key=lambda t: t[0])
    for it in range(iters):
        a, b = sorted(rng.sample(range(n), 2)); b = min(b, a + 12); new = cur[:]
        x = new.pop(a); new.insert(b, x); nk = key(profile(A, new))
        if nk <= ck: cur, ck = new, nk
    return cur, ck

if __name__ == "__main__":
    n, ng, it = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); seed = int(sys.argv[4]) if len(sys.argv) > 4 else n; tag = f"{n}" if len(sys.argv) <= 4 else f"{n}_s{seed}"
    write_qasm(n, 1000 + n, f"work/fi{n}.qasm")
    nq, T, tail = clifford_body(f"work/fi{n}.qasm"); assert nq == n
    A, _ = graph_form(T)
    order, k = search(A, seed, ng, it)
    am = Amplitudes(n, T, tail, order); R = max(max(s["r"], s["rn"]) for s in am.steps)
    steps = plan_inplace(A, order, R); write_inplace_plan(f"work/fi{tag}_inplace_plan.txt", steps, R, am.M, n)
    json.dump(order, open(f"work/fi{tag}_order.json", "w")); json.dump([s["v"] for s in am.steps], open(f"work/fi{tag}_step_qubits.json", "w"))
    work = sum(2 ** s["rn"] for s in am.steps)
    print(json.dumps(dict(n=n, edges=int(A.sum()) // 2, width=R, natural=max(profile(A, list(range(n)))), log2_sum_2r=round(math.log2(work), 2),
                          inplace_GB=round(8 * 2 ** R / 1e9, 1), order_file=f"work/fi{tag}_order.json")))
