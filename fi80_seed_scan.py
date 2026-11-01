"""Width of 80-qubit instances of the #151 construction for several circuit seeds (frontier_instance.py fixes the
circuit seed at 1000+n; this scans others to see whether any instance of the family has width 34, which fits a
256 GiB node unsliced). Each line: circuit seed, width found by search(), order file. Usage:
  python fi80_seed_scan.py <n> <ngreedy> <iters> <cseed1> [<cseed2> ...]"""
import sys, json, math
from frontier_instance import write_qasm, search
from tableau import clifford_body
from graphform import graph_form
from order import profile
n, ng, it = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
for cs in map(int, sys.argv[4:]):
    write_qasm(n, cs, f"work/fi{n}_c{cs}.qasm"); nq, T, tail = clifford_body(f"work/fi{n}_c{cs}.qasm"); A, _ = graph_form(T)
    order, k = search(A, 1, ng, it); w = max(profile(A, order)); json.dump(order, open(f"work/fi{n}_c{cs}_order.json", "w"))
    print(json.dumps(dict(n=n, cseed=cs, width=w, natural=max(profile(A, list(range(n)))), edges=int(A.sum()) // 2, order_file=f"work/fi{n}_c{cs}_order.json")), flush=True)
