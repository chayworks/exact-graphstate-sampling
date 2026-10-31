"""Minimum stabilizer weight of generated #151-style instances versus n (how the smallest informative marginal grows).
For each n: build the instance (frontier.ibm_like, seed 1000+n, as in frontier.py), take its graph form, write the
adjacency, and run minweight_bin (exhaustive over popcount(v) <= K); the run is a proof when its minimum is <= K+1.
Usage: python minweight_scan.py "30:14,40:14,50:13,60:12" <path to minweight_bin, run via docker>  -> results/minweight_scan.txt"""
import sys, subprocess, json
from frontier import ibm_like
from graphform import graph_form
spec = [(int(a), int(b)) for a, b in (x.split(":") for x in sys.argv[1].split(","))]
with open("results/minweight_scan.txt", "a") as out:
    for n, K in spec:
        A, _ = graph_form(ibm_like(n, 1000 + n)); path = f"work/scan_A{n}.txt"
        with open(path, "w") as f:
            for i in range(n): f.write("".join(str(int(x)) for x in A[i]) + "\n")
        cmd = ["docker", "run", "--rm", "-v", "D:/GIT/Innovation/graphsamp:/w", "-w", "/w", "gcc:14", "bash", "-c",
               f"OMP_NUM_THREADS=16 ./work/minweight_bin {path} {K}"]
        r = subprocess.run(cmd, capture_output=True, text=True, env={**__import__('os').environ, "MSYS_NO_PATHCONV": "1"})
        last = [l for l in r.stdout.strip().splitlines() if l.strip()][-2:]
        rec = {"n": n, "edges": int(A.sum()) // 2, "K": K, "result": last}
        print(json.dumps(rec), flush=True); out.write(json.dumps(rec) + "\n"); out.flush()
