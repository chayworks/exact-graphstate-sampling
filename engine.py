"""Write plans / jobs for dp.cpp, with vertex slicing, and run it in the gcc:14 container."""
import os, subprocess, itertools, numpy as np
from plan import plan

def write_plan(path, steps, m):
    maxr = max(max(s["r"], s["rn"]) for s in steps)
    with open(path, "w") as f:
        f.write(f"{m} {len(steps)} {maxr}\n")
        for s in steps:
            pm = sum(1 << p for p in s["piv"])
            f.write(f"{s['v']} {s['r']} {s['rn']} {s['rho']} {pm} {len(s['ker'])} " + " ".join(map(str, s["ker"])) + "\n")
    return maxr

class Sliced:
    """Z(w) for graph A with the vertices in S summed outside the DP (graph A - S inside)."""
    def __init__(self, A, order, S):
        n = len(A); self.S = list(S); self.keep = [q for q in order if q not in set(S)]
        idx = {q: i for i, q in enumerate(self.keep)}; self.idx = idx; m = len(self.keep)
        Ak = np.zeros((m, m), np.uint8)
        for i, p in enumerate(self.keep):
            for j, q in enumerate(self.keep): Ak[i, j] = A[p, q]
        self.A = A; self.m = m; self.steps = plan(Ak, list(range(m)))
    def slices(self, w):
        """Yield (prefactor, weights for kept vertices) for every assignment of the sliced vertices."""
        A, S = self.A, self.S
        for xs in itertools.product((0, 1), repeat=len(S)):
            pre = 1 + 0j; ph = 0
            for a, (v, xv) in enumerate(zip(S, xs)):
                pre *= w[v][xv]
                for v2, x2 in list(zip(S, xs))[a + 1:]: ph ^= int(A[v, v2]) & xv & x2
            pre *= (-1) ** ph
            ww = []
            for q in self.keep:
                f = sum(int(A[q, v]) & xv for v, xv in zip(S, xs)) & 1
                ww.append((w[q][0], -w[q][1] if f else w[q][1]))
            yield pre, ww

def write_jobs(path, sl, ws):
    with open(path, "w") as f:
        f.write(f"{len(ws)}\n")
        for w in ws:
            sls = list(sl.slices(w)); f.write(f"{len(sls)}\n")
            for pre, ww in sls:
                f.write("%.17g %.17g\n" % (pre.real, pre.imag))
                for a, b in ww: f.write("%.17g %.17g %.17g %.17g\n" % (a.real, a.imag, b.real, b.imag))

def run_dp(workdir, plan_file, jobs_file, threads=16, mem="8g"):
    host = os.path.abspath(workdir).replace("\\", "/")
    cmd = ["docker", "run", "--rm", "--memory", mem, "-v", f"{host}:/w", "-w", "/w", "-e", f"OMP_NUM_THREADS={threads}",
           "gcc:14", "./dp", plan_file, jobs_file]
    out = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, "MSYS_NO_PATHCONV": "1"})
    if out.returncode: raise RuntimeError(out.stderr)
    return [complex(float(a), float(b)) for a, b in (l.split() for l in out.stdout.strip().splitlines())]

def write_sampler_plan(path, am):
    """Plan + per-step measurement matrices M_v (rows b, columns x) for gsamp.cpp."""
    maxr = max(max(s["r"], s["rn"]) for s in am.steps)
    with open(path, "w") as f:
        f.write(f"{am.n} {len(am.steps)} {maxr}\n")
        for s in am.steps:
            pm = sum(1 << p for p in s["piv"]); M = am.M[s["v"]]
            f.write(f"{s['v']} {s['r']} {s['rn']} {s['rho']} {pm} {len(s['ker'])} " + " ".join(map(str, s["ker"])) + " "
                    + " ".join("%.17g %.17g" % (complex(M[b, x]).real, complex(M[b, x]).imag) for b in (0, 1) for x in (0, 1)) + "\n")
    return maxr
