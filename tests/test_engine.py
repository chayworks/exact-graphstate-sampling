"""dp.cpp (real binary, in the container) with vertex slicing == brute force.
Kills: wrong slice prefactor sign (edges inside S), missing neighbour Z-phase from sliced vertices,
pdep/pivmask errors, kernel-span errors in C++. Negative control: dropping the neighbour phase must disagree."""
import random, numpy as np
from ref import brute
import engine
from engine import Sliced, write_plan, write_jobs, run_dp

rng = random.Random(5); cases = []
for t in range(40):
    n = rng.randint(3, 12); p = rng.random()
    A = np.zeros((n, n), np.uint8)
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < p: A[i, j] = A[j, i] = 1
    w = [(complex(rng.gauss(0, 1), rng.gauss(0, 1)), complex(rng.gauss(0, 1), rng.gauss(0, 1))) for _ in range(n)]
    if t % 2: w = [(np.complex128(a), np.complex128(b)) for a, b in w]   # numpy scalars, as Amplitudes.weights returns
    order = list(range(n)); rng.shuffle(order); S =rng.sample(range(n), rng.randint(0, min(3, n - 1)))
    cases.append((A, w, order, S))
worst = 0
for i, (A, w, order, S) in enumerate(cases):
    sl = Sliced(A, order, S); write_plan(f"work/p{i}.txt", sl.steps, sl.m); write_jobs(f"work/j{i}.txt", sl, [w])
z = []
import subprocess, os
script = " && ".join(f"./dp work/p{i}.txt work/j{i}.txt" for i in range(len(cases)))
out = subprocess.run(["docker", "run", "--rm", "-v", os.path.abspath(".").replace("\\", "/") + ":/w", "-w", "/w", "-e", "OMP_NUM_THREADS=4",
                      "gcc:14", "bash", "-c", script], capture_output=True, text=True, env={**os.environ, "MSYS_NO_PATHCONV": "1"})
assert out.returncode == 0, out.stderr
vals = [complex(*map(float, l.split())) for l in out.stdout.strip().splitlines()]
assert len(vals) == len(cases)
for (A, w, order, S), zd in zip(cases, vals):
    zb = brute(A, w); worst = max(worst, abs(zb - zd) / max(1, abs(zb)))
# negative control: slices without the neighbour phase
A, w, order, S = next(c for c in cases if c[3] and any(c[0][c[3][0], q] for q in range(len(c[0])) if q not in c[3]))
sl = Sliced(A, order, S); bad = 0
from ref import run
for pre, ww in sl.slices(w):
    raw = [w[q] for q in sl.keep]; bad += pre * run(sl.steps, raw)
print(f"40 cases (C++ binary, slicing 0-3 vertices): max rel err {worst:.2e}; negative control rel diff {abs(bad - brute(A, w)) / max(1, abs(brute(A, w))):.3f}")
assert worst < 1e-10 and abs(bad - brute(A, w)) > 1e-6
