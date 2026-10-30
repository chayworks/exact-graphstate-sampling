"""gsamp.cpp (real binary, container): 'prob' mode must give log2 |<b|psi>|^2 for every string, against a
dense statevector; 'sample' mode histogram must match (chi-square-style check on one larger case).
Kills: wrong M row/column order in the plan file, sign of rho term, float storage overflow,
off-by-one in bit order. Negative control: swapping M rows (b=0 <-> b=1) must disagree."""
import os, random, subprocess, numpy as np
from tableau import Tableau
from amp import Amplitudes, tail_matrix
from engine import write_sampler_plan
from test_tableau import G1, apply1, cz
TH = 0.9553166181245092

def build(n, rng, nops):
    ops = []
    for _ in range(nops):
        g = rng.choice(["h", "s", "sx", "x", "z", "cz", "cz", "cz"])
        ops.append((g, rng.sample(range(n), 2)) if g == "cz" else (g, [rng.randrange(n)]))
    tail = []
    for q in range(n): tail += [("tdg", [q]), ("h", [q]), (f"rz({TH})", [q]), ("h", [q])]
    T = Tableau(n); psi = np.zeros(2 ** n, complex); psi[0] = 1
    for name, qs in ops:
        T.apply(name, qs); psi = cz(psi, *qs, n) if name == "cz" else apply1(psi, G1[name], qs[0], n)
    for name, qs in tail: psi = apply1(psi, tail_matrix([(name, qs)]), qs[0], n)
    return T, tail, psi

def docker(cmd):
    out = subprocess.run(["docker", "run", "--rm", "-v", os.path.abspath(".").replace("\\", "/") + ":/w", "-w", "/w",
                          "-e", "OMP_NUM_THREADS=4", "gcc:14", "bash", "-c", cmd], capture_output=True, text=True,
                         env={**os.environ, "MSYS_NO_PATHCONV": "1"})
    assert out.returncode == 0, out.stderr; return out.stdout

rng = random.Random(9); cmds = []; cases = []
for i in range(25):
    n = rng.randint(2, 9); T, tail, psi = build(n, rng, rng.randint(10, 80))
    order = list(range(n)); rng.shuffle(order); am = Amplitudes(n, T, tail, order)
    write_sampler_plan(f"work/gs{i}.txt", am)
    if i < 5:  # negative control plan: rows of M swapped
        bad = Amplitudes(n, T, tail, order); bad.M = [m[::-1].copy() for m in bad.M]; write_sampler_plan(f"work/gsbad{i}.txt", bad)
    steps_v = [s["v"] for s in am.steps]
    with open(f"work/gs{i}_str.txt", "w") as f:
        for idx in range(2 ** n):
            b = [(idx >> (n - 1 - q)) & 1 for q in range(n)]; f.write("".join(str(b[v]) for v in steps_v) + "\n")
    cases.append((n, psi)); cmds.append(f"echo CASE{i}; ./gsamp work/gs{i}.txt prob work/gs{i}_str.txt")
    if i < 5: cmds.append(f"echo BAD{i}; ./gsamp work/gsbad{i}.txt prob work/gs{i}_str.txt")
out = docker("; ".join(cmds)).split("\n"); blocks = {}; key = None
for line in out:
    if line.startswith(("CASE", "BAD")): key = line; blocks[key] = []
    elif line.strip(): blocks[key].append(float(line))
worst = 0; ctrl = 0
for i, (n, psi) in enumerate(cases):
    p = np.exp2(np.array(blocks[f"CASE{i}"])); assert len(p) == 2 ** n
    worst = max(worst, np.abs(p - np.abs(psi) ** 2).max())
    if i < 5: ctrl = max(ctrl, np.abs(np.exp2(np.array(blocks[f"BAD{i}"])) - np.abs(psi) ** 2).max())
# sample mode: 20000 samples on a 9-qubit case, compare with exact
n = 9; T, tail, psi = build(n, random.Random(77), 120); order = list(range(n)); random.Random(1).shuffle(order)
am = Amplitudes(n, T, tail, order); write_sampler_plan("work/gs_s.txt", am); steps_v = [s["v"] for s in am.steps]
lines = docker("./gsamp work/gs_s.txt sample 20000 1000").strip().split("\n"); cnt = np.zeros(2 ** n)
maxlp = 0
for l in lines:
    bits, lp, _ = l.split(); b = [0] * n
    for k, v in enumerate(steps_v): b[v] = int(bits[k])
    idx = int("".join(map(str, b)), 2); cnt[idx] += 1; maxlp = max(maxlp, abs(np.exp2(float(lp)) - abs(psi[idx]) ** 2))
emp = cnt / cnt.sum(); tv = 0.5 * np.abs(emp - np.abs(psi) ** 2).sum()
tv_null = np.mean([0.5 * np.abs(np.bincount(np.random.default_rng(s).choice(2 ** n, 20000, p=np.abs(psi) ** 2 / np.sum(np.abs(psi) ** 2)), minlength=2 ** n) / 20000 - np.abs(psi) ** 2).sum() for s in range(20)])
print(f"prob mode: 25 circuits, every string, max |p - p_sv| = {worst:.2e}; negative control (M rows swapped) = {ctrl:.3f}")
print(f"sample mode: 20000 samples, TV to exact = {tv:.4f} (exact-sampler finite-size TV, 20 draws: {tv_null:.4f}); "
      f"reported Pr of each sample vs exact: max diff {maxlp:.2e}")
assert worst < 1e-6 and ctrl > 1e-2 and tv < 1.5 * tv_null and maxlp < 1e-6
