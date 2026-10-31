"""gsamp2.cpp (in-place arena, real binary): prob mode == statevector on every string; sample-mode TV vs
an ideal sampler. Kills: wrong base tables, wrong g slot, rhoc/par sign split, tau mixups.
Negative control: M rows swapped in the plan must disagree."""
import os, random, subprocess, numpy as np
from amp import Amplitudes
from inplace import plan_inplace, write_inplace_plan
from test_gsamp import build, docker

def plan_for(n, T, tail, order, bad=False):
    am = Amplitudes(n, T, tail, order); R = max(1, max(max(s["r"], s["rn"]) for s in am.steps))
    st = plan_inplace(am.A, order, R); M = [m[::-1].copy() for m in am.M] if bad else am.M
    return am, st, R, M

rng = random.Random(9); cmds = []; cases = []
for i in range(25):
    n = rng.randint(2, 9); T, tail, psi = build(n, rng, rng.randint(10, 80))
    order = list(range(n)); rng.shuffle(order)
    am, st, R, M = plan_for(n, T, tail, order); write_inplace_plan(f"work/gi{i}.txt", st, R, M, n)
    if i < 5: _, stb, Rb, Mb = plan_for(n, T, tail, order, bad=True); write_inplace_plan(f"work/gibad{i}.txt", stb, Rb, Mb, n)
    with open(f"work/gi{i}_str.txt", "w") as f:
        for idx in range(2 ** n):
            b = [(idx >> (n - 1 - q)) & 1 for q in range(n)]; f.write("".join(str(b[s["v"]]) for s in st) + "\n")
    cases.append((n, psi)); cmds.append(f"echo CASE{i}; ./gsamp2 work/gi{i}.txt prob work/gi{i}_str.txt")
    if i < 5: cmds.append(f"echo BAD{i}; ./gsamp2 work/gibad{i}.txt prob work/gi{i}_str.txt")
blocks = {}; key = None
for line in docker("; ".join(cmds)).split("\n"):
    if line.startswith(("CASE", "BAD")): key = line; blocks[key] = []
    elif line.strip(): blocks[key].append(float(line))
worst = 0; ctrl = 0
for i, (n, psi) in enumerate(cases):
    p = np.exp2(np.array(blocks[f"CASE{i}"])); assert len(p) == 2 ** n
    worst = max(worst, np.abs(p - np.abs(psi) ** 2).max())
    if i < 5: ctrl = max(ctrl, np.abs(np.exp2(np.array(blocks[f"BAD{i}"])) - np.abs(psi) ** 2).max())
n = 9; T, tail, psi = build(n, random.Random(77), 120); order = list(range(n)); random.Random(1).shuffle(order)
am, st, R, M = plan_for(n, T, tail, order); write_inplace_plan("work/gi_s.txt", st, R, M, n)
cnt = np.zeros(2 ** n)
for l in docker("./gsamp2 work/gi_s.txt sample 20000 1000").strip().split("\n"):
    bits = l.split()[0]; b = [0] * n
    for k, s in enumerate(st): b[s["v"]] = int(bits[k])
    cnt[int("".join(map(str, b)), 2)] += 1
pe = np.abs(psi) ** 2; tv = 0.5 * np.abs(cnt / cnt.sum() - pe).sum()
tv_null = np.mean([0.5 * np.abs(np.bincount(np.random.default_rng(s).choice(2 ** n, 20000, p=pe / pe.sum()), minlength=2 ** n) / 20000 - pe).sum() for s in range(20)])
print(f"prob mode: 25 circuits, every string, max |p - p_sv| = {worst:.2e}; negative control = {ctrl:.3f}")
print(f"sample mode: 20000 samples, TV = {tv:.4f} (ideal sampler, same size: {tv_null:.4f})")
assert worst < 1e-6 and ctrl > 1e-2 and tv < 1.5 * tv_null
