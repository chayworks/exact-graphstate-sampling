"""The real in-place sampler binary (gsamp2) on T-gadget circuits: prob mode must equal the dense statevector of
the original T circuit, |amp_T|^2 = 2^k * Pr_gadget(b, ancillas = 0), on every data string of random small circuits
with mid-circuit T gates. Negative control: plans with the M rows swapped must disagree."""
import os, random, itertools, numpy as np
from doped_gadget import run_dense
from doped_to_qasm import convert
from tableau import clifford_body
from amp import Amplitudes
from inplace import plan_inplace, write_inplace_plan
from test_gsamp import docker
rng = np.random.default_rng(11); cmds = []; cases = []
for i in range(12):
    n = int(rng.integers(2, 6)); ops = [("h", [q]) for q in range(n)]
    for _ in range(int(rng.integers(8, 20))):
        g = rng.choice(["s", "sx", "sxdg", "h", "cz", "rz(pi/4)"], p=[.15, .2, .1, .15, .3, .1])
        if g == "cz": a, b = rng.choice(n, 2, replace=False); ops.append(("cz", [int(a), int(b)]))
        else: ops.append((str(g), [int(rng.integers(n))]))
    for _ in range(int(rng.integers(1, 5))): ops.insert(int(rng.integers(1, len(ops))), ("rz(pi/4)", [int(rng.integers(n))]))
    dense = run_dense([("t", qs) if g == "rz(pi/4)" else (g, qs) for g, qs in ops], n)
    N, body, tail, anc = convert(n, ops); k = len(anc)
    path = f"work/gd{i}.qasm"
    with open(path, "w") as f:
        f.write(f"OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[{N}];\n")
        for g, qs in body + tail: f.write(g + " " + ",".join(f"q[{q}]" for q in qs) + ";\n")
    nn, T, tl = clifford_body(path)
    order = list(range(N)); random.Random(i).shuffle(order)
    am = Amplitudes(N, T, tl, order); R = max(1, max(max(s["r"], s["rn"]) for s in am.steps))
    st = plan_inplace(am.A, order, R); write_inplace_plan(f"work/gd{i}.txt", st, R, am.M, N)
    if i < 4: write_inplace_plan(f"work/gdbad{i}.txt", st, R, [m[::-1].copy() for m in am.M], N)
    with open(f"work/gd{i}_str.txt", "w") as f:
        for bits in itertools.product([0, 1], repeat=n):
            b = list(bits) + [0] * k; f.write("".join(str(b[s["v"]]) for s in st) + "\n")
    cases.append((n, k, dense)); cmds.append(f"echo CASE{i}; ./gsamp2 work/gd{i}.txt prob work/gd{i}_str.txt")
    if i < 4: cmds.append(f"echo BAD{i}; ./gsamp2 work/gdbad{i}.txt prob work/gd{i}_str.txt")
blocks = {}; key = None
for line in docker("; ".join(cmds)).split("\n"):
    if line.startswith(("CASE", "BAD")): key = line; blocks[key] = []
    elif line.strip(): blocks[key].append(float(line))
worst = 0.0; ctrl = 0.0
for i, (n, k, dense) in enumerate(cases):
    p = np.exp2(np.array(blocks[f"CASE{i}"]) + k); assert len(p) == 2 ** n
    worst = max(worst, np.abs(p - np.abs(dense) ** 2).max())
    if i < 4: ctrl = max(ctrl, np.abs(np.exp2(np.array(blocks[f"BAD{i}"]) + k) - np.abs(dense) ** 2).max())
print(f"gsamp2 prob mode on {len(cases)} T-gadget circuits (2-5 data qubits, 1-8 T gates), every data string: max |2^k Pr_gadget - |amp_T|^2| = {worst:.2e}; negative control {ctrl:.3f}")
assert worst < 1e-6 and ctrl > 1e-2
