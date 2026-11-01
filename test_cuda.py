"""The CUDA port (gsamp_cuda.cu) against dense statevectors and against the CPU binary gsamp2.
 1. small T-gadget circuits: prob mode on every data string vs dense (as test_gsamp2_doped.py);
 2. hybrid, prefix+resume (one dump) and sliced prefix+resume (two dumps) vs the CPU binary on the hy21 case;
 3. a generated 50-qubit instance (fstate construction, width about 22): 24 random strings, GPU prob and chain
    conditionals vs CPU gsamp2, plus 200 samples' reported log2 Pr re-scored, plus timing of one pass.
Builds the GPU binary inside the CUDA image (nvcc -arch=native). Negative control: M rows swapped must disagree."""
import os, random, itertools, subprocess, json, time, numpy as np
from doped_gadget import run_dense
from doped_to_qasm import convert
from tableau import clifford_body
from amp import Amplitudes
from inplace import plan_inplace, write_inplace_plan
from test_gsamp import docker as docker_cpu
IMG = "nvidia/cuda:12.8.1-devel-ubuntu24.04"
def docker_gpu(cmd):
    out = subprocess.run(["docker", "run", "--rm", "--gpus", "all", "-v", os.path.abspath(".").replace("\\", "/") + ":/w", "-w", "/w", "--entrypoint", "bash", IMG, "-c", cmd],
                         capture_output=True, text=True, env={**os.environ, "MSYS_NO_PATHCONV": "1"})
    assert out.returncode == 0, out.stderr[-2000:]; return out.stdout
print(docker_gpu("nvcc -O3 -arch=native -o gsamp_cuda gsamp_cuda.cu 2>&1 | tail -3; echo built"))
# --- 1 ---
rng = np.random.default_rng(5); cmds = []; cases = []
for i in range(10):
    n = int(rng.integers(2, 6)); ops = [("h", [q]) for q in range(n)]
    for _ in range(int(rng.integers(8, 20))):
        g = rng.choice(["s", "sx", "sxdg", "h", "cz", "rz(pi/4)"], p=[.15, .2, .1, .15, .3, .1])
        if g == "cz": a, b = rng.choice(n, 2, replace=False); ops.append(("cz", [int(a), int(b)]))
        else: ops.append((str(g), [int(rng.integers(n))]))
    for _ in range(int(rng.integers(1, 5))): ops.insert(int(rng.integers(1, len(ops))), ("rz(pi/4)", [int(rng.integers(n))]))
    dense = run_dense([("t", qs) if g == "rz(pi/4)" else (g, qs) for g, qs in ops], n)
    N, body, tail, anc = convert(n, ops); k = len(anc); path = f"work/cu{i}.qasm"
    with open(path, "w") as f:
        f.write(f"OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[{N}];\n")
        for g, qs in body + tail: f.write(g + " " + ",".join(f"q[{q}]" for q in qs) + ";\n")
    nn, T, tl = clifford_body(path); order = list(range(N)); random.Random(i).shuffle(order)
    am = Amplitudes(N, T, tl, order); R = max(1, max(max(s["r"], s["rn"]) for s in am.steps)); st = plan_inplace(am.A, order, R)
    write_inplace_plan(f"work/cu{i}.txt", st, R, am.M, N)
    if i < 3: write_inplace_plan(f"work/cubad{i}.txt", st, R, [m[::-1].copy() for m in am.M], N)
    with open(f"work/cu{i}_str.txt", "w") as f:
        for bits in itertools.product([0, 1], repeat=n):
            b = list(bits) + [0] * k; f.write("".join(str(b[s["v"]]) for s in st) + "\n")
    cases.append((n, k, dense)); cmds.append(f"echo CASE{i}; ./gsamp_cuda work/cu{i}.txt prob work/cu{i}_str.txt")
    if i < 3: cmds.append(f"echo BAD{i}; ./gsamp_cuda work/cubad{i}.txt prob work/cu{i}_str.txt")
blocks = {}; key = None
for line in docker_gpu("; ".join(cmds)).split("\n"):
    if line.startswith(("CASE", "BAD")): key = line; blocks[key] = []
    elif line.strip(): blocks[key].append(float(line))
worst = 0.0; ctrl = 0.0
for i, (n, k, dense) in enumerate(cases):
    p = np.exp2(np.array(blocks[f"CASE{i}"]) + k); worst = max(worst, np.abs(p - np.abs(dense) ** 2).max())
    if i < 3: ctrl = max(ctrl, np.abs(np.exp2(np.array(blocks[f"BAD{i}"]) + k) - np.abs(dense) ** 2).max())
print(f"1. GPU prob mode on 10 T-gadget circuits, every data string: max |2^k Pr - |amp_T|^2| = {worst:.2e}; negative control {ctrl:.3f}")
assert worst < 1e-6 and ctrl > 1e-2
# --- 2: hybrid / prefix+resume / sliced, GPU vs CPU on the hy21 case (files from test_hybrid.py) ---
assert os.path.exists("work/hy21.txt") and os.path.exists("work/hy21_sl0.txt"), "run test_hybrid.py first"
N21 = len(open("work/hy21_full.txt").readline().strip()); K = int(open("work/hy21_K.txt").read()) if os.path.exists("work/hy21_K.txt") else None
if K is None:
    st21 = json.load(open("work/hy21_steps.json")); K = max(i for i, v in enumerate(st21) if v < 3) + 1           # as in test_hybrid: suffix = data vertices 3, 4
cpu = docker_cpu(f"./gsamp2 work/hy21.txt hybrid work/hy21_full.txt 5").strip().split("\n")
gpu = docker_gpu(f"./gsamp_cuda work/hy21.txt hybrid work/hy21_full.txt 5").strip().split("\n")
d_h = max(abs(float(a.split()[1]) - float(b.split()[1])) for a, b in zip(cpu, gpu))
d_c = max(abs(float(x) - float(y)) for a, b in zip(cpu, gpu) for x, y in zip(a.split()[2:], b.split()[2:]))
gpu2 = docker_gpu(f"./gsamp_cuda work/hy21_sl0.txt prefix work/hy21_full.txt {K} work/cu_d0; ./gsamp_cuda work/hy21_sl1.txt prefix work/hy21_full.txt {K} work/cu_d1; ./gsamp_cuda work/hy21_sl0.txt resume work/hy21_full.txt {K} 3 work/cu_d0 work/cu_d1")
res = [l.split() for l in gpu2.strip().split("\n") if len(l.split()[0]) == N21]
d_s = max(abs(float(r[1]) - float(h.split()[1])) for r, h in zip(res, cpu))
print(f"2. hybrid GPU vs CPU: log2 Pr max diff {d_h:.1e}, conditionals {d_c:.1e}; sliced prefix+resume on GPU vs CPU hybrid: {d_s:.1e}")
assert d_h < 1e-6 and d_c < 1e-6 and d_s < 1e-6
# --- 3: 50-qubit fstate-construction instance ---
from frontier_instance import write_qasm
import frontier_instance
if not os.path.exists("work/cu50.qasm"): write_qasm(50, 1050, "work/cu50.qasm")
nn, T, tl = clifford_body("work/cu50.qasm"); order = list(range(nn)); am = Amplitudes(nn, T, tl, order); R = max(max(s["r"], s["rn"]) for s in am.steps)
st = plan_inplace(am.A, order, R); write_inplace_plan("work/cu50.txt", st, R, am.M, nn); json.dump([s["v"] for s in st], open("work/cu50_steps.json", "w"))
rs = random.Random(3); open("work/cu50_str.txt", "w").write("\n".join("".join(rs.choice("01") for _ in range(nn)) for _ in range(24)) + "\n")
t0 = time.time(); c_prob = docker_cpu("OMP_NUM_THREADS=16 ./gsamp2 work/cu50.txt chain work/cu50_str.txt").strip().split("\n"); tc = time.time() - t0
t0 = time.time(); g_prob = docker_gpu("./gsamp_cuda work/cu50.txt chain work/cu50_str.txt").strip().split("\n"); tg = time.time() - t0
d50 = max(abs(float(x) - float(y)) for a, b in zip(c_prob, g_prob) for x, y in zip(a.split(), b.split()))
samp = docker_gpu("./gsamp_cuda work/cu50.txt sample 200 777").strip().split("\n")
open("work/cu50_samp_str.txt", "w").write("\n".join(l.split()[0] for l in samp) + "\n")
resc = docker_cpu("OMP_NUM_THREADS=16 ./gsamp2 work/cu50.txt prob work/cu50_samp_str.txt").strip().split("\n")
d_samp = max(abs(float(a.split()[1]) - float(b)) for a, b in zip(samp, resc))
z = np.exp2(np.array([float(l.split()[1]) for l in samp]) + nn)
print(f"3. 50-qubit instance (width {R}): chain mode GPU vs CPU on 24 strings: max diff {d50:.1e} (CPU {tc:.1f} s, GPU {tg:.1f} s incl. container start); "
      f"200 GPU samples re-scored on CPU: max diff {d_samp:.1e}; samples' XEB {z.mean() - 1:+.3f} +/- {z.std(ddof=1) / 200 ** 0.5:.3f}")
assert d50 < 1e-5 and d_samp < 1e-5
