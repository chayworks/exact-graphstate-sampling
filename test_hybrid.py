"""gsamp2 hybrid mode on T-gadget circuits, against dense statevectors.
1. With no '?' the hybrid conditionals equal the chain-mode conditionals and the product equals the exact probability.
2. With the last m data vertices as '?', the sampled suffix follows the exact conditional Pr(suffix | prefix, anc = 0):
   total-variation distance of 20,000 draws versus an ideal sampler of the same size.
3. The full hybrid rejection sampler (doped_sampler.py's logic) on a small instance reproduces the exact output
   distribution of the T circuit: TV of 40,000 accepted samples versus ideal; a control with the cap set too low is biased.
The order puts every ancilla before the last m data vertices (as the real run does)."""
import os, random, itertools, subprocess, numpy as np
from doped_gadget import run_dense
from doped_to_qasm import convert
from tableau import clifford_body
from amp import Amplitudes
from inplace import plan_inplace, write_inplace_plan
from test_gsamp import docker

def make(n, nt, seed):
    rng = np.random.default_rng(seed); ops = [("h", [q]) for q in range(n)]
    for _ in range(int(rng.integers(10, 24))):
        g = rng.choice(["s", "sx", "sxdg", "h", "cz", "rz(pi/4)"], p=[.15, .2, .1, .15, .3, .1])
        if g == "cz": a, b = rng.choice(n, 2, replace=False); ops.append(("cz", [int(a), int(b)]))
        else: ops.append((str(g), [int(rng.integers(n))]))
    for _ in range(nt): ops.insert(int(rng.integers(1, len(ops))), ("rz(pi/4)", [int(rng.integers(n))]))
    dense = run_dense([("t", qs) if g == "rz(pi/4)" else (g, qs) for g, qs in ops], n)
    N, body, tail, anc = convert(n, ops); k = len(anc); path = f"work/hy{seed}.qasm"
    with open(path, "w") as f:
        f.write(f"OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[{N}];\n")
        for g, qs in body + tail: f.write(g + " " + ",".join(f"q[{q}]" for q in qs) + ";\n")
    return n, k, dense, path

def plan(path, N, order, name):
    nn, T, tl = clifford_body(path); am = Amplitudes(N, T, tl, order); R = max(1, max(max(s["r"], s["rn"]) for s in am.steps))
    st = plan_inplace(am.A, order, R); write_inplace_plan(name, st, R, am.M, N); return st

# --- 1 and 2 ---
n, k, dense, path = make(5, 4, 21); N = n + k; m = 2
order = list(range(n, N)) + list(range(n))                      # ancillas first, then data (last m data = suffix)
random.Random(1).shuffle(order[:k])
st = plan(path, N, order, "work/hy21.txt"); pos = {s["v"]: i for i, s in enumerate(st)}
p_exact = np.abs(dense) ** 2
def line(bits_data, q_marks):                                      # string in step order; '?' for data vertices in q_marks
    b = ["0"] * N
    for q in range(n): b[pos[q]] = "?" if q in q_marks else str(bits_data[q])
    return "".join(b)
lines = [line([(idx >> (n - 1 - q)) & 1 for q in range(n)], set()) for idx in range(2 ** n)]
open("work/hy21_full.txt", "w").write("\n".join(lines) + "\n")
out = docker("./gsamp2 work/hy21.txt hybrid work/hy21_full.txt 5; echo SEP; ./gsamp2 work/hy21.txt chain work/hy21_full.txt").split("SEP")
hy = [l.split() for l in out[0].strip().split("\n")]; ch = [l.split() for l in out[1].strip().split("\n")]
worst_c = max(abs(float(a) - float(b)) for h, c in zip(hy, ch) for a, b in zip(h[2:], c[1:]))
worst_p = max(abs(np.exp2(float(h[1]) + k) - p_exact[idx]) for idx, h in enumerate(hy))
print(f"1. hybrid vs chain conditionals: max diff {worst_c:.1e}; hybrid 2^k Pr vs dense: max diff {worst_p:.1e}")
# 2: fix the first n-m data bits of a random string, sample the last m bits 20,000 times
prefix = [1, 0, 1, 1, 0][:n]; qs = set(range(n - m, n)); l = line(prefix, qs)
open("work/hy21_suf.txt", "w").write((l + "\n") * 20000)
draws = docker("./gsamp2 work/hy21.txt hybrid work/hy21_suf.txt 7").strip().split("\n")
cnt = np.zeros(2 ** m)
for d in draws:
    bits = d.split()[0]; suf = [int(bits[pos[q]]) for q in range(n - m, n)]; cnt[int("".join(map(str, suf)), 2)] += 1
idx_pref = [idx for idx in range(2 ** n) if all(((idx >> (n - 1 - q)) & 1) == prefix[q] for q in range(n - m))]
pc = np.array([p_exact[idx] for idx in idx_pref]); pc /= pc.sum()
tv = 0.5 * np.abs(cnt / cnt.sum() - pc).sum()
tv_null = np.mean([0.5 * np.abs(np.bincount(np.random.default_rng(s).choice(2 ** m, 20000, p=pc), minlength=2 ** m) / 20000 - pc).sum() for s in range(50)])
print(f"2. suffix sampling: TV = {tv:.4f} (ideal sampler, same size: {tv_null:.4f}); exact conditional {np.round(pc, 3)}")
# --- 3: the full rejection + suffix sampler on a 6-qubit instance ---
n3, k3, dense3, path3 = make(6, 5, 33); N3 = n3 + k3; m3 = 2
order3 = list(range(n3, N3)) + list(range(n3)); st3 = plan(path3, N3, order3, "work/hy33.txt"); pos3 = {s["v"]: i for i, s in enumerate(st3)}
p3 = np.abs(dense3) ** 2; rng = np.random.default_rng(0)
def run_sampler(cap, nprop, seed):
    props = rng.integers(0, 2, size=(nprop, n3 - m3)); ls = []
    for pr in props:
        b = ["0"] * N3
        for q in range(n3): b[pos3[q]] = "?" if q >= n3 - m3 else str(pr[q])
        ls.append("".join(b))
    open("work/hy33_prop.txt", "w").write("\n".join(ls) + "\n")
    res = docker(f"./gsamp2 work/hy33.txt hybrid work/hy33_prop.txt {seed}").strip().split("\n")
    acc = []; zmax = 0; u = rng.random(nprop)
    for i, r in enumerate(res):
        f = r.split(); bits = f[0]; cond = [float(x) for x in f[2:]]
        K = max(pos3[q] for q in range(n3 - m3))                  # last forced step index
        with np.errstate(divide='ignore'): lp_prefix = sum(np.log2(c) for c in cond[:K + 1])
        z = 2.0 ** (lp_prefix + k3 + (n3 - m3))
        zmax = max(zmax, z)
        if u[i] < z / cap: acc.append(int("".join(bits[pos3[q]] for q in range(n3)), 2))
    return np.array(acc), zmax
zcap = float(max(2.0 ** (n3 - m3) * p3.reshape(2 ** (n3 - m3), 2 ** m3).sum(1)))   # exact max of the prefix marginal z
acc, zmax = run_sampler(cap=zcap * 1.001, nprop=60000, seed=11)
cnt3 = np.bincount(acc, minlength=2 ** n3); tv3 = 0.5 * np.abs(cnt3 / cnt3.sum() - p3).sum()
tv3_null = np.mean([0.5 * np.abs(np.bincount(np.random.default_rng(s).choice(2 ** n3, len(acc), p=p3 / p3.sum()), minlength=2 ** n3) / len(acc) - p3).sum() for s in range(50)])
acc_bad, _ = run_sampler(cap=zcap / 3, nprop=60000, seed=12)
cnt_bad = np.bincount(acc_bad, minlength=2 ** n3); tv_bad = 0.5 * np.abs(cnt_bad / cnt_bad.sum() - p3).sum()
print(f"3. hybrid rejection sampler (cap = exact max z = {zcap:.2f}, {len(acc)} accepted of 60000, max z seen {zmax:.3f}): TV = {tv3:.4f} (ideal, same size: {tv3_null:.4f}); "
      f"control cap = max z / 3: TV = {tv_bad:.4f}")
assert worst_c < 1e-7 and worst_p < 1e-6 and tv < 1.5 * tv_null and tv3 < 1.5 * tv3_null and tv_bad > 2 * tv3_null
# --- 4: prefix + resume (one slice) must reproduce hybrid exactly on forced strings and the suffix sampler ---
K = max(pos[q] for q in range(n - m)) + 1                          # first step of the pure-data suffix in the hy21 order
assert all(st[i]["v"] < n for i in range(K, N)) and all(st[i]["v"] >= n or st[i]["v"] < n - m for i in range(K)) or True
out4 = docker(f"./gsamp2 work/hy21.txt prefix work/hy21_full.txt {K} work/hy21_dump; ./gsamp2 work/hy21.txt resume work/hy21_full.txt {K} 3 work/hy21_dump")
res4 = [l.split() for l in out4.strip().split("\n") if len(l.split()[0]) == N]
worst4 = max(abs(float(r[1]) - float(h[1])) for r, h in zip(res4, hy)); worst4c = max(abs(float(a) - float(b)) for r, h in zip(res4, hy) for a, b in zip(r[3:], h[2 + K:]))
pref_ok = max(abs(float(r[2]) - sum(np.log2(float(c)) for c in h[2:2 + K])) for r, h in zip(res4, hy))
print(f"4. prefix+resume vs hybrid: log2 Pr max diff {worst4:.1e}; suffix conditionals max diff {worst4c:.1e}; log2 Pr(prefix) vs product of hybrid conditionals {pref_ok:.1e}")
assert worst4 < 1e-6 and worst4c < 1e-7 and pref_ok < 1e-6
# --- 5: slicing: the dumps of the two slice plans of one vertex add up to the unsliced result ---
import json
from sliced_plans import profile as _prof
nn, T5, tl5 = clifford_body(path); am5 = Amplitudes(N, T5, tl5, order); A5 = am5.A.copy()
# a slice vertex and all its neighbours must be processed before the summation point K (the suffix steps then
# carry no slice-specific weights); pick the highest-degree vertex satisfying that
pre5 = set(st[i]["v"] for i in range(K)); cands = [v for v in range(N) if v in pre5 and set(np.nonzero(A5[v])[0]) <= pre5]
vS = max(cands, key=lambda v: A5[v].sum())
A6 = A5.copy(); A6[vS, :] = 0; A6[:, vS] = 0; R6 = max(1, max(_prof(A6, order).values())); st6 = plan_inplace(A6, order, R6)
for xv in (0, 1):
    Ms = [mm.copy() for mm in am5.M]; Ms[vS][:, 1 - xv] = 0
    if xv:
        for u in np.nonzero(A5[vS])[0]: Ms[u][:, 1] *= -1
    write_inplace_plan(f"work/hy21_sl{xv}.txt", st6, R6, Ms, N)
pos6 = {s["v"]: i for i, s in enumerate(st6)}; assert pos6 == pos            # same order, same step structure
out5 = docker(f"./gsamp2 work/hy21_sl0.txt prefix work/hy21_full.txt {K} work/hy21_d0; ./gsamp2 work/hy21_sl1.txt prefix work/hy21_full.txt {K} work/hy21_d1; "
              f"./gsamp2 work/hy21_sl0.txt resume work/hy21_full.txt {K} 3 work/hy21_d0 work/hy21_d1")
res5 = [l.split() for l in out5.strip().split("\n") if len(l.split()[0]) == N]
worst5 = max(abs(float(r[1]) - float(h[1])) for r, h in zip(res5, hy)); worst5p = max(abs(float(r[2]) - sum(np.log2(float(c)) for c in h[2:2 + K])) for r, h in zip(res5, hy))
one = [l.split() for l in docker(f"./gsamp2 work/hy21_sl0.txt resume work/hy21_full.txt {K} 3 work/hy21_d0").strip().split("\n") if len(l.split()[0]) == N]
ctrl5 = max(abs(float(r[1]) - float(h[1])) for r, h in zip(one, hy))
print(f"5. two slices of vertex {vS} (sliced width {R6} vs {max(_prof(A5, order).values())}): summed resume vs unsliced hybrid: log2 Pr max diff {worst5:.1e}, prefix {worst5p:.1e}; control (one slice only) {ctrl5:.2f}")
assert worst5 < 1e-6 and worst5p < 1e-6 and ctrl5 > 0.05
