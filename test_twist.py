"""Slicing with suffix neighbours (the 'twist' of sliced_plans.py / gsamp2 resume --twist).
1. Toy T-gadget circuit (the hy21 case of test_hybrid.py, 5 data qubits + 4 ancillas, suffix = last 2 data vertices):
   slice a vertex that has a neighbour in the suffix; the twisted sum of the two slice dumps must reproduce the
   unsliced hybrid result on all 32 forced strings (log2 Pr, prefix marginal) and the dense statevector; the untwisted
   sum (control) must not. A two-vertex slice checks the XOR of the twist vectors.
2. Production scale: the 70-qubit fstate instance with the width-30 order, K = 50 (20 bits sampled): the hybrid run on
   3 proposals equals the twisted two-slice run bit for bit (same seed) and in log2 Pr; control untwisted differs.
Run test_hybrid.py's prerequisites are not needed; this builds its own files (work/tw21*)."""
import json, random, numpy as np
from doped_gadget import run_dense
from doped_to_qasm import convert
from tableau import clifford_body
from amp import Amplitudes
from inplace import plan_inplace, write_inplace_plan
from sliced_plans import make_sliced_plans, gf2_solve, profile
import os, subprocess
def docker(cmd):
    out = subprocess.run(["docker", "run", "--rm", "-v", os.path.abspath(".").replace("\\", "/") + ":/w", "-w", "/w",
                          "-e", "OMP_NUM_THREADS=8", "gcc:14", "bash", "-c", cmd], capture_output=True, text=True,
                         env={**os.environ, "MSYS_NO_PATHCONV": "1"})
    assert out.returncode == 0, out.stderr; return out.stdout

def build_toy(n, nt, seed):
    rng = np.random.default_rng(seed); ops = [("h", [q]) for q in range(n)]
    for _ in range(int(rng.integers(10, 24))):
        g = rng.choice(["s", "sx", "sxdg", "h", "cz", "rz(pi/4)"], p=[.15, .2, .1, .15, .3, .1])
        if g == "cz": a, b = rng.choice(n, 2, replace=False); ops.append(("cz", [int(a), int(b)]))
        else: ops.append((str(g), [int(rng.integers(n))]))
    for _ in range(nt): ops.insert(int(rng.integers(1, len(ops))), ("rz(pi/4)", [int(rng.integers(n))]))
    dense = run_dense([("t", qs) if g == "rz(pi/4)" else (g, qs) for g, qs in ops], n)
    N, body, tail, anc = convert(n, ops); path = f"work/tw{seed}.qasm"
    with open(path, "w") as f:
        f.write(f"OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[{N}];\n")
        for g, qs in body + tail: f.write(g + " " + ",".join(f"q[{q}]" for q in qs) + ";\n")
    return n, len(anc), dense, path

def run_sliced(base, nsl, K, strings, seed, twist):
    cmd = "; ".join(f"./gsamp2 {base}_s{format(j, f'0{nsl.bit_length() - 1}b') if nsl > 1 else ''}.txt prefix {strings} {K} {base}_d{j}" for j in range(nsl))
    dumps = " ".join(f"{base}_d{j}" for j in range(nsl)); tw = f" --twist {base}_twist.txt" if twist else ""
    plan0 = f"{base}_s{'0' * (nsl.bit_length() - 1)}.txt"
    out = docker(f"{cmd}; ./gsamp2 {plan0} resume {strings} {K} {seed} {dumps}{tw}")
    ns = int(open(plan0).readline().split()[1])                      # resume lines start with the full step string
    return [l.split() for l in out.strip().split("\n") if l.split() and len(l.split()[0]) == ns]

# ---------------- 1. toy ----------------
def twist_ok(A, order, K, S):
    """True if every slice vertex's suffix neighbourhood lies in the span of the sliced plan's cut basis at step K."""
    A2 = A.copy()
    for v in S: A2[v, :] = 0; A2[:, v] = 0
    st2 = plan_inplace(A2, order, max(1, max(profile(A2, order).values()))); B = st2[K]["B"]; sufs = set(order[K:])
    return all(gf2_solve(B, sum(1 << int(u) for u in np.nonzero(A[v])[0] if u in sufs)) is not None for v in S)
# a toy whose suffix (last m data vertices) has a prefix vertex with suffix neighbours and a valid twist
for seed in range(21, 60):
    n, k, dense, path = build_toy(6, 4, seed); N = n + k; m = 3
    order = list(range(n, N)) + list(range(n)); random.Random(1).shuffle(order[:k])
    nn, T, tl = clifford_body(path); am = Amplitudes(N, T, tl, order); A = am.A.copy(); M = [x.copy() for x in am.M]
    R = max(1, max(max(s["r"], s["rn"]) for s in am.steps)); st = plan_inplace(A, order, R)
    pos = {s["v"]: i for i, s in enumerate(st)}; K = max(pos[q] for q in range(n - m)) + 1; suf = set(order[K:]); pre = set(order[:K])
    with_suf = [v for v in pre if any(u in suf for u in np.nonzero(A[v])[0]) and twist_ok(A, order, K, [v])]
    if with_suf: break
print(f"toy seed {seed}: {n} data + {k} ancillas, K = {K}, suffix {sorted(suf)}, slice candidates with suffix neighbours and a valid twist: {with_suf}")
write_inplace_plan("work/tw21.txt", st, R, M, N)
lines = []
for idx in range(2 ** n):
    b = ["0"] * N
    for q in range(n): b[pos[q]] = str((idx >> (n - 1 - q)) & 1)
    lines.append("".join(b))
open("work/tw21_full.txt", "w").write("\n".join(lines) + "\n")
hy = [l.split() for l in docker("./gsamp2 work/tw21.txt hybrid work/tw21_full.txt 5").strip().split("\n")]
p_exact = np.abs(dense) ** 2
wd = max(abs(np.exp2(float(h[1]) + k) - p_exact[i]) for i, h in enumerate(hy)); print(f"toy: hybrid 2^k Pr vs dense max diff {wd:.1e}; k = {k}, N = {N}, K = {K}"); assert wd < 1e-7
assert with_suf, "no toy with a prefix vertex with a suffix neighbour and a valid twist"
vS = max(with_suf, key=lambda v: A[v].sum())
info = make_sliced_plans(A, M, order, [vS], "work/tw21a", K=K, ND=n); assert info["twist"] and info["step_vertices"] == [s["v"] for s in st]
res = run_sliced("work/tw21a", 2, K, "work/tw21_full.txt", 3, True)
d_lp = max(abs(float(r[1]) - float(h[1])) for r, h in zip(res, hy)); d_pre = max(abs(float(r[2]) - sum(np.log2(float(c)) for c in h[2:2 + K])) for r, h in zip(res, hy))
ctrl = run_sliced("work/tw21a", 2, K, "work/tw21_full.txt", 3, False); d_ctrl = max(abs(float(r[1]) - float(h[1])) for r, h in zip(ctrl, hy))
print(f"1a. one slice vertex {vS} (suffix neighbours {sorted(u for u in np.nonzero(A[vS])[0] if u in suf)}): twisted sum vs unsliced hybrid: log2 Pr {d_lp:.1e}, prefix {d_pre:.1e}; untwisted control {d_ctrl:.2f}")
assert d_lp < 1e-6 and d_pre < 1e-6 and d_ctrl > 0.05
# two slice vertices (both with suffix neighbours if possible)
others = [v for v in pre if v != vS and twist_ok(A, order, K, [vS, v])]
others.sort(key=lambda v: (-any(u in suf for u in np.nonzero(A[v])[0]), -A[v].sum()))      # prefer a second vertex with suffix neighbours
v2 = others[0]; info2 = make_sliced_plans(A, M, order, [vS, v2], "work/tw21b", K=K, ND=n)
res2 = run_sliced("work/tw21b", 4, K, "work/tw21_full.txt", 3, info2["twist"])
d2 = max(abs(float(r[1]) - float(h[1])) for r, h in zip(res2, hy)); d2p = max(abs(float(r[2]) - sum(np.log2(float(c)) for c in h[2:2 + K])) for r, h in zip(res2, hy))
print(f"1b. two slice vertices {vS},{v2}: four twisted slices vs unsliced: log2 Pr {d2:.1e}, prefix {d2p:.1e} (twist {info2['twist']})")
assert d2 < 1e-6 and d2p < 1e-6

# ---------------- 2. 70-qubit fstate, width 30, K = 50 ----------------
order70 = json.load(open("cloud2/w30_step_qubits.json")); n70 = 70; K70 = 50
nn, T70, tl70 = clifford_body("cloud7/nq70.qasm"); am70 = Amplitudes(n70, T70, tl70, order70); A70 = am70.A.copy(); M70 = [x.copy() for x in am70.M]
w70 = max(profile(A70, order70).values()); st70 = plan_inplace(A70, order70, w70)
pre70 = set(order70[:K70]); suf70 = set(order70[K70:]); B70 = st70[K70]["B"]
best = None
for v in order70[:K70]:
    A2 = A70.copy(); A2[v, :] = 0; A2[:, v] = 0; wv = max(profile(A2, order70).values())
    s_v = sum(1 << int(u) for u in np.nonzero(A70[v])[0] if u in suf70)
    if gf2_solve(B70, s_v) is not None and (best is None or wv < best[0]): best = (wv, v)
print(f"2.  fstate order width {w70}; best single slice with a valid twist: vertex {best[1]} -> width {best[0]}")
rng = random.Random(11); props = []
for _ in range(3):
    b = ["?"] * n70
    for i in range(K70): b[i] = str(rng.randrange(2))
    props.append("".join(b))
open("work/tw70_props.txt", "w").write("\n".join(props) + "\n")
hy70 = [l.split() for l in docker("./gsamp2 cloud2/nq70_inplace_plan.txt hybrid work/tw70_props.txt 17").strip().split("\n")]
info70 = make_sliced_plans(A70, M70, order70, [best[1]], "work/tw70", K=K70, ND=n70); assert info70["R"] == best[0]
res70 = run_sliced("work/tw70", 2, K70, "work/tw70_props.txt", 17, info70["twist"])
same_bits = all(r[0] == h[0] for r, h in zip(res70, hy70)); d70 = max(abs(float(r[1]) - float(h[1])) for r, h in zip(res70, hy70))
d70p = max(abs(float(r[2]) - sum(np.log2(float(c)) for c in h[2:2 + K70])) for r, h in zip(res70, hy70))
ctrl70 = run_sliced("work/tw70", 2, K70, "work/tw70_props.txt", 17, False); c70 = max(abs(float(r[3]) - float(h[2 + K70])) for r, h in zip(ctrl70, hy70)); t70 = max(abs(float(r[3]) - float(h[2 + K70])) for r, h in zip(res70, hy70))
print(f"2.  70 qubits, slice vertex {best[1]} (width {info70['R']}, twist {info70['twist']}): same sampled bits {same_bits}; log2 Pr {d70:.1e}, prefix {d70p:.1e}; first suffix conditional: twisted diff {t70:.1e}, untwisted control diff {c70:.3f} (the prefix norm alone does not see the shift on a stabilizer state)")
print("    log2 Pr of the 3 proposals:", [h[1] for h in hy70])
assert same_bits and d70 < 1e-5 and d70p < 1e-5 and t70 < 1e-6 and c70 > 1e-4   # stabilizer conditionals move little; the toy's control is the strong one
