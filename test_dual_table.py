"""The in-place table is the dual (Walsh) representation: T(u) is the summed weight of the prefix assignments whose
linear functional on the suffix is phi(u) = XOR_i u_i B_i (B = the cut basis stored by plan_inplace). This is the
mechanism behind the slicing twist (a phase (-1)^{<s, x_suffix>} translates the table by a with s = sum a_i B_i), so it
is checked directly: gsamp2 'prefix' dumps of a toy T-gadget circuit (unsliced and both slices of a vertex with a suffix
neighbour) must equal a brute-force dual table entry by entry (after undoing the dump's logscale), and the brute-force
slice identity F = F_0 + (-1)^{<s,x>} F_1 must hold on every suffix assignment. Uses the toy of test_twist.py (seed 23)."""
import struct, itertools, random, os, subprocess, numpy as np
from tableau import clifford_body
from amp import Amplitudes
from inplace import plan_inplace, write_inplace_plan
from sliced_plans import make_sliced_plans, gf2_solve, profile
from doped_gadget import run_dense
from doped_to_qasm import convert

def build_toy(n, nt, seed):
    rng = np.random.default_rng(seed); ops = [("h", [q]) for q in range(n)]
    for _ in range(int(rng.integers(10, 24))):
        g = rng.choice(["s", "sx", "sxdg", "h", "cz", "rz(pi/4)"], p=[.15, .2, .1, .15, .3, .1])
        if g == "cz": a, bq = rng.choice(n, 2, replace=False); ops.append(("cz", [int(a), int(bq)]))
        else: ops.append((str(g), [int(rng.integers(n))]))
    for _ in range(nt): ops.insert(int(rng.integers(1, len(ops))), ("rz(pi/4)", [int(rng.integers(n))]))
    N, body, tail, anc = convert(n, ops); path = f"work/dual{seed}.qasm"
    with open(path, "w") as f:
        f.write(f"OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[{N}];\n")
        for g, qs in body + tail: f.write(g + " " + ",".join(f"q[{q}]" for q in qs) + ";\n")
    return n, len(anc), path
def docker(cmd):
    out = subprocess.run(["docker", "run", "--rm", "-v", os.path.abspath(".").replace("\\", "/") + ":/w", "-w", "/w",
                          "-e", "OMP_NUM_THREADS=4", "gcc:14", "bash", "-c", cmd], capture_output=True, text=True, env={**os.environ, "MSYS_NO_PATHCONV": "1"})
    assert out.returncode == 0, out.stderr; return out.stdout

n, k, path = build_toy(6, 4, 23); N = n + k; m = 3
order = list(range(n, N)) + list(range(n)); random.Random(1).shuffle(order[:k])
nn, T, tl = clifford_body(path); am = Amplitudes(N, T, tl, order); A = am.A.copy(); M = [x.copy() for x in am.M]
R = max(1, max(max(s["r"], s["rn"]) for s in am.steps)); st = plan_inplace(A, order, R); write_inplace_plan("work/dual_full.txt", st, R, M, N)
pos = {s["v"]: i for i, s in enumerate(st)}; K = max(pos[q] for q in range(n - m)) + 1; P = order[:K]; S = order[K:]
# two forced strings: all zeros, and alternating bits
strings = ["0" * N, "".join("01"[i % 2] for i in range(N))]; open("work/dual_str.txt", "w").write("\n".join(strings) + "\n")

def read_dump(fn):
    d = open(fn, "rb").read(); nk, ntq = struct.unpack("ii", d[:8]); ls, = struct.unpack("d", d[8:16])
    return nk, ntq, ls, np.frombuffer(d[16:], dtype=np.complex64).reshape(1 << nk, ntq)
def brute(A2, Ms, b):
    Td = {}
    for xs in itertools.product([0, 1], repeat=K):
        x = dict(zip(P, xs)); w = 1.0 + 0j; phi = 0
        for v in P:
            w *= Ms[v][b[v], x[v]]
            if x[v]: phi ^= sum(1 << int(u) for u in np.nonzero(A2[v])[0] if u in S)
        for i, v in enumerate(P):
            for u in P[i + 1:]:
                if A2[v, u] and x[v] and x[u]: w = -w
        Td[phi] = Td.get(phi, 0) + w
    return Td
def compare(plan, A2, Ms, tag, stK):
    docker(f"./gsamp2 {plan} prefix work/dual_str.txt {K} work/dual_{tag}"); worst = 0; nonzero = 0
    B = stK["B"]; keep = stK["keep"]; Tq = stK["Tq"]; inv = {v: q for q, v in Tq.items()}; tqs = sorted(Tq.values())
    for li, s in enumerate(strings):
        b = [int(ch) for ch in s]; b = {st[i]["v"]: b[i] for i in range(N)}         # outcome per vertex (string is in step order)
        nk, ntq, ls, vals = read_dump(f"work/dual_{tag}_{li}.bin"); Td = brute(A2, Ms, b)
        for c in range(1 << nk):
            for ti, tq in enumerate(tqs):
                u = sum(((c >> j) & 1) << keep[j] for j in range(len(keep))) ^ inv[tq]; phi = 0
                for i in range(len(B)):
                    if (u >> i) & 1: phi ^= int(B[i])
                dv = complex(vals[c, ti]) * 2 ** ls; bv = Td.get(phi, 0); worst = max(worst, abs(dv - bv)); nonzero += abs(bv) > 1e-9
    return worst, nonzero
w_full, nz = compare("work/dual_full.txt", A, M, "full", st[K])
print(f"unsliced: dump vs brute-force dual table, max |diff| {w_full:.1e} over 2 strings ({nz} nonzero entries)")
vS = next(v for v in sorted(P, key=lambda v: -int(A[v].sum())) if any(u in S for u in np.nonzero(A[v])[0])
          and gf2_solve(plan_inplace(np.where(np.arange(N)[:, None] == v, 0, np.where(np.arange(N)[None, :] == v, 0, A)), order, R)[K]["B"],
                        sum(1 << int(u) for u in np.nonzero(A[v])[0] if u in S)) is not None)
info = make_sliced_plans(A, M, order, [vS], "work/dual_sl", K=K, ND=n); A2 = A.copy(); A2[vS, :] = 0; A2[:, vS] = 0; st2 = plan_inplace(A2, order, info["R"])
s_v = sum(1 << int(u) for u in np.nonzero(A[vS])[0] if u in S); Msl = []
for xv in (0, 1):
    Ms = [x.copy() for x in M]; Ms[vS][:, 1 - xv] = 0
    if xv:
        for u in np.nonzero(A[vS])[0]:
            if u in set(P): Ms[u][:, 1] *= -1
    Msl.append(Ms); w, nz = compare(f"work/dual_sl_s{xv}.txt", A2, Ms, f"s{xv}", st2[K]); print(f"slice x_{vS} = {xv}: dump vs brute-force dual table, max |diff| {w:.1e} ({nz} nonzero entries)"); w_full = max(w_full, w)
assert w_full < 1e-6
def F(Td, xs): return sum(w * (-1) ** bin(phi & xs).count("1") for phi, w in Td.items())
worst = 0; ctrl = 0
for s in strings:
    b = {st[i]["v"]: int(s[i]) for i in range(N)}; Tf = brute(A, M, b); T0 = brute(A2, Msl[0], b); T1 = brute(A2, Msl[1], b)
    for bits in itertools.product([0, 1], repeat=len(S)):
        xs = sum(bit << int(v) for bit, v in zip(bits, S)); sg = (-1) ** bin(s_v & xs).count("1")
        worst = max(worst, abs(F(Tf, xs) - (F(T0, xs) + sg * F(T1, xs)))); ctrl = max(ctrl, abs(F(Tf, xs) - (F(T0, xs) + F(T1, xs))))
print(f"slice identity F = F0 + (-1)^<s,x> F1 on all suffix assignments: max |diff| {worst:.1e}; without the sign (control) {ctrl:.3f}")
assert worst < 1e-9 and ctrl > 1e-3
