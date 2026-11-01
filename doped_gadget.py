"""The doped instance as a graph state with projected ancillas.

Every mid-circuit T = rz(pi/4) on qubit q is replaced by a fresh ancilla a in |0>, CNOT(q -> a), and the projection
of a onto <phi| with |phi> = (|0> + e^{-i pi/4} |1>)/sqrt2:  <phi|_a CNOT(q->a) |psi>_q |0>_a = T|psi> / sqrt2.
The projection is a product measurement: apply T then H to a and take outcome 0.  The whole doped circuit is then a
Clifford circuit on 70 + 468 qubits followed by a product measurement (data in Z, ancillas as above), i.e. exactly
the setting of the cut-rank sampler/amplitude engine, with amplitude(doped) = 2^{468/2} * amplitude(gadget).

This script (1) checks the gadget identity on random small T-circuits against dense statevectors, (2) builds the
538-qubit tableau and its graph form, writes the adjacency and the per-vertex tails, and (3) reports cut-rank
widths for several vertex orders (natural, time-ordered ancillas first, spatial sweep by data qubit).
Usage: python doped_gadget.py <doped.qasm> <out_prefix>"""
import sys, re, json, itertools, numpy as np
from tableau import Tableau, gf2_rank
from graphform import graph_form

# ---------- dense helpers for the small-case test ----------
I2 = np.eye(2, dtype=complex); H = np.array([[1, 1], [1, -1]]) / np.sqrt(2); S = np.diag([1, 1j]); SDG = S.conj().T
SX = np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]]) / 2; SXDG = SX.conj().T; T = np.diag([1, np.exp(1j * np.pi / 4)])
ONE = {"h": H, "s": S, "sdg": SDG, "sx": SX, "sxdg": SXDG, "t": T, "x": np.array([[0, 1], [1, 0]], complex), "z": np.diag([1, -1]).astype(complex)}
def apply1(psi, G, q, n): psi = psi.reshape([2] * n); psi = np.moveaxis(np.tensordot(G, psi, axes=([1], [q])), 0, q); return psi.reshape(-1)
def apply_cz(psi, a, b, n):
    psi = psi.reshape([2] * n).copy(); idx = [slice(None)] * n; idx[a] = 1; idx[b] = 1; psi[tuple(idx)] *= -1; return psi.reshape(-1)
def apply_cnot(psi, c, t, n):
    psi = psi.reshape([2] * n).copy(); idx0 = [slice(None)] * n; idx0[c] = 1; sub = psi[tuple(idx0)]; sub = np.flip(sub, axis=t - (1 if t > c else 0)); psi[tuple(idx0)] = sub; return psi.reshape(-1)
def run_dense(ops, n):
    psi = np.zeros(2 ** n, complex); psi[0] = 1
    for g, qs in ops:
        if g == "cz": psi = apply_cz(psi, qs[0], qs[1], n)
        elif g == "cnot": psi = apply_cnot(psi, qs[0], qs[1], n)
        else: psi = apply1(psi, ONE[g], qs[0], n)
    return psi

def gadgetize(ops, n):
    """Replace each ('t', [q]) by an ancilla: returns (clifford ops on n + k qubits, k, ancilla list [(a, q, position)])."""
    out = []; anc = []; k = 0
    for pos, (g, qs) in enumerate(ops):
        if g == "t":
            a = n + k; k += 1; anc.append((a, qs[0], pos)); out.append(("cnot", [qs[0], a]))
        else: out.append((g, qs))
    return out, k, anc

def test_gadget(trials=20, seed=0):
    rng = np.random.default_rng(seed)
    for tr in range(trials):
        n = int(rng.integers(2, 5)); nt = int(rng.integers(1, 4)); ops = [("h", [q]) for q in range(n)]
        for _ in range(int(rng.integers(6, 14))):
            g = rng.choice(["s", "sx", "sxdg", "h", "cz", "t"], p=[.15, .2, .1, .15, .3, .1])
            if g == "cz" and n >= 2: a, b = rng.choice(n, 2, replace=False); ops.append(("cz", [int(a), int(b)]))
            elif g != "cz": ops.append((str(g), [int(rng.integers(n))]))
        for _ in range(nt): ops.insert(int(rng.integers(1, len(ops))), ("t", [int(rng.integers(n))]))
        psi = run_dense(ops, n)                                                   # the T circuit
        gops, k, anc = gadgetize(ops, n); phi = run_dense(gops, n + k)            # the gadget circuit, ancillas unmeasured
        for a, _, _ in anc: phi = apply1(apply1(phi, T, a, n + k), H, a, n + k)   # tails on ancillas: T then H, outcome 0
        phi = phi.reshape([2] * (n + k))
        proj = phi[(slice(None),) * n + (0,) * k].reshape(-1)                      # ancilla outcomes all 0
        assert np.allclose(proj * 2 ** (k / 2), psi, atol=1e-10), f"gadget identity fails on trial {tr}"
        # control: wrong tail (H only, no T) must not reproduce the T circuit
        phi2 = run_dense(gops, n + k)
        for a, _, _ in anc: phi2 = apply1(phi2, H, a, n + k)
        proj2 = phi2.reshape([2] * (n + k))[(slice(None),) * n + (0,) * k].reshape(-1)
        assert not np.allclose(proj2 * 2 ** (k / 2), psi, atol=1e-6), "control (missing T in the tail) unexpectedly agreed"
    print(f"gadget identity: {trials} random circuits with 1-3 T gates agree with dense statevectors to 1e-10; control fails as intended")

def parse(path):
    n = None; ops = []
    for line in open(path):
        line = line.strip().rstrip(";")
        if not line or line.startswith(("OPENQASM", "include")): continue
        m = re.match(r"qreg q\[(\d+)\]", line)
        if m: n = int(m.group(1)); continue
        name, args = line.split(" ", 1); qs = [int(a) for a in re.findall(r"q\[(\d+)\]", args)]
        if name == "rz(pi/4)": name = "t"
        ops.append((name, qs))
    return n, ops

def cut_profile(A, order):
    n = len(A); r = []
    for k in range(1, n):
        P = order[:k]; F = order[k:]; r.append(gf2_rank(A[np.ix_(F, P)]))
    return r

if __name__ == "__main__":
    test_gadget()
    n, ops = parse(sys.argv[1]); out = sys.argv[2]
    gops, k, anc = gadgetize(ops, n); N = n + k
    print(f"{sys.argv[1]}: {n} data qubits, {len(ops)} gates, {k} T gates -> {N} vertices")
    Tb = Tableau(N)
    for g, qs in gops:
        if g == "cnot": Tb.h(qs[1]); Tb.cz(qs[0], qs[1]); Tb.h(qs[1])
        elif g == "sxdg": Tb.h(qs[0]); Tb.sdg(qs[0]); Tb.h(qs[0])
        else: Tb.apply(g, qs)
    A, U = graph_form(Tb)
    with open(out + "_A.txt", "w") as f:
        for i in range(N): f.write("".join(str(int(x)) for x in A[i]) + "\n")
    json.dump({"n_data": n, "n_anc": k, "ancillas": [{"vertex": a, "qubit": q, "gate_index": p} for a, q, p in anc]}, open(out + "_meta.json", "w"))
    deg = A.sum(1); print(f"graph: {int(A.sum()) // 2} edges; degree min/mean/max {deg.min()}/{deg.mean():.1f}/{deg.max()}; "
                          f"data-data edges {int(A[:n, :n].sum()) // 2}, data-anc {int(A[:n, n:].sum())}, anc-anc {int(A[n:, n:].sum()) // 2}")
    print(f"rank A[anc, data] = {gf2_rank(A[n:, :n])}  (cut rank between all ancillas and all data)")
    orders = {}
    orders["natural (data 0..69, then ancillas in time order)"] = list(range(N))
    orders["ancillas in time order, then data 0..69"] = [a for a, _, _ in anc] + list(range(n))
    by_q = {q: [] for q in range(n)}
    for a, q, p in anc: by_q[q].append(a)
    sweep = []
    for q in range(n): sweep += by_q[q] + [q]
    orders["spatial sweep: qubit 0..69, each qubit's ancillas (time order) then its data vertex"] = sweep
    sweep2 = []
    for q in range(n): sweep2 += [q] + by_q[q]
    orders["spatial sweep, data vertex first in each group"] = sweep2
    res = {}
    for name, o in orders.items():
        pr = cut_profile(A, o); w = max(pr); res[name] = {"width": w, "log2_sum_2r": float(np.log2(np.sum(2.0 ** np.array(pr))))}
        print(f"{name}: width {w}, log2 sum 2^r = {res[name]['log2_sum_2r']:.2f}, profile max at step {int(np.argmax(pr)) + 1}")
    json.dump(res, open(out + "_widths.json", "w"), indent=1)
