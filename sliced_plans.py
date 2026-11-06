"""Sliced in-place plans. Slicing vertex v fixes x_v: the Ising sum splits into the two values of x_v, each a sum on the
graph with v's edges removed, with weights w_v(x_v) kept and the phases (-1)^{x_v x_u} folded into M_u[:, 1] for every
neighbour u. All 2^{|S|} slice assignments share the same graph A' and therefore the same plan structure (addresses,
step sizes); only the M matrices differ, so their dumped tables add up (gsamp2 prefix/resume).

Neighbours in the suffix (steps >= K, the sampled part of a hybrid run) cannot carry the folded phase, because the
suffix steps run once on the summed table. Instead the phase (-1)^{x_v sum_{u in N(v), u suffix} x_u} is applied at
the summation. The table at step K is indexed by logical coordinates u whose field vector phi(u) = sum_i u_i B_i (B_i
the cut basis of the old coordinates, vectors over the suffix) is the linear functional the processed part imposes on
x_suffix: T(u) collects the weight of all prefix assignments with functional phi(u), and the partial sum is
F(x_S) = sum_u T(u) (-1)^{<phi(u), x_S>}. Multiplying F by (-1)^{<s_v, x_S>} therefore translates the table:
T(u) -> T(u + a) with s_v = sum_i a_i B_i, which exists when s_v lies in span(B). The 'twist' file gives resume that
translation per slice (as a shift of the kept-coordinate index and of the coset offset). At the summation cut of a
hybrid run (rank <= |suffix|) the span is usually everything, so dense graphs can be sliced.

Usage: python sliced_plans.py <circuit.qasm> <order.json> <out_prefix> [v1,v2,...] [--nd ND] [--K K]
With no slice list: scores candidate vertices and prints the best; with a list: writes <out_prefix>_s<assignment>.txt
plans, <out_prefix>_slices.json and, when a slice vertex has suffix neighbours, <out_prefix>_twist.txt.
--nd: number of data vertices (default 70; vertices >= ND are ancillas). --K: summation step (default: first step
of the pure-data suffix; for a hybrid run on an undoped instance give K = steps minus the sampled suffix)."""
import sys, json, itertools, numpy as np
from tableau import clifford_body, gf2_rank
from amp import Amplitudes
from inplace import plan_inplace, write_inplace_plan

def profile(A, order, ks=None):
    n = len(A); ks = ks or range(1, n); return {k: gf2_rank(A[np.ix_(order[k:], order[:k])]) for k in ks}

def gf2_solve(B, s):
    """a with XOR_i a_i B_i = s (python ints as bit vectors), or None."""
    s = int(s); rows = [(int(b), 1 << i) for i, b in enumerate(B)]; basis = []            # (vector, combination)
    for v, comb in rows:
        for bv, bc in basis:
            if (v >> (bv.bit_length() - 1)) & 1: v ^= bv; comb ^= bc
        if v: basis.append((v, comb)); basis.sort(key=lambda p: -p[0])
    comb = 0
    for bv, bc in basis:
        if (s >> (bv.bit_length() - 1)) & 1: s ^= bv; comb ^= bc
    return comb if s == 0 else None

def make_sliced_plans(A, M, order, S, out, K=None, ND=70):
    """Writes the plans for every assignment of the slice vertices S. Returns dict(R, K, twist, step_vertices)."""
    N = len(A); A2 = A.copy()
    for v in S: A2[v, :] = 0; A2[:, v] = 0
    R = max(1, max(profile(A2, order).values())); st = plan_inplace(A2, order, R); sv = [s["v"] for s in st]
    if K is None:
        K = next(i for i in range(len(sv)) if all(v < ND for v in sv[i:]))
        if K == 0: K = len(sv)
    pre = set(sv[:K]); suf = set(sv[K:])
    for v in S: assert v in pre, f"slice vertex {v} lies in the suffix"
    # per slice vertex: the suffix part of its neighbourhood as a vertex bitmask, and its expansion in the cut basis at K
    avec = {}
    if K < len(st):
        B = st[K]["B"]; keep = st[K]["keep"]; Tq = st[K]["Tq"]; Qall = st[K]["Qall"]
    for v in S:
        s_v = sum(1 << int(u) for u in np.nonzero(A[v])[0] if u in suf)       # python ints: vertices can exceed 63
        if s_v:
            assert K < len(st), "suffix neighbours but K is past the end"
            a = gf2_solve(B, s_v); assert a is not None, f"slice vertex {v}: suffix neighbourhood not in the span of the cut at step {K}"
            avec[v] = a
    twist = []
    for bits in itertools.product([0, 1], repeat=len(S)):
        Ms = [m.copy() for m in M]; a = 0
        Sset = set(S)
        for v, xv in zip(S, bits):
            Ms[v][:, 1 - xv] = 0                                   # only x_v = xv contributes
            if xv:
                for u in np.nonzero(A[v])[0]:
                    if u in pre and u not in Sset: Ms[u][:, 1] *= -1   # (-1)^{x_v x_u} for prefix neighbours
                a ^= avec.get(v, 0)                                # suffix neighbours: sign at the summation
        # an edge between two slice vertices contributes (-1)^{x_i x_j} once (applying it from both ends cancels it)
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                if A[S[i], S[j]] and bits[i] and bits[j]: Ms[S[j]][:, 1] *= -1
        tag = "".join(map(str, bits)); write_inplace_plan(f"{out}_s{tag}.txt", st, R, Ms, N)
        if avec:
            # translation by a in logical coordinates: a = a_keep + a_Q (F_2^r = span(e_keep) + Q); the dump entry (c, tq)
            # of this slice is read from (c ^ cshift, tq ^ tqshift)
            keepmask = sum(1 << i for i in keep)
            aq = next(qv for qv in Qall if (a ^ qv) & ~keepmask == 0); ak = a ^ aq
            cshift = sum(((ak >> keep[j]) & 1) << j for j in range(len(keep))); twist.append((cshift, Tq[aq]))
    if avec:
        with open(f"{out}_twist.txt", "w") as f:
            for cshift, tqshift in twist: f.write(f"{cshift} {tqshift}\n")
    json.dump({"slices": S, "R": R, "K_suffix_start": K, "twist": bool(avec), "step_vertices": sv}, open(f"{out}_slices.json", "w"))
    return dict(R=R, K=K, twist=bool(avec), step_vertices=sv)

def _opt(flag, default):
    if flag in sys.argv: i = sys.argv.index(flag); v = sys.argv[i + 1]; del sys.argv[i:i + 2]; return int(v)
    return default

if __name__ == "__main__":
    ND = _opt("--nd", 70); KOPT = _opt("--K", None)
    qasm, oj, out = sys.argv[1:4]; order = json.load(open(oj))
    N, T, tail = clifford_body(qasm); am = Amplitudes(N, T, tail, order); A = am.A.copy(); M = [m.copy() for m in am.M]
    pr = [s["r"] for s in am.steps]; w = max(pr); hi = [k for k, r in enumerate(pr) if r >= w - 1]
    print(f"N = {N}, width {w}, widest cuts at steps {hi[0]}..{hi[-1]} ({len(hi)} steps)")
    if len(sys.argv) < 5:
        cand = sorted(set(order[max(0, hi[0] - 40):min(N, hi[-1] + 40)])); ks = list(range(max(1, hi[0] - 5), min(N, hi[-1] + 6)))
        best = []
        for v in cand:
            A2 = A.copy(); A2[v, :] = 0; A2[:, v] = 0; p2 = profile(A2, order, ks); best.append((max(p2.values()), v))
        best.sort(); print("best single slices (new max over the widest region, vertex):", best[:12])
        v1 = best[0][1]; A2 = A.copy(); A2[v1, :] = 0; A2[:, v1] = 0; best2 = []
        for v in cand:
            if v == v1: continue
            A3 = A2.copy(); A3[v, :] = 0; A3[:, v] = 0; p3 = profile(A3, order, ks); best2.append((max(p3.values()), v))
        best2.sort(); print(f"best second slice with {v1}:", best2[:8])
        raise SystemExit
    S = [int(x) for x in sys.argv[4].split(",")]
    info = make_sliced_plans(A, M, order, S, out, KOPT, ND)
    print(f"sliced graph: width {info['R']} (was {w}); wrote {2 ** len(S)} slice plans; summation at step K = {info['K']} of {len(info['step_vertices'])}; "
          f"twist {'written' if info['twist'] else 'not needed'}; table 2^{info['R']} complex64 = {2 ** info['R'] * 8 / 1e9:.1f} GB")
