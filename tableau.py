"""Stabilizer tableau of the Clifford body of a QASM circuit (h, s, sdg, sx, x, z, cz, t-free).

Rows are the n stabilizer generators of C|0>; xs, zs are (n, n) uint8 arrays, r the sign bits.
Conjugation rules (Aaronson-Gottesman), phases tracked; sx = Rx(pi/2) up to global phase.
"""
import re
import numpy as np

def parse(path):
    n = None; ops = []
    for line in open(path):
        line = line.strip().rstrip(";")
        if not line or line.startswith(("OPENQASM", "include")): continue
        m = re.match(r"qreg q\[(\d+)\]", line)
        if m: n = int(m.group(1)); continue
        name, args = line.split(" ", 1)
        qs = [int(a) for a in re.findall(r"q\[(\d+)\]", args)]
        ops.append((name, qs))
    return n, ops

class Tableau:
    def __init__(s, n):
        s.n = n; s.xs = np.zeros((n, n), np.uint8); s.zs = np.eye(n, dtype=np.uint8); s.r = np.zeros(n, np.uint8)
    def h(s, q):
        s.r ^= s.xs[:, q] & s.zs[:, q]; s.xs[:, q], s.zs[:, q] = s.zs[:, q].copy(), s.xs[:, q].copy()
    def s_(s, q):
        s.r ^= s.xs[:, q] & s.zs[:, q]; s.zs[:, q] ^= s.xs[:, q]
    def z(s, q): s.r ^= s.xs[:, q]
    def x(s, q): s.r ^= s.zs[:, q]
    def sdg(s, q): s.z(q); s.s_(q)
    def sx(s, q): s.h(q); s.s_(q); s.h(q)          # H S H = sqrt(X) up to global phase
    def cz(s, a, b):
        s.r ^= s.xs[:, a] & s.xs[:, b] & (s.zs[:, a] ^ s.zs[:, b])
        s.zs[:, a] ^= s.xs[:, b]; s.zs[:, b] ^= s.xs[:, a]
    def apply(s, name, qs):
        {"h": s.h, "s": s.s_, "sdg": s.sdg, "sx": s.sx, "sxdg": s.sxdg, "x": s.x, "z": s.z}[name](*qs) if name != "cz" else s.cz(*qs)
    def sxdg(s, q): s.h(q); s.sdg(q); s.h(q)      # H S^dag H = sqrt(X)^dag up to global phase

def clifford_body(path):
    """Return (n, tableau of the leading Clifford part, remaining non-Clifford tail ops)."""
    n, ops = parse(path)
    cl = {"h", "s", "sdg", "sx", "sxdg", "x", "z", "cz"}
    T = Tableau(n); i = 0
    for i, (name, qs) in enumerate(ops):
        if name not in cl: break
        T.apply(name, qs)
    else:
        i = len(ops)
    return n, T, ops[i:]

def gf2_rank(M):
    M = (np.asarray(M) & 1).astype(np.uint8).copy(); r = 0; rows, cols = M.shape
    for c in range(cols):
        p = np.nonzero(M[r:, c])[0]
        if len(p) == 0: continue
        p = p[0] + r; M[[r, p]] = M[[p, r]]
        nz = np.nonzero(M[:, c])[0]; nz = nz[nz != r]; M[nz] ^= M[r]
        r += 1
        if r == rows: break
    return r

def cut_rank(T, A):
    """Entanglement (in ebits) of the stabilizer state across A | complement: rank(G|_B) - |B|."""
    B = np.array(sorted(set(range(T.n)) - set(A)))
    G = np.concatenate([T.xs[:, B], T.zs[:, B]], axis=1)
    return gf2_rank(G) - len(B)
