"""Chunked version of xcheck70: writes, for each target string, its 2^s slices split into `nch` job files
(work/x70c_<target>_<chunk>.txt), so each dp.cpp call is a few minutes in the foreground. Targets: the first
`k` cloud samples, plus sample 1 with step-bit 35 flipped (negative control). `combine` sums the chunk
outputs and prints gsamp's log2 Pr next to dp.cpp's |Z|^2/2^70.

Usage: python xcheck70c.py make <qasm> <samples> <k> <s> <nch>
       python xcheck70c.py combine <samples> <k>
"""
import sys, json, math, glob, itertools, numpy as np
from engine import Sliced, write_plan

if sys.argv[1] == "make":
    from tableau import clifford_body
    from amp import Amplitudes
    F, SF, k, s, nch = sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    n, T, tail = clifford_body(F); order = json.load(open("order_best2.json")); am = Amplitudes(n, T, tail, order)
    steps_v = [st["v"] for st in am.steps]
    assert steps_v == json.load(open("cloud/nq70_step_qubits.json")), "step order differs from the cloud plan"
    S = [4, 60, 23, 0, 41, 2, 18, 3, 24, 32][:s]          # greedy slice list printed by xcheck70.py for this order
    sl = Sliced(am.A, order, S); maxr = write_plan("work/x70c_plan.txt", sl.steps, sl.m)
    rows = [l.split() for l in open(SF) if l.strip()][:k]
    targets = [r[0] for r in rows]; f = list(rows[0][0]); f[35] = "1" if f[35] == "0" else "0"; targets.append("".join(f))
    for ti, bits in enumerate(targets):
        b = [0] * n
        for kk, q in enumerate(steps_v): b[q] = int(bits[kk])
        sls = list(sl.slices(am.weights(b))); per = len(sls) // nch
        for c in range(nch):
            part = sls[c * per:(c + 1) * per] if c < nch - 1 else sls[c * per:]
            with open(f"work/x70c_{ti}_{c}.txt", "w") as fh:
                fh.write(f"1\n{len(part)}\n")
                for pre, ww in part:
                    fh.write("%.17g %.17g\n" % (pre.real, pre.imag))
                    for a, d in ww: fh.write("%.17g %.17g %.17g %.17g\n" % (a.real, a.imag, d.real, d.imag))
    print(f"sliced {S}: width {maxr}; {len(targets)} targets x {nch} chunks of ~{2 ** s // nch} slices")
else:
    SF, k = sys.argv[2], int(sys.argv[3]); rows = [l.split() for l in open(SF) if l.strip()][:k]
    for ti in range(k + 1):
        parts = sorted(glob.glob(f"work/x70c_{ti}_*.out"))
        zs = [complex(*map(float, open(p).read().split())) for p in parts]
        z = sum(zs); lab = f"seed {rows[ti][2]}: gsamp log2 Pr = {float(rows[ti][1]):.7f}" if ti < k else "negative control (sample 1, step-bit 35 flipped):"
        print(f"{lab}   dp.cpp log2 |Z|^2/2^70 = {math.log2(abs(z) ** 2 / 2 ** 70):.7f}   ({len(parts)} chunks)")
