"""Independent recomputation of more nq70 samples (chunked, resumable), files named by target label.
Targets: cloud samples with indices given; control = sample 1 with 10 fixed step-bits flipped (a large,
deterministic change, unlike the 1-bit control in xcheck70c, which happened to move Pr by only 1.6%).

make:    python xcheck70d.py make <qasm> <samples> <idx,idx,...> <nch>
run:     python xcheck70d.py run            (runs every chunk without a .out, in order)
combine: python xcheck70d.py combine <samples>
"""
import sys, json, math, glob, os, subprocess
FLIPS = [3, 11, 19, 27, 35, 43, 51, 59, 64, 68]
if sys.argv[1] == "make":
    from tableau import clifford_body
    from amp import Amplitudes
    from engine import Sliced, write_plan
    F, SF, idxs, nch = sys.argv[2], sys.argv[3], [int(i) for i in sys.argv[4].split(",")], int(sys.argv[5])
    n, T, tail = clifford_body(F); order = json.load(open("order_best2.json")); am = Amplitudes(n, T, tail, order)
    steps_v = [st["v"] for st in am.steps]; assert steps_v == json.load(open("cloud/nq70_step_qubits.json"))
    sl = Sliced(am.A, order, [4, 60, 23, 0, 41, 2, 18, 3, 24]); write_plan("work/x70e_plan.txt", sl.steps, sl.m)
    rows = [l.split() for l in open(SF) if l.strip()]
    targets = {f"s{i}": rows[i][0] for i in idxs}
    c = list(rows[0][0])
    for p in FLIPS: c[p] = "1" if c[p] == "0" else "0"
    pass
    for lab, bits in targets.items():
        b = [0] * n
        for kk, q in enumerate(steps_v): b[q] = int(bits[kk])
        sls = list(sl.slices(am.weights(b))); per = len(sls) // nch
        for ch in range(nch):
            part = sls[ch * per:(ch + 1) * per] if ch < nch - 1 else sls[ch * per:]
            with open(f"work/x70e_{lab}_{ch}.txt", "w") as fh:
                fh.write(f"1\n{len(part)}\n")
                for pre, ww in part:
                    fh.write("%.17g %.17g\n" % (pre.real, pre.imag))
                    for a, d in ww: fh.write("%.17g %.17g %.17g %.17g\n" % (a.real, a.imag, d.real, d.imag))
    json.dump(targets, open("work/x70e_targets.json", "w")); print("targets:", list(targets))
elif sys.argv[1] == "run":
    host = os.path.abspath(".").replace("\\", "/")
    for job in sorted(glob.glob("work/x70e_*_[0-9].txt")):
        out = job[:-4] + ".out"
        if os.path.exists(out): continue
        r = subprocess.run(["docker", "run", "--rm", "--memory", "3g", "-v", f"{host}:/w", "-w", "/w", "-e", "OMP_NUM_THREADS=16",
                            "gcc:14", "./dp", "work/x70e_plan.txt", job.replace("\\", "/")], capture_output=True, text=True,
                           env={**os.environ, "MSYS_NO_PATHCONV": "1"})
        assert r.returncode == 0, r.stderr
        open(out + ".tmp", "w").write(r.stdout); os.replace(out + ".tmp", out); print("done", job, flush=True)
else:
    rows = [l.split() for l in open(sys.argv[2]) if l.strip()]; targets = json.load(open("work/x70e_targets.json"))
    for lab in targets:
        parts = sorted(glob.glob(f"work/x70e_{lab}_*.out")); z = sum(complex(*map(float, open(p).read().split())) for p in parts)
        rep = f"gsamp log2 Pr = {float(rows[int(lab[1:])][1]):.7f}" if lab.startswith("s") else "control: sample 0 with 10 step-bits flipped"
        print(f"{lab:7s} {rep}   dp.cpp log2 |Z|^2/2^70 = {math.log2(abs(z) ** 2 / 2 ** 70):.7f}   ({len(parts)} chunks)")
