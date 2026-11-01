"""Driver for the doped instance: hybrid rejection sampling and exact scoring of given strings with gsamp2.

Conventions. Plan: an in-place plan of the 538-vertex gadget circuit whose last m steps are data vertices (the
deferred suffix); step_vertices.json gives the vertex of each step. Data strings are 70 bits indexed by qubit; IBM's
published shots are Qiskit little-endian (character k is qubit 69-k), the flag --qiskit converts them.
Ancilla bits are always 0. Pr(all 468 ancillas = 0) = 2^-468 exactly, so for a data string b:
    log2 Pr_doped(b) = log2 Pr_gadget(b, anc = 0) + 468.

sample: proposals draw the first 70-m data bits uniformly; one hybrid pass per proposal gives the exact marginal
Pr(prefix bits, anc = 0) (the suffix is traced out by orthonormality) and an exact draw of the last m bits. The
proposal is accepted with probability z / cap, z = 2^(70-m) * 2^468 * Pr(prefix, anc = 0) (the marginal of the
doped distribution in units of the uniform one). Any proposal with z > cap is reported; with m = 20 the marginal
averages 2^20 output probabilities and z stays within a fraction of a percent of 1.
Usage:
  python doped_sampler.py sample  <plan> <step_vertices.json> <nprop> <seed> <cap> <out.txt> [--threads T] [--m M]
  python doped_sampler.py score   <plan> <step_vertices.json> <strings.txt> <out.txt> [--qiskit] [--threads T]
Outputs are appended line by line so a killed run can be resumed by passing a new seed."""
import sys, json, os, subprocess, time, math, random
NA = 468; ND = 70
def argval(flag, default=None):
    if flag in sys.argv: i = sys.argv.index(flag); v = sys.argv[i + 1]; del sys.argv[i:i + 2]; return v
    return default
def run(cmd, threads):
    env = {**os.environ, "OMP_NUM_THREADS": str(threads)}; return subprocess.run(cmd, capture_output=True, text=True, env=env, check=True).stdout
if __name__ == "__main__":
    threads = int(argval("--threads", os.cpu_count())); BIN = argval("--bin", "./gsamp2"); mcap = argval("--m", None); qiskit = "--qiskit" in sys.argv
    if qiskit: sys.argv.remove("--qiskit")
    mode, plan, svj = sys.argv[1:4]; sv = json.load(open(svj)); ns = len(sv); pos = {v: i for i, v in enumerate(sv)}
    data_steps = [i for i, v in enumerate(sv) if v < ND]; m_suffix = 0
    while m_suffix < ns and sv[ns - 1 - m_suffix] < ND: m_suffix += 1
    if mcap is not None: m_suffix = min(m_suffix, int(mcap))      # --m: sample at most this many trailing data bits
    def to_step_string(bits70, free_suffix):
        s = ["0"] * ns
        for q in range(ND): s[pos[q]] = str(bits70[q])
        if free_suffix:
            for i in range(ns - m_suffix, ns): s[i] = "?"
        return "".join(s)
    if mode == "score":
        strings, out = sys.argv[4], sys.argv[5]; lines = [l.strip() for l in open(strings) if l.strip()]
        done = sum(1 for _ in open(out)) if os.path.exists(out) else 0
        with open(out, "a") as f:
            for li, s in enumerate(lines[done:], start=done):
                assert len(s) == ND and set(s) <= {"0", "1"}
                b = [int(s[ND - 1 - q]) for q in range(ND)] if qiskit else [int(c) for c in s]
                tmp = f"{out}.cur"; open(tmp, "w").write(to_step_string(b, False) + "\n"); t0 = time.time()
                lp = float(run([BIN, plan, "prob", tmp], threads).split()[0]) + NA
                f.write(f"{s} {lp:.9f} {time.time() - t0:.1f}\n"); f.flush()
    elif mode == "sample":
        nprop, seed, cap, out = int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), sys.argv[7]
        rng = random.Random(seed); props = [[rng.randrange(2) for _ in range(ND)] for _ in range(nprop)]; u = [rng.random() for _ in range(nprop)]
        K = ns - m_suffix; zmax = 0.0; acc = 0
        with open(out, "a") as f:
            f.write(f"# seed {seed} cap {cap} m_suffix {m_suffix} K {K} nprop {nprop}\n")
            for i in range(nprop):
                tmp = f"{out}.cur"; open(tmp, "w").write(to_step_string(props[i], True) + "\n"); t0 = time.time()
                r = run([BIN, plan, "hybrid", tmp, str(seed * 1000003 + i)], threads).split()
                bits_step, lp_full, cond = r[0], float(r[1]) + NA, [float(x) for x in r[2:]]
                lp_prefix = (sum(math.log2(c) for c in cond[:K]) + NA) if all(c > 0 for c in cond[:K]) else -math.inf
                z = 2.0 ** (lp_prefix + (ND - m_suffix)); zmax = max(zmax, z); ok = u[i] < z / cap
                b = "".join(bits_step[pos[q]] for q in range(ND))
                f.write(f"{'ACC' if ok else 'rej'} {b} {lp_full:.9f} z={z:.6f} {time.time() - t0:.1f}s{' Z_ABOVE_CAP' if z > cap else ''}\n"); f.flush(); acc += ok
            f.write(f"# done: {acc} accepted of {nprop}; max z {zmax:.6f}; cap {cap}\n")
        print(f"{acc} accepted of {nprop}; max z {zmax:.6f} (cap {cap})")
    else: raise SystemExit("mode must be sample or score")
