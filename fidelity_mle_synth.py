"""Does ibm_fidelity_mle.py recover a known fidelity? Synthetic shots: with probability F a shot is drawn from the
ideal (Porter-Thomas) law, z ~ Gamma(2, 1) with density z e^{-z}, else from the uniform law, z ~ Exp(1); log2 p =
log2 z - 70 is written in the format of results/ibmxeb/ibm_out.txt and the estimator is run unchanged with the real
exact samples supplying E[z]. Usage: python fidelity_mle_synth.py results/samples_nq70.txt  -> results/fidelity_mle_synthetic.txt"""
import sys, os, subprocess, re, numpy as np
samples = sys.argv[1]; rng = np.random.default_rng(20261031); N = 1726; out = []
for F in (0.0, 0.4, 1.0):
    ideal = rng.random(N) < F; z = np.where(ideal, rng.gamma(2.0, 1.0, N), rng.exponential(1.0, N))
    os.makedirs("work", exist_ok=True); path = f"work/synth_{F}.txt"
    with open(path, "w") as f:
        for v in z: f.write(f"{np.log2(v) - 70:.12f}\n")
    r = subprocess.run([sys.executable, "ibm_fidelity_mle.py", path, samples], capture_output=True, text=True)
    if "F_mle" not in r.stdout: print("estimator output:", r.stdout[-500:], r.stderr[-800:]); raise SystemExit(1)
    xeb = re.search(r"F_xeb = XEB/\(E\[z\]-1\)\s+([-\d.]+) \+/- ([-\d.]+)", r.stdout); mle = re.search(r"F_mle \(profile 95%\)\s+([-\d.]+)\s+\[([-\d.]+), ([-\d.]+)\]", r.stdout)
    if xeb is None or mle is None: print("UNPARSED estimator output:", r.stdout); raise SystemExit(1)
    line = f"true F = {F:.1f}: F_xeb = {xeb.group(1)} +/- {xeb.group(2)}; F_mle = {mle.group(1)} [{mle.group(2)}, {mle.group(3)}]"
    print(line); out.append(line)
with open("results/fidelity_mle_synthetic.txt", "w") as f:
    f.write("# python fidelity_mle_synth.py results/samples_nq70.txt   (N = 1726 synthetic shots per F, seed 20261031)\n" + "\n".join(out) + "\n")
