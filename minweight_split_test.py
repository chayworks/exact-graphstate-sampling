"""Validate the three-way split of minweight_exact.py on instances whose minimum weight is known from the direct
enumeration (results/minweight_scan.txt): run the split with a small t (so the minimiser is NOT in the directly
enumerated part), combine the three cases, and check the minimum equals the direct answer. Also reports which case
attains it, so the mutation "drop that case" is seen to change the answer.
Usage: python minweight_split_test.py <A.txt> <t> <W> <direct_dmin>   (runs the C++ enumerators through docker)"""
import sys, os, subprocess, re
A, t, W, direct = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); base = A.replace("_A.txt", "")
env = {**os.environ, "MSYS_NO_PATHCONV": "1", "PYTHONIOENCODING": "utf-8"}
def sh(cmd): return subprocess.run(cmd, capture_output=True, text=True, env=env, encoding="utf-8").stdout
py = sh([sys.executable, "-P", "minweight_exact.py", A, str(t), str(W)])
dk = ["docker", "run", "--rm", "-v", "D:/GIT/Innovation/graphsamp:/w", "-w", "/w", "gcc:14", "bash", "-c"]
part1 = sh(dk + [f"OMP_NUM_THREADS=16 ./work/minweight_bin {A} {t}"])
case1 = sh(dk + [f"OMP_NUM_THREADS=16 ./work/minweight_case1_bin {base}_G.txt {base}_ker.txt {t}"])
m1 = int(re.search(r"minimum over popcount\(v\) <= \d+: weight (\d+)", part1).group(1))
mc1 = int(re.search(r"case \(1\): minimum .* is weight (\d+)", case1).group(1))
mk = re.search(r"elements are stabilizers X\^k of weight \[([^\]]*)\]", py); kw = [int(x) for x in mk.group(1).split(",")] if mk and mk.group(1).strip() else []
mc2 = int(re.search(r"case \(2\): .* weight (\d+) at", py).group(1))
split = min([m1, mc1, mc2] + kw)
print(f"{A}: t = {t}, W = {W}; |v| <= t part: {m1}; case (1) |Av| <= t: {mc1}; kernel elements: {kw}; case (2) |e| <= {2 * (W - t - 1)}: {mc2}")
print(f"split minimum {split} vs direct d_min {direct}: {'AGREE' if split == direct else 'DISAGREE'}; attained by: "
      + ", ".join(nm for nm, v in [("|v|<=t", m1), ("case 1", mc1), ("case 2", mc2)] + [("kernel", min(kw)) if kw else ("kernel", 10**9)] if v == split))
print(f"mutation check: without case 2 the split would report {min([m1, mc1] + kw)}; without case 1 {min([m1, mc2])}; without the |v|<=t part {min([mc1, mc2] + kw)}")
