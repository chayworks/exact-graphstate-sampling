"""Linear XEB of our ring-instance (#245) samples against IBM's hardware shots on the same circuit.

Usage: python ring64_xeb_compare.py <samples.txt> <experiment1_probabilities.csv>
samples.txt: run-14 sampler output (lines "ACC <bits> <log2 Pr> z=... <secs>"; bits indexed by qubit, character k = qubit k,
  i.e. the reverse of Qiskit order; rejected proposals are "rej" lines and are not samples).
csv: Renaud's Zenodo 22882048 file, column two_n_p = 2^64 Pr(shot) for each of IBM's 1,389 hardware shots.
Linear XEB = mean(2^n Pr) - 1; standard error = sd / sqrt(N). Both sides use exact probabilities of the same circuit
(ours from our sampler's own evaluation, the hardware's from Renaud's independent code; the two agree on the three gate
shots to 3.8e-6 in log2)."""
import sys, csv, math, statistics as st

def summary(v):
    n = len(v); m = st.mean(v); se = st.stdev(v) / math.sqrt(n); return n, m - 1, se

sp, cp = sys.argv[1:3]
ours = [2 ** (64 + float(l.split()[2])) for l in open(sp) if l.startswith("ACC")]
hw = [float(r["two_n_p"]) for r in csv.DictReader(open(cp))]
n1, x1, s1 = summary(ours); n2, x2, s2 = summary(hw)
d = x1 - x2; sd = math.hypot(s1, s2)
print(f"classical samples (ours): N = {n1}, XEB = {x1:.3f} +/- {s1:.3f}")
print(f"IBM hardware shots      : N = {n2}, XEB = {x2:.3f} +/- {s2:.3f}")
print(f"difference {d:.3f} +/- {sd:.3f} ({d / sd:.1f} sigma)")
