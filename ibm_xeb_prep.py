"""Prepare IBM's 1,726 hardware samples (#151) for exact scoring with gsamp2 (prob mode).
Qiskit strings are little-endian: character k is qubit 69-k. Writes three string files in the plan's step order:
  ibm_qiskit.txt   (correct order), ibm_reversed.txt (bit-reversed control, first 200),
  uniform.txt      (128 uniform random strings, seed 20261030, control)."""
import sys, json, random, numpy as np
from inplace import plan_inplace
a = np.load(sys.argv[1], allow_pickle=False); steps_v = json.load(open(sys.argv[2])); out = sys.argv[3]
n = 70; assert a.shape == (1726,) and all(len(s) == n and set(s) <= {"0", "1"} for s in a)
def write(path, strings, mapper):
    with open(path, "w") as f:
        for s in strings:
            b = mapper(s); f.write("".join(str(b[q]) for q in steps_v) + "\n")
write(f"{out}/ibm_qiskit.txt", a, lambda s: [int(s[n - 1 - q]) for q in range(n)])
write(f"{out}/ibm_reversed.txt", a[:200], lambda s: [int(s[q]) for q in range(n)])
rng = random.Random(20261030); uni = ["".join(rng.choice("01") for _ in range(n)) for _ in range(128)]
write(f"{out}/uniform.txt", uni, lambda s: [int(c) for c in s])
print("wrote", len(a), 200, 128)
