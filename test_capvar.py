"""The grouped second-moment formula of capvar.py against brute force over every prefix string, on random small
circuits with the #151 tail (same generator as test_gsamp), with random prefixes of every size. Kills: a wrong pairing
rule (supports must be equal), a dropped sign, a wrong Bloch component, a wrong normalisation. Negative control: one
element's sign flipped changes the grouped value whenever that element shares its support with another element or
is the only one (its square is unchanged only if its group sum is 0 or it is alone), so the control is reported as the
largest change seen. Also exports one case for capvar.cpp and prints its expected value."""
import random, itertools, numpy as np
from test_gsamp import build
from stabprefix import prefix_prob_factors, stab_group_on
from capvar import grouped_Ez2, export_group
rng = random.Random(11); worst = 0.0; ctrl = 0.0; ncase = 0; nontriv = 0; exported = None
for case in range(40):
    n = rng.randint(3, 8); T, tail, psi = build(n, rng, rng.randint(5, 60))
    per = [[] for _ in range(n)]
    for name, qs in tail: per[qs[0]].append((name, qs))
    from stabprefix import bloch; from amp import tail_matrix
    blochs = [bloch(tail_matrix(per[q])) for q in range(n)]
    if not all(np.allclose(np.abs(b), 1 / np.sqrt(3), atol=1e-9) for b in blochs): continue   # the formula's 3^-|U| needs the face basis
    K = rng.randint(1, n); P = rng.sample(range(n), K)
    B = np.array([[(i >> k) & 1 for k in range(K)] for i in range(2 ** K)], np.uint8)
    f, ng = prefix_prob_factors(T, blochs, P, B); brute = float(np.mean(f ** 2))      # f = 2^K Pr(b_P) = z(b)
    grouped, ng2, groups = grouped_Ez2(T, blochs, P); assert ng == ng2
    worst = max(worst, abs(brute - grouped)); ncase += 1; nontriv += ng > 1
    if ng > 1:
        flipped = max(abs(grouped_Ez2(T, blochs, P, flip=k)[0] - brute) for k in range(1, min(ng, 8)))
        ctrl = max(ctrl, flipped)
    if exported is None and ng >= 16: exported = (n, K, ng, brute); export_group(T, blochs, P, "work/capvar_test_gens.txt")
print(f"{ncase} random circuits (3-8 qubits, random prefixes), {nontriv} with a nontrivial prefix group: "
      f"max |E[z^2] brute - grouped| = {worst:.2e}; largest change from a flipped sign (control) = {ctrl:.3e}")
if exported: print(f"exported work/capvar_test_gens.txt: n = {exported[0]}, K = {exported[1]}, group size {exported[2]}, expected E[z^2] = {exported[3]:.12f}")
assert worst < 1e-10 and ctrl > 1e-3
