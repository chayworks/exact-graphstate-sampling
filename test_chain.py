"""gsamp2 chain mode: per-step conditionals Pr(b_k | prefix) and log2 Pr must equal the Python in-place
reference (which computes both branch norms explicitly) on random strings of small circuits, and log2 Pr
must equal gsamp2's two-pass prob mode. Kills: the single-pass identity n0 + n1 = 2 * previous norm (if it
were false, the conditionals would disagree), the printed step order, and the unnormalised-storage scaling.
Negative control: a plan with the M rows swapped must disagree."""
import random, subprocess, numpy as np
from inplace import sample_chain_probs
from test_gsamp import build, docker
from test_gsamp2 import plan_for
from inplace import write_inplace_plan

rng = random.Random(31); cmds = []; cases = []
for i in range(12):
    n = rng.randint(2, 9); T, tail, psi = build(n, rng, rng.randint(10, 80))
    order = list(range(n)); rng.shuffle(order)
    am, st, R, M = plan_for(n, T, tail, order); write_inplace_plan(f"work/gc{i}.txt", st, R, M, n)
    strings = ["".join(str(rng.randint(0, 1)) for _ in range(n)) for _ in range(6)]
    open(f"work/gc{i}_str.txt", "w").write("\n".join(strings) + "\n")
    if i < 3: _, stb, Rb, Mb = plan_for(n, T, tail, order, bad=True); write_inplace_plan(f"work/gcbad{i}.txt", stb, Rb, Mb, n)
    cases.append((st, R, M, strings))
    cmds += [f"echo CHAIN{i}; ./gsamp2 work/gc{i}.txt chain work/gc{i}_str.txt", f"echo PROB{i}; ./gsamp2 work/gc{i}.txt prob work/gc{i}_str.txt"]
    if i < 3: cmds.append(f"echo BAD{i}; ./gsamp2 work/gcbad{i}.txt chain work/gc{i}_str.txt")
blocks = {}; key = None
for line in docker("; ".join(cmds)).split("\n"):
    if line[:1].isalpha(): key = line; blocks[key] = []
    elif line.strip(): blocks[key].append([float(x) for x in line.split()])
worst_c = worst_p = 0; ctrl = 0; nfinite = 0
for i, (st, R, M, strings) in enumerate(cases):
    for j, s in enumerate(strings):
        bits = [int(ch) for ch in s]; ref = sample_chain_probs(st, max(R, 1), M, bits)
        got = blocks[f"CHAIN{i}"][j]; lp, cond = got[0], got[1:]
        assert len(cond) == len(bits)
        worst_c = max(worst_c, np.abs(np.array(cond) - np.array(ref)).max())
        pref = float(np.prod(ref)); lp_prob = blocks[f"PROB{i}"][j][0]
        if pref > 1e-12:
            nfinite += 1; worst_p = max(worst_p, abs(lp - np.log2(pref)), abs(lp - lp_prob))
        else: assert lp == -np.inf or lp < -35, (lp, pref)
        if i < 3: ctrl = max(ctrl, np.abs(np.array(blocks[f"BAD{i}"][j][1:]) - np.array(ref)).max())
print(f"chain mode: 12 circuits x 6 strings; max |cond - ref| = {worst_c:.2e}; max |log2P - ref|, |log2P - prob mode| = {worst_p:.2e} "
      f"over {nfinite} strings with Pr > 1e-12; negative control (M rows swapped) = {ctrl:.3f}")
assert worst_c < 1e-6 and worst_p < 1e-6 and ctrl > 1e-2 and nfinite >= 40
