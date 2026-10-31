"""gsamp (out-of-place) chain mode: per-step conditionals and log2 Pr must equal the Python in-place
reference (which computes both branch norms explicitly) and gsamp's two-pass prob mode, on random strings
of small circuits. Kills: the single-pass identity n0 + n1 = 2 * previous norm, the printed step order and
the unnormalised-storage scaling. Negative control: M rows swapped in the plan must disagree."""
import random, numpy as np
from amp import Amplitudes
from engine import write_sampler_plan
from inplace import plan_inplace, sample_chain_probs
from test_gsamp import build, docker

rng = random.Random(47); cmds = []; cases = []
for i in range(12):
    n = rng.randint(2, 9); T, tail, psi = build(n, rng, rng.randint(10, 80))
    order = list(range(n)); rng.shuffle(order)
    am = Amplitudes(n, T, tail, order); write_sampler_plan(f"work/gsc{i}.txt", am)
    R = max(1, max(max(s["r"], s["rn"]) for s in am.steps)); st = plan_inplace(am.A, order, R)
    strings = ["".join(str(rng.randint(0, 1)) for _ in range(n)) for _ in range(6)]
    open(f"work/gsc{i}_str.txt", "w").write("\n".join(strings) + "\n")
    if i < 3:
        amb = Amplitudes(n, T, tail, order); amb.M = [m[::-1].copy() for m in am.M]; write_sampler_plan(f"work/gscbad{i}.txt", amb)
    cases.append((st, R, am.M, strings))
    cmds += [f"echo CHAIN{i}; ./gsamp work/gsc{i}.txt chain work/gsc{i}_str.txt", f"echo PROB{i}; ./gsamp work/gsc{i}.txt prob work/gsc{i}_str.txt"]
    if i < 3: cmds.append(f"echo BAD{i}; ./gsamp work/gscbad{i}.txt chain work/gsc{i}_str.txt")
blocks = {}; key = None
for line in docker("; ".join(cmds)).split("\n"):
    if line[:1].isalpha(): key = line; blocks[key] = []
    elif line.strip(): blocks[key].append([float(x) for x in line.split()])
worst_c = worst_p = ctrl = 0; nfinite = 0
for i, (st, R, M, strings) in enumerate(cases):
    for j, s in enumerate(strings):
        bits = [int(ch) for ch in s]; ref = sample_chain_probs(st, R, M, bits)
        got = blocks[f"CHAIN{i}"][j]; lp, cond = got[0], got[1:]; assert len(cond) == len(bits)
        worst_c = max(worst_c, np.abs(np.array(cond) - np.array(ref)).max())
        pref = float(np.prod(ref)); lp_prob = blocks[f"PROB{i}"][j][0]
        if pref > 1e-12: nfinite += 1; worst_p = max(worst_p, abs(lp - np.log2(pref)), abs(lp - lp_prob))
        else: assert lp == -np.inf or lp < -35, (lp, pref)
        if i < 3: ctrl = max(ctrl, np.abs(np.array(blocks[f"BAD{i}"][j][1:]) - np.array(ref)).max())
print(f"gsamp chain mode: 12 circuits x 6 strings; max |cond - ref| = {worst_c:.2e}; max |log2P - ref|, |log2P - prob| = {worst_p:.2e} "
      f"over {nfinite} strings with Pr > 1e-12; negative control (M rows swapped) = {ctrl:.3f}")
# float32 tables: log2 Pr carries ~1e-6 of accumulated rounding (the two-pass prob mode itself agrees with the
# cloud run to 1e-6); a wrong scale or step order gives O(1) differences, as the control shows.
assert worst_c < 1e-6 and worst_p < 5e-6 and ctrl > 1e-2 and nfinite >= 40
