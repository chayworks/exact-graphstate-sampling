# Exact Born sampling of a 70-qubit random-graph-state advantage experiment on one CPU node

*Draft technical note, not submitted. Numbers in braces {…} are pending; every other number comes from a
committed artifact listed in the appendix.*

## Abstract

IBM and the University of Chicago proposed sampling random graph states in a non-stabilizer product
basis as a verifiable quantum-advantage task. They ran a 70-qubit, depth-70 instance with an estimated
fidelity of 0.37 and extrapolated classical costs of 10²⁵ s (matrix product states) and 10⁷ s per sample
(extended stabilizer simulation). We show that this instance admits **exact sampling from the ideal Born
distribution** on a single 16-core CPU node: **13.0 s per sample**, with 1,012 samples drawn and two
complex64 tables of 2³¹ entries.

The method combines two standard ingredients:
- reduction to a graph state with a dynamic program whose state space is set by the cut rank (the
  rank-width simulability of graph-state measurements);
- a marginal identity that makes the right environment trivial: after the first t qubits are measured, the
  remaining states Z^v|G⟩ are orthonormal, so Pr(b₁…b_t) = 2^{-t} Σ|P_t|².

Each sample is therefore drawn in two passes, with no rejection and no truncation. This differs from
capped rejection sampling (#260 on the tracker), whose output distribution departs from the Born law.

We verify the sampler against dense statevectors on small instances, and on the full instance against an
independent amplitude algorithm. An in-place variant draws an exact sample in about 17 s on a laptop with
an 8.6 GB table.

We also score IBM's 1,726 published hardware shots exactly: linear XEB **0.466 ± 0.031** (controls at
zero), with every per-shot probability agreeing with #260's independent code to 10⁻⁵ in log₂. Two fidelity
estimators built on the exact probabilities give 0.46 ± 0.04 and 0.48 [0.43, 0.54], above the 0.37 that
IBM reports from direct fidelity estimation of the post-selected shots (95% lower bound 0.092), and far
narrower.

Finally, we map the **hardness frontier** of the construction. The width that governs the attack grows as n/2 − 3.5, so memory doubles every two qubits. At depth n the
attack needs one server at n = 70, a large supercomputer at n ≈ 90, and more than any existing machine
at n ≳ 100. That locates the size at which this construction escapes rank-width methods.

## 1. The task

Circuit `nq70_depth70_checks27_basis_fstate`, in the Quantum Advantage Tracker (issue #151):
1. Hadamards on 70 qubits in a line.
2. Brickwork of nearest-neighbour CZ gates to CZ depth 70, with a random √X or S√X on every qubit
   between layers.
3. On every qubit, T†·H·Rz(θ)·H with θ = arccos(1/√3), a rotation into the "face-state" basis that maximizes
   stabilizer extent.
4. Measurement in Z.

The hardware run (IBM Boston, Heron R3, with 27 ancillas for spacetime Pauli checks) reports a post-selected
fidelity of 0.37 from direct fidelity estimation of the underlying graph state, with a 95% lower bound of
0.092. 1,726 post-selected shots are published with the instance.

## 2. Method

**Graph form.** The Clifford part maps to a stabilizer tableau, and then, by Gaussian elimination and
local Cliffords, to |S⟩ = ⊗U_j† |G_A⟩. Then

  ⟨b|ψ⟩ = 2^{-n/2} Σ_x Π_j M_j[b_j, x_j] (−1)^{Σ_{i<j} A_ij x_i x_j},  where M_j = W_j U_j†.

**Cut-rank dynamic program.** Fix a vertex order. After the prefix P is processed, the future depends
only on s = A[F,P]·x_P ∈ GF(2)^{|F|}. That vector lives in a space of dimension r = rank A[F,P], so the
prefix is a table P[c] of 2^r complex numbers in a chosen basis of that space. Each step is a surjective
GF(2)-linear map (c, x_v) ↦ c′. We evaluate it in gather form: a `pdep` scatter plus a kernel of
dimension at most 2.

**Exact sampling.** With prefix P measured with outcomes b_P, the post-measurement state is

  2^{-|P|/2} Σ_c P[c] Z^{B c} |G_F⟩.

For v ≠ 0, Z^v has no X part, so it anticommutes with some graph stabilizer and ⟨G|Z^v|G⟩ = 0. The
states {Z^v|G_F⟩} are therefore orthonormal, and Pr(b_P) = 2^{-|P|} Σ_c |P[c]|². Measuring qubit by qubit
in the DP order:
1. compute both candidate extensions' norms (one read-only pass);
2. draw the outcome;
3. store the chosen extension (a second pass).

This is the graph-state analogue of perfect MPS sampling (Ferris & Vidal 2012). The output is the exact
Born distribution, up to floating-point error.

**Cost for this instance.** A greedy-plus-local-search order (`order2.py`) reaches maximum cut
rank 30, the value #260 also reports; the natural order gives 32. The production sampling run used an earlier
width-31 order and cost 4·Σ_k 2^{r_k} = 6.0·10¹⁰ table updates per sample with two complex64 tables of 2³¹
entries; the later runs (laptop, hardware-shot scoring) use the width-30 order.

## 3. Results

1,012 samples in 13,202 s on an AWS r7a.4xlarge (16 cores, AMD EPYC 9R14), peak RSS 33.6 GB, about $4.50
of on-demand time.
- All samples are distinct.
- Linear XEB against their own exact probabilities is 1.020 ± 0.047, i.e. 2⁷⁰ Σ p² ≈ 2.02,
  consistent with Porter–Thomas-like output.

**Laptop.** With an in-place variant of the program, the arena is a single complex64 table of 2^30
entries (8.6 GB) for a width-30 order. On a laptop (AMD Ryzen 9 9955HX, 16 threads, container memory
capped at 10 GB), 3 exact samples took 55 s including allocation, i.e. **about 17 s per exact sample**
(`results/laptop_samples_w30.txt`). The in-place sampler was validated on every bitstring of small
circuits. On the full instance it reproduces the cloud samples' log₂ Pr to 3·10⁻⁷. The independent amplitude engine confirms the first laptop sample:
−68.5019327 against −68.5019326 reported (`results/xcheck70e_laptop.txt`).

**Scoring IBM's hardware shots.** The sampler's probability mode evaluates the exact Born probability of a
given bitstring through the same code path as sampling. We scored all 1,726 hardware shots published with
#151 on an r7a.4xlarge: 6.9 s per shot, 3 h 19 min in all, 16.8 GB peak (`results/ibmxeb/`). The
run first reproduced two earlier cloud samples' log₂ Pr to 5·10⁻⁷ as a gate.

| strings | N | linear XEB |
|---|---|---|
| IBM hardware shots | 1,726 | **+0.466 ± 0.031** (bootstrap 95% [+0.41, +0.53]) |
| uniform random strings (control) | 64 | +0.16 ± 0.12 |
| IBM shots with the bit order reversed (control) | 64 | −0.15 ± 0.10 |

Our log₂ Pr agrees with #260's independently written FP64 code on every one of the 1,726 strings, to
9.4·10⁻⁶ at worst (median 2.8·10⁻⁷), and the XEB of +0.4660 is the same to four digits. Two programs
written independently assign the same probability to each of IBM's shots.

Normalising by the ideal value 2⁷⁰Σp² − 1 = 1.020 ± 0.047, measured on our own exact samples, gives a
fidelity estimate F_XEB = 0.457 ± 0.037. A second estimator uses the whole distribution of z = 2⁷⁰·p(x)
over the shots rather than its mean. If each shot is drawn from the ideal distribution with probability F
and uniformly otherwise, the density of z over shots is (Fz + 1 − F)e^{−z}, and the maximum-likelihood
fit gives F_MLE = 0.48 [0.43, 0.54] (profile likelihood, 95%). The mixture fits the observed shot
distribution (Kolmogorov–Smirnov distance 0.018, 5% threshold 0.033; `ibm_fidelity_mle.py`,
which recovers F = 0, 0.4 and 1 on synthetic data).

The instance README (tracker repository, `random_graph_sampling/README.md`) states that "XEB requires
classical simulations for the outcome probabilities, hence can only be extrapolated for the full-depth
circuit". Here it is computed exactly for every published shot. The shots are post-selected by spacetime
Pauli checks on 27 ancillas (IBM Boston, Heron R3), and the README reports a post-selected state fidelity
of 0.37 by direct fidelity estimation with 80 random stabilizers, "bounded above 0.092 with 95%
confidence". Our two estimates, 0.457 ± 0.037 and 0.48 [0.43, 0.54], are consistent with that bound and
higher than the point estimate, with a much narrower interval. The quantities differ (the full circuit
including the non-Clifford measurement layer, versus the stabilizers of the graph state; the README argues
they agree to first order under Pauli twirling), and linear XEB is known to over-estimate fidelity under
non-depolarising noise, so we report the gap without interpreting it.

## 4. Verification

1. **Small instances.** On circuits of the same form with up to 9 qubits, the distribution implied by
   the sampler matches the statevector on every bitstring to 2·10⁻¹⁵. The C++ sampler's 20,000-shot
   histogram has the same total-variation distance from the exact distribution as an ideal sampler's.
   Every test has a negative control that fails as intended.
2. **Full instance.** A sliced amplitude engine recomputes log₂ Pr of sampled bitstrings; it is a
   different algorithm from the sampler.

   | sample | sampler log₂ Pr | independent engine log₂ Pr |
   |---|---|---|
   | 1 | −69.5085302 | −69.5085294 |
   | 2 | −68.8761651 | −68.8761651 |
   | 3 | −68.3318209 | −68.3318207 |
   | 4 | −69.0978566 | −69.0978571 |
   | 5 | −70.2946875 | −70.2946876 |

   Control: sample 1 with 10 bits flipped gives −72.8808316. That is 3.37 lower in log₂ (10× less likely),
   against agreement within 1e-6 on the real samples.
3. **Reproducibility.** The production run first reproduced, bit for bit, three samples generated on a
   different machine, compiler and thread count.
4. **Engine-independent check.** While the cut rank of the measured prefix equals its length, the reduced
   state is maximally mixed, so the prefix marginal is exactly uniform; for the production order that holds
   for the first 30 steps. Over the 1,012 samples the first 30 step-bits show pair parities (435 tests) and
   triple parities (4,060 tests) at chance (4.8% and 4.4% beyond 2σ, 4.6% expected; largest |z| 3.2 and
   3.8). Of the 30 single-bit means, 4 fall beyond 2σ (p ≈ 0.05), the largest at 2.4σ. A control in which
   one bit is replaced by a copy of another fails the test (`check_uniform_prefix.py`,
   `results/uniform_prefix_check.txt`). This uses only the rank profile and the sample file.

## 5. Hardness frontier

We generated circuits with the #151 construction (linear chain, alternating CZ brickwork with random √X or
S√X layers, CZ depth n) for n = 30–120, with one instance per n, and ran the same order search
(`frontier.py`, output `results/frontier.jsonl`). At n = 70 the generated instance
matches the real one: 1,192 vs 1,150 edges, natural-order width 33 vs 32, best width 31 vs 31 at the
same search budget (a longer search on the real instance reaches 30).

| n | best width found | memory (2 complex64 tables) | table updates per exact sample |
|---|---|---|---|
| 30 | 12 | < 1 MB | 2^16.7 |
| 50 | 22 | 67 MB | 2^26.4 |
| 60 | 26 | 1.1 GB | 2^30.8 |
| **70** | **31** | **34 GB** | **2^36.0** |
| 80 | 36 | 1.1 TB | 2^41.2 |
| 90 | 42 | 70 TB | 2^46.3 |
| 100 | 47 | 2.3 PB | 2^51.4 |
| 110 | 51 | 36 PB | 2^55.9 |
| 120 | 56 | 1.2 EB | 2^61.2 |

The width fits 0.492·n − 2.9; at every n it is between n/2 − 4 and n/2 − 3. So memory doubles with
every two qubits. At depth n, this method moves from one server (n = 70) to one large-memory node
(n ≈ 80), to a leadership-class supercomputer (n ≈ 90), and past any existing machine (n ≳ 100).

Caveats. The widths are upper bounds from heuristic orders; better orders or tree decompositions would
lower them. Slicing trades memory for time at a measured cost of about 1.5× per sliced vertex on this
instance, not 2×. The table has one instance per n. For a fixed classical budget, then, the safe
direction for the experiment is more qubits at full depth. A width of ~n/2 is what random graphs
give, and the 70-qubit instance sits about 30 qubits below the point where this attack stops being
practical.

## 6. Related work

- Rank-width simulability of graph-state measurements: Van den Nest, Dür, Vidal, Briegel, PRA 75, 012337 (2007).
- ZX-calculus rank-width contraction: Kuyanov & Kissinger, arXiv:2603.06764.
- Exact amplitudes of IBM's doped circuits by T-injection and rank-width: Renaud, `dcs-rankwidth`; tracker #249.
- Capped rejection sampling for this instance on GPUs: tracker #260.
- Exact tensor-network amplitudes for the doped 70-qubit instance: Manabe, Gu, Pan, arXiv:2608.13110.

## Appendix: artifacts

- `results/samples_nq70.txt`, `run.log`, `time_nq70.txt`, `inputs.sha256`, `check55_result.txt`: the
  production sampling run and its reproduction check.
- `results/xcheck70_independent.txt`, `xcheck70d.txt`, `xcheck70e_laptop.txt`: independent
  recomputation of five cloud samples, the 10-bit control, and one laptop sample.
- `results/laptop_samples_w30.txt`: the laptop run (3 samples, timing).
- `results/ibmxeb/`: exact log₂ Pr of the 1,726 hardware shots and the controls, timing, gate,
  input hashes, and `analysis.txt` (XEB, comparison with #260, fidelity estimators).
- `results/uniform_prefix_check.txt`: the engine-independent uniformity check.
- `results/frontier.jsonl`: the frontier search output.
- `README.md`: test table and commands.
