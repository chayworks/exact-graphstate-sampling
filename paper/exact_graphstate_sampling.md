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
independent amplitude algorithm.

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

The hardware run reports a fidelity of 0.37 from direct fidelity estimation of the underlying graph state.

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

**Cost for this instance.** Our best order has maximum cut rank 31 (#260 reports 30). Each sample then
costs 4·Σ_k 2^{r_k} = 6.0·10¹⁰ table updates, with two complex64 tables of 2³¹ entries.

## 3. Results

1,012 samples in 13,202 s on an AWS r7a.4xlarge (16 cores, AMD EPYC 9R14), peak RSS 33.6 GB.
- All samples are distinct.
- Linear XEB against their own exact probabilities is 1.020 ± 0.047, i.e. 2⁷⁰ Σ p² ≈ 2.02,
  consistent with Porter–Thomas-like output.

## 4. Verification

1. **Small instances.** On circuits of the same form with up to 9 qubits, the distribution implied by
   the sampler matches the statevector on every bitstring to 2·10⁻¹⁵. The C++ sampler's 20,000-shot
   histogram has the same total-variation distance from the exact distribution as an ideal sampler's.
   Every test has a negative control that fails as intended.
2. **Full instance.** A sliced amplitude engine recomputes log₂ Pr of sampled bitstrings; it is a
   different algorithm from the sampler. Results {table: samples 1–5, control}.
3. **Reproducibility.** The production run first reproduced, bit for bit, three samples generated on a
   different machine, compiler and thread count.

## 5. Hardness frontier

We generated circuits with the #151 construction (linear chain, alternating CZ brickwork with random √X or
S√X layers, CZ depth n) for n = 30–120, with one instance per n, and ran the same order search
(`graphsamp/frontier.py`, output `graphsamp/results/frontier.jsonl`). At n = 70 the generated instance
matches the real one: 1,192 vs 1,150 edges, natural-order width 33 vs 32, best width 31 vs 31.

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

- `graphsamp/results/`: samples, run log, timing, input hashes, reproduction check.
- `graphsamp/results/xcheck70_independent.txt`: independent recomputation of samples 1 and 2.
- `graphsamp/README.md`: test table and commands.
