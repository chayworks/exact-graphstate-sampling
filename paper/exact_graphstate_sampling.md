# When are graph-state advantage experiments classically easy? Width ceilings, gate-cut bounds and sampler certificates, with IBM's 70-qubit instances as case studies

**Chaitanya Talapalli** — independent researcher — ORCID 0009-0000-8811-3349 — chayt128@gmail.com

Code, samples and raw outputs: https://github.com/chayworks/exact-graphstate-sampling (Apache-2.0).

*Draft of 6 November 2026. Every number is read from a committed artifact listed in the appendix, and every
verification step has a negative control that fails as intended.*

## Abstract

Graph-state sampling, a Clifford circuit followed by a product measurement, has been proposed as a verifiable
quantum-advantage task. For the cut-rank (linear rank-width) dynamic program that samples such states exactly we give
two general bounds on its width. Every graph on n qubits has width at most ⌊n/2⌋. A one-dimensional CZ brickwork of
depth d has width at most ⌈d/2⌉, whatever n, so on a chain the lever is min(n, d), not the qubit count. For IBM's
70-qubit instance the output is exactly uniform on every set of at most 14 qubits, because the graph state is a
[[70, 0, 15]] code, which limits what any low-order statistic of the hardware shots can show. For capped hybrid
samplers we give a total-variation bound from the exact second moment of the prefix marginal (1.1·10⁻³ at 70 qubits)
and a statistical certificate from the sampler's own proposals that assumes nothing about the output distribution.

As case studies, IBM's 70-qubit instance is sampled exactly on a laptop (11 s per sample, 8.6 GB), and its 1,726
hardware shots are scored exactly (XEB 0.466 ± 0.031). Instances at 76 and 80 qubits and IBM's two doped instances
run on one node, cross-checked against independent tensor-network and rank-width computations. Several case-study
results were obtained concurrently or earlier by others, credited in §9.

## 1. Introduction

In 2026 IBM and the University of Chicago posted a 70-qubit, CZ-depth-70 random-graph-state sampling experiment as a
candidate for verifiable quantum advantage (Quantum Advantage Tracker #151), with a post-selected state fidelity of
0.37 ± 0.02 (95% lower bound 0.092) and extrapolated classical costs of 10²⁵ s for matrix product states and 10⁷ s per
sample for extended stabilizer simulation; a doped sibling with 468 T gates (#228) and a 64-qubit ring instance with
314 T gates (#245) followed. Several groups have since attacked these instances on the tracker: capped and exact
sampling of the 70-qubit instance on GPUs (#260, #261), exact tensor-network amplitudes of the doped instance (#242,
[9]), low-precision samples of it (#258), and rank-width probabilities and samples of the ring instance (#249, #263).
This note is about why these instances are easy for one family of exact methods and what would make them hard.

**General results.**

1. *A ceiling for every graph (§5).* Any n-qubit graph state measured in any product basis is exactly samplable with
   2^{⌊n/2⌋} table entries along any linear order: 275 GB at 70 qubits, 8.8 TB at 80, 281 TB at 90. The bound is for
   linear orders; tree decompositions reach rank-width about n/3 on random graphs [14] but cost about 4^{width} per
   step with known algorithms [5, 15, 16], which is more than 2^{n/2} at that width.
2. *A gate-cut bound (§5).* For a one-dimensional brickwork of CZ depth d the spatial order has width at most ⌈d/2⌉
   for every n, measured equal to it at n = 200, d = 40 (width 20 against a ceiling of 100). With each T gate written
   as a projected ancilla the bound becomes ⌈d/2⌉ + ⌊g/2⌋, with g the largest number of vertices of one qubit's group.
   On a chain, adding qubits at fixed depth does not raise the width.
3. *Where the signal is (§6).* The 70-qubit graph state is a [[70, 0, 15]] stabilizer code, distance proven by
   exhaustive enumeration, so every marginal of the ideal output on at most 14 qubits is exactly uniform. Pauli-type
   noise cannot bias those marginals either, so they test readout asymmetry and non-Pauli channels only; IBM's shots
   pass all 974,120 parity tests up to weight four.
4. *Certificates for capped samplers (§7).* When a prefix is proposed uniformly and accepted against its exact
   marginal z with cap c, the output law is proportional to min(z, c) times the exact suffix, and its total-variation
   distance to the Born law is at most E_uniform[(z − c)⁺]. We bound this by the exact second moment of z, enumerated
   over the 2^{30} prefix stabilizers at 70 qubits (1.1·10⁻³), and, at any size, from the proposals themselves with an
   empirical Bernstein bound (fidelity ≥ 0.887–0.940 at 95% on the four instances run).

**Case studies.** IBM's 70-qubit instance is sampled exactly from the Born distribution, 13.0 s per sample on a
16-core node and 11 s on a laptop with an 8.6 GB table, and the exact probability of every hardware shot gives
linear XEB 0.466 ± 0.031, which the instance README says can only be extrapolated (#261 scores the shots at the same
value with a different method). Exact samples at 76 qubits and hybrid samples at 80 run on one node, the latter by a
translation of the dual table that makes slicing possible when every vertex touches the sampled suffix (§3.4). IBM's
doped instances become graph states in a product basis once each T gate is a projected ancilla: the 70-qubit one at
width 34, agreeing with #242's tensor-network amplitudes shot by shot up to a common factor of 1 − 1.8·10⁻⁴, and the
ring at width 34 after slicing, agreeing with #249 to 4·10⁻⁶ in log₂.

The method joins two standard ingredients by one identity: rank-width simulability of graph-state measurements [1]
and perfect sampling in the style of Ferris and Vidal [4], joined by the Schmidt decomposition of a graph state across
a cut [2]. The post-measurement states Z^v|G⟩ are orthonormal, so every prefix marginal is a plain sum of squares,
Pr(b₁…b_t) = 2^{−t} Σ|P_t|². We verify the sampler against dense statevectors on small instances, against an
independent amplitude algorithm on the full instance, against a closed-form stabilizer-group formula for its per-step
conditionals (3·10⁻⁸ on all 1,726 shots), and across three dynamic-program paths.

The note is organised as follows: §1.1 states the task; §2 the method; §3 the case studies on the 70-, 76- and
80-qubit instances and on IBM's shots; §4 the verification; §5 the ceiling, the gate-cut bound and the frontier; §6
the stabilizer structure of the output and the simulation-free test of the shots; §7 the doped instances and the
sampler certificates; §8 what the work does and does not show; §9 related work; the appendix lists the artifact
behind every number.

### 1.1 The task

Circuit `nq70_depth70_checks27_basis_fstate`, Quantum Advantage Tracker issue #151:
1. Hadamards on 70 qubits in a line.
2. Brickwork of nearest-neighbour CZ gates to CZ depth 70, with a random √X or S√X on every qubit between
   layers.
3. On every qubit, T†·H·Rz(θ)·H with θ = arccos(1/√3), a rotation into the "face-state" basis that maximises
   stabilizer extent.
4. Measurement in Z.

The hardware run (IBM Boston, Heron R3, with 27 ancillas for spacetime Pauli checks) reports a post-selected
fidelity of 0.37 ± 0.02 from direct fidelity estimation of the underlying graph state, with a 95% lower bound of
0.092. 1,726 post-selected shots are published with the instance. The construction was chosen so that
stabilizer-extent methods cost 10⁷ s per sample and matrix-product-state methods 10²⁵ s, according to the
instance README.

## 2. Method

**Graph form.** The Clifford part maps to a stabilizer tableau, and then, by Gaussian elimination and
local Cliffords, to |S⟩ = ⊗U_j† |G_A⟩. Then

  ⟨b|ψ⟩ = 2^{-n/2} Σ_x Π_j M_j[b_j, x_j] (−1)^{Σ_{i<j} A_ij x_i x_j},  where M_j = W_j U_j†.

**Cut-rank dynamic program.** Fix a vertex order. After the prefix P is processed, the future depends
only on s = A[F,P]·x_P ∈ GF(2)^{|F|}. That vector lives in a space of dimension r = rank A[F,P], so the
prefix is a table P[c] of 2^r complex numbers in a chosen basis of that space. Each step is a surjective
GF(2)-linear map (c, x_v) ↦ c′. We evaluate it in gather form: a `pdep` scatter plus a kernel of
dimension at most 2. A generic matrix product state of the same cut needs bond dimension 2^r and therefore
cost polynomial in 2^r; here the cost is linear in 2^r.

**Exact sampling.** With prefix P measured with outcomes b_P, the post-measurement state is

  2^{-|P|/2} Σ_c P[c] Z^{B c} |G_F⟩.

For v ≠ 0, Z^v has no X part, so it anticommutes with some graph stabilizer and ⟨G|Z^v|G⟩ = 0. The
states Z^v|G_F⟩, one for each v, are therefore orthonormal, and Pr(b_P) = 2^{-|P|} Σ_c |P[c]|². Measuring qubit by qubit
in the DP order:
1. compute both candidate extensions' norms (one read-only pass);
2. draw the outcome;
3. store the chosen extension (a second pass).

This is the graph-state analogue of perfect MPS sampling (Ferris & Vidal 2012). The output of this chain-rule
sampler, used for the 70- and 76-qubit instances, is the exact Born distribution up to floating-point error, and
unlike rejection sampling its cost does not depend on the shape of the output distribution. Where the chain rule
is unavailable (projected ancillas, §7) or too costly (a sliced table, §3.4) we use a hybrid sampler that
re-introduces a capped acceptance step on a prefix marginal, and we label its output accordingly rather than as exact.
Its ingredients are not new: Manabe, Gu and Pan [9, §§4.1–4.2] describe batched sampling (a uniform prefix with the
suffix drawn exactly from one batch of amplitudes) and frugal rejection sampling, analyse both under a Porter–Thomas
assumption, and note that combining them reduces the batched sampler's error. The hybrid sampler is that combination;
what we add is a bound on its excluded mass that does not assume Porter–Thomas statistics (§7).

**In-place layout.** A second program keeps a single table and updates it in place, so the memory is
2^{width} complex64 entries rather than two tables of 2^{width+1}: 8.6 GB for the width-30 order of this
instance, which is why a laptop can run it.

**Cost for this instance.** A greedy-plus-local-search order (`order2.py`) reaches maximum cut rank 30, the
value #260 also reports; the natural order gives 32. The production sampling run used an earlier width-31
order and cost 4·Σ_k 2^{r_k} = 6.0·10¹⁰ table updates per sample with two complex64 tables of 2³¹ entries;
the later runs (laptop, hardware-shot scoring, chain mode) use the width-30 order.

## 3. Results

### 3.1 Sampling IBM's instance

1,012 samples in 13,202 s on an AWS r7a.4xlarge (16 cores, AMD EPYC 9R14), peak RSS 34.4 GB, about $4.50
of on-demand time.
- All samples are distinct.
- Linear XEB against their own exact probabilities is 1.020 ± 0.047, i.e. 2⁷⁰ Σ p² ≈ 2.02, the Porter–Thomas
  value.
- The distribution of z = 2⁷⁰·p over the samples follows the Porter–Thomas size-biased law z·e^{−z}: KS distance
  0.016 (5% level 0.043), mean z² 6.3 ± 0.3 against 6, mean log₂ Pr −69.398 ± 0.038 against the predicted −69.390
  (`porter_thomas_test.py`, `results/porter_thomas_nq70.txt`). So the ideal distribution over strings is Exp(1) as
  far as 1,012 samples can tell, which the fidelity estimators of §3.2 assume.

**Laptop.** With the in-place program, the arena is a single complex64 table of 2^30 entries (8.6 GB) for a
width-30 order. On a laptop (AMD Ryzen 9 9955HX, 16 threads, container memory capped at 10 GB), 3 exact
samples took 55 s including allocation, i.e. about 18 s per exact sample; the 1,012-sample rerun of §4 item 7 on
the same laptop and plan, using all 32 hardware threads, took 190 min, i.e. **11.3 s per exact sample**
(`results/time_nq70_reseeded.txt`). The first run's timing is
(`results/laptop_samples_w30.txt`). The in-place sampler was validated on every bitstring of small circuits.
On the full instance it reproduces the cloud samples' log₂ Pr to 3·10⁻⁷, and the independent amplitude
engine confirms the first laptop sample: −68.5019327 against −68.5019326 reported
(`results/xcheck70e_laptop.txt`).

### 3.2 IBM's hardware shots, scored exactly

The sampler's probability mode evaluates the exact Born probability of a given bitstring through the same
code path as sampling. We scored all 1,726 hardware shots published with #151 on an r7a.4xlarge: 6.9 s per
shot, 3 h 19 min in all, 17.2 GB peak (`results/ibmxeb/`). The run first reproduced two earlier cloud
samples' log₂ Pr to 5·10⁻⁷ as a gate.

| strings | N | linear XEB |
|---|---|---|
| IBM hardware shots | 1,726 | **+0.466 ± 0.031** (bootstrap 95% [+0.41, +0.53]) |
| uniform random strings (control) | 64 | +0.16 ± 0.12 |
| IBM shots with the bit order reversed (control) | 64 | −0.15 ± 0.10 |

Our log₂ Pr agrees with #260's independently written FP64 code on every one of the 1,726 strings, to
9.4·10⁻⁶ at worst (median 2.8·10⁻⁷), and the XEB of +0.4660 is the same to four digits. Two programs
written independently assign the same probability to each of IBM's shots. The uniform control is small (64 strings,
1.3σ from zero); it checks the sign and scale of the statistic rather than its tail, and a control as large as the
shot set would cost two hours of the same node.

Normalising by the ideal value 2⁷⁰Σp² − 1 = 1.020 ± 0.047, measured on our own exact samples, gives a
fidelity estimate F_XEB = 0.457 ± 0.037. A second estimator uses the whole distribution of z = 2⁷⁰·p(x)
over the shots rather than its mean. If each shot is drawn from the ideal distribution with probability F
and uniformly otherwise, the density of z over shots is (Fz + 1 − F)e^{−z}, and the maximum-likelihood
fit gives F_MLE = 0.48 [0.43, 0.54] (profile likelihood, 95%). This model takes the ideal z over strings to be
Exp(1), which §3.1 supports for this instance; the synthetic F = 1 case below (0.98 [0.95, 1.00]) shows the size of
the residual model bias. The mixture fits the observed shot
distribution (Kolmogorov–Smirnov distance 0.018, 5% threshold 0.033; `ibm_fidelity_mle.py`). On synthetic
shots of known fidelity the two estimators return 0.00 ± 0.02 and 0.00 [0.00, 0.05] for F = 0, 0.42 ± 0.04 and
0.41 [0.35, 0.46] for F = 0.4, and 0.94 ± 0.05 and 0.98 [0.95, 1.00] for F = 1 (`fidelity_mle_synth.py`,
`results/fidelity_mle_synthetic.txt`).

![Exact probabilities of IBM's shots](fig/fig2_shots.png)

*Figure 1. Distribution of z = 2⁷⁰·p(x) over IBM's 1,726 post-selected shots (exact p), our 1,012 exact samples, and the fitted clean/uniform mixture (`paper/make_figures.py`).*

The instance README (tracker repository, `random_graph_sampling/README.md`) states that "XEB requires
classical simulations for the outcome probabilities, hence can only be extrapolated for the full-depth
circuit". Here it is computed exactly for every published shot. The shots are post-selected by spacetime
Pauli checks on 27 ancillas, and the README reports a post-selected state fidelity of 0.37 by direct
fidelity estimation with 80 random stabilizers, 0.37 ± 0.02 and "bounded above 0.092 with 95% confidence". Our
two estimates, 0.457 ± 0.037 and 0.48 [0.43, 0.54], are consistent with the lower bound and about 2σ above the point
estimate if ±0.02 is its full uncertainty; the README's own 95% bound suggests it is not. The quantities also differ
(the full circuit including the non-Clifford measurement layer, versus the stabilizers of the graph state; the README
argues they agree to first order under Pauli twirling), and linear XEB can over-estimate fidelity under
non-depolarising noise, so we do not interpret the difference. What the exact probabilities add is an interval about
four times narrower than the one IBM's lower bound implies, from a per-shot benchmark at full depth.

### 3.3 A 76-qubit instance, sampled

To test the frontier with a measurement rather than a projection, we generated a 76-qubit, CZ-depth-76
instance of the same construction (`frontier_instance.py`, seed 1076, 1,465 edges), found an order of width
33 (natural order 38), and ran the in-place sampler on the same 16-core node class (r7a.4xlarge, 128 GiB):
**19 exact samples in 3,604 s, 190 s per sample, 69 GB peak** (`results/nq76/`). The node first reproduced
the two 70-qubit reference probabilities as a gate. All 19 samples are distinct; their self-consistency XEB
is 0.83 ± 0.22. The cross-check scheduled in that run produced nothing (its complex128 tables, 137 GB, were
killed by the kernel; `results/nq76/run.log`), so it was repeated in a second run. The independent sliced
amplitude engine, rebuilt with complex64 tables and itself gated on the two known 70-qubit values (−69.5085289 and −68.8761655 against −69.5085302 and −68.8761651), recomputes
the first two samples' log₂ Pr as −74.1389442 and −75.3402128 against the sampler's −74.1389441 and
−75.3402127 (`results/nq76/analysis.txt`).

### 3.4 An 80-qubit instance, sampled on one node

The 80-qubit, depth-80 instance (seed 1080, 1,573 edges) has an order of width 35, a 275 GB in-place table that no
single node in our quota holds. Slicing one vertex (fixing its variable and summing the two half-problems) lowers the
width to 34, but every vertex of this dense graph (degrees 30–50) has neighbours among the 20 sampled suffix qubits,
which the slicing of §7 forbids: the slice's phase on a suffix neighbour cannot be folded into that neighbour's weight,
because the suffix steps run once on the summed table. The way out is a property of the table itself. The in-place
table is the dual representation: its entry at logical coordinate u is the summed weight of all prefix assignments
whose linear functional on the suffix is φ(u) = Σ_i u_i B_i, where B is the cut basis the planner stores, and the
partial amplitude is F(x_S) = Σ_u T(u)(−1)^{⟨φ(u),x_S⟩}. A phase (−1)^{⟨s,x_S⟩} on the suffix is therefore a
translation of the table, T(u) → T(u + a) with s = Σ a_i B_i (the shift theorem of the Walsh–Hadamard transform:
multiplying by a character translates the transform). The translation exists whenever s lies in the span of B, i.e.
whenever the slice vertex's row on the suffix is a combination of the remaining prefix rows; at the summation cut of
a hybrid run (rank ≤ 20, spanned by 59 prefix rows) that held for every slice we planned, and the planner refuses a
slice for which it does not (no fallback was needed or built). The planner writes the translation for each slice and
the engine applies it when it sums the dumps (`sliced_plans.py`, `gsamp2 resume --twist`). We checked the
representation entry by entry against a brute-force sum on a toy T-gadget circuit (`test_dual_table.py`, 2.4·10⁻⁸ at worst, with the
unsigned sum off by 1.41 as a control), the twisted two-slice sampler against the unsliced one on that toy (9.1·10⁻⁸
in log₂ Pr, with the untwisted sum off by 2.54 as a control) and at full scale on the 70-qubit instance (three
proposals, identical sampled bits, 1.8·10⁻⁶ in log₂ Pr and 3.2·10⁻⁷ on the prefix marginal; `test_twist.py`); the
outputs of both scripts are in `results/tests/`. The same checks ran as gates on the cloud node before the 80-qubit sampling
(`results/nq80/run.log`).

With vertex 34 sliced (width 34, 137 GB per slice) the sampler is the hybrid one of §7, not the chain rule: a
chain-rule sample would need both slice passes at every one of the 80 steps, since the in-place table cannot be
branched, about 80 times the cost of a proposal. Each proposal draws the first K = 60 bits uniformly; one pass per
slice gives their exact marginal z (in units of the uniform distribution) together with exact draws of the last 20
bits, and the proposal is accepted with probability z/c, c = 1.03. One r7a.8xlarge (32 cores, 256 GiB) produced
**79 samples of the 80-qubit instance at 246 s per proposal** (two passes of 2^{34} per proposal, peak resident set
128 GiB, 5.5 h of sampling; `results/nq80/`): 79 of 80 proposals accepted, z between 0.9970 and 1.0018 (the spread
expected when z averages 2^{20} Porter–Thomas probabilities is 2^{−10}), none above the cap; mean log₂ Pr
−79.29 ± 0.12 against the Porter–Thomas value −79.39 and −80 for the uniform distribution; XEB against their own
exact probabilities 1.02 ± 0.14, KS distance to z·e^{−z} 0.095 (5% level 0.153; `results/nq80/porter_thomas_nq80.txt`);
all 79 strings distinct. Because a proposal with z > c would always be accepted, the output law is proportional to
min(z, c) times the exact suffix; it equals the Born law exactly if no prefix has z > c, and in general its total-variation
distance to it is at most E_uniform[(z − c)⁺] (§7). The 80 proposals bound that by the certificate of §7 at TV ≤ 0.113,
classical fidelity ≥ 0.887 at 95% confidence (`results/certificates_nq80_fstate.txt`); a rigorous second-moment bound
like the one at 70 qubits would need the 2^{40} prefix stabilizers, so we do not call these exact samples. The first two were re-scored with a second slicing of the same
engine (vertices 10 and 40, four slice passes of 472 s), which agrees with the sampler to 5.2·10⁻⁷ and 4.4·10⁻⁷ in
log₂ Pr; no independent engine has evaluated an 80-qubit probability (the independent engine's checks are at 70 and
76 qubits). Before sampling, the node passed two gates: the two known 70-qubit values to
2.6·10⁻⁷, and the twisted two-slice run against the unsliced hybrid sampler on three 70-qubit proposals (identical
bits, 8.2·10⁻⁷ in log₂ Pr). A scan of 16 other circuit seeds at n = 80 found widths 35–36 only
(`results/frontier_instances/fi80_cseed_scan.log`), so the slicing, not a lucky instance, is what put 80 qubits on one node.

## 4. Verification

1. **Small instances.** On circuits of the same form with up to 9 qubits, the distribution implied by
   the sampler matches the statevector on every bitstring to 2·10⁻¹⁵. The C++ sampler's 20,000-shot
   histogram has the same total-variation distance from the exact distribution as an ideal sampler's.
   Every test has a negative control that fails as intended (`test_*.py`; outputs in `results/tests/`).
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
3. **Reproducibility.** The production run first reproduced, bit for bit, three samples of a 55-qubit
   sub-instance (width 26) generated on a different machine, compiler and thread count
   (`results/check55_result.txt`).
4. **Engine-independent check.** While the cut rank of the measured prefix equals its length, the reduced
   state is maximally mixed, so the prefix marginal is exactly uniform; for the production order that holds
   for the first 30 steps. Over the 1,012 samples the first 30 step-bits show pair parities (435 tests) and
   triple parities (4,060 tests) at chance (4.8% and 4.4% beyond 2σ, 4.6% expected; largest |z| 3.2 and
   3.8). Of the 30 single-bit means, 4 fall beyond 2σ (p ≈ 0.05), the largest at 2.4σ. A control in which
   one bit is replaced by a copy of another fails the test (`check_uniform_prefix.py`,
   `results/uniform_prefix_check.txt`). This uses only the rank profile and the sample file.
5. **Stabilizer-only conditionals.** Once the cut rank falls below the prefix length, the prefix marginal has
   a closed form from the stabilizer group alone: Pr(b_P) = 2^{-t} Σ_S sign(S) Π_{j∈supp S} (−1)^{b_j} n_j[S_j],
   summed over the 2^{t−rank} stabilizers supported on the prefix P, with n_j the Bloch vector of the measured
   operator. We enumerate that sum exactly while the deficit is at most 16 (steps 1–42 of the width-30
   order) and compare the resulting conditionals Pr(b_k | b_1…b_{k−1}) with the sampler's chain output for
   all 1,726 hardware shots: the two agree to 3.3·10⁻⁸ (half a float32 ulp at ½, the sampler's storage
   precision), while the conditionals
   themselves depart from ½ by up to 8·10⁻⁵ at those steps (`stabprefix.py`,
   `results/chain/stabprefix_ibm_fwd.txt`). On small circuits the same formula reproduces the
   statevector marginals to 6·10⁻¹⁵ on 162 prefixes, and a flipped generator sign fails it. This is a
   check of the dynamic program's fine structure by a method that never builds the DP table.

   The same formula explains where the information in a string lives: a stabilizer of weight w contributes
   at most max|n|^w, so the conditionals stay within 10⁻⁴ of ½ until the prefix supports many low-weight
   stabilizers, which in this order happens only in the last ~10 steps; the per-step spread of the chain
   conditionals grows by about 1.4× per step from 7·10⁻⁴ at step 50 to 0.56 at step 70, and the last 8 steps
   carry 2.83 of the 2.85 variance of log₂ Pr across the hardware shots
   (`results/chain/analysis_profile.txt`). For our 1,012 exact samples the picture is the same: the
   spread of the conditionals first exceeds 10⁻³ at step 51 (52 in the reversed order), and the last 4 and
   8 steps account for 90% and 100% of the variance of log₂ Pr (`results/chain/analysis_profile.txt`).
6. **Three dynamic-program paths agree.** Scoring every string along the width-30 order, along its reverse
   (the same width, a different table layout and step sequence), and, for our samples, along the width-31
   production order gives the same log₂ Pr: the forward and reversed orders agree to 1.4·10⁻⁵ at worst
   (median 4·10⁻⁷) on the 1,726 hardware shots and to 3.5·10⁻⁶ on the 1,012 samples, and the width-31 run
   agrees with the width-30 order to 4.3·10⁻⁶ (`results/chain/analysis_crossorder.txt`). The
   stabilizer-only check of item 5 also holds along the reversed order and on the ideal samples (3.3·10⁻⁸,
   `results/chain/stabprefix_rev_ideal.txt`).
7. **Simulation-free test, and the sampler's random stream.** By §6 every marginal of the ideal distribution
   on at most 14 qubits is exactly uniform, so all single-bit means, pair parities and triple parities are
   exactly zero (for ±1 bits), and their sums of squared z-scores are χ² with as many degrees of freedom as
   tests. Our 1,012 samples: pairs p = 0.76, triples p = 0.97, quadruples p = 0.12, single bits χ² = 100 on 70
   (p = 0.010); the same run scores truly uniform random bits as a calibration and a control with one bit
   biased by 15% fails it (p = 1.8·10⁻⁷) (`hw_marginals.py`, `results/hw_marginals.txt`).
   The single-bit p = 0.010 is unremarkable among the dozens of tests run, but we checked it: the generator was
   std::mt19937_64 seeded with consecutive integers 10⁶ + i, a known weak initialisation [18], and
   a direct test of the drawing step against the exact conditionals passes where they carry the information (χ² 23.1
   on 20, p = 0.29; a control drawing with p^{1.1} fails at p = 7·10⁻⁵⁰; `draw_check.py`). A second run of 1,012
   samples seeded through `std::seed_seq` (11.3 s per sample on the laptop; `results/samples_nq70_reseeded.txt`)
   shows no excess: single bits p = 0.65, pairs 0.31, triples 0.64, quadruples 0.64 (`results/hw_marginals_reseeded.txt`).

![Where the information lives](fig/fig3_information.png)

*Figure 2. Spread of the exact per-step conditionals along the width-30 order for the hardware shots and the ideal samples, and the mean of 2c_k − 1 over ideal samples where it is resolved (the per-step XEB signal). The first 22 conditionals are exactly ½; the information sits in the last ~8 steps.*

**Per-step fidelity of the hardware shots.** Dividing the hardware mean of 2c_k − 1 by the ideal mean gives a
fidelity-like ratio for each informative step. Because the information sits in the last steps, only four
steps per order are informative (ideal mean above 0.05): in the forward order, qubits 9, 8, 13 and 50 give
0.45 ± 0.12, 0.54 ± 0.08, 0.52 ± 0.07 and 0.48 ± 0.04; in the reversed order, qubits 45, 68, 53 and 21 give
0.47 ± 0.10, 0.43 ± 0.07, 0.43 ± 0.06 and 0.53 ± 0.05. All eight are consistent with the global fidelity
estimates of 0.46–0.48 and none stands out, so the shots show no sign of a qubit that is much worse than
the rest among those eight. A per-qubit map of all 70 qubits would need an order ending in each group of
about four qubits, roughly 18 further width-30 passes; we did not run them.

## 5. Hardness frontier

We generated circuits with the #151 construction (linear chain, alternating CZ brickwork with random √X or
S√X layers, CZ depth n) for n = 30–120, with one instance per n, and ran the same order search
(`frontier.py`, output `results/frontier.jsonl`). At n = 70 the generated instance matches the real one:
1,192 vs 1,150 edges, natural-order width 33 vs 32, best width 31 vs 31 at the same search budget (a
longer search on the real instance reaches 30).

| n | best width found | memory (2 complex64 tables) | table updates per exact sample |
|---|---|---|---|
| 30 | 12 | < 1 MB | 2^16.7 |
| 50 | 22 | 67 MB | 2^26.4 |
| 60 | 26 | 1.1 GB | 2^30.8 |
| **70** | **31** | **34 GB** | **2^36.0** |
| **76 (sampled)** | **33** | **69 GB in-place** | **2^38.2** |
| **80 (sampled, hybrid)** | **35, sliced to 34** | **275 GB in-place; 137 GB per slice** | **2^40.3** |
| 90 | 42 | 70 TB | 2^46.3 |
| 100 | 47 | 2.3 PB | 2^51.4 |
| 110 | 51 | 36 PB | 2^55.9 |
| 120 | 56 | 1.2 EB | 2^61.2 |

![Hardness frontier](fig/fig1_frontier.png)

*Figure 3. Best width found versus n for the construction, the ⌊n/2⌋ ceiling that holds for every graph, the three instances sampled (70 and 76 exactly, 80 with the hybrid sampler), and the memory of the in-place table. Shaded bands mark what a laptop, one large cloud node (up to 32 TB), a supercomputer, and no existing machine can hold; the largest single shared-memory server sold, 64 TB [21], sits inside the second band. The 76- and 80-qubit order searches are in `results/frontier_instances/` (fi76_s2: width 33, log₂ Σ_k 2^{r_k} = 36.2, i.e. 2^{38.2} updates; fi80_s5: width 35, log₂ Σ_k 2^{r_k} = 38.3, i.e. 2^{40.3} updates, run at width 34 with one vertex sliced, §3.4).*

![Cut-rank profiles](fig/fig5_rankprofile.png)

*Figure 4. Cut rank r_k after step k along the orders used in this note, against the fraction of the order processed: IBM's 70-qubit instance in the production order (width 31) and the width-30 order, the 80-qubit instance in its sliced order (width 34), and the 538-vertex gadget graph of the doped instance in the sweep order (width 34). For the three curves read from in-place plans the table at step k holds 2^{r_{k−1}+1} entries, so the plan's width is one more than the largest r_k shown (29 → 30, 33 → 34, 33 → 34); the production-order curve is computed from the adjacency as the rank after each step, and its largest value, 31, is the two-table width. The dotted line is the ceiling min(k, n − k) for n = 70 (`paper/make_figures.py`, from the committed plans and `results/nq70_graph_adjacency.txt`).*

At this search budget the width is n/2 − 4 to n/2 − 3 at every n (least-squares fit 0.492·n − 2.9), and longer
searches on the 70-, 76- and 80-qubit instances reach n/2 − 5 (widths 30, 33, 35); the true linear rank-width of these
graphs is unknown. Because the construction sets d = n, the table cannot separate the two levers; the gate-cut bound
below does. At d = n, memory doubles with every two qubits, and the method moves from one 128 GB server (n = 70) to one
256 GiB node with slicing (n = 80), to the largest shared-memory server sold, 64 TB [21] (n ≈ 93 at the fitted width),
and to a leadership-class machine at n ≈ 100 (2.3 PB with two tables). Those are memory statements. The step is a GF(2)
linear shuffle of the whole table, which is bound by communication once the table is distributed, so a supercomputer
run would be slower per update than the single-node rates measured here, by a factor we have not measured.

**Proposition (a ceiling for every graph).** *Let |G⟩ be a graph state on n qubits and let each qubit be
measured in an arbitrary single-qubit basis. For every linear order of the qubits, the cut rank of every
prefix is at most ⌊n/2⌋. Consequently the sampler above draws an exact Born sample with at most 2^{⌊n/2⌋}
complex table entries and at most 4·n·2^{⌊n/2⌋} table updates, whatever the graph and whatever the
measurement bases.*

*Proof.* The cut rank of a prefix P is the GF(2) rank of the |F| × |P| adjacency block A[F,P], which is at
most min(|P|, |F|) = min(t, n − t) ≤ ⌊n/2⌋. The in-place sampler's table has 2^{r_k} entries at step k and
each step makes at most four passes over it. ∎

The proposition is elementary, but its consequence is not widely stated: at n = 70 the ceiling is 275 GB and
about 9.6·10¹² updates per sample; at the rates measured in this work (6.0·10¹⁰ updates in 13.0 s for the
two-table sampler, 2^{38.2} updates in 190 s for the in-place sampler at 76 qubits, i.e. 4.6·10⁹ and
1.7·10⁹ per second on 16 cores) that is between 35 minutes and 1.6 hours per sample. The ceiling's memory is
2.2 TB at n = 76, 8.8 TB at 80, 17.6 TB at 82, 281 TB at 90 and 9 PB at 100. Against the machines that exist: a
2 TB cloud node holds the ceiling to n = 75, a 64 TB server [21] to n = 85, and 281 TB of aggregate memory, the scale
of a leadership machine, to n = 90. So no graph-state sampling experiment with product measurements needs more than
one shared-memory server below about 86 qubits, whatever the graph. The random brickwork sits 3–5 below this ceiling
at our search budgets (Figure 3), so at d = n the choice of graph buys at most a factor of about 30 in memory.

**Proposition 2 (a gate-cut bound).** *Let a circuit on n qubits in a line consist of single-qubit Clifford gates and
CZ gates between neighbours, with at most D CZ gates on any bond, followed by a product measurement, and let |G⟩ be the
graph state it prepares up to local Cliffords. Along the spatial order 1, …, n every cut rank is at most D. If each
non-Clifford gate is replaced by a projected ancilla attached to its qubit (§7) and the order lists each qubit's
ancillas next to it, every cut rank is at most D + ⌊g/2⌋, with g the largest number of vertices of one qubit's group.*

*Proof.* The cut rank across a bipartition equals the entanglement entropy of |G⟩ across it in bits [2], local
Cliffords leave it unchanged, and each CZ gate raises the entanglement across a cut it crosses by at most one bit while
gates on one side leave it unchanged; a cut between qubits k and k + 1 is crossed only by the CZs on that bond. With
ancillas, a cut between two qubit groups is crossed only by the CZs on one bond, since each ancilla's CNOT stays inside
its qubit's group. A cut inside a group of g vertices differs from the cuts at either end of the group by moving at
most ⌊g/2⌋ vertices, and moving one vertex changes a cut rank by at most one. ∎

For the brickwork of depth d, D = ⌈d/2⌉. Measured along the spatial order (`gatecut_check.py`,
`results/gatecut_check.txt`):

| n | d | spatial-order width | ⌈d/2⌉ | ⌊n/2⌋ |
|---|---|---|---|---|
| 40 | 40 | 18 | 20 | 20 |
| 70 | 30 | 15 | 15 | 35 |
| 70 | 70 | 35 | 35 | 35 |
| 70 | 100 | 35 | 50 | 35 |
| 100 | 40 | 20 | 20 | 50 |
| 150 | 70 | 35 | 35 | 75 |
| 200 | 40 | 20 | 20 | 100 |

The bound is attained whenever d < n, so on a chain the width is min(⌈d/2⌉, ⌊n/2⌋) up to what an order search finds
below it (30 rather than 35 at n = d = 70): a 200-qubit chain of depth 40 has width 20. The
same bound explains the width of 35 that Manabe, Gu and Pan's contraction reaches on the doped instance [9] and the
spatial-sweep width of 34 in §7 (D = 35). The ring instance is not explained by it: a cut of a ring crosses two bonds,
so the bound is 74 at d = 73, while Renaud's order reaches 36 [6].

**Is the ceiling tight?** Running the same order search (four randomised greedy orders then 3,000 local-search
moves, two instances per cell; the command is recorded in `results/ceiling_tightness.cmd.txt`) on Erdős–Rényi
graphs G(n, ½) and G(n, 0.3) and on the brickwork construction gives
(`ceiling_tightness.py`, `results/ceiling_tightness.jsonl`):

| n | ceiling ⌊n/2⌋ | G(n, ½) | G(n, 0.3) | brickwork |
|---|---|---|---|---|
| 40 | 20 | 16, 16 | 16, 16 | 16, 17 |
| 50 | 25 | 21, 21 | 22, 21 | 22, 21 |
| 60 | 30 | 26, 26 | 26, 26 | 26, 26 |
| 70 | 35 | 32, 32 | 32, 31 | 31, 31 |

Random graphs of either density sit 3–4 below the ceiling and the brickwork within one of them, so for this
method the construction is as hard as a random graph. Whether the ceiling is tight is open: the gap of 3–5 is what
this heuristic finds at this budget, and it grew by one when the search on the real instances was lengthened.

**Caveats.** The widths are upper bounds from heuristic orders. Both bounds are for linear orders: the rank-width of
a random graph G(n, p) is ⌈n/3⌉ − O(1) [14], below the ⌊n/2⌋ ceiling, and tree-shaped rank decompositions can be used
for simulation [5, 15, 16], but the known dynamic programs over them cost about 4^{width} per step, 2^{2n/3} at width
n/3, against 2^{n/2} here; a cheaper combining step would change the design advice of §8. Slicing trades memory for time; we have not measured the exchange rate on one machine (the sliced and
unsliced doped runs of §7 also differ in engine). The table has one instance per n. The ceiling bounds this
family of methods, not all classical methods, and applies to a circuit only once it is a graph state followed by a
product measurement: for the doped sibling instance that form is reached through the T-gadget of §7, on a
538-vertex graph whose ceiling (2^{269}) is vacuous; what makes it tractable there is the measured width of 34, not
the ceiling.

## 6. Where the signal is: no marginal on fewer than 15 qubits carries any

In the graph frame every stabilizer is S_v = X^v Z^{Av} with support v ∪ Av, and by the formula of §4 item 5
the marginal of the output distribution on a qubit set Q is Pr(b_Q) = 2^{-|Q|} Σ_{S supported on Q} sign(S)
Π_j (−1)^{b_j} n_j[S_j]. The measured Bloch vectors have no zero component (all three components are ±1/√3
in the graph frame), so:

**Proposition.** *For the face-state basis (every Bloch component ±1/√3), the marginal on Q is exactly uniform
if and only if no nontrivial stabilizer is supported inside Q. In particular every marginal on fewer than d_min
qubits is exactly uniform, where d_min is the minimum weight of a stabilizer of the graph state, and a
stabilizer of weight w changes a marginal by at most the fraction 3^{−w/2} of its uniform value.*

*Proof.* "If" is the formula with the identity as the only term. For "only if", let U be a minimal nonempty
support among the stabilizers inside Q and H_U the subgroup of stabilizers with support inside U; by
minimality every nontrivial element of H_U has support exactly U. For any j ∈ U the map S ↦ S_j (the Pauli letter at j, modulo phase) is a
homomorphism from H_U to the one-qubit Pauli group modulo phases whose kernel consists of elements with identity at j, hence
of support smaller than U, hence trivial; so H_U embeds in the Klein group and has 1 or 3 nontrivial
elements. The Fourier coefficient of the marginal at the character Π_{j∈U}(−1)^{b_j} is the sum over exactly
those elements of sign(S) Π_j n_j[S_j] = ±3^{−|U|/2}, an odd number of terms of equal magnitude, which cannot
vanish. ∎ (For a general product basis with nonzero Bloch components "if" still holds and "only if" holds
unless three such terms cancel exactly.)

In coding terms |G⟩ is a [[n, 0, d_min]] stabilizer code and d_min is its distance; the uniformity of every marginal
on fewer than d_min qubits in every product basis is the statement that the state is (d_min − 1)-uniform, the known
equivalence between k-uniform states and pure codes of distance k + 1 (Scott 2004). The "if" direction of the
proposition is that fact; the "only if" for the face-state basis is the part specific to this measurement; and the
minimum-weight search below is a code-distance computation.

For IBM's instance **d_min = 15 exactly**, by a three-way split that is exhaustive (`minweight_exact.py`,
`results/minweight_exact.txt`): a stabilizer of weight ≤ 14 has either at most 11 X-support vertices (all
2.6·10¹² such v enumerated, lightest 15; `minweight.cpp`, `results/minweight_k11.txt`), or at most 11
Z-support vertices c = Av (all such c enumerated through a generalised inverse of the rank-68 adjacency and
its 4-element kernel coset, lightest 15; `minweight_case1.cpp`, `results/minweight_case1_k11.txt`), or at least
12 of each, in which case v and Av differ in at most 4 positions and v is recovered from (A+I)v = e for the
974,120 vectors e with |e| ≤ 4 (lightest 16; A+I has full rank, so e = 0 contributes nothing and each e has
one solution). The first enumeration alone cannot certify the minimum, which is why its output file ends
"NOT PROVEN"; the certificate is the union of the three. The lightest stabilizer has v = {0, 5, 27, 42, 48, 49, 52, 58, 61}
and support of 15 qubits. The split was validated on the 30- and 40-qubit generated instances, whose minima
(6 and 8) are known from direct enumeration: it returns the same minima, and on the 30-qubit instance, run
with t = 2 so that the minimiser lies outside the directly enumerated part, only case (2) finds the weight-6
stabilizer, so dropping that case would report 8 (at 40 qubits, t = 4, the direct part finds the minimiser and
the other two cases return 9) (`minweight_split_test.py`, `results/minweight_split_test.txt`). So **every marginal on at most 14 qubits is exactly uniform**, and even the first
non-uniform marginals deviate by at most about 3^{−7.5} ≈ 3·10⁻⁴ per stabilizer. This is the stabilizer-side
statement of Figure 2: no test restricted to 14 or fewer qubits, however many shots, can tell the ideal output
from coin flips; the signal lives in correlations of weight ≥ 15, and verifying such an experiment is global by
nature.

The threshold grows with the size of the experiment. For generated instances of the same construction
(`frontier.py` seeds, exhaustive enumeration with the proof condition met at every size; `minweight_scan.py`,
`results/minweight_scan.txt`):

| n | 30 | 40 | 50 | 60 | 70 (IBM's instance) |
|---|---|---|---|---|---|
| minimum stabilizer weight d_min | 6 | 8 | 9 | 12 | 15 |

So the smallest marginal that carries any signal grows roughly like n/5: at 100 qubits one would expect every
marginal on about 20 qubits to be exactly uniform. The ratio d_min/n ≈ 0.19–0.21 at n = 60–70 is close to the quantum
Gilbert–Varshamov rate for random stabilizer codes of rate zero (δ ≈ 0.19), as one would expect of a random-like
construction; we note this as an observation.

The same fact turns the hardware shots into a simulation-free diagnostic, but of a restricted kind. A Pauli error
anywhere in the circuit, propagated through the Clifford part, acts on |G⟩ as a Pauli operator and leaves a
stabilizer state with the same group up to signs; by the proposition its marginals on fewer than 15 qubits are still
exactly uniform, whatever the error rate, and so are those of any mixture of such states. Pauli-type noise
(depolarising, dephasing, Pauli errors on the gates, and the frames selected by the Pauli checks) therefore cannot
bias these marginals at all. A coherent mis-rotation of the final layer cannot either, since the "if" direction holds
in every product basis. What can is asymmetric readout (a bit that reads 1 more often than 0 regardless of the state)
and non-Pauli channels such as amplitude damping or leakage. For IBM's 1,726 post-selected shots all 70 single-bit
means (χ² 61 on 70, p = 0.78), 2,415 pair parities (p = 0.36, largest |z| 4.4), 54,740 triple parities (p = 0.28) and
916,895 quadruple parities (p = 0.95, largest |z| 5.6, which is p ≈ 0.02 after correcting for the number of tests) are
at chance (`results/hw_marginals.txt`). What 1,726 shots can detect in one test is a marginal departing from ½ by about 3.6 percentage points (3σ); across
all 974,120 tests the family-wide 5% threshold is about 5.4σ, or 6.6 points. So the statement is: after post-selection,
readout asymmetry and other non-Pauli effects on IBM Boston are below about 6.6 points on every qubit and every parity
of up to four qubits, and below about 3.6 points on any single one chosen in advance. The χ² values treat the z-scores as
independent, which overlapping parities are not, so their p-values are approximate. The test says nothing about
the Pauli-type noise that presumably accounts for most of the gap between fidelity 0.46 and 1, which is why it is
consistent with, but no evidence for, global depolarisation (Figure 5).

![Low-order parities of the hardware shots](fig/fig4_marginals.png)

*Figure 5. Distribution of z-scores of all 1-, 2-, 3- and 4-bit parity means over IBM's 1,726 shots and over our 1,012 exact samples, against the standard normal that an exactly uniform marginal implies. The χ² p-values are those quoted in this section for the shots and in §4 item 7 for the samples (`hw_marginals.py`).*

## 7. The doped instance

IBM's flagship variant, `nq70_depth70_checks27_doped`, inserts 468 T gates throughout the same 70-qubit brickwork and
measures in Z; its README calls it adversarial for both tensor-network and stabilizer methods (best known stabilizer
rank 2^{185.5}) and reports a hardware fidelity of 0.26 ± 0.02 after post-selection. The only classical entry on the
tracker computes exact amplitudes by tensor-network contraction (#242: contraction width 35, 2^{48.1} complex FMAs per
batch of 256 amplitudes, all 2,051 IBM shots in 37.3 minutes on 32 nodes of eight H100 GPUs, log-XEB 0.350 ± 0.027);
#258 sampled it at reduced precision on one desktop GPU [19], and #261 states that its half/half method does not
apply to it.

**T gates as projected ancillas.** For a data qubit q, a fresh ancilla a in |0⟩, CNOT(q → a), and the projection of a
onto |φ⟩ = (|0⟩ + e^{−iπ/4}|1⟩)/√2 give ⟨φ|_a CNOT |ψ⟩_q|0⟩_a = T|ψ⟩/√2 for every |ψ⟩ (the orthogonal outcome gives
ZT|ψ⟩/√2, so each projection succeeds with probability exactly ½ whatever the state). The projection is a product
measurement: apply T then H to the ancilla and keep outcome 0. Hence the doped circuit is a Clifford circuit on
70 + 468 = 538 qubits followed by a product measurement, a graph state measured in a product basis, and
Pr(all 468 projections succeed) = 2^{−468} exactly, so Pr_doped(b) = 2^{468}·Pr_gadget(b, anc = 0). This gadget is the
T-injection that Renaud's `dcs-rankwidth` (tracker #249) applies with a linear rank-width order to a different doped
instance (64 qubits, 314 T gates, width 36, probabilities on 192-core, 3 TB nodes); what is added here is its
application to the 70-qubit flagship with an order of width 34 on one 128 GiB node, the string-by-string comparison
with #242, and the hybrid sampler below. We checked the
identity on 20 random circuits against dense statevectors (1e-10), the whole pipeline on 25 more (1e-15), and the
sampler binary itself on 12 (1e-8), each with a control that drops the ancilla's T (`doped_gadget.py`, `test_doped.py`,
`test_gsamp2_doped.py`).

**Width.** On the 538-vertex graph the natural order has width 137, but a spatial sweep (qubit 0 to 69, each qubit's
ancillas in time order then its data vertex) has width **34** (`results/doped/doped_widths.json`): the cut between the
left and right halves of the chain is crossed only by the CZ gates on one bond, at most 35 of them, and the sweep
realises nearly that bound; simulated annealing from the sweep and single- and double-vertex slice searches found
nothing lower. One exact amplitude therefore costs a table of 2^{34} complex64 (128 GiB) and 2^{40.9} table updates.
For comparison, #242's contraction has width 35 and costs 2^{48.1} complex multiply-adds per batch of 256 amplitudes,
2^{40.1} per amplitude: the two operation counts per amplitude (table updates here, complex multiply-adds there; not
a like-for-like FLOP count) are within a factor of two, and the gap in wall-clock time (their 37 minutes for all 2,051
shots on 256 H100 GPUs against our 19.4 minutes per shot on one 32-core node, run 11) is mostly hardware.

**Why the chain rule fails here, and what replaces it.** The orthonormality identity traces out the unmeasured suffix,
but the ancillas must be projected, not traced (a traced ancilla dephases its T gate). Putting all ancillas before any
sampled data qubit makes the prefix–suffix cut contain the whole ancilla|data cut, whose rank is 70; exact chain-rule
sampling would need 2^{70} entries. Instead we move the last m = 20 data vertices to the end of the sweep order, which
leaves the width at 34 and makes the final cuts 20, 20, 19, …, 1, and sample in two stages: a proposal draws the
first 50 data bits uniformly; one pass gives their exact marginal Pr(prefix, anc = 0) (the 20 unmeasured data qubits
are traced by orthonormality) together with exact draws of the last 20 bits from their conditionals; the proposal is
accepted with probability z/c where z = 2^{50}·2^{468}·Pr(prefix, anc = 0) is the marginal in units of the uniform
distribution. Because z averages 2^{20} output probabilities it is concentrated near 1 (relative spread about
2^{−10}), so with c = 1.03 almost every proposal is accepted. A proposal with z > c would be accepted with
probability one, so the output law is q ∝ min(z, c) on the prefix times the exact suffix: the Born law if no prefix has
z > c, and in general within total variation ε = E_uniform[(z − c)⁺] of it (with w = min(z, c) and m = E[w] = 1 − ε,
TV(q, p) ≤ ½(E|w/m − w| + E|w − z|) = ε), with classical fidelity at least m, since √(wz) ≥ w. We call the output capped
rather than exact. The difference from #260 is one of degree: there the cap sits at O(1) times the spread of the capped
quantity, here at about 30 standard deviations.

**A certificate from the proposals.** Each proposal is a uniform draw of w = min(z, c) ∈ [0, c], so a one-sided
confidence bound on its mean bounds ε without any model of z. With the empirical Bernstein inequality [17] at 95%,
the proposals give fidelity ≥ 0.940 for the 70-qubit hybrid check (150 proposals), 0.887 at 80 qubits (80), 0.936 for
the doped instance (139) and 0.912 for the ring (102) (`hybrid_certificate.py`, which also checks the inequalities on
exact small laws and the coverage of the bound; `results/certificates_nq80_fstate.txt`,
`results/doped/certificate_run11_13.txt`, `results/ring64/run14/certificate.txt`). The bound is limited by the number
of proposals, not by z: its range term is about 7c·ln(40)/(3n).

For the 70-qubit instance the excluded mass is now bounded without a model. The second moment of z over uniformly
random prefixes is a sum over the stabilizers supported on the prefix: E[z²] = Σ_U 3^{−|U|} (Σ_{S: supp S = U} v(S))²,
where U runs over supports, v(S) = ±1 is the sign of S times the signs of the measured Bloch components on its
support, and pairs of stabilizers with different supports cancel in the sum over the prefix bits (`capvar.py`,
checked against brute force over every prefix string on 38 random circuits to 4·10⁻¹⁵ with a flipped-sign control,
`results/tests/test_capvar.txt`). For the first 50 qubits of the width-30 order that group has 2^{30} elements;
enumerating all of them by a Gray code with phase tracking (`capvar.cpp`, 24 minutes on the laptop) gives
E[z²] = 1.000000970065, i.e. Var(z) = 9.70·10⁻⁷ against the 1.00·10⁻⁶ observed over the 150 proposals; 98.4% of the
elements have a support of their own and the cross terms sum to −2·10⁻¹¹ (`results/capvar_nq70_K50.txt`). Chebyshev
then gives Pr_uniform(z > 1.03) ≤ 1.08·10⁻³, and by Cauchy–Schwarz the Born mass of those prefixes is at most
Pr + √(Var · Pr) = 1.11·10⁻³, so the accepted output of the hybrid sampler on this instance is within total-variation
distance 1.1·10⁻³ of the Born distribution, proven. Since the bound scales as Var/(c − 1)² while the acceptance rate
is about 1/c, a larger cap is cheap: c = 1.3 would give about 1.1·10⁻⁵ at about 30% more proposals. At 80 qubits the
same enumeration has 2^{40} elements with almost as many distinct supports, which an in-memory hash cannot hold; an
external sort over about 10³ disk buckets would, at roughly a day on one node, but we have not run it. For the doped
order it has 2^{485}; there only the proposal certificate above applies. On small T-circuits the hybrid sampler reproduces the exact output distribution
(total variation at the level of an ideal sampler of the same size) and a cap set below the true maximum of z is
detected as a biased output (`test_hybrid.py`). At full scale we ran it on the 70-qubit fstate instance, where the
exact sampler is the ground truth: 150 proposals with the width-30 plan (50 bits proposed, 20 drawn exactly), 147
accepted, z between 0.9970 and 1.0025 (spread 0.0010, as the 2^{−10} estimate predicts), no proposal above the cap;
the accepted samples' XEB against their own exact probabilities is 1.018 ± 0.128 against the exact sampler's
1.020 ± 0.047, their mean log₂ Pr −69.37 ± 0.09 against −69.40 ± 0.04, and three accepted samples re-scored in
prob mode agree with the sampler's reported log₂ Pr to 9·10⁻⁷ (`results/fstate_hybrid_check.txt`,
`results/fstate_hybrid_samples.txt`). Each proposal costs one pass of the width-30 plan.

![Hybrid sampler marginals](fig/fig6_cap.png)

*Figure 6. The prefix marginal z of every proposal of the hybrid sampler on the three instances it was run on (70-qubit Clifford instance, 150 proposals; 80-qubit instance, 80; doped instance, 13), in units of the uniform distribution. The spread is 10⁻³, as the average of 2^{20} Porter–Thomas probabilities predicts, and the cap c = 1.03 lies about 30 standard deviations to the right of every observed value; the mass beyond it is what separates these samples from exact ones (`results/fstate_hybrid_samples.txt`, `results/nq80/samples_nq80.txt`, `results/doped/run11/samples.txt`).*

**Full scale: the first IBM doped shot, seven ways.** The first of IBM's 2,051 published doped shots was scored on
five rented nodes with four different arrangements of the same computation (`graphsamp/results/doped/`, each value
read from that run's output file):

| run | node | plan | log₂ Pr(shot 1) | seconds per pass |
|---|---|---|---|---|
| 7 | r5a.8xlarge, 256 GiB | width 34, m = 20, unsliced | −68.295330880 | 9,872 |
| 7 | same | width 34, m = 10 (different order) | −68.295331335 | — |
| 9 | r7iz.8xlarge, huge pages | m = 20 / m = 10 | −68.295331925 / −68.295330954 | 4,919 |
| 10 | r7a.4xlarge, 128 GiB | two slices at ancillas 135, 104 (4 plans of width 33) | −68.295331107 | 5,888 (4 slice passes) |
| 11 (spot) | r7iz.8xlarge, two-table engine | same two slices, two tables of 2^{33} | −68.295331101 | 1,592 (4 slice passes) |
| 11 (on-demand) | r7a.8xlarge, two-table engine | same | −68.295330552 | 1,168 (4 slice passes) |

The seven values span 1.4·10⁻⁶ in log₂, the rounding level of complex64 tables accumulated in different orders. The
two-table engine's sliced evaluation (run 11, 1,592 s on r7iz.8xlarge) against the in-place engine's unsliced pass on
the same machine type (run 9, 4,919 s) is 3.1 times faster, but engine and slicing both changed between the two, so it
is not an engine-only figure. Runs 9
and 10 then scored the same uniformly random control string (not a hardware shot) at −72.910036757 and −72.910036839,
an agreement of 8·10⁻⁸ between a different machine, layout and slicing; both runs ended before reaching further strings
(run 9 was reclaimed by the spot market, run 10 was stopped to free its quota). An
eighth value comes from a different method altogether: the tensor-network entry #242 published its amplitudes
(Zenodo 10.5281/zenodo.21912448; q0-first strings, which are our Qiskit strings reversed, all 2,051 matching), and its
shot 1 has log₂ p = −68.295074. The relative difference between the two methods is 1.8·10⁻⁴ (`doped_tn242_compare.py`,
`results/doped/tn242_crosscheck_run7.txt`); the 16 further shots of run 11 below show that it is a common factor
across strings (15 of 16 within 2.2·10⁻⁵ in log₂ of the mean offset).

**What the published exact probabilities say about the hardware run.** Since our values agree with #242's, their
2,051 probabilities can be analysed with the estimators of §3.2 (`ibm_fidelity_mle.py`, `doped_tn242_analysis.py`;
outputs in `results/doped/tn242_fidelity.txt` and `tn242_porter_thomas.txt`). Each of their batches also holds the
255 sibling strings that differ from the shot on qubits 0–7, so the ideal distribution itself can be tested: over
523,005 siblings z = 2^{70}p has mean 0.9991 ± 0.0014, mean z² 2.0013 and KS distance 0.0016 to Exp(1) (5% level
0.0019), so the doped circuit's output is Porter–Thomas over strings, as the Clifford instance's is over our samples
(§3.1), and E_ideal[z] = 2 is justified.
For the shots: linear XEB fidelity 0.342 ± 0.028, mixture maximum likelihood 0.343 with 95% interval [0.293, 0.393],
and the KS distance of the shots' z to the fitted mixture (F z + 1 − F)e^{−z} is 0.014 (5% level 0.030), so the
depolarizing mixture describes the data; the log-XEB of #242 is 0.350 ± 0.027. A local check uses the conditional
distribution of qubits 0–7 given the other 62 bits: the shots' mean of 256·p(x)/Σ_siblings p − 1 is 0.339 ± 0.028,
the same fidelity read off eight qubits at a time. The three agree with each other; the linear-XEB value 0.342 ± 0.028
was first quoted by #258 [19] from #242's amplitudes. They are above the README's 0.26 ± 0.02, whose own 95% lower bound
of 0.036 suggests a larger uncertainty, and which was obtained on the device by other means.

**Sixteen more shots, and hybrid samples (run 11).** Run 11 (r7a.8xlarge on-demand, 32 cores, 256 GiB;
`results/doped/run11/`) first reproduced shot 1 at −68.295330552 (1,168 s, peak resident set 134,262,660 kB), then
scored the last 16 of IBM's 2,051 shots at 1,165 s each with the two-table engine. Against #242 the 16 differ by
almost the same amount: log₂ difference −2.56·10⁻⁴ on average, 15 of 16 within [−2.76, −2.54]·10⁻⁴, and the one
exception, the least probable of the 16 (log₂ p = −74.73, z ≈ 0.04), at −1.63·10⁻⁴ (`tn242_compare.txt`). The mean
offset is therefore a common normalisation factor of 1 − 1.8·10⁻⁴ between the two pipelines rather than per-string
rounding; the spread of 2.2·10⁻⁵ among the 15, and the larger deviation on the smallest amplitude, are what absolute
rounding at the level of a 2^{48} complex64 contraction would produce on top of it. The factor's origin is not identified: float32 rounding of the repeated gate constants (350 H at 1/√2,
468 rz(π/4)) predicts a factor of 4·10⁻⁵ of the opposite sign, and our seven in-house values of shot 1 agree to
1.4·10⁻⁶. It shifts #242's log-XEB by 1.8·10⁻⁴. Recomputing two amplitudes with complex128 tables (run 15, the same
two-slice plans; `results/doped/run15_c128/`) moves our values by 2.5·10⁻⁷ (shot 1, −68.295331132) and 1.4·10⁻⁶ (the
last shot, −69.167144551) and leaves the offset from #242 at −2.57·10⁻⁴ and −2.60·10⁻⁴ in log₂, so the factor is not
complex64 rounding here; it lies on #242's side or in a convention between the pipelines. The same pipeline reproduces
Renaud's independent probabilities of the ring instance to 3.8·10⁻⁶ (below). The 16 shots' own linear XEB is 0.33 ± 0.25 (bootstrap 95%
[−0.14, 0.83]), consistent with the 0.342 ± 0.028 from all 2,051 published amplitudes (`analysis.txt`). The run then
drew samples of the doped instance with the hybrid sampler: 13 proposals in 15,183 s, 13 accepted, prefix
marginals z between 0.9981 and 1.0007 against the cap of 1.03 (mean 0.9996, sd 0.0008), 1,166 s per proposal
(`samples.txt`); their log₂ Pr averages −69.59 ± 0.36 (Porter–Thomas −69.39) and their self-consistency XEB is
0.81 ± 0.36, consistent with 1 at this sample size (`porter_thomas_samples.txt`). As at 80 qubits these follow the
capped law, not the exact one. A dedicated run on the same node class (run 13,
`results/doped/run13/`; gate on shot 1 to 10⁻⁶) added 126 proposals at 1,194 s each, 122 accepted, none above the cap,
z between 0.99783 and 1.00269. Together runs 11 and 13 give 135 distinct samples from 139 proposals (z sd 0.00093);
their self-consistency XEB is 1.15 ± 0.13 and their mean log₂ Pr −69.28 ± 0.10 against the Porter–Thomas −69.39. They
are not the first samples of this instance: tracker entry #258 (8 October 2026) drew 74 samples with a chain-sweep
contraction in 5-bit block floating point on one desktop GPU, at a predicted fidelity of about 0.87 set by the
rounding. Ours come from complex64 amplitudes, and their 139 proposals certify the sampler at fidelity ≥ 0.936 (95%).

**IBM's ring instance (#245, run 14).** The current doped candidate has 64 qubits on a ring, depth 73 and 314 T gates
(378-vertex gadget graph). Renaud's published order (#249) gives width 36 on our graph with identical vertex numbering;
slicing ancillas 317 and 211 brings it to 34 (four plans, no twist; `cloud_ring64/`). On one 256 GiB node the gate
re-scored three of IBM's hardware shots against #249's exact probabilities to 3.8·10⁻⁶ in log₂ (2,425 s for three);
an earlier attempt failed this gate and exposed a sign error for an edge between two slice vertices, fixed and checked
against a dense simulation on an 8-qubit piece of the circuit (`results/tests/test_slice_adjacent.txt`). The run then
drew 102 hybrid samples from 102 proposals (817 s each, z ≤ 1.00236), all distinct, Porter–Thomas-consistent (mean
2⁶⁴Pr 1.91 ± 0.14), with linear XEB 0.91 ± 0.14 against 0.390 ± 0.034 for IBM's 1,389 shots, and with TV ≤ 0.088,
classical fidelity ≥ 0.912 at 95% from the proposals alone (`results/ring64/run14/`). Tracker entry #263 [20] sampled
the same instance first, while this run was in progress: 1,389 samples in 864 s on 8 H200 GPUs with a uniform prefix,
at certified fidelity ≥ 0.908. Run 14 is an independent CPU implementation that agrees with both #249 and #263.

## 8. What this does and does not show

It shows:
- the 70-qubit instance is sampled exactly, from the Born distribution, on a laptop, so it does not
  demonstrate a task beyond classical reach;
- the published hardware data can be benchmarked exactly at full depth, which the experiment's README says can only
  be extrapolated, and the resulting fidelity estimate is consistent with IBM's and about four times narrower;
- at depth equal to the qubit count the construction stays within one node's reach to about 80 qubits (exact samples
  at 76; hybrid samples at 80 on a 256 GiB node; a 64 TB server would hold the fitted width to about 93), and no choice
  of graph changes this by more than a few qubits; on a chain, depth rather than qubit count sets the width
  (Proposition 2);
- the output distribution is exactly uniform on every set of at most 14 qubits (the graph state is a [[70, 0, 15]]
  code), and IBM's shots show no readout-type bias above the 4-point level that 1,726 shots can detect (Pauli noise
  cannot show in these marginals, so this is not a statement about the noise that sets the fidelity);
- the doped flagship's amplitudes are computed exactly on one node after the T-gadget (width 34), agreeing with the
  tensor-network entry #242 string by string up to a common factor of 1 − 1.8·10⁻⁴ that complex128 recomputation shows
  is not rounding on our side.

It does not show:
- a new complexity-theoretic result; the rank-width reduction is from 2007 and the ceiling is elementary;
- priority on attacking this instance; #260 applied the same reduction with capped rejection sampling the
  day before our first run, and we cite it throughout; nor priority in sampling the doped instances, which #258
  (doped 70-qubit, low precision) and #263 (ring) did first;
- exact samples of the doped instances or at 80 qubits: the hybrid sampler is capped, its error is proven below
  1.1·10⁻³ only at 70 qubits, and elsewhere it is certified statistically (fidelity ≥ 0.887–0.936 at 95%);
- a bound for tree decompositions, whose rank-width is lower than the linear one on random graphs (§5), or anything
  about graph-state experiments with entangled measurements.

**What we would suggest the experimenters run next.** On a chain, width is set by depth: Proposition 2 caps it at
⌈d/2⌉ whatever the qubit count, so escaping a 64 TB server needs both d and n above about 90, and adding
qubits at fixed depth buys nothing. Layouts whose cuts are crossed by more gates per layer (two-dimensional or
long-range couplings) raise D, and with it the bound; at n ≥ 100 and d ≥ n the ceiling then puts any graph beyond
linear rank-width methods (9 PB), with the caveat of §5 about tree decompositions. Doping with a few hundred T gates
does not help by itself: the T-gadget keeps the width at 34 for the doped instance (§7), as Proposition 2 predicts. The exact scorer here
computes full-depth XEB for any graph-state run up to about 80 qubits on one node (§3.4), so an experiment at that
scale can report an exact benchmark instead of an extrapolated one.

## 9. Related work and references

The method is the linear-order case of the rank-width simulability of graph-state measurements [1], with the Schmidt
structure of graph states across a cut [2] supplying the prefix marginals and chain-rule sampling in the manner of
perfect sampling for tensor networks [4]. The stabilizer-code framing of §6 is the equivalence between k-uniform states
and pure codes of distance k + 1 [3]. Rank-width contraction of ZX diagrams [5] is the general form of the same idea.
On the tracker, Renaud [6] applied the T-gadget with a linear rank-width order to a different doped instance (64 qubits,
314 T gates, width 36, probabilities only); Qiu [7, 8] attacked the 70-qubit instance of §1.1 with capped rejection
sampling and then with exact half/half sampling on 8 H200 GPUs (1,726 samples in 1,265 s; the half/half split is the
ceiling of §5 with no order search, 32 times the laptop's table, and its issue states that it does not apply to the
doped instances); Manabe, Gu and Pan [9] computed exact amplitudes of the doped instance by tensor-network contraction
on 256 H100 GPUs and published them [10], which §7 compares string by string. The experiment itself is described in
[11] and [12]; the random-graph-state hardness argument behind the construction is [13]. Tree-shaped rank
decompositions, and their cost, are treated in [5, 15, 16], and the rank-width of random graphs in [14]. On the doped
instances, #258 [19] drew low-precision samples of the 70-qubit one and #263 [20] certified samples of the ring. Our
first exact samples of the 70-qubit instance are dated 29 October 2026 (`results/run.log`); [8] was posted first.

1. M. Van den Nest, W. Dür, G. Vidal and H. J. Briegel, Classical simulation versus universality in measurement-based
   quantum computation, Phys. Rev. A 75, 012337 (2007).
2. M. Hein, J. Eisert and H. J. Briegel, Multiparty entanglement in graph states, Phys. Rev. A 69, 062311 (2004).
3. A. J. Scott, Multipartite entanglement, quantum-error-correcting codes, and entangling power of quantum evolutions,
   Phys. Rev. A 69, 052330 (2004).
4. A. J. Ferris and G. Vidal, Perfect sampling with unitary tensor networks, Phys. Rev. B 85, 165146 (2012).
5. F. Kuyanov and A. Kissinger, Efficient classical simulation of low-rank-width quantum circuits using ZX-calculus,
   arXiv:2603.06764 (2026).
6. K. Renaud, injected_graph_scan_nq64_depth73 (code `dcs-rankwidth`), Quantum Advantage Tracker issue #249,
   https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/249 (22 September 2026).
7. Z. Qiu, capped rejection sampling of random_graph_sampling_nq70_depth70_checks27_basis_fstate, Quantum Advantage
   Tracker issue #260 (30 October 2026; closed by the author on 2 November 2026 in favour of [8]).
8. Z. Qiu, Exact sampling of random_graph_sampling_nq70_depth70_checks27_basis_fstate: 1,726 samples in 1,265 s on
   8 H200 GPUs, Quantum Advantage Tracker issue #261,
   https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/issues/261 (2 November 2026).
9. H. Manabe, H. Gu and F. Pan, Classical simulation and design frontiers for IBM's doped Clifford sampling experiment,
   arXiv:2608.13110 (2026); tracker issue #242.
10. H. Manabe, H. Gu and F. Pan, Data for exact tensor-network simulation of IBM's doped-Clifford circuit, Zenodo,
    doi:10.5281/zenodo.21912448 (2026), CC BY 4.0.
11. Quantum Advantage Tracker, random_graph_sampling_nq70_depth70_checks27_basis_fstate (issue #151, 29 April 2026) and
    doped_random_graph_sampling_nq70_depth70_checks27 (issue #228, 9 July 2026),
    https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io.
12. S. Martiel, J.-U. Chung, A. Seif, S. Ghosh, I. Hincks, A. Deshpande, H. Manabe, H. Gu, F. Pan, B. Fefferman,
    J. M. Gambetta and A. Javadi-Abhari, Sampling hard circuits with verifiably high fidelity, arXiv:2607.25941 (2026).
13. S. Ghosh, D. Hangleiter and J. Helsen, Random regular graph states are complex at almost any depth,
    arXiv:2412.07058 (2024).
14. C. Lee, J. Lee and S. Oum, Rank-width of random graphs, J. Graph Theory 70, 339 (2012), arXiv:1001.0461.
15. J. Codsi and T. Laakkonen, Unifying graph measures and stabilizer decompositions for the classical simulation of
    quantum circuits, arXiv:2603.06377 (2026).
16. A. de Colnet, F. Geerts, R. Hai, A. Laarman, J. H. Lee and G. A. Pérez, Quadratic sums-of-powers for
    fixed-parameter tractable quantum-circuit simulation, arXiv:2605.29944 (2026).
17. A. Maurer and M. Pontil, Empirical Bernstein bounds and sample variance penalization, COLT (2009).
18. M. Matsumoto, I. Wada, A. Kuramoto and H. Ashihara, Common defects in initialization of pseudorandom number
    generators, ACM Trans. Model. Comput. Simul. 17, 15 (2007).
19. GitHub user dylanneve1, qsimlab_packed_chainsweep_doped70, Quantum Advantage Tracker issue #258 (8 October 2026).
20. Z. Qiu, Classical sampling of ibm_boston_dcs_nq64_depth73_checks12 (#245), Quantum Advantage Tracker issue #263
    (5 November 2026).
21. Hewlett Packard Enterprise, HPE introduces industry-first 64 TB memory server (press release, May 2026).

## Appendix: the artifact behind each claim

All paths are relative to the repository root (`graphsamp/` in the working repository; the root of the public
repository). Commands run from that root; `results/run.log` and the per-run `README.md` files record the exact command
lines of the cloud runs and the machines they ran on.

| Claim (section) | Artifact | Produced by |
|------|-----------|--------|
| 1,012 exact samples, 13,202 s, peak 34.4 GB (§3.1) | `results/samples_nq70.txt`, `run.log`, `time_nq70.txt`, `inputs.sha256` | `gsamp` sample mode on `cloud/nq70_plan.txt` (command in `run.log`) |
| Bit-for-bit reproduction of 55-qubit samples (§4.3) | `results/check55_result.txt` | `xcheck55.py` |
| Laptop: 3 samples in 55 s, 8.6 GB (§3.1) | `results/laptop_samples_w30.txt` | `gsamp2 cloud2/nq70_inplace_plan.txt sample 3 2000000` |
| Independent engine: 5 cloud samples to 10⁻⁶, 10-bit control, 1 laptop sample (§4.2) | `results/xcheck70_independent.txt`, `xcheck70d.txt`, `xcheck70e_laptop.txt` | `xcheck70c.py`, `xcheck70d.py`, `xcheck70e.py` |
| Self-XEB 1.020 ± 0.047 (§3.1) | `results/ibmxeb/analysis.txt` | `analyze.py results/samples_nq70.txt cloud/nq70_step_qubits.json` |
| Porter–Thomas test of the samples (§3.1) | `results/porter_thomas_nq70.txt` | `porter_thomas_test.py results/samples_nq70.txt 70` |
| Reseeded rerun: 1,012 samples, marginals at chance (§4.7) | `results/samples_nq70_reseeded.txt`, `time_nq70_reseeded.txt`, `hw_marginals_reseeded.txt`, `porter_thomas_nq70_reseeded.txt` | `gsamp2 ... sample2 1012 3000000` (Docker, laptop); `hw_marginals.py`; `porter_thomas_test.py` |
| IBM's 1,726 shots: XEB 0.466 ± 0.031, controls, #260 agreement, 6.9 s/shot (§3.2) | `results/ibmxeb/*_out.txt`, `time_ibm.txt`, `gate_result.txt`, `analysis.txt` | `ibm_xeb_prep.py`; `gsamp2 ... prob`; `ibm_xeb_analyze.py` |
| Fidelity estimators and their synthetic calibration (§3.2) | `results/ibmxeb/analysis.txt`, `results/fidelity_mle_synthetic.txt` | `ibm_fidelity_mle.py`, `fidelity_mle_synth.py` |
| Per-step conditionals, stabilizer-only check, cross-order check, information profile (§4.5, §4.6, Figure 2) | `results/chain/` (`stabprefix_*.txt`, `analysis_crossorder.txt`, `analysis_profile.txt`) | `gsamp2 ... chain`; `stabprefix.py`; `chain_profile.py` |
| Engine-independent uniform-prefix check (§4.4) | `results/uniform_prefix_check.txt` | `check_uniform_prefix.py` |
| Low-order marginals of shots and samples, biased-bit control (§4.7, §6, Figure 5) | `results/hw_marginals.txt`, `hw_marginals.json`, `hw_marginals_z.npz` | `hw_marginals.py` |
| Random-stream identity and scan; drawing-step test (§4.7) | `results/rng_check.txt`, `rng_scan.txt`, `rngbits_seeds_*.txt`, `draw_check.txt` | `rngbits.cpp`, `rng_check.py`, `draw_check.py` |
| 76-qubit run: 19 samples, 190 s each, independent check (§3.3) | `results/nq76/` (`samples_nq76.txt`, `time_sample.txt`, `xcheck_out.txt`, `analysis.txt`) | cloud run 5/6 (`run.log`); `report_run.py` (`xcheck_out.txt` holds re, im of the amplitude scaled by 2^{n/2}; log₂ Pr = log₂(re² + im²) − n) |
| 80-qubit run: gates, 79 hybrid samples, re-scored strings (§3.4) | `results/nq80/` (`samples_nq80.txt`, `xcheck_nq80.txt`, `gate*_result.txt`, `time_sample.txt`, `analysis.txt`, `porter_thomas_nq80.txt`) | cloud run 12 (`cloud_nq80/run_hybrid.sh`); `nq80_analyze.py` |
| Dual-table translation and twist checks (§3.4) | `results/tests/test_dual_table.txt`, `test_twist.txt`; `results/nq80/gate_tw_result.txt` | `test_dual_table.py`, `test_twist.py`; run 12 gate |
| Frontier widths n = 30–120 (§5, Figure 3) | `results/frontier.jsonl`; `results/frontier_instances/*.log` | `frontier.py`; `frontier_instance.py` + `order_search.py` |
| Ceiling tightness on random graphs (§5) | `results/ceiling_tightness.jsonl`, `ceiling_tightness.cmd.txt` | `ceiling_tightness.py` (command in the `.cmd.txt`) |
| Cut-rank profiles (Figure 4) | `cloud/nq70_step_qubits.json`, `results/nq70_graph_adjacency.txt`, the in-place plans | `paper/make_figures.py` |
| Minimum stabilizer weight 15, three-way certificate, validation on 30/40 qubits, scan (§6) | `results/minweight_exact.txt`, `minweight_k11.txt`, `minweight_case1_k11.txt`, `minweight_split_test.txt`, `minweight_scan.txt`, `nq70_graph_geninverse.txt`, `nq70_graph_kernel.txt` | `minweight_exact.py`, `minweight.cpp`, `minweight_case1.cpp`, `minweight_split_test.py`, `minweight_scan.py` |
| Doped instance: gadget identity, widths (§7) | `results/tests/doped_gadget.txt`, `test_doped.txt`, `test_gsamp2_doped.txt`; `results/doped/doped_widths.json` | `doped_gadget.py`, `test_doped.py`, `test_gsamp2_doped.py`; `doped_to_qasm.py` + order sweep |
| Doped shot 1 seven ways; 16 further shots; 135 samples (§7) | `results/doped/run7/`, `run9/`, `run10/`, `run11/`, `run13/` (`shots_scored.txt`, `samples.txt`, `tn242_compare.txt`, `analysis.txt`, `porter_thomas_samples.txt`) | `cloud_doped/run*.sh`; `doped_sampler.py`; `doped_tn242_compare.py`; `doped_analyze.py` |
| Exact second moment of the prefix marginal at 70 qubits, K = 50: Var(z) = 9.70·10⁻⁷, TV bound 1.1·10⁻³ (§7) | `results/capvar_nq70_K50.txt`, `capvar_nq70_K50_gens.txt`, `results/tests/test_capvar.txt` | `capvar.py export ...; capvar` (C++), `test_capvar.py` |
| Gate-cut bound, widths along the spatial order for n ≠ d (§5, Proposition 2) | `results/gatecut_check.txt` | `gatecut_check.py` |
| Proposal certificates for the hybrid sampler (§7) | `results/certificates_nq80_fstate.txt`, `results/doped/certificate_run11_13.txt`, `results/ring64/run14/certificate.txt` | `hybrid_certificate.py` (`--selftest`) |
| complex128 recomputation of two doped amplitudes (§7) | `results/doped/run15_c128/` | `cloud_doped/run_c128.sh`, `gsamp2.cpp -DGS_DOUBLE` |
| Ring instance #245, run 14 (§7) | `results/ring64/run14/` | `cloud_ring64/run_ring.sh`, `ring64_xeb_compare.py`, `porter_thomas_test.py` |
| Hybrid sampler at 70 qubits against the exact sampler (§7, Figure 6) | `results/fstate_hybrid_check.txt`, `fstate_hybrid_samples.txt` | `fstate_hybrid_check.py` |
| #242's amplitudes: Porter–Thomas, fidelity three ways, local 8-qubit estimator (§7) | `results/doped/tn242_log2p.txt`, `tn242_fidelity.txt`, `tn242_porter_thomas.txt`, `tn242_log_xeb.json` | `doped_tn242_analysis.py` on the Zenodo release |
| Small-circuit correctness tests with negative controls (§4.1) | `results/tests/*.txt` (each ends with the script's exit code) | `test_*.py` |
| Test table, reproduction commands, sampler modes | `README.md` | |
