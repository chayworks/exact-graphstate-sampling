# graphsamp: exact sampling of graph states measured in a product basis

Target: `random_graph_sampling/nq70_depth70_checks27_basis_fstate` on the Quantum Advantage Tracker
(IBM / UChicago, issue #151). It reports hardware fidelity 0.37 and extrapolated classical costs of 10^25 s
(MPS) and 10^7 s per sample (extended stabilizer). *Update 2026-10-30:* #260 (posted the same day as this
note) attacks the instance with the same rank-width reduction (max rank 30) plus capped rejection sampling,
which is not an exact Born sampler. What is distinctive here is exact Born sampling on CPU.

## Idea
The circuit is a Clifford body followed by one layer of single-qubit rotations. The Clifford state is a
graph state |G_A> up to local Cliffords, so

    <b|psi> = 2^{-n/2} sum_x prod_j M_j[b_j, x_j] (-1)^{x^T A x}  (upper triangle),

where M_j combines the rotation with the local Clifford. Processing the vertices in some order, the
processed ones affect the rest only through the GF(2) field A[F,P] x_P. That field takes 2^{r_k} values,
where r_k is the cut rank, so the dynamic program holds 2^{r_k} complex numbers, not 2^{2 r_k} as a generic MPS does.

Sampling: after measuring the first t vertices, the remaining state is a combination of
Z^v|G_F>, and these states are orthonormal. So **Pr(b_1..b_t) = 2^{-t} sum_c |P_t[c]|^2**. One exact
sample therefore costs two DP passes, with no right environment.

## Status: what is proven, what is projected

Proven by running it (each test has a negative control that is caught):

| Test | What it checks | Result |
|---|---|---|
| `test_tableau.py` | stabilizer tableau vs statevector | 200 circuits, 2.6e-16 |
| `test_graphform.py` | local-Clifford graph form | 150 circuits, overlap 1 - 2e-15 |
| `test_plan.py` | DP vs brute force | 120 cases, 3e-15 |
| `test_amp.py` | Born probabilities, end to end vs statevector | 60 circuits, 1e-15 |
| `test_sampler.py` | reference sampler: distribution implied by the chain, every string | 2e-15 |
| `test_engine.py` | C++ amplitude engine with slicing | 1e-14 |
| `test_gsamp.py` | C++ sampler, exact probabilities | 6e-9 |
| `test_gsamp.py` | C++ sampler, 20,000 samples | TV 0.0570; ideal sampler at that size: 0.0570 |
| `xcheck55.py` | 55-qubit real-data subgraph: gsamp Pr vs dp.cpp \|amplitude\|^2 | agree to ~4e-7 in log2 |

nq70 structure (`cutprofile.py`, `order.py`):
- linear cut-rank width is 32 in natural order and 31 after search;
- the best balanced cut found has rank 29 (the README says 2^30 is the best known).

Measured throughput (`bench_sampler.py`): on the 55-vertex subgraph (width 26), ~2.5 s per exact sample.
That is 6.4e8 visits/s on 16 threads, sharing the CPU with another job.

**Full instance, run (2026-10-29, `results/`):**
- **Machine and run:** AWS r7a.4xlarge (16 cores, AMD EPYC 9R14, 123 GB). Order `order_best2.json`: width 31,
  6.0e10 visits per sample. **1,012 exact samples in 13,202 s, i.e. 13.0 s per sample**
  (`results/run.log`, `results/time_nq70.txt`). Peak RSS 33.6 GB.
- **Inputs:** hash-identical to this repository's files (`results/inputs.sha256`). Instance QASM SHA-256
  `e0189bdd2e67d1456cb92f5c750b880cbee935220e539c04934906a94535ea5f`.
- **Reproduction gate before sampling:** the three local 55-qubit samples were reproduced bit for bit on a
  different CPU, compiler and thread count (`results/check55_result.txt`).
- **Independent check** (`xcheck70c.py`, `results/xcheck70_independent.txt`): dp.cpp sums sliced amplitudes,
  a different algorithm from gsamp's chain of conditionals. For cloud samples 1 and 2 it gives log2 Pr =
  -69.5085294 and -68.8761651; gsamp reported -69.5085302 and -68.8761651. Control: sample 1 with step-bit 35
  flipped gives -69.5312514. That differs by 0.023, against agreement of about 1e-6. The control is weaker
  than ideal, because this particular flip changes Pr by only 1.6%.
- **Self-consistency statistic:** linear XEB of our own samples is 1.020 +/- 0.047, i.e. 2^70 sum p^2 ≈ 2.02
  (Porter–Thomas would give 2). This is not an independent check.

Not yet done: a check on more samples, and a check of the samples' marginals against something independent
of both engines (for example, Pauli expectations computable with stabilizer methods).

**Unverified:** novelty. Rank-width-based simulation of graph-state measurements is known in principle
(Van den Nest, Dür, Vidal, Briegel, PRA 75, 012337 (2007)). The literature has not been checked for this
exact sampler. The claim we would make is about this instance, not a new complexity result.
