# Circuit instance description (Doped Random Graph State Sampling):

`nq70_depth70_checks27_doped.qasm`: This prepares a T-doped random graph state on a 70 X 70 circuit (70 qubits with CZ-depth 70). The 70 logical qubits are arranged on a 1D lattice (LNN) and 27 ancilla qubits are used for error detection. A maximum of 468 nontrivial T gates can be placed throughout the circuit. The best known Schmidt rank is $2^{30}$ and the best known stabilizer rank is $2^{185.5}$. After post-selection, 0.038% of the shots remain, and the state fidelity is 0.26 $\pm$ 0.02 (bounded above 0.036 with 95% confidence).

- `nq70_depth70_checks27_doped_checks.qasm`: This includes the ancillas and spacetime Pauli checks for the above circuit.

## Doped Random Graph State Sampling

We present circuits for doped random graph state sampling. This generalizes [random graph state sampling (RGS)](https://github.com/quantum-advantage-tracker/quantum-advantage-tracker.github.io/tree/main/data/classically-verifiable-problems/circuit-models/random_graph_sampling) by inserting non-Clifford gates throughout the circuit, instead of reserving them for a final basis rotation. This increases the amount of magic and therefore the complexity of extended stabilizer simulations.

The arguments for computational complexity largely follow that of RGS, although we highlight differences in the proof of verifiability (see section C).

## A. Computational Complexity

Our work, like [random circuit sampling (RCS)](https://www.nature.com/articles/s41586-019-1666-5), relies on the asymptotic hardness of classical algorithms for random quantum state sampling. However, instead of Haar-random states, we focus on random graph states which permit both error detection and verifiability in the experiment.

It [has been shown](https://arxiv.org/pdf/2412.07058) that sampling random _regular_ graph states in a random product basis is average-case #P-hard. This implies that this task is hard for _any_ classical algorithm, unless there is a collapse in the polynomial hierarchy. The problem can also be reduced to the simulation of measurement-based quantum computing (MBQC).

In this work, we dope the graph state preparation circuits with T gates, distributing them throughout the bulk of the circuit. T-doped Clifford circuits have been used to study [quantum chaos](https://arxiv.org/pdf/2102.08406), converging to relative error unitary k-designs (and therefore approximating Haar random states) when using a linear number of non-Clifford gates. We conjecture that sampling from these states in the Z-basis is also #P-hard - and, as noted above, the larger amount of available magic resources implies an even greater complexity for classical simulations.

## B. Hardness of classical simulation

For the specific problem size that we ran, we quantify the hardness of simulation through two quantities: entanglement and magic. The former leads to the difficulty of tensor network-based methods, and the latter the difficulty of Clifford-based methods.
In random graph state sampling, both of these properties scale with the system size, implying that classical simulations using known methods will become effectively intractable past some finite circuit size.

### `Quantifying Entanglement`

We quantify the entanglement of the graph state by estimating the Schmidt rank, the $GF(2)$ rank of its adjacency matrix, across random bipartitions. For an $N$-qubit graph state, the Schmidt rank is upper bounded by $2^{N/2}$. As this corresponds exactly to the bond dimension, the contraction cost of tensor network simulations will scale exponentially with this quantity.

To maximize the entanglement in our graph state, we use an ansatz of repeating layers of odd/even CZ gates (brickwork layout), followed by random $\sqrt{X}$ or $S;\ \sqrt{X}$ rotations on each qubit. We do this until the circuit reaches a CZ-depth of $N$. It is known that this depth is enough to prepare _any_ graph state on an LNN architecture (to appear).

By taking the minimum over 100 million random bipartitions, we numerically verify that our circuit nearly saturates the entanglement upper bound with a Schmidt rank of $2^{30}$.

<p align="center">
  <img width="590" height="455" alt="image" src="https://github.com/user-attachments/assets/53d47789-5109-4678-995d-26dd1b5368dc" />

<em>Figure 1. Quimb matrix product state (MPS) contraction times for the 70 X 70 circuit with increasing depth (maximum depth 24). Linear extrapolation from a logarithmic plot of the data (R<sup>2</sup> > 0.98) yields a predicted contraction time of 10<sup>25</sup> seconds.</em>

</p>

Note that we do not directly measure the entanglement for the doped random graph state. However, T gates can be modeled as random multi-qubit Pauli rotations when pulled to the end of the circuit, which are unlikely to disentangle the initial Clifford block.

### `Quantifying Non-stabilizerness`

Alternatively, stabilizer rank algorithms can be used, which are oblivious to the amount of entanglement, but scale exponentially with the amount of magic (non-stabilizerness) in the circuit. These methods decompose quantum states into a sum of stabilizers which can be efficiently simulated.

The complexity therefore scales with the minimum number of states used in this decomposition, i. e. the stabilizer rank. For circuits consisting only of Clifford and T gates, the [best known upper bound](https://arxiv.org/pdf/2106.07740v1) is given by

$$stabilizer\ rank \leq 2^{0.3963 \cdot (T\ count)}$$

So, as claimed, the nonstabilizerness increases exponentially with the system size.

(Recall that for random graph state sampling, we instead characterized the nonstabilizerness with the [stabilizer extent](https://arxiv.org/pdf/1808.00128), the minimum $\lVert c \rVert^2$ over all stabilizer decompositions. We do not use this metric here as the stabilizer rank lower bounds this quantity for approximate error simulations.)

<p align="center">
  <img width="570" height="476" alt="image" src="https://github.com/user-attachments/assets/1dea42a2-5c21-4ae4-b0d4-72f78d859e55" />

<em>Figure 2. Extended stabilizer simulation timing in quizx for 70 X 70 circuits with various T-counts (maximum 85). Linear extrapolation from a logarithmic plot of the data (R<sup>2</sup> > 0.98) yields a predicted time of 10<sup>42</sup> seconds to find the probability of one bitstring.</em>

</p>

Furthermore, by doping on a square $N$ X $N$ circuit, non-Clifford gates can be placed across depth $N$, yielding a maximum of $O(N^2)$ (instead of O(N)) magic resources - which further increases the complexity of extended stabilizer simulations.

We also make note of [Clifford Augmented Matrix Product State (CAMPS) simulators](https://arxiv.org/pdf/2412.17209), which combine tensor networks with Clifford tableau simulators. These algorithms reduce the bond dimension necessary to represent the state by propagating magic gates to the front of the circuit, using a tableau for the entangled Clifford bulk and a smaller bond dimension MPS for the magic layer. After the first $N$ magic gates however, the upper bound on the bond dimension necessary for the MPS increases exponentially, as well as the maximum bond dimension necessary to sample bitstring probabilities.

To the best of our knowledge, then, the large Schmidt rank and stabilizer extent of our T-doped graph states will be adversarial for exact or naïve approximate simulations.

## C. Verifiability

For sampling-based experiments, it suffices to show that samples can be drawn from a quantum computer with greater fidelity than through classical means. With doped random graph state sampling, we claim to have a more direct measure of the state fidelity that, relative to prior RCS-style experiments, requires fewer assumptions about the noise in our circuit.

For comparison, RCS uses the fact that the outcome probabilties of Haar-random states is described by the Porter-Thomas distribution, a distinctly non-classical model. The closeness of the sample distribution to the Porter-Thomas is quantified with the linear cross-entropy benchmarking (XEB) score, which compares the output distribution (in the Z-basis) of $M$ sample outputs $x$ to their quantum probabilities:

$$F_{XEB} =  \frac{2^N}{M} \sum_i^M | \bra{0} C \ket{x_i} |^2 - 1$$

As samples drawn noiselessly from the Porter-Thomas distribution have $F_{XEB} = 1$ and samples drawn from the (classical) uniform distribution have $F_{XEB} = 0$, a non-zero fidelity is evidence of quantum behavior.

A key detail is that XEB requires classical simulations for the outcome probabilities, hence can only be extrapolated for the full-depth circuit. The authors present arguments for why this can be extrapolated from smaller depth or less entangled circuits. XEB also is only close to state fidelity under assumptions of weak noise.

[Quantinuum's RCS experiments](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.15.021052) also measure a proxy for fidelity, running a mirrored version of their circuit with half the depth and calculating the probability of return to the expected state (mirror fidelity).

### `Error detection`

To avoid the problem of vanishing fidelities with large circuits, we use the [spacetime error detection protocol](https://arxiv.org/pdf/2504.15725) and mitigate noise in our samples. As the graph state is prepared with Clifford gates, it is possible to augment the circuit with ancilla qubits and insert spacetime Pauli checks on their support at various depths. While these Pauli checks collectively stabilize the circuit, errors will not generally commute with these checks, manifesting as a "syndrome" error on the ancillas. These errored shots can be detected and post-selected out which, as illustrated in the reference above, can result in order of magnitude improvements over bare state fidelity, all with a significantly lower sampling overhead than PEC.

With these spacetime Pauli checks, T gates can not be arbitrarily placed, being restricted to locations in which they will simultaenously commute with all the stabilizers of the checks. In spite of this restriction, we find that a large number of T-gates remain available for the circuit sizes presented.

Our circuits are mapped onto a one dimensional chain with ancillas attached to exactly one data qubit. On IBM's Heron devices, which have the heavy-hex architecture, it is likely that every other qubit will have degree 3. This dense placement of ancillas enables effective error detection, and guarantees that the number of ancillas can scale with the size of the circuit.

<p align="center">
  <img width="594" height="483" alt="image" src="https://github.com/user-attachments/assets/1624d6ee-c943-4fb5-801a-9ad1fb3ac94c" />

<em>Figure 3. Physical layout for the 70 X 70 circuit (with 27 ancilla qubits) on IBM Boston, a Heron R3 device with heavy-hex architecture. The logical qubits (green) are arranged on a 1D chain with dangling ancilla qubits (red).</em>

</p>

### `Measuring Fidelity`

As the graph states can be prepared with high fidelity with error detection, it is efficient to use [direct fidelity estimation](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.106.230501). The fidelity for target state $\sigma$ and noisy output state $\rho$ can be approximated by randomly sampling the expectation values of $M$ random Paulis $P$:

$$F \approx \frac{1}{M} \sum_k^M \frac{\braket{P_k}_{\rho}}{\braket{P_k}_{\sigma}} $$

For stabilizer states, whose expectation values can be bounded, this requires only a constant shot overhead. Given the $O(1/\sqrt{M})$ scaling in uncertainty, we choose enough random stabilizers to bound the fidelity above 1% with 95% confidence.

For non-stabilizer states, expectation values can become arbitrarily small, which requires a commensurately large shot overhead. We argue, however, that the fidelity of the graph state is equivalent to the doped random graph state.

### `Verfiability with Fidelity`

As in RGS, the equality between undoped and doped graph state fidelities requires that:

- The undoped and doped states share the same gates, the only exception being the angle of Z rotations, which toggles between I and T gates.
- On IBM's devices, Z rotations are [virtually implemented](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.96.022330) as framechanges, hence are noiseless.
- We employ Pauli twirling, or randomized compiling, on the graph state preparation circuit. With the assumption of independent error on the single qubit gates, this is enough to ensure that the [process fidelity](https://arxiv.org/pdf/2503.05943) of the doped graph state is equal (up to first order) to the average graph state process fidelity.

We highlight a subtle point when doping our circuits: inserting T gates into the circuit, even if done perfectly, can change the interference between coherent errors within the circuit. When implementing Pauli twirling then, we first separate the circuit into layers of Clifford $C_i$ and T gate $T_i$ layers: $C_L...T_1C_1T_0C_0$. Twirling these Clifford layers (including the ancilla qubits) suppresses the buildup of coherent errors, ensuring that the angle of the T-layer gates does not affect the underlying noise.

To build confidence in this experimentally, we measure the state fidelity (or a proxy) in a few different regimes, aiming to demonstrate that this unaffected by the T-count:

- Regime 1: Zero T gates, for which we use direct fidelity estimation.
- Regime 2: Small constant number $(\ll O(N))$ of T gates, for which we also use direct fidelity estimation.
  - In this regime, there is still relatively low overhead in measuring the fidelity. Classical simulations are done to sample random Pauli observables in an unbiased manner and calculate the corresponding expectation values.
- Regime 3: $O(N)$ T gates, for which we use XEB as a proxy for fidelity.
  - After $O(N)$ T gates, the second moment converges closely to a Porter-Thomas distribution, making XEB a reasonable metric. The number of T gates is chosen such that classical simulations are still tractable.

Provided that classical simulations fail to faithfully sample from our state, we can certify quantum advantage by showing that the fidelity of the graph state - and consequently the T-doped graph state - is bounded above zero.

<p align="center">
  <img width="752" height="530" alt="image" src="https://github.com/user-attachments/assets/b4d5a8bf-b659-4863-8cbd-2747ac1af9db" />

<em>Figure 4. Fidelity for the 70 X 70 circuit on IBM Boston at various T-counts. For T-count=0,5 we utilize direct fidelity estimation with 120 randomly drawn Paulis, and for T-count=70,80 we measure the linear cross entropy of the resultant samples. After rescaling for readout error, we note no negative trend in fidelity is observed as T-count increases. </em>

</p>

Addendum: In Fig. 4 the fidelities are rescaled by readout error. For points using direct fidelity estimation, the observables are rescaled according to [readout error mitigation](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.105.032620). For XEB, we assume that the probability distribution with readout error is uncorrelated with the ideal probability distribution, separate circuit fidelity $F$ into

$$F = F_{ideal}F_{readout}$$

where $F_{ideal}$ is the circuit fidelity with perfect readout error, and $F_{readout}$ is the fidelity of the readout (see supplementary by [Arute et. al](https://arxiv.org/pdf/1910.11333)). After readout error mitigation, each qubit should have symmetric 0 and 1 state error $p_i$, which yields

$$F_{readout} = \Pi_{i=1}^{N} (1 - p_i)$$

which can be used to rescaled the XEB scores.

### `Verifiability with Logical Fault Rate`

The spacetime Pauli checks, in addition to allowing for post-selection of errors and improvements in fidelity, can be used to characterize the rate of logical faults occurring in the circuit. We aim to relate this rate to the fidelity of the circuit and to bound the difference in fidelity between the undoped and doped graph states.

For simplicity, consider the stochastic Pauli fault model. Note that T gates do not disturb the spacetime Pauli checks, implying that the set of logical faults that can be accepted/rejected during post-selection is identical in the undoped and doped cases, and are perfect, implying that the probability of a fault occurring is independent of the T-count.

For the set of accepted shots, faults can be separated into two categories: harmless if they commute with stabilizers, and harmful if they do not and result in a logical error. The fidelity $F$ can be expressed in terms of these probabilities:

$$ F = P(id\ | accept) + P(harmless\ | accept) $$

where $P(accept), P(id),\ P(harmless)$ respectively are the probabilities of acceptance, no fault occurring, and a harmless fault occurring.

As stated above, T gates do not affect the set of accepted faults nor the probability of a fault occurring, hence $P(id\ | accept)$ is equal for the undoped and doped cases.

<p align="center">
  <img width="1011" height="511" alt="image" src="https://github.com/user-attachments/assets/90562100-1b92-4274-a100-8c4607ad7d0e" />

<em>Figure 5. Post-selection rate for the 70 X 70 circuit on IBM Boston shows no dependence on T-count across the measured values. </em>

</p>

However, T gates can affect that a harmless fault will occur, either by converting a previously harmful fault to a harmless one, or, in an adversarial case, converting a previously harmless fault to a harmful one. The difference in fidelity between the graph and doped graph states can therefore be upper bounded

$$ F_{graph} - F_{doped} = P_{graph}(harmless\ | accept) - P_{doped}(harmless\ | accept) \leq max(P(harmless\ | accept))$$

by the maximum probability of a harmless fault occurring.

We conjecture that the harmless faults should not occur often, a fact that we numerically verify below.

<p align="center">
<img width="3680" height="1050" alt="image" src="https://github.com/user-attachments/assets/7d0666d7-6c3d-4245-b0e6-d5cf7686f670" />

<em>Figure 6. Stochastic Pauli fault simulations for the 70 X 70 circuit. The model includes global depolarizing noise, decoherence error (200 μs), and readout error. Left: The probability of acceptance decreases exponentially with error rate. Middle: For accepted shots, the probability of no fault occurring is exponentially suppressed with increasing noise, while the probability of harmless and harmful faults occurring increases. Right: Only a small systematic difference exists between the fidelities of the graph (Clifford) and doped graph states across all error rates.</em>

</p>

The graph state fidelity, then, is close to a lower bound for the doped graph state fidelity, which concludes our argument for verifiability.

We note ongoing work to extend this argument to coherent errors, which would eliminate the requirement of perfect Pauli twirling.

## Institutions

IBM, University of Chicago
