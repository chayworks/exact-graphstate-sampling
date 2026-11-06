# Exact classical sampling of the IQP generative models run on IBM's 156-qubit Aachen processor

Target: O. Balló-Gimbernat et al., *Shallow instantaneous quantum polynomial-time circuits for generative modeling on
noisy intermediate-scale quantum hardware*, arXiv:2511.05267 (v3, March 2026), Phys. Rev. A (April 2026). The abstract
says the approach keeps "the guarantee of sampling hardness"; models are trained classically and sampled on
ibm_aachen at 28, 45, 91 and 153 qubits. Code and data: Zenodo 10.5281/zenodo.18983816 (CC BY 4.0 / MIT), which
contains the trained parameters (`results/trained_params_iqp/`) and the hardware samples
(`results/evaluation_results_iqp/`). The data is not in this repository (`data_iqpnisq/` is git-ignored); MD5 of
`IQP-NISQ.zip`: 6a3eb30b933f7757d32a2b461555daed.

## Result

Their circuits are exp(i Σ θ_j X_{g_j})|0⟩ with one X_i per qubit and one X_iX_j per coupler of the heavy-hex
chip (`src/utils/hardware.py`). All generators commute and act on a graph whose rows are joined by four bridges, so
along the chip's row order the output distribution is a transfer-matrix product with bond dimension 2^6 = 64. The
width depends only on the coupling map: it is 6 for every model they ran and for the full 156-qubit chip with any
angles and any number of layers (`results/full_chip.txt`). Exact samples and exact probabilities therefore cost
milliseconds on a laptop (Python and NumPy, one process):

| qubits | trained models with hardware shots | width | exact sampler, ms per sample |
|---|---|---|---|
| 28 | 18 | 4 | 0.7 |
| 45 | 18 | 6 | 3.7 |
| 91 | 18 | 6 | 7.6 |
| 153 | 18 | 6 | 8.3 |
| 156, full chip, random angles, 3 layers | n/a | 6 | 6.0 |

So the instances run on hardware are not hard to sample classically, and neither is any model of this ansatz on
this chip. This does not contradict IQP hardness results, which concern circuits whose interaction graphs are dense
or wide (and anticoncentrating outputs); a nearest-neighbour IQP circuit on a graph of bounded width is a standard
classically easy case. The trained distributions are also far from anticoncentrated: exact samples have mean
log₂ p of −21.2, −32.8, −60.5 and −92.4 at 28, 45, 91 and 153 qubits (means over the 18 models), against −28, −45, −91 and −153 for uniform.

## How it was checked

- `iqp_exact.py --selftest` (`results/selftest.txt`): probabilities of every string against two independent dense
  references (Walsh–Hadamard of the phase vector, and gate-by-gate statevector) on 40 random graphs and vertex orders
  up to 11 qubits, max |ΔP| 2.2·10⁻¹⁶; a control with one edge moved is detected (|ΔP| 0.10); 40,000 samples have
  TV 0.0289 to the exact law against 0.0283 for ideal i.i.d. draws.
- `convention_check.py` (`results/convention_check.txt`): the angle convention exp(iθX) and the bit order (column k is
  qubit k) are fixed by IQPopt's own simulated samples of six trained 28-qubit models: mean log₂ p −20.0 under
  exp(iθX) against −27 to −29 for θ/2 and −30 to −36 for 2θ, and −35.0 with the bits reversed; our exact samples of the
  same model give −20.02 ± 0.06 against IQPopt's −20.00 ± 0.11.
- `check153.py` (`results/check153.txt`): at 153 qubits, where no statevector exists, 153 single-qubit and 300
  two-qubit expectation values from 4,000 exact samples agree with IQPopt's training formula (Monte Carlo over 2^16
  z, a separate code path, itself checked against the statevector to 9·10⁻¹⁶): pull sd 0.948, 0.935 and 1.016 over three
  seeds, mean pulls +0.34, −0.15, +0.12 (correlated through shared draws). Control: scaling the largest coupler's
  angle by 1.5 in the sampler is detected at 7.8σ and 7.3σ on its two qubits.

## The hardware shots, scored exactly (`score_hardware.py`, `summarize.py`, `results/score_hardware.txt`)

Every one of the 72 hardware sample files (512 shots each) was scored under its trained model. Because these laws
are strongly peaked, linear XEB is dominated by rare high-probability strings (values from 0.000 to 649 at 153
qubits) and the log-likelihood gap can be negative (noise that drives bits to 0 can make shots of a sparse-graph
model more likely than exact samples). Neither is a fidelity here; they are reported as measured:

| qubits | log-likelihood gap of the hardware shots, bits per shot: mean [min, max] over 18 models | median linear XEB |
|---|---|---|
| 28 | 0.67 [−1.06, 2.42] | 0.92 |
| 45 | 1.67 [−1.67, 3.71] | 0.70 |
| 91 | 3.66 [−1.55, 10.12] | 0.81 |
| 153 | 7.25 [−4.48, 20.14] | 0.41 |

The exact sampler reproduces each trained model with no gap by construction, in milliseconds, so for these models
classical sampling is both cheaper and more faithful than the hardware.

## The paper's own graph metrics, hardware against exact classical samples (`graph_metrics.py`)

`graph_metrics.py` reimplements the metrics of `src/analysis_scripts/analyze_nisq_results.py` (edge probability,
percentage of bipartite graphs, TVD of the degree distribution from Binomial(n − 1, p_target)). On the hardware
shots it reproduces all 72 rows of the paper's `nisq_analysis_results.csv` to the printed precision (max differences
0.0005, 0.05 %, 0.0005; `results/paper_nisq_analysis_results.csv`). The same metrics on 512 exact classical samples
of each trained model (`results/graph_metrics.txt`), averaged over models:

| qubits | data | target p | hardware: p, bipartite %, degree TVD | exact classical: p, bipartite %, degree TVD |
|---|---|---|---|---|
| 28 | Bipartite | 0.275 | 0.259, 50.8, 0.071 | 0.266, 49.0, 0.070 |
| 28 | ER | 0.475 | 0.451, 23.0, 0.057 | 0.476, 19.6, 0.018 |
| 45 | Bipartite | 0.252 | 0.302, 17.2, 0.137 | 0.308, 17.9, 0.151 |
| 45 | ER | 0.473 | 0.453, 15.7, 0.057 | 0.474, 14.6, 0.013 |
| 91 | Bipartite | 0.228 | 0.201, 20.1, 0.214 | 0.204, 23.2, 0.215 |
| 91 | ER | 0.474 | 0.446, 7.9, 0.097 | 0.474, 7.5, 0.018 |
| 153 | Bipartite | 0.214 | 0.153, 29.8, 0.255 | 0.154, 35.9, 0.260 |
| 153 | ER | 0.473 | 0.438, 2.9, 0.137 | 0.473, 1.3, 0.014 |

Two readings follow. For the Erdős–Rényi models the trained model itself hits the target density and degree law at
every size (TVD 0.013–0.018), and the hardware's shortfall (density 0.438 against 0.473 at 153 qubits, TVD up to 0.137)
is device noise. For the bipartite models the exact samples are no closer to the target than the hardware: the loss of
bipartite structure with size (100 % in the data, about 20–36 % in samples) is already in the trained model, so the
paper's "global structural features exhibit significant degradation beyond 91 qubits" is, for this family, a
limitation of the ansatz and its training rather than of the hardware.

## Not done

No noise model is fitted, the graph-level metrics of the paper (edge density, bipartiteness) are not recomputed,
and nothing has been sent to the authors or submitted.
