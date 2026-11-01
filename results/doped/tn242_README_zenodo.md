# Data for exact tensor-network simulation of IBM's doped-Clifford circuit

This is the dataset for *Classical Simulation and Design Frontiers for IBM's Doped Clifford Sampling Experiment*.
It contains the amplitudes used to calculate log-XEB and the contraction paths for the circuit.

## Contents

- `data/amplitude_batches.npz`: 2,051 batches of 256 complex64 amplitudes.
- `analysis/log_xeb.json`: the log-XEB result and estimator definition.
- `path/constructive/` and `path/bm4/`: the constructive path and its BM4 variant used for GPU execution. Each directory contains `compiled_program.npz`, `step_info.csv`, and `metadata.json`.
- `circuit/nq70_depth70_checks27_doped.qasm`: the unmodified public circuit from the IBM dataset.
- `scripts/verify.py`: a verification script that recomputes log-XEB and verifies the contraction widths and costs.

## Verify the release

Python 3.10 or later and NumPy are sufficient:

```bash
python -m pip install -r requirements.txt
python scripts/verify.py
```

## Amplitude batches

`data/amplitude_batches.npz` contains the following arrays:

| member | shape | dtype | meaning |
|---|---:|---|---|
| `raw_vectors` | `(2051, 256)` | `<c8` | Raw production output vectors |
| `row_indices` | `(2051,)` | `<i8` | Zero-based IBM row indices |
| `bitstrings_q0_first` | `(2051,)` | Unicode | Requested 70-bit strings |
| `assignments_q8_q69` | `(2051,)` | Unicode | Fixed values on `q8,...,q69` |
| `selected_vector_indices` | `(2051,)` | `uint16` | Entry matching each full bitstring |
| `recovery_factor` | scalar | `<f8` | Scale converting raw to physical amplitudes |
| `metadata_json` | scalar | Unicode | Schema and convention metadata |

The calculation leaves `q0,...,q7` open.
The vector entry corresponding to a requested bitstring is therefore

```python
index = int(bitstring_q0_first[:8], 2)
amplitude = raw_vectors[row, index] * recovery_factor
probability = abs(amplitude) ** 2
```

## Contraction paths

The constructive path follows Algorithm 1 of the paper, with `q0,...,q7` left open instead of the lower qubits used there.
The BM4 path is a branch-merge transform of the constructive path designed to accelerate GPU execution.

| path | role | width | cost (complex FMAs) |
|---|---|---:|---:|
| Constructive | Base path before execution blocking | 35 | 303,946,252,222,144 |
| BM4 | Blocked GPU production schedule | 35 | 566,903,624,260,592 |

Both `compiled_program.npz` files contain the same 4,830 complex64 input tensors and projector locations, with different sequences of 4,829 binary contractions.
`step_info.csv` lists each contraction, and `metadata.json` records hashes, widths, and costs.

## AI-assisted development

AI coding agents, including OpenAI Codex, were used to assist with software development, documentation, and construction of verification checks.
All AI-assisted outputs included in this release were reviewed by the author, who remains responsible for the released data, software, and scientific claims.

## Third-party data

The circuit and IBM bitstrings originate from Chung *et al.*, *Dataset for Sampling hard circuits with verifiably high fidelity*, Zenodo [10.5281/zenodo.21633064](https://doi.org/10.5281/zenodo.21633064), licensed under CC BY 4.0.
See `ATTRIBUTION.md` for checksums and full attribution.

## License

Generated data, contraction programs, analysis records, metadata, and documentation are licensed under CC BY 4.0.
The verification script is licensed under the MIT License.
The redistributed IBM circuit and bitstrings remain under their original CC BY 4.0 license.
See `LICENSE.md` for the exact scope and `ATTRIBUTION.md` for third-party attribution.
