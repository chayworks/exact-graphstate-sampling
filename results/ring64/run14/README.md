# Cloud run 14: hybrid samples of IBM's 64-qubit ring instance (#245)

Circuit `dcs_nq64_depth73_checks12_t314.qasm` (MD5 054441020510c19c628284de146f5f34, identical to IBM's
`v2/loop64_314t.qasm`); 314 T gates as gadgets, 378-vertex graph. Plan: Renaud's published order (#249) with the last
20 data vertices moved to the end (width 36), slice vertices 317 and 211 (four plans, width 34, no twist),
`cloud_ring64/`. Node: r7a.8xlarge on-demand (32 vCPU AMD EPYC 9R14, 246 GB), launch template version 6, booted
21:03 UTC 2026-11-04, terminated itself after 20:53 UTC 2026-11-05. Inputs and SHA-256 in `inputs.sha256`.
A first attempt failed its gate and exposed the slice-edge sign bug fixed in commit f2e0bf7.

- Gate (`gate.txt`, `gate_result.txt`, `time_gate.txt`): IBM hardware shots 0-2 scored against Renaud's exact
  probabilities (Zenodo 22882048): ours -67.331000392, -64.672098827, -64.550927311; differences -3.8e-6, -1.9e-6,
  -2.1e-7 in log2: PASS; 2,425 s for three evaluations, peak resident set 134,226,068 kB.
- Sampling (`samples.txt`, `time_sample.txt`): seed 20261104, cap 1.03, m = 20 (44 bits proposed, 20 drawn exactly);
  102 proposals, 102 accepted, none above the cap, max z 1.002356; 83,350 s (about 817 s per sample); all 102
  strings distinct. Strings are indexed by qubit (character k = qubit k), the reverse of Qiskit order.
- `porter_thomas_samples.txt`: mean 2^64 Pr 1.910 +/- 0.141 (Porter-Thomas 2), mean square 5.65 +/- 0.91 (6),
  KS distance to z e^-z 0.086 (5% level 0.135).
- `xeb_vs_hardware.txt` (`ring64_xeb_compare.py`): linear XEB 0.910 +/- 0.141 for these samples against
  0.390 +/- 0.034 for IBM's 1,389 hardware shots (probabilities from #249), difference 3.6 sigma.
- `certificate.txt` (`hybrid_certificate.py`): from the 102 uniform proposals alone, without a Porter-Thomas
  assumption, TV(sampled, ideal) <= 0.088 and classical fidelity >= 0.912 at 95% confidence (empirical Bernstein on
  min(z, c); floating-point error not included).

Context: tracker #263 (posted 06:41 UTC 2026-11-05, while this run was sampling) drew 1,389 samples of the same
instance in 864 s on 8 H200 GPUs with a uniform prefix and no acceptance step, certified fidelity >= 0.908. This run
is an independent CPU implementation of the same family; it is not the first sampling of #245.
