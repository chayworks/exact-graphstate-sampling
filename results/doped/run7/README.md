# Cloud run 7 (2026-11-01): first exact full-depth probability of an IBM doped-circuit hardware shot

Node: spot r5a.8xlarge (32 vCPU AMD EPYC 7571, 249 GB), us-east-1d, after a first r7a.8xlarge spot was reclaimed at 16 min.
In-place sampler, width 34 (table 2^34 complex64; peak RSS 134,226,760 KiB = 137 GB).

Gate (one shot through two plans with different step structure): PASS
- IBM doped shot #1 (`1010100101110001...`): log2 Pr = -68.295330880 (m20 plan) vs -68.295331335 (m10 plan): agree to 4.6e-7.
  That is z = 2^70 Pr = 3.26 times the uniform value for this shot.
- One pass took 9,872 s (2 h 44 min wall, 3015% CPU) on this node: about 2e8 table updates per second, 8x slower than
  the r7a rate measured at 76 qubits. The scoring phase (timeout 2.5 h) therefore finished no further string and the
  sampling phase none either; the run ended at its budget. Cost about $6.

Files: `gate_m20.txt`, `gate_m10.txt`, `gate_result.txt`, `time_gate.txt`, `inputs.sha256`, `run.log`.
Relaunched as run 8 on an r7iz.8xlarge spot (newer cores) with the same inputs.
