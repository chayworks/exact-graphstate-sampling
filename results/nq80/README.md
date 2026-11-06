# Cloud run 12: hybrid-sampler samples of the 80-qubit, depth-80 instance on one node

(The output law is proportional to min(z, 1.03) on the prefix times the exact suffix, equal to the Born law if no
prefix has z > 1.03; the 80 proposals certify fidelity ≥ 0.887 at 95% (`results/certificates_nq80_fstate.txt`), but the
mass above the cap is not bounded rigorously, so these are not labelled exact samples. See the note, §3.4.)

Instance: seed 1080 of the construction (`cloud_nq80/fi80_order.json` order, width 35, 1,573 edges, degrees 30–50).
Node: on-demand r7a.8xlarge (32 vCPU, 256 GiB), us-east-1a, launch template version 4 = `cloud_nq80/run_hybrid.sh`;
booted 23:13 UTC 2026-11-01, finished and self-terminated 04:59 UTC 2026-11-02. Inputs and their SHA-256 are in
`inputs.sha256` (the same files as `cloud_nq80/`).

Method: hybrid sampler (`doped_sampler.py sample ... --nd 80 --na 0 --m 20`, cap 1.03, seed 20261102): 60 bits proposed
uniformly and accepted on their exact one-pass marginal, last 20 bits drawn exactly. Vertex 34 is sliced (plans
`fi80h_s0/s1.txt`, width 34, K = 60) and its phase on suffix neighbours is applied as a translation of the dual table
at the summation (`fi80h_twist.txt`: `0 0` / `40146 35739773`; `gsamp2 resume --twist`). Two passes of 2^34 complex64
per proposal.

Files (all written by the node):

- `gate70_out.txt`, `gate70_result.txt`: the two known 70-qubit values, PASS (−69.508529982592, −68.876164845958 against
  −69.508530214929, −68.876165103181; 2.3e-7 and 2.6e-7).
- `gate_tw_hybrid.txt`, `gate_tw_resume.txt`, `gate_tw_result.txt`, `tw70_props.txt`: unsliced hybrid vs twisted two-slice
  resume on three 70-qubit proposals (seed 17, K = 50), PASS: identical sampled bits, log₂ Pr differences 8.2e-7, 1e-8, 2.8e-7.
- `samples_nq80.txt`: 80 proposals (79 `ACC`, 1 `rej`), each with log₂ Pr, z and the time of its two passes.
- `xcheck_strings.txt`, `xcheck_nq80.txt`: the first two accepted samples re-scored with the (10, 40) pair slicing
  (plans `fi80p_s00..s11.txt`, four passes, `fi80p_twist.txt`), 475 s and 470 s.
- `run.log`: the driver's timeline; `time_sample.txt`: GNU time of the sampling step (budget 19,700 s reached, exit 124,
  peak resident set 134,222,824 kB).
- `analysis.txt`: output of `python nq80_analyze.py results/nq80/samples_nq80.txt results/nq80/xcheck_nq80.txt --n 80`.

Result (from `analysis.txt`): 79 of 80 accepted (98.8%), none above the cap; z mean 0.99983, sd 0.00096, min 0.9970,
max 1.0018; linear XEB of the samples against their own exact probabilities +1.020 ± 0.138; mean log₂ Pr −79.293 ± 0.118
(min −82.888, max −77.304); 245.9 s per proposal (238.6–274.8); the two re-scored samples differ from the sampler by
5.2e-7 and 4.4e-7 in log₂ Pr. All 79 strings are distinct and 80 bits long.
