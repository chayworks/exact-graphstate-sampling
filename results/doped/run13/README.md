# Cloud run 13: hybrid samples of IBM's doped 70-qubit instance (#228)

Node: r7a.8xlarge on-demand (32 vCPU AMD EPYC 9R14, 246 GB), us-east-1a, launch template version 5 =
`cloud_doped/run_sample.sh`; booted 01:16:58 UTC 2026-11-03. Stopped by us at 19:42 UTC 2026-11-04, about two hours
before its 44-hour budget ended, to free the node for run 14; the sampler uploads `samples.txt` every five minutes,
and the last upload (19:39 UTC) contains all 126 completed proposals. Because the run was stopped, `time_sample.txt`
is empty and `run.log` has no "sampling ended" line. Inputs and SHA-256 in `inputs.sha256` (same files as run 11).

- Gate (`gate.txt`, `gate_result.txt`, `time_gate.txt`): IBM's doped shot 1 re-scored at −68.295331011 against the
  reference −68.295331 of runs 7, 9, 10, 11: PASS; 1,207.6 s for one sliced evaluation (four slice passes), peak
  resident set 134,262,424 kB.
- Sampling (`samples.txt`): hybrid sampler, seed 20261103, cap 1.03, m = 20 (50 data bits proposed, 20 drawn exactly):
  126 proposals, 122 accepted, none above the cap; prefix marginal z mean 1.00013, sd 0.00092, range
  [0.99783, 1.00269]; 1,194 s per proposal; all 122 strings distinct.
- `porter_thomas_samples.txt` (`python porter_thomas_test.py results/doped/run13/samples.txt 70 --acc`): mean z
  2.181 ± 0.136 (Porter–Thomas 2), mean z² 6.99 ± 0.90 (6), mean log₂ Pr −69.245 ± 0.102 (−69.390), KS distance to
  z·e^{−z} 0.108 (5% level 0.123); self-consistency XEB 1.181 ± 0.136.

Runs 11 and 13 together: 139 proposals, 135 accepted samples, all distinct; z in [0.99783, 1.00269], sd 0.00093;
self-consistency XEB 1.145 ± 0.127; mean log₂ Pr −69.278 ± 0.098. The samples follow the capped
law, proportional to min(z, 1.03) on the prefix (hybrid sampler), not exact samples; with run 11 the 139 proposals
certify fidelity ≥ 0.936 at 95% (`results/doped/certificate_run11_13.txt`); see the note, §7.
