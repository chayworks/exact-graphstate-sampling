# Cloud run 10 (doped instance, in-place engine, two-vertex slicing)

r7a.4xlarge on-demand (AMD EPYC 9R14, 16 vCPU, 123 GiB), booted 2026-11-01 19:40 UTC, terminated by us at 23:05 UTC
to give its quota to run 12. Inputs: the width-33 slice plans doped_m20sl_s{00,01,10,11} (slice vertices 135 and 104;
sha256 in `inputs.sha256`), `doped_sampler.py --slices --K 518`, `gsamp2.cpp`.

- `gate.txt`: IBM doped shot #1 = -68.295331107 (run 7: -68.295330880; `gate_result.txt`: PASS); 5,888 s for one
  sliced evaluation = 4 slice passes of a 2^33 table (`time_gate.txt`: peak RSS 67,118,812 kB).
- `shots_scored.txt`: the first string of the score list (a uniform random control): -72.910036839 in 5,938 s;
  run 9 scored the same string unsliced on an r7iz.8xlarge: -72.910036757.
- `run.log`: abridged (the progress lines were identical); `time_gate.txt`: abridged to the informative fields.
