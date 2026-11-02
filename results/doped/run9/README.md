# Cloud run 9 (doped instance, in-place engine with huge pages)

r7iz.8xlarge spot (Intel Xeon Gold 6455B, 32 vCPU, 247 GiB), us-east-1b, booted 2026-11-01 18:33 UTC, reclaimed by
spot ("Server.SpotInstanceTermination") at 23:50 UTC. Inputs: the m20 and m10 in-place plans of width 34 (sha256 in
`inputs.sha256`), `doped_sampler.py`, `gsamp2.cpp` (2 MB-aligned arena, MADV_HUGEPAGE).

- `gate_m20.txt`, `gate_m10.txt`: IBM doped shot #1 through the two orders: -68.295331925 and -68.295330954
  (`gate_result.txt`: PASS); 4,919 s for the m20 pass (`time_gate.txt`: peak RSS 134,227,596 kB).
- `shots_scored.txt`: the first string of the score list (a uniform random control, not an IBM shot):
  -72.910036757, 4,936 s. Run 10 scored the same string on an r7a.4xlarge with the sliced plans: -72.910036839.
- No samples: the sampling phase was never reached.
- `run.log`: the node's log; `dnf.log`: package installation.
