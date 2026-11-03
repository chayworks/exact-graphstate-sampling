# Cloud run 11, first attempt (doped instance, two-table engine gsamp.cpp, two-vertex slicing)

r7iz.8xlarge spot (Intel Xeon Gold 6455B, 32 vCPU, 247 GiB), us-east-1b, booted 2026-11-02 00:23 UTC, reclaimed by
spot 00:56 UTC, six minutes after the gate. Inputs (sha256 in `inputs.sha256`): the four width-33 two-table slice plans
doped_m20g_s{00,01,10,11} (slice vertices 135, 104; `sliced_plans_gsamp.py`), `doped_sampler.py --slices --K 518 --bin ./gsamp`.

- `gate.txt`: IBM doped shot #1 = -68.295331101 in **1,591.8 s** (one sliced evaluation = 4 slice passes of two 2^33
  tables; `time_gate.txt`: 26:32 wall, peak RSS 134,262,440 kB). Runs 7/9/10 gave -68.295330880 / -68.295331925,
  -68.295330954 / -68.295331107; run 9 needed 4,919 s on the same machine type with the in-place engine: 3.1x per shot.
- No IBM shots were scored before the termination (`run.log`). The run was relaunched with the same inputs (run 11b).
