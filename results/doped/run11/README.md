# Cloud run 11 (second attempt, on-demand): IBM's doped shots scored and the first hybrid samples

Node: r7a.8xlarge on-demand (32 vCPU AMD EPYC 9R14, 246 GB), us-east-1a, launch template version 3 =
`cloud_doped/run_gsamp.sh`; booted 14:58:52 UTC 2026-11-02, finished and self-terminated 00:50:53 UTC 2026-11-03.
Engine: two-table sampler `gsamp.cpp` on the four width-33 slice plans (`cloud_doped/doped_m20g_s00..s11.txt`, slices
at ancilla vertices 135 and 104, K = 518); inputs and their SHA-256 in `inputs.sha256`. The first attempt of this run
(spot, reclaimed after its gate) is archived in `attempt1_spot/`.

Timeline (`run.log`):
- 15:18:44 gate: shot 1 of IBM's 2,051 doped shots re-scored at −68.295330552 (reference −68.295331 from runs 7, 9, 10):
  PASS; 1,168 s for one sliced evaluation (4 slice passes), peak resident set 134,262,660 kB (`gate.txt`,
  `gate_result.txt`, `time_gate.txt`).
- 15:18 to 20:37 scoring: IBM's shots taken from the end of the published list (shot 2051, 2050, ...; `score_list.txt`,
  two uniform controls last and not reached): 16 shots in 19,145 s, 1,165 s each (`shots_scored.txt`, lines
  `<Qiskit string> <log2 Pr> <seconds>`).
- 20:37 to 00:50 sampling: hybrid sampler, seed 20261101, cap 1.03, m = 20 (50 data bits proposed, 20 drawn exactly):
  13 proposals in 15,183 s, 13 accepted, z in [0.9981, 1.0007], 1,166 s per proposal (`samples.txt`, lines
  `ACC|rej <70 data bits, qubit order> <log2 Pr_doped> z=<...> <seconds>`).

Analyses:
- `tn242_compare.txt` (`python doped_tn242_compare.py results/doped/tn242_log2p.txt results/doped/run11/shots_scored.txt`):
  all 16 shots match strings in #242's Zenodo release; log₂ difference ours − TN: mean −2.56e-4, max |.| 2.76e-4; 15 of 16
  within [−2.76, −2.54]e-4 and the least probable shot (log₂ p −74.73) at −1.63e-4, i.e. a common factor of 1 − 1.8e-4
  between the two pipelines (origin unidentified) plus rounding that is largest on the smallest amplitude.
- `analysis.txt` (`python doped_analyze.py results/doped/run11 0`): the 16 shots' linear XEB +0.332 ± 0.253
  (bootstrap 95% [−0.136, +0.829]); hybrid sampler 13 of 13 accepted, z mean 0.9996, sd 0.0008, none above the cap;
  accepted samples' self-consistency XEB +0.807 ± 0.358; mean time per proposal 1,166 s.
- `porter_thomas_samples.txt` (`python porter_thomas_test.py results/doped/run11/samples.txt 70 --acc`): mean z
  1.81 ± 0.36 (Porter–Thomas 2), mean log₂ Pr −69.59 ± 0.36 (Porter–Thomas −69.39), KS 0.14 (5% level 0.38).

The samples are the Born law conditioned on z ≤ 1.03 (hybrid sampler), not exact samples; see the note, §7.
