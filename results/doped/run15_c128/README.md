# Cloud run 15: two doped amplitudes in complex128, to locate the 1 − 1.8·10⁻⁴ factor against #242

Node: r7a.8xlarge on-demand (32 vCPU AMD EPYC 9R14, 246 GB), launch template version 7 = `cloud_doped/run_c128.sh`,
booted 04:18:55 UTC 2026-11-06, finished 06:13:52 UTC and terminated itself. gsamp2 built with `-DGS_DOUBLE`
(complex128 tables and dumps; locally `test_gsamp2_doped.py` gives 1.2e-13 against dense statevectors, 8.6e-9 in
complex64, and the sliced driver `results/tests/test_slice_adjacent_c128.txt` 5e-10, 6.7e-6 in complex64). Run 10's
two-slice plans (ancillas 135, 104; four plans of width 33; one 2^33 complex128 table, peak 134,227,620 kB).
Inputs and SHA-256 in `inputs.sha256`.

| IBM doped shot | complex128 (this run) | our complex64 | #242 (tensor network) | c128 − c64 | c128 − #242 |
|---|---|---|---|---|---|
| first (shot 1) | −68.295331132 | −68.295330880 (run 7; seven runs within 1.4e-6) | −68.295074 | −2.5e-7 | −2.57e-4 |
| last (shot 2,051) | −69.167144551 | −69.167143120 (run 11) | −69.166885 | −1.4e-6 | −2.60e-4 |

(log₂ Pr; 3,441 s and 3,429 s per amplitude.) Double precision moves our values by at most 1.4e-6, the complex64
rounding level, and leaves the offset from #242 unchanged: the common factor 1 − 1.8·10⁻⁴ is not complex64 rounding in
this pipeline. Its source is on #242's side or in a convention between the two pipelines; this run cannot say which.
Independently, the same pipeline reproduces Renaud's #249 probabilities of the ring instance to 3.8e-6 (run 14).
