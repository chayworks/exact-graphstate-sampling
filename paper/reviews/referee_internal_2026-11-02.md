# Internal referee report, 2026-11-02

Independent review of `paper/exact_graphstate_sampling.md` (draft of 31 Oct 2026, 622 lines) by a reviewer with no
knowledge of the project history, briefed to check every number against `graphsamp/results/` and to referee the
mathematics, statistics, novelty and presentation. Read-only; one script recomputed statistics from committed sample
files. Kept verbatim; the author's responses and the resulting changes are tracked in `paper/reviews/response_2026-11-02.md`.

---

## Summary

The note claims that IBM/UChicago's 70-qubit, depth-70 graph-state sampling instance (tracker #151) is sampled exactly from the Born distribution by a linear-order cut-rank dynamic program at 13 s/sample on 16 cores and 18 s on a laptop, that IBM's 1,726 hardware shots can be scored exactly (XEB 0.466 ± 0.031, fidelity 0.46–0.48), and that 76- and 80-qubit instances of the same construction were sampled exactly on one cloud node, the latter via a "dual-table translation" slicing. It further asserts a universal ⌊n/2⌋ cut-rank ceiling for any graph state in any product basis, that the instance's output is exactly uniform on every ≤14-qubit marginal (minimum stabilizer weight 15, proven), and, in §7, that the doped sibling instance reduces to a 538-vertex graph state of width 34 whose first shot probability agrees with the tensor-network entry #242 to 1.8·10⁻⁴. The artifact trail is unusually good and nearly every number I checked matches; the problems are in what the numbers are claimed to mean.

## Major issues

**M1. The 80-qubit samples are not "exact" in the sense the note uses elsewhere; they are capped rejection samples, the very method the note distinguishes itself from.**
Lines 30–32 (abstract), 69 (results table), 222–226, 579 say "79 exact samples". Line 116–118 says: *"Unlike rejection sampling, the cost does not depend on the shape of the output distribution, and unlike capped rejection sampling (#260) the output law is the Born law."* But §3.4 and §7 (lines 500–508) describe the 80-qubit sampler as: propose 60 bits uniformly, accept with probability z/c, c = 1.03, with the output *"exactly except for the mass of prefixes with z > c, which is ... under the Porter–Thomas model, ... below 10⁻³⁰."* That is capped rejection sampling. The improvement over #260 is real (z averages 2²⁰ probabilities, so the observed spread is 10⁻³ rather than O(1)), but the exactness claim rests on an unproven tail model for a distribution the note itself says is *not* Porter–Thomas (line 545: "unlike the Clifford instance the doped circuit's output is Porter–Thomas"). From the artifact (`results/nq80/analysis.txt`: 80 proposals, max z 1.0018) the only data-driven statement is that the cap-exceedance probability is below roughly 4% (rule of three), which bounds the total-variation error at the percent level, not 10⁻³⁰.
What would fix it: Var_uniform(z) = 2^K Σ_b Pr(b_P)² − 1 = Σ_{S≠I supported in P} 3^{−|S|} is exactly computable by enumerating the 2^{K−r} = 2⁴⁰ stabilizers inside the proposed prefix (the note already enumerates 2.6·10¹² vectors in `minweight.cpp`); Chebyshev then gives a rigorous bound on P(z > c) and on the truncated mass, which would turn "exact under a model" into "exact to a proven TV ≤ ε". Alternatively, retitle the 80-qubit result as what it is and fix line 118. Note also that no independent engine checked any 80-qubit probability: `results/nq80/xcheck_nq80.txt` uses the same `gsamp2` code with a different slicing (the abstract's "unrelated slicing" is accurate but should say "same engine").

**M2. Internal contradiction on the doped instance.**
Line 397–400 (§5 caveats): *"the doped sibling instance on the tracker, with 468 T gates spread through the circuit, is outside its reach (its amplitudes were computed by tensor-network methods ...)."* Line 573 (§8, "It does not show"): *"anything about the doped instance."* But §7 (lines 477–551) shows precisely that the doped instance is a 538-vertex graph state in a product basis with a width-34 order, computes its first shot's probability five ways on single nodes (`results/doped/run7..run11`), matches #242 to 1.8·10⁻⁴ (`tn242_crosscheck_run7.txt`), and re-analyses all 2,051 #242 amplitudes. Both stale sentences must go, and the abstract should state the §7 result, which is arguably the strongest technical content in the note (see "underclaimed" below).

**M3. "Every number in this note is read from a committed artifact" (line 7) is false for several numbers.**
I searched `results/` with `grep --text` for every figure quoted from the toy and slicing tests and found none of:
- line 216–219: `test_dual_table.py` 2·10⁻⁸; `test_twist.py` 9·10⁻⁸, control "off by 2.5", 1.7·10⁻⁶, 3.5·10⁻⁷ (the scripts print these at runtime; no output is committed);
- line 395–396: *"Slicing trades memory for time at a measured cost of about 1.5× per sliced vertex on this instance, not 2×"* — no artifact, no pointer, and the only sliced-vs-unsliced timings in the repository (run 9: 4,919 s unsliced on r7iz; run 10: 5,888 s sliced on a different, 16-core r7a) are across machines and do not support a per-vertex factor;
- line 483–485: doped-gadget checks on "20 random circuits (1e-10)", "25 more (1e-15)", "12 (1e-8)" — test scripts only;
- line 235–237: small-instance 2·10⁻¹⁵ and the 20,000-shot histogram — in `graphsamp/README.md`'s test table, not in `results/`.
Either commit the test outputs or weaken line 7.

**M4. The low-order-marginal null results on the hardware shots are over-interpreted.**
Lines 46–48, 458–461 and the Figure 4 title: *"The errors that bring the fidelity to about 0.46 therefore leave no low-order trace at all: no biased qubit, no correlated pair ... the residual noise acts on the output like a global depolarisation."* With N = 1,726, the standard error of a single-bit mean is 0.024, so a readout asymmetry or single-qubit bias of 2–3% (z < 1.3) is invisible, and a fidelity of 0.46 over 70 qubits corresponds to ~1% error per qubit. The note's own power demonstration (line 293–294, `hw_marginals.txt`) uses a 15% bias. The honest statement is "no bias above ~7% on any qubit or ~7% on any parity is detectable at this N", which is what one would expect under almost any noise model and is not evidence for global depolarisation. The sentence in the abstract ("whatever brings the fidelity to 0.46 leaves no trace below weight 5") should be qualified the same way.

**M5. The frontier/tightness claims are contradicted by the note's own longer searches.**
Line 355: *"at every n it is between n/2 − 4 and n/2 − 3"*; line 37–38 and 375: "within 3–4 of this ceiling". But: the real 70-qubit instance reaches width 30 (gap 5, line 124–125); the sampled 76-qubit instance has width 33 (gap 5, `fi76_s2.log`); the sampled 80-qubit instance has width 35 (gap 5, `fi80_s5.log`), while `frontier.jsonl` lists 36 for the same 1,573-edge instance because of a shorter search. The gap grows with search effort, so "tight to within the heuristic's few bits" (line 392–393) is a statement about the heuristic's budget. State the gap as 3–5 at the budgets used and say the true linear rank-width is unknown. Also the §5 table row for n = 80 (line 345) quotes 2^{41.2} updates, which is `frontier.jsonl`'s width-36 order (log2_work 41.21), not the width-35 order that was sampled (`fi80_s5.log`: log2 Σ2^{r_k} = 38.29, i.e. 2^{40.3}).

**M6. "One large-memory node below about 76 qubits" (abstract line 36, lines 372–375) depends on an unstated node size and is inconsistent with Figure 1.**
Line 372 calls 2.2 TB "the largest single cloud nodes"; Figure 1 labels the node band "up to 512 GB"; AWS u-24tb1 / u7in-32tb nodes hold 24–32 TiB, at which the ceiling fits to n ≈ 82–84 and the measured construction width (0.49n − 2.9) fits to n ≈ 90. The thresholds "76 / 80 / 90" should be stated as functions of an explicit memory figure, and the abstract's "or a supercomputer below about 90" revised accordingly.

**M7. Fidelity comparison framing.** Line 24–26: estimates "above IBM's direct-fidelity estimate of 0.37 ... and much narrower". The note later (line 184–187) concedes the two quantities differ and that linear XEB over-estimates fidelity under non-depolarising noise. The abstract should carry that caveat, and the mixture-MLE model (density (Fz + 1 − F)e^{−z}) assumes the ideal z over strings is Exp(1); line 136 says the instance is only "Porter–Thomas-like" (2⁷⁰Σp² = 2.02) and line 545 says it is not Porter–Thomas. No test of the 1,012 samples' z against z·e^{−z} is reported. Either show that test or state the model dependence of F_MLE.

## Minor issues

1. Line 530 vs 534: "The six values span 1.1·10⁻⁶" then "The sixth value comes from a different method altogether" — the table has six in-house values; the TN value is the seventh. (Span of the six is 1.05·10⁻⁶, fine.)
2. Line 494: "our 82 minutes per shot on one 32-core node" is stale; run 11 in the same table gives 1,592 s = 26.5 min (`run11/README.md`).
3. Line 492–495: "the two methods do the same amount of arithmetic ... the gap in wall-clock time ... is hardware, not algorithm" compares "table updates" to "complex FMAs" and uses #242's per-amplitude cost amortised over batches of 256. This is not a FLOP-normalised comparison; soften.
4. Line 214: "at the summation cut ... it always did" — the span condition is presented as luck; explain when it can fail (sliced column outside the column space of the remaining prefix) and that the fallback (extending the small summation-cut table by one rank) is cheap.
5. The 80-qubit section never says why exact chain-rule sampling was not used (unlike §7, which explains it for ancillas). The reason appears to be that the two slice tables cannot coexist in 256 GiB at every step; say so.
6. Line 225: "mean log₂ Pr −79.29 ± 0.12 against −80 for the uniform distribution": the relevant comparison is the Porter–Thomas value −80 + 0.61 = −79.39, which it matches.
7. Line 292–315 (single-bit χ² p = 0.010 on the sampler's own output): recorded as a flag; a second 1,012-sample run with a different seed base costs ~$5 and would settle it. Seeding `std::mt19937_64` with consecutive integers is also poor practice and worth a sentence.
8. Line 156–158: the uniform control has 64 strings (+0.16 ± 0.12, 1.3σ from zero). At 6.9 s/string, 1,000 controls cost two hours; the control should be at least as powerful as the claim it supports.
9. Notation: F denotes the unmeasured suffix in §2 but the measured set in §4 item 5 and §6.
10. Line 50–56 credits "the observation that the post-measurement states Z^v|G⟩ are orthonormal" as the joining identity; this is the Schmidt decomposition of a graph state across a cut (Hein, Eisert, Briegel, PRA 69, 062311 (2004)), and should be cited as such.
11. Line 213: the "dual table translation" is the Walsh–Hadamard shift theorem (a character multiplying F is a translation of its transform); name it.
12. Line 199 lists the 76-qubit engine gate values "−69.5085289 and −68.8761655"; `nq76/analysis.txt` has −69.50852886805679 and −68.8761654749801. Consistent, fine.
13. `results/nq76/xcheck_out.txt` contains two unlabeled number pairs (1.7897… 0.6554…; −0.7132… 1.0349…) that are not the log₂ Pr values quoted; `analysis.txt` reads them through `xcheck_make_run6.txt`. A reader following the appendix pointer will not be able to interpret the file. Label it.
14. §7 placeholder (line 553–554): known pending; not counted.

## Numerical checks performed

| Claim (line) | Artifact | Match |
|---|---|---|
| 1,012 samples, 13,202 s, 13.0 s each, peak 34.4 GB (62, 133) | `time_nq70.txt`: 3:40:01, RSS 33,558,388 kB; `samples_nq70.txt`: 1,012 lines, all distinct (recomputed) | yes |
| Self-XEB 1.020 ± 0.047 (136) | recomputed from `samples_nq70.txt`: 1.0200 ± 0.0471 | yes |
| Production order cost 6.0·10¹⁰ updates (126) | `cloud/nq70_plan.txt`: 4·Σ2^{rn} = 6.01·10¹⁰ (recomputed) | yes |
| Laptop 3 samples, 55 s, 8.6 GB (139–142) | `laptop_samples_w30.txt`: real 54.999 s; 2³⁰·8 B = 8.59 GB | yes |
| Independent engine 5 samples ≤10⁻⁶, control −72.8808316, 3.37 lower (242–251) | `xcheck70_independent.txt`, `xcheck70d.txt` | yes |
| Laptop sample −68.5019327 vs −68.5019326 (144) | `xcheck70e_laptop.txt` | yes |
| IBM XEB +0.466 ± 0.031, bootstrap [0.41, 0.53], controls +0.16 ± 0.12 / −0.15 ± 0.10 (156–158) | `ibmxeb/analysis.txt`; XEB recomputed from `ibm_out.txt`: 0.4660 ± 0.0311 | yes |
| vs #260: 9.4·10⁻⁶ max, 2.8·10⁻⁷ median (161) | `ibmxeb/analysis.txt` | yes (the #260 file itself: not checked) |
| 6.9 s/shot, 3 h 19 min, 17.2 GB (150–151) | `ibmxeb/time_ibm.txt`: 3:19:28, 16,781,328 kB | yes |
| Gate 5·10⁻⁷ (152) | `gate_result.txt` vs `samples_nq70.txt`: 4.6·10⁻⁷, 3.8·10⁻⁷ | yes |
| F_XEB 0.457 ± 0.037; F_MLE 0.48 [0.43, 0.54]; KS 0.018 vs 0.033 (165–169) | `ibmxeb/analysis.txt` | yes |
| Synthetic estimator table (170–172) | `fidelity_mle_synthetic.txt` | yes |
| 76 qubits: 1,465 edges, width 33/natural 38, 19 samples, 3,604 s, 69 GB, XEB 0.83 ± 0.22, dp check 10⁻⁷ (193–201) | `fi76_s2.log`, `nq76/run.log`, `time_sample.txt` (67,113,728 kB), `analysis.txt` | yes |
| 80 qubits: 1,573 edges, width 35, 275 GB; 79/80 accepted; z 0.9970–1.0018; 246 s; 128 GiB; 5.5 h; XEB 1.02 ± 0.14; mean log₂ Pr −79.29 ± 0.12; xcheck 5.2/4.4·10⁻⁷; gates 2.6·10⁻⁷ and 8.2·10⁻⁷ (205–230) | `fi80_s5.log`, `nq80/samples_nq80.txt` (recomputed: 79 ACC, 1 rej, 79 distinct, z 0.9970–1.0018, 245.9 s), `time_sample.txt` (134,222,824 kB, 5:28:20), `analysis.txt`, `gate70_result.txt`, `gate_tw_result.txt` | yes |
| 16 other seeds give widths 35–36 (230) | `fi80_cseed_scan.log` | yes |
| Toy twist/dual-table numbers (216–219) | none in `results/` | **no artifact** |
| 55-qubit bit-for-bit reproduction (252–254) | `check55_result.txt` | yes |
| Uniform-prefix check (257–261) | `uniform_prefix_check.txt` | yes |
| Stabilizer-only conditionals 3.3·10⁻⁸, deviations ≤8·10⁻⁵, steps 1–42 (267–269) | `chain/stabprefix_ibm_fwd.txt`, `stabprefix_rev_ideal.txt` | yes |
| Information profile numbers (276–281) | `chain/analysis_profile.txt` | yes |
| Cross-order 1.4·10⁻⁵ / 3.5·10⁻⁶ / 4.3·10⁻⁶ (284–286) | `chain/analysis_crossorder.txt` | yes |
| hw_marginals, rng_check, draw_check, rng_scan figures (292–315, 456–458) | respective files | yes |
| Per-step fidelities (323–325) | `chain/analysis_profile.txt` | yes |
| Frontier table and fit 0.492n − 2.9 (338–355) | `frontier.jsonl`; fit recomputed: 0.492n − 2.91 | yes, except n = 80 row's 2^{41.2} is the width-36 order (M5) |
| Ceiling memory 275 GB / 2.2 TB / 281 TB / 9 PB; 9.6·10¹² updates; 35 min–1.6 h (369–373) | arithmetic | yes |
| Ceiling-tightness table (384–389) | `ceiling_tightness.jsonl`, `.cmd.txt` | yes |
| "1.5× per sliced vertex" (396) | none | **no artifact** |
| d_min = 15, three-way split, 974,120 vectors, weight-16 case-2 minimum, lightest v (424–432) | `minweight_exact.txt`, `minweight_k11.txt`, `minweight_case1_k11.txt` | yes |
| Split test on n = 30, 40 (433–437) | `minweight_split_test.txt` | yes |
| d_min scan 6, 8, 9, 12 (449) | `minweight_scan.txt` | yes |
| Doped widths 137 / 34 / 2^{40.9} (487–491) | `doped/doped_widths.json` | yes |
| Doped shot-1 table, seconds per pass, 3.1× (522–530) | `run7/9/10/11 README, gate_result, time_gate` | yes |
| Control string −72.910036757 / −72.910036839 (531) | `run9/run10 shots_scored.txt` | yes |
| TN #242 −68.295074, rel. diff 1.8·10⁻⁴ (536) | `tn242_crosscheck_run7.txt` | yes (Zenodo file: not checked) |
| #242 analysis: 0.342 ± 0.028, 0.343 [0.293, 0.393], KS 0.014/0.030, siblings 0.9991 ± 0.0014, z² 2.0013, KS 0.0016/0.0019, local 0.339 ± 0.028, log-XEB 0.350 ± 0.027 (544–548) | `tn242_fidelity.txt`, `tn242_porter_thomas.txt`, `tn242_log_xeb.json` | yes |
| Small-circuit and gadget test tolerances (235–237, 483–485) | README test table / scripts only | **no committed output** |
| IBM README quotations, DFE 0.37/0.092, #242 run details, #260 reported width | external | not checked |

## Mathematics

- Prefix-marginal identity (lines 105–110): correct. Z^v (v ≠ 0) anticommutes with K_u for any u ∈ v, so ⟨G|Z^v|G⟩ = 0 and ⟨G|Z^v Z^w|G⟩ = δ_{vw}; the full state is 2^{−|P|/2}Σ_{x_P}|x_P⟩Z^{A[F,P]x_P}|G_F⟩, and c ↦ Bc is a bijection onto the column space. This is the graph-state Schmidt decomposition (Hein–Eisert–Briegel 2004) and should be cited rather than presented as an observation.
- ⌊n/2⌋ ceiling (359–367): trivially correct; "kernel of dimension at most 2" per step follows from the rank changing by at most ±1 when a vertex moves across the cut. The "at most four passes" constant is asserted, not shown, but is plausible from §2.
- T-gadget (477–482): ⟨φ|_a CNOT|ψ⟩|0⟩ = T|ψ⟩/√2 with |φ⟩ = T†|+⟩, orthogonal outcome ZT|ψ⟩/√2, each with probability exactly ½ independent of |ψ⟩, hence Pr(all succeed) = 2^{−468} and Pr_doped = 2^{468}Pr_gadget. Correct. "Apply T then H, keep 0" is correct (⟨0|HT = ⟨φ|).
- §6 proposition (409–422): correct. H_U embeds in Z₂×Z₂ via the letter at any j ∈ U; the Fourier coefficient at the character on U is a sum of 1 or 3 terms of magnitude 3^{−|U|/2}, hence nonzero. The general-basis caveat is correctly stated.
- Hybrid rejection sampler (500–508): the accept-with-z/c argument is correct *if* z ≤ c for every prefix. That premise is asserted via a Porter–Thomas model, not shown (M1). The formula z = 2^{50}·2^{468}·Pr(prefix, anc = 0) is right for n = 70; for n = 80 the prefix is 60 bits (K = 60 in `nq80/README.md`), which the note does not say.
- "The 80-qubit samples are exact": as a claim about the algorithm, it is exact-up-to-cap; as a claim verified at 80 qubits, the evidence is the sampler's self-consistency plus a same-engine re-slicing. The DP itself was independently verified only at 55, 70 and 76 qubits.

## Statistics

- XEB error bars are standard errors of the mean of z with bootstrap confirmation; correct. Error propagation into F_XEB (0.457 ± 0.037) checks.
- The mixture MLE assumes an Exp(1) ideal z-distribution; the synthetic test (which resamples the real 1,012 samples) shows a mild bias at F = 1 (0.98 [0.95, 1.00]) and 0.94 ± 0.05 for F_XEB, which is consistent with the model being slightly off; this should be said.
- "Porter–Thomas" is used inconsistently (lines 137, 545) with no test of the Clifford instance's z-distribution reported.
- The χ² parity tests are correctly constructed and the 15%-bias control passes, but the power is not stated (M4). The quadruple-parity Bonferroni remark (line 457) is appropriate.
- 79 samples support "consistent with Porter–Thomas statistics" (XEB 1.02 ± 0.14) and nothing stronger; they cannot detect a 10% TV error in the sampler. The note does not claim more than consistency here, but the abstract's "exact" is carried by the algorithm argument, which M1 addresses.
- The p = 0.010 single-bit anomaly in the sampler's own output is handled carefully but remains unexplained; an independent repeat is cheap.

## Novelty assessment

- *Rank-width sampling.* The algorithm is the linear-order (path) special case of Van den Nest et al. (2007) with the Schmidt structure of graph states supplying marginals; chain-rule sampling from exact marginals is the standard "perfect sampling" pattern (Ferris–Vidal). The note says this itself (§8). The engineering contributions (2^r rather than 2^{2r} storage, in-place single table, pdep-gather kernel) are real but incremental, and Kuyanov–Kissinger's ZX rank-width contraction and #260 do the same reduction. Not overclaimed in §8; the abstract's "The method is two standard ingredients joined by one identity" is fair.
- *Dual-table translation (§3.4).* Not new as mathematics: it is the shift theorem for the Walsh–Hadamard transform (a Z-phase on the suffix relabels the Schmidt basis Z^s|G_F⟩ → Z^{s+s'}|G_F⟩, i.e. translates the coefficient table). Its use to make slicing commute with suffix sampling is a sensible engineering step; present it as such.
- *Exact scoring of IBM's shots, 76/80-qubit runs, ceiling observation, d_min = 15 and the ≤14-qubit uniformity, hardware parity test:* genuine empirical contributions; the exact full-depth XEB on the published shots is the most useful item for the experimenters, and #260 computed the same probabilities independently, which strengthens rather than weakens it.
- *Underclaimed:* the doped instance (§7). Reducing a 468-T-gate circuit to a width-34 graph state via projected ancillas, computing an exact amplitude on one node in 26 minutes, and matching an independent tensor-network result string by string is the note's most interesting result and is listed under "does not show" (M2). The #242 sibling-string Porter–Thomas analysis and the local 8-qubit fidelity estimator are also new and buried.
- *Overclaimed:* "exact" at 80 qubits (M1); "no low-order trace at all"/"global depolarisation" (M4); "within 3–4 of the ceiling" and "as hard as a random graph" as statements about the construction rather than about the heuristic (M5); "hardware, not algorithm" (minor 3).

## Clarity and structure

The note is long and the strongest results are not where a reader looks first: the abstract leads with timings, while the exact per-shot benchmark, the stabilizer-only fine-structure verification, d_min = 15, and the doped reduction are deeper in. §4 item 7 (RNG and drawing-step forensics) is 27 lines of material that belongs in an appendix. §7's five-way table with three stale sentences around it reads as a work log. Phrases such as "A ceiling that holds for every graph", "the only lever left to the experimenters", and "leaves no trace" are stronger than the evidence. A busy expert will (a) read "exact" at 80 qubits and later discover a cap, (b) read "the doped instance is outside its reach" and then find §7, and (c) read "no low-order bias at all" and ask for the power.

## What a referee would require before acceptance, in priority order

1. Resolve M1: either a rigorous bound on the cap-exceedance mass (exact Var(z) via enumeration of the 2⁴⁰ in-prefix stabilizers, then Chebyshev), or rename the 80-qubit result as capped rejection sampling with a measured bound, and fix line 118.
2. Remove the contradictions on the doped instance (lines 397–400, 573) and move §7's result into the abstract.
3. Commit outputs for every number currently backed only by a test script (M3), or delete line 7's guarantee and the "1.5× per sliced vertex" figure.
4. Restate the low-order-marginal null results with their detection power (M4) and drop "global depolarisation".
5. Correct the frontier statements to the gaps actually found (3–5), fix the n = 80 row's update count, and state node memory explicitly in all threshold claims (M5, M6).
6. Put the fidelity-model caveat in the abstract and report a test of the Clifford instance's z-distribution or remove the "unlike the Clifford instance" claim (M7).
7. Fix the stale §7 sentences (minor 1–2), cite Hein–Eisert–Briegel for the Schmidt structure, name the Walsh shift theorem, and complete the pending §7 placeholder.
8. Optional but cheap: a second 1,012-sample run with a new seed base; 1,000 uniform control strings for the hardware XEB.

## Recommendation

**Major revision.** The computational work is careful and the artifact trail checks out line by line; I found no number that contradicts its artifact except the n = 80 update count. But the headline "exact" at 80 qubits is model-dependent and is in direct tension with the note's own distinction from #260; the text contradicts itself on the doped instance; several quoted numbers have no committed artifact despite an explicit promise that all do; and the physical interpretations placed on the null results and on the heuristic frontier go beyond what the data support. All of this is fixable without new large computations, with the exception of item 1, which needs either a new (modest) enumeration or a change of wording.
