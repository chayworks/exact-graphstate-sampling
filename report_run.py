"""Report for an exact-sampling run of an n-qubit instance.
Usage: python report_run.py <samples.txt> <n> [<dp_out.txt> <xcheck_make.txt>]
samples.txt: "<bits in step order> <log2Pr> <seed>" per line (gsamp/gsamp2 sample output).
dp_out.txt: per checked string "re im" of the amplitude from dp.cpp (|amp|^2 * 2^-n is the probability);
xcheck_make.txt: the JSON line written by xcheck_make.py (which strings were checked, with the sampler's log2 Pr)."""
import sys, json, math, numpy as np
SF, n = sys.argv[1], int(sys.argv[2])
rows = [l.split() for l in open(SF) if len(l.split()) >= 2 and len(l.split()[0]) == n]
lp = np.array([float(r[1]) for r in rows]); N = len(lp); z = np.exp2(lp + n)
print(f"{N} exact samples of the {n}-qubit instance; all distinct: {len(set(r[0] for r in rows)) == N}")
print(f"mean log2 Pr {lp.mean():.3f} (uniform would be {-n}); linear XEB against own exact probabilities "
      f"{z.mean() - 1:+.3f} +/- {z.std(ddof=1) / math.sqrt(N):.3f} (Porter-Thomas: +1; estimates 2^n sum p^2 - 1)")
if len(sys.argv) > 4:
    amps = [[float(x) for x in l.split()] for l in open(sys.argv[3]) if l.strip()]
    meta = json.loads([l for l in open(sys.argv[4]) if l.startswith("{")][-1])
    print(f"independent engine (dp.cpp, {len(meta['sliced'])} sliced vertices, width {meta['dp_width']}):")
    for (s, lps), a in zip(meta["strings"], amps):
        lpd = math.log2(a[0] ** 2 + a[1] ** 2) - n
        print(f"  {s[:16]}...  sampler {lps:.7f}  dp {lpd:.7f}  diff {lpd - lps:+.1e}")
