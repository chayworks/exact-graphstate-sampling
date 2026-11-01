"""String-by-string comparison of our exact log2-probabilities of IBM's doped shots with the exact tensor-network
amplitudes of tracker entry #242 (Manabe, Gu, Pan; Zenodo 10.5281/zenodo.21912448, CC BY 4.0). Their strings are
q0-first; ours are the Qiskit little-endian strings of the IBM file (character k = qubit 69-k), so a string matches when
reversed. Usage:
  python doped_tn242_compare.py <amplitude_batches.npz | tn242_log2p.txt> <scored.txt> [<scored2.txt> ...]
Scored files are doped_sampler.py outputs: '<bits> <log2Pr> <seconds>' (score mode) or 'ACC|rej <bits> <log2Pr> ...'
(sample mode; sampled strings are not IBM shots and are skipped unless they happen to coincide). With an .npz the
extracted table is also written next to it as tn242_log2p.txt ('<q0-first string> <log2 p>')."""
import sys, math

def load_tn(path):
    if path.endswith(".npz"):
        import numpy as np
        z = np.load(path); rf = float(z["recovery_factor"]); bs = z["bitstrings_q0_first"]
        amp = z["raw_vectors"][np.arange(len(bs)), z["selected_vector_indices"]] * rf
        lp = np.log2(np.abs(amp) ** 2)
        out = path[: -len(".npz")].rsplit("/", 1)[0] + "/tn242_log2p.txt" if "/" in path else "tn242_log2p.txt"
        with open(out, "w") as f:
            for b, v in zip(bs, lp): f.write(f"{b} {v:.6f}\n")
        return {str(b): float(v) for b, v in zip(bs, lp)}
    return {l.split()[0]: float(l.split()[1]) for l in open(path) if l.strip()}

if __name__ == "__main__":
    tn = load_tn(sys.argv[1]); rows = []
    for fn in sys.argv[2:]:
        for l in open(fn):
            t = l.split()
            if not t or t[0].startswith("#"): continue
            if t[0] in ("ACC", "rej"): t = t[1:]
            if len(t) < 2 or len(t[0]) != 70: continue
            s, lp = t[0], float(t[1]); key = s[::-1]
            if key in tn: rows.append((fn, s, lp, tn[key]))
    print(f"{len(tn)} TN strings; {len(rows)} of our scored strings match (reversed)")
    print(f"{'file':28s} {'string (Qiskit order)':70s} {'ours log2 p':>14s} {'TN log2 p':>12s} {'diff':>10s} {'rel. prob':>10s}")
    for fn, s, lp, tl in rows:
        print(f"{fn[-28:]:28s} {s} {lp:14.7f} {tl:12.6f} {lp - tl:10.2e} {2 ** (lp - tl) - 1:10.2e}")
    if rows:
        d = [r[2] - r[3] for r in rows]; m = sum(d) / len(d)
        print(f"log2 difference ours - TN: mean {m:.2e}, max |.| {max(abs(x) for x in d):.2e}, n = {len(d)}; "
              f"TN values are complex64 outputs of a 2^48 FMA contraction, so a relative spread of 1e-4 is within their precision")
