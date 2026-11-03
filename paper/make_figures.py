"""Figures for the note, drawn from committed artifacts only.
  fig1_frontier.png     width vs n for the construction (frontier.jsonl), the n/2 ceiling, measured points, memory axis
  fig2_shots.png        distribution of z = 2^70 p over IBM's 1,726 shots vs the fitted clean/uniform mixture and the ideal
  fig3_information.png  spread of the exact per-step conditionals along the order: hardware shots and ideal samples
Colours: a three-slot categorical palette validated for colour-vision deficiency (blue, orange, aqua).
Run from the repository root: python paper/make_figures.py"""
import json, math, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
R = "graphsamp/results/"; OUT = "paper/fig/"
import os; os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, AQUA, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 8.5, "axes.titlesize": 10, "axes.labelsize": 9, "axes.edgecolor": MUTED,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})

# ---- fig 1: frontier -------------------------------------------------------------------------------------------
rec = [json.loads(l) for l in open(R + "frontier.jsonl") if l.strip()]
n = np.array([r["n"] for r in rec]); w = np.array([r["width"] for r in rec])
fig, ax = plt.subplots(figsize=(6.6, 4.4))
# machine bands, in units of the in-place table (2^width complex64 entries = 2^(width+3) bytes)
bands = [(0, math.log2(8.6e9) - 3, "laptop", "#f3f3f1"), (math.log2(8.6e9) - 3, math.log2(512e9) - 3, "one cloud node (up to 512 GB)", "#e8e8e5"),
         (math.log2(512e9) - 3, math.log2(2.8e14) - 3, "supercomputer-scale memory", "#dcdcd8"), (math.log2(2.8e14) - 3, 70, "beyond any existing machine", "#cfcfca")]
for lo, hi, lab, col in bands:
    ax.axhspan(lo, hi, color=col, lw=0, zorder=0); ax.text(121, (max(lo, 8) + min(hi, 62)) / 2, lab, fontsize=7.5, color=MUTED, ha="right", va="center")
ax.plot(n, n // 2, ls="--", lw=1.2, color=INK, label="ceiling for every graph: ⌊n/2⌋")
ax.plot(n, w, "o-", ms=4, lw=1.6, color=BLUE, label="IBM's construction at depth n: best order found")
ax.plot([70], [30], "s", color=ORANGE, ms=8, mec="white", mew=1, zorder=5, label="sampled exactly: IBM's 70-qubit instance (#151), width 30")
ax.plot([76], [33], "D", color=ORANGE, ms=7, mec="white", mew=1, zorder=5, label="sampled exactly: 76-qubit instance, width 33")
ax.plot([80], [35], "^", color=ORANGE, ms=8, mec="white", mew=1, zorder=5, label="sampled (hybrid, capped): 80-qubit instance, width 35 sliced to 34")
ax.annotate("18 s / sample, laptop", (70, 30), (44, 33.5), fontsize=7.5, color=INK, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7))
ax.annotate("190 s / sample, 16 cores", (76, 33), (82, 27), fontsize=7.5, color=INK, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7))
ax.annotate("246 s / proposal, 32 cores, 128 GiB", (80, 35), (42, 39.5), fontsize=7.5, color=INK, arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7))
ax.set_xlabel("qubits n (CZ depth n)"); ax.set_ylabel("cut-rank width of the best linear order")
ax2 = ax.secondary_yaxis("right", functions=(lambda y: y + 3, lambda y: y - 3))
ax2.set_ylabel("log₂ bytes of the sampler's table", color=MUTED); ax2.spines["right"].set_visible(True)
ax2.set_yticks([33, 37, 40, 43, 47, 50, 53]); ax2.set_yticklabels(["8 GB", "128 GB", "1 TB", "8 TB", "128 TB", "1 PB", "8 PB"], fontsize=7.5, color=MUTED)
ax.legend(fontsize=7.2, loc="upper left"); ax.set_xlim(28, 122); ax.set_ylim(8, 62); ax.grid(alpha=0.25, lw=0.5)
ax.set_title("Where random-graph-state sampling becomes hard: memory doubles every two qubits")
fig.tight_layout(); fig.savefig(OUT + "fig1_frontier.png", dpi=220); plt.close(fig)

# ---- fig 2: hardware shots -------------------------------------------------------------------------------------
lp = np.array([float(l) for l in open(R + "ibmxeb/ibm_out.txt") if l.strip()]); z = np.exp2(lp + 70)
zo = np.array([float(l.split()[1]) for l in open(R + "samples_nq70.txt") if len(l.split()) > 1 and len(l.split()[0]) == 70]); zo = np.exp2(zo + 70)
F = 0.4815
fig, ax = plt.subplots(figsize=(6.6, 4.2))
bins = np.linspace(0, 6, 49); x = np.linspace(0, 6, 400)
ax.hist(z, bins=bins, density=True, alpha=0.6, color=BLUE, lw=0, label="IBM's hardware shots, N = 1,726 (exact Born probability of each)")
ax.hist(zo, bins=bins, density=True, histtype="step", color=ORANGE, lw=1.6, label="our exact samples of the ideal distribution, N = 1,012")
ax.plot(x, (F * x + 1 - F) * np.exp(-x), "-", color=INK, lw=1.4, label=f"fit to the shots: ideal with probability F = {F:.2f}, uniform otherwise")
ax.plot(x, x * np.exp(-x), ":", color=MUTED, lw=1.1, label="ideal device (Porter–Thomas): z·e⁻ᶻ")
ax.plot(x, np.exp(-x), "--", color=MUTED, lw=1.1, label="fully depolarised device: e⁻ᶻ")
ax.set_xlabel("z = 2⁷⁰ · p(x)   (1 = a uniformly random string)"); ax.set_ylabel("density over shots"); ax.legend(fontsize=7.2); ax.grid(alpha=0.25, lw=0.5)
ax.set_title("IBM's 70-qubit shots, scored exactly: linear XEB 0.466 ± 0.031, fidelity ≈ 0.46–0.48")
fig.tight_layout(); fig.savefig(OUT + "fig2_shots.png", dpi=220); plt.close(fig)

# ---- fig 3: information along the order ------------------------------------------------------------------------
def load(p):
    A = np.array([[float(v) for v in l.split()] for l in open(p) if len(l.split()) == 71]); return A[:, 1:]
ch, ci = load(R + "chain/ibm_fwd_chain.txt"), load(R + "chain/ideal_fwd_chain.txt")
k = np.arange(1, 71)
fig, ax = plt.subplots(figsize=(6.6, 4.2))
ax.axvspan(62.5, 70.5, color="#f3f3f1", lw=0, zorder=0); ax.text(66.5, 1e-5, "last 8 steps:\n100% of\nvar(log₂ Pr)", fontsize=7.5, color=MUTED, ha="center")
ax.semilogy(k, np.maximum((2 * ci - 1).std(0), 1e-9), "o-", ms=3, lw=1.2, color=ORANGE, label="ideal samples: spread (std) of 2c_k − 1")
ax.semilogy(k, np.maximum((2 * ch - 1).std(0), 1e-9), "s-", ms=3, lw=1.2, color=BLUE, label="hardware shots: spread (std) of 2c_k − 1")
mi = (2 * ci - 1).mean(0); se = (2 * ci - 1).std(0, ddof=1) / math.sqrt(ci.shape[0]); sig = (mi > 2 * se) & (mi > 1e-6)   # below 1e-6 the mean is at the print resolution of the conditionals, not signal
ax.semilogy(k[sig], mi[sig], "^", ms=4, color=AQUA, label="ideal samples: mean of 2c_k − 1 where above 2 SE and 10⁻⁶ (per-step XEB signal)")
ax.axvline(22.5, color=MUTED, ls=":", lw=1); ax.text(1, 1.5e-3, "first 22 conditionals are exactly ½\n(cut rank = prefix length)", fontsize=7.5, color=MUTED)
ax.set_xlabel("step k along the width-30 order (one qubit per step)"); ax.set_ylabel("|2·Pr(b_k | b₁…b_{k−1}) − 1|")
ax.set_ylim(1e-9, 2); ax.legend(fontsize=7.2, loc="upper left"); ax.grid(alpha=0.25, lw=0.5, which="both")
ax.set_title("Where the information in a 70-bit string lives: the last ~8 conditionals")
fig.tight_layout(); fig.savefig(OUT + "fig3_information.png", dpi=220); plt.close(fig)

# ---- fig 4: low-order parities of the hardware shots (simulation-free) ----------------------------------------
import json as _json
Z = np.load(R + "hw_marginals_z.npz"); M = _json.load(open(R + "hw_marginals.json"))
fig, axs = plt.subplots(1, 4, figsize=(9.6, 2.9), sharey=True)
xx = np.linspace(-5, 5, 300)
for i, k in enumerate((1, 2, 3, 4)):
    ax = axs[i]; zh, zi = Z[f"hw_k{k}"], Z[f"ideal_k{k}"]; bins = np.linspace(-5, 5, 26 if k == 1 else 51)
    ax.hist(zh, bins=bins, density=True, alpha=0.6, color=BLUE, lw=0, label="IBM hardware shots")
    ax.hist(zi, bins=bins, density=True, histtype="step", color=ORANGE, lw=1.3, label="our exact samples")
    ax.plot(xx, np.exp(-xx ** 2 / 2) / np.sqrt(2 * np.pi), "-", color=INK, lw=1.1, label="standard normal (exactly uniform marginal)")
    ph, pi = M["hardware"][str(k)]["chi2_p"], M["ideal"][str(k)]["chi2_p"]; nt = M["hardware"][str(k)]["tests"]
    ax.set_title(f"{k}-bit parities, {nt:,} tests", fontsize=9); ax.text(0.97, 0.95, f"χ² p: shots {ph:.2f}\nsamples {pi:.2f}", transform=ax.transAxes, fontsize=7.5, va="top", ha="right", color=MUTED)
    ax.set_xlabel("z-score of the parity mean"); ax.grid(alpha=0.25, lw=0.5); ax.set_xlim(-5, 5)
axs[0].set_ylabel("density"); h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=3, fontsize=7.5, bbox_to_anchor=(0.5, 0.0))
fig.suptitle("IBM's shots carry no low-order bias: every marginal on ≤ 14 qubits of the ideal output is exactly uniform", fontsize=9.5)
fig.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig(OUT + "fig4_marginals.png", dpi=220); plt.close(fig)
print("wrote", OUT)

# ---- fig 5: cut-rank profiles along the orders ---------------------------------------------------------------------
def gf2_rank(M):
    M = M.copy() % 2; r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]: M[i, :] ^= M[r, :]
        r += 1
        if r == rows: break
    return r
def profile_from_order(adj_path, order):
    A = np.array([[int(c) for c in l.strip()] for l in open(adj_path) if l.strip()], np.uint8); n = len(A)
    out = []
    for k in range(1, n + 1):
        P = order[:k]; F = order[k:]
        out.append(gf2_rank(A[np.ix_(F, P)]) if F else 0)
    return out
def profile_from_plan(path):
    L = open(path).read().split("\n"); n, steps, R = map(int, L[0].split())
    return [int(l.split()[1]) for l in L[1:1 + steps]]
w31 = profile_from_order(R + "nq70_graph_adjacency.txt", json.load(open("graphsamp/cloud/nq70_step_qubits.json")))
w30 = profile_from_plan("graphsamp/cloud2/nq70_inplace_plan.txt")
d34 = profile_from_plan("graphsamp/cloud_doped/doped_m20_inplace_plan.txt")
f80 = profile_from_plan("graphsamp/cloud_nq80/fi80h_s0.txt")
fig, ax = plt.subplots(figsize=(6.6, 4.0))
for prof, lab, col, ls in [(w31, "IBM's 70-qubit instance, production order (width 31)", BLUE, "-"), (w30, "same instance, width-30 order", AQUA, "-"),
                           (f80, "80-qubit instance, sliced order (width 34)", ORANGE, "-"), (d34, "doped instance: 538-vertex gadget graph, sweep order (width 34)", INK, "--")]:
    n = len(prof); ax.plot(np.arange(1, n + 1) / n, prof, ls, lw=1.4, color=col, label=lab)
x = np.linspace(0, 1, 200); ax.plot(x, 70 * np.minimum(x, 1 - x), ":", color=MUTED, lw=1, label="ceiling for 70 qubits: min(k, n − k)")
ax.set_xlabel("fraction of the order processed, k / n"); ax.set_ylabel("cut rank r_k after step k"); ax.set_ylim(0, 40)
ax.legend(fontsize=7.2, loc="upper right"); ax.grid(alpha=0.25, lw=0.5)
ax.set_title("Cut-rank profiles along the orders: the table at step k has 2^{r_k} entries")
fig.tight_layout(); fig.savefig(OUT + "fig5_rankprofile.png", dpi=220); plt.close(fig)
print("fig5 maxima:", max(w31), max(w30), max(f80), max(d34))

# ---- fig 6: hybrid-sampler prefix marginals against the cap --------------------------------------------------------
def zs(path):
    return np.array([float(l.split("z=")[1].split()[0]) for l in open(path) if l.startswith(("ACC", "rej"))])
runs = [(zs(R + "fstate_hybrid_samples.txt"), "70-qubit Clifford instance, 150 proposals (K = 50)", BLUE),
        (zs(R + "nq80/samples_nq80.txt"), "80-qubit instance, 80 proposals (K = 60)", ORANGE),
        (zs(R + "doped/run11/samples.txt"), "doped instance, 13 proposals (K = 50, run 11)", AQUA)]
fig, ax = plt.subplots(figsize=(6.6, 3.6))
bins = np.linspace(0.994, 1.006, 25)
for z, lab, col in runs:
    ax.hist(z, bins=bins, histtype="step", lw=1.6, color=col, label=f"{lab}: mean {z.mean():.4f}, sd {z.std(ddof=1):.4f}")
ax.text(1.0058, ax.get_ylim()[1] * 0.55, "cap c = 1.03:\n30 sd to the right →", ha="right", fontsize=8, color=INK)
ax.set_xlabel("prefix marginal z = 2^K Pr(first K bits) in units of the uniform distribution"); ax.set_ylabel("proposals")
ax.set_xlim(0.994, 1.006); ax.legend(fontsize=7.2, loc="upper left"); ax.grid(alpha=0.25, lw=0.5)
ax.set_title("Hybrid sampler: prefix marginals of every proposal, three instances")
fig.tight_layout(); fig.savefig(OUT + "fig6_cap.png", dpi=220); plt.close(fig)
print("fig6 z ranges:", [(float(z.min()), float(z.max())) for z, _, _ in runs])
