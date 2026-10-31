"""Figures for the note, drawn from committed artifacts only.
  fig1_frontier.png   width vs n for the construction (frontier.jsonl), the n/2 ceiling, measured points, memory axis
  fig2_shots.png      distribution of z = 2^70 p over IBM's 1,726 shots vs the fitted clean/uniform mixture and the ideal
  fig3_information.png  spread of the exact per-step conditionals along the order: hardware shots and ideal samples
Run from the repository root: python paper/make_figures.py"""
import json, math, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
R = "results/"; OUT = "paper/fig/"
import os; os.makedirs(OUT, exist_ok=True)

# ---- fig 1: frontier -------------------------------------------------------------------------------------------
rec = [json.loads(l) for l in open(R + "frontier.jsonl") if l.strip()]
n = np.array([r["n"] for r in rec]); w = np.array([r["width"] for r in rec])
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.plot(n, n // 2, "k--", lw=1, label="ceiling for any graph: ⌊n/2⌋")
ax.plot(n, w, "o-", color="#1f77b4", label="best order found, generated instances (seed 1000+n)")
ax.plot([70], [30], "s", color="#d62728", ms=8, label="IBM's instance (#151): width 30, sampled exactly")
ax.plot([76], [33], "D", color="#d62728", ms=7, label="76-qubit generated instance: width 33, sampled exactly")
ax.plot([80], [35], "^", color="#2ca02c", ms=8, label="80-qubit generated instance: width 35 (order found)")
ax.set_xlabel("qubits n (CZ depth n)"); ax.set_ylabel("linear cut-rank width")
ax2 = ax.secondary_yaxis("right", functions=(lambda y: y + math.log2(8), lambda y: y - math.log2(8)))
ax2.set_ylabel("log₂ bytes of the in-place table")
for y, lab in [(math.log2(8.6e9), "laptop 8.6 GB"), (math.log2(128e9), "128 GB node"), (math.log2(512e9), "512 GB node"), (math.log2(2.8e14), "281 TB")]:
    ax.axhline(y - 3, color="grey", lw=0.5, ls=":"); ax.text(30.5, y - 3 + 0.3, lab, fontsize=7, color="grey")
ax.legend(fontsize=7, loc="upper left"); ax.set_xlim(28, 122); ax.set_ylim(8, 62); ax.grid(alpha=0.3)
ax.set_title("Hardness frontier of random-graph-state sampling", fontsize=10)
fig.tight_layout(); fig.savefig(OUT + "fig1_frontier.png", dpi=200); plt.close(fig)

# ---- fig 2: hardware shots -------------------------------------------------------------------------------------
lp = np.array([float(l) for l in open(R + "ibmxeb/ibm_out.txt") if l.strip()]); z = np.exp2(lp + 70)
zo = np.array([float(l.split()[1]) for l in open(R + "samples_nq70.txt") if len(l.split()) > 1 and len(l.split()[0]) == 70]); zo = np.exp2(zo + 70)
F = 0.4815
fig, ax = plt.subplots(figsize=(6.4, 4.2))
bins = np.linspace(0, 6, 49); x = np.linspace(0, 6, 400)
ax.hist(z, bins=bins, density=True, alpha=0.55, color="#1f77b4", label="IBM hardware shots, N = 1,726 (exact p)")
ax.hist(zo, bins=bins, density=True, histtype="step", color="#d62728", lw=1.2, label="our exact samples, N = 1,012 (ideal)")
ax.plot(x, (F * x + 1 - F) * np.exp(-x), "k-", lw=1.2, label=f"mixture fit: clean with prob F = {F:.2f}, else uniform")
ax.plot(x, x * np.exp(-x), "k:", lw=1, label="ideal (Porter–Thomas): z·e⁻ᶻ")
ax.plot(x, np.exp(-x), "k--", lw=1, label="uniform strings: e⁻ᶻ")
ax.set_xlabel("z = 2⁷⁰ · p(x)"); ax.set_ylabel("density"); ax.legend(fontsize=7); ax.grid(alpha=0.3)
ax.set_title("Exact probabilities of IBM's post-selected shots (XEB 0.466 ± 0.031)", fontsize=10)
fig.tight_layout(); fig.savefig(OUT + "fig2_shots.png", dpi=200); plt.close(fig)

# ---- fig 3: information along the order ------------------------------------------------------------------------
def load(p):
    A = np.array([[float(v) for v in l.split()] for l in open(p) if len(l.split()) == 71]); return A[:, 1:]
ch, ci = load(R + "chain/ibm_fwd_chain.txt"), load(R + "chain/ideal_fwd_chain.txt")
k = np.arange(1, 71)
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.semilogy(k, np.maximum((2 * ci - 1).std(0), 1e-9), "o-", ms=3, color="#d62728", label="ideal samples: std of 2c_k − 1")
ax.semilogy(k, np.maximum((2 * ch - 1).std(0), 1e-9), "s-", ms=3, color="#1f77b4", label="hardware shots: std of 2c_k − 1")
ax.semilogy(k, np.maximum((2 * ci - 1).mean(0), 1e-9), "^--", ms=3, color="#2ca02c", label="ideal samples: mean of 2c_k − 1 (the XEB signal)")
ax.axvline(22.5, color="grey", ls=":", lw=1); ax.text(23, 2e-9, "cut rank < prefix length from step 23", fontsize=7, color="grey")
ax.set_xlabel("step k along the width-30 order"); ax.set_ylabel("|2·Pr(b_k | b₁…b_{k−1}) − 1|")
ax.set_ylim(1e-9, 2); ax.legend(fontsize=7, loc="upper left"); ax.grid(alpha=0.3, which="both")
ax.set_title("Where the information in a string lives: the last ~8 conditionals", fontsize=10)
fig.tight_layout(); fig.savefig(OUT + "fig3_information.png", dpi=200); plt.close(fig)
print("wrote", OUT)
