# Supplementary Fig. S3 | Gene-set null test for the year-5 PC1 x time interaction
# Data: random66_null_results.csv (1,000 random 66-gene sets, from null66_and_pdhc_test.R)
# Observed: beta = -0.1831, empirical p = 0.007 (6/1000 random sets with |beta| >= |beta_obs|)

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42   # TrueType embedding: keep whole words in PDF text
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["axes.linewidth"] = 0.8
matplotlib.rcParams["xtick.major.width"] = 0.8
matplotlib.rcParams["ytick.major.width"] = 0.8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

BASE = r"E:\PPMI帕金森数据库专用\结果\新结果\实验\movementdisorders\new\revision_v2\2026-08-31-第四版-评审修订"

df = pd.read_csv(BASE + r"\random66_null_results.csv")
null_beta = df["beta"].to_numpy()
beta_obs = -0.1831

n = len(null_beta)
n_exceed = int(np.sum(np.abs(null_beta) >= abs(beta_obs)))
emp_p = (n_exceed + 1) / (n + 1)  # add-one correction, same as null66_and_pdhc_test.R

q95 = np.quantile(np.abs(null_beta), 0.95)
print(f"n = {n}, |beta| >= |obs|: {n_exceed}, empirical p = {emp_p:.3f}")
print(f"null |beta| 95th percentile = {q95:.4f}, max = {np.abs(null_beta).max():.4f}")

fig, ax = plt.subplots(figsize=(3.54, 3.54), dpi=300)  # square, ~90 mm single column

# Histogram of signed null estimates
bins = np.linspace(-0.28, 0.28, 29)
counts, edges = np.histogram(null_beta, bins=bins)
centers = (edges[:-1] + edges[1:]) / 2
bar_colors = np.where(np.abs(centers) >= abs(beta_obs), "#C1272D", "#BDBDBD")
ax.bar(centers, counts, width=(edges[1] - edges[0]) * 0.92, color=bar_colors,
       edgecolor="white", linewidth=0.4)

# Observed estimate
ax.axvline(beta_obs, color="#C1272D", lw=1.4, ls="--", zorder=5)

ax.set_xlabel("Year-5 × gene-set interaction (β)", fontsize=8)
ax.set_ylabel("Number of random sets", fontsize=8)
ax.tick_params(labelsize=7.5)
ax.spines[["top", "right"]].set_visible(False)
ax.set_xlim(-0.29, 0.29)

fig.tight_layout(pad=0.4)
for ext in ("pdf", "png"):
    fig.savefig(f"{BASE}\\FigureS3_null_test.{ext}", dpi=300, bbox_inches="tight")
print("saved FigureS3_null_test.pdf / .png")

# Verify PDF text extraction (whole words, no letter splitting)
from pypdf import PdfReader
r = PdfReader(f"{BASE}\\FigureS3_null_test.pdf")
txt = " ".join((p.extract_text() or "") for p in r.pages)
print("PDF text sample:", txt[:150].replace("\n", " "))
