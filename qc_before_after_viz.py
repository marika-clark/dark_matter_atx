#!/usr/bin/env python3
"""
qc_before_after_viz.py
======================

Before-vs-after QC diagnostic visualizations for the GSE206283 scRNA pipeline,
using the finalized TEMPERED GMM filter (qc_reports/tempered/cell_qc_flags_tempered.csv.gz).

- Violins / scatters / summary table use ALL 472,830 cells (metrics only; no matrix load).
- PCA(50)+UMAP embeddings use a seeded ~60k stratified subsample (proportional per
  sample) because this machine has ~8.6 GB RAM. Subsampling is labeled in the figures.
- No batch integration: this is a QC diagnostic, so per-sample structure is expected.

Outputs -> qc_reports/visualizations/  (PNG + PDF)
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

BASE    = "/Users/marikaclark/dark_matter_atx/data/Giroux_2022/GSE206283_scRNA_RAW"
SAMPLES = os.path.join(BASE, "samples")
FLAGS   = os.path.join(BASE, "qc_reports", "tempered", "cell_qc_flags_tempered.csv.gz")
OUTDIR  = os.path.join(BASE, "qc_reports", "visualizations")
os.makedirs(OUTDIR, exist_ok=True)

TARGET_SUB = 60000
SEED = 0
sc.settings.verbosity = 1


def savefig(fig, stem, dpi=200):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUTDIR, f"{stem}.{ext}"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def pipeline(a):
    """Standard scanpy QC-diagnostic embedding on a raw-counts AnnData (in place-ish)."""
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=2000)
    a = a[:, a.var.highly_variable].copy()
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=50, svd_solver="arpack", random_state=SEED)
    sc.pp.neighbors(a, n_pcs=50, n_neighbors=15, random_state=SEED)
    sc.tl.umap(a, random_state=SEED)
    return a


def umap_scatter(a, color, stem, title, cmap="viridis", categorical=False):
    xy = a.obsm["X_umap"]
    fig, ax = plt.subplots(figsize=(7, 6))
    if categorical:
        cats = {"kept (pass)": ("#4C9F70", a.obs[color].values), }
        keep = a.obs[color].values.astype(bool)
        ax.scatter(xy[~keep, 0], xy[~keep, 1], s=3, alpha=0.4, rasterized=True,
                   color="#C1352B", edgecolors="none", label=f"removed (n={int((~keep).sum()):,})")
        ax.scatter(xy[keep, 0], xy[keep, 1], s=3, alpha=0.4, rasterized=True,
                   color="#4C9F70", edgecolors="none", label=f"kept (n={int(keep.sum()):,})")
        lg = ax.legend(markerscale=4, fontsize=9, loc="best")
        for lh in lg.legend_handles:
            lh.set_alpha(1)
    else:
        vals = a.obs[color].values
        sctr = ax.scatter(xy[:, 0], xy[:, 1], s=3, alpha=0.5, rasterized=True,
                          c=vals, cmap=cmap, vmin=0, vmax=np.percentile(vals, 99),
                          edgecolors="none")
        fig.colorbar(sctr, ax=ax, label=color)
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2"); ax.set_title(title)
    ax.set_xticks([]); ax.set_yticks([])
    savefig(fig, stem)


def main():
    # ============ PART A: full-cohort metrics (no matrix load) ============
    f = pd.read_csv(FLAGS)
    f = f.rename(columns={"percent.mt": "pct_mt"})
    before = f.copy(); before["stage"] = "before (raw)"
    after = f[f.QC_pass].copy(); after["stage"] = "after (tempered GMM)"
    both = pd.concat([before, after], ignore_index=True)
    order = ["before (raw)", "after (tempered GMM)"]

    # --- violins before vs after ---
    metrics = [("nFeature_RNA", "genes / cell", True),
               ("nCount_RNA", "UMIs / cell", True),
               ("pct_mt", "percent.mt", False)]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, (col, lab, logy) in zip(axes, metrics):
        sns.violinplot(data=both, x="stage", y=col, order=order, hue="stage",
                       palette={"before (raw)": "#9AA7B0", "after (tempered GMM)": "#4C9F70"},
                       legend=False, cut=0, density_norm="width", inner="quartile", ax=ax)
        for c in ax.collections:
            c.set_rasterized(True)
        if logy:
            ax.set_yscale("log")
        ax.set_title(lab); ax.set_xlabel("")
    fig.suptitle("QC metric distributions before vs after tempered GMM filtering (all 472,830 cells)", y=1.02)
    savefig(fig, "violin_before_after")

    # --- boundary scatters (colored by QC_pass) ---
    keep = f.QC_pass.values.astype(bool)
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    for ax, (xcol, ycol, ylog, ylab) in zip(
            axes, [("nCount_RNA", "pct_mt", False, "percent.mt"),
                   ("nCount_RNA", "nFeature_RNA", True, "nFeature_RNA (log)")]):
        ax.scatter(f[xcol][~keep], f[ycol][~keep], s=2, alpha=0.2, rasterized=True,
                   color="#C1352B", edgecolors="none", label=f"removed (n={int((~keep).sum()):,})")
        ax.scatter(f[xcol][keep], f[ycol][keep], s=2, alpha=0.2, rasterized=True,
                   color="#4C9F70", edgecolors="none", label=f"kept (n={int(keep.sum()):,})")
        ax.set_xscale("log")
        if ylog:
            ax.set_yscale("log")
        ax.set_xlabel("nCount_RNA (log)"); ax.set_ylabel(ylab)
        lg = ax.legend(markerscale=6, fontsize=9)
        for lh in lg.legend_handles:
            lh.set_alpha(1)
    fig.suptitle("Tempered GMM decision boundary (kept vs removed)", y=1.00)
    savefig(fig, "scatter_boundary_by_qc")

    # --- summary markdown table ---
    def med(d, c): return float(np.median(d[c]))
    tbl = [
        "# QC before vs after — tempered GMM (GSE206283, 41 samples)", "",
        "| Metric | Before (raw) | After (tempered GMM) |",
        "|---|---|---|",
        f"| Total cells | {len(before):,} | {len(after):,} ({100*len(after)/len(before):.1f}% retained) |",
        f"| Median genes / cell | {med(before,'nFeature_RNA'):.0f} | {med(after,'nFeature_RNA'):.0f} |",
        f"| Median UMIs / cell | {med(before,'nCount_RNA'):.0f} | {med(after,'nCount_RNA'):.0f} |",
        f"| Median percent.mt | {med(before,'pct_mt'):.2f}% | {med(after,'pct_mt'):.2f}% |",
        "",
        "_Embeddings below use a seeded ~60k stratified subsample (8.6 GB RAM); "
        "metrics above use all cells._", ""]
    with open(os.path.join(OUTDIR, "qc_before_after_summary.md"), "w") as fh:
        fh.write("\n".join(tbl))
    print("\n".join(tbl))

    # ============ PART B: ~60k stratified subsample embeddings ============
    meta = pd.read_csv(os.path.join(SAMPLES, "metadata.csv"))
    rng = np.random.default_rng(SEED)
    frac = TARGET_SUB / len(f)
    subs = []
    for i, row in meta.iterrows():
        s = str(row["sample"])
        a = sc.read_10x_mtx(os.path.join(SAMPLES, s), var_names="gene_symbols", gex_only=True)
        a.var_names_make_unique()
        n_i = max(1, int(round(a.n_obs * frac)))
        idx = rng.choice(a.n_obs, size=min(n_i, a.n_obs), replace=False)
        sub = a[idx].copy()
        sub.obs["sample_id"] = s
        sub.obs["barcode"] = sub.obs_names.values
        subs.append(sub)
        print(f"[{i+1:2d}/41] {s}: sampled {sub.n_obs}/{a.n_obs}", flush=True)
        del a
    counts = ad.concat(subs, join="outer", fill_value=0.0, index_unique="_")
    counts.obs["key"] = counts.obs.sample_id.astype(str) + "|" + counts.obs.barcode.astype(str)
    fl = f.assign(key=f.sample_id.astype(str) + "|" + f.barcode.astype(str)).set_index("key")
    assert counts.obs["key"].isin(fl.index).all(), "subsample cell missing from flags"
    j = fl.loc[counts.obs["key"].values]
    counts.obs["pct_mt"] = j["pct_mt"].values
    counts.obs["QC_pass"] = j["QC_pass"].values.astype(bool)
    print(f"*** subsample: {counts.n_obs:,} cells x {counts.n_vars} genes "
          f"({int(counts.obs.QC_pass.sum()):,} pass)")

    # raw embedding
    raw = pipeline(counts.copy())
    umap_scatter(raw, "pct_mt", "umap_raw_percent_mt",
                 f"RAW UMAP (~{raw.n_obs//1000}k subsample) colored by percent.mt", cmap="inferno")
    umap_scatter(raw, "QC_pass", "umap_raw_qc_pass",
                 f"RAW UMAP colored by tempered GMM QC_pass", categorical=True)
    del raw

    # filtered embedding (same subsample, kept cells only)
    filt = pipeline(counts[counts.obs.QC_pass].copy())
    umap_scatter(filt, "pct_mt", "umap_filtered_percent_mt",
                 f"FILTERED UMAP (~{filt.n_obs//1000}k kept) colored by percent.mt", cmap="inferno")

    print("\n*** Visualizations complete ->", OUTDIR)
    for fn in sorted(os.listdir(OUTDIR)):
        print("   ", fn)


if __name__ == "__main__":
    main()
