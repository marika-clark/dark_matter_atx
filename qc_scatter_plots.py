#!/usr/bin/env python3
"""
qc_scatter_plots.py
===================

Standalone cell-level QC visualization for the GSE206283 (Giroux 2022) scRNA
cohort. CellExpress `--only-qc` renders density plots but does not persist an
h5ad and does not show *cell-level* joint relationships, so this script
recomputes per-cell QC metrics directly from the raw 10X trios with Scanpy and
produces the standard QC diagnostic plots.

Outputs (into --outdir, default: <dataset>/qc_figures/):
  1. qc_scatter_counts_vs_genes.png  -- n_counts vs n_genes, colored by pct_counts_mt
  2. qc_scatter_counts_vs_mt.png     -- n_counts vs pct_counts_mt
  3. qc_violins_by_sample.png        -- n_genes / n_counts / pct_counts_mt violins per sample_id
  qc_metrics.csv.gz                  -- cached per-cell metrics (so re-plotting is instant)

Design notes:
  * Only obs-level metrics are retained per sample; the full ~470k x 32k matrix
    is never concatenated, keeping memory low.
  * Scatter points are rasterized (dpi-bounded) so the PNGs stay small despite
    ~470k points; optional --hexbin swaps to hexbin density aggregation.
  * Proposed QC thresholds are drawn as reference lines so the cloud can be
    judged against them.
"""

from __future__ import annotations
import argparse
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

DEFAULT_INPUT = "/Users/marikaclark/dark_matter_atx/data/Giroux_2022/GSE206283_scRNA_RAW/samples"

# Proposed (convention-based) thresholds, drawn as reference lines only.
THRESH = {"min_umi": 500, "min_genes": 250, "max_mt": 15.0}

# Consistent colors for the cohort's severity groups.
SEV_PALETTE = {
    "Mild": "#4C9F70", "Severe": "#C1352B",
    "Exposed": "#E0A030", "Healthy": "#3B6EA5", "NA": "#999999",
}


def compute_metrics(input_dir: str) -> pd.DataFrame:
    """Load each sample's 10X trio, compute per-cell QC, return combined obs."""
    meta = pd.read_csv(os.path.join(input_dir, "metadata.csv"))
    if "sample" not in meta.columns:
        sys.exit("metadata.csv must contain a 'sample' column")
    frames = []
    for i, row in meta.reset_index(drop=True).iterrows():
        s = str(row["sample"])
        d = os.path.join(input_dir, s)
        a = sc.read_10x_mtx(d, var_names="gene_symbols", gex_only=True)
        a.var_names_make_unique()
        a.var["mt"] = a.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            a, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
        )
        obs = a.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
        obs.columns = ["n_counts", "n_genes", "pct_counts_mt"]
        obs["sample_id"] = s
        obs["severity"] = str(row.get("severity", "NA"))
        frames.append(obs.reset_index(drop=True))
        print(f"[{i + 1:2d}/{len(meta)}] {s}: {a.n_obs:>6,} cells, "
              f"{int(a.var['mt'].sum())} MT genes", flush=True)
    df = pd.concat(frames, ignore_index=True)
    print(f"*** combined: {len(df):,} cells across {df['sample_id'].nunique()} samples")
    return df, meta


def scatter_counts_vs_genes(df, outpath, hexbin, dpi):
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    if hexbin:
        hb = ax.hexbin(df.n_counts, df.n_genes, C=df.pct_counts_mt,
                       reduce_C_function=np.mean, gridsize=120, xscale="log",
                       yscale="log", cmap="viridis", mincnt=1, vmin=0, vmax=20)
        cbar = fig.colorbar(hb, ax=ax)
    else:
        scat = ax.scatter(df.n_counts, df.n_genes, c=df.pct_counts_mt, s=2,
                          alpha=0.3, cmap="viridis", vmin=0, vmax=20,
                          edgecolors="none", rasterized=True)
        ax.set_xscale("log"); ax.set_yscale("log")
        cbar = fig.colorbar(scat, ax=ax)
    cbar.set_label("pct_counts_mt (clipped at 20%)")
    ax.axvline(THRESH["min_umi"], ls="--", lw=1, color="crimson")
    ax.axhline(THRESH["min_genes"], ls="--", lw=1, color="crimson")
    ax.set_xlabel("n_counts (UMIs per cell, log)")
    ax.set_ylabel("n_genes per cell (log)")
    ax.set_title("UMI vs gene count, colored by mitochondrial %\n"
                 f"ref lines: min_umi={THRESH['min_umi']}, min_genes={THRESH['min_genes']}")
    fig.tight_layout(); fig.savefig(outpath, dpi=dpi); plt.close(fig)
    print("wrote", outpath)


def scatter_counts_vs_mt(df, outpath, hexbin, dpi):
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    if hexbin:
        hb = ax.hexbin(df.n_counts, df.pct_counts_mt, gridsize=120, xscale="log",
                       cmap="mako_r", mincnt=1, bins="log")
        cbar = fig.colorbar(hb, ax=ax); cbar.set_label("log10(cell count)")
    else:
        ax.scatter(df.n_counts, df.pct_counts_mt, s=2, alpha=0.3,
                   color="#3B6EA5", edgecolors="none", rasterized=True)
        ax.set_xscale("log")
    ax.axvline(THRESH["min_umi"], ls="--", lw=1, color="crimson")
    ax.axhline(THRESH["max_mt"], ls="--", lw=1, color="crimson")
    ax.set_xlabel("n_counts (UMIs per cell, log)")
    ax.set_ylabel("pct_counts_mt")
    ax.set_title("UMI count vs mitochondrial %\n"
                 f"ref lines: min_umi={THRESH['min_umi']}, max_mt={THRESH['max_mt']}%")
    fig.tight_layout(); fig.savefig(outpath, dpi=dpi); plt.close(fig)
    print("wrote", outpath)


def _violin(ax, df, col, order, colors, logy):
    # hue=x with legend off is the cross-version-safe way to color per-group.
    try:
        sns.violinplot(data=df, x="sample_id", y=col, order=order, hue="sample_id",
                       palette=colors, legend=False, cut=0, inner="box",
                       linewidth=0.4, density_norm="width", ax=ax)
    except TypeError:
        # older seaborn: `scale` instead of `density_norm`, no hue/legend kwargs
        sns.violinplot(data=df, x="sample_id", y=col, order=order, palette=colors,
                       cut=0, inner="box", linewidth=0.4, scale="width", ax=ax)
    for coll in ax.collections:
        coll.set_rasterized(True)
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel("")
    ax.set_ylabel(col)
    ax.tick_params(axis="x", labelrotation=90, labelsize=6)


def violins(df, meta, outpath, dpi):
    order = [str(s) for s in meta["sample"].tolist()]
    sev = {str(r["sample"]): str(r.get("severity", "NA")) for _, r in meta.iterrows()}
    colors = [SEV_PALETTE.get(sev.get(s, "NA"), "#999999") for s in order]
    fig, axes = plt.subplots(3, 1, figsize=(20, 13), sharex=True)
    _violin(axes[0], df, "n_genes", order, colors, logy=True)
    _violin(axes[1], df, "n_counts", order, colors, logy=True)
    _violin(axes[2], df, "pct_counts_mt", order, colors, logy=False)
    axes[0].axhline(THRESH["min_genes"], ls="--", lw=1, color="crimson")
    axes[1].axhline(THRESH["min_umi"], ls="--", lw=1, color="crimson")
    axes[2].axhline(THRESH["max_mt"], ls="--", lw=1, color="crimson")
    handles = [Line2D([0], [0], marker="s", ls="", markerfacecolor=c,
                      markeredgecolor="none", label=k)
               for k, c in SEV_PALETTE.items() if k != "NA"]
    axes[0].legend(handles=handles, title="severity", ncol=4, loc="upper right",
                   fontsize=8, framealpha=0.9)
    axes[0].set_title("Per-sample QC distributions (violins colored by severity; "
                      "dashed = proposed thresholds)")
    fig.tight_layout(); fig.savefig(outpath, dpi=dpi); plt.close(fig)
    print("wrote", outpath)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", default=DEFAULT_INPUT,
                   help="dir with per-sample 10X subfolders + metadata.csv")
    p.add_argument("--outdir", default=None,
                   help="output dir (default: <parent-of-input>/qc_figures)")
    p.add_argument("--hexbin", action="store_true",
                   help="use hexbin density aggregation instead of rasterized points")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--cache", default=None,
                   help="path to a cached qc_metrics.csv.gz to skip recomputation")
    args = p.parse_args()

    outdir = args.outdir or os.path.join(
        os.path.dirname(args.input.rstrip("/")), "qc_figures")
    os.makedirs(outdir, exist_ok=True)

    if args.cache and os.path.exists(args.cache):
        print("loading cached metrics:", args.cache)
        df = pd.read_csv(args.cache)
        meta = pd.read_csv(os.path.join(args.input, "metadata.csv"))
    else:
        df, meta = compute_metrics(args.input)
        cache_path = os.path.join(outdir, "qc_metrics.csv.gz")
        df.to_csv(cache_path, index=False, compression="gzip")
        print("wrote", cache_path)

    scatter_counts_vs_genes(df, os.path.join(outdir, "qc_scatter_counts_vs_genes.png"),
                            args.hexbin, args.dpi)
    scatter_counts_vs_mt(df, os.path.join(outdir, "qc_scatter_counts_vs_mt.png"),
                         args.hexbin, args.dpi)
    violins(df, meta, os.path.join(outdir, "qc_violins_by_sample.png"), args.dpi)
    print("*** done ->", outdir)


if __name__ == "__main__":
    main()
