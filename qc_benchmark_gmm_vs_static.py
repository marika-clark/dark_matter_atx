#!/usr/bin/env python3
"""
qc_benchmark_gmm_vs_static.py
=============================

Benchmark a GMM QC result (per-cell flags CSV with a QC_pass column) against
standard STATIC literature thresholds (min_umi=500, max_mt=15%) across the 41
GSE206283 scRNA samples.

Central question: is the GMM aggressively discarding biologically intact cells
(e.g. low-count / quiescent immune subsets) that the static filter keeps? We
answer it by profiling the four agreement partitions and by measuring canonical
PBMC lineage-marker expression in the GMM-discarded-but-static-kept pool.

Use --flags/--tag to benchmark either the baseline or the tempered GMM:
    python qc_benchmark_gmm_vs_static.py                       # baseline
    python qc_benchmark_gmm_vs_static.py \
        --flags .../tempered/cell_qc_flags_tempered.csv.gz --tag tempered

Outputs (qc_reports/ and qc_reports/benchmark_figures/), suffixed by --tag:
  - qc_benchmark_summary[_tag].csv        method retention + per-partition stats
  - qc_benchmark_markers[_tag].csv        per-partition marker mean expr + %positive
  - qc_benchmark_per_sample[_tag].csv     per-sample static vs gmm retention
  - benchmark_figures/confusion_matrix[_tag].png
  - benchmark_figures/scatter_counts_vs_mt_by_partition[_tag].png
  - benchmark_figures/marker_violins_by_partition[_tag].png
  - benchmark_figures/marker_dotplot_by_partition[_tag].png
"""
from __future__ import annotations
import argparse
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
FLAGS   = os.path.join(BASE, "qc_reports", "cell_qc_flags.csv.gz")
OUTDIR  = os.path.join(BASE, "qc_reports")
FIGDIR  = os.path.join(OUTDIR, "benchmark_figures")
os.makedirs(FIGDIR, exist_ok=True)

MIN_UMI, MAX_MT = 500, 15.0
MARKERS = ["CD3D", "CD4", "CD8A", "MS4A1", "CD14", "FCGR3A", "NCAM1"]
LINEAGE = {"CD3D": "T", "CD8A": "CD8 T", "MS4A1": "B", "CD14": "Monocyte",
           "FCGR3A": "NK/Mono", "NCAM1": "NK", "CD4": "CD4 T"}
PARTS = ["kept_by_both", "kept_static_only", "kept_gmm_only", "removed_by_both"]
PAL = {"kept_by_both": "#4C9F70", "kept_static_only": "#C1352B",
       "kept_gmm_only": "#E0A030", "removed_by_both": "#B0B0B0"}


def load_marker_expr():
    """Per sample: normalize (CP10k+log1p) on full genes, then keep only markers."""
    meta = pd.read_csv(os.path.join(SAMPLES, "metadata.csv"))
    subs = []
    for i, row in meta.iterrows():
        s = str(row["sample"])
        a = sc.read_10x_mtx(os.path.join(SAMPLES, s), var_names="gene_symbols", gex_only=True)
        a.var_names_make_unique()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        present = [g for g in MARKERS if g in a.var_names]
        sub = a[:, present].copy()
        sub.obs["sample_id"] = s
        sub.obs["barcode"] = sub.obs_names.values
        subs.append(sub)
        print(f"[{i+1:2d}/41] {s}: {a.n_obs:>6,} cells | markers {len(present)}/7", flush=True)
        del a
    adata = ad.concat(subs, join="outer", fill_value=0.0, index_unique="_")
    adata = adata[:, [m for m in MARKERS if m in adata.var_names]].copy()
    return adata


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flags", default=FLAGS, help="per-cell flags CSV with a QC_pass column")
    ap.add_argument("--tag", default="", help="suffix for outputs, e.g. 'tempered'")
    args = ap.parse_args()
    flags_path = args.flags
    tag = ("_" + args.tag) if args.tag else ""

    def op(fname):
        if tag:
            stem, _, ext = fname.partition(".")
            fname = f"{stem}{tag}.{ext}"
        return os.path.join(OUTDIR, fname)

    def fp(fname):
        if tag:
            stem, _, ext = fname.partition(".")
            fname = f"{stem}{tag}.{ext}"
        return os.path.join(FIGDIR, fname)

    # ---- marker expression AnnData ----
    adata = load_marker_expr()
    print(f"*** marker AnnData: {adata.n_obs:,} cells x {adata.n_vars} markers")

    # ---- join GMM flags on (sample_id, barcode) ----
    flags = pd.read_csv(flags_path)
    flags["key"] = flags["sample_id"].astype(str) + "|" + flags["barcode"].astype(str)
    adata.obs["key"] = adata.obs["sample_id"].astype(str) + "|" + adata.obs["barcode"].astype(str)
    assert adata.obs["key"].is_unique, "non-unique (sample_id,barcode) keys in adata"
    assert flags["key"].is_unique, "non-unique keys in flags file"
    fl = flags.set_index("key")
    missing = (~adata.obs["key"].isin(fl.index)).sum()
    assert missing == 0, f"{missing} cells have no matching GMM flag"
    j = fl.loc[adata.obs["key"].values]
    adata.obs["n_counts"]      = j["nCount_RNA"].values
    adata.obs["n_genes"]       = j["nFeature_RNA"].values
    adata.obs["pct_counts_mt"] = j["percent.mt"].values
    adata.obs["prob_damaged"]  = j["prob_damaged"].values
    adata.obs["severity"]      = j["severity"].values

    # ---- classification flags (booleans in adata.obs) ----
    adata.obs["pass_gmm"]    = j["QC_pass"].values.astype(bool)
    adata.obs["pass_static"] = ((adata.obs["n_counts"] >= MIN_UMI) &
                                (adata.obs["pct_counts_mt"] <= MAX_MT)).values
    ps, pg = adata.obs["pass_static"], adata.obs["pass_gmm"]
    adata.obs["kept_by_both"]     = ps & pg
    adata.obs["removed_by_both"]  = (~ps) & (~pg)
    adata.obs["kept_static_only"] = ps & (~pg)      # static KEEPS, GMM REMOVES  <- key pool
    adata.obs["kept_gmm_only"]    = (~ps) & pg
    adata.obs["partition"] = np.select(
        [adata.obs.kept_by_both, adata.obs.kept_static_only,
         adata.obs.kept_gmm_only, adata.obs.removed_by_both], PARTS, default="NA")

    obs = adata.obs
    N = adata.n_obs

    # ---- summary CSV: method retention + per-partition metric stats ----
    rows = []
    for name, mask in [("STATIC_kept", ps), ("GMM_kept", pg)]:
        d = obs[mask.values]
        r = {"group": name, "n_cells": int(mask.sum()),
             "pct_of_total": round(100 * mask.mean(), 2)}
        for c in ["n_counts", "n_genes", "pct_counts_mt"]:
            r[f"{c}_mean"] = round(d[c].mean(), 3)
            r[f"{c}_median"] = round(d[c].median(), 3)
            r[f"{c}_p95"] = round(np.percentile(d[c], 95), 3)
        rows.append(r)
    for p in PARTS:
        d = obs[obs.partition == p]
        r = {"group": p, "n_cells": len(d), "pct_of_total": round(100 * len(d) / N, 2)}
        for c in ["n_counts", "n_genes", "pct_counts_mt"]:
            r[f"{c}_mean"] = round(d[c].mean(), 3) if len(d) else np.nan
            r[f"{c}_median"] = round(d[c].median(), 3) if len(d) else np.nan
            r[f"{c}_p95"] = round(np.percentile(d[c], 95), 3) if len(d) else np.nan
        rows.append(r)
    summary = pd.DataFrame(rows)
    summary.to_csv(op("qc_benchmark_summary.csv"), index=False)

    # ---- per-sample retention ----
    ps_df = obs.groupby("sample_id").agg(
        raw=("pass_static", "size"),
        static_kept=("pass_static", "sum"),
        gmm_kept=("pass_gmm", "sum")).reset_index()
    ps_df["static_pct"] = (100 * ps_df.static_kept / ps_df.raw).round(2)
    ps_df["gmm_pct"] = (100 * ps_df.gmm_kept / ps_df.raw).round(2)
    ps_df["gmm_minus_static_pp"] = (ps_df.gmm_pct - ps_df.static_pct).round(2)
    ps_df.to_csv(op("qc_benchmark_per_sample.csv"), index=False)

    # ---- marker integrity ----
    mk = MARKERS.copy()
    X = adata[:, mk].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    expr = pd.DataFrame(X, columns=mk, index=obs.index)
    expr["partition"] = obs["partition"].values
    mean_by = expr.groupby("partition")[mk].mean()
    pct_by = expr.groupby("partition")[mk].apply(lambda df: (df > 0).mean() * 100)
    any_pos = (expr[mk] > 0).any(axis=1)
    frac_pos = any_pos.groupby(expr["partition"]).mean().mul(100).round(2)
    mk_out = mean_by.round(4).add_suffix("_meanExpr").join(pct_by.round(2).add_suffix("_pctPos"))
    mk_out["pct_any_marker_pos"] = frac_pos
    mk_out.to_csv(op("qc_benchmark_markers.csv"))
    ref = mean_by.loc["kept_by_both"].replace(0, np.nan)
    mis = float((mean_by.loc["kept_static_only"] / ref).mean()) if "kept_static_only" in mean_by.index else np.nan

    # ===================== FIGURES =====================
    # (a) confusion matrix
    ct = pd.crosstab(obs.pass_static, obs.pass_gmm)
    ct = ct.reindex(index=[True, False], columns=[True, False])
    fig, axx = plt.subplots(figsize=(5.5, 4.6))
    ann = np.empty(ct.shape, dtype=object)
    for r in range(2):
        for c in range(2):
            ann[r, c] = f"{ct.values[r,c]:,}\n({100*ct.values[r,c]/N:.1f}%)"
    sns.heatmap(ct, annot=ann, fmt="", cmap="Blues", cbar=False, ax=axx,
                xticklabels=["pass", "fail"], yticklabels=["pass", "fail"])
    axx.set_xlabel("GMM"); axx.set_ylabel("Static (min_umi=500, max_mt=15%)")
    axx.set_title(f"Static vs GMM cell agreement{(' ['+args.tag+']') if args.tag else ''}")
    fig.tight_layout(); fig.savefig(fp("confusion_matrix.png"), dpi=200); plt.close(fig)

    # (b) scatter n_counts vs pct_counts_mt by partition
    fig, axx = plt.subplots(figsize=(8, 6.5))
    for p in ["kept_by_both", "kept_gmm_only", "removed_by_both", "kept_static_only"]:
        d = obs[obs.partition == p]
        axx.scatter(d.n_counts, d.pct_counts_mt, s=2, alpha=0.25, edgecolors="none",
                    rasterized=True, color=PAL[p], label=f"{p} (n={len(d):,})")
    axx.set_xscale("log"); axx.axvline(MIN_UMI, ls="--", lw=1, color="k")
    axx.axhline(MAX_MT, ls="--", lw=1, color="k")
    axx.set_xlabel("n_counts (log)"); axx.set_ylabel("pct_counts_mt")
    axx.set_title(f"Static vs GMM partitions{(' ['+args.tag+']') if args.tag else ''} (dashed = static thresholds)")
    lg = axx.legend(markerscale=6, fontsize=8, framealpha=0.9)
    for lh in lg.legend_handles:
        lh.set_alpha(1)
    fig.tight_layout(); fig.savefig(fp("scatter_counts_vs_mt_by_partition.png"), dpi=180); plt.close(fig)

    # (c) marker violins across partitions
    long = expr.melt(id_vars="partition", value_vars=mk, var_name="marker", value_name="expr")
    g = sns.catplot(data=long, x="partition", y="expr", col="marker", col_wrap=4,
                    kind="violin", order=PARTS, palette=PAL, cut=0, density_norm="width",
                    height=3, aspect=1.1, inner=None, sharey=False)
    for ax in g.axes.flat:
        for coll in ax.collections:
            coll.set_rasterized(True)
        ax.tick_params(axis="x", labelrotation=90, labelsize=7)
        ax.set_xlabel("")
    g.figure.suptitle(f"Canonical PBMC marker expression by partition{(' ['+args.tag+']') if args.tag else ''}", y=1.02)
    g.savefig(fp("marker_violins_by_partition.png"), dpi=150, bbox_inches="tight")
    plt.close(g.figure)

    # (d) marker dotplot: size=%pos, color=mean expr
    fig, axx = plt.subplots(figsize=(8, 4.2))
    xs = np.arange(len(mk)); ys = np.arange(len(PARTS))
    for yi, p in enumerate(PARTS):
        for xi, gene in enumerate(mk):
            axx.scatter(xi, yi, s=pct_by.loc[p, gene] * 6,
                        c=[mean_by.loc[p, gene]], cmap="Reds",
                        vmin=0, vmax=float(mean_by.values.max()), edgecolors="grey", linewidths=.3)
    axx.set_xticks(xs); axx.set_xticklabels([f"{m}\n{LINEAGE[m]}" for m in mk], fontsize=8)
    axx.set_yticks(ys); axx.set_yticklabels(PARTS)
    axx.set_title(f"Marker mean expr (color) & %positive (size){(' ['+args.tag+']') if args.tag else ''}")
    sm = plt.cm.ScalarMappable(cmap="Reds", norm=plt.Normalize(0, float(mean_by.values.max())))
    fig.colorbar(sm, ax=axx, label="mean log-norm expr")
    fig.tight_layout(); fig.savefig(fp("marker_dotplot_by_partition.png"), dpi=200); plt.close(fig)

    # persist obs flags
    obs.reset_index().rename(columns={"index": "cell"}).to_csv(
        op("qc_benchmark_cell_partitions.csv.gz"), index=False, compression="gzip")

    # ===================== console report =====================
    stat_pct = 100 * ps.mean(); gmm_pct = 100 * pg.mean()
    kso = obs[obs.kept_static_only]
    print("\n" + "=" * 68)
    print(f"QC BENCHMARK{(' ['+args.tag+']') if args.tag else ''} — GMM vs STATIC (min_umi=500, max_mt=15%)")
    print("=" * 68)
    print(f"Total cells: {N:,} across {obs.sample_id.nunique()} samples")
    print(f"STATIC retention: {int(ps.sum()):,} ({stat_pct:.1f}%)  | attrition {100-stat_pct:.1f}%")
    print(f"GMM    retention: {int(pg.sum()):,} ({gmm_pct:.1f}%)  | attrition {100-gmm_pct:.1f}%")
    print("\nPartitions:")
    for p in PARTS:
        n = int((obs.partition == p).sum())
        print(f"  {p:18} {n:>8,} ({100*n/N:5.1f}%)")
    print(f"\nKEY POOL — kept_static_only (static keeps, GMM discards): {len(kso):,} cells")
    print(f"  median n_counts={kso.n_counts.median():.0f}  median n_genes={kso.n_genes.median():.0f}  median mt={kso['pct_counts_mt'].median():.1f}%")
    print(f"  % expressing >=1 lineage marker: {frac_pos.get('kept_static_only', float('nan')):.1f}%")
    print(f"  (kept_by_both reference: {frac_pos.get('kept_by_both', float('nan')):.1f}%)")
    print(f"  Marker Integrity Score (mean marker expr ratio vs kept_by_both): {mis:.2f}")
    print("\nPer-marker mean log-norm expr (kept_static_only vs kept_by_both):")
    for m in mk:
        print(f"  {m:7} static_only={mean_by.loc['kept_static_only', m]:.3f}  both={mean_by.loc['kept_by_both', m]:.3f}")
    print(f"\nOutputs -> qc_reports/  (tag='{args.tag}')")


if __name__ == "__main__":
    main()
