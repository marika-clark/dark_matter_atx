#!/usr/bin/env python3
"""
run_cellexpress.py
===================

A thin, auditable wrapper around AstraZeneca's CellExpress pipeline
(https://github.com/AstraZeneca/cellatria/tree/main/cellexpress).

WHY THIS EXISTS
---------------
CellExpress is a vetted, published, standalone scRNA-seq processing
pipeline (QC -> normalization -> clustering -> annotation -> report).
The point of orchestrating it through Claude Code is to let Claude
reason about *which parameters to use* (a scientific judgment call
that deserves review), while NOT letting Claude re-implement or
freelance the actual analysis logic (that's what CellExpress is for).

This script draws a hard line between those two responsibilities:
  - Claude/you decide the parameters -> this script just assembles
    them into the JSON config CellExpress expects.
  - CellExpress does the actual science -> this script never touches
    the count matrices itself, it only invokes the pipeline and reads
    back its own reported outputs.

That means Claude Code can be told "use run_cellexpress.py to run the
pipeline" instead of being given raw bash access to call main.py with
an arbitrary, unreviewed command line.

WORKFLOW THIS SUPPORTS
-----------------------
1. `--only-qc` pass: get QC metrics without committing to thresholds.
2. `--dry-run`: build and print the config WITHOUT executing anything,
   so a human (or the CLAUDE.md gatekeeping rule) can review the
   parameter choices before a real run happens.
3. Real run: executes CellExpress, streams its output live, then
   locates the output folder it produced and prints a structured
   summary (cell counts, cluster counts, annotation coverage) so the
   "verify outputs before calling it done" step has something
   concrete to check against, rather than requiring someone to
   manually dig through the HTML report.

USAGE
-----
    # Step 1: look at QC metrics only, no processing
    python run_cellexpress.py \\
        --input /data/my_project --project my_project \\
        --species hs --tissue lung --disease asthma \\
        --only-qc

    # Step 2: review proposed parameters without running anything
    python run_cellexpress.py \\
        --input /data/my_project --project my_project \\
        --species hs --tissue lung --disease asthma \\
        --min-umi-per-cell 1000 --max-mt-percent 12.5 \\
        --resolution 0.8 --dry-run

    # Step 3: the real run, once the config above is approved
    python run_cellexpress.py \\
        --input /data/my_project --project my_project \\
        --species hs --tissue lung --disease asthma \\
        --min-umi-per-cell 1000 --max-mt-percent 12.5 \\
        --resolution 0.8

    # Or reuse a config you already reviewed/saved:
    python run_cellexpress.py --config-in /path/to/reviewed_config.json
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Where to find CellExpress itself.
#
# Why configurable rather than hardcoded: CellExpress might be a local git
# clone (e.g. ~/cellatria/cellexpress/main.py) or the copy baked into the
# CellAtria Docker image (/opt/cellatria/cellexpress/main.py). Rather than
# guess wrong silently, we require it to be findable and fail loudly if not.
# ---------------------------------------------------------------------------
DEFAULT_CELLEXPRESS_CANDIDATES = [
    Path("/opt/cellatria/cellexpress/main.py"),          # inside the Docker image
    Path.home() / "cellatria" / "cellexpress" / "main.py",  # a typical local clone
]


def find_cellexpress_main(explicit_path: str | None) -> Path:
    """Resolve the path to CellExpress's main.py.

    We fail loudly here rather than silently falling back to a guess,
    because running the wrong copy of the pipeline (e.g. an out-of-date
    local clone vs. the Docker image's copy) would silently undermine the
    reproducibility guarantees CellExpress is supposed to provide.
    """
    if explicit_path:
        p = Path(explicit_path).expanduser().resolve()
        if not p.is_file():
            raise FileNotFoundError(f"--cellexpress-path does not exist: {p}")
        return p

    for candidate in DEFAULT_CELLEXPRESS_CANDIDATES:
        if candidate.is_file():
            return candidate

    raise FileNotFoundError(
        "Could not locate cellexpress/main.py. Pass it explicitly with "
        "--cellexpress-path /path/to/cellexpress/main.py."
    )


# ---------------------------------------------------------------------------
# Config assembly
#
# WHY we only include parameters the caller actually specified (not every
# CellExpress default): the point of the review step is to show exactly
# which decisions were made deliberately vs. left to CellExpress's own
# defaults. Padding the config with every default value would make it
# harder to tell "chosen" from "inherited" at a glance during review.
# ---------------------------------------------------------------------------

# Maps our argparse dest names -> the exact key CellExpress expects in its
# JSON config (per cellexpress/README.md "Pipeline Arguments" section).
CONFIG_KEY_MAP: dict[str, str] = {
    # General (required)
    "input": "input",
    "project": "project",
    "species": "species",
    "tissue": "tissue",
    "disease": "disease",
    # QC
    "min_umi_per_cell": "min_umi_per_cell",
    "max_umi_per_cell": "max_umi_per_cell",
    "min_genes_per_cell": "min_genes_per_cell",
    "max_genes_per_cell": "max_genes_per_cell",
    "min_cell": "min_cell",
    "max_mt_percent": "max_mt_percent",
    "doublet_method": "doublet_method",
    "scrublet_cutoff": "scrublet_cutoff",
    # Analysis
    "norm_target_sum": "norm_target_sum",
    "n_top_genes": "n_top_genes",
    "regress_out": "regress_out",
    "scale_max_value": "scale_max_value",
    "n_pcs": "n_pcs",
    "batch_correction": "batch_correction",
    "batch_vars": "batch_vars",
    "n_neighbors": "n_neighbors",
    "resolution": "resolution",
    "compute_tsne": "compute_tsne",
    # Annotation
    "annotation_method": "annotation_method",
    "sci_model_path": "sci_model_path",
    "cty_model_path": "cty_model_path",
    "cty_model_name": "cty_model_name",
    # Differential expression
    "pval_threshold": "pval_threshold",
    "logfc_threshold": "logfc_threshold",
    "dea_method": "dea_method",
    "pts_threshold": "pts_threshold",
    "top_n_deg_leidn": "top_n_deg_leidn",
    "top_n_deg_scim": "top_n_deg_scim",
    "top_n_deg_cltpst": "top_n_deg_cltpst",
    # Additional
    "doc_url": "doc_url",
    "data_url": "data_url",
    "only_qc": "only_qc",
    "fix_gene_names": "fix_gene_names",
    "plot_alpha": "plot_alpha",
}

REQUIRED_KEYS = {"input", "project", "species", "tissue", "disease"}


def build_config(args: argparse.Namespace) -> dict[str, Any]:
    """Turn the parsed CLI args into a CellExpress-shaped config dict,
    keeping only the parameters that were actually supplied.
    """
    config: dict[str, Any] = {}
    for dest, config_key in CONFIG_KEY_MAP.items():
        value = getattr(args, dest, None)
        if value is not None:
            config[config_key] = value

    missing = REQUIRED_KEYS - config.keys()
    if missing:
        raise ValueError(
            f"Missing required parameter(s): {', '.join(sorted(missing))}. "
            "These are required by CellExpress and have no safe default."
        )

    # --input must be an absolute path with no trailing slash, per the
    # CellExpress README's data-governance requirement. We normalize it
    # here rather than trusting the caller, since a subtle path error
    # would otherwise surface as an opaque CellExpress failure later.
    input_path = Path(config["input"]).expanduser().resolve()
    if not input_path.is_dir():
        raise FileNotFoundError(f"--input directory does not exist: {input_path}")
    config["input"] = str(input_path)

    return config


# ---------------------------------------------------------------------------
# Running the pipeline
# ---------------------------------------------------------------------------

def run_cellexpress(
    main_py: Path,
    config_path: Path,
    python_exe: str,
) -> int:
    """Invoke CellExpress as a subprocess and stream its output live.

    We stream rather than capture-then-print so a long-running pipeline
    (these can take a while on real datasets) doesn't look "hung" -- you
    see progress as it happens, the same as running it by hand.
    """
    cmd = [python_exe, str(main_py), "--config", str(config_path)]
    print(f"[run_cellexpress] Executing: {' '.join(cmd)}\n", flush=True)

    process = subprocess.run(cmd)
    return process.returncode


# ---------------------------------------------------------------------------
# Output discovery + summarization
#
# WHY: a subprocess exiting 0 tells you CellExpress didn't crash. It does
# NOT tell you the biology came out sane (empty clusters, everything
# annotated "unknown", QC that dropped 95% of cells, etc). This is the
# "verify outputs before calling the stage done" step from CLAUDE.md,
# automated so it happens the same way on every run.
# ---------------------------------------------------------------------------

def find_latest_output_dir(search_root: Path) -> Path | None:
    """Find the most recently created outputs_cellexpress_v*_* folder.

    CellExpress writes its output folder inside (or relative to) the run
    location with a version + run-UID suffix we can't predict in advance,
    so we search for it by name pattern and pick the newest match.
    """
    candidates = sorted(
        search_root.glob("outputs_cellexpress_v*"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    return candidates[0] if candidates else None


def summarize_output(output_dir: Path) -> dict[str, Any]:
    """Read back what CellExpress actually produced and summarize it.

    Kept dependency-light: only pulls in anndata if it's available, and
    degrades gracefully (skips the h5ad-derived stats) if not, since the
    config/report files alone still convey useful information.
    """
    summary: dict[str, Any] = {"output_dir": str(output_dir)}

    html_reports = list(output_dir.glob("report_cellexpress_v*.html"))
    config_snapshots = list(output_dir.glob("config_cellexpress_v*.json"))
    qced_files = list(output_dir.glob("counts-qced_cellexpress_v*.h5ad"))
    adata_files = list(output_dir.glob("adata_cellexpress_v*.h5ad"))

    summary["html_report"] = str(html_reports[0]) if html_reports else None
    summary["config_snapshot"] = str(config_snapshots[0]) if config_snapshots else None

    if config_snapshots:
        with open(config_snapshots[0]) as f:
            summary["config_used"] = json.load(f)

    try:
        import anndata as ad  # local import: keep this an optional dependency
    except ImportError:
        summary["h5ad_summary"] = (
            "anndata not installed in this environment -- skipped reading "
            ".h5ad outputs. Install anndata to get cluster/annotation stats."
        )
        return summary

    h5ad_summary: dict[str, Any] = {}

    if qced_files:
        # backed='r' avoids loading the full matrix into memory just to
        # read shape/obs -- important since these files can be large.
        qced = ad.read_h5ad(qced_files[0], backed="r")
        h5ad_summary["qced_counts"] = {
            "file": str(qced_files[0]),
            "n_cells": qced.n_obs,
            "n_genes": qced.n_vars,
        }

    if adata_files:
        adata = ad.read_h5ad(adata_files[0], backed="r")
        adata_info: dict[str, Any] = {
            "file": str(adata_files[0]),
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
        }

        # Cluster sizes, if Leiden clustering ran.
        for col in ("leiden", "cluster"):
            if col in adata.obs.columns:
                adata_info["cluster_sizes"] = (
                    adata.obs[col].value_counts().to_dict()
                )
                break

        # Annotation coverage: fraction of cells with a non-null label,
        # for whichever annotation method(s) were requested.
        annotation_cols = [
            c for c in adata.obs.columns
            if "scimilarity" in c.lower() or "celltypist" in c.lower()
            or "cell_type" in c.lower()
        ]
        if annotation_cols:
            coverage = {}
            for col in annotation_cols:
                non_null = adata.obs[col].notna().sum()
                coverage[col] = round(non_null / adata.n_obs, 4) if adata.n_obs else None
            adata_info["annotation_coverage"] = coverage

        h5ad_summary["annotated_adata"] = adata_info

        # QC attrition: how many cells survived QC relative to the fully
        # processed object. Both files come out of the same pipeline run
        # so a mismatch here usually points to a downstream filtering
        # step worth double-checking, not a bug in this script.
        if qced_files:
            before = h5ad_summary["qced_counts"]["n_cells"]
            after = adata_info["n_cells"]
            h5ad_summary["qc_attrition"] = {
                "cells_before_further_processing": before,
                "cells_after_processing": after,
                "retained_fraction": round(after / before, 4) if before else None,
            }

    summary["h5ad_summary"] = h5ad_summary
    return summary


def print_report(summary: dict[str, Any]) -> None:
    print("\n" + "=" * 70)
    print("CellExpress run summary")
    print("=" * 70)
    print(f"Output directory: {summary.get('output_dir')}")
    print(f"HTML report:      {summary.get('html_report')}")
    print(f"Config snapshot:  {summary.get('config_snapshot')}")

    h5ad = summary.get("h5ad_summary")
    if isinstance(h5ad, str):
        print(f"\n{h5ad}")
        return
    if not h5ad:
        print("\nNo .h5ad outputs found to summarize.")
        return

    if "qced_counts" in h5ad:
        qc = h5ad["qced_counts"]
        print(f"\nQC'd counts matrix: {qc['n_cells']} cells x {qc['n_genes']} genes")

    if "annotated_adata" in h5ad:
        info = h5ad["annotated_adata"]
        print(f"Annotated object:   {info['n_cells']} cells x {info['n_genes']} genes")

        if "cluster_sizes" in info:
            print("\nCluster sizes:")
            for cluster, n in sorted(info["cluster_sizes"].items(), key=lambda kv: str(kv[0])):
                print(f"  {cluster}: {n}")

        if "annotation_coverage" in info:
            print("\nAnnotation coverage (fraction of cells labeled):")
            for col, frac in info["annotation_coverage"].items():
                print(f"  {col}: {frac}")

    if "qc_attrition" in h5ad:
        attr = h5ad["qc_attrition"]
        print(
            f"\nCells retained through processing: "
            f"{attr['cells_after_processing']} / {attr['cells_before_further_processing']} "
            f"({attr['retained_fraction']:.1%})"
            if attr["retained_fraction"] is not None
            else "\nCells retained through processing: unavailable"
        )

    print("=" * 70 + "\n")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Wrapper: build a CellExpress config, run it, and summarize outputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- Where CellExpress lives / how to invoke it -----------------------
    p.add_argument("--cellexpress-path", default=None,
                    help="Path to cellexpress/main.py. Auto-detected if omitted.")
    p.add_argument("--python-exe", default=sys.executable,
                    help="Python executable to run CellExpress with.")
    p.add_argument("--config-in", default=None,
                    help="Path to an existing, already-reviewed config JSON. "
                         "If given, all other parameter flags are ignored.")
    p.add_argument("--config-out", default=None,
                    help="Where to write the assembled config JSON. "
                         "Defaults to a timestamped file in the input directory.")

    # --- Modes --------------------------------------------------------------
    p.add_argument("--dry-run", action="store_true",
                    help="Build and print the config, but do not execute CellExpress. "
                         "Use this for the human-approval step before a real run.")
    p.add_argument("--only-qc", dest="only_qc", action="store_const", const="yes", default=None,
                    help="Run CellExpress in QC-only mode (no downstream processing).")

    # --- General (required unless --config-in is used) ----------------------
    p.add_argument("--input", help="Absolute path to input directory (sample folders + metadata.csv).")
    p.add_argument("--project", help="Project name, used to tag outputs.")
    p.add_argument("--species", choices=["hs", "mm"], help="hs (human) or mm (mouse).")
    p.add_argument("--tissue", help="Tissue name.")
    p.add_argument("--disease", help="Disease name.")

    # --- QC -------------------------------------------------------------
    p.add_argument("--min-umi-per-cell", dest="min_umi_per_cell", type=int)
    p.add_argument("--max-umi-per-cell", dest="max_umi_per_cell", type=int)
    p.add_argument("--min-genes-per-cell", dest="min_genes_per_cell", type=int)
    p.add_argument("--max-genes-per-cell", dest="max_genes_per_cell", type=int)
    p.add_argument("--min-cell", dest="min_cell", type=int)
    p.add_argument("--max-mt-percent", dest="max_mt_percent", type=float)
    p.add_argument("--doublet-method", dest="doublet_method", choices=["scrublet"])
    p.add_argument("--scrublet-cutoff", dest="scrublet_cutoff", type=float)

    # --- Analysis ---------------------------------------------------------
    p.add_argument("--norm-target-sum", dest="norm_target_sum", type=float)
    p.add_argument("--n-top-genes", dest="n_top_genes", type=int)
    p.add_argument("--regress-out", dest="regress_out", choices=["yes", "no"])
    p.add_argument("--scale-max-value", dest="scale_max_value", type=float)
    p.add_argument("--n-pcs", dest="n_pcs", type=int)
    p.add_argument("--batch-correction", dest="batch_correction", choices=["harmony", "scvi"])
    p.add_argument("--batch-vars", dest="batch_vars",
                    help="Comma-separated obs columns (e.g. donor_id,sample_id).")
    p.add_argument("--n-neighbors", dest="n_neighbors", type=int)
    p.add_argument("--resolution", dest="resolution", type=float)
    p.add_argument("--compute-tsne", dest="compute_tsne", choices=["yes", "no"])

    # --- Annotation ---------------------------------------------------------
    p.add_argument("--annotation-method", dest="annotation_method",
                    help="Comma-separated: scimilarity,celltypist")
    p.add_argument("--sci-model-path", dest="sci_model_path")
    p.add_argument("--cty-model-path", dest="cty_model_path")
    p.add_argument("--cty-model-name", dest="cty_model_name")

    # --- Differential expression --------------------------------------------
    p.add_argument("--pval-threshold", dest="pval_threshold", type=float)
    p.add_argument("--logfc-threshold", dest="logfc_threshold", type=float)
    p.add_argument("--dea-method", dest="dea_method", choices=["wilcoxon", "t-test", "logreg"])
    p.add_argument("--pts-threshold", dest="pts_threshold", type=float)
    p.add_argument("--top-n-deg-leidn", dest="top_n_deg_leidn", type=int)
    p.add_argument("--top-n-deg-scim", dest="top_n_deg_scim", type=int)
    p.add_argument("--top-n-deg-cltpst", dest="top_n_deg_cltpst", type=int)

    # --- Additional ---------------------------------------------------------
    p.add_argument("--doc-url", dest="doc_url")
    p.add_argument("--data-url", dest="data_url")
    p.add_argument("--fix-gene-names", dest="fix_gene_names")
    p.add_argument("--plot-alpha", dest="plot_alpha", type=float)

    return p


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    # ---- Assemble (or load) the config -----------------------------------
    if args.config_in:
        config_path = Path(args.config_in).expanduser().resolve()
        if not config_path.is_file():
            print(f"error: --config-in file not found: {config_path}", file=sys.stderr)
            return 2
        with open(config_path) as f:
            config = json.load(f)
        print(f"[run_cellexpress] Using existing config: {config_path}")
    else:
        try:
            config = build_config(args)
        except (ValueError, FileNotFoundError) as e:
            print(f"error: {e}", file=sys.stderr)
            return 2

        config_out = (
            Path(args.config_out).expanduser().resolve()
            if args.config_out
            else Path(config["input"]) / f"config_review_{int(time.time())}.json"
        )
        config_out.write_text(json.dumps(config, indent=2))
        config_path = config_out
        print(f"[run_cellexpress] Wrote config for review: {config_path}")

    # ---- Always show the config that will be used --------------------------
    print("\nConfig to be used:")
    print(json.dumps(config, indent=2))

    if args.dry_run:
        print(
            "\n[run_cellexpress] --dry-run: stopping here. "
            "Nothing was executed. Review the config above/in the file, "
            "then re-run without --dry-run to proceed."
        )
        return 0

    # ---- Resolve CellExpress location and run -------------------------------
    try:
        main_py = find_cellexpress_main(args.cellexpress_path)
    except FileNotFoundError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    input_dir = Path(config["input"])
    returncode = run_cellexpress(main_py, config_path, args.python_exe)

    if returncode != 0:
        print(f"\n[run_cellexpress] CellExpress exited with code {returncode}.", file=sys.stderr)
        return returncode

    # ---- Find and summarize outputs -----------------------------------------
    # CellExpress writes its output folder relative to where it's run / the
    # input directory; we search the input directory first since that's the
    # documented convention, and fall back to the current working directory.
    output_dir = find_latest_output_dir(input_dir) or find_latest_output_dir(Path.cwd())
    if output_dir is None:
        print(
            "\n[run_cellexpress] CellExpress finished successfully, but no "
            "'outputs_cellexpress_v*' folder was found to summarize. Check "
            "the input directory and CellExpress's own logs manually.",
            file=sys.stderr,
        )
        return 0

    summary = summarize_output(output_dir)
    print_report(summary)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
