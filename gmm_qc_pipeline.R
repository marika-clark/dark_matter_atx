#!/usr/bin/env Rscript
# =============================================================================
# gmm_qc_pipeline.R
# Adaptive 3D Gaussian Mixture Model (GMM) QC filtering across the 41
# GSE206283 (Giroux 2022) scRNA-seq samples.
#
# Strategy: NO fixed global thresholds. Per sample, fit a 3D GMM (mclust) on
# transformed QC metrics and flag the dynamically-identified "damaged" mixture
# component by posterior probability.
#
# Metrics (per cell):
#   x1 = log10(nCount_RNA + 1)
#   x2 = log10(nFeature_RNA + 1)
#   x3 = logit(clip(percent.mt/100, [EPS, 1-EPS]))
#
# Damaged component = argmax over components of  z(mean logit_mito) - z(mean log_genes)
# prob_damaged = posterior of that component ; QC_pass = prob_damaged < 0.5
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(mclust); library(ggplot2)
  library(cowplot); library(boot); library(data.table)
})
set.seed(42)

INPUT_DIR <- "/Users/marikaclark/dark_matter_atx/data/Giroux_2022/GSE206283_scRNA_RAW/samples"
OUT_DIR   <- "/Users/marikaclark/dark_matter_atx/data/Giroux_2022/GSE206283_scRNA_RAW/qc_reports"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

EPS      <- 1e-4     # logit clip bound
PROB_CUT <- 0.5      # posterior threshold for "damaged"
G_RANGE  <- 2:4      # mixture components tested
BOOT_R   <- 1000     # bootstrap replicates for retention CI

REP_SAMPLES <- c("GSM6249236_26-0", "GSM6249271_145-14",
                 "GSM6249239_27-2", "GSM6249274_2213")

logit <- function(p) log(p / (1 - p))

boot_ci <- function(pass_vec, R = BOOT_R) {
  if (length(pass_vec) < 20 || length(unique(pass_vec)) < 2) return(c(NA, NA))
  b  <- boot(pass_vec, statistic = function(d, idx) mean(d[idx]), R = R)
  ci <- tryCatch(boot.ci(b, type = "perc")$percent[4:5],
                 error = function(e) c(NA, NA))
  round(100 * ci, 2)
}

meta <- read.csv(file.path(INPUT_DIR, "metadata.csv"), stringsAsFactors = FALSE)
stopifnot("sample" %in% colnames(meta), nrow(meta) == 41)

all_meta <- vector("list", nrow(meta))
summ     <- vector("list", nrow(meta))

for (i in seq_len(nrow(meta))) {
  s   <- meta$sample[i]
  sev <- if ("severity" %in% colnames(meta)) meta$severity[i] else NA
  counts <- tryCatch(Read10X(file.path(INPUT_DIR, s)), error = function(e) NULL)
  if (is.null(counts)) {                       # robustness: skip a bad folder, don't abort
    cat(sprintf("[%2d/41] %-18s  IO_FAILED — skipped\n", i, s))
    summ[[i]] <- data.table(sample_id = s, severity = sev, G_selected = NA_integer_,
      damaged_component = NA_character_, note = "io_failed",
      raw_cells = NA_integer_, kept_cells = NA_integer_, pct_kept = NA_real_,
      pct_kept_lo = NA_real_, pct_kept_hi = NA_real_,
      mean_mt_before = NA_real_, mean_mt_after = NA_real_,
      mean_counts_before = NA_real_, mean_counts_after = NA_real_,
      mean_genes_before = NA_real_, mean_genes_after = NA_real_)
    next
  }
  obj <- CreateSeuratObject(counts = counts, project = s,
                            min.cells = 0, min.features = 0)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  n_mt <- sum(grepl("^MT-", rownames(obj)))    # guard: Ensembl-ID features would match 0
  if (n_mt == 0) cat(sprintf("    WARNING %s: no ^MT- genes matched — percent.mt degenerate\n", s))
  md <- obj@meta.data

  # --- transformed 3D matrix ---
  logc  <- log10(md$nCount_RNA + 1)
  logg  <- log10(md$nFeature_RNA + 1)
  p     <- pmin(pmax(md$percent.mt / 100, EPS), 1 - EPS)     # clip to avoid +/-Inf
  lmito <- logit(p)
  X <- cbind(log_counts = logc, log_genes = logg, logit_mito = lmito)

  # --- fit GMM with graceful fallback ---
  note <- "ok"
  fit  <- tryCatch(Mclust(X, G = G_RANGE, verbose = FALSE), error = function(e) NULL)
  if (is.null(fit)) { fit <- tryCatch(Mclust(X, G = 2, verbose = FALSE),
                                      error = function(e) NULL); note <- "fallback_G2" }

  if (is.null(fit)) {                          # ultimate fallback: keep all cells
    prob_damaged <- rep(0, nrow(md)); comp <- rep(1L, nrow(md))
    damaged_str <- NA_character_; Gsel <- 1L;  note <- "gmm_failed_all_pass"
  } else {
    Gsel  <- fit$G
    comp  <- fit$classification
    means <- fit$parameters$mean               # d x G
    zsc   <- function(v) if (length(v) > 1 && sd(v) > 0) (v - mean(v)) / sd(v) else rep(0, length(v))
    score <- zsc(means["logit_mito", ]) - zsc(means["log_genes", ])
    # damaged = ALL components with above-average damage score; sum their posteriors
    damaged_comps <- which(score > 0)
    if (length(damaged_comps) == 0) damaged_comps <- which.max(score)  # degenerate guard
    damaged_comps <- as.integer(damaged_comps)
    prob_damaged  <- if (length(damaged_comps) == 1) fit$z[, damaged_comps]
                     else rowSums(fit$z[, damaged_comps, drop = FALSE])
    damaged_str   <- paste(damaged_comps, collapse = ";")
    if (Gsel < 2) note <- "selected_1_component"
  }
  qc_pass <- prob_damaged < PROB_CUT
  if (n_mt == 0) note <- paste0(note, "|no_MT")

  # --- write flags back into the Seurat object's meta.data (per spec) ---
  obj$gmm_component <- comp
  obj$prob_damaged  <- prob_damaged
  obj$QC_pass       <- qc_pass

  rec <- data.table(barcode = rownames(md), sample_id = s, severity = sev,
                    nCount_RNA = md$nCount_RNA, nFeature_RNA = md$nFeature_RNA,
                    percent.mt = md$percent.mt, gmm_component = comp,
                    prob_damaged = prob_damaged, QC_pass = qc_pass)
  all_meta[[i]] <- rec

  ci <- boot_ci(as.integer(qc_pass))
  kept <- sum(qc_pass); raw <- nrow(md)
  summ[[i]] <- data.table(
    sample_id = s, severity = sev, G_selected = Gsel,
    damaged_component = damaged_str, note = note,
    raw_cells = raw, kept_cells = kept, pct_kept = round(100 * kept / raw, 2),
    pct_kept_lo = ci[1], pct_kept_hi = ci[2],
    mean_mt_before = round(mean(md$percent.mt), 3),
    mean_mt_after  = if (kept > 0) round(mean(md$percent.mt[qc_pass]), 3) else NA_real_,
    mean_counts_before = round(mean(md$nCount_RNA), 1),
    mean_counts_after  = if (kept > 0) round(mean(md$nCount_RNA[qc_pass]), 1) else NA_real_,
    mean_genes_before  = round(mean(md$nFeature_RNA), 1),
    mean_genes_after   = if (kept > 0) round(mean(md$nFeature_RNA[qc_pass]), 1) else NA_real_)

  # --- diagnostic plots for representative samples ---
  if (s %in% REP_SAMPLES) {
    dfp <- as.data.frame(rec); dfp$component <- factor(dfp$gmm_component)
    base <- function() list(geom_point(size = .3, alpha = .4), scale_x_log10(),
                            theme_cowplot(12), xlab("nFeature_RNA (log10)"), ylab("percent.mt"))
    p1 <- ggplot(dfp, aes(nFeature_RNA, percent.mt, color = component)) + base() +
          ggtitle(sprintf("%s  mclust components (G=%d)", s, Gsel))
    p2 <- ggplot(dfp, aes(nFeature_RNA, percent.mt, color = prob_damaged)) + base() +
          scale_color_viridis_c() + ggtitle("posterior prob_damaged")
    p3 <- ggplot(dfp, aes(nFeature_RNA, percent.mt, color = QC_pass)) + base() +
          scale_color_manual(values = c(`FALSE` = "#C1352B", `TRUE` = "#4C9F70")) +
          ggtitle(sprintf("kept=%d / removed=%d", kept, raw - kept))
    ggsave(file.path(OUT_DIR, paste0("gmm_qc_", s, ".png")),
           plot_grid(p1, p2, p3, ncol = 3), width = 18, height = 5.2, dpi = 150)
  }

  cat(sprintf("[%2d/41] %-18s G=%d damaged={%s} kept=%d/%d (%.1f%%) [%s]\n",
              i, s, Gsel, ifelse(is.na(damaged_str), "NA", damaged_str),
              kept, raw, 100 * kept / raw, note))
  rm(obj, counts, X); invisible(gc(verbose = FALSE))
}

# ===================== write tables =====================
full <- rbindlist(all_meta)
fwrite(full, file.path(OUT_DIR, "cell_qc_flags.csv.gz"))
sdt <- rbindlist(summ)
tot <- data.table(sample_id = "ALL", severity = "", G_selected = NA_integer_,
                  damaged_component = NA_character_, note = "",
                  raw_cells = sum(sdt$raw_cells), kept_cells = sum(sdt$kept_cells),
                  pct_kept = round(100 * sum(sdt$kept_cells) / sum(sdt$raw_cells), 2),
                  pct_kept_lo = NA, pct_kept_hi = NA,
                  mean_mt_before = NA, mean_mt_after = NA,
                  mean_counts_before = NA, mean_counts_after = NA,
                  mean_genes_before = NA, mean_genes_after = NA)
fwrite(rbind(sdt, tot), file.path(OUT_DIR, "sample_qc_summary.csv"))

# ===================== cohort bootstrap CI =====================
cohort_ci <- boot_ci(as.integer(full$QC_pass), R = BOOT_R)

# ===================== markdown report =====================
notes_tab <- sdt[, .N, by = note][order(-N)]
low_ret   <- sdt[pct_kept < 50][order(pct_kept)]
rep <- c(
"# Adaptive 3D GMM QC Filtering — GSE206283 (Giroux 2022)", "",
sprintf("Generated: %s | R %s | mclust %s | Seurat %s",
        as.character(Sys.Date()), getRversion(),
        packageVersion("mclust"), packageVersion("Seurat")), "",
"## 1. Mathematical transformation",
"Per cell, three metrics are modeled jointly:",
"- `x1 = log10(nCount_RNA + 1)`",
"- `x2 = log10(nFeature_RNA + 1)`",
sprintf("- `x3 = logit(p)`, `p = clip(percent.mt/100, [%.0e, 1-%.0e])` (clip avoids +/-Inf)", EPS, EPS),
"`percent.mt` computed with Seurat `PercentageFeatureSet(pattern='^MT-')`.", "",
"## 2. GMM fitting & damaged-component identification",
sprintf("Per sample: `Mclust(X, G=%d:%d)` selects G and covariance model by BIC.",
        min(G_RANGE), max(G_RANGE)),
"Damaged components = ALL k with `score_k = z(mean logit_mito_k) - z(mean log_genes_k) > 0`",
"(z standardized across components; falls back to argmax if none exceed 0),",
"i.e. above-average mito and/or below-average genes.",
sprintf("`prob_damaged` = SUM of posteriors over damaged components; `QC_pass = prob_damaged < %.2f`.", PROB_CUT), "",
"### Component-selection outcomes across 41 samples",
paste(capture.output(print(notes_tab)), collapse = "\n"), "",
sprintf("G selected: %s", paste(sprintf("G=%d:%d samples",
        sort(unique(sdt$G_selected)),
        as.integer(table(factor(sdt$G_selected, sort(unique(sdt$G_selected)))))), collapse = ", ")), "",
"## 3. Overall retention",
sprintf("- Raw cells: **%s**", format(sum(sdt$raw_cells), big.mark=",")),
sprintf("- Cells passing QC: **%s** (%.2f%%)",
        format(sum(sdt$kept_cells), big.mark=","),
        100 * sum(sdt$kept_cells) / sum(sdt$raw_cells)),
sprintf("- Bootstrap 95%% CI on cohort retention (R=%d): [%.2f%%, %.2f%%]",
        BOOT_R, cohort_ci[1], cohort_ci[2]),
sprintf("- Per-sample retention range: %.1f%% – %.1f%% (median %.1f%%)",
        min(sdt$pct_kept), max(sdt$pct_kept), median(sdt$pct_kept)), "",
if (nrow(low_ret) > 0)
  paste0("**Samples with <50% retention (review):**\n",
         paste(sprintf("- %s: %.1f%% kept", low_ret$sample_id, low_ret$pct_kept), collapse = "\n"))
else "All 41 samples retained >=50% of cells.", "",
"## 4. Outputs",
sprintf("- `sample_qc_summary.csv` — per-sample raw/kept counts, %%, bootstrap CI, mean mito/counts/genes before vs after"),
"- `cell_qc_flags.csv.gz` — per-cell gmm_component, prob_damaged, QC_pass",
sprintf("- `gmm_qc_<sample>.png` — 3-panel diagnostics for: %s", paste(REP_SAMPLES, collapse = ", ")),
sprintf("- All under: `%s`", OUT_DIR), "")
writeLines(rep, file.path(OUT_DIR, "GMM_QC_report.md"))
cat("\n", paste(rep, collapse = "\n"), "\n", sep = "")
cat("\n*** GMM QC pipeline complete. Outputs ->", OUT_DIR, "\n")
