#!/usr/bin/env Rscript
# =============================================================================
# gmm_qc_pipeline_tempered.R
# TEMPERED variant of the 3D GMM QC (see gmm_qc_pipeline.R for the baseline).
#
# Change vs baseline: the "damaged" component is no longer the relative
# argmax/score>0 rule (which over-removed low-gene but low-mito quiescent
# immune cells). Instead a component is flagged damaged by ABSOLUTE biology:
#     mean(percent.mt) > MITO_FLOOR   OR   mean(nFeature_RNA) < GENE_FLOOR
# GMM still fits per sample; prob_damaged = summed posterior over damaged
# components; QC_pass = prob_damaged < 0.5.  Outputs go to qc_reports/tempered/
# so the baseline results are preserved for benchmarking.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(mclust); library(ggplot2)
  library(cowplot); library(boot); library(data.table)
})
set.seed(42)

INPUT_DIR <- "/Users/marikaclark/dark_matter_atx/data/Giroux_2022/GSE206283_scRNA_RAW/samples"
OUT_DIR   <- "/Users/marikaclark/dark_matter_atx/data/Giroux_2022/GSE206283_scRNA_RAW/qc_reports/tempered"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

EPS        <- 1e-4
PROB_CUT   <- 0.5
G_RANGE    <- 2:4
BOOT_R     <- 1000
MITO_FLOOR <- 10     # component mean percent.mt above this => damaged
GENE_FLOOR <- 300    # component mean nFeature_RNA below this => debris

REP_SAMPLES <- c("GSM6249236_26-0", "GSM6249271_145-14",
                 "GSM6249239_27-2", "GSM6249274_2213")

logit <- function(p) log(p / (1 - p))
boot_ci <- function(pass_vec, R = BOOT_R) {
  if (length(pass_vec) < 20 || length(unique(pass_vec)) < 2) return(c(NA, NA))
  b  <- boot(pass_vec, statistic = function(d, idx) mean(d[idx]), R = R)
  ci <- tryCatch(boot.ci(b, type = "perc")$percent[4:5], error = function(e) c(NA, NA))
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
  if (is.null(counts)) {
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
  obj <- CreateSeuratObject(counts = counts, project = s, min.cells = 0, min.features = 0)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  n_mt <- sum(grepl("^MT-", rownames(obj)))
  if (n_mt == 0) cat(sprintf("    WARNING %s: no ^MT- genes matched\n", s))
  md <- obj@meta.data

  logc  <- log10(md$nCount_RNA + 1)
  logg  <- log10(md$nFeature_RNA + 1)
  p     <- pmin(pmax(md$percent.mt / 100, EPS), 1 - EPS)
  lmito <- logit(p)
  X <- cbind(log_counts = logc, log_genes = logg, logit_mito = lmito)

  note <- "ok"
  fit  <- tryCatch(Mclust(X, G = G_RANGE, verbose = FALSE), error = function(e) NULL)
  if (is.null(fit)) { fit <- tryCatch(Mclust(X, G = 2, verbose = FALSE),
                                      error = function(e) NULL); note <- "fallback_G2" }

  if (is.null(fit)) {
    prob_damaged <- rep(0, nrow(md)); comp <- rep(1L, nrow(md))
    damaged_str <- NA_character_; Gsel <- 1L; note <- "gmm_failed_all_pass"
  } else {
    Gsel <- fit$G
    comp <- fit$classification
    # ---- TEMPERED: absolute biological gates on component means ----
    comp_mt   <- tapply(md$percent.mt,   comp, mean)
    comp_gene <- tapply(md$nFeature_RNA, comp, mean)
    dmg <- which(comp_mt > MITO_FLOOR | comp_gene < GENE_FLOOR)
    damaged_comps <- as.integer(names(dmg))
    if (length(damaged_comps) == 0) {
      prob_damaged <- rep(0, nrow(md)); damaged_str <- "none"; note <- "no_damaged_component"
    } else {
      prob_damaged <- if (length(damaged_comps) == 1) fit$z[, damaged_comps]
                      else rowSums(fit$z[, damaged_comps, drop = FALSE])
      damaged_str  <- paste(damaged_comps, collapse = ";")
    }
    if (Gsel < 2) note <- "selected_1_component"
  }
  qc_pass <- prob_damaged < PROB_CUT

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

  if (s %in% REP_SAMPLES) {
    dfp <- as.data.frame(rec); dfp$component <- factor(dfp$gmm_component)
    base <- function() list(geom_point(size = .3, alpha = .4), scale_x_log10(),
                            theme_cowplot(12), xlab("nFeature_RNA (log10)"), ylab("percent.mt"))
    p1 <- ggplot(dfp, aes(nFeature_RNA, percent.mt, color = component)) + base() +
          ggtitle(sprintf("%s  mclust components (G=%d)", s, Gsel))
    p2 <- ggplot(dfp, aes(nFeature_RNA, percent.mt, color = prob_damaged)) + base() +
          scale_color_viridis_c() + ggtitle("posterior prob_damaged (tempered)")
    p3 <- ggplot(dfp, aes(nFeature_RNA, percent.mt, color = QC_pass)) + base() +
          scale_color_manual(values = c(`FALSE` = "#C1352B", `TRUE` = "#4C9F70")) +
          ggtitle(sprintf("kept=%d / removed=%d", kept, raw - kept))
    ggsave(file.path(OUT_DIR, paste0("gmm_qc_tempered_", s, ".png")),
           plot_grid(p1, p2, p3, ncol = 3), width = 18, height = 5.2, dpi = 150)
  }

  cat(sprintf("[%2d/41] %-18s G=%d damaged={%s} kept=%d/%d (%.1f%%) [%s]\n",
              i, s, Gsel, ifelse(is.na(damaged_str), "NA", damaged_str),
              kept, raw, 100 * kept / raw, note))
  rm(obj, counts, X); invisible(gc(verbose = FALSE))
}

full <- rbindlist(all_meta)
fwrite(full, file.path(OUT_DIR, "cell_qc_flags_tempered.csv.gz"))
sdt <- rbindlist(summ)
tot <- data.table(sample_id = "ALL", severity = "", G_selected = NA,
                  damaged_component = NA_character_, note = "",
                  raw_cells = sum(sdt$raw_cells, na.rm = TRUE),
                  kept_cells = sum(sdt$kept_cells, na.rm = TRUE),
                  pct_kept = round(100 * sum(sdt$kept_cells, na.rm = TRUE) / sum(sdt$raw_cells, na.rm = TRUE), 2),
                  pct_kept_lo = NA_real_, pct_kept_hi = NA_real_,
                  mean_mt_before = NA_real_, mean_mt_after = NA_real_,
                  mean_counts_before = NA_real_, mean_counts_after = NA_real_,
                  mean_genes_before = NA_real_, mean_genes_after = NA_real_)
fwrite(rbind(sdt, tot), file.path(OUT_DIR, "sample_qc_summary_tempered.csv"))

cohort_ci <- boot_ci(as.integer(full$QC_pass), R = BOOT_R)
notes_tab <- sdt[, .N, by = note][order(-N)]
rep <- c(
"# Tempered 3D GMM QC — GSE206283 (Giroux 2022)", "",
sprintf("Generated: %s | R %s | mclust %s", as.character(Sys.Date()), getRversion(), packageVersion("mclust")), "",
sprintf("Damaged component rule (ABSOLUTE): mean percent.mt > %g OR mean nFeature_RNA < %g.", MITO_FLOOR, GENE_FLOOR),
sprintf("prob_damaged = summed posterior over damaged components; QC_pass = prob_damaged < %.2f.", PROB_CUT), "",
"## Component-note outcomes",
paste(capture.output(print(notes_tab)), collapse = "\n"), "",
"## Overall retention",
sprintf("- Raw cells: %s", format(sum(sdt$raw_cells, na.rm = TRUE), big.mark = ",")),
sprintf("- Passing QC: %s (%.2f%%)", format(sum(sdt$kept_cells, na.rm = TRUE), big.mark = ","),
        100 * sum(sdt$kept_cells, na.rm = TRUE) / sum(sdt$raw_cells, na.rm = TRUE)),
sprintf("- Bootstrap 95%% CI on cohort retention: [%.2f%%, %.2f%%]", cohort_ci[1], cohort_ci[2]),
sprintf("- Per-sample retention: %.1f%%–%.1f%% (median %.1f%%)",
        min(sdt$pct_kept, na.rm = TRUE), max(sdt$pct_kept, na.rm = TRUE), median(sdt$pct_kept, na.rm = TRUE)), "",
"## Outputs",
"- sample_qc_summary_tempered.csv, cell_qc_flags_tempered.csv.gz, gmm_qc_tempered_<sample>.png",
sprintf("- Dir: %s", OUT_DIR), "")
writeLines(rep, file.path(OUT_DIR, "GMM_QC_report_tempered.md"))
cat("\n", paste(rep, collapse = "\n"), "\n", sep = "")
cat("\n*** Tempered GMM QC complete. Outputs ->", OUT_DIR, "\n")
