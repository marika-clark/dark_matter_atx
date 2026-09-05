# Adaptive 3D GMM QC Filtering — GSE206283 (Giroux 2022)

Generated: 2026-09-05 | R 4.5.1 | mclust 6.1.2 | Seurat 5.3.1

## 1. Mathematical transformation
Per cell, three metrics are modeled jointly:
- `x1 = log10(nCount_RNA + 1)`
- `x2 = log10(nFeature_RNA + 1)`
- `x3 = logit(p)`, `p = clip(percent.mt/100, [1e-04, 1-1e-04])` (clip avoids +/-Inf)
`percent.mt` computed with Seurat `PercentageFeatureSet(pattern='^MT-')`.

## 2. GMM fitting & damaged-component identification
Per sample: `Mclust(X, G=2:4)` selects G and covariance model by BIC.
Damaged components = ALL k with `score_k = z(mean logit_mito_k) - z(mean log_genes_k) > 0`
(z standardized across components; falls back to argmax if none exceed 0),
i.e. above-average mito and/or below-average genes.
`prob_damaged` = SUM of posteriors over damaged components; `QC_pass = prob_damaged < 0.50`.

### Component-selection outcomes across 41 samples
     note     N
   <char> <int>
1:     ok    41

G selected: G=4:41 samples

## 3. Overall retention
- Raw cells: **472,830**
- Cells passing QC: **355,297** (75.14%)
- Bootstrap 95% CI on cohort retention (R=1000): [75.02%, 75.26%]
- Per-sample retention range: 52.4% – 99.6% (median 75.9%)

All 41 samples retained >=50% of cells.

## 4. Outputs
- `sample_qc_summary.csv` — per-sample raw/kept counts, %, bootstrap CI, mean mito/counts/genes before vs after
- `cell_qc_flags.csv.gz` — per-cell gmm_component, prob_damaged, QC_pass
- `gmm_qc_<sample>.png` — 3-panel diagnostics for: GSM6249236_26-0, GSM6249271_145-14, GSM6249239_27-2, GSM6249274_2213
- All under: `/Users/marikaclark/dark_matter_atx/data/Giroux_2022/GSE206283_scRNA_RAW/qc_reports`

