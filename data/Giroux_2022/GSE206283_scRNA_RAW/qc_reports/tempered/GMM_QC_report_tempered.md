# Tempered 3D GMM QC — GSE206283 (Giroux 2022)

Generated: 2026-09-21 | R 4.5.1 | mclust 6.1.2

Damaged component rule (ABSOLUTE): mean percent.mt > 10 OR mean nFeature_RNA < 300.
prob_damaged = summed posterior over damaged components; QC_pass = prob_damaged < 0.50.

## Component-note outcomes
                   note     N
                 <char> <int>
1:                   ok    38
2: no_damaged_component     3

## Overall retention
- Raw cells: 472,830
- Passing QC: 410,863 (86.89%)
- Bootstrap 95% CI on cohort retention: [86.80%, 86.99%]
- Per-sample retention: 59.0%–100.0% (median 88.9%)

## Outputs
- sample_qc_summary_tempered.csv, cell_qc_flags_tempered.csv.gz, gmm_qc_tempered_<sample>.png
- Dir: /Users/marikaclark/dark_matter_atx/data/Giroux_2022/GSE206283_scRNA_RAW/qc_reports/tempered

