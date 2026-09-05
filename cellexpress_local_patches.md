# Local patches to the CellExpress clone

The CellExpress pipeline is cloned at `~/cellatria` (AstraZeneca/cellatria,
`--depth 1` of `main`). The following **local divergences from upstream** are
required to run it on this machine. Re-apply them after any fresh clone / pull.

## 1. pandoc stack-size off-by-one (blocks all HTML report rendering)

**Files:**
- `~/cellatria/cellexpress/report_onlyqc_wrapper.py:80`
- `~/cellatria/cellexpress/report_wrapper.py:200`

**Change:** `options(pandoc.stack.size = '4096m')` → `'2048m'`

**Why:** Upstream hardcodes `'4096m'` = 4,294,967,296 bytes, which is exactly
**1 byte over** the max RTS `-K` value accepted by the local pandoc 2.12
(`4,294,967,295` = 2^32 − 1). pandoc aborts with
`error in RTS option -K4096m: size outside allowed range (0 - 4294967295)`,
which fails the R Markdown render and makes CellExpress exit 1 — even though
all QC/analysis computation succeeded. `2048m` is ample for report rendering
and well under the limit. Worth reporting upstream to AstraZeneca.

**Environment:** pandoc 2.12 (anaconda), R 4.5.1, macOS arm64.

## Reproducibility notes
- Python deps for the pipeline were installed into the `thesis_agent_1` conda
  env; a full freeze is in `cellexpress_requirements.txt` (repo root).
- Extra R package installed beyond a base R: `visNetwork` (see
  `cellexpress_R_extra.txt`).
