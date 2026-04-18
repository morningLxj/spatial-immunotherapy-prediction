# Environment Setup

This document summarizes a practical setup path for running the curated public repository workflows.

## Python

Use Python `3.10+` if possible.

```bash
python -m venv .venv
.venv\Scripts\activate
pip install --upgrade pip
pip install -r requirements.txt
```

If you prefer Conda:

```bash
conda create -n spatial-immuno python=3.10
conda activate spatial-immuno
pip install -r requirements.txt
```

## R

Install R `4.2+` and then add the MR-related packages:

```r
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
install.packages(c("tidyverse", "data.table"))
```

## Recommended Python Packages

The public repository mainly relies on:

- `numpy`
- `pandas`
- `scipy`
- `scikit-learn`
- `xgboost`
- `scanpy`
- `anndata`
- `lifelines`
- `matplotlib`
- `seaborn`
- `python-docx`

## Workspace Assumptions

Several scripts expect a local workspace root that contains source tables, expression matrices, and intermediate analysis files.

Typical examples include:

- `integrated_data.csv`
- `multi_omics_feature_importance.csv`
- `gse31210_clinical.csv`
- `GSE31210_series_matrix.txt.gz`
- `gpl570_probe_mapping.csv`
- local manuscript table source files

When possible, run scripts with explicit arguments such as `--root` and `--out` rather than editing paths inside the code.

## Output Policy

Do not commit large regenerated outputs by default.

Write local outputs into:

- `results/figures/`
- `results/tables/`
- `results/audits/`

or into a separate private workspace directory.

## Verification

Before running the reporting pipeline, confirm:

- Python imports succeed without errors
- R can load `TwoSampleMR`
- required local input files are present
- the repository root is the current working directory when relative commands are used
