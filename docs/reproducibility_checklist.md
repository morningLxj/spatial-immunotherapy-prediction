# Reproducibility Checklist

This checklist summarizes the minimum items needed to understand and rerun the public repository workflows.

## Environment

- Python environment created and dependencies installed from `requirements.txt`
- R environment available for MR-related scripts
- `TwoSampleMR` installed in R
- Local file paths adjusted if scripts are run outside the original workspace layout

## Required Input Data

- `integrated_data.csv`
- `multi_omics_feature_importance.csv`
- `gse31210_clinical.csv`
- `GSE31210_series_matrix.txt.gz`
- `gpl570_probe_mapping.csv`
- local table source files under the manuscript workspace
- spatial objects or exported matrices required by Visium and CosMx analysis steps

## Core Reproducible Outputs

- Main and supplementary tables via `code/07_reporting/rebuild_tables_and_docs.py`
- repeated-CV XGBoost AUC distribution via `code/05_validation/append_xgboost_auc_distribution.py`
- figure rebuild outputs for Figures 1-7 via `code/06_plotting/rebuild_figure*.py`

## Reporting Consistency Checks

- Table 1 smoking explanation appears in the table note logic
- Table 3 uses the updated "Integrated Genetic, Spatial, and Clinical Evidence" title
- Supplementary Table 10 uses journal-style column names
- figure wording is aligned with the latest conservative manuscript framing

## Before Sharing Results

- confirm no patient-level raw data are staged in Git
- confirm generated manuscript binaries are excluded
- confirm exported figures and tables are written to local output folders rather than committed by default
- confirm manuscript text and repository summary use the same model, table, and figure terminology
