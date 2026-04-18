# Data Notes

This repository does not include raw patient-level datasets or large generated outputs.

## Not Tracked in Git

- TCGA downloads
- GEO series matrices
- eQTL resources
- Visium or CosMx raw objects
- generated PDFs, TIFFs, PNGs, and DOCX files
- intermediate CSV exports from local reporting runs

## Expected Local Inputs

Typical local workflows depend on files such as:

- `integrated_data.csv`
- `multi_omics_feature_importance.csv`
- `gse31210_clinical.csv`
- `GSE31210_series_matrix.txt.gz`
- `gpl570_probe_mapping.csv`
- local `spring模板/` table sources
- local spatial objects used for Visium and CosMx analysis

These assets should be stored outside the public repository or regenerated locally from approved sources.

## Public Data Sources

- TCGA NSCLC cohorts from the NCI GDC portal
- GEO datasets including `GSE31210`, `GSE126044`, `GSE135222`, and `GSE91061`
- eQTL resources from the eQTL Catalogue
- public spatial transcriptomics references compatible with Visium and CosMx workflows

## Recommended Practice

- Keep a private local workspace for large data and generated manuscript outputs.
- Use this repository for reusable code, documentation, and lightweight configuration only.
- If a collaborator needs to rerun the pipeline, provide a short environment guide plus a private data manifest rather than committing data snapshots into Git.
