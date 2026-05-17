# Public Release Layout

This repository preserves its implementation-oriented directory tree while exposing a manuscript-facing public release view aligned to the current PLOS ONE submission.

## Public Release Sections

### 01_data_accession_manifest

- data source references
- accession identifiers
- public data notes
- Primary links: `docs/data_notes.md`, `01_data_accession_manifest.md`

### 02_processing_and_qc

- preprocessing and QC entry points
- representative processing scripts
- Primary links: `02_processing_and_qc.md`, `code/04_spatial_analysis/prepare_visium.py`

### 03_feature_selection_and_genetic_prioritization

- nested cross-validation and feature selection
- MR-based genetic prioritization
- Primary links: `03_feature_selection_and_genetic_prioritization.md`, `code/02_feature_selection/`, `code/03_mendelian_randomization/`

### 04_spatial_and_consistency_analysis

- spatial pattern analysis
- computational consistency analysis
- Primary links: `04_spatial_and_consistency_analysis.md`, `code/04_spatial_analysis/`, `code/06_plotting/`

### 05_reporting_assets

- figure rebuild scripts
- table rebuild scripts
- reproducibility-facing reporting docs
- Primary links: `05_reporting_assets.md`, `code/06_plotting/`, `code/07_reporting/`

### 06_session_info

- environment notes
- requirements
- reproducibility checklist
- Primary links: `06_session_info.md`, `requirements.txt`, `docs/environment_setup.md`, `docs/reproducibility_checklist.md`

## Physical Repository Structure

The implementation-oriented tree remains:

- `code/`
- `data/`
- `docs/`
- `results/`

This avoids breaking existing script paths while keeping the public navigation consistent with the manuscript.
