# A Causal Inference-Guided Spatial Immune Framework Centered on C1Q in Non-Small Cell Lung Cancer

This repository contains the curated analysis scripts, figure rebuild utilities, and reporting pipeline for our NSCLC study on C1Q-centered spatial immune organization.

The current repository state has been updated to match the latest submission-ready analysis package:

- `95` robust features retained under nested cross-validation
- `198` genes prioritized by Mendelian randomization
- final prognostic model updated to `XGBoost`
- final main and supplementary tables generated from a unified reporting script
- figure wording and table labels synchronized with the latest manuscript language

## Study Overview

We developed a conservative multi-layer framework that integrates:

- nested cross-validation and stability-driven feature discovery
- Mendelian randomization for directionally informative genetic support
- spatial transcriptomics from Visium and CosMx datasets
- in silico perturbation and mediation analysis
- prognostic and exploratory external validation

Rather than relying on a single analysis layer, the project emphasizes convergence across genetic, spatial, and clinical evidence.

## Current Key Results

- **Feature robustness**: `95` features were retained across nested resampling.
- **Genetic prioritization**: `198` genes showed putative immune-related MR support.
- **Spatial organization**: C1Q family genes form structured immune hotspots, while `SPP1` marks a distinct exclusion-associated niche.
- **Perturbation consistency**: computational attenuation of the C1Q axis weakens hotspot organization in a graded manner.
- **Clinical modeling**: `XGBoost` was selected as the final model based on the best cross-validated AUC in the TCGA training cohort.
- **External validation**: the harmonized risk score remained prognostically relevant in `GSE31210`, with exploratory assessment in immunotherapy-related cohorts.

## Repository Layout

```text
spatial-immunotherapy-prediction/
|-- code/
|   |-- 02_feature_selection/
|   |-- 03_mendelian_randomization/
|   |-- 04_spatial_analysis/
|   |-- 05_validation/
|   |-- 06_plotting/
|   `-- 07_reporting/
|-- data/
|   |-- raw/
|   `-- processed/
|-- docs/
|-- results/
|   |-- figures/
|   |-- tables/
|   `-- audits/
|-- .gitignore
|-- LICENSE
|-- README.md
`-- requirements.txt
```

See [docs/repository_layout.md](docs/repository_layout.md) for the mapping between the repository and the latest manuscript-ready workspace assets.

## Included Code Tracks

- `code/02_feature_selection/`: robust feature selection and nested CV utilities
- `code/03_mendelian_randomization/`: MR analysis scripts
- `code/04_spatial_analysis/`: spatial preprocessing and Visium-related analysis
- `code/05_validation/`: external validation and pan-cancer summary scripts
- `code/06_plotting/`: final figure rebuild scripts used to align manuscript figures with the latest wording
- `code/07_reporting/`: final table generation and manuscript formatting scripts

## Quick Start

### Python

```bash
pip install -r requirements.txt
```

### R

The MR workflow depends on `TwoSampleMR` and common tidyverse tooling:

```r
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
install.packages(c("tidyverse", "data.table"))
```

## Reproducing Final Reporting Assets

### Rebuild the unified tables

```bash
python code/07_reporting/rebuild_tables_and_docs.py --root <workspace_root> --out <output_dir>
```

### Recompute the XGBoost repeated-CV AUC distribution

```bash
python code/05_validation/append_xgboost_auc_distribution.py --root <workspace_root> --out-dir <output_dir>
```

### Rebuild manuscript figures

```bash
python code/06_plotting/rebuild_figure1.py
python code/06_plotting/rebuild_figure2.py
python code/06_plotting/rebuild_figure3.py
python code/06_plotting/rebuild_figure4.py
python code/06_plotting/rebuild_figure5.py
python code/06_plotting/rebuild_figure6.py
python code/06_plotting/rebuild_figure7.py
```

## Data Availability

Large source datasets and manuscript output files are not tracked in this repository.

- **TCGA NSCLC**: [NCI GDC Data Portal](https://portal.gdc.cancer.gov)
- **GEO cohorts**: `GSE31210`, `GSE126044`, `GSE135222`, `GSE91061`
- **eQTL resources**: [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/)
- **Spatial transcriptomics**: 10x Genomics Visium and CosMx-compatible resources used in the local analysis workspace

See [docs/data_notes.md](docs/data_notes.md) for expected local paths and tracking policy.

## Status Notes

- This repository is a curated code-and-documentation layer, not a mirror of the entire local workspace.
- Intermediate logs, large result files, manuscript binaries, and private/raw datasets are intentionally excluded.
- Final table naming and wording now match the latest submission package, including the updated Table 3 title and Supplementary Table 10 column labels.
- See [docs/latest_submission_sync.md](docs/latest_submission_sync.md) for a concise summary of the repository refresh.

## Repository Metadata

- Citation metadata: [`CITATION.cff`](CITATION.cff)
- Reproducibility checklist: [docs/reproducibility_checklist.md](docs/reproducibility_checklist.md)

## Citation

If you use this repository, please cite the manuscript version associated with the latest submission package.

> Li X, Zhang F, Zheng X, Xu X, Luo C. A Causal Inference-Guided Spatial Immune Framework Centered on C1Q in Non-Small Cell Lung Cancer. Under review.

## License

This project is licensed under the MIT License.
