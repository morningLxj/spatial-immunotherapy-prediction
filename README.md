# Public Reproducibility Release for a Spatial Immune Pattern Study in NSCLC

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-under_review-blue)](README.md)
[![Release](https://img.shields.io/badge/public_release-PLOS%20ONE%20aligned-4c956c)](README.md)

This repository is a curated public release aligned to the current PLOS ONE manuscript submission. It is intended to support transparency, data-source tracing, analytical navigation, and reproducibility-oriented review of the study workflow.

The repository focuses on:

- public data accession and processing notes
- robust feature selection and genetic prioritization
- spatial pattern analysis
- computational consistency analysis
- reporting assets and reproducibility support

This repository should be read as a manuscript support package rather than a clinical deployment or mechanistic proof package.

## Public Release Navigation

- [01 Data Accession Manifest](01_data_accession_manifest.md)
- [02 Processing and QC](02_processing_and_qc.md)
- [03 Feature Selection and Genetic Prioritization](03_feature_selection_and_genetic_prioritization.md)
- [04 Spatial and Consistency Analysis](04_spatial_and_consistency_analysis.md)
- [05 Reporting Assets](05_reporting_assets.md)
- [06 Session Info](06_session_info.md)
- [Public Release Layout Mapping](docs/public_release_layout.md)

## Quick Links

- Environment setup: [docs/environment_setup.md](docs/environment_setup.md)
- Recommended run order: [docs/run_order.md](docs/run_order.md)
- Repository layout: [docs/repository_layout.md](docs/repository_layout.md)
- Data notes: [docs/data_notes.md](docs/data_notes.md)
- Submission sync summary: [docs/latest_submission_sync.md](docs/latest_submission_sync.md)
- Reproducibility checklist: [docs/reproducibility_checklist.md](docs/reproducibility_checklist.md)
- Citation metadata: [`CITATION.cff`](CITATION.cff)

## Study Overview

This public release documents a multi-layer analysis workflow spanning data accession, feature selection, genetic prioritization, spatial pattern analysis, computational consistency analysis, and manuscript-facing reporting.

The current public release highlights:

- `95` robust features retained under nested cross-validation
- `198` genes carried forward as Mendelian-randomization-based genetic prioritization results
- spatial pattern analyses centered on C1Q-related and comparator myeloid signals
- computational attenuation analyses used as internal consistency checks
- synchronized figure and table generation for the current manuscript-facing package

Mendelian randomization is presented here as supportive genetic prioritization rather than definitive proof of biological causality. External validation results are presented as supportive transportability evidence rather than clinical deployment claims.

## Public Release View

The repository exposes a manuscript-facing public navigation layer built around six sections:

1. data accession manifest
2. processing and QC
3. feature selection and genetic prioritization
4. spatial and consistency analysis
5. reporting assets
6. session info

These sections are documented through the top-level entry files listed above, while the implementation-oriented repository tree remains unchanged.

## Repository Notes

- The physical repository structure remains organized as `code/`, `data/`, `docs/`, and `results/`.
- The six public-release entry files provide the manuscript-facing navigation layer used for editorial and reviewer access.
- Large source datasets, manuscript binaries, and private workspace outputs are intentionally excluded.
- This repository is a curated code-and-documentation layer, not a mirror of the entire private analysis workspace.

## Implementation-Oriented Structure

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

See [docs/public_release_layout.md](docs/public_release_layout.md) for the manuscript-facing mapping between this tree and the public release sections.

## Included Code Tracks

- `code/02_feature_selection/`: nested cross-validation and feature-selection utilities
- `code/03_mendelian_randomization/`: scripts supporting genetic prioritization analyses
- `code/04_spatial_analysis/`: spatial preprocessing and analysis utilities
- `code/05_validation/`: transportability-oriented validation and related summary scripts
- `code/06_plotting/`: figure rebuild scripts aligned to the current manuscript wording
- `code/07_reporting/`: main and supplementary table generation plus manuscript-facing formatting helpers

## Quick Start

### Python

```bash
pip install -r requirements.txt
```

### R

The MR-related workflow depends on `TwoSampleMR` and common tidyverse tooling:

```r
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
install.packages(c("tidyverse", "data.table"))
```

For a practical setup guide, see [docs/environment_setup.md](docs/environment_setup.md).

## Reproducibility-Oriented Reporting

For the recommended execution order, see [docs/run_order.md](docs/run_order.md).

### Rebuild the unified tables

```bash
python code/07_reporting/rebuild_tables_and_docs.py --root <workspace_root> --out <output_dir>
```

### Recompute the repeated-CV XGBoost AUC distribution

```bash
python code/05_validation/append_xgboost_auc_distribution.py --root <workspace_root> --out-dir <output_dir>
```

### Rebuild manuscript-facing figures

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
- **Spatial transcriptomics**: public resources compatible with Visium and CosMx-style workflows used in the local analysis workspace

See [docs/data_notes.md](docs/data_notes.md) for expected local inputs and tracking policy.

## Status Notes

- Final table naming and wording are synchronized with the current manuscript-facing submission package.
- Public-facing repository text is being aligned to a conservative PLOS ONE submission framing.
- The preferred entry path for editorial or reviewer browsing is: `README` -> `01` to `06` public entry files -> linked implementation paths.

## Repository Metadata

- Citation metadata: [`CITATION.cff`](CITATION.cff)
- Reproducibility checklist: [docs/reproducibility_checklist.md](docs/reproducibility_checklist.md)

## Citation

If you use this repository, please cite the manuscript version associated with the latest public release and submission package.

## License

This project is licensed under the MIT License.
