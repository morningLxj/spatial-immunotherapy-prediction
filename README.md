# Convergent public-data evidence for C1Q-associated macrophage programs in non-small cell lung cancer

This repository contains reproducibility materials for a retrospective public-data integration study of C1Q-associated macrophage programs in non-small cell lung cancer (NSCLC). The study evaluates spatial organization, immune-context associations, public proteomic context, comparator-marker analyses and cross-cohort survival traceability using publicly available datasets.

The C1Q-associated score is presented as a biological traceability measure. It is **not** presented as a clinically deployable prognostic model, a treatment-response biomarker, or mechanistic proof. Immunotherapy-associated datasets are used only as exploratory immune-context resources.

Repository version: PLOS ONE public-data integration release, 2026-05-18.

## Repository contents

```text
01_data_accession_manifest/     Public dataset accessions and intended analytical use
02_processed_source_tables/     Machine-readable processed source tables used in the manuscript
03_analysis_scripts/            Lightweight scripts for table checks and source-data summaries
04_figure_source_data/          Figure-source data and source-data manifest
05_reproducibility_log/         Repository notes, limitations and audit report
06_session_info/                Package requirements and session information
```

## How to use this repository

1. Review `01_data_accession_manifest/data_accession_manifest.csv` for the public data resources and their roles.
2. Use the processed tables in `02_processed_source_tables/` and `04_figure_source_data/` to reproduce the reported source-data summaries.
3. Run the lightweight validation scripts:

```bash
python 03_analysis_scripts/01_validate_source_tables.py
python 03_analysis_scripts/02_summarize_external_survival_traceability.py
python 03_analysis_scripts/03_check_repository_language.py
```

The scripts are designed to verify key reported values and repository consistency from the uploaded source tables. They do not re-download raw public datasets.

## Interpretation boundary

This repository follows the conservative interpretation used in the PLOS ONE submission:

- public-data integration study;
- spatial organization and immune-context characterization;
- cross-cohort survival traceability rather than clinical validation;
- exploratory immunotherapy-associated immune-context analyses only;
- no claim of clinical deployment, treatment-response prediction, experimental perturbation, or mechanistic proof.

## Data availability

All source datasets are public. Processed cohort-level derived tables and figure-source files needed to support the reported findings are included here and in the manuscript supporting information.

## Citation

If you use this repository, please cite the associated manuscript when available:

> Convergent public-data evidence for C1Q-associated macrophage programs in non-small cell lung cancer: spatial organization and survival traceability.

## Contact

For questions about this repository, please contact the corresponding author listed in the manuscript.
