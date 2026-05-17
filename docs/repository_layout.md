# Repository Layout

This repository is a curated publication-facing subset of the broader local analysis workspace.

## Design Principles

- Keep code that is needed to understand or reproduce the final analysis logic.
- Exclude large raw datasets, manuscript binaries, generated figures, and one-off troubleshooting files.
- Group scripts by analytical stage instead of by local working-history order.
- Preserve the final reporting path used for the submission-ready tables and manuscript assets.

## Directory Map

### `code/02_feature_selection/`

Core feature-selection and nested cross-validation scripts.

### `code/03_mendelian_randomization/`

MR scripts used to prioritize genes with directionally informative genetic support.

### `code/04_spatial_analysis/`

Spatial preprocessing and Visium-related analysis utilities.

### `code/05_validation/`

External validation and pan-cancer summary scripts, including the targeted XGBoost AUC distribution updater added during final reporting cleanup.

### `code/06_plotting/`

Final figure rebuild scripts used to synchronize manuscript figures with the latest wording and panel logic.

### `code/07_reporting/`

Final reporting scripts, including:

- unified main/supplementary table generation
- export to CSV, XLSX, and DOCX
- submission-style manuscript formatting helpers

### `docs/`

Repository-facing documentation for readers and collaborators.

### `results/`

Optional local output folder for regenerated figures, tables, and audit files. This directory is scaffolded for structure but generated outputs are not committed by default.

## Public Release Alignment

For manuscript-facing navigation, this repository is also exposed through six top-level public release entry files:

- `01_data_accession_manifest.md`
- `02_processing_and_qc.md`
- `03_feature_selection_and_genetic_prioritization.md`
- `04_spatial_and_consistency_analysis.md`
- `05_reporting_assets.md`
- `06_session_info.md`

See [public_release_layout.md](public_release_layout.md) for the mapping.

## Local Workspace Mapping

The repository was synchronized against a larger local workspace containing:

- manuscript drafting files
- regenerated journal figures
- submission package exports
- audit reports and validation notes

Only the stable and reusable pieces are intentionally mirrored here.

## Notes for Maintainers

- If a script hardcodes local paths, prefer converting it to CLI arguments before treating it as a reusable workflow entry point.
- If a script was used only for a one-time forensic check, keep it out of the repository unless it documents a recurring validation step.
- When manuscript wording changes, update both the figure/table scripts and the README summary so the public repository stays aligned with the final submission package.
