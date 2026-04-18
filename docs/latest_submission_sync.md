# Latest Submission Sync

This repository refresh aligns the public codebase with the latest submission-ready local analysis package.

## Main Content Updates

- README updated from the older "Spatial-Causal" framing to the current C1Q-centered, multi-layer framework language.
- Final prognostic model description updated from earlier AdaBoost-centered wording to the current XGBoost-centered reporting.
- Key counts updated to the latest values:
  - `95` robust features
  - `198` MR-prioritized genes
- Final reporting language now reflects the current main-table and supplementary-table architecture.

## New Curated Additions

- final figure rebuild scripts for Figures 1-7
- shared plotting helper for figure rebuild workflows
- unified table generator for main and supplementary tables
- XGBoost repeated-CV AUC distribution updater
- submission-style DOCX formatter

## Reporting Alignment

The refreshed reporting code is synchronized with the latest local submission package, including:

- Table 1 smoking explanation moved to table note logic
- softened Table 3 title using "Genetic, Spatial, and Clinical" wording
- journal-style Supplementary Table 10 column names

## Repository Policy

This sync intentionally does not add:

- raw patient-level datasets
- generated manuscript binaries
- large figure exports
- one-off logs from local debugging sessions

The goal is to keep the public repository readable, defensible, and close to the manuscript logic without exposing the entire private workspace.
