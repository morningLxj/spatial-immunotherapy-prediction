# Run Order

This document gives a practical execution order for the curated repository workflows.

## 1. Prepare the environment

Follow [environment_setup.md](environment_setup.md) and confirm that both Python and R dependencies are available.

## 2. Confirm local inputs

Before launching any major script, verify that your local workspace contains:

- integrated feature matrix and response/survival data
- ranked feature list
- external validation cohort files
- local table source files
- spatial analysis inputs used by the plotting and reporting scripts

## 3. Rebuild reporting tables

This is the most useful starting point if your goal is to regenerate the latest manuscript-facing tables.

```bash
python code/07_reporting/rebuild_tables_and_docs.py --root <workspace_root> --out <output_dir>
```

Expected outputs include:

- main table CSV files
- supplementary table CSV files
- `Main_Tables.xlsx`
- `Supplementary_Tables.xlsx`
- `Main_Tables.docx`
- `Supplementary_Tables.docx`

## 4. Refresh the XGBoost AUC distribution

Use this when you need the repeated-CV distribution and summary files for the final selected model.

```bash
python code/05_validation/append_xgboost_auc_distribution.py --root <workspace_root> --out-dir <output_dir>
```

## 5. Rebuild manuscript figures

Run the figure scripts individually as needed:

```bash
python code/06_plotting/rebuild_figure1.py
python code/06_plotting/rebuild_figure2.py
python code/06_plotting/rebuild_figure3.py
python code/06_plotting/rebuild_figure4.py
python code/06_plotting/rebuild_figure5.py
python code/06_plotting/rebuild_figure6.py
python code/06_plotting/rebuild_figure7.py
```

These scripts are most useful when you want the final wording-aligned versions used in the latest submission package.

## 6. Run MR-specific or legacy analysis scripts only if needed

The repository still contains earlier-stage or component-specific scripts under:

- `code/02_feature_selection/`
- `code/03_mendelian_randomization/`
- `code/04_spatial_analysis/`
- `code/05_validation/`

Use these when you need deeper analysis components rather than just the submission-facing rebuild path.

## 7. Apply manuscript formatting

If you have a local DOCX manuscript and want to standardize the formatting:

```bash
python code/07_reporting/format_submission_docx.py <path_to_docx>
```

## 8. Final consistency check

Before sharing outputs, confirm:

- figure wording matches the latest manuscript language
- Table 1 note logic is preserved
- Table 3 title uses the updated safer wording
- Supplementary Table 10 uses the journal-style column names
- local output directories, not Git, hold regenerated binaries
