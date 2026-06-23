# PLOS ONE Submission Repository

Public repository for the PLOS ONE manuscript version:

**Convergent public-data evidence for C1Q-associated macrophage programs in non-small cell lung cancer: spatial organization and survival traceability**

This branch is the PLOS ONE submission version of the repository. It is maintained separately from the Cancer Informatics version so that journal-specific manuscript files, figure files, statements, and supporting information can evolve without overwriting the other submission track.

## Branch Boundary

- PLOS ONE version: `plos-one-submission`
- Cancer Informatics version: `cancer-informatics-submission`
- Default public line at the time of branching: `main`

The two submission branches share the same project origin but should be edited independently.

## PLOS ONE Upload Package

The directory `PLOS_ONE_Final_Upload_Slim/` contains the PLOS ONE upload-ready files supplied for this branch:

| Path | Contents |
|---|---|
| `PLOS_ONE_Final_Upload_Slim/Manuscript.docx` | PLOS ONE manuscript file |
| `PLOS_ONE_Final_Upload_Slim/Cover_Letter.docx` | PLOS ONE cover letter |
| `PLOS_ONE_Final_Upload_Slim/Author_Statements.docx` | Author contributions, ethics, funding, competing interests, and AI-use statement |
| `PLOS_ONE_Final_Upload_Slim/Data_Availability_Statement.docx` | PLOS ONE data availability statement |
| `PLOS_ONE_Final_Upload_Slim/Code_Availability_Statement.docx` | Repository and reproducibility statement |
| `PLOS_ONE_Final_Upload_Slim/PLOS_Human_Participants_Research_Checklist.docx` | PLOS human participants checklist |
| `PLOS_ONE_Final_Upload_Slim/Figures/` | Main and supplementary figure files in TIFF format |
| `PLOS_ONE_Final_Upload_Slim/Supplementary Figures.docx` | Supplementary figure document |
| `PLOS_ONE_Final_Upload_Slim/Supplementary_Tables.docx` | Supplementary table document |
| `PLOS_ONE_Final_Upload_Slim/S1_File_Source_Data_and_Machine_Readable_Tables.zip` | Source data and machine-readable tables |

## Repository Structure for PLOS ONE

The PLOS ONE code-availability statement refers to a reproducibility-oriented repository structure. This branch provides that structure as follows:

| Path | Purpose |
|---|---|
| `01_data_accession_manifest/` | Dataset and evidence-source manifests |
| `02_processed_source_tables/` | Machine-readable processed source-data archive |
| `03_analysis_scripts/` | Lightweight verification and reproducibility scripts |
| `04_figure_source_data/` | PLOS ONE figure files and figure-source references |
| `05_reproducibility_log/` | Release notes and source-data inventory |
| `06_session_info/` | Branch provenance and session information notes |

The earlier `manuscript/`, `supplementary/`, `figures/`, `data/`, `code/`, and `docs/` directories remain available for traceability and historical manuscript-source context.

## Evidence Boundary

This PLOS ONE version is framed as a retrospective public-data integration study. The C1Q-associated score is interpreted as a biological traceability measure, not as a deployable clinical prediction tool, treatment-selection rule, therapeutic target, or proof of causal mechanism.

Key interpretation limits:

- The study uses retrospective public datasets.
- No new human participant data were generated.
- No in-house wet-lab validation is claimed.
- External survival findings are treated as traceability and robustness assessment.
- Immunotherapy-associated cohorts are exploratory context only.
- Spatial analyses describe context and organization, not direct cell-cell signaling proof.

## Data Sources

The PLOS ONE manuscript uses public resources including TCGA-LUAD, TCGA-LUSC, GEO survival cohorts, exploratory immunotherapy-associated cohorts, public spatial transcriptomics resources, and public protein-level resources. Dataset provenance and source-data details are provided in the upload package and repository manifests.

## Reproducibility Check

The repository includes a lightweight source-data integrity check:

```bash
Rscript code/check_source_data.R
```

This script verifies the repository source-data archive and writes `release/source_data_inventory.csv`. It is an integrity check for the processed source-data package, not a full re-analysis from raw public repositories.

## Citation

If this repository is used before journal publication, cite the repository URL and PLOS ONE manuscript title. After publication, replace this section with the final DOI and archived repository DOI.

## License

Code and repository documentation are released under the MIT License. Manuscript files, figure files, and processed public-data derivatives are provided for scholarly review and should be reused according to the final journal license and original data-provider terms. See `LICENSES.md`.
