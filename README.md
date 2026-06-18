# C1Q/SPP1 Macrophage Spatial Context in NSCLC

Public repository for the manuscript:

**C1Q/SPP1 Macrophage Spatial Context as a Reproducible Traceability Feature in Non-Small Cell Lung Cancer: A Public Multi-Cohort Cancer-Informatics Evidence Synthesis**

This repository supports a conservative cancer-informatics evidence synthesis. The central claim is that the C1Q/SPP1 axis marks a reproducible macrophage-associated spatial context in public NSCLC resources. The repository is not intended to support clinical deployment, treatment selection, therapeutic targeting, or definitive causal mechanism claims.

## Repository Contents

| Path | Contents |
|---|---|
| `manuscript/` | Strictly revised manuscript in LaTeX, PDF, and Word format |
| `supplementary/` | Strictly revised supplementary material in LaTeX, PDF, and Word format |
| `figures/` | Main and supplementary figures used by the manuscript |
| `data/` | Source-data manifest, evidence-chain map, and machine-readable source-data archive |
| `code/` | Lightweight checks for source-data completeness and manuscript consistency |
| `docs/` | Data/code availability, reporting boundary notes, and reviewer-facing documentation |

## Evidence Boundary

The analyses are interpreted as hypothesis-generating and traceability-oriented:

- MR analyses are used for gene prioritization under instrumental-variable assumptions, not as proof of mechanism.
- Counterfactual C1QA attenuation is a scoring-sensitivity analysis, not an experimental perturbation.
- External survival and immunotherapy-response cohorts are retrospective and limited; they are not clinical validation sets.
- Spatial statistics describe context and association; they do not prove direct cell-cell signaling.
- No locked clinical threshold, clinical assay, therapeutic target, or treatment recommendation is provided.

## Data Sources

The manuscript analyzes public resources including:

- TCGA NSCLC data from the NCI Genomic Data Commons.
- GEO datasets including GSE31210, GSE91061, GSE126044, and GSE135222.
- eQTL resources including eQTL Catalogue/eQTLGen and GTEx context.
- Human Protein Atlas immunohistochemistry images.
- Public spatial transcriptomics resources summarized in the manuscript and source-data tables.

See `data/source_data_manifest.csv` and `data/evidence_chain_map.csv` for claim-to-source mapping.

## Reproducibility

The file `data/S1_source_data_and_machine_readable_tables.zip` contains the machine-readable tables used for source-data checks. To inspect the archive and generate a source-data inventory:

```bash
Rscript code/check_source_data.R
```

This script verifies the expected archive exists, lists its contents, and writes `release/source_data_inventory.csv`. It is a repository integrity check, not a full re-analysis pipeline from raw public data.

## Manuscript Builds

The manuscript and supplement were compiled locally with TeX Live 2025:

```bash
cd manuscript
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error main_manuscript.tex

cd ../supplementary
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error supplementary_material.tex
```

The LaTeX source uses relative paths to `figures/`. If compiling from a different working directory, preserve the repository structure.

## Citation

If this repository is used before journal publication, cite the repository URL and manuscript title. After publication, replace this section with the final DOI and archived repository DOI.

## License

Code and repository documentation are released under the MIT License. Manuscript files, figures, and processed public-data derivatives are provided for scholarly review and should be reused according to the final journal license and original data-provider terms. See `LICENSES.md`.
