# Data and Code Availability

All primary datasets analyzed in the manuscript are public. This repository provides processed source-data tables, figures, manuscript files, and traceability documentation needed for reviewer inspection.

## Public Data Sources

- TCGA NSCLC: NCI Genomic Data Commons.
- GEO: GSE31210, GSE91061, GSE126044, and GSE135222.
- eQTL summary statistics: eQTL Catalogue and eQTLGen.
- GTEx v8: GTEx Portal.
- Protein-expression context: Human Protein Atlas. Figure 8 uses public HPA context only; image-source and license notes are tracked in `data/hpa_image_sources.csv`.

## Repository Data Files

- `data/source_data_manifest.csv`: package-level source-data manifest.
- `data/evidence_chain_map.csv`: claim-to-source mapping with allowed conclusions and inference boundaries.
- `data/hpa_image_sources.csv`: HPA source-page and license notes for public protein-expression context.
- `data/S1_source_data_and_machine_readable_tables.zip`: machine-readable source-data tables.

## Code Status

This repository currently provides source-data integrity checks rather than a complete raw-data-to-figure pipeline. The check script is:

- `code/check_source_data.R`

The public-data acquisition and complete re-analysis workflow should be archived as a versioned release before final publication if required by the journal or reviewers.

## Inference Boundary

The repository does not provide a clinical model package, locked clinical cutoff, treatment-selection tool, therapeutic-target validation workflow, or author-generated IHC validation dataset.
