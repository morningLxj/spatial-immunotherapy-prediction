# Spatial-Causal Integration Reveals Complement C1q Hotspots as Drivers of Immunotherapy Response in Non-Small Cell Lung Cancer

## 📖 Abstract

This repository contains the source code and analysis pipeline for the paper **"Spatial-Causal Integration Reveals Complement C1q Hotspots as Drivers of Immunotherapy Response in Non-Small Cell Lung Cancer"**.

We established a comprehensive "Spatial-Causal" framework integrating multi-omics data from 1,100 NSCLC patients (TCGA) with rigorous statistical and spatial validation. To ensure reproducibility and combat overfitting, we employed a LASSO-based feature selection within a 5-fold nested cross-validation (nested CV) framework to screen for robust molecular features. Causal regulatory networks were constructed using Mendelian randomization (MR) with immune cell-specific eQTL data. We further resolved the spatial architecture of identified targets using 10x Genomics Visium spatial transcriptomics and the GigaTIME deep learning model.

## 🚀 Key Findings

- **Robust Biomarkers**: Identified 91 robust genes (including *ACP5* and *KCNAB2*) selected in 100% of nested CV folds.
- **Causal Drivers**: Validated 198 causal genes via Mendelian randomization, linking *KCNAB2* and *C1QA* to immune phenotypes.
- **Spatial Hotspots**: Revealed that C1q family members form distinct "immune hotspots" (Moran's I > 0.35) that co-localize with macrophages.
- **Clinical Prediction**: The derived AdaBoost model significantly outperforms standard clinical baselines (AUC=0.743).

## 📂 Project Structure

```
spatial-immunotherapy-prediction/
├── code/
│   ├── 01_preprocessing/          # Data cleaning and normalization scripts
│   ├── 02_feature_selection/      # Nested Cross-Validation & LASSO implementation
│   ├── 03_mendelian_randomization/# MR analysis using TwoSampleMR (R)
│   ├── 04_spatial_analysis/       # Spatial transcriptomics (Visium) & GigaTIME integration
│   ├── 05_validation/             # External validation (GSE31210) & Pan-cancer analysis
│   └── 06_plotting/               # Scripts to reproduce manuscript figures (Fig 1-8)
├── data/                          # Placeholder for raw/processed data (see Data Availability)
├── results/                       # Output tables and intermediate files
├── requirements.txt               # Python dependencies
└── README.md                      # This file
```

## 🛠️ Installation & Usage

### 1. Prerequisites
- Python 3.8+
- R 4.0+ (for MR analysis)

### 2. Python Environment
```bash
pip install -r requirements.txt
```

### 3. R Dependencies
The Mendelian randomization analysis requires the `TwoSampleMR` package:
```R
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
install.packages("tidyverse")
```

### 4. Running the Analysis

**Feature Selection (Nested CV):**
```bash
# Run LASSO feature selection with nested CV
python code/02_feature_selection/run_nested_cv.py
```

**Mendelian Randomization:**
```R
# Run MR analysis (requires R)
Rscript code/03_mendelian_randomization/mr_full_analysis.R
```

**Spatial Analysis:**
```bash
# Process Visium data and calculate Moran's I
python code/04_spatial_analysis/prepare_visium.py
```

**Reproduce Figures:**
```bash
python code/06_plotting/generate_figure2.py
python code/06_plotting/create_figure7.py
# ... and other scripts in 06_plotting/
```

## 📊 Data Availability

- **TCGA NSCLC**: [NCI GDC Data Portal](https://portal.gdc.cancer.gov)
- **GEO Datasets**: GSE31210, GSE126044, GSE135222 ([NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/))
- **eQTL Data**: [eQTL Catalogue](https://www.ebi.ac.uk/eqtl/)
- **Spatial Data**: 10x Genomics Demonstration Data

## 📝 Citation

If you use this code or findings in your research, please cite:

> **Li X, Zhang F, Zheng X, Xu X, Luo C.** Spatial-Causal Integration Reveals Complement C1q Hotspots as Drivers of Immunotherapy Response in Non-Small Cell Lung Cancer. *Cancer Cell International* (Under Review).

## 📄 License

This project is licensed under the MIT License.
