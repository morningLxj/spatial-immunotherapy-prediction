import pandas as pd
import scanpy as sc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import scipy.stats
import scipy.sparse

def get_significance_label(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return ''

def main():
    # Paths
    data_dir = r"D:\ZhouFX\Data"
    h5_path = os.path.join(data_dir, "filtered_feature_bc_matrix.h5")
    existing_csv_path = r"D:\ZhouFX\GigaTIME_Results\Final_Analysis_Correlation_Data.csv"
    output_matrix_path = r"D:\ZhouFX\GigaTIME_Results\Full_Correlation_Matrix.csv"
    # Old paths
    # figure_pdf = r"D:\ZhouFX\论文图表\精选\FigureS6.pdf"
    # figure_png = r"D:\ZhouFX\论文图表\精选\FigureS6.png"

    # New paths consistent with other figures
    os.makedirs("Final_Analysis/Plots", exist_ok=True)
    figure_pdf = "Final_Analysis/Plots/Figure_S6_Hub_Gene_Immune_Correlation.pdf"
    figure_png = "Final_Analysis/Plots/Figure_S6_Hub_Gene_Immune_Correlation.png"

    # 1. Load existing correlation data (contains GigaTIME predictions and barcodes)
    print(f"Loading {existing_csv_path}...")
    df = pd.read_csv(existing_csv_path)
    barcodes = df['Barcode'].tolist()

    # Identify GigaTIME marker columns
    marker_cols = [c for c in df.columns if c.startswith('GigaTIME_')]
    print(f"Found {len(marker_cols)} GigaTIME markers.")

    # 2. Load Raw 10X Data to get Hub Genes expression
    print(f"Loading 10X data from {h5_path}...")
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()

    # Target Genes: Updated to include core genes KCNAB2 and C1QA
    # We strictly focus on the core genes to align with Figure 6 and S5
    target_genes = ['KCNAB2', 'C1QA']

    # Check which genes are present
    available_genes = [g for g in target_genes if g in adata.var_names]
    missing_genes = [g for g in target_genes if g not in adata.var_names]
    print(f"Available genes: {available_genes}")
    print(f"Missing genes: {missing_genes}")

    # 3. Extract expression data
    # Check intersection
    valid_barcodes = [b for b in barcodes if b in adata.obs_names]
    print(f"Barcodes intersection: {len(valid_barcodes)} / {len(barcodes)}")

    if len(valid_barcodes) < len(barcodes) * 0.9:
        print("Warning: Low barcode overlap. Checking format...")

    # Subset adata
    adata = adata[valid_barcodes, available_genes]

    # Create a DataFrame for expression
    if scipy.sparse.issparse(adata.X):
        expr_data = adata.X.toarray()
    else:
        expr_data = adata.X

    df_expr = pd.DataFrame(expr_data, index=valid_barcodes, columns=available_genes)
    df_expr['Barcode'] = valid_barcodes

    # Merge with original df
    cols_to_drop = [g for g in available_genes if g in df.columns]
    df = df.drop(columns=cols_to_drop)

    df_merged = pd.merge(df, df_expr, on='Barcode', how='inner')
    print(f"Merged Data Shape: {df_merged.shape}")

    # 4. Calculate Correlations and P-values
    print("Calculating correlations...")
    correlation_results = []

    for gene in available_genes:
        for marker in marker_cols:
            marker_name = marker.replace('GigaTIME_', '')

            # Pearson correlation and p-value
            r, p = scipy.stats.pearsonr(df_merged[gene], df_merged[marker])

            correlation_results.append({
                'Gene': gene,
                'Marker': marker_name,
                'Correlation': r,
                'P-value': p
            })

    df_corr = pd.DataFrame(correlation_results)

    # Save full matrix data
    df_corr.to_csv(output_matrix_path, index=False)
    print(f"Saved correlation data to {output_matrix_path}")

    # 5. Plot Clustered Heatmap (SCI Style)
    print("Generating SCI-style Heatmap...")

    # Pivot data
    heatmap_data = df_corr.pivot(index='Gene', columns='Marker', values='Correlation')
    p_data = df_corr.pivot(index='Gene', columns='Marker', values='P-value')

    # Ensure specific order for rows (Core genes at top)
    desired_order = ['KCNAB2', 'C1QA']
    existing_order = [g for g in desired_order if g in heatmap_data.index]
    remaining = [g for g in heatmap_data.index if g not in existing_order]
    final_order = existing_order + remaining

    heatmap_data = heatmap_data.reindex(final_order)
    p_data = p_data.reindex(final_order)

    # Create annotation matrix (r value + stars)
    annot_labels = heatmap_data.copy().astype(str)
    for gene in heatmap_data.index:
        for marker in heatmap_data.columns:
            r = heatmap_data.loc[gene, marker]
            p = p_data.loc[gene, marker]
            stars = get_significance_label(p)
            annot_labels.loc[gene, marker] = f"{r:.2f}\n{stars}"

    # Set style
    sns.set(style='white', font='Arial', font_scale=1.2)

    # Use Clustermap for grouping similar markers
    # col_cluster=True will group markers that behave similarly
    # row_cluster=False (since we want to keep our manual priority order)

    # Determine figure size
    figsize = (16, 5) # Compact height for 2 genes

    g = sns.clustermap(heatmap_data,
                       cmap='RdBu_r',
                       center=0,
                       annot=annot_labels,
                       fmt="",  # Raw string format
                       figsize=figsize,
                       dendrogram_ratio=(0.1, 0.2), # Control dendrogram size
                       cbar_pos=(0.02, 0.8, 0.03, 0.15), # Position colorbar
                       cbar_kws={'label': 'Pearson r'},
                       linewidths=1.0,
                       linecolor='white',
                       col_cluster=True,
                       row_cluster=False, # Keep Hub genes manual order
                       yticklabels=True,
                       xticklabels=True)

    # Adjust axes
    g.ax_heatmap.set_xlabel('')
    g.ax_heatmap.set_ylabel('')

    # Rotate labels
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    # Highlight Core Genes labels
    for label in g.ax_heatmap.get_yticklabels():
        if label.get_text() in ['KCNAB2', 'C1QA']:
            label.set_color('#d73027') # Red
            label.set_fontweight('bold')
            label.set_fontsize(14)

    # Title
    g.fig.suptitle('Correlation between Core Genes (KCNAB2, C1QA) and Immune Markers', y=0.98, fontsize=18, fontweight='bold')

    # Save
    plt.savefig(figure_pdf, dpi=600, bbox_inches='tight')
    plt.savefig(figure_png, dpi=600, bbox_inches='tight')
    print(f"Figure saved to {figure_pdf} and {figure_png}")

if __name__ == "__main__":
    main()
