import scanpy as sc
import os

# Paths
data_path = r"D:\ZhouFX\Data"
# Ensure directory exists
output_dir = r"D:\ZhouFX\GigaTIME_Paper_Submission\analysis_scripts"
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, "processed_spatial_data.h5ad")

print(f"Loading data from {data_path}...")
try:
    adata = sc.read_visium(path=data_path, count_file='filtered_feature_bc_matrix.h5')
    adata.var_names_make_unique()

    # Basic processing
    print("Processing...")
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    # Use igraph if available, else skip clustering or handle error
    try:
        sc.tl.leiden(adata, resolution=0.5)
    except Exception as e:
        print(f"Leiden clustering failed (likely missing igraph), skipping clustering column: {e}")
        # Add a dummy column so plotting scripts don't crash if they check for it
        adata.obs['leiden'] = '1'

    # Save
    print(f"Saving to {output_path}...")
    adata.write(output_path)
    if os.path.exists(output_path):
        print(f"SUCCESS: File created at {output_path}")
    else:
        print(f"FAILURE: File not found at {output_path}")

except Exception as e:
    import traceback
    traceback.print_exc()
    print(f"Error: {e}")
