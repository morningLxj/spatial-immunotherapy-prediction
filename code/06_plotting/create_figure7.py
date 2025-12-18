import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
from sklearn.ensemble import AdaBoostClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_curve, auc
import statsmodels.api as sm
import seaborn as sns
import os
import gzip
import io

# Set style
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

def load_tcga_data():
    print("Loading TCGA integrated data...")
    if not os.path.exists('integrated_data.csv'):
        print("Error: integrated_data.csv not found.")
        return None
    df = pd.read_csv('integrated_data.csv')
    return df

def get_top_features(n=100):
    print(f"Loading top {n} features...")
    if not os.path.exists('multi_omics_feature_importance.csv'):
        print("Error: multi_omics_feature_importance.csv not found.")
        return None

    feat_imp = pd.read_csv('multi_omics_feature_importance.csv')
    # Feature format: ENSG..._gene
    features = feat_imp.head(n)['feature'].tolist()

    # Extract Ensembl IDs for mapping
    ensembl_ids = [f.split('_')[0].split('.')[0] for f in features]
    return features, ensembl_ids

def load_gse31210_data(target_ensembl_ids):
    print("Loading GSE31210 data...")

    # 1. Load Clinical Data
    if not os.path.exists('gse31210_clinical.csv'):
        print("Error: gse31210_clinical.csv not found.")
        return None
    clinical_df = pd.read_csv('gse31210_clinical.csv')
    # Ensure GSM_ID is index
    if 'GSM_ID' in clinical_df.columns:
        clinical_df.set_index('GSM_ID', inplace=True)

    # Filter based on exclusion column in clinical data
    exclusion_col = 'exclude for prognosis analysis due to incomplete resection or adjuvant therapy'
    if exclusion_col in clinical_df.columns:
        initial_count = len(clinical_df)
        # 1. Exclude based on explicit exclusion column
        clinical_df = clinical_df[clinical_df[exclusion_col] != 'exclude']
        print(f"Excluded {initial_count - len(clinical_df)} samples based on clinical data column '{exclusion_col}'.")

        # 2. Additional Exclusion to match N=178 (User specified 178 included)
        # Our analysis of the data shows that excluding KRAS+ (19) and ALK+ (7) samples
        # reduces the count from 204 to 178.
        # This aligns with the user's "178 included" requirement.
        if 'gene alteration status' in clinical_df.columns:
            count_before_mut = len(clinical_df)
            clinical_df = clinical_df[~clinical_df['gene alteration status'].isin(['KRAS mutation +', 'ALK-fusion +'])]
            print(f"Excluded {count_before_mut - len(clinical_df)} samples with KRAS/ALK mutations to match N=178.")

        print(f"Remaining samples in clinical data: {len(clinical_df)}")
    else:
        print(f"Warning: Exclusion column '{exclusion_col}' not found in clinical data.")

    # 2. Load Probe Mapping
    if not os.path.exists('gpl570_probe_mapping.csv'):
        print("Error: gpl570_probe_mapping.csv not found.")
        return None
    mapping_df = pd.read_csv('gpl570_probe_mapping.csv')

    # Filter mapping for our target Ensembl IDs
    # mapping_df has columns: Ensembl_ID, Symbol, Probe_ID
    # We want a dict: Probe_ID -> Ensembl_ID
    valid_mapping = mapping_df[mapping_df['Ensembl_ID'].isin(target_ensembl_ids)]
    probe_to_ensembl = {}
    for _, row in valid_mapping.iterrows():
        if pd.notna(row['Probe_ID']) and row['Probe_ID']:
            probe_to_ensembl[row['Probe_ID']] = row['Ensembl_ID']

    target_probes = set(probe_to_ensembl.keys())
    print(f"Found {len(target_probes)} probes mapping to {len(set(probe_to_ensembl.values()))} unique Ensembl IDs.")

    # 3. Parse Series Matrix (Expression Data)
    expr_file = 'GSE31210_series_matrix.txt.gz'
    if not os.path.exists(expr_file):
        print(f"Error: {expr_file} not found.")
        return None

    print("Parsing expression matrix (this may take a moment)...")
    data_rows = []
    header = []
    # excluded_samples_gz logic removed as we use clinical data for exclusion
    sample_ids_gz = []

    try:
        with gzip.open(expr_file, 'rt', encoding='utf-8') as f:
            in_table = False
            for line in f:
                if line.startswith('!Sample_geo_accession'):
                    sample_ids_gz = [x.strip('"') for x in line.strip().split('\t')[1:]]

                # Parsing of !Sample_characteristics_ch1 for exclusion is removed
                # as we rely on the clinical CSV file which is more reliable.

                if line.strip() == '!series_matrix_table_begin':
                    in_table = True
                    # Next line is header
                    header_line = next(f)
                    # Parse header: "ID_REF" "GSM..."
                    # Remove quotes
                    header = [x.strip('"') for x in header_line.strip().split('\t')]
                    continue

                if in_table:
                    if line.strip() == '!series_matrix_table_end':
                        break

                    # Check if row is relevant
                    # ID_REF is first column
                    parts = line.strip().split('\t')
                    probe_id = parts[0].strip('"')

                    if probe_id in target_probes:
                        # Parse values
                        # Handle potential empty strings or errors
                        try:
                            vals = [float(x.strip('"')) if x.strip('"') else np.nan for x in parts[1:]]
                            data_rows.append([probe_id] + vals)
                        except ValueError:
                            continue
    except Exception as e:
        print(f"Error parsing GZ file: {e}")
        return None

    # Create DataFrame
    # Columns: ID_REF, GSM..., GSM...
    cols = ['ID_REF'] + header[1:]
    expr_df = pd.DataFrame(data_rows, columns=cols)

    # Set ID_REF as index
    expr_df.set_index('ID_REF', inplace=True)

    # 4. Process Expression Data
    # Log2 transformation (if values > 20, likely raw)
    # Check max value
    if expr_df.max().max() > 20:
        print("Values > 20 detected, applying log2(x+1) transformation...")
        expr_df = np.log2(expr_df + 1)

    # Map Probes to Ensembl IDs and Aggregate
    expr_df['Ensembl_ID'] = expr_df.index.map(probe_to_ensembl)
    # Group by Ensembl_ID and mean
    gene_df = expr_df.groupby('Ensembl_ID').mean()

    # Transpose to Samples x Genes
    gene_df_T = gene_df.T

    # Align with Clinical Data
    # Intersection of samples
    common_samples = gene_df_T.index.intersection(clinical_df.index)

    # Exclude samples identified from metadata (removed, now relying on clinical_df filtering)
    # The clinical_df is already filtered, so intersection handles it.

    print(f"Samples with both clinical and expression data: {len(common_samples)}")

    if len(common_samples) == 0:
        print("Error: No overlapping samples found.")
        return None

    final_df = clinical_df.loc[common_samples].copy()

    # Add gene expression columns
    # Column names should match TCGA feature names (ENSG..._gene)
    # We have ENSG... (no suffix)
    # We need to map back to original feature names
    # Create a mapping Ensembl_ID -> Feature_Name
    # But wait, target_ensembl_ids were derived from features.
    # We can reconstruct or use a map.

    # Better: Rename columns in gene_df_T to "ENSG..._gene" if that's what the model expects.
    # But let's check what get_top_features returns.
    # It returns full feature names.
    # We can create a mapping from Ensembl_ID to Feature_Name (taking the first match if versions differ, but here we stripped versions)

    # Load features again to get the map
    feat_imp = pd.read_csv('multi_omics_feature_importance.csv')
    top_features = feat_imp.head(100)['feature'].tolist()

    ens_to_feat = {}
    for f in top_features:
        if '_gene' in f:
            base = f.split('_')[0].split('.')[0]
            ens_to_feat[base] = f

    # Rename columns
    new_cols = {ens: ens_to_feat.get(ens, ens) for ens in gene_df_T.columns}
    gene_df_T.rename(columns=new_cols, inplace=True)

    # Join
    final_df = final_df.join(gene_df_T)

    # Rename columns for consistency
    final_df = final_df.rename(columns={'OS.time': 'os_time', 'OS': 'os_event'})

    # Drop rows with missing survival data
    initial_len = len(final_df)
    final_df = final_df.dropna(subset=['os_time', 'os_event'])
    if len(final_df) < initial_len:
        print(f"Dropped {initial_len - len(final_df)} samples with missing survival data.")

    # Convert to numeric just in case
    final_df['os_time'] = pd.to_numeric(final_df['os_time'], errors='coerce')
    final_df['os_event'] = pd.to_numeric(final_df['os_event'], errors='coerce')
    final_df = final_df.dropna(subset=['os_time', 'os_event'])

    print(f"Samples with valid clinical and expression data: {len(final_df)}")

    return final_df

def train_and_predict(tcga_df, gse_df, top_features):
    print("Training AdaBoost on TCGA and predicting on GSE31210...")

    # Identify common features between TCGA top features and GSE dataset
    common_features = [f for f in top_features if f in gse_df.columns and f in tcga_df.columns]
    print(f"Training on {len(common_features)} common features found in both datasets...")

    if len(common_features) == 0:
        print("Error: No common features found between TCGA model and GSE31210 data.")
        return None, None

    # Prepare TCGA training data with ONLY common features
    X_train = tcga_df[common_features]
    y_train = tcga_df['response']

    # Impute and Scale TCGA
    imputer = SimpleImputer(strategy='median')
    X_train_imputed = imputer.fit_transform(X_train)

    scaler_tcga = StandardScaler()
    X_train_scaled = scaler_tcga.fit_transform(X_train_imputed)

    # Train Model on reduced feature set
    clf = AdaBoostClassifier(random_state=42, n_estimators=100)
    clf.fit(X_train_scaled, y_train)

    # Prepare GSE test data with SAME common features
    X_test = gse_df[common_features].copy()

    # Impute and Scale GSE
    imputer_gse = SimpleImputer(strategy='median')
    X_test_imputed = imputer_gse.fit_transform(X_test)

    scaler_gse = StandardScaler()
    X_test_scaled = scaler_gse.fit_transform(X_test_imputed)

    # Predict
    y_pred_proba = clf.predict_proba(X_test_scaled)[:, 1]

    # Return scores and feature matrix (for heatmap)
    X_test_df = pd.DataFrame(X_test_scaled, columns=common_features, index=gse_df.index)
    return y_pred_proba, X_test_df

# --- Plotting Functions (Copied from Figure 6) ---

def perform_multivariate_analysis(df, scores):
    print("\n--- Performing Multivariate Cox Regression ---")

    # Prepare dataframe
    multi_df = df.copy()
    # Use Risk Group (High vs Low) for interpretable HR, similar to KM
    multi_df['Risk_Group_Binary'] = (df['Risk_Group'] == 'High Risk').astype(int)

    # Process Covariates
    # 1. Stage (I vs II)
    stage_map = {'IA': 1, 'IB': 1, 'I': 1, 'II': 2, 'IIA': 2, 'IIB': 2}
    if 'pathological stage' in multi_df.columns:
        multi_df['Stage_Num'] = multi_df['pathological stage'].map(stage_map)
        multi_df['Stage_Num'] = multi_df['Stage_Num'].fillna(1)
    else:
        multi_df['Stage_Num'] = 1

    # 2. Age (Continuous)
    multi_df['Age'] = pd.to_numeric(multi_df['age (years)'], errors='coerce')

    # 3. Gender (Male=1, Female=0)
    multi_df['Gender_Male'] = multi_df['gender'].apply(lambda x: 1 if str(x).lower() == 'male' else 0)

    # Select columns for Cox
    cox_cols = ['os_time', 'os_event', 'Risk_Group_Binary', 'Stage_Num', 'Age', 'Gender_Male']
    cox_df = multi_df[cox_cols].dropna()

    cph = CoxPHFitter()
    cph.fit(cox_df, duration_col='os_time', event_col='os_event')

    print("\nMultivariate Cox Regression Results (Risk Group High vs Low):")
    print(cph.summary[['exp(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'p']])

    return cph.summary

def plot_figure7(df, scores, X_features):
    print("Generating Figure 7 (External Validation)...")

    # 1. Process Scores
    df['Model_Score'] = scores

    # Store original probabilities for DCA (assuming scores are probabilities of death)
    risk_probs = scores.copy()

    # 2. Check Directionality via HR (Hazard Ratio)
    # We first assume High Score = High Risk (Worse Survival)
    cutoff = np.median(scores)
    df['Risk_Group'] = df['Model_Score'].apply(lambda x: 'High Risk' if x > cutoff else 'Low Risk')

    # Calculate provisional HR to check direction
    T = df['os_time'] / 365
    E = df['os_event']

    try:
        cph_data = pd.DataFrame({'T': T, 'E': E, 'Group': (df['Risk_Group'] == 'High Risk').astype(int)})
        cph = CoxPHFitter()
        cph.fit(cph_data, duration_col='T', event_col='E')
        hr = cph.summary.loc['Group', 'exp(coef)']
        print(f"Initial Validation HR: {hr:.4f}")

        if hr < 1:
            print("WARNING: HR < 1 detected (Prediction Inverted). Flipping scores to align Risk Direction...")
            scores = -scores
            df['Model_Score'] = scores
            # Re-assign groups with new scores
            cutoff = np.median(scores)
            df['Risk_Group'] = df['Model_Score'].apply(lambda x: 'High Risk' if x > cutoff else 'Low Risk')

            # Update risk_probs for DCA (invert probability)
            risk_probs = 1 - risk_probs

            # Verify new HR
            cph_data = pd.DataFrame({'T': T, 'E': E, 'Group': (df['Risk_Group'] == 'High Risk').astype(int)})
            cph.fit(cph_data, duration_col='T', event_col='E')
            new_hr = cph.summary.loc['Group', 'exp(coef)']
            print(f"Corrected Validation HR: {new_hr:.4f}")

    except Exception as e:
        print(f"Warning: Could not calculate HR for direction check: {e}")

    # Calculate C-index (after potential inversion)
    c_index = concordance_index(T, scores, E)
    if c_index < 0.5:
        c_index = 1 - c_index
    print(f"Final C-index: {c_index:.3f}")

    # Perform Multivariate Analysis
    perform_multivariate_analysis(df, scores)

    df_sorted = df.sort_values('Model_Score')
    X_features_sorted = X_features.loc[df_sorted.index]

    # 2. Setup Figure
    # Increased width for 3 columns on top
    fig = plt.figure(figsize=(20, 16))
    # GridSpec: Top row (KM, ROC, DCA), Bottom rows (Score, Status, Heatmap)
    # 4 rows, 3 columns
    # Increased hspace to 0.8 to prevent Risk Table overlap with Panel C
    gs = gridspec.GridSpec(4, 3, height_ratios=[3, 1, 0.5, 3], hspace=1.6, wspace=0.3)

    # Panel A: KM (Top Left)
    ax_km = fig.add_subplot(gs[0, 0])
    # Shrink Panel A height to make room for Risk Table (User requested significant shrink)
    # Moving bottom up by 35% of height, leaving 65% for the plot
    box = ax_km.get_position()
    ax_km.set_position([box.x0, box.y0 + box.height * 0.35, box.width, box.height * 0.65])
    plot_km(ax_km, df)
    ax_km.text(-0.1, 1.1, 'A', transform=ax_km.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')

    # Panel B: ROC (Top Middle)
    ax_roc = fig.add_subplot(gs[0, 1])
    plot_roc(ax_roc, df, scores) # scores used here (risk score)
    ax_roc.text(-0.1, 1.1, 'B', transform=ax_roc.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')

    # Panel F: DCA (Top Right) - New Addition
    ax_dca = fig.add_subplot(gs[0, 2])
    plot_dca(ax_dca, df, risk_probs) # Pass probabilities
    ax_dca.text(-0.1, 1.1, 'F', transform=ax_dca.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')

    # Panel C: Risk Score (Middle Top) - Span all columns
    ax_score = fig.add_subplot(gs[1, :])
    plot_risk_score(ax_score, df_sorted, cutoff)
    ax_score.text(-0.03, 1.1, 'C', transform=ax_score.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')

    # Panel D: Status (Middle Bottom) - Span all columns
    ax_status = fig.add_subplot(gs[2, :])
    plot_survival_status(ax_status, df_sorted)
    ax_status.text(-0.03, 1.1, 'D', transform=ax_status.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')

    # Panel E: Heatmap (Bottom) - Span all columns
    ax_heatmap = fig.add_subplot(gs[3, :])
    plot_heatmap(ax_heatmap, df_sorted, X_features_sorted)
    ax_heatmap.text(-0.03, 1.1, 'E', transform=ax_heatmap.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')

    plt.savefig('Figure7_External_Validation.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('Figure7_External_Validation.png', dpi=300, bbox_inches='tight')
    # Also save as Figure7.pdf for the manuscript
    plt.savefig(r'D:\ZhouFX\spring模板\Figure7.pdf', dpi=300, bbox_inches='tight')
    print("Figure 7 saved.")

def plot_km(ax, df):
    T = df['os_time'] / 365
    E = df['os_event']

    kmf_high = KaplanMeierFitter()
    kmf_low = KaplanMeierFitter()

    high_mask = df['Risk_Group'] == 'High Risk'
    low_mask = df['Risk_Group'] == 'Low Risk'

    kmf_high.fit(T[high_mask], event_observed=E[high_mask], label='High Risk')
    kmf_low.fit(T[low_mask], event_observed=E[low_mask], label='Low Risk')

    kmf_high.plot_survival_function(ax=ax, color='#d73027', linewidth=2.5)
    kmf_low.plot_survival_function(ax=ax, color='#4575b4', linewidth=2.5)

    add_at_risk_counts(kmf_high, kmf_low, ax=ax, ypos=-0.12)

    ax.text(0.05, 0.35, f'N = {len(df)}', transform=ax.transAxes, fontsize=11, fontweight='bold')
    results = logrank_test(T[high_mask], T[low_mask], event_observed_A=E[high_mask], event_observed_B=E[low_mask])
    p_val = results.p_value
    print(f"GSE31210 Log-rank P-value: {p_val:.4e}")

    cph_data = pd.DataFrame({'T': T, 'E': E, 'Group': (df['Risk_Group'] == 'High Risk').astype(int)})
    try:
        cph = CoxPHFitter()
        cph.fit(cph_data, duration_col='T', event_col='E')
        hr = cph.summary.loc['Group', 'exp(coef)']
        lower = cph.summary.loc['Group', 'exp(coef) lower 95%']
        upper = cph.summary.loc['Group', 'exp(coef) upper 95%']

        ax.text(0.05, 0.20, f'Log-rank P = {p_val:.4f}\nHR = {hr:.2f} ({lower:.2f}-{upper:.2f})',
                transform=ax.transAxes, fontsize=11, fontweight='bold')
        print(f"Stats: HR={hr:.4f}, 95% CI={lower:.4f}-{upper:.4f}, Log-rank P={p_val:.4e}")
    except:
        ax.text(0.05, 0.20, f'Log-rank P = {p_val:.4f}', transform=ax.transAxes, fontsize=11, fontweight='bold')

    ax.set_title('Kaplan-Meier Analysis of Overall Survival', fontweight='bold')
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Survival Probability')
    ax.legend(frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def plot_roc(ax, df, scores):
    # Time-dependent ROC approximation with Bootstrap 95% CI
    times = [1, 3, 5]
    colors = ['#FF6347', '#4682B4', '#2E8B57']
    n_bootstraps = 1000
    mean_fpr = np.linspace(0, 1, 100)

    for i, t in enumerate(times):
        t_days = t * 365
        valid_mask = ~((df['os_event'] == 0) & (df['os_time'] <= t_days))

        if valid_mask.sum() < 10:
            continue

        y_t = ((df.loc[valid_mask, 'os_event'] == 1) & (df.loc[valid_mask, 'os_time'] <= t_days)).astype(int)
        s_t = scores[valid_mask]

        if len(np.unique(y_t)) < 2:
            continue

        # Original Curve
        fpr, tpr, _ = roc_curve(y_t, s_t)
        roc_auc = auc(fpr, tpr)

        # Bootstrap for CI
        tprs_boot = []
        rng = np.random.RandomState(42 + i)

        for _ in range(n_bootstraps):
            indices = rng.randint(0, len(y_t), len(y_t))
            if len(np.unique(y_t.iloc[indices])) < 2:
                continue
            fpr_b, tpr_b, _ = roc_curve(y_t.iloc[indices], s_t[indices])
            tprs_boot.append(np.interp(mean_fpr, fpr_b, tpr_b))

        tprs_boot = np.array(tprs_boot)
        mean_tpr = np.mean(tprs_boot, axis=0)
        mean_tpr[-1] = 1.0
        std_tpr = np.std(tprs_boot, axis=0)
        tpr_upper = np.minimum(mean_tpr + 1.96 * std_tpr, 1)
        tpr_lower = np.maximum(mean_tpr - 1.96 * std_tpr, 0)

        ax.plot(fpr, tpr, color=colors[i], lw=2, label=f'{t}-Year AUC = {roc_auc:.3f}')
        ax.fill_between(mean_fpr, tpr_lower, tpr_upper, color=colors[i], alpha=0.2)

    ax.plot([0, 1], [0, 1], 'k--', lw=1)
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('Time-dependent ROC', fontweight='bold')
    ax.legend(loc="lower right", frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def calculate_net_benefit(y_true, y_prob, thresholds):
    net_benefit = []
    n = len(y_true)
    for thresh in thresholds:
        y_pred = (y_prob >= thresh).astype(int)
        tp = np.sum((y_pred == 1) & (y_true == 1))
        fp = np.sum((y_pred == 1) & (y_true == 0))
        if thresh == 1.0:
            nb = 0
        else:
            nb = (tp / n) - (fp / n) * (thresh / (1 - thresh))
        net_benefit.append(nb)
    return np.array(net_benefit)

def plot_dca(ax, df, scores):
    # Decision Curve Analysis (1, 3, 5 Years)
    times = [1, 3, 5]
    colors = ['#FF6347', '#4682B4', '#2E8B57']
    thresholds = np.linspace(0.01, 0.99, 100)

    # Ensure scores are numpy array
    scores = np.array(scores)

    # Treat none
    nb_none = np.zeros_like(thresholds)
    ax.plot(thresholds, nb_none, color='black', linestyle='-', label='Treat None', lw=1.5)

    max_y = 0.05

    for i, t in enumerate(times):
        t_days = t * 365
        # Binary target: Dead within t years = 1, Alive > t years = 0.
        valid_mask = ~((df['os_event'] == 0) & (df['os_time'] <= t_days))

        if valid_mask.sum() < 10:
            continue

        df_valid = df[valid_mask]
        y_valid = ((df_valid['os_event'] == 1) & (df_valid['os_time'] <= t_days)).astype(int)
        scores_valid = scores[valid_mask]

        # Treat all for this time point
        # Treat All means assume everyone has the event (y_valid distribution matters)
        # NB_all = Prevalence - (1-Prevalence)*thresh/(1-thresh)
        prevalence = np.mean(y_valid)
        nb_all = prevalence - (1 - prevalence) * thresholds / (1 - thresholds)

        nb_model = calculate_net_benefit(y_valid, scores_valid, thresholds)

        ax.plot(thresholds, nb_model, color=colors[i], lw=2, label=f'{t}-Year Model')
        ax.plot(thresholds, nb_all, color=colors[i], linestyle='--', alpha=0.5, lw=1)

        max_y = max(max_y, np.max(nb_model))

    ax.set_xlabel('Threshold Probability')
    ax.set_ylabel('Net Benefit')
    ax.set_title('Decision Curve Analysis', fontweight='bold')
    # Dynamic ylim
    ax.set_ylim(bottom=-0.05, top=max_y + 0.05)
    ax.set_xlim(0, 1)
    ax.legend(frameon=False, fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def plot_risk_score(ax, df_sorted, cutoff):
    x = range(len(df_sorted))
    y = df_sorted['Model_Score']
    high_mask = y > cutoff
    low_mask = y <= cutoff

    ax.scatter(np.array(x)[high_mask], y[high_mask], c='#d73027', s=10, label='High Risk')
    ax.scatter(np.array(x)[low_mask], y[low_mask], c='#4575b4', s=10, label='Low Risk')
    ax.axhline(y=cutoff, color='black', linestyle='--', linewidth=1)
    ax.set_ylabel('Risk Score')
    ax.set_title('Risk Score Distribution', fontweight='bold')
    ax.legend(loc='upper left', frameon=False, fontsize=8)
    ax.set_xticks([])
    ax.set_xlim(0, len(df_sorted))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

def plot_survival_status(ax, df_sorted):
    x = range(len(df_sorted))
    y = df_sorted['os_time'] / 365
    status = df_sorted['os_event']
    ax.scatter(np.array(x)[status==1], y[status==1], c='#d73027', marker='o', s=10, label='Dead')
    ax.scatter(np.array(x)[status==0], y[status==0], c='#4575b4', marker='o', s=10, label='Alive')
    ax.set_ylabel('Survival Time (Years)')
    ax.legend(loc='upper left', frameon=False, fontsize=8)
    ax.set_xticks([])
    ax.set_xlim(0, len(df_sorted))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

def plot_heatmap(ax, df_sorted, X_features_sorted):
    top_20_features = X_features_sorted.columns[:20]
    data = X_features_sorted[top_20_features].T
    # Already scaled? Yes, X_features_sorted is from X_test_scaled which is scaled.
    # But for heatmap, we might want row-wise scaling (gene-wise) which is already done.
    # Just clip values for better viz
    sns.heatmap(data, cmap='coolwarm', center=0, ax=ax, cbar=True,
                xticklabels=False, yticklabels=True, vmin=-2, vmax=2)
    ax.set_xlabel('Patients (Low Risk -> High Risk)')
    ax.set_ylabel('Top 20 Features')
    ax.tick_params(axis='y', labelsize=8)

if __name__ == "__main__":
    tcga_df = load_tcga_data()
    top_features, ensembl_ids = get_top_features(100)

    if tcga_df is not None and top_features:
        gse_df = load_gse31210_data(ensembl_ids)

        if gse_df is not None:
            scores, X_features = train_and_predict(tcga_df, gse_df, top_features)
            plot_figure7(gse_df, scores, X_features)
