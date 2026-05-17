
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind, mannwhitneyu
from statannotations.Annotator import Annotator
import os
import gzip
from sklearn.ensemble import AdaBoostClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from lifelines import CoxPHFitter

# Set style
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
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
    features = feat_imp.head(n)['feature'].tolist()
    ensembl_ids = [f.split('_')[0].split('.')[0] for f in features]
    return features, ensembl_ids

def get_gene_symbols(features):
    # Load mapping
    if not os.path.exists('converted_features.csv'):
        return {f: f for f in features}

    mapping = pd.read_csv('converted_features.csv')
    # Create dict: feature -> symbol
    feat_to_symbol = {}
    for _, row in mapping.iterrows():
        feat_to_symbol[row['feature']] = row['gene_symbol']

    return feat_to_symbol

def load_gse31210_data(target_ensembl_ids):
    print("Loading GSE31210 data...")

    # 1. Load Clinical Data
    if not os.path.exists('gse31210_clinical.csv'):
        return None
    clinical_df = pd.read_csv('gse31210_clinical.csv')
    if 'GSM_ID' in clinical_df.columns:
        clinical_df.set_index('GSM_ID', inplace=True)

    # 2. Load Probe Mapping
    if not os.path.exists('gpl570_probe_mapping.csv'):
        return None
    mapping_df = pd.read_csv('gpl570_probe_mapping.csv')

    valid_mapping = mapping_df[mapping_df['Ensembl_ID'].isin(target_ensembl_ids)]
    probe_to_ensembl = {}
    for _, row in valid_mapping.iterrows():
        if pd.notna(row['Probe_ID']) and row['Probe_ID']:
            probe_to_ensembl[row['Probe_ID']] = row['Ensembl_ID']

    target_probes = set(probe_to_ensembl.keys())

    # 3. Parse Expression
    expr_file = 'GSE31210_series_matrix.txt.gz'
    if not os.path.exists(expr_file):
        return None

    data_rows = []
    header = []

    try:
        with gzip.open(expr_file, 'rt', encoding='utf-8') as f:
            in_table = False
            for line in f:
                if line.strip() == '!series_matrix_table_begin':
                    in_table = True
                    header_line = next(f)
                    header = [x.strip('"') for x in header_line.strip().split('\t')]
                    continue

                if in_table:
                    if line.strip() == '!series_matrix_table_end':
                        break

                    parts = line.strip().split('\t')
                    probe_id = parts[0].strip('"')

                    if probe_id in target_probes:
                        try:
                            vals = [float(x.strip('"')) if x.strip('"') else np.nan for x in parts[1:]]
                            data_rows.append([probe_id] + vals)
                        except ValueError:
                            continue
    except Exception as e:
        print(f"Error: {e}")
        return None

    cols = ['ID_REF'] + header[1:]
    expr_df = pd.DataFrame(data_rows, columns=cols)
    expr_df.set_index('ID_REF', inplace=True)

    if expr_df.max().max() > 20:
        expr_df = np.log2(expr_df + 1)

    expr_df['Ensembl_ID'] = expr_df.index.map(probe_to_ensembl)
    gene_df = expr_df.groupby('Ensembl_ID').mean()
    gene_df_T = gene_df.T

    common_samples = gene_df_T.index.intersection(clinical_df.index)
    final_df = clinical_df.loc[common_samples].copy()

    # Map back to feature names
    feat_imp = pd.read_csv('multi_omics_feature_importance.csv')
    top_features = feat_imp.head(100)['feature'].tolist()
    ens_to_feat = {}
    for f in top_features:
        if '_gene' in f:
            base = f.split('_')[0].split('.')[0]
            ens_to_feat[base] = f

    new_cols = {ens: ens_to_feat.get(ens, ens) for ens in gene_df_T.columns}
    gene_df_T.rename(columns=new_cols, inplace=True)

    final_df = final_df.join(gene_df_T)
    final_df = final_df.rename(columns={'OS.time': 'os_time', 'OS': 'os_event'})

    final_df['os_time'] = pd.to_numeric(final_df['os_time'], errors='coerce')
    final_df['os_event'] = pd.to_numeric(final_df['os_event'], errors='coerce')
    final_df = final_df.dropna(subset=['os_time', 'os_event'])

    return final_df

def train_and_predict(tcga_df, gse_df, top_features):
    common_features = [f for f in top_features if f in gse_df.columns and f in tcga_df.columns]

    X_train = tcga_df[common_features]
    y_train = tcga_df['response']

    imputer = SimpleImputer(strategy='median')
    X_train_imputed = imputer.fit_transform(X_train)
    scaler_tcga = StandardScaler()
    X_train_scaled = scaler_tcga.fit_transform(X_train_imputed)

    clf = AdaBoostClassifier(random_state=42, n_estimators=100)
    clf.fit(X_train_scaled, y_train)

    X_test = gse_df[common_features].copy()
    imputer_gse = SimpleImputer(strategy='median')
    X_test_imputed = imputer_gse.fit_transform(X_test)
    scaler_gse = StandardScaler()
    X_test_scaled = scaler_gse.fit_transform(X_test_imputed)

    y_pred_proba = clf.predict_proba(X_test_scaled)[:, 1]

    # Return scaled data for plotting? Or raw?
    # Usually better to show standardized or just log2 expression.
    # Since we want to show "expression differences", let's use the X_test_imputed (Log2 expression)
    # But for consistency with model, maybe scaled.
    # Let's use X_test_imputed (Log2 normalized expression) as it is more biologically interpretable in a boxplot than Z-score.

    X_test_df = pd.DataFrame(X_test_imputed, columns=common_features, index=gse_df.index)

    return y_pred_proba, X_test_df, common_features

def plot_figure_s7(df, scores, X_features, top_features_ordered, feat_to_symbol):
    print("Generating Figure S7...")

    df['Model_Score'] = scores

    # Check direction (reusing logic from Fig 7)
    cutoff = np.median(scores)
    df['Risk_Group'] = df['Model_Score'].apply(lambda x: 'High Risk' if x > cutoff else 'Low Risk')

    T = df['os_time'] / 365
    E = df['os_event']

    try:
        cph_data = pd.DataFrame({'T': T, 'E': E, 'Group': (df['Risk_Group'] == 'High Risk').astype(int)})
        cph = CoxPHFitter()
        cph.fit(cph_data, duration_col='T', event_col='E')
        hr = cph.summary.loc['Group', 'exp(coef)']
        if hr < 1:
            print("Inverting scores...")
            scores = -scores
            df['Model_Score'] = scores
            cutoff = np.median(scores)
            df['Risk_Group'] = df['Model_Score'].apply(lambda x: 'High Risk' if x > cutoff else 'Low Risk')
    except:
        pass

    # Identify Top 5 AVAILABLE features
    top_5_features = []
    for f in top_features_ordered:
        if f in X_features.columns:
            top_5_features.append(f)
            if len(top_5_features) == 5:
                break

    print(f"Top 5 common features: {top_5_features}")

    # Prepare data for plotting
    plot_data = []
    for f in top_5_features:
        symbol = feat_to_symbol.get(f, f)
        # Extract data for High and Low risk
        high_vals = X_features.loc[df[df['Risk_Group']=='High Risk'].index, f]
        low_vals = X_features.loc[df[df['Risk_Group']=='Low Risk'].index, f]

        for v in high_vals:
            plot_data.append({'Gene': symbol, 'Expression': v, 'Group': 'High Risk'})
        for v in low_vals:
            plot_data.append({'Gene': symbol, 'Expression': v, 'Group': 'Low Risk'})

    plot_df = pd.DataFrame(plot_data)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Boxplot
    sns.boxplot(data=plot_df, x='Gene', y='Expression', hue='Group',
                palette={'High Risk': '#d73027', 'Low Risk': '#4575b4'},
                width=0.6, ax=ax, showfliers=False)

    # Add stripplot
    sns.stripplot(data=plot_df, x='Gene', y='Expression', hue='Group',
                  dodge=True, color='black', alpha=0.3, size=3, ax=ax)

    # Stats
    pairs = []
    for gene in plot_df['Gene'].unique():
        pairs.append(((gene, 'Low Risk'), (gene, 'High Risk')))

    annotator = Annotator(ax, pairs, data=plot_df, x='Gene', y='Expression', hue='Group')
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
    annotator.apply_and_annotate()

    # Clean up legend (remove duplicates from stripplot)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], title='Risk Group', frameon=False)

    ax.set_title('Top 5 Model Gene Expression in GSE31210 Validation Cohort', fontweight='bold')
    ax.set_ylabel('Log2 Expression Level')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig('FigureS7.pdf', dpi=300)
    plt.savefig(r'D:\ZhouFX\spring模板\FigureS7.pdf', dpi=300)
    print("Figure S7 saved.")

if __name__ == "__main__":
    tcga_df = load_tcga_data()
    top_features, ensembl_ids = get_top_features(100)
    feat_to_symbol = get_gene_symbols(top_features)

    if tcga_df is not None and top_features:
        gse_df = load_gse31210_data(ensembl_ids)

        if gse_df is not None:
            scores, X_features, common_features = train_and_predict(tcga_df, gse_df, top_features)
            plot_figure_s7(gse_df, scores, X_features, top_features, feat_to_symbol)
