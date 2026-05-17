import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Set font to Arial
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

def generate_pdf():
    # Load data
    df = pd.read_csv('Table_S1_Full.csv')
    
    # Format numbers
    for col in df.columns:
        if 'Train_' in col or 'Test_' in col:
            df[col] = df[col].apply(lambda x: f"{x:.4f}")
    
    # Rename columns for brevity
    rename_map = {
        'Hyperparameters': 'Hyperparameters',
        'Train_Accuracy': 'Train Acc',
        'Train_Precision': 'Train Prec',
        'Train_Recall': 'Train Rec',
        'Train_F1': 'Train F1',
        'Train_AUC': 'Train AUC',
        'Val_Accuracy': 'Val Acc',
        'Val_Precision': 'Val Prec',
        'Val_Recall': 'Val Rec',
        'Val_F1': 'Val F1',
        'Val_AUC': 'Val AUC',
        'Test_Accuracy': 'Test Acc',
        'Test_Precision': 'Test Prec',
        'Test_Recall': 'Test Rec',
        'Test_F1': 'Test F1',
        'Test_AUC': 'Test AUC',
    }
    df_display = df.rename(columns=rename_map)
    
    # Wrap hyperparameters text
    def wrap_text(text, width=30):
        import textwrap
        return '\n'.join(textwrap.wrap(text, width))
    
    df_display['Hyperparameters'] = df_display['Hyperparameters'].apply(lambda x: wrap_text(str(x), 25))
    
    # Create figure
    # 17 columns. Model + Params + 15 metrics.
    # Params needs more space.
    # Ratios: Model(1), Params(3), Metrics(1)*15
    col_widths = [0.12, 0.25] + [0.06] * 15
    
    fig_width = 24
    fig_height = len(df) * 0.8 + 2
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.axis('off')
    
    # Create table
    table = ax.table(
        cellText=df_display.values,
        colLabels=df_display.columns,
        cellLoc='center',
        loc='center',
        colWidths=col_widths
    )
    
    # Style table
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    
    # Header styling
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight='bold')
            cell.set_facecolor('#f0f0f0')
    
    # Add title
    plt.title('Table S1: Detailed Hyperparameters and Performance Metrics of Machine Learning Models', 
              fontsize=16, fontweight='bold', pad=20)
    
    # Save
    plt.savefig('Table_S1.pdf', dpi=600, bbox_inches='tight')
    plt.close()
    print("Table S1 saved to Table_S1.pdf")

if __name__ == "__main__":
    generate_pdf()
