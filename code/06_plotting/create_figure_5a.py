import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

# 读取Moran's I结果
df = pd.read_csv('d:/ZhouFX/Spatial_Analysis_Results/morans_i_results_manual.csv')
df.columns = ['Gene', 'Morans_I']

# 创建Figure 5A: Moran's I分布图
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# 左图：Moran's I值分布直方图
ax1.hist(df['Morans_I'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
ax1.axvline(x=0, color='red', linestyle='--', alpha=0.7, label='零自相关线')
ax1.set_xlabel('Moran\'s I 值', fontsize=12)
ax1.set_ylabel('基因数量', fontsize=12)
ax1.set_title('空间自相关系数分布', fontsize=14, fontweight='bold')
ax1.legend()
ax1.grid(True, alpha=0.3)

# 右图：Top 20基因的Moran's I值
top_20 = df.head(20)
colors = ['red' if gene in ['C1QA', 'C1QC', 'C1QB'] else 'steelblue' for gene in top_20['Gene']]
bars = ax2.barh(range(len(top_20)), top_20['Morans_I'], color=colors, alpha=0.7)
ax2.set_yticks(range(len(top_20)))
ax2.set_yticklabels(top_20['Gene'])
ax2.set_xlabel('Moran\'s I 值', fontsize=12)
ax2.set_title('Top 20 基因空间自相关系数', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)

# 添加数值标签
for i, (bar, value) in enumerate(zip(bars, top_20['Morans_I'])):
    ax2.text(value + 0.01, i, f'{value:.3f}', va='center', fontsize=9)

# 添加图例说明
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='red', alpha=0.7, label='C1Q家族基因'),
    Patch(facecolor='steelblue', alpha=0.7, label='其他基因')
]
ax2.legend(handles=legend_elements, loc='lower right')

plt.tight_layout()
plt.savefig('d:/ZhouFX/Figure_5A_Morans_I_Distribution.png', dpi=300, bbox_inches='tight')
plt.show()

# 打印关键统计信息
print("=== Moran's I 分析结果 ===")
print(f"总基因数: {len(df)}")
print(f"正空间自相关基因数 (I > 0): {sum(df['Morans_I'] > 0)}")
print(f"负空间自相关基因数 (I < 0): {sum(df['Morans_I'] < 0)}")
print(f"平均Moran's I值: {df['Morans_I'].mean():.4f}")
print(f"Moran's I值标准差: {df['Morans_I'].std():.4f}")

print("\n=== C1Q家族基因排名 ===")
c1q_genes = df[df['Gene'].str.contains('C1Q')]
print(c1q_genes.to_string(index=False))

print("\n=== Top 10 高空间自相关基因 ===")
print(df.head(10).to_string(index=False))