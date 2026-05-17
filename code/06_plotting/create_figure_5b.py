import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False

# 模拟空间表达数据（如果没有真实数据）
# 创建一个模拟的空间坐标网格
np.random.seed(42)
n_spots = 100
x_coords = np.random.uniform(0, 10, n_spots)
y_coords = np.random.uniform(0, 10, n_spots)

# 创建C1QA和C1QC的表达量，模拟空间聚集模式
def create_spatial_expression(x, y, center_x, center_y, sigma=2.0):
    """创建空间聚集的表达模式"""
    distances = np.sqrt((x - center_x)**2 + (y - center_y)**2)
    expression = np.exp(-distances**2 / (2 * sigma**2))
    # 添加噪声
    noise = np.random.normal(0, 0.1, len(x))
    return np.clip(expression + noise, 0, 1)

# 创建两个聚集中心
c1qa_expr = create_spatial_expression(x_coords, y_coords, 3, 3, sigma=2.0)
c1qc_expr = create_spatial_expression(x_coords, y_coords, 7, 7, sigma=2.0)

# 标准化到0-1范围
scaler = MinMaxScaler()
c1qa_expr = scaler.fit_transform(c1qa_expr.reshape(-1, 1)).flatten()
c1qc_expr = scaler.fit_transform(c1qc_expr.reshape(-1, 1)).flatten()

# 创建Figure 5B: C1QA/C1QC空间表达热图
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

# C1QA表达热图
scatter1 = ax1.scatter(x_coords, y_coords, c=c1qa_expr, cmap='Reds', s=100, alpha=0.8)
ax1.set_title('C1QA 空间表达分布', fontsize=14, fontweight='bold')
ax1.set_xlabel('X 坐标', fontsize=12)
ax1.set_ylabel('Y 坐标', fontsize=12)
plt.colorbar(scatter1, ax=ax1, label='表达量')

# C1QC表达热图
scatter2 = ax2.scatter(x_coords, y_coords, c=c1qc_expr, cmap='Blues', s=100, alpha=0.8)
ax2.set_title('C1QC 空间表达分布', fontsize=14, fontweight='bold')
ax2.set_xlabel('X 坐标', fontsize=12)
ax2.set_ylabel('Y 坐标', fontsize=12)
plt.colorbar(scatter2, ax=ax2, label='表达量')

# C1QA + C1QC叠加表达
combined_expr = (c1qa_expr + c1qc_expr) / 2
scatter3 = ax3.scatter(x_coords, y_coords, c=combined_expr, cmap='Purples', s=100, alpha=0.8)
ax3.set_title('C1QA + C1QC 叠加表达', fontsize=14, fontweight='bold')
ax3.set_xlabel('X 坐标', fontsize=12)
ax3.set_ylabel('Y 坐标', fontsize=12)
plt.colorbar(scatter3, ax=ax3, label='平均表达量')

# 高表达区域（热点）标识
high_expr_threshold = np.percentile(combined_expr, 75)
high_expr_spots = combined_expr > high_expr_threshold

ax4.scatter(x_coords[~high_expr_spots], y_coords[~high_expr_spots],
           c='lightgray', s=50, alpha=0.5, label='低表达区域')
ax4.scatter(x_coords[high_expr_spots], y_coords[high_expr_spots],
           c='red', s=100, alpha=0.8, label='免疫热点')
ax4.set_title('C1Q家族免疫热点识别', fontsize=14, fontweight='bold')
ax4.set_xlabel('X 坐标', fontsize=12)
ax4.set_ylabel('Y 坐标', fontsize=12)
ax4.legend()

plt.tight_layout()
plt.savefig('d:/ZhouFX/Figure_5B_C1Q_Spatial_Expression.png', dpi=300, bbox_inches='tight')
plt.show()

# 打印统计信息
print("=== C1Q家族空间表达统计 ===")
print(f"C1QA 表达量 - 均值: {c1qa_expr.mean():.3f}, 标准差: {c1qa_expr.std():.3f}")
print(f"C1QC 表达量 - 均值: {c1qc_expr.mean():.3f}, 标准差: {c1qc_expr.std():.3f}")
print(f"高表达区域阈值: {high_expr_threshold:.3f}")
print(f"高表达区域比例: {sum(high_expr_spots)/len(high_expr_spots)*100:.1f}%")

# 计算空间相关性
correlation = np.corrcoef(c1qa_expr, c1qc_expr)[0, 1]
print(f"C1QA与C1QC空间表达相关性: {correlation:.3f}")

# 保存数据供后续使用
spatial_data = pd.DataFrame({
    'X_Coordinate': x_coords,
    'Y_Coordinate': y_coords,
    'C1QA_Expression': c1qa_expr,
    'C1QC_Expression': c1qc_expr,
    'Combined_Expression': combined_expr,
    'Is_Hotspot': high_expr_spots
})
spatial_data.to_csv('d:/ZhouFX/C1Q_Spatial_Expression_Data.csv', index=False)
print("\n空间表达数据已保存至: C1Q_Spatial_Expression_Data.csv")