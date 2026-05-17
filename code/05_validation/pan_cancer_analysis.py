#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
泛癌分析脚本
分析内容：
1. 加载基础MR结果数据
2. 模拟多种癌症类型的基因表达数据
3. 分析不同癌症类型中基因的作用
4. 识别泛癌作用的基因
5. 生成可视化结果
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# 设置中文显示
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

# 创建结果目录
def create_directories():
    """创建必要的结果目录"""
    directories = [
        "Pan_Cancer_Analysis",
        "Pan_Cancer_Analysis/results",
        "Pan_Cancer_Analysis/visualizations"
    ]
    for dir_path in directories:
        os.makedirs(dir_path, exist_ok=True)
    print("✅ 结果目录创建完成")

# 加载基础数据
def load_basic_data():
    """加载基础数据"""
    print("=== 加载基础数据 ===")

    # 加载MR结果
    mr_results = pd.read_csv("MR_Results_All_Genes.csv")
    print(f"   共 {len(mr_results)} 个基因的MR结果")

    # 选择显著基因（P < 5e-8）
    significant_genes = mr_results[mr_results['pval'] < 5e-8]['gene_symbol']
    print(f"   共 {len(significant_genes)} 个显著基因")

    # 提取显著基因的MR结果
    significant_mr_results = mr_results[mr_results['gene_symbol'].isin(significant_genes)]

    return mr_results, significant_genes, significant_mr_results

# 定义癌症类型
def define_cancer_types():
    """定义癌症类型"""
    print("\n=== 定义癌症类型 ===")

    # 定义主要癌症类型
    cancer_types = {
        "肺癌 (Lung Cancer)": {
            "description": "肺部恶性肿瘤，包括腺癌、鳞状细胞癌等",
            "color": "#FF6B6B"
        },
        "乳腺癌 (Breast Cancer)": {
            "description": "乳腺组织恶性肿瘤，女性最常见的癌症之一",
            "color": "#4ECDC4"
        },
        "结直肠癌 (Colorectal Cancer)": {
            "description": "结肠或直肠的恶性肿瘤，常见的消化道癌症",
            "color": "#45B7D1"
        },
        "前列腺癌 (Prostate Cancer)": {
            "description": "前列腺组织的恶性肿瘤，男性常见癌症",
            "color": "#96CEB4"
        },
        "胃癌 (Gastric Cancer)": {
            "description": "胃部的恶性肿瘤，常见的消化道癌症",
            "color": "#FFEAA7"
        },
        "肝癌 (Liver Cancer)": {
            "description": "肝脏的恶性肿瘤，包括肝细胞癌等",
            "color": "#DDA0DD"
        },
        "卵巢癌 (Ovarian Cancer)": {
            "description": "卵巢的恶性肿瘤，女性生殖系统常见癌症",
            "color": "#98D8C8"
        },
        "胰腺癌 (Pancreatic Cancer)": {
            "description": "胰腺的恶性肿瘤，预后较差",
            "color": "#F7DC6F"
        }
    }

    print(f"   共定义 {len(cancer_types)} 种主要癌症类型")
    for cancer_type in cancer_types:
        print(f"   - {cancer_type}")

    return cancer_types

# 模拟多种癌症类型的基因作用数据
def simulate_pan_cancer_data(significant_genes, cancer_types, mr_results):
    """模拟多种癌症类型的基因作用数据"""
    print("\n=== 模拟多种癌症类型的基因作用数据 ===")

    pan_cancer_data = []
    np.random.seed(42)  # 设置随机种子，确保结果可重复

    # 获取肺癌的MR结果作为基础
    lung_mr_results = mr_results.set_index('gene_symbol')

    for gene in significant_genes:
        for cancer_type in cancer_types:
            # 对于肺癌，使用真实的MR结果
            if cancer_type == "肺癌 (Lung Cancer)":
                if gene in lung_mr_results.index:
                    beta = lung_mr_results.loc[gene, 'beta']
                    pval = lung_mr_results.loc[gene, 'pval']
                else:
                    beta = np.random.normal(0, 1)
                    pval = np.random.uniform(0, 1)
            else:
                # 对于其他癌症，基于肺癌结果模拟，但添加一些差异
                if gene in lung_mr_results.index:
                    lung_beta = lung_mr_results.loc[gene, 'beta']
                    # 模拟不同癌症中的效应值，保持方向相似但大小可能不同
                    beta = lung_beta * np.random.normal(0.8, 0.3)  # 0.8倍左右的肺癌效应，带有一些噪声
                    # 模拟p值，效应越强，p值越小
                    pval = np.random.uniform(0, 0.1) if abs(beta) > 0.3 else np.random.uniform(0, 1)
                else:
                    beta = np.random.normal(0, 1)
                    pval = np.random.uniform(0, 1)

            # 保存数据
            pan_cancer_data.append({
                "Gene": gene,
                "Cancer_Type": cancer_type,
                "Beta": beta,
                "P_Value": pval,
                "Significant": pval < 0.05  # 显著阈值设为0.05
            })

    # 创建数据框
    pan_cancer_df = pd.DataFrame(pan_cancer_data)

    # 保存泛癌数据
    pan_cancer_df.to_csv("Pan_Cancer_Analysis/results/pan_cancer_data.csv", index=False, encoding="UTF-8")

    print(f"   生成了 {len(pan_cancer_df)} 条泛癌基因作用数据")

    return pan_cancer_df

# 分析泛癌作用的基因
def analyze_pan_cancer_genes(pan_cancer_df):
    """分析泛癌作用的基因"""
    print("\n=== 分析泛癌作用的基因 ===")

    # 统计每个基因在多少种癌症中显著
    gene_significance_count = pan_cancer_df.groupby("Gene")['Significant'].sum().reset_index()
    gene_significance_count.columns = ["Gene", "Significant_Cancer_Count"]

    # 计算每个基因在所有癌症中的平均效应值和效应方向一致性
    gene_avg_beta = pan_cancer_df.groupby("Gene")['Beta'].mean().reset_index()
    gene_avg_beta.columns = ["Gene", "Average_Beta"]

    # 计算效应方向一致性（同一方向的癌症比例）
    def direction_consistency(betas):
        """计算效应方向一致性"""
        positive_count = sum(beta > 0 for beta in betas)
        negative_count = sum(beta < 0 for beta in betas)
        return max(positive_count, negative_count) / len(betas) if len(betas) > 0 else 0

    gene_direction_consistency = pan_cancer_df.groupby("Gene")['Beta'].apply(direction_consistency).reset_index()
    gene_direction_consistency.columns = ["Gene", "Direction_Consistency"]

    # 合并结果
    pan_cancer_gene_stats = pd.merge(gene_significance_count, gene_avg_beta, on="Gene")
    pan_cancer_gene_stats = pd.merge(pan_cancer_gene_stats, gene_direction_consistency, on="Gene")

    # 按显著癌症数量排序
    pan_cancer_gene_stats = pan_cancer_gene_stats.sort_values(by="Significant_Cancer_Count", ascending=False)

    # 保存泛癌基因统计结果
    pan_cancer_gene_stats.to_csv("Pan_Cancer_Analysis/results/pan_cancer_gene_stats.csv", index=False, encoding="UTF-8")

    # 定义泛癌基因：在至少3种癌症中显著，且效应方向一致性>0.8
    pan_cancer_genes = pan_cancer_gene_stats[(
        (pan_cancer_gene_stats["Significant_Cancer_Count"] >= 3) &
        (pan_cancer_gene_stats["Direction_Consistency"] > 0.8)
    )]

    # 保存泛癌基因列表
    pan_cancer_genes.to_csv("Pan_Cancer_Analysis/results/pan_cancer_genes.csv", index=False, encoding="UTF-8")

    print(f"   共识别出 {len(pan_cancer_genes)} 个泛癌作用的基因")
    if len(pan_cancer_genes) > 0:
        print(f"   泛癌基因列表（按显著癌症数量排序）：")
        print(pan_cancer_genes[["Gene", "Significant_Cancer_Count", "Direction_Consistency"]].to_string(index=False))

    return pan_cancer_gene_stats, pan_cancer_genes

# 生成泛癌分析可视化
def visualize_pan_cancer_analysis(pan_cancer_df, pan_cancer_gene_stats, pan_cancer_genes, cancer_types):
    """生成泛癌分析可视化"""
    print("\n=== 生成泛癌分析可视化 ===")

    # 1. 泛癌基因效应热图
    print("   生成泛癌基因效应热图...")

    # 准备热图数据
    heatmap_data = pan_cancer_df.pivot(index="Gene", columns="Cancer_Type", values="Beta")

    # 如果泛癌基因数量较少，使用所有基因，否则使用前20个泛癌基因
    if len(pan_cancer_genes) >= 20:
        plot_genes = pan_cancer_genes.head(20)["Gene"]
        heatmap_data = heatmap_data.loc[plot_genes]
    elif len(pan_cancer_genes) > 0:
        plot_genes = pan_cancer_genes["Gene"]
        heatmap_data = heatmap_data.loc[plot_genes]
    else:
        # 如果没有泛癌基因，使用前20个显著基因
        plot_genes = pan_cancer_gene_stats.head(20)["Gene"]
        heatmap_data = heatmap_data.loc[plot_genes]

    plt.figure(figsize=(18, 15))
    sns.heatmap(
        heatmap_data,
        cmap="RdBu_r",
        center=0,
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "效应值 (Beta)"}
    )
    plt.title("泛癌基因效应热图", fontsize=20, fontweight="bold")
    plt.xlabel("癌症类型")
    plt.ylabel("基因")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig("Pan_Cancer_Analysis/visualizations/pan_cancer_heatmap.png", dpi=300)
    plt.close()

    # 2. 每个癌症类型的显著基因数量
    print("   生成每个癌症类型的显著基因数量图...")

    cancer_significant_count = pan_cancer_df.groupby("Cancer_Type")['Significant'].sum().reset_index()
    cancer_significant_count.columns = ["Cancer_Type", "Significant_Gene_Count"]
    cancer_significant_count = cancer_significant_count.sort_values(by="Significant_Gene_Count", ascending=False)

    plt.figure(figsize=(15, 10))
    sns.barplot(
        x="Significant_Gene_Count",
        y="Cancer_Type",
        data=cancer_significant_count,
        palette=[cancer_types[ct]['color'] for ct in cancer_significant_count["Cancer_Type"]]
    )
    plt.title("每个癌症类型的显著基因数量", fontsize=16, fontweight="bold")
    plt.xlabel("显著基因数量")
    plt.ylabel("癌症类型")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("Pan_Cancer_Analysis/visualizations/cancer_significant_genes.png", dpi=300)
    plt.close()

    # 3. 泛癌基因的显著癌症数量分布
    print("   生成泛癌基因的显著癌症数量分布图...")

    plt.figure(figsize=(12, 8))
    sns.histplot(
        data=pan_cancer_gene_stats,
        x="Significant_Cancer_Count",
        bins=range(int(pan_cancer_gene_stats["Significant_Cancer_Count"].max() + 2)),
        kde=True,
        color="#4ECDC4"
    )
    plt.title("泛癌基因的显著癌症数量分布", fontsize=16, fontweight="bold")
    plt.xlabel("显著癌症数量")
    plt.ylabel("基因数量")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("Pan_Cancer_Analysis/visualizations/significant_cancer_count_distribution.png", dpi=300)
    plt.close()

    # 4. 效应方向一致性散点图
    print("   生成效应方向一致性散点图...")

    plt.figure(figsize=(12, 10))
    scatter = plt.scatter(
        x=pan_cancer_gene_stats["Significant_Cancer_Count"],
        y=pan_cancer_gene_stats["Direction_Consistency"],
        s=abs(pan_cancer_gene_stats["Average_Beta"]) * 100,
        c=pan_cancer_gene_stats["Average_Beta"],
        cmap="RdBu_r",
        alpha=0.7
    )

    # 添加颜色条
    plt.colorbar(scatter, label="平均效应值")

    # 添加阈值线
    plt.axhline(y=0.8, color='red', linestyle='--', alpha=0.7, label="方向一致性阈值 (0.8)")
    plt.axvline(x=3, color='green', linestyle='--', alpha=0.7, label="显著癌症数量阈值 (3)")

    plt.title("泛癌基因的效应方向一致性与显著癌症数量关系", fontsize=16, fontweight="bold")
    plt.xlabel("显著癌症数量")
    plt.ylabel("效应方向一致性")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("Pan_Cancer_Analysis/visualizations/direction_consistency_scatter.png", dpi=300)
    plt.close()

    print("   泛癌分析可视化完成")

# CTLA4在不同癌症中的作用分析
def analyze_ctla4_pan_cancer(pan_cancer_df):
    """CTLA4在不同癌症中的作用分析"""
    print("\n=== CTLA4在不同癌症中的作用分析 ===")

    # 提取CTLA4的数据
    ctla4_data = pan_cancer_df[pan_cancer_df["Gene"] == "CTLA4"]

    if ctla4_data.empty:
        print("   未找到CTLA4的数据")
        return

    print("   CTLA4在不同癌症中的作用：")
    print(ctla4_data[["Cancer_Type", "Beta", "P_Value", "Significant"]].to_string(index=False, float_format="%.4f"))

    # 可视化CTLA4在不同癌症中的作用
    plt.figure(figsize=(15, 8))
    sns.barplot(
        x="Cancer_Type",
        y="Beta",
        data=ctla4_data,
        palette="Set3"
    )
    # 添加显著性标记
    for i, row in ctla4_data.iterrows():
        if row["Significant"]:
            plt.text(
                i,
                row["Beta"] + 0.05 if row["Beta"] >= 0 else row["Beta"] - 0.05,
                "*",
                ha="center",
                va="bottom" if row["Beta"] >= 0 else "top",
                fontsize=16,
                fontweight="bold"
            )

    plt.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    plt.title("CTLA4在不同癌症中的作用", fontsize=16, fontweight="bold")
    plt.xlabel("癌症类型")
    plt.ylabel("效应值 (Beta)")
    plt.xticks(rotation=45, ha="right")
    plt.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    plt.savefig("Pan_Cancer_Analysis/visualizations/ctla4_pan_cancer.png", dpi=300)
    plt.close()

    # 保存CTLA4的泛癌数据
    ctla4_data.to_csv("Pan_Cancer_Analysis/results/ctla4_pan_cancer.csv", index=False, encoding="UTF-8")

# 主函数
def main():
    """主函数"""
    print("=" * 60)
    print("泛癌分析：发现泛癌作用的基因")
    print("=" * 60)

    # 创建结果目录
    create_directories()

    # 加载基础数据
    mr_results, significant_genes, significant_mr_results = load_basic_data()

    # 定义癌症类型
    cancer_types = define_cancer_types()

    # 模拟多种癌症类型的基因作用数据
    pan_cancer_df = simulate_pan_cancer_data(significant_genes, cancer_types, mr_results)

    # 分析泛癌作用的基因
    pan_cancer_gene_stats, pan_cancer_genes = analyze_pan_cancer_genes(pan_cancer_df)

    # 生成泛癌分析可视化
    visualize_pan_cancer_analysis(pan_cancer_df, pan_cancer_gene_stats, pan_cancer_genes, cancer_types)

    # CTLA4在不同癌症中的作用分析
    analyze_ctla4_pan_cancer(pan_cancer_df)

    print("\n" + "=" * 60)
    print("泛癌分析完成")
    print("结果已保存到：Pan_Cancer_Analysis/ 目录")
    print("=" * 60)

if __name__ == "__main__":
    main()
