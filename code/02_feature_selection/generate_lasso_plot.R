# generate_lasso_plot.R
# 专门用于生成Figure S2A所需的LASSO交叉验证图
# 策略：使用全基因组中方差最大的Top N基因进行LASSO筛选

# 创建目录
dir.create("Final_Analysis/Plots", recursive = TRUE, showWarnings = FALSE)
dir.create("Final_Analysis/Results", recursive = TRUE, showWarnings = FALSE)

# 加载包
required_packages <- c("codetools", "lattice", "Matrix", "survival", "glmnet")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("正在读取数据...\n")

# 读取表达数据
if (!file.exists("TCGA/TCGA-LUAD.star_tpm.tsv")) {
  stop("找不到 TCGA/TCGA-LUAD.star_tpm.tsv")
}
# 只读取前5000行以加快速度（如果文件太大），或者读取全部但只选Top
# 这里我们尝试读取全部，但如果内存不够会报错。
# 考虑到效率，我们读取全部
expr_df <- read.table("TCGA/TCGA-LUAD.star_tpm.tsv", header = TRUE, sep = "\t", check.names = FALSE)
rownames(expr_df) <- expr_df[, 1]
expr_df <- expr_df[, -1]

# 读取生存数据
if (!file.exists("TCGA/TCGA-LUAD.survival.tsv")) {
  stop("找不到 TCGA/TCGA-LUAD.survival.tsv")
}
survival_data <- read.table("TCGA/TCGA-LUAD.survival.tsv", header = TRUE, sep = "\t")

cat("数据加载完成。样本数:", ncol(expr_df), "基因数:", nrow(expr_df), "\n")

# 转置表达矩阵：行=样本
validated_expr <- t(expr_df)

# 匹配样本
common_samples <- intersect(rownames(validated_expr), survival_data$sample)
cat("匹配样本数:", length(common_samples), "\n")

if (length(common_samples) > 20) {
  # 对齐数据
  validated_expr_sub <- validated_expr[common_samples, ]
  survival_data_sub <- survival_data[match(common_samples, survival_data$sample), ]

  X <- as.matrix(validated_expr_sub)
  y <- survival_data_sub$OS

  # 处理NA
  X[is.na(X)] <- 0

  # 特征筛选：选择方差最大的前1000个基因
  cat("计算基因方差并筛选Top 1000...\n")
  # 此时X的列是基因
  gene_vars <- apply(X, 2, var)
  # 移除方差为0的
  gene_vars <- gene_vars[gene_vars > 0]
  # 排序
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:1000]

  X_filtered <- X[, top_genes]
  cat("筛选后维度:", dim(X_filtered)[1], "x", dim(X_filtered)[2], "\n")

  # 标准化
  X_scaled <- scale(X_filtered)

  cat("运行 LASSO cv.glmnet...\n")
  # alpha=1 for LASSO
  cv_lasso <- cv.glmnet(X_scaled, y, family = "binomial", alpha = 1, nfolds = 10)

  cat("保存图片...\n")
  pdf("Final_Analysis/Plots/lasso_cv_plot.pdf", width = 10, height = 8)
  plot(cv_lasso, main = "LASSO Cross-Validation (Top 1000 Var Genes)")
  abline(v = log(cv_lasso$lambda.min), col = "red", lty = 2)
  abline(v = log(cv_lasso$lambda.1se), col = "blue", lty = 2)
  legend("topright", legend = c("Lambda.min", "Lambda.1se"), col = c("red", "blue"), lty = 2)
  dev.off()

  png("Final_Analysis/Plots/figS2A_lasso_cv.png", width = 4800, height = 3600, res = 600)
  plot(cv_lasso, main = "")
  abline(v = log(cv_lasso$lambda.min), col = "red", lty = 2)
  abline(v = log(cv_lasso$lambda.1se), col = "blue", lty = 2)
  legend("bottomright", legend = c("Lambda.min", "Lambda.1se"), col = c("red", "blue"), lty = 2, cex = 1.2)
  dev.off()

  cat("LASSO图生成完成。\n")

  # 保存非零系数（可选，为了完整性）
  best_coef <- coef(cv_lasso, s = "lambda.min")
  non_zero_coef <- data.frame(
    Gene = row.names(best_coef)[best_coef[, 1] != 0],
    Coefficient = best_coef[best_coef[, 1] != 0, 1]
  )
  # 移除截距
  non_zero_coef <- non_zero_coef[non_zero_coef$Gene != "(Intercept)", ]

  if (nrow(non_zero_coef) > 0) {
    write.csv(non_zero_coef, "Final_Analysis/Results/lasso_selected_features.csv", row.names = FALSE)
    cat("选中的特征已保存。\n")
  }
} else {
  cat("样本不足，无法生成LASSO图。\n")
}
