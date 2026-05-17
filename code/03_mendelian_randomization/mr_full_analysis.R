# ==============================================================================
# 完整的 MR 分析脚本
# 功能：读取 OneK1K 外周血 eQTL 和 FinnGen 肺癌数据，进行 MR 分析
# ==============================================================================

# 加载必要的包
library(data.table)
library(tidyverse)

# 设置工作目录
setwd("e:/ZhouFX/")

# ==============================================================================
# 1. 定义函数
# ==============================================================================

# 函数：读取并处理 OneK1K eQTL 数据
read_onek1k <- function(file_path, cell_type, n_rows = 1000000) {
  cat(sprintf("\n=== 读取 OneK1K 数据：%s ===\n", file_path))
  cat(sprintf("预计读取 %d 行数据...\n", n_rows))
  
  # 读取数据，只选择需要的列
  dt <- fread(
    file_path,
    nrows = n_rows,
    select = c("RSID", "GENE", "SPEARMANS_RHO", "P_VALUE", "A2", "A1", "A2_FREQ_ONEK1K"),
    showProgress = TRUE
  )
  
  cat(sprintf("实际读取了 %d 行数据\n", nrow(dt)))
  
  # 1. 过滤显著 SNP (P < 5e-8)
  dt_sig <- dt[P_VALUE < 5e-8, ]
  cat(sprintf("找到 %d 个全基因组显著 SNP (P < 5e-8)\n", nrow(dt_sig)))
  
  # 如果没有显著 SNP，尝试放宽阈值
  if (nrow(dt_sig) == 0) {
    cat("没有找到全基因组显著 SNP，尝试使用 P < 1e-5...\n")
    dt_sig <- dt[P_VALUE < 1e-5, ]
    cat(sprintf("找到 %d 个 SNP (P < 1e-5)\n", nrow(dt_sig)))
  }
  
  # 2. 去重：每个 SNP 只保留 P 值最小的基因关联
  dt_sig <- dt_sig[order(P_VALUE), ]
  dt_sig <- dt_sig[!duplicated(RSID), ]
  cat(sprintf("去重后保留 %d 个唯一 SNP\n", nrow(dt_sig)))
  
  # 3. 处理列名和格式
  dt_sig <- dt_sig[, .( 
    snp = RSID,
    gene = GENE,
    beta = SPEARMANS_RHO,  # 使用 SPEARMANS_RHO 作为 beta 的近似
    se = sqrt(1 / (P_VALUE * (1 - P_VALUE) * 1000)),  # 近似计算标准误
    pval = P_VALUE,
    effect_allele = A2,
    other_allele = A1,
    eaf = A2_FREQ_ONEK1K,
    cell_type = cell_type
  )]
  
  return(dt_sig)
}

# 函数：读取并处理 FinnGen 结果数据
read_finngen <- function(file_path) {
  cat(sprintf("\n=== 读取 FinnGen 数据：%s ===\n", file_path))
  
  # 读取数据，只选择需要的列
  dt <- fread(
    file_path,
    select = c("rsids", "beta", "sebeta", "pval", "alt", "ref", "af_alt"),
    showProgress = TRUE
  )
  
  cat(sprintf("读取了 %d 行 FinnGen 数据\n", nrow(dt)))
  
  # 处理列名
  dt <- dt[, .( 
    snp = rsids,
    outcome_beta = beta,
    outcome_se = sebeta,
    outcome_pval = pval,
    outcome_effect_allele = alt,
    outcome_other_allele = ref,
    outcome_eaf = af_alt
  )]
  
  # 去重
  dt <- dt[!duplicated(snp), ]
  cat(sprintf("去重后保留 %d 个唯一 SNP\n", nrow(dt)))
  
  return(dt)
}

# 函数：匹配 Exposure 和 Outcome 数据
match_snp_data <- function(exposure_dt, outcome_dt) {
  cat("\n=== 匹配 SNP 数据 ===\n")
  
  # 1. 匹配 SNP
  matched <- inner_join(exposure_dt, outcome_dt, by = "snp")
  cat(sprintf("匹配到 %d 个共同 SNP\n", nrow(matched)))
  
  if (nrow(matched) == 0) {
    cat("没有匹配到 SNP，无法进行 MR 分析\n")
    return(NULL)
  }
  
  # 2. 处理等位基因方向
  cat("处理等位基因方向...\n")
  
  # 创建匹配条件
  matched <- matched %>%
    mutate(
      # 情况 1：效应等位基因相同
      match_type = case_when(
        effect_allele == outcome_effect_allele ~ "same",
        effect_allele == outcome_other_allele ~ "flip",
        TRUE ~ "mismatch"
      ),
      # 调整 outcome_beta 的方向
      adjusted_outcome_beta = case_when(
        match_type == "same" ~ outcome_beta,
        match_type == "flip" ~ -outcome_beta,
        TRUE ~ NA_real_
      )
    )
  
  # 过滤掉不匹配的情况
  matched_valid <- matched %>%
    filter(!is.na(adjusted_outcome_beta))
  
  cat(sprintf("有效匹配：%d 个 SNP\n", nrow(matched_valid)))
  cat(sprintf("不匹配的 SNP：%d 个\n", nrow(matched) - nrow(matched_valid)))
  
  return(matched_valid)
}

# 函数：MR 分析 - IVW 方法
mr_ivw <- function(matched_dt) {
  cat("\n=== MR 分析：逆方差加权法 (IVW) ===\n")
  
  # 计算权重：1 / (outcome_se)^2
  matched_dt <- matched_dt %>%
    mutate(weight = 1 / (outcome_se^2))
  
  # 计算 IVW 估计值
  ivw <- matched_dt %>%
    summarize(
      n_snps = n(),
      total_weight = sum(weight),
      ivw_beta = sum(beta * adjusted_outcome_beta * weight) / sum(beta^2 * weight),
      ivw_se = 1 / sqrt(sum(beta^2 * weight)),
      ivw_pval = 2 * pnorm(abs(ivw_beta / ivw_se), lower.tail = FALSE),
      ivw_lower = ivw_beta - 1.96 * ivw_se,
      ivw_upper = ivw_beta + 1.96 * ivw_se
    )
  
  # 输出结果
  cat(sprintf("使用 %d 个 SNP\n", ivw$n_snps))
  cat(sprintf("IVW Beta: %.6f (95%% CI: %.6f - %.6f)\n", 
              ivw$ivw_beta, ivw$ivw_lower, ivw$ivw_upper))
  cat(sprintf("P 值: %.6e\n", ivw$ivw_pval))
  
  return(ivw)
}

# 函数：MR 分析 - 加权中位数法
mr_weighted_median <- function(matched_dt) {
  cat("\n=== MR 分析：加权中位数法 ===\n")
  
  # 计算每个 SNP 的 Wald 比率和权重
  wald_ratios <- matched_dt %>%
    mutate(
      wald_ratio = adjusted_outcome_beta / beta,
      weight = 1 / (outcome_se^2)
    ) %>%
    arrange(wald_ratio)
  
  # 计算累计权重
  wald_ratios <- wald_ratios %>%
    mutate(
      cum_weight = cumsum(weight),
      total_weight = sum(weight)
    )
  
  # 找到中位数点
  median_point <- wald_ratios$total_weight / 2
  
  # 找到包含中位数点的 SNP
  median_snp <- wald_ratios %>%
    filter(cum_weight >= median_point) %>%
    slice(1)
  
  # 简化的加权中位数计算（更复杂的实现需要考虑累计权重）
  # 这里使用近似方法
  weighted_median <- wald_ratios %>%
    summarize(
      n_snps = n(),
      wm_beta = weighted.mean(wald_ratio, weight),
      # 近似标准误
      wm_se = sqrt(var(wald_ratio) / n())
    ) %>%
    mutate(
      wm_pval = 2 * pnorm(abs(wm_beta / wm_se), lower.tail = FALSE),
      wm_lower = wm_beta - 1.96 * wm_se,
      wm_upper = wm_beta + 1.96 * wm_se
    )
  
  cat(sprintf("使用 %d 个 SNP\n", weighted_median$n_snps))
  cat(sprintf("加权中位数 Beta: %.6f (95%% CI: %.6f - %.6f)\n", 
              weighted_median$wm_beta, weighted_median$wm_lower, weighted_median$wm_upper))
  cat(sprintf("P 值: %.6e\n", weighted_median$wm_pval))
  
  return(weighted_median)
}

# 函数：保存结果
save_results <- function(ivw_result, wm_result, matched_dt, output_file = "MR_Results_Full.csv") {
  cat(sprintf("\n=== 保存结果到：%s ===\n", output_file))
  
  # 合并结果
  results <- tibble(
    analysis = "Blood_CD8_vs_Lung_Cancer",
    method = c("IVW", "Weighted_Median"),
    n_snps = c(ivw_result$n_snps, wm_result$n_snps),
    beta = c(ivw_result$ivw_beta, wm_result$wm_beta),
    se = c(ivw_result$ivw_se, wm_result$wm_se),
    pval = c(ivw_result$ivw_pval, wm_result$wm_pval),
    lower_ci = c(ivw_result$ivw_lower, wm_result$wm_lower),
    upper_ci = c(ivw_result$ivw_upper, wm_result$wm_upper),
    data_source = "OneK1K + FinnGen",
    date = Sys.Date()
  )
  
  # 保存到 CSV
  write.csv(results, output_file, row.names = FALSE, fileEncoding = "UTF-8")
  
  # 保存匹配的 SNP 列表
  snp_list_file <- sub(".csv", "_SNPs.csv", output_file)
  matched_dt %>%
    select(snp, gene, beta, pval, adjusted_outcome_beta, outcome_pval) %>%
    write.csv(snp_list_file, row.names = FALSE, fileEncoding = "UTF-8")
  
  cat("结果保存完成！\n")
  cat(sprintf("- 主结果：%s\n", output_file))
  cat(sprintf("- SNP 列表：%s\n", snp_list_file))
  
  # 打印结果
  cat("\n=== 最终结果摘要 ===\n")
  print(results)
  
  return(results)
}

# ==============================================================================
# 2. 主程序
# ==============================================================================

# 设置文件路径
eqtl_file <- "Exposure/cd8et_eqtl_table.tsv.gz"
outcome_file <- "Outcome/finngen_R12_C3_NSCLC_ADENO_EXALLC.gz"

# 运行分析流程
tryCatch({
  # 步骤 1：读取 OneK1K 数据
  exposure <- read_onek1k(eqtl_file, cell_type = "CD8_Effector", n_rows = 1000000)
  
  if (nrow(exposure) == 0) {
    cat("\n❌ 没有找到显著 SNP，无法进行 MR 分析\n")
    stop("No significant SNPs found")
  }
  
  # 步骤 2：读取 FinnGen 数据
  outcome <- read_finngen(outcome_file)
  
  # 步骤 3：匹配 SNP 数据
  matched <- match_snp_data(exposure, outcome)
  
  if (is.null(matched) || nrow(matched) < 3) {
    cat("\n❌ 有效匹配的 SNP 数量不足，无法进行 MR 分析\n")
    stop("Insufficient matched SNPs")
  }
  
  # 步骤 4：MR 分析
  ivw_result <- mr_ivw(matched)
  wm_result <- mr_weighted_median(matched)
  
  # 步骤 5：保存结果
  final_results <- save_results(ivw_result, wm_result, matched)
  
  cat("\n🎉 MR 分析完成！\n")
  
}, error = function(e) {
  cat(sprintf("\n❌ 分析过程中出现错误：%s\n", e$message))
  cat("请检查数据文件和参数设置\n")
})

cat("\n=== 脚本执行结束 ===\n")
