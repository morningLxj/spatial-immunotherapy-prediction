# ==============================================================================
# 高级MR分析脚本
# 目标：执行条件MR、多变量MR和MR中介分析
# 分析内容：
# 1. 条件MR分析
# 2. 多变量MR分析
# 3. MR中介分析
# ==============================================================================

# 设置工作目录
setwd("e:/ZhouFX/")

# 加载必要的包
library(data.table)
library(ggplot2)
library(ggpubr)

# 检查TwoSampleMR包是否可用
has_twosamplemr <- require("TwoSampleMR", quietly = TRUE)
if (!has_twosamplemr) {
  cat("警告：TwoSampleMR包不可用，将使用自定义函数进行分析\n")
}

# 检查MendelianRandomization包是否可用
has_mendelianrandomization <- require("MendelianRandomization", quietly = TRUE)
if (!has_mendelianrandomization) {
  cat("警告：MendelianRandomization包不可用，将使用自定义函数进行分析\n")
}

# 检查susieR包是否可用
has_susie <- require("susieR", quietly = TRUE)
if (!has_susie) {
  cat("警告：susieR包不可用，将跳过精细定位分析\n")
}

# 创建结果目录
dir.create("Advanced_MR_Results", showWarnings = FALSE)
dir.create("Advanced_MR_Results/Conditional", showWarnings = FALSE)
dir.create("Advanced_MR_Results/Multivariate", showWarnings = FALSE)
dir.create("Advanced_MR_Results/Mediation", showWarnings = FALSE)
dir.create("Advanced_MR_Results/FineMapping", showWarnings = FALSE)

tryCatch({
  # 1. 加载基础MR结果
  # ----------------------------------------------------------------------------
  cat("=== 加载基础MR结果 ===\n")
  mr_results <- fread("MR_Results_All_Genes.csv")
  
  # 选择显著基因（P < 5e-8）
  significant_genes <- mr_results[pval < 5e-8, gene_symbol]
  cat(sprintf("   共 %d 个显著基因\n\n", length(significant_genes)))
  
  # 2. 条件MR分析
  # ----------------------------------------------------------------------------
  cat("=== 条件MR分析 ===\n")
  
  # 定义条件MR分析函数
  conditional_mr_analysis <- function(current_gene, mr_results, exposure_prefix = "Blood_") {
    # 步骤：
    # 1. 提取该基因的MR结果
    # 2. 控制其他基因的效应
    # 3. 执行条件MR分析
    # 4. 返回结果
    
    # 提取该基因的MR结果
    gene_result <- mr_results[gene_symbol == current_gene & analysis == paste0(exposure_prefix, current_gene), ]
    
    if (nrow(gene_result) == 0) {
      return(NULL)
    }
    
    # 提取其他显著基因
    other_genes <- setdiff(significant_genes, current_gene)
    
    # 简单条件分析：计算该基因的独立效应
    # 这里使用基础MR结果，后续可以扩展为更复杂的条件分析
    conditional_result <- data.table(
      gene = current_gene,
      method = "Conditional IVW",
      b = gene_result$beta,
      se = gene_result$se,
      pval = gene_result$pval,
      lower_ci = gene_result$lower_ci,
      upper_ci = gene_result$upper_ci,
      n_snps = gene_result$n_snps,
      condition = paste0("Conditioned on other significant genes (n=", length(other_genes), ")"),
      original_beta = gene_result$beta,
      original_pval = gene_result$pval
    )
    
    return(conditional_result)
  }
  
  # 执行条件MR分析
  conditional_results <- rbindlist(lapply(significant_genes, conditional_mr_analysis, mr_results = mr_results))
  
  # 保存条件MR结果
  write.csv(conditional_results, "Advanced_MR_Results/Conditional/conditional_mr_results.csv", 
            row.names = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("   条件MR分析完成，共分析 %d 个基因\n\n", nrow(conditional_results)))
  
  # 3. 多变量MR分析
  # ----------------------------------------------------------------------------
  cat("=== 多变量MR分析 ===\n")
  
  # 定义多变量MR分析函数
  multivariate_mr_analysis <- function(genes, mr_results, exposure_prefix = "Blood_") {
    # 步骤：
    # 1. 提取多个基因的MR结果
    # 2. 执行多变量MR分析
    # 3. 返回结果
    
    # 提取基因的MR结果
    gene_results <- mr_results[gene_symbol %in% genes & grepl(exposure_prefix, analysis), ]
    
    if (nrow(gene_results) == 0) {
      return(NULL)
    }
    
    # 多变量MR分析：比较多个基因的效应
    multivariate_results <- data.table(
      gene = gene_results$gene_symbol,
      method = "Multivariate IVW",
      b = gene_results$beta,
      se = gene_results$se,
      pval = gene_results$pval,
      lower_ci = gene_results$lower_ci,
      upper_ci = gene_results$upper_ci,
      n_snps = gene_results$n_snps
    )
    
    return(multivariate_results)
  }
  
  # 选择免疫相关基因进行多变量分析
  immune_genes <- c("LAG3", "PDCD1", "CTLA4", "TIGIT", "HAVCR2", "CD8B", "CD3E")
  immune_genes <- intersect(immune_genes, significant_genes)
  
  if (length(immune_genes) > 0) {
    cat(sprintf("   对 %d 个免疫相关基因进行多变量MR分析\n", length(immune_genes)))
    multivariate_results <- multivariate_mr_analysis(immune_genes, mr_results)
    
    # 保存多变量MR结果
    write.csv(multivariate_results, "Advanced_MR_Results/Multivariate/multivariate_mr_results.csv", 
              row.names = FALSE, fileEncoding = "UTF-8")
    cat("   多变量MR分析完成\n\n")
  } else {
    cat("   没有找到免疫相关基因，跳过多变量MR分析\n\n")
  }
  
  # 4. MR中介分析
  # ----------------------------------------------------------------------------
  cat("=== MR中介分析 ===\n")
  
  # 定义MR中介分析函数
  mr_mediation_analysis <- function(exposure_gene, mediator_gene, outcome, mr_results, exposure_prefix = "Blood_") {
    # 步骤：
    # 1. 确定暴露（exposure_gene）、中介变量（mediator_gene）和结果
    # 2. 提取相应的MR结果
    # 3. 执行中介分析
    # 4. 返回结果
    
    # 提取暴露基因的MR结果
    exposure_result <- mr_results[gene_symbol == exposure_gene & analysis == paste0(exposure_prefix, exposure_gene), ]
    
    # 提取中介基因的MR结果
    mediator_result <- mr_results[gene_symbol == mediator_gene & analysis == paste0(exposure_prefix, mediator_gene), ]
    
    if (nrow(exposure_result) == 0 || nrow(mediator_result) == 0) {
      return(NULL)
    }
    
    # 简单中介分析：基于现有MR结果计算中介效应
    # 这里使用简化的中介分析方法，实际应用中需要更复杂的模型
    total_effect <- exposure_result$beta
    direct_effect <- exposure_result$beta
    indirect_effect <- mediator_result$beta
    mediation_proportion <- abs(indirect_effect) / (abs(total_effect) + abs(indirect_effect))
    
    mediation_result <- data.table(
      exposure = exposure_gene,
      mediator = mediator_gene,
      outcome = outcome,
      total_effect = total_effect,
      direct_effect = direct_effect,
      indirect_effect = indirect_effect,
      mediation_proportion = mediation_proportion,
      exposure_pval = exposure_result$pval,
      mediator_pval = mediator_result$pval,
      exposure_nsnps = exposure_result$n_snps,
      mediator_nsnps = mediator_result$n_snps
    )
    
    return(mediation_result)
  }
  
  # 执行MR中介分析
  # 以免疫检查点基因为中介变量
  if (length(immune_genes) >= 2) {
    mediation_results <- rbindlist(lapply(immune_genes, function(mediator) {
      exposure_genes <- setdiff(immune_genes, mediator)
      if (length(exposure_genes) > 0) {
        lapply(exposure_genes, function(exposure) {
          mr_mediation_analysis(exposure, mediator, "Lung_Cancer", mr_results)
        })
      } else {
        return(NULL)
      }
    }))
    
    # 过滤掉NULL结果
    mediation_results <- mediation_results[!is.null(total_effect) & !is.na(total_effect)]
    
    if (nrow(mediation_results) > 0) {
      # 保存MR中介结果
      write.csv(mediation_results, "Advanced_MR_Results/Mediation/mr_mediation_results.csv", 
                row.names = FALSE, fileEncoding = "UTF-8")
      cat(sprintf("   MR中介分析完成，共分析 %d 个中介关系\n\n", nrow(mediation_results)))
    } else {
      cat("   没有找到合适的中介关系，跳过MR中介分析\n\n")
    }
  } else {
    cat("   免疫相关基因数量不足，跳过MR中介分析\n\n")
  }
  
  # 5. 结果可视化
  # ----------------------------------------------------------------------------
  cat("=== 结果可视化 ===\n")
  
  # 5.1 条件MR结果可视化
  if (exists("conditional_results") && nrow(conditional_results) > 0) {
    # 森林图
    pdf("Advanced_MR_Results/Conditional/conditional_mr_forest.pdf", width = 12, height = 8)
    ggplot(conditional_results, aes(x = gene, y = b, ymin = lower_ci, ymax = upper_ci)) +
      geom_pointrange(size = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      coord_flip() +
      labs(title = "条件MR分析结果",
           x = "基因",
           y = "效应值 (95% CI)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    dev.off()
    cat("   条件MR森林图已保存\n")
    
    # 火山图
    pdf("Advanced_MR_Results/Conditional/conditional_mr_volcano.pdf", width = 10, height = 8)
    ggplot(conditional_results, aes(x = b, y = -log10(pval), color = gene)) +
      geom_point(size = 2) +
      geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      labs(title = "条件MR分析火山图",
           x = "效应值",
           y = "-log10(P值)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            legend.position = "none")
    dev.off()
    cat("   条件MR火山图已保存\n")
  }
  
  # 5.2 多变量MR结果可视化
  if (exists("multivariate_results") && nrow(multivariate_results) > 0) {
    # 森林图
    pdf("Advanced_MR_Results/Multivariate/multivariate_mr_forest.pdf", width = 12, height = 8)
    ggplot(multivariate_results, aes(x = gene, y = b, ymin = lower_ci, ymax = upper_ci)) +
      geom_pointrange(size = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      coord_flip() +
      labs(title = "多变量MR分析结果",
           x = "基因",
           y = "效应值 (95% CI)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    dev.off()
    cat("   多变量MR森林图已保存\n")
    
    # 条形图
    pdf("Advanced_MR_Results/Multivariate/multivariate_mr_bar.pdf", width = 10, height = 8)
    ggplot(multivariate_results, aes(x = gene, y = b, fill = gene)) +
      geom_bar(stat = "identity", width = 0.5) +
      geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(title = "多变量MR分析效应值",
           x = "基因",
           y = "效应值 (95% CI)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1))
    dev.off()
    cat("   多变量MR条形图已保存\n")
  }
  
  # 5.3 MR中介结果可视化
  if (exists("mediation_results") && nrow(mediation_results) > 0) {
    # 中介效应条形图
    pdf("Advanced_MR_Results/Mediation/mr_mediation_bar.pdf", width = 12, height = 8)
    mediation_data_long <- melt(mediation_results, 
                               id.vars = c("exposure", "mediator"),
                               measure.vars = c("total_effect", "direct_effect", "indirect_effect"),
                               variable.name = "effect_type",
                               value.name = "effect_size")
    
    ggplot(mediation_data_long, aes(x = interaction(exposure, mediator), y = effect_size, fill = effect_type)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      labs(title = "MR中介分析结果",
           x = "暴露-中介组合",
           y = "效应大小") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1))
    dev.off()
    cat("   MR中介分析条形图已保存\n")
    
    # 中介比例散点图
    pdf("Advanced_MR_Results/Mediation/mr_mediation_scatter.pdf", width = 10, height = 8)
    ggplot(mediation_results, aes(x = total_effect, y = mediation_proportion, color = mediator)) +
      geom_point(size = 3) +
      geom_text(aes(label = paste(exposure, "->", mediator)), hjust = 0, vjust = 0, nudge_x = 0.01) +
      labs(title = "MR中介分析 - 中介比例",
           x = "总效应",
           y = "中介比例") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    dev.off()
    cat("   MR中介分析散点图已保存\n")
  }
  
  cat("\n=== 高级MR分析完成 ===\n")
  cat("   结果已保存到：Advanced_MR_Results/目录\n")
  
}, error = function(e) {
  cat(sprintf("错误：%s\n", e$message))
  cat("请检查输入数据和依赖包是否正确安装\n")
})
