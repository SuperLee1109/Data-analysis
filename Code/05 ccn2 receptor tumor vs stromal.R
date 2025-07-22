# CCN2受体基因在敏感及耐药组之间的差异分析
# 专注于肿瘤细胞和间质细胞

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)
library(ggpubr)
library(rstatix)

# CCN2受体基因列表
ccn2_receptors <- c("ITGAV", "ITGB3", "ITGA5", "LRP1", "DDR2")

# 1. 检查数据和必要信息
# 请根据你的实际数据修改以下变量名：
# - 数据对象名：scRNA
# - 敏感/耐药信息列名：通常是 "response", "sensitivity", "drug_response" 等
# - 细胞类型列名：通常是 "cell_type", "celltype", "annotation" 等

# 检查可用的meta.data列
cat("可用的meta.data列：\n")
print(colnames(scRNA@meta.data))

# 根据用户信息指定列名
response_col <- "group"  # 敏感/耐药信息列
cell_type_col <- "cell_type"  # 细胞类型列名（如果不同请修改）

# 检查敏感/耐药分组
if(response_col %in% colnames(scRNA@meta.data)) {
  cat("\n敏感/耐药分组情况：\n")
  print(table(scRNA@meta.data[[response_col]]))
} else {
  cat("错误：找不到敏感/耐药信息列，请检查列名\n")
  stop("请修改 response_col 变量为正确的列名")
}

# 检查细胞类型
if(cell_type_col %in% colnames(scRNA@meta.data)) {
  cat("\n细胞类型分布：\n")
  print(table(scRNA@meta.data[[cell_type_col]]))
} else {
  cat("错误：找不到细胞类型列，请检查列名\n")
  stop("请修改 cell_type_col 变量为正确的列名")
}

# 2. 标准化敏感/耐药标签
# 根据用户提供的信息指定标签
sensitive_labels <- "S"  # 敏感标签
resistant_labels <- "R"  # 耐药标签

cat("\n敏感标签：S\n")
cat("耐药标签：R\n")

# 创建标准化的响应分组
scRNA$drug_response <- ifelse(scRNA@meta.data[[response_col]] == "S", 
                               "Sensitive",
                               ifelse(scRNA@meta.data[[response_col]] == "R", 
                                      "Resistant", "Other"))

cat("\n标准化后的分组：\n")
print(table(scRNA$drug_response))

# 3. 识别肿瘤细胞和间质细胞
# 根据用户提供的信息指定细胞类型
tumor_cells <- "tumor"  # 肿瘤细胞类型
stromal_cells <- "fibroblast"  # 间质细胞类型

cat("\n指定的肿瘤细胞类型：tumor\n")
cat("指定的间质细胞类型：fibroblast\n")

# 创建细胞大类分组
scRNA$cell_major_type <- ifelse(scRNA@meta.data[[cell_type_col]] == "tumor", 
                                 "Tumor",
                                 ifelse(scRNA@meta.data[[cell_type_col]] == "fibroblast", 
                                        "Stromal", "Other"))

# 4. 检查基因存在性
genes_present <- ccn2_receptors[ccn2_receptors %in% rownames(scRNA)]
genes_missing <- ccn2_receptors[!ccn2_receptors %in% rownames(scRNA)]

cat("\n存在的CCN2受体基因:", paste(genes_present, collapse = ", "), "\n")
if(length(genes_missing) > 0) {
  cat("缺失的基因:", paste(genes_missing, collapse = ", "), "\n")
}

# 5. 创建分析用的子集数据
# 只保留有明确敏感/耐药信息和肿瘤/间质细胞的数据
analysis_data <- subset(scRNA, 
                        drug_response %in% c("Sensitive", "Resistant") & 
                          cell_major_type %in% c("Tumor", "Stromal"))

cat("\n分析数据概况：\n")
cat("总细胞数：", ncol(analysis_data), "\n")
print(table(analysis_data$drug_response, analysis_data$cell_major_type))

# 6. 敏感vs耐药的表达比较分析

# 6.1 在肿瘤细胞中的比较
cat("\n=== 肿瘤细胞中CCN2受体基因敏感vs耐药分析 ===\n")

tumor_data <- subset(analysis_data, cell_major_type == "Tumor")
tumor_comparison_plots <- list()
tumor_stats <- data.frame()

for(gene in genes_present) {
  # 获取表达数据
  sensitive_expr <- GetAssayData(tumor_data)[gene, tumor_data$drug_response == "Sensitive"]
  resistant_expr <- GetAssayData(tumor_data)[gene, tumor_data$drug_response == "Resistant"]
  
  # 统计分析
  sensitive_mean <- mean(sensitive_expr)
  resistant_mean <- mean(resistant_expr)
  fold_change <- resistant_mean / (sensitive_mean + 1e-10)  # 避免除零
  
  # Wilcoxon检验
  if(length(sensitive_expr) > 10 && length(resistant_expr) > 10) {
    test_result <- wilcox.test(sensitive_expr, resistant_expr)
    p_value <- test_result$p.value
  } else {
    p_value <- NA
  }
  
  # 存储统计结果
  tumor_stats <- rbind(tumor_stats, data.frame(
    Gene = gene,
    Cell_Type = "Tumor",
    Sensitive_Mean = sensitive_mean,
    Resistant_Mean = resistant_mean,
    Fold_Change = fold_change,
    P_Value = p_value
  ))
  
  # 创建小提琴图
  plot_data <- data.frame(
    Expression = c(sensitive_expr, resistant_expr),
    Response = rep(c("Sensitive", "Resistant"), 
                   c(length(sensitive_expr), length(resistant_expr)))
  )
  
  p <- ggplot(plot_data, aes(x = Response, y = Expression, fill = Response)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
    scale_fill_manual(values = c("Sensitive" = "#4CAF50", "Resistant" = "#F44336")) +
    labs(title = paste(gene, "in Tumor Cells"),
         subtitle = paste("p =", ifelse(is.na(p_value), "N/A", 
                                        sprintf("%.2e", p_value))),
         x = "Drug Response", y = "Expression Level") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # 添加统计显著性标记
  if(!is.na(p_value) && p_value < 0.05) {
    significance <- ifelse(p_value < 0.001, "***", 
                           ifelse(p_value < 0.01, "**", "*"))
    p <- p + annotate("text", x = 1.5, y = max(plot_data$Expression) * 0.9, 
                      label = significance, size = 6)
  }
  
  tumor_comparison_plots[[gene]] <- p
}

# 组合肿瘤细胞比较图
tumor_combined_plot <- wrap_plots(tumor_comparison_plots, ncol = 3)
print(tumor_combined_plot)

# 6.2 在间质细胞中的比较
cat("\n=== 间质细胞中CCN2受体基因敏感vs耐药分析 ===\n")

stromal_data <- subset(analysis_data, cell_major_type == "Stromal")
stromal_comparison_plots <- list()
stromal_stats <- data.frame()

for(gene in genes_present) {
  # 获取表达数据
  sensitive_expr <- GetAssayData(stromal_data)[gene, stromal_data$drug_response == "Sensitive"]
  resistant_expr <- GetAssayData(stromal_data)[gene, stromal_data$drug_response == "Resistant"]
  
  # 统计分析
  sensitive_mean <- mean(sensitive_expr)
  resistant_mean <- mean(resistant_expr)
  fold_change <- resistant_mean / (sensitive_mean + 1e-10)
  
  # Wilcoxon检验
  if(length(sensitive_expr) > 10 && length(resistant_expr) > 10) {
    test_result <- wilcox.test(sensitive_expr, resistant_expr)
    p_value <- test_result$p.value
  } else {
    p_value <- NA
  }
  
  # 存储统计结果
  stromal_stats <- rbind(stromal_stats, data.frame(
    Gene = gene,
    Cell_Type = "Stromal",
    Sensitive_Mean = sensitive_mean,
    Resistant_Mean = resistant_mean,
    Fold_Change = fold_change,
    P_Value = p_value
  ))
  
  # 创建小提琴图
  plot_data <- data.frame(
    Expression = c(sensitive_expr, resistant_expr),
    Response = rep(c("Sensitive", "Resistant"), 
                   c(length(sensitive_expr), length(resistant_expr)))
  )
  
  p <- ggplot(plot_data, aes(x = Response, y = Expression, fill = Response)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
    scale_fill_manual(values = c("Sensitive" = "#4CAF50", "Resistant" = "#F44336")) +
    labs(title = paste(gene, "in Stromal Cells"),
         subtitle = paste("p =", ifelse(is.na(p_value), "N/A", 
                                        sprintf("%.2e", p_value))),
         x = "Drug Response", y = "Expression Level") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # 添加统计显著性标记
  if(!is.na(p_value) && p_value < 0.05) {
    significance <- ifelse(p_value < 0.001, "***", 
                           ifelse(p_value < 0.01, "**", "*"))
    p <- p + annotate("text", x = 1.5, y = max(plot_data$Expression) * 0.9, 
                      label = significance, size = 6)
  }
  
  stromal_comparison_plots[[gene]] <- p
}

# 组合间质细胞比较图
stromal_combined_plot <- wrap_plots(stromal_comparison_plots, ncol = 3)
print(stromal_combined_plot)

# 7. 综合统计结果
all_stats <- rbind(tumor_stats, stromal_stats)

# 添加FDR校正
all_stats$FDR <- p.adjust(all_stats$P_Value, method = "fdr")

# 添加显著性标记
all_stats$Significance <- ifelse(is.na(all_stats$FDR), "N/A",
                                 ifelse(all_stats$FDR < 0.001, "***",
                                        ifelse(all_stats$FDR < 0.01, "**",
                                               ifelse(all_stats$FDR < 0.05, "*", "ns"))))

cat("\n=== 综合统计结果 ===\n")
print(all_stats)

# 8. 创建综合热图显示fold change
library(reshape2)
library(ComplexHeatmap)

# 准备热图数据
heatmap_data <- all_stats %>%
  select(Gene, Cell_Type, Fold_Change) %>%
  reshape2::dcast(Gene ~ Cell_Type, value.var = "Fold_Change")

rownames(heatmap_data) <- heatmap_data$Gene
heatmap_data$Gene <- NULL

# 创建显著性注释
sig_data <- all_stats %>%
  select(Gene, Cell_Type, Significance) %>%
  reshape2::dcast(Gene ~ Cell_Type, value.var = "Significance")

rownames(sig_data) <- sig_data$Gene
sig_data$Gene <- NULL

# 绘制热图
if(require(ComplexHeatmap, quietly = TRUE)) {
  library(ComplexHeatmap)
  library(circlize)
  
  # 创建颜色映射
  col_fun <- colorRamp2(c(0.5, 1, 2), c("blue", "white", "red"))
  
  # 创建热图
  ht <- Heatmap(as.matrix(heatmap_data),
                name = "Fold Change\n(Resistant/Sensitive)",
                col = col_fun,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  if(!is.na(sig_data[i, j]) && sig_data[i, j] != "ns" && sig_data[i, j] != "N/A") {
                    grid.text(sig_data[i, j], x, y, gp = gpar(fontsize = 10, fontface = "bold"))
                  }
                },
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                column_title = "CCN2 Receptors: Resistant vs Sensitive Expression",
                row_title = "Genes",
                heatmap_legend_param = list(title_position = "topcenter"))
  
  print(ht)
} else {
  cat("ComplexHeatmap包未安装，跳过热图绘制\n")
}

# 9. UMAP图显示敏感耐药差异
# 在UMAP上显示敏感/耐药分组
response_umap <- DimPlot(analysis_data, 
                         group.by = "drug_response",
                         split.by = "cell_major_type",
                         cols = c("Sensitive" = "#4CAF50", "Resistant" = "#F44336")) +
  ggtitle("Drug Response in Tumor vs Stromal Cells")

print(response_umap)

# 10. 保存结果
# 创建输出目录
if(!dir.exists("CCN2_sensitivity_analysis")) {
  dir.create("CCN2_sensitivity_analysis")
}

# 保存统计结果
write.csv(all_stats, "CCN2_sensitivity_analysis/CCN2_receptors_sensitivity_stats.csv", 
          row.names = FALSE)

# 保存图片
ggsave("CCN2_sensitivity_analysis/Tumor_Sensitivity_Comparison.pdf", 
       tumor_combined_plot, width = 15, height = 10, dpi = 300)

ggsave("CCN2_sensitivity_analysis/Stromal_Sensitivity_Comparison.pdf", 
       stromal_combined_plot, width = 15, height = 10, dpi = 300)

ggsave("CCN2_sensitivity_analysis/Response_UMAP.pdf", 
       response_umap, width = 12, height = 6, dpi = 300)

cat("\n=== 分析完成 ===\n")
cat("结果已保存到 CCN2_sensitivity_analysis/ 目录\n")
cat("主要发现：\n")

# 总结显著差异
significant_results <- all_stats[!is.na(all_stats$FDR) & all_stats$FDR < 0.05, ]
if(nrow(significant_results) > 0) {
  cat("显著差异的基因:\n")
  for(i in 1:nrow(significant_results)) {
    gene <- significant_results[i, "Gene"]
    cell_type <- significant_results[i, "Cell_Type"]
    fc <- round(significant_results[i, "Fold_Change"], 2)
    direction <- ifelse(fc > 1, "上调", "下调")
    cat(sprintf("- %s在%s细胞中%s (FC=%.2f, FDR=%.2e)\n", 
                gene, cell_type, direction, fc, significant_results[i, "FDR"]))
  }
} else {
  cat("未发现显著差异的基因\n")
}

