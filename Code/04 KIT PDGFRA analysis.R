# GIST肿瘤PDGFRA vs KIT突变细胞亚群分布比较分析

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(rstatix)

# 1. 检查数据和样本信息
# 根据用户提供的信息，对象名为 scRNA
cat("数据概况：\n")
cat("总细胞数：", ncol(scRNA), "\n")
cat("总基因数：", nrow(scRNA), "\n")

# 检查meta.data列
cat("\n可用的meta.data列：\n")
print(colnames(scRNA@meta.data))

# 根据用户提供的信息指定列名
sample_col <- "orig.ident"  # 用户的样本ID列
cell_type_col <- "cell_type"  # 请根据实际情况修改细胞类型列名

# 检查样本ID列是否存在
if(!sample_col %in% colnames(scRNA@meta.data)) {
  cat("错误：找不到样本ID列，请检查列名\n")
  cat("可选的列名：", paste(colnames(scRNA@meta.data), collapse = ", "), "\n")
  stop("请修改 sample_col 变量为正确的列名")
}

# 检查细胞类型列是否存在
if(!cell_type_col %in% colnames(scRNA@meta.data)) {
  cat("错误：找不到细胞类型列，请检查列名\n")
  stop("请修改 cell_type_col 变量为正确的列名")
}

# 2. 根据样本ID创建突变类型分组
# 根据用户实际数据更新样本列表（使用下划线格式）
pdgfra_samples <- c("P_G06")
kit_samples <- c("M_P01", "M_L01", "M_P02", "M_P03", "P_G04", "P_G05", "P_C07", "M_L07")

# 检查样本是否存在于数据中
all_samples <- unique(scRNA@meta.data[[sample_col]])
cat("\n数据中的样本：", paste(all_samples, collapse = ", "), "\n")

# 检查PDGFRA样本
pdgfra_found <- intersect(pdgfra_samples, all_samples)
pdgfra_missing <- setdiff(pdgfra_samples, all_samples)
cat("\n找到的PDGFRA样本：", paste(pdgfra_found, collapse = ", "), "\n")
if(length(pdgfra_missing) > 0) {
  cat("缺失的PDGFRA样本：", paste(pdgfra_missing, collapse = ", "), "\n")
}

# 检查KIT样本
kit_found <- intersect(kit_samples, all_samples)
kit_missing <- setdiff(kit_samples, all_samples)
cat("找到的KIT样本：", paste(kit_found, collapse = ", "), "\n")
if(length(kit_missing) > 0) {
  cat("缺失的KIT样本：", paste(kit_missing, collapse = ", "), "\n")
}

# 创建突变类型分组
scRNA$mutation_type <- ifelse(scRNA@meta.data[[sample_col]] %in% pdgfra_found, 
                              "PDGFRA_mutant",
                              ifelse(scRNA@meta.data[[sample_col]] %in% kit_found, 
                                     "KIT_mutant", "Other"))

# 显示分组结果
cat("\n突变类型分组结果：\n")
print(table(scRNA$mutation_type))

# 按样本和突变类型的分布
cat("\n各样本的突变类型：\n")
sample_mutation_table <- table(scRNA@meta.data[[sample_col]], scRNA$mutation_type)
print(sample_mutation_table)

# 3. 只保留PDGFRA和KIT突变样本进行分析
analysis_data <- subset(scRNA, mutation_type %in% c("PDGFRA_mutant", "KIT_mutant"))

cat("\n分析数据概况：\n")
cat("保留的细胞数：", ncol(analysis_data), "\n")
print(table(analysis_data$mutation_type))

# 4. 计算各细胞亚群在两组中的分布
# 每个突变类型中各细胞亚群的绝对数量
cell_counts <- table(analysis_data@meta.data[[cell_type_col]], analysis_data$mutation_type)
cat("\n各细胞亚群绝对数量：\n")
print(cell_counts)

# 每个突变类型中各细胞亚群的比例
cell_proportions <- prop.table(cell_counts, margin = 2) * 100  # 按列计算比例
cat("\n各细胞亚群比例(%)：\n")
print(round(cell_proportions, 2))

# 转换为数据框便于分析
cell_counts_df <- as.data.frame(cell_counts)
colnames(cell_counts_df) <- c("Cell_Type", "Mutation_Type", "Count")

cell_prop_df <- as.data.frame(cell_proportions)
colnames(cell_prop_df) <- c("Cell_Type", "Mutation_Type", "Proportion")

# 5. 统计检验
# 由于PDGFRA组只有一个样本，我们主要进行描述性比较
cat("\n=== 细胞亚群分布比较 ===\n")

# 计算每个细胞类型在两组间的差异
comparison_results <- data.frame()
for(cell_type in unique(cell_prop_df$Cell_Type)) {
  pdgfra_prop <- cell_prop_df[cell_prop_df$Cell_Type == cell_type & 
                                cell_prop_df$Mutation_Type == "PDGFRA_mutant", "Proportion"]
  kit_prop <- cell_prop_df[cell_prop_df$Cell_Type == cell_type & 
                             cell_prop_df$Mutation_Type == "KIT_mutant", "Proportion"]
  
  # 处理缺失值
  if(length(pdgfra_prop) == 0) pdgfra_prop <- 0
  if(length(kit_prop) == 0) kit_prop <- 0
  
  fold_change <- ifelse(kit_prop == 0, Inf, pdgfra_prop / kit_prop)
  difference <- pdgfra_prop - kit_prop
  
  comparison_results <- rbind(comparison_results, data.frame(
    Cell_Type = cell_type,
    PDGFRA_Proportion = pdgfra_prop,
    KIT_Proportion = kit_prop,
    Difference = difference,
    Fold_Change = fold_change
  ))
}

print(comparison_results)

# 6. 可视化分析

# 6.1 堆积柱状图显示细胞亚群分布
stacked_plot <- ggplot(cell_prop_df, aes(x = Mutation_Type, y = Proportion, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  labs(title = "Cell Type Distribution by Mutation Type",
       x = "Mutation Type", 
       y = "Proportion (%)",
       fill = "Cell Type") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

print(stacked_plot)

# 6.2 并排柱状图显示细胞亚群比例
grouped_plot <- ggplot(cell_prop_df, aes(x = Cell_Type, y = Proportion, fill = Mutation_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("PDGFRA_mutant" = "#FF6B6B", "KIT_mutant" = "#4ECDC4")) +
  labs(title = "Cell Type Proportions: PDGFRA vs KIT Mutant",
       x = "Cell Type", 
       y = "Proportion (%)",
       fill = "Mutation Type") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  geom_text(aes(label = paste0(round(Proportion, 1), "%")), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3)

print(grouped_plot)

# 6.3 差异热图
library(reshape2)
# 重塑数据用于热图
heatmap_data <- dcast(cell_prop_df, Cell_Type ~ Mutation_Type, value.var = "Proportion")
rownames(heatmap_data) <- heatmap_data$Cell_Type
heatmap_data$Cell_Type <- NULL

# 如果有ComplexHeatmap包，创建热图
if(require(ComplexHeatmap, quietly = TRUE)) {
  library(ComplexHeatmap)
  library(circlize)
  
  col_fun <- colorRamp2(c(0, max(heatmap_data, na.rm = TRUE)/2, max(heatmap_data, na.rm = TRUE)), 
                        c("white", "lightblue", "darkblue"))
  
  ht <- Heatmap(as.matrix(heatmap_data),
                name = "Proportion (%)",
                col = col_fun,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%.1f", heatmap_data[i, j]), x, y, 
                            gp = gpar(fontsize = 10))
                },
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                column_title = "Cell Type Proportions by Mutation Type",
                row_names_side = "left",
                column_names_side = "top")
  
  print(ht)
}

# 6.4 UMAP图显示突变类型分布
umap_mutation <- DimPlot(analysis_data, 
                         group.by = "mutation_type",
                         cols = c("PDGFRA_mutant" = "#FF6B6B", "KIT_mutant" = "#4ECDC4")) +
  ggtitle("UMAP: PDGFRA vs KIT Mutant Samples") +
  theme(plot.title = element_text(hjust = 0.5))

print(umap_mutation)

# 6.5 按细胞类型分面的UMAP图
umap_celltype <- DimPlot(analysis_data, 
                         group.by = "mutation_type",
                         split.by = cell_type_col,
                         cols = c("PDGFRA_mutant" = "#FF6B6B", "KIT_mutant" = "#4ECDC4")) +
  ggtitle("UMAP by Cell Type and Mutation Status")

print(umap_celltype)

# 7. 按样本分析每个细胞类型的分布
# 计算每个样本中各细胞类型的比例
sample_cell_prop <- analysis_data@meta.data %>%
  group_by(.data[[sample_col]], .data[[cell_type_col]], mutation_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(.data[[sample_col]]) %>%
  mutate(total = sum(count),
         proportion = count / total * 100)

colnames(sample_cell_prop)[1:2] <- c("Sample", "Cell_Type")

# 样本水平的箱线图（如果有足够样本）
if(length(kit_found) > 1) {
  box_plot <- ggplot(sample_cell_prop, aes(x = Cell_Type, y = proportion, fill = mutation_type)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.8) +
    scale_fill_manual(values = c("PDGFRA_mutant" = "#FF6B6B", "KIT_mutant" = "#4ECDC4")) +
    labs(title = "Cell Type Proportions by Sample",
         x = "Cell Type", 
         y = "Proportion (%)",
         fill = "Mutation Type") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(box_plot)
}

# 8. 保存结果
# 创建输出目录
if(!dir.exists("GIST_mutation_analysis")) {
  dir.create("GIST_mutation_analysis")
}

# 保存统计结果
write.csv(comparison_results, "GIST_mutation_analysis/cell_type_comparison.csv", row.names = FALSE)
write.csv(cell_prop_df, "GIST_mutation_analysis/cell_proportions.csv", row.names = FALSE)
write.csv(sample_cell_prop, "GIST_mutation_analysis/sample_level_proportions.csv", row.names = FALSE)

# 保存图片
ggsave("GIST_mutation_analysis/stacked_distribution.pdf", 
       stacked_plot, width = 10, height = 8, dpi = 300)

ggsave("GIST_mutation_analysis/grouped_comparison.pdf", 
       grouped_plot, width = 12, height = 8, dpi = 300)

ggsave("GIST_mutation_analysis/umap_mutation_type.pdf", 
       umap_mutation, width = 10, height = 8, dpi = 300)

if(exists("box_plot")) {
  ggsave("GIST_mutation_analysis/sample_boxplot.pdf", 
         box_plot, width = 12, height = 8, dpi = 300)
}

# 9. 总结报告
cat("\n=== 分析总结 ===\n")
cat("PDGFRA突变组样本数：", length(pdgfra_found), "\n")
cat("KIT突变组样本数：", length(kit_found), "\n")
cat("总分析细胞数：", ncol(analysis_data), "\n")
cat("细胞类型数量：", length(unique(analysis_data@meta.data[[cell_type_col]])), "\n")

cat("\n主要发现：\n")
# 找出比例差异最大的细胞类型
max_diff_idx <- which.max(abs(comparison_results$Difference))
max_diff_cell <- comparison_results[max_diff_idx, "Cell_Type"]
max_diff_value <- comparison_results[max_diff_idx, "Difference"]

cat("1. 比例差异最大的细胞类型：", max_diff_cell, 
    sprintf("(差异: %.2f%%)\n", max_diff_value))

# 找出比例最高的细胞类型
pdgfra_dominant <- comparison_results[which.max(comparison_results$PDGFRA_Proportion), "Cell_Type"]
kit_dominant <- comparison_results[which.max(comparison_results$KIT_Proportion), "Cell_Type"]

cat("2. PDGFRA组主要细胞类型：", pdgfra_dominant, "\n")
cat("3. KIT组主要细胞类型：", kit_dominant, "\n")

cat("\n结果文件已保存到 GIST_mutation_analysis/ 目录\n")
cat("分析完成！\n")