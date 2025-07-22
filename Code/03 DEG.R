library(Seurat)
library(dior)
library(reticulate)
library(SingleCellExperiment)
library(ggplot2)
###
scRNA <- read_h5(file = "scRNA.h5",
                 assay.name = 'RNA', 
                 target.object = 'seurat')
###创建分组
sample_names <- scRNA$orig.ident
# 创建分组信息向量
group <- ifelse(sample_names %in% c("M_P01", "M_L01", "M_P03"), "R",
                ifelse(sample_names %in% c("M_P02", "P_G04", "P_G05", "P_C07", "M_L07"), "S", NA))
# 将分组信息添加到 Seurat 对象的元数据中
scRNA$group <- group
##提取细胞进行二次聚类
table(scRNA$group)
Idents(scRNA) <- "cell_type"
Idents(scRNA) <- "leiden"
####
table(scRNA$cell_type)
scRNA@meta.data
library(future)
library(tidyverse)
options(future.globals.maxSize = 20 * 1024^3)
plan(multisession, workers = 4) #开启多核运算 (4个核)
logFCfilter=0.5
adjPvalFilter=0.05
scRNA.markers <- FindAllMarkers(object = scRNA,
                                only.pos = FALSE,
                                min.pct = 0.25,
                                logfc.threshold = logFCfilter)
sig.markers=scRNA.markers[(abs(as.numeric(as.vector(scRNA.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(scRNA.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="markers_by_leiden.xls",sep="\t",row.names=F,quote=F)
top10 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
plan('sequential') #终止多核运算
saveRDS(scRNA,file = "scRNA_sub.rds")
###
library(Seurat)
library(dplyr)
library(openxlsx)
library(ggplot2)
Idents(scRNA) <- "cell_type"
# 3. 对每个细胞亚群进行差异基因分析
perform_deg_analysis <- function(seurat_obj, cluster_id, group_col = "group", 
                                 control_group = "Control", treatment_group = "Treatment") {
  
  # 提取指定细胞亚群
  cluster_cells <- subset(seurat_obj, idents = cluster_id)
  
  # 设置分组信息
  Idents(cluster_cells) <- cluster_cells[[group_col]]
  
  # 进行差异基因分析
  deg_results <- FindMarkers(
    cluster_cells,
    ident.1 = treatment_group,
    ident.2 = control_group,
    min.pct = 0.1,          # 至少在10%的细胞中表达
    logfc.threshold = 0.25, # log2FC阈值
    test.use = "wilcox",    # 使用Wilcoxon秩和检验
    verbose = FALSE
  )
  
  # 添加基因名称列
  deg_results$gene <- rownames(deg_results)
  
  # 添加细胞亚群信息
  deg_results$cluster <- paste0("Cluster_", cluster_id)
  
  # 添加显著性标记
  deg_results$significance <- ifelse(deg_results$p_val_adj < 0.05, 
                                     ifelse(deg_results$avg_log2FC > 0, "Up", "Down"), 
                                     "NS")
  
  # 重新排列列顺序
  deg_results <- deg_results[, c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", 
                                 "p_val", "p_val_adj", "significance")]
  
  return(deg_results)
}

# 4. 批量分析所有细胞亚群
analyze_all_clusters <- function(seurat_obj, group_col = "group", 
                                 control_group = "Control", treatment_group = "Treatment") {
  
  # 获取所有细胞亚群
  all_clusters <- levels(Idents(seurat_obj))
  
  # 存储所有结果
  all_deg_results <- list()
  
  # 对每个细胞亚群进行分析
  for (cluster in all_clusters) {
    cat("正在分析细胞亚群:", cluster, "\n")
    
    # 检查该亚群是否有足够的细胞
    cluster_cells <- subset(seurat_obj, idents = cluster)
    group_counts <- table(cluster_cells[[group_col]])
    
    if (all(group_counts >= 3)) {  # 每组至少3个细胞
      deg_result <- perform_deg_analysis(seurat_obj, cluster, group_col, 
                                         control_group, treatment_group)
      all_deg_results[[paste0("Cluster_", cluster)]] <- deg_result
    } else {
      cat("警告: 细胞亚群", cluster, "的细胞数量不足，跳过分析\n")
    }
  }
  
  return(all_deg_results)
}
# 5. 导出结果到Excel
export_to_excel <- function(deg_results_list, filename = "DEG_results.xlsx") {
  
  # 创建工作簿
  wb <- createWorkbook()
  
  # 为每个细胞亚群创建工作表
  for (cluster_name in names(deg_results_list)) {
    
    # 添加工作表
    addWorksheet(wb, cluster_name)
    
    # 写入数据
    writeData(wb, cluster_name, deg_results_list[[cluster_name]])
    
    # 设置表头样式
    headerStyle <- createStyle(
      fontSize = 12,
      fontColour = "white",
      halign = "center",
      fgFill = "#4F81BD",
      border = "TopBottom",
      borderStyle = "medium"
    )
    
    # 应用表头样式
    addStyle(wb, cluster_name, headerStyle, rows = 1, cols = 1:ncol(deg_results_list[[cluster_name]]))
    
    # 设置列宽
    setColWidths(wb, cluster_name, cols = 1:ncol(deg_results_list[[cluster_name]]), 
                 widths = c(15, 12, 12, 10, 10, 12, 12, 12))
  }
  
  # 创建汇总表
  summary_data <- data.frame()
  for (cluster_name in names(deg_results_list)) {
    cluster_data <- deg_results_list[[cluster_name]]
    summary_row <- data.frame(
      Cluster = cluster_name,
      Total_genes = nrow(cluster_data),
      Upregulated = sum(cluster_data$significance == "Up"),
      Downregulated = sum(cluster_data$significance == "Down"),
      Non_significant = sum(cluster_data$significance == "NS")
    )
    summary_data <- rbind(summary_data, summary_row)
  }
  
  # 添加汇总表
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", summary_data)
  addStyle(wb, "Summary", headerStyle, rows = 1, cols = 1:ncol(summary_data))
  setColWidths(wb, "Summary", cols = 1:ncol(summary_data), widths = 15)
  
  # 保存文件
  saveWorkbook(wb, filename, overwrite = TRUE)
  cat("结果已保存到:", filename, "\n")
}


# 示例1：基本用法
deg_results <- analyze_all_clusters(
  seurat_obj = scRNA,
  group_col = "group",      # 你的分组列名
  control_group = "R",    # 对照组名称
  treatment_group = "S" # 实验组名称
)

export_to_excel(deg_results, filename = "DEG_results.xlsx")

