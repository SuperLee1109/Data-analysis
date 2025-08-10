# --------------------------------
# 0. ENVIRONMENT CONFIGURATION
# --------------------------------
# Set seed for reproducibility
set.seed(2024)

# Load required packages
required_packages <- c("Seurat", "dior", "reticulate", "SingleCellExperiment", 
                       "ggplot2", "harmony", "clustree", "patchwork", "future", "tidyverse")

invisible(suppressPackageStartupMessages({
  lapply(required_packages, library, character.only = TRUE)
}))

# Record computational environment
cat("Analysis environment:\n")
cat("R version:", R.version$version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("Packages loaded:\n")
print(sessionInfo()$otherPkgs)

# Create output directories
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# --------------------------------
# 1. ANALYSIS PARAMETERS
# --------------------------------
params <- list(
  # Input/output files
  input_h5 = "scRNA.h5",
  output_h5 = "scRNA_sub.h5",
  output_rds = "scRNA_sub.rds",
  
  # Sample grouping definition
  resistant_groups = c("M_P01", "M_L01", "M_P03"),
  sensitive_groups = c("M_P02", "P_G04", "P_G05", "P_C07", "M_L07"),
  
  # Processing parameters
  target_cell_type = "tumor",
  normalization_method = "LogNormalize",
  scale_factor = 10000,
  n_variable_features = 2000,
  
  # Dimensionality reduction
  max_pcs = 50,
  elbow_cutoff_method = "co1",  # Options: "co1", "co2", "manual"
  manual_pc_num = 13,           # Used if method = "manual"
  
  # Clustering
  resolution_test = seq(0.1, 1.2, by = 0.2),
  selected_resolution = 0.3,
  
  # Differential expression
  logFC_threshold = 0.5,
  adj_pval_threshold = 0.05,
  min_pct = 0.25,
  top_n_markers = 10
)

# Configure parallel processing
options(future.globals.maxSize = 20 * 1024^3)
plan(multisession, workers = 4)
cat("Parallel processing enabled with 4 workers\n")

# --------------------------------
# 2. DATA LOADING & PREPROCESSING
# --------------------------------
cat("\n>> Loading integrated dataset from:", params$input_h5, "\n")

# Load Seurat object from H5 file
scRNA <- dior::read_h5(
  file = params$input_h5,
  assay.name = 'RNA',
  target.object = 'seurat'
)

# Add clinical response groups
scRNA$group <- ifelse(
  scRNA$orig.ident %in% params$responder_groups, "R",
  ifelse(scRNA$orig.ident %in% params$stable_groups, "S", NA)
)

# Report group distribution
cat("\nClinical group distribution:\n")
print(table(scRNA$group))

# --------------------------------
# 3. TUMOR CELL SUBSETTING
# --------------------------------
cat("\n>> Subsetting tumor cells...\n")

# Create tumor subset
Idents(scRNA) <- "cell_type"
scRNA_tumor <- subset(scRNA, cell_type == params$target_cell_type)

# Visualize subset
tumor_umap <- DimPlot(scRNA_tumor, group.by = "cell_type", 
                      reduction = "umap", label = TRUE) +
  ggtitle("Tumor Cell Subset") +
  theme_minimal()

ggsave("results/figures/Fig1_tumor_subset.pdf", tumor_umap,
       width = 8, height = 6, dpi = 300)

# --------------------------------
# 4. NORMALIZATION & FEATURE SELECTION
# --------------------------------
cat("\n>> Processing tumor subset...\n")

scRNA_tumor <- scRNA_tumor %>%
  NormalizeData(
    normalization.method = params$normalization_method,
    scale.factor = params$scale_factor
  ) %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = params$n_variable_features
  ) %>%
  ScaleData(features = VariableFeatures(.))

# Visualize variable features
var_feat_plot <- VariableFeaturePlot(scRNA_tumor) +
  theme_minimal() +
  ggtitle(paste("Top", params$n_variable_features, "Variable Features"))

ggsave("results/figures/Fig2_variable_features.pdf", var_feat_plot,
       width = 8, height = 6, dpi = 300)

# --------------------------------
# 5. DIMENSIONALITY REDUCTION
# --------------------------------
cat("\n>> Performing PCA...\n")

scRNA_tumor <- RunPCA(
  scRNA_tumor,
  npcs = params$max_pcs,
  features = VariableFeatures(scRNA_tumor)
)

# Automated PC selection
pct <- scRNA_tumor[["pca"]]@stdev / sum(scRNA_tumor[["pca"]]@stdev) * 100
cumu <- cumsum(pct)

if(params$elbow_cutoff_method == "co1") {
  co1 <- which(cumu > 90 & pct < 5)[1]
  pc.num <- 1:co1
} else if(params$elbow_cutoff_method == "co2") {
  co2 <- sort(which((pct[1:(length(pct)-1] - pct[2:length(pct)]) > 0.1), 
               decreasing = TRUE)[1] + 1
  pc.num <- 1:co2
} else {
  pc.num <- 1:params$manual_pc_num
}

cat("Selected PCs:", length(pc.num), "\n")

# Create elbow plot
plot_df <- data.frame(PC = 1:params$max_pcs, 
                      Variance = pct, 
                      Cumulative = cumu)

elbow_plot <- ggplot(plot_df, aes(Cumulative, Variance, label = PC)) +
  geom_point() +
  geom_text(vjust = -1, size = 3) +
  geom_vline(xintercept = 90, color = "grey50", linetype = 2) +
  geom_hline(yintercept = 5, color = "grey50", linetype = 2) +
  geom_vline(xintercept = cumu[max(pc.num)], color = "red", linetype = 2) +
  labs(title = "PCA Elbow Plot",
       subtitle = paste("Selected", length(pc.num), "PCs"),
       x = "Cumulative Variance (%)",
       y = "Individual Variance (%)") +
  theme_minimal()

ggsave("results/figures/Fig3_PCA_elbow.pdf", elbow_plot,
       width = 8, height = 6, dpi = 300)

# --------------------------------
# 6. BATCH CORRECTION & INTEGRATION
# --------------------------------
cat("\n>> Running Harmony integration...\n")

scRNA_tumor <- RunHarmony(
  scRNA_tumor,
  group.by.vars = "orig.ident",
  max.iter.harmony = 20,
  reduction = "pca",
  dims.use = pc.num
)

# --------------------------------
# 7. CLUSTERING ANALYSIS
# --------------------------------
cat("\n>> Performing clustering...\n")

# Multi-resolution clustering
scRNA_tumor <- FindNeighbors(
  scRNA_tumor,
  reduction = "harmony",
  dims = pc.num
) %>%
  FindClusters(resolution = params$resolution_test)

# Clustering tree visualization
clustree_plot <- clustree(scRNA_tumor@meta.data, prefix = "RNA_snn_res.") +
  ggtitle("Clustering Resolution Selection")

ggsave("results/figures/Fig4_clustree.pdf", clustree_plot,
       width = 10, height = 8, dpi = 300)

# Apply selected resolution
scRNA_tumor <- FindClusters(
  scRNA_tumor,
  resolution = params$selected_resolution
)

# Run UMAP
scRNA_tumor <- RunUMAP(
  scRNA_tumor,
  reduction = "harmony",
  dims = pc.num
)

# --------------------------------
# 8. VISUALIZATION & ANNOTATION
# --------------------------------
cat("\n>> Generating final visualizations...\n")

# UMAP by cluster
cluster_umap <- DimPlot(scRNA_tumor, reduction = "umap", label = TRUE) +
  ggtitle("Tumor Cell Clusters") +
  theme_minimal()

# UMAP by clinical group
group_umap <- DimPlot(scRNA_tumor, group.by = "group", reduction = "umap") +
  ggtitle("Clinical Response Groups") +
  theme_minimal()

# Combined plot
combined_plot <- cluster_umap + group_umap +
  plot_annotation(tag_levels = "A")

ggsave("results/figures/Fig5_umap_clusters.pdf", combined_plot,
       width = 14, height = 6, dpi = 300)

# Marker visualization (example: CTGF)
ctgf_plot <- FeaturePlot(scRNA_tumor, features = "CTGF", split.by = "group",
                         pt.size = 1.5, cols = c("lightgrey", "blue")) &
  theme_minimal()

ggsave("results/figures/SupFig1_CTGF_expression.pdf", ctgf_plot,
       width = 12, height = 6, dpi = 300)

# --------------------------------
# 9. DIFFERENTIAL EXPRESSION ANALYSIS
# --------------------------------
cat("\n>> Identifying cluster markers...\n")

scRNA_markers <- FindAllMarkers(
  object = scRNA_tumor,
  only.pos = FALSE,
  min.pct = params$min_pct,
  logfc.threshold = params$logFC_threshold
)

# Filter significant markers
sig_markers <- scRNA_markers %>%
  filter(abs(avg_log2FC) > params$logFC_threshold & p_val_adj < params$adj_pval_threshold)

# Save full marker table
write.csv(sig_markers, "results/tables/all_significant_markers.csv", row.names = FALSE)

# Save top markers per cluster
top_markers <- sig_markers %>%
  group_by(cluster) %>%
  top_n(n = params$top_n_markers, wt = abs(avg_log2FC))

write.csv(top_markers, "results/tables/top_markers_per_cluster.csv", row.names = FALSE)

# --------------------------------
# 10. DATA EXPORT & SESSION INFO
# --------------------------------
cat("\n>> Saving results...\n")

# Save processed data
dior::write_h5(scRNA_tumor, file = params$output_h5, object.type = "seurat")
saveRDS(scRNA_tumor, file = params$output_rds)

# Save session information
writeLines(capture.output(sessionInfo()), "results/session_info.txt")

# Disable parallel processing
plan(sequential)
cat("\nAnalysis completed successfully. Outputs saved to results/ directory\n")
