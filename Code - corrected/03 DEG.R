# --------------------------------
# 0. ENVIRONMENT CONFIGURATION
# --------------------------------
set.seed(2024)  # Fix random seed for reproducibility

# Load required packages
required_packages <- c("Seurat", "dior", "reticulate", "SingleCellExperiment",
                      "ggplot2", "tidyverse", "future", "openxlsx", "dplyr")

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
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# Configure parallel processing
options(future.globals.maxSize = 20 * 1024^3)
plan(multisession, workers = 4)
cat("Parallel processing enabled with 4 workers\n")

# --------------------------------
# 1. ANALYSIS PARAMETERS
# --------------------------------
params <- list(
  # Input/output files
  input_h5 = "scRNA.h5",
  output_rds = "scRNA_sub.rds",
  
  # Group definitions
  responder_groups = c("M_P01", "M_L01", "M_P03"),
  stable_groups = c("M_P02", "P_G04", "P_G05", "P_C07", "M_L07"),
  
  # Differential expression settings
  group_column = "group",
  control_group = "R",          # Responder group
  treatment_group = "S",        # Stable group
  min_pct = 0.25,               # Minimum expression percentage
  logFC_threshold = 0.5,        # Log fold-change threshold
  adj_pval_threshold = 0.05,    # Adjusted p-value threshold
  min_cells_per_group = 10,     # Minimum cells per group for DE analysis
  test_method = "wilcox"        # Differential expression test method
)

# --------------------------------
# 2. DATA LOADING & GROUP ASSIGNMENT
# --------------------------------
cat("\n>> Loading integrated dataset from:", params$input_h5, "\n")

# Load Seurat object
scRNA <- dior::read_h5(
  file = params$input_h5,
  assay.name = 'RNA',
  target.object = 'seurat'
)

# Add clinical response groups
scRNA$group <- ifelse(
  scRNA$orig.ident %in% params$responder_groups, params$control_group,
  ifelse(scRNA$orig.ident %in% params$stable_groups, params$treatment_group, NA)
)

# Report group distribution
cat("\nClinical group distribution:\n")
print(table(scRNA$group))

# --------------------------------
# 3. PANEL-WIDE DIFFERENTIAL EXPRESSION
# --------------------------------
cat("\n>> Performing differential expression across all cell types...\n")

# Set identity to cell type
Idents(scRNA) <- "cell_type"

# Run differential expression
scRNA.markers <- FindAllMarkers(
  object = scRNA,
  only.pos = FALSE,
  min.pct = params$min_pct,
  logfc.threshold = params$logFC_threshold,
  test.use = params$test_method
)

# Filter significant markers
sig.markers <- scRNA.markers %>%
  filter(abs(avg_log2FC) > params$logFC_threshold & 
         p_val_adj < params$adj_pval_threshold)

# Save results
write.table(sig.markers, file = "results/tables/pan_celltype_markers.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Extract top markers
top10 <- sig.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = abs(avg_log2FC))

# --------------------------------
# 4. CELL-TYPE SPECIFIC DIFFERENTIAL EXPRESSION
# --------------------------------
cat("\n>> Performing cell-type specific differential expression...\n")

#' Perform differential expression analysis for a specific cell type
#' 
#' @param seurat_obj Seurat object
#' @param cluster_id Cell type/cluster identifier
#' @param group_col Metadata column containing group information
#' @param control_group Name of control group
#' @param treatment_group Name of treatment group
#' @return Dataframe with differential expression results
perform_cluster_deg <- function(seurat_obj, cluster_id, group_col, 
                                control_group, treatment_group) {
  
  # Extract cell subset
  cluster_cells <- subset(seurat_obj, idents = cluster_id)
  
  # Check cell count thresholds
  cell_counts <- table(cluster_cells[[group_col]])
  if (any(cell_counts < params$min_cells_per_group)) {
    warning(paste("Skipping cluster", cluster_id, 
                  "- insufficient cells (min:", params$min_cells_per_group, ")"))
    return(NULL)
  }
  
  # Set group identities
  Idents(cluster_cells) <- group_col
  
  # Run differential expression
  deg_results <- FindMarkers(
    object = cluster_cells,
    ident.1 = treatment_group,
    ident.2 = control_group,
    min.pct = params$min_pct,
    logfc.threshold = params$logFC_threshold,
    test.use = params$test_method,
    verbose = FALSE
  )
  
  # Add gene and cluster information
  deg_results$gene <- rownames(deg_results)
  deg_results$cluster <- cluster_id
  
  # Add significance classification
  deg_results <- deg_results %>%
    mutate(
      significance = case_when(
        p_val_adj >= 0.05 ~ "NS",
        p_val_adj < 0.05 & avg_log2FC > 0 ~ "Up",
        p_val_adj < 0.05 & avg_log2FC < 0 ~ "Down"
      )
    )
  
  # Reorder columns
  deg_results <- deg_results %>%
    select(gene, cluster, avg_log2FC, pct.1, pct.2, 
           p_val, p_val_adj, significance)
  
  return(deg_results)
}

# Batch DE analysis across all cell types
all_cluster_deg <- lapply(levels(Idents(scRNA)), function(cluster_id) {
  cat("Processing cluster:", cluster_id, "\n")
  deg_df <- perform_cluster_deg(
    seurat_obj = scRNA,
    cluster_id = cluster_id,
    group_col = params$group_column,
    control_group = params$control_group,
    treatment_group = params$treatment_group
  )
  return(deg_df)
}) %>% bind_rows()

# Filter out NULL results
all_cluster_deg <- all_cluster_deg %>% filter(!is.na(gene))

# --------------------------------
# 5. RESULT EXPORT & VISUALIZATION
# --------------------------------
cat("\n>> Exporting results...\n")

#' Export differential expression results to Excel with formatting
#' 
#' @param deg_results Differential expression results dataframe
#' @param filename Output filename
export_deg_results <- function(deg_results, filename) {
  
  # Create workbook
  wb <- createWorkbook()
  
  # Add summary sheet
  summary_data <- deg_results %>%
    group_by(cluster) %>%
    summarise(
      Total_genes = n(),
      Upregulated = sum(significance == "Up"),
      Downregulated = sum(significance == "Down"),
      Non_significant = sum(significance == "NS")
    )
  
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", summary_data)
  
  # Apply formatting to summary sheet
  headerStyle <- createStyle(
    textDecoration = "bold",
    fgFill = "#4F81BD",
    fontColour = "white",
    halign = "center"
  )
  
  addStyle(wb, "Summary", headerStyle, rows = 1, cols = 1:5)
  
  # Add sheets for each cluster
  clusters <- unique(deg_results$cluster)
  for (cluster in clusters) {
    cluster_data <- deg_results %>% filter(cluster == !!cluster)
    
    # Sort by significance and fold-change
    cluster_data <- cluster_data %>%
      arrange(significance, desc(abs(avg_log2FC)))
    
    # Create worksheet
    addWorksheet(wb, sheetName = cluster)
    writeData(wb, cluster, cluster_data)
    
    # Apply conditional formatting
    posStyle <- createStyle(fontColour = "#E41A1C")  # Red for upregulated
    negStyle <- createStyle(fontColour = "#377EB8")  # Blue for downregulated
    
    conditionalFormatting(wb, cluster, 
                          cols = 8, rows = 2:(nrow(cluster_data)+1,
                          rule = '=="Up"', style = posStyle)
    
    conditionalFormatting(wb, cluster, 
                          cols = 8, rows = 2:(nrow(cluster_data)+1,
                          rule = '=="Down"', style = negStyle)
    
    # Add header style
    addStyle(wb, cluster, headerStyle, rows = 1, cols = 1:ncol(cluster_data))
  }
  
  # Save workbook
  saveWorkbook(wb, file.path("results/tables", filename), overwrite = TRUE)
  cat("Results exported to:", file.path("results/tables", filename), "\n")
}

# Export results
export_deg_results(all_cluster_deg, "celltype_specific_deg_results.xlsx")

# Save significant genes
write.csv(all_cluster_deg, "results/tables/all_celltype_specific_deg.csv", 
          row.names = FALSE)

# --------------------------------
# 6. SESSION CLEANUP
# --------------------------------
# Save processed data
saveRDS(scRNA, file = params$output_rds)

# Save session information
writeLines(capture.output(sessionInfo()), "results/session_info.txt")

# Disable parallel processing
plan(sequential)
cat("\nDifferential expression analysis completed successfully\n")
