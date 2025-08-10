# ========================================================
# CCN2 RECEPTOR GENE EXPRESSION ANALYSIS IN GIST
# Comparing Sensitive vs Resistant Groups in Tumor and Stromal Cells

# --------------------------------
# 0. ENVIRONMENT CONFIGURATION
# --------------------------------
# Load required packages
required_packages <- c("Seurat", "ggplot2", "dplyr", "patchwork", "viridis", 
                      "ggpubr", "rstatix", "ComplexHeatmap", "circlize")

invisible(suppressPackageStartupMessages({
  lapply(required_packages, library, character.only = TRUE)
}))

# Set plotting theme
theme_set(theme_minimal(base_size = 12) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5),
                  legend.position = "right"))

# Create output directory
output_dir <- "CCN2_sensitivity_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------------------
# 1. ANALYSIS PARAMETERS
# --------------------------------
params <- list(
  # Metadata columns
  response_col = "group",        # Column with sensitivity information
  cell_type_col = "cell_type",   # Column with cell type annotations
  
  # Group definitions
  sensitive_labels = "S",        # Labels for sensitive group
  resistant_labels = "R",        # Labels for resistant group
  
  # Cell type definitions
  tumor_cells = "tumor",         # Tumor cell type identifier
  stromal_cells = "fibroblast",  # Stromal cell type identifier
  
  # Genes of interest
  ccn2_receptors = c("ITGAV", "ITGB3", "ITGA5", "LRP1", "DDR2"),
  
  # Visualization parameters
  sensitive_color = "#4CAF50",   # Green for sensitive
  resistant_color = "#F44336",   # Red for resistant
  tumor_color = "#1F77B4",       # Blue for tumor
  stromal_color = "#FF7F0E",     # Orange for stromal
  
  # Output settings
  figure_width = 15,             # Figure width in inches
  figure_height = 10,            # Figure height in inches
  figure_dpi = 300               # Figure resolution
)

# --------------------------------
# 2. DATA VALIDATION & PREPARATION
# --------------------------------
cat("\n>> Validating data structure...\n")

# Validate metadata columns
if(!params$response_col %in% colnames(scRNA@meta.data)) {
  stop(paste("Response column not found. Available columns:",
             paste(colnames(scRNA@meta.data), collapse = ", ")))
}

if(!params$cell_type_col %in% colnames(scRNA@meta.data)) {
  stop(paste("Cell type column not found. Available columns:",
             paste(colnames(scRNA@meta.data), collapse = ", ")))
}

# Create standardized response groups
scRNA$drug_response <- case_when(
  scRNA@meta.data[[params$response_col]] %in% params$sensitive_labels ~ "Sensitive",
  scRNA@meta.data[[params$response_col]] %in% params$resistant_labels ~ "Resistant",
  TRUE ~ "Other"
)

# Create major cell type groups
scRNA$cell_major_type <- case_when(
  scRNA@meta.data[[params$cell_type_col]] == params$tumor_cells ~ "Tumor",
  scRNA@meta.data[[params$cell_type_col]] == params$stromal_cells ~ "Stromal",
  TRUE ~ "Other"
)

# Filter for analysis
analysis_data <- subset(scRNA, 
                        drug_response %in% c("Sensitive", "Resistant") & 
                          cell_major_type %in% c("Tumor", "Stromal"))

# Report analysis dataset
cat("\nAnalysis dataset summary:\n")
cat("Total cells:", ncol(analysis_data), "\n")
print(table(analysis_data$drug_response, analysis_data$cell_major_type))

# Validate genes of interest
genes_present <- intersect(params$ccn2_receptors, rownames(analysis_data))
genes_missing <- setdiff(params$ccn2_receptors, rownames(analysis_data))

if(length(genes_present) == 0) stop("No CCN2 receptor genes found in dataset")
cat("\nCCN2 receptors present:", paste(genes_present, collapse = ", "), "\n")
if(length(genes_missing) > 0) {
  cat("Missing genes:", paste(genes_missing, collapse = ", "), "\n")
}

# --------------------------------
# 3. DIFFERENTIAL EXPRESSION ANALYSIS
# --------------------------------
cat("\n>> Performing differential expression analysis...\n")

#' Perform DE analysis for a specific cell type
#' 
#' @param seurat_obj Seurat object
#' @param cell_type Cell type to analyze ("Tumor" or "Stromal")
#' @param min_cells Minimum cells per group for analysis
#' @return List containing plots and statistics
analyze_cell_type <- function(seurat_obj, cell_type, min_cells = 20) {
  # Subset data
  cell_data <- subset(seurat_obj, cell_major_type == cell_type)
  
  # Check cell counts
  cell_counts <- table(cell_data$drug_response)
  if(any(cell_counts < min_cells)) {
    cat("Skipping", cell_type, "- insufficient cells (min:", min_cells, "required)\n")
    return(list(plots = NULL, stats = NULL))
  }
  
  # Initialize storage
  plots <- list()
  stats <- data.frame()
  
  # Analyze each gene
  for(gene in genes_present) {
    # Extract expression data
    expr_data <- FetchData(cell_data, vars = gene)
    expr_data$Response <- cell_data$drug_response
    
    # Calculate summary statistics
    sum_stats <- expr_data %>%
      group_by(Response) %>%
      summarise(
        Mean = mean(.data[[gene]]),
        SD = sd(.data[[gene]]),
        Cells = n(),
        .groups = "drop"
      )
    
    # Wilcoxon test
    if(sum_stats$Cells[sum_stats$Response == "Sensitive"] >= min_cells &&
       sum_stats$Cells[sum_stats$Response == "Resistant"] >= min_cells) {
      test_res <- wilcox.test(as.formula(paste(gene, "~ Response")), data = expr_data)
      p_value <- test_res$p.value
    } else {
      p_value <- NA
    }
    
    # Fold change
    sensitive_mean <- sum_stats$Mean[sum_stats$Response == "Sensitive"]
    resistant_mean <- sum_stats$Mean[sum_stats$Response == "Resistant"]
    fold_change <- resistant_mean / (sensitive_mean + 1e-10)
    
    # Store statistics
    stats <- rbind(stats, data.frame(
      Gene = gene,
      Cell_Type = cell_type,
      Sensitive_Mean = sensitive_mean,
      Resistant_Mean = resistant_mean,
      Fold_Change = fold_change,
      P_Value = p_value
    ))
    
    # Create violin plot
    p <- ggplot(expr_data, aes(x = Response, y = .data[[gene]], fill = Response)) +
      geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
      geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.shape = NA) +
      geom_jitter(width = 0.1, size = 0.5, alpha = 0.3) +
      scale_fill_manual(values = c("Sensitive" = params$sensitive_color, 
                                   "Resistant" = params$resistant_color)) +
      labs(title = paste(gene, "in", cell_type, "Cells"),
           x = "Drug Response", y = "Expression Level") +
      theme(legend.position = "none")
    
    # Add significance annotation
    if(!is.na(p_value) && p_value < 0.05) {
      sig_label <- case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      )
      
      y_max <- max(expr_data[[gene]]) * 1.1
      p <- p + 
        geom_segment(x = 1, xend = 2, y = y_max * 0.95, yend = y_max * 0.95) +
        annotate("text", x = 1.5, y = y_max, 
                 label = sig_label, size = 5, vjust = 0.5) +
        annotate("text", x = 1.5, y = y_max * 0.9, 
                 label = sprintf("p = %.2e", p_value), size = 3)
    }
    
    plots[[gene]] <- p
  }
  
  return(list(plots = plots, stats = stats))
}

# Run analysis for tumor and stromal cells
tumor_results <- analyze_cell_type(analysis_data, "Tumor")
stromal_results <- analyze_cell_type(analysis_data, "Stromal")

# Combine statistics
all_stats <- rbind(tumor_results$stats, stromal_results$stats) %>%
  mutate(FDR = p.adjust(P_Value, method = "fdr"),
         Significance = case_when(
           FDR < 0.001 ~ "***",
           FDR < 0.01 ~ "**",
           FDR < 0.05 ~ "*",
           is.na(FDR) ~ "NA",
           TRUE ~ "ns"
         ))

# --------------------------------
# 4. VISUALIZATION
# --------------------------------
cat("\n>> Generating visualizations...\n")

## 4.1 Combined expression plots
if(!is.null(tumor_results$plots)) {
  tumor_combined <- wrap_plots(tumor_results$plots, ncol = 3) +
    plot_annotation(title = "CCN2 Receptor Expression in Tumor Cells",
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
}

if(!is.null(stromal_results$plots)) {
  stromal_combined <- wrap_plots(stromal_results$plots, ncol = 3) +
    plot_annotation(title = "CCN2 Receptor Expression in Stromal Cells",
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
}

## 4.2 Heatmap of fold changes
# Prepare heatmap data
heatmap_data <- all_stats %>%
  select(Gene, Cell_Type, Fold_Change) %>%
  tidyr::pivot_wider(names_from = Cell_Type, values_from = Fold_Change) %>%
  column_to_rownames("Gene")

# Prepare significance data
sig_data <- all_stats %>%
  select(Gene, Cell_Type, Significance) %>%
  tidyr::pivot_wider(names_from = Cell_Type, values_from = Significance) %>%
  column_to_rownames("Gene")

# Create heatmap
heatmap_plot <- Heatmap(
  as.matrix(log2(heatmap_data + 1e-5)),  # Log-transform for better visualization
  name = "log2(FC)",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  rect_gp = gpar(col = "gray80", lwd = 0.5),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(!is.na(sig_data[i, j]) && sig_data[i, j] != "ns") {
      grid.text(sig_data[i, j], x, y, gp = gpar(fontsize = 10, fontface = "bold"))
    }
  },
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "top",
  column_title = "CCN2 Receptors: Resistant vs Sensitive (log2 Fold Change)",
  heatmap_legend_param = list(title_position = "topcenter")
)

## 4.3 UMAP visualization
umap_plot <- DimPlot(analysis_data, 
                     group.by = "drug_response",
                     split.by = "cell_major_type",
                     cols = c("Sensitive" = params$sensitive_color, 
                              "Resistant" = params$resistant_color)) +
  ggtitle("Drug Response in Tumor vs Stromal Cells") +
  theme(plot.title = element_text(hjust = 0.5))

# --------------------------------
# 5. RESULTS EXPORT
# --------------------------------
cat("\n>> Saving results...\n")

# Save statistics
write.csv(all_stats, file.path(output_dir, "CCN2_receptors_sensitivity_stats.csv"), 
          row.names = FALSE)

# Save visualizations
if(exists("tumor_combined")) {
  ggsave(file.path(output_dir, "Tumor_CCN2_expression.pdf"), tumor_combined,
         width = params$figure_width, height = params$figure_height, dpi = params$figure_dpi)
}

if(exists("stromal_combined")) {
  ggsave(file.path(output_dir, "Stromal_CCN2_expression.pdf"), stromal_combined,
         width = params$figure_width, height = params$figure_height, dpi = params$figure_dpi)
}

ggsave(file.path(output_dir, "Response_UMAP.pdf"), umap_plot,
       width = 12, height = 6, dpi = params$figure_dpi)

pdf(file.path(output_dir, "Fold_change_heatmap.pdf"), width = 6, height = 5)
draw(heatmap_plot)
dev.off()

# --------------------------------
# 6. ANALYSIS SUMMARY
# --------------------------------
cat("\n===== ANALYSIS SUMMARY =====\n")
cat("Analysis performed on", ncol(analysis_data), "cells\n")
cat("Tumor cells (Sensitive/Resistant):", 
    sum(analysis_data$cell_major_type == "Tumor" & analysis_data$drug_response == "Sensitive"), "/",
    sum(analysis_data$cell_major_type == "Tumor" & analysis_data$drug_response == "Resistant"), "\n")
cat("Stromal cells (Sensitive/Resistant):", 
    sum(analysis_data$cell_major_type == "Stromal" & analysis_data$drug_response == "Sensitive"), "/",
    sum(analysis_data$cell_major_type == "Stromal" & analysis_data$drug_response == "Resistant"), "\n\n")

# Identify significant results
sig_results <- all_stats %>% filter(FDR < 0.05)

if(nrow(sig_results) > 0) {
  cat("Significant differential expression results (FDR < 0.05):\n")
  for(i in 1:nrow(sig_results)) {
    res <- sig_results[i, ]
    direction <- ifelse(res$Fold_Change > 1, "upregulated", "downregulated")
    cat(sprintf("- %s in %s cells: %s (FC=%.2f, FDR=%.2e)\n",
                res$Gene, res$Cell_Type, direction, res$Fold_Change, res$FDR))
  }
} else {
  cat("No significant differential expression found at FDR < 0.05\n")
}

# Identify most significant result
if(nrow(sig_results) > 0) {
  top_hit <- sig_results %>% arrange(FDR) %>% slice(1)
  cat(sprintf("\nMost significant result: %s in %s cells (FDR=%.2e)\n",
              top_hit$Gene, top_hit$Cell_Type, top_hit$FDR))
}

cat("\nResults saved to:", output_dir, "\n")
cat("Analysis completed successfully.\n")
