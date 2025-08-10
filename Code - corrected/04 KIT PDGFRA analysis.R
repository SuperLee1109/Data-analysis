# ====================================================
# GIST TUMOR CELL SUBPOPULATION DISTRIBUTION ANALYSIS
# Comparing PDGFRA vs KIT Mutation Profiles

# --------------------------------
# 0. ENVIRONMENT CONFIGURATION
# --------------------------------
# Load required packages
required_packages <- c("Seurat", "ggplot2", "dplyr", "patchwork", "reshape2", 
                      "RColorBrewer", "ggpubr", "rstatix", "ComplexHeatmap", "circlize")

invisible(suppressPackageStartupMessages({
  lapply(required_packages, library, character.only = TRUE)
}))

# Set plotting theme
theme_set(theme_minimal(base_size = 12) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  legend.position = "right"))

# Create output directory
output_dir <- "GIST_mutation_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------------------
# 1. ANALYSIS PARAMETERS
# --------------------------------
params <- list(
  # Metadata column names
  sample_col = "orig.ident",       # Column containing sample IDs
  cell_type_col = "cell_type",     # Column containing cell type annotations
  
  # Mutation groups
  pdgfra_samples = c("P_G06"),     # PDGFRA mutant samples
  kit_samples = c("M_P01", "M_L01", "M_P02", "M_P03", 
                  "P_G04", "P_G05", "P_C07", "M_L07"), # KIT mutant samples
  
  # Visualization parameters
  pdgfra_color = "#E64B35",        # PDGFRA mutant color (red)
  kit_color = "#4DBBD5",           # KIT mutant color (blue)
  cell_type_palette = "Set3",      # ColorBrewer palette for cell types
  
  # Output settings
  figure_width = 10,               # Figure width in inches
  figure_height = 8,               # Figure height in inches
  figure_dpi = 300                 # Figure resolution
)

# --------------------------------
# 2. DATA VALIDATION & PREPARATION
# --------------------------------
cat("\n>> Validating data structure...\n")

# Check data dimensions
cat("Total cells:", ncol(scRNA), "\n")
cat("Total genes:", nrow(scRNA), "\n")

# Validate metadata columns
if(!params$sample_col %in% colnames(scRNA@meta.data)) {
  stop(paste("Sample ID column not found. Available columns:",
             paste(colnames(scRNA@meta.data), collapse = ", ")))
}

if(!params$cell_type_col %in% colnames(scRNA@meta.data)) {
  stop(paste("Cell type column not found. Available columns:",
             paste(colnames(scRNA@meta.data), collapse = ", ")))
}

# Get sample list
all_samples <- unique(scRNA@meta.data[[params$sample_col]])
cat("\nSamples in dataset:", paste(all_samples, collapse = ", "), "\n")

# Validate mutation groups
pdgfra_found <- intersect(params$pdgfra_samples, all_samples)
kit_found <- intersect(params$kit_samples, all_samples)

if(length(pdgfra_found) == 0) stop("No PDGFRA mutant samples found in dataset")
if(length(kit_found) == 0) stop("No KIT mutant samples found in dataset")

# Create mutation type annotation
scRNA$mutation_type <- case_when(
  scRNA@meta.data[[params$sample_col]] %in% pdgfra_found ~ "PDGFRA_mutant",
  scRNA@meta.data[[params$sample_col]] %in% kit_found ~ "KIT_mutant",
  TRUE ~ "Other"
)

# Filter for mutation analysis
analysis_data <- subset(scRNA, mutation_type %in% c("PDGFRA_mutant", "KIT_mutant"))

# Report analysis dataset
cat("\nAnalysis dataset summary:\n")
cat("Cells retained:", ncol(analysis_data), "\n")
print(table(analysis_data$mutation_type))

# --------------------------------
# 3. CELL TYPE PROPORTION ANALYSIS
# --------------------------------
cat("\n>> Analyzing cell type proportions...\n")

# Calculate cell type proportions
cell_counts <- table(analysis_data@meta.data[[params$cell_type_col]], 
                     analysis_data$mutation_type)
cell_proportions <- prop.table(cell_counts, margin = 2) * 100

# Convert to data frames
cell_counts_df <- as.data.frame(cell_counts) %>%
  rename(Cell_Type = Var1, Mutation_Type = Var2, Count = Freq)

cell_prop_df <- as.data.frame(cell_proportions) %>%
  rename(Cell_Type = Var1, Mutation_Type = Var2, Proportion = Freq) %>%
  mutate(Proportion = round(Proportion, 2))

# Calculate comparative metrics
comparison_df <- cell_prop_df %>%
  pivot_wider(names_from = Mutation_Type, values_from = Proportion) %>%
  mutate(
    Difference = PDGFRA_mutant - KIT_mutant,
    Fold_Change = ifelse(KIT_mutant == 0, Inf, PDGFRA_mutant / KIT_mutant)
  ) %>%
  arrange(desc(abs(Difference)))

# --------------------------------
# 4. VISUALIZATION
# --------------------------------
cat("\n>> Generating visualizations...\n")

## 4.1 Stacked barplot of cell type distribution
stacked_plot <- ggplot(cell_prop_df, 
                       aes(x = Mutation_Type, y = Proportion, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = params$cell_type_palette) +
  labs(title = "Cell Type Distribution by Mutation Status",
       x = "Mutation Type", 
       y = "Proportion (%)",
       fill = "Cell Type") +
  theme(legend.key.size = unit(0.5, "cm"))

## 4.2 Grouped barplot of cell type proportions
grouped_plot <- ggplot(cell_prop_df, 
                       aes(x = Cell_Type, y = Proportion, fill = Mutation_Type)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  scale_fill_manual(values = c("PDGFRA_mutant" = params$pdgfra_color, 
                               "KIT_mutant" = params$kit_color)) +
  labs(title = "Cell Type Proportions by Mutation Status",
       x = "Cell Type", 
       y = "Proportion (%)",
       fill = "Mutation Type") +
  geom_text(aes(label = sprintf("%.1f%%", Proportion)), 
            position = position_dodge(0.9), vjust = -0.5, size = 3) +
  ylim(0, max(cell_prop_df$Proportion) * 1.2)

## 4.3 Heatmap of cell type proportions
heatmap_data <- dcast(cell_prop_df, Cell_Type ~ Mutation_Type, value.var = "Proportion")
rownames(heatmap_data) <- heatmap_data$Cell_Type
heatmap_data$Cell_Type <- NULL

heatmap_plot <- Heatmap(
  as.matrix(heatmap_data),
  name = "Proportion (%)",
  col = colorRamp2(c(0, max(heatmap_data)/2, max(heatmap_data)), 
                   c("white", "lightblue", "darkblue")),
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", heatmap_data[i, j]), x, y, gp = gpar(fontsize = 10))
  },
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = "Cell Type Proportions by Mutation Status",
  row_title = "Cell Type",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10)
)

## 4.4 UMAP visualizations
# UMAP by mutation type
umap_mutation <- DimPlot(analysis_data, group.by = "mutation_type",
                         cols = c("PDGFRA_mutant" = params$pdgfra_color, 
                                  "KIT_mutant" = params$kit_color)) +
  ggtitle("Mutation Status Distribution") +
  theme(legend.position = "right")

# UMAP split by cell type
umap_celltype <- DimPlot(analysis_data, group.by = "mutation_type",
                         split.by = params$cell_type_col,
                         cols = c("PDGFRA_mutant" = params$pdgfra_color, 
                                  "KIT_mutant" = params$kit_color)) +
  ggtitle("Mutation Status by Cell Type") +
  theme(legend.position = "bottom")

# --------------------------------
# 5. SAMPLE-LEVEL ANALYSIS
# --------------------------------
cat("\n>> Performing sample-level analysis...\n")

# Calculate sample-level proportions
sample_prop_df <- analysis_data@meta.data %>%
  group_by(.data[[params$sample_col]], .data[[params$cell_type_col]], mutation_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(.data[[params$sample_col]]) %>%
  mutate(total = sum(count),
         proportion = count / total * 100) %>%
  rename(Sample = 1, Cell_Type = 2)

# Generate boxplot (if sufficient samples)
if(length(unique(sample_prop_df$Sample)) > 3) {
  box_plot <- ggplot(sample_prop_df, 
                     aes(x = Cell_Type, y = proportion, fill = mutation_type)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 2) +
    scale_fill_manual(values = c("PDGFRA_mutant" = params$pdgfra_color, 
                                 "KIT_mutant" = params$kit_color)) +
    labs(title = "Cell Type Proportions by Sample",
         x = "Cell Type", 
         y = "Proportion (%)",
         fill = "Mutation Type") +
    stat_compare_means(aes(group = mutation_type), 
                       method = "wilcox", label = "p.signif", hide.ns = TRUE)
}

# --------------------------------
# 6. RESULTS EXPORT
# --------------------------------
cat("\n>> Saving results...\n")

# Save statistical results
write.csv(comparison_df, file.path(output_dir, "cell_type_comparison.csv"), row.names = FALSE)
write.csv(cell_prop_df, file.path(output_dir, "cell_proportions.csv"), row.names = FALSE)
write.csv(sample_prop_df, file.path(output_dir, "sample_level_proportions.csv"), row.names = FALSE)

# Save visualizations
ggsave(file.path(output_dir, "stacked_distribution.pdf"), stacked_plot,
       width = params$figure_width, height = params$figure_height, dpi = params$figure_dpi)

ggsave(file.path(output_dir, "grouped_comparison.pdf"), grouped_plot,
       width = params$figure_width + 2, height = params$figure_height, dpi = params$figure_dpi)

ggsave(file.path(output_dir, "umap_mutation_type.pdf"), umap_mutation,
       width = params$figure_width, height = params$figure_height, dpi = params$figure_dpi)

pdf(file.path(output_dir, "cell_proportion_heatmap.pdf"), 
    width = 6, height = 8)
draw(heatmap_plot)
dev.off()

if(exists("box_plot")) {
  ggsave(file.path(output_dir, "sample_boxplot.pdf"), box_plot,
         width = params$figure_width + 2, height = params$figure_height, dpi = params$figure_dpi)
}

# --------------------------------
# 7. ANALYSIS SUMMARY
# --------------------------------
cat("\n===== ANALYSIS SUMMARY =====\n")
cat("PDGFRA mutant samples:", length(pdgfra_found), "(", paste(pdgfra_found, collapse = ", "), ")\n")
cat("KIT mutant samples:", length(kit_found), "(", paste(kit_found, collapse = ", "), ")\n")
cat("Total cells analyzed:", ncol(analysis_data), "\n")
cat("Cell types identified:", nlevels(factor(analysis_data@meta.data[[params$cell_type_col]])), "\n\n")

# Identify key findings
top_diff <- comparison_df %>% slice_max(abs(Difference), n = 1)
pdgfra_dominant <- comparison_df %>% slice_max(PDGFRA_mutant, n = 1)
kit_dominant <- comparison_df %>% slice_max(KIT_mutant, n = 1)

cat("Key Findings:\n")
cat("1. Most differentially distributed cell type:", top_diff$Cell_Type, 
    sprintf("(Δ = %.1f%%)\n", top_diff$Difference))
cat("2. Most abundant in PDGFRA mutants:", pdgfra_dominant$Cell_Type, 
    sprintf("(%.1f%%)\n", pdgfra_dominant$PDGFRA_mutant))
cat("3. Most abundant in KIT mutants:", kit_dominant$Cell_Type, 
    sprintf("(%.1f%%)\n", kit_dominant$KIT_mutant))

cat("\nResults saved to:", output_dir, "\n")
cat("Analysis completed successfully.\n")
