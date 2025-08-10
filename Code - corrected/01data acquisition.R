# --------------------------------
# 0. ENVIRONMENT CONFIGURATION
# --------------------------------
set.seed(2024)  # Fix random seed for reproducibility

# Load required packages with version control
required_packages <- c("Seurat", "harmony", "tidyverse", "patchwork", "dior")
invisible(lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}))

# Record computational environment
cat("Analysis environment:\n")
cat("R version:", R.version$version.string, "\n")
cat("Platform:", R.version$platform, "\n\n")
cat("Package versions:\n")
print(sapply(required_packages, function(pkg) paste(pkg, packageVersion(pkg))))

# Create standardized directory structure
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("processed_data", showWarnings = FALSE)

# --------------------------------
# 1. ANALYSIS PARAMETERS
# --------------------------------
params <- list(
  # Data paths
  raw_data_dir = "GSE254762_RAW/",
  
  # Sample metadata (aligned with GEO)
  samples = c("M_P01", "M_L01", "M_P02", "M_P03", 
              "P_G04", "P_G05", "P_G06", "P_C07", "M_L07"),
  
  # Quality control thresholds
  min_cells = 3,       # Minimum cells expressing a gene
  min_features = 200,  # Minimum features per cell
  mt_threshold = 15,   # Maximum mitochondrial percentage
  rb_threshold = 40,   # Maximum ribosomal percentage
  hb_threshold = 1,    # Maximum hemoglobin percentage
  
  # Data export formats
  export_formats = c("rds", "h5")  # Multiple format support
)

# --------------------------------
# 2. DATA LOADING & PREPROCESSING
# --------------------------------
cat("\n>> Loading raw data from:", params$raw_data_dir, "\n")

# Construct sample paths
sample_paths <- file.path(params$raw_data_dir, dir(params$raw_data_dir))

# Validate file integrity
missing_files <- sample_paths[!file.exists(sample_paths)]
if (length(missing_files) > 0) {
  stop(paste("Missing data files:", paste(missing_files, collapse = ", ")))
}

# Initialize data storage
scRNA_list <- vector("list", length(params$samples))

for(i in seq_along(sample_paths)) {
  cat("Processing sample:", params$samples[i], "\n")
  
  # Read count matrix
  counts <- Seurat::Read10X(data.dir = sample_paths[i])
  
  # Create Seurat object with QC thresholds
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = params$samples[i],
    min.cells = params$min_cells,
    min.features = params$min_features
  )
  
  # Calculate QC metrics
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  
  # Detect hemoglobin genes
  hb_genes <- c("HBA1","HBA2","HBB","HBD","HBE1",
                "HBG1","HBG2","HBM","HBQ1","HBZ")
  hb_genes <- CaseMatch(hb_genes, rownames(seurat_obj))
  if(length(hb_genes) > 0) {
    seurat_obj[["percent.HB"]] <- PercentageFeatureSet(seurat_obj, features = hb_genes)
  }
  
  # Add sample prefix to cell barcodes
  seurat_obj <- RenameCells(seurat_obj, add.cell.id = params$samples[i])
  
  scRNA_list[[i]] <- seurat_obj
}
names(scRNA_list) <- params$samples

# --------------------------------
# 3. DATA INTEGRATION & EXPORT
# --------------------------------
cat("\n>> Merging datasets...\n")
scRNA <- merge(scRNA_list[[1]], scRNA_list[2:length(scRNA_list)])

# Save merged data in RDS format
if("rds" %in% params$export_formats) {
  saveRDS(scRNA, file = "processed_data/scRNA_unfiltered.rds")
}

# Export H5 format (compatible with SingleCellExperiment)
if("h5" %in% params$export_formats) {
  dior::write_h5(
    scRNA, 
    file = "processed_data/scRNA_unfiltered.h5",
    object.type = "seurat",
    overwrite = TRUE
  )
  cat("H5 file exported: processed_data/scRNA_unfiltered.h5\n")
}

# Report dataset statistics
cat("\nDataset composition:\n")
sample_counts <- table(scRNA$orig.ident)
print(sample_counts)

# Export sample distribution table
write.csv(data.frame(Sample = names(sample_counts), 
          Cells = as.integer(sample_counts)),
          file = "results/tables/sample_distribution.csv",
          row.names = FALSE)

# --------------------------------
# 4. QUALITY CONTROL VISUALIZATION
# --------------------------------
cat("\n>> Generating QC diagnostic plots...\n")

# Dynamic feature detection for QC plots
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")
if("percent.HB" %in% colnames(scRNA@meta.data)) {
  qc_features <- c(qc_features, "percent.HB")
}

# Generate publication-quality QC plot
qc_plot <- VlnPlot(scRNA, features = qc_features,
                   pt.size = 0.1, ncol = 3) +
  theme_minimal() +
  ggtitle("Pre-filtering Quality Control Metrics")

ggsave("results/figures/Fig1_QC_metrics.pdf", qc_plot, 
       width = 12, height = 8, dpi = 300)

# Save complete session information
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
cat("\nAnalysis completed successfully. Outputs saved to processed_data/ and results/\n")
