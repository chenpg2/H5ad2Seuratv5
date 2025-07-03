#!/usr/bin/env Rscript
# Convert H5AD to Seurat RDS
# Usage: Rscript convert_h5ad_to_rds.R input.h5ad output_prefix

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript convert_h5ad_to_rds.R input.h5ad output_prefix\n")
  cat("Example: Rscript convert_h5ad_to_rds.R data.h5ad testis_data\n")
  quit(status = 1)
}

input_file <- args[1]
output_prefix <- args[2]

cat("=== H5AD to Seurat RDS Converter ===\n")
cat("Input file:", input_file, "\n")
cat("Output prefix:", output_prefix, "\n\n")

# Check if input file exists
if (!file.exists(input_file)) {
  cat("Error: Input file does not exist:", input_file, "\n")
  quit(status = 1)
}

# Try to load with SeuratDisk first
if (requireNamespace("SeuratDisk", quietly = TRUE)) {
  library(SeuratDisk)
  
  cat("Attempting to load with SeuratDisk...\n")
  tryCatch({
    seurat_obj <- LoadH5Seurat(input_file)
    cat("✅ Successfully loaded with SeuratDisk\n")
    
    # Display object information
    cat("\n=== Seurat Object Summary ===\n")
    print(seurat_obj)
    cat("Cells:", ncol(seurat_obj), "\n")
    cat("Genes:", nrow(seurat_obj), "\n")
    cat("Assays:", names(seurat_obj@assays), "\n")
    
    # Check metadata
    cat("Metadata columns:", ncol(seurat_obj@meta.data), "\n")
    if ("celltype" %in% colnames(seurat_obj@meta.data)) {
      cat("Cell types found:", length(unique(seurat_obj@meta.data$celltype)), "\n")
      cat("Cell type distribution:\n")
      print(table(seurat_obj@meta.data$celltype))
    }
    
    # Check reductions
    if (length(seurat_obj@reductions) > 0) {
      cat("Available reductions:", names(seurat_obj@reductions), "\n")
    }
    
    # Save as RDS
    rds_file <- paste0(output_prefix, "_seurat.rds")
    saveRDS(seurat_obj, file = rds_file)
    cat("\n✅ Seurat object saved as:", rds_file, "\n")
    
    # Create a summary file
    summary_file <- paste0(output_prefix, "_summary.txt")
    sink(summary_file)
    cat("Seurat Object Summary\n")
    cat("====================\n")
    cat("Input file:", input_file, "\n")
    cat("Output file:", rds_file, "\n")
    cat("Conversion date:", Sys.time(), "\n\n")
    print(seurat_obj)
    cat("\nMetadata columns:\n")
    cat(paste(colnames(seurat_obj@meta.data), collapse = ", "), "\n")
    if (length(seurat_obj@reductions) > 0) {
      cat("\nReductions:\n")
      cat(paste(names(seurat_obj@reductions), collapse = ", "), "\n")
    }
    if ("celltype" %in% colnames(seurat_obj@meta.data)) {
      cat("\nCell type distribution:\n")
      print(table(seurat_obj@meta.data$celltype))
    }
    sink()
    
    cat("Summary saved as:", summary_file, "\n")
    cat("\n🎉 Conversion completed successfully!\n")
    
  }, error = function(e) {
    cat("❌ SeuratDisk loading failed:", e$message, "\n")
    cat("This might be due to h5ad format compatibility issues.\n")
    cat("Please try converting the h5ad file using the Python script first.\n")
    quit(status = 1)
  })
  
} else {
  cat("❌ SeuratDisk not available\n")
  cat("Please install SeuratDisk:\n")
  cat("remotes::install_github('mojaveazure/seurat-disk')\n")
  quit(status = 1)
}

cat("\n=== Conversion Complete ===\n")
cat("To load in R:\n")
cat("seurat_obj <- readRDS('", rds_file, "')\n", sep = "")
cat("print(seurat_obj)\n")
cat("DimPlot(seurat_obj, reduction = 'umap')\n")
