#!/usr/bin/env Rscript
# Create Seurat RDS from CSV files
# Usage: Rscript create_seurat_from_csv.R metadata.csv umap.csv pca.csv output_prefix

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  cat("Usage: Rscript create_seurat_from_csv.R metadata.csv umap.csv pca.csv output_prefix\n")
  cat("Example: Rscript create_seurat_from_csv.R metadata.csv umap.csv pca.csv testis_data\n")
  quit(status = 1)
}

metadata_file <- args[1]
umap_file <- args[2]
pca_file <- args[3]
output_prefix <- args[4]

cat("=== CSV to Seurat RDS Converter ===\n")
cat("Metadata file:", metadata_file, "\n")
cat("UMAP file:", umap_file, "\n")
cat("PCA file:", pca_file, "\n")
cat("Output prefix:", output_prefix, "\n\n")

# Check if files exist
files_to_check <- c(metadata_file, umap_file, pca_file)
for (file in files_to_check) {
  if (!file.exists(file)) {
    cat("Error: File does not exist:", file, "\n")
    quit(status = 1)
  }
}

# Load metadata
cat("Loading metadata...\n")
metadata <- read.csv(metadata_file, row.names = 1)
cat("✅ Loaded metadata:", nrow(metadata), "cells ×", ncol(metadata), "columns\n")

# Display metadata summary
cat("Metadata columns:", paste(colnames(metadata), collapse = ", "), "\n")

# Check for cell type information
if ("celltype" %in% colnames(metadata)) {
  cat("Cell types found:", length(unique(metadata$celltype)), "\n")
  celltype_counts <- table(metadata$celltype)
  cat("Cell type distribution:\n")
  print(celltype_counts)
}

# Check for sample information
if ("sample" %in% colnames(metadata)) {
  cat("Samples found:", length(unique(metadata$sample)), "\n")
  sample_counts <- table(metadata$sample)
  cat("Sample distribution:\n")
  print(sample_counts)
}

# Create a minimal count matrix for demonstration
# In a real scenario, you would load actual count data
cat("\nCreating minimal Seurat object...\n")
n_cells <- nrow(metadata)
n_genes <- min(2000, n_cells)  # Use reasonable number of genes

# Create a sparse random count matrix as placeholder
set.seed(42)  # For reproducibility
counts <- Matrix::rsparsematrix(
  nrow = n_genes,
  ncol = n_cells,
  density = 0.1,  # 10% non-zero values
  rand.x = function(n) rpois(n, lambda = 5)  # Poisson distributed counts
)

# Set row and column names
rownames(counts) <- paste0("Gene_", 1:n_genes)
colnames(counts) <- rownames(metadata)

# Create Seurat object
cat("Creating Seurat object...\n")
seurat_obj <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata,
  project = output_prefix,
  min.cells = 0,
  min.features = 0
)

cat("✅ Created Seurat object with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes\n")

# Load and add UMAP coordinates
if (file.exists(umap_file)) {
  cat("Loading UMAP coordinates...\n")
  umap_coords <- read.csv(umap_file, row.names = 1)
  cat("✅ Loaded UMAP:", nrow(umap_coords), "cells ×", ncol(umap_coords), "dimensions\n")
  
  # Ensure cell names match
  common_cells <- intersect(colnames(seurat_obj), rownames(umap_coords))
  if (length(common_cells) > 0) {
    umap_coords <- umap_coords[common_cells, ]
    
    # Add UMAP reduction
    seurat_obj[["umap"]] <- CreateDimReducObject(
      embeddings = as.matrix(umap_coords),
      key = "UMAP_",
      assay = DefaultAssay(seurat_obj)
    )
    cat("✅ Added UMAP reduction\n")
  } else {
    cat("⚠️ No matching cells between Seurat object and UMAP coordinates\n")
  }
}

# Load and add PCA coordinates
if (file.exists(pca_file)) {
  cat("Loading PCA coordinates...\n")
  pca_coords <- read.csv(pca_file, row.names = 1)
  cat("✅ Loaded PCA:", nrow(pca_coords), "cells ×", ncol(pca_coords), "dimensions\n")
  
  # Ensure cell names match
  common_cells <- intersect(colnames(seurat_obj), rownames(pca_coords))
  if (length(common_cells) > 0) {
    pca_coords <- pca_coords[common_cells, ]
    
    # Add PCA reduction
    seurat_obj[["pca"]] <- CreateDimReducObject(
      embeddings = as.matrix(pca_coords),
      key = "PC_",
      assay = DefaultAssay(seurat_obj)
    )
    cat("✅ Added PCA reduction\n")
  } else {
    cat("⚠️ No matching cells between Seurat object and PCA coordinates\n")
  }
}

# Display final object information
cat("\n=== Final Seurat Object Summary ===\n")
print(seurat_obj)
cat("Cells:", ncol(seurat_obj), "\n")
cat("Genes:", nrow(seurat_obj), "\n")
cat("Assays:", names(seurat_obj@assays), "\n")
cat("Metadata columns:", ncol(seurat_obj@meta.data), "\n")

if (length(seurat_obj@reductions) > 0) {
  cat("Reductions:", names(seurat_obj@reductions), "\n")
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
cat("Created from CSV files\n")
cat("Metadata file:", metadata_file, "\n")
cat("UMAP file:", umap_file, "\n")
cat("PCA file:", pca_file, "\n")
cat("Output file:", rds_file, "\n")
cat("Creation date:", Sys.time(), "\n\n")
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
if ("sample" %in% colnames(seurat_obj@meta.data)) {
  cat("\nSample distribution:\n")
  print(table(seurat_obj@meta.data$sample))
}
sink()

cat("Summary saved as:", summary_file, "\n")

# Create a simple R script to load and visualize
load_script <- paste0(output_prefix, "_load_and_plot.R")
cat("# Load and visualize Seurat object\n", file = load_script)
cat("library(Seurat)\n", file = load_script, append = TRUE)
cat("library(ggplot2)\n", file = load_script, append = TRUE)
cat("\n", file = load_script, append = TRUE)
cat("# Load the Seurat object\n", file = load_script, append = TRUE)
cat(paste0("seurat_obj <- readRDS('", rds_file, "')\n"), file = load_script, append = TRUE)
cat("\n", file = load_script, append = TRUE)
cat("# Display object information\n", file = load_script, append = TRUE)
cat("print(seurat_obj)\n", file = load_script, append = TRUE)
cat("\n", file = load_script, append = TRUE)
cat("# Create UMAP plot\n", file = load_script, append = TRUE)
cat("if ('umap' %in% names(seurat_obj@reductions)) {\n", file = load_script, append = TRUE)
cat("  p1 <- DimPlot(seurat_obj, reduction = 'umap', group.by = 'celltype') + \n", file = load_script, append = TRUE)
cat("    ggtitle('UMAP by Cell Type')\n", file = load_script, append = TRUE)
cat("  print(p1)\n", file = load_script, append = TRUE)
cat("  \n", file = load_script, append = TRUE)
cat("  # Save plot\n", file = load_script, append = TRUE)
cat(paste0("  ggsave('", output_prefix, "_umap_celltype.pdf', p1, width = 12, height = 8)\n"), file = load_script, append = TRUE)
cat("}\n", file = load_script, append = TRUE)

cat("Load script created:", load_script, "\n")

cat("\n🎉 Conversion completed successfully!\n")
cat("\n=== Next Steps ===\n")
cat("1. Load in R: seurat_obj <- readRDS('", rds_file, "')\n", sep = "")
cat("2. Run visualization script: source('", load_script, "')\n", sep = "")
cat("3. Check summary: cat(readLines('", summary_file, "'), sep = '\\n')\n", sep = "")
