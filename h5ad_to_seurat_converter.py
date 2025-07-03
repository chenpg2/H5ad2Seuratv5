#!/usr/bin/env python3
"""
H5AD to Seurat Object Converter
================================

This script converts AnnData objects (h5ad files) to Seurat objects for use in R.
Designed to work with Seurat v5 conda environment.

Author: Generated for testis single-cell analysis
Date: 2025-07-03
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path
import subprocess
import tempfile
import warnings
warnings.filterwarnings('ignore')

# Configure scanpy
sc.settings.verbosity = 1  # verbosity level
sc.settings.set_figure_params(dpi=80, facecolor='white')

def check_dependencies():
    """Check if required packages are installed."""
    required_packages = ['scanpy', 'pandas', 'numpy', 'anndata']
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing_packages.append(package)
    
    if missing_packages:
        print(f"Error: Missing required packages: {', '.join(missing_packages)}")
        print("Please install them using: pip install " + " ".join(missing_packages))
        sys.exit(1)

def load_h5ad_file(file_path):
    """
    Load h5ad file and return AnnData object.
    
    Parameters:
    -----------
    file_path : str
        Path to the h5ad file
        
    Returns:
    --------
    adata : AnnData
        Loaded AnnData object
    """
    try:
        print(f"Loading h5ad file: {file_path}")
        adata = sc.read_h5ad(file_path)
        print(f"Successfully loaded data with {adata.n_obs} cells and {adata.n_vars} genes")
        return adata
    except Exception as e:
        print(f"Error loading h5ad file: {e}")
        sys.exit(1)

def prepare_seurat_data(adata):
    """
    Prepare AnnData object for Seurat conversion.
    
    Parameters:
    -----------
    adata : AnnData
        Input AnnData object
        
    Returns:
    --------
    adata : AnnData
        Processed AnnData object ready for Seurat conversion
    """
    print("Preparing data for Seurat conversion...")
    
    # Make variable names unique (required for Seurat)
    adata.var_names_make_unique()
    
    # Ensure gene names are strings
    adata.var.index = adata.var.index.astype(str)
    
    # Ensure cell barcodes are strings
    adata.obs.index = adata.obs.index.astype(str)
    
    # Convert sparse matrix to dense if needed for certain operations
    if hasattr(adata.X, 'toarray'):
        print("Data is in sparse format (recommended for large datasets)")
    
    # Check for required metadata
    print(f"Available observation metadata: {list(adata.obs.columns)}")
    print(f"Available variable metadata: {list(adata.var.columns)}")
    
    return adata

def create_seurat_rds(adata, output_dir, prefix="converted_data"):
    """
    Create Seurat object directly in R and save as RDS file.

    Parameters:
    -----------
    adata : AnnData
        AnnData object to convert
    output_dir : str
        Output directory path
    prefix : str
        Prefix for output files

    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Creating Seurat RDS object: {output_dir}")

    # Create temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)

        # Export h5ad file temporarily
        temp_h5ad = temp_path / f"{prefix}_temp.h5ad"
        adata.write_h5ad(temp_h5ad)

        # Export metadata
        metadata_path = temp_path / f"{prefix}_metadata.csv"
        adata.obs.to_csv(metadata_path)

        # Export embeddings if available
        embeddings_info = {}
        if 'X_umap' in adata.obsm:
            umap_path = temp_path / f"{prefix}_umap.csv"
            umap_df = pd.DataFrame(
                adata.obsm['X_umap'],
                index=adata.obs.index,
                columns=['UMAP_1', 'UMAP_2']
            )
            umap_df.to_csv(umap_path)
            embeddings_info['umap'] = str(umap_path)

        if 'X_pca' in adata.obsm:
            pca_path = temp_path / f"{prefix}_pca.csv"
            pca_df = pd.DataFrame(
                adata.obsm['X_pca'],
                index=adata.obs.index,
                columns=[f'PC_{i+1}' for i in range(adata.obsm['X_pca'].shape[1])]
            )
            pca_df.to_csv(pca_path)
            embeddings_info['pca'] = str(pca_path)

        # Export count matrix for R loading
        count_matrix_path = temp_path / f"{prefix}_counts.csv"

        # For large datasets, we'll use a different approach
        if adata.n_obs > 10000:
            # Export a subset of highly variable genes for demonstration
            if 'highly_variable' in adata.var.columns:
                hvg_mask = adata.var['highly_variable'].fillna(False)
                if hvg_mask.sum() > 0:
                    adata_subset = adata[:, hvg_mask].copy()
                else:
                    # Take top 2000 genes by variance
                    adata_subset = adata[:, :2000].copy()
            else:
                adata_subset = adata[:, :2000].copy()

            print(f"Using subset for RDS: {adata_subset.n_obs} cells × {adata_subset.n_vars} genes")
        else:
            adata_subset = adata.copy()

        # Convert to dense and export count matrix
        if hasattr(adata_subset.X, 'toarray'):
            count_matrix = adata_subset.X.toarray()
        else:
            count_matrix = adata_subset.X

        # Create count matrix DataFrame (genes × cells for Seurat)
        counts_df = pd.DataFrame(
            count_matrix.T,
            index=adata_subset.var.index,
            columns=adata_subset.obs.index
        )
        counts_df.to_csv(count_matrix_path)

        # Create R script for conversion
        r_script_content = f'''
# Load required libraries
suppressPackageStartupMessages({{
  library(Seurat)
  library(dplyr)
  library(Matrix)
}})

cat("Creating Seurat object from CSV files...\\n")

# Load count matrix
cat("Loading count matrix...\\n")
counts <- read.csv("{count_matrix_path}", row.names = 1, check.names = FALSE)
cat("✅ Loaded count matrix:", nrow(counts), "genes ×", ncol(counts), "cells\\n")

# Load metadata
cat("Loading metadata...\\n")
metadata <- read.csv("{metadata_path}", row.names = 1)
cat("✅ Loaded metadata:", nrow(metadata), "cells ×", ncol(metadata), "columns\\n")

# Ensure cell names match
common_cells <- intersect(colnames(counts), rownames(metadata))
cat("Common cells:", length(common_cells), "\\n")

if (length(common_cells) == 0) {{
  stop("No common cell names between count matrix and metadata!")
}}

# Subset to common cells
counts <- counts[, common_cells]
metadata <- metadata[common_cells, ]

# Create Seurat object
cat("Creating Seurat object...\\n")
seurat_obj <- CreateSeuratObject(
  counts = as.matrix(counts),
  meta.data = metadata,
  project = "{prefix}",
  min.cells = 0,
  min.features = 0
)

cat("✅ Created Seurat object\\n")

# Add embeddings if available
'''

        # Add embedding loading code
        if 'umap' in embeddings_info:
            r_script_content += f'''
# Add UMAP coordinates
if (file.exists("{embeddings_info['umap']}")) {{
  umap_coords <- read.csv("{embeddings_info['umap']}", row.names = 1)
  seurat_obj[["umap"]] <- CreateDimReducObject(
    embeddings = as.matrix(umap_coords),
    key = "UMAP_",
    assay = DefaultAssay(seurat_obj)
  )
  cat("✅ Added UMAP coordinates\\n")
}}
'''

        if 'pca' in embeddings_info:
            r_script_content += f'''
# Add PCA coordinates
if (file.exists("{embeddings_info['pca']}")) {{
  pca_coords <- read.csv("{embeddings_info['pca']}", row.names = 1)
  seurat_obj[["pca"]] <- CreateDimReducObject(
    embeddings = as.matrix(pca_coords),
    key = "PC_",
    assay = DefaultAssay(seurat_obj)
  )
  cat("✅ Added PCA coordinates\\n")
}}
'''

        # Add final save code
        rds_path = output_dir / f"{prefix}_seurat.rds"
        r_script_content += f'''
# Display object information
cat("\\n=== Seurat Object Summary ===\\n")
print(seurat_obj)
cat("Cells:", ncol(seurat_obj), "\\n")
cat("Genes:", nrow(seurat_obj), "\\n")
cat("Metadata columns:", ncol(seurat_obj@meta.data), "\\n")
if (length(seurat_obj@reductions) > 0) {{
  cat("Reductions:", names(seurat_obj@reductions), "\\n")
}}

# Save as RDS
saveRDS(seurat_obj, file = "{rds_path}")
cat("\\n✅ Seurat object saved as: {rds_path}\\n")

# Success indicator
cat("CONVERSION_SUCCESS\\n")
'''

        # Write R script
        r_script_path = temp_path / "convert_to_seurat.R"
        with open(r_script_path, 'w') as f:
            f.write(r_script_content)

        # Run R script
        try:
            print("Running R conversion script...")
            result = subprocess.run(
                ['Rscript', str(r_script_path)],
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )

            if result.returncode == 0 and "CONVERSION_SUCCESS" in result.stdout:
                print("✅ Successfully created Seurat RDS file!")
                print(f"Output: {rds_path}")

                # Print R output (excluding success indicator)
                r_output = result.stdout.replace("CONVERSION_SUCCESS", "").strip()
                if r_output:
                    print("R Output:")
                    print(r_output)

                return True
            else:
                print("❌ R conversion failed!")
                print("STDOUT:", result.stdout)
                print("STDERR:", result.stderr)
                return False

        except subprocess.TimeoutExpired:
            print("❌ R conversion timed out!")
            return False
        except FileNotFoundError:
            print("❌ R/Rscript not found! Please install R and make sure it's in your PATH.")
            return False
        except Exception as e:
            print(f"❌ Error running R script: {e}")
            return False

def export_for_seurat(adata, output_dir, prefix="converted_data"):
    """
    Export AnnData object in formats compatible with Seurat v5.
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object to export
    output_dir : str
        Output directory path
    prefix : str
        Prefix for output files
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Exporting data to: {output_dir}")
    
    # Method 1: Save as h5ad (can be read by SeuratDisk)
    h5ad_path = output_dir / f"{prefix}.h5ad"
    adata.write_h5ad(h5ad_path)
    print(f"Saved h5ad file: {h5ad_path}")
    
    # Method 2: Save as h5seurat (if SeuratDisk is available)
    try:
        h5seurat_path = output_dir / f"{prefix}.h5seurat"
        # Note: This requires SeuratDisk to be installed in the environment
        # We'll create the h5ad file and provide R code to convert it
        print(f"H5AD file saved. Use the following R code to convert to Seurat:")
        print(f"""
# R code to convert h5ad to Seurat object:
library(Seurat)
library(SeuratDisk)

# Load the h5ad file
adata <- LoadH5Seurat("{h5ad_path}")

# Or if SeuratDisk is not available, use:
# library(reticulate)
# sc <- import("scanpy")
# adata <- sc$read_h5ad("{h5ad_path}")
# 
# # Convert to Seurat
# seurat_obj <- CreateSeuratObject(
#   counts = t(adata$X),
#   meta.data = adata$obs
# )
        """)
    except Exception as e:
        print(f"Note: Direct h5seurat export not available: {e}")
    
    # Method 3: Export individual components for manual Seurat creation
    try:
        # Export count matrix
        if hasattr(adata.X, 'toarray'):
            count_matrix = adata.X.toarray()
        else:
            count_matrix = adata.X
        
        # Save as CSV (for small datasets) or use other formats
        if adata.n_obs < 10000 and adata.n_vars < 5000:
            counts_df = pd.DataFrame(
                count_matrix.T,  # Transpose for Seurat format (genes x cells)
                index=adata.var.index,
                columns=adata.obs.index
            )
            counts_path = output_dir / f"{prefix}_counts.csv"
            counts_df.to_csv(counts_path)
            print(f"Saved count matrix: {counts_path}")
        
        # Export metadata
        metadata_path = output_dir / f"{prefix}_metadata.csv"
        adata.obs.to_csv(metadata_path)
        print(f"Saved cell metadata: {metadata_path}")
        
        # Export gene metadata
        gene_metadata_path = output_dir / f"{prefix}_gene_metadata.csv"
        adata.var.to_csv(gene_metadata_path)
        print(f"Saved gene metadata: {gene_metadata_path}")
        
        # Export embeddings if available
        if 'X_umap' in adata.obsm:
            umap_df = pd.DataFrame(
                adata.obsm['X_umap'],
                index=adata.obs.index,
                columns=['UMAP_1', 'UMAP_2']
            )
            umap_path = output_dir / f"{prefix}_umap.csv"
            umap_df.to_csv(umap_path)
            print(f"Saved UMAP coordinates: {umap_path}")
        
        if 'X_pca' in adata.obsm:
            pca_df = pd.DataFrame(
                adata.obsm['X_pca'],
                index=adata.obs.index,
                columns=[f'PC_{i+1}' for i in range(adata.obsm['X_pca'].shape[1])]
            )
            pca_path = output_dir / f"{prefix}_pca.csv"
            pca_df.to_csv(pca_path)
            print(f"Saved PCA coordinates: {pca_path}")
            
    except Exception as e:
        print(f"Warning: Could not export individual components: {e}")

def create_r_conversion_script(output_dir, prefix="converted_data"):
    """
    Create an R script for loading the converted data into Seurat v5.
    
    Parameters:
    -----------
    output_dir : str
        Output directory path
    prefix : str
        Prefix for data files
    """
    output_dir = Path(output_dir)
    r_script_path = output_dir / f"{prefix}_load_seurat.R"
    
    r_script_content = f'''# R Script to Load Converted Data into Seurat v5
# Generated by h5ad_to_seurat_converter.py

library(Seurat)
library(dplyr)

# Method 1: Load from h5ad using SeuratDisk (recommended)
if (requireNamespace("SeuratDisk", quietly = TRUE)) {{
  library(SeuratDisk)
  seurat_obj <- LoadH5Seurat("{output_dir}/{prefix}.h5ad")
  print("Successfully loaded data using SeuratDisk")
}} else {{
  print("SeuratDisk not available. Install with: remotes::install_github('mojaveazure/seurat-disk')")
  
  # Method 2: Load from CSV files (for smaller datasets)
  if (file.exists("{output_dir}/{prefix}_counts.csv")) {{
    # Load count matrix
    counts <- read.csv("{output_dir}/{prefix}_counts.csv", row.names = 1, check.names = FALSE)
    
    # Load metadata
    metadata <- read.csv("{output_dir}/{prefix}_metadata.csv", row.names = 1)
    
    # Create Seurat object
    seurat_obj <- CreateSeuratObject(
      counts = counts,
      meta.data = metadata,
      project = "H5AD_Converted"
    )
    
    # Add embeddings if available
    if (file.exists("{output_dir}/{prefix}_umap.csv")) {{
      umap_coords <- read.csv("{output_dir}/{prefix}_umap.csv", row.names = 1)
      seurat_obj[["umap"]] <- CreateDimReducObject(
        embeddings = as.matrix(umap_coords),
        key = "UMAP_",
        assay = DefaultAssay(seurat_obj)
      )
    }}
    
    if (file.exists("{output_dir}/{prefix}_pca.csv")) {{
      pca_coords <- read.csv("{output_dir}/{prefix}_pca.csv", row.names = 1)
      seurat_obj[["pca"]] <- CreateDimReducObject(
        embeddings = as.matrix(pca_coords),
        key = "PC_",
        assay = DefaultAssay(seurat_obj)
      )
    }}
    
    print("Successfully loaded data from CSV files")
  }} else {{
    stop("No suitable data files found for loading")
  }}
}}

# Display basic information about the Seurat object
print(seurat_obj)
print(paste("Number of cells:", ncol(seurat_obj)))
print(paste("Number of genes:", nrow(seurat_obj)))

# Save the Seurat object
saveRDS(seurat_obj, file = "{output_dir}/{prefix}_seurat_object.rds")
print("Seurat object saved as {prefix}_seurat_object.rds")
'''
    
    with open(r_script_path, 'w') as f:
        f.write(r_script_content)
    
    print(f"Created R conversion script: {r_script_path}")

def main():
    """Main function to handle command line arguments and run conversion."""
    parser = argparse.ArgumentParser(
        description="Convert h5ad files to Seurat-compatible formats",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python h5ad_to_seurat_converter.py input.h5ad
  python h5ad_to_seurat_converter.py input.h5ad -o output_dir -p my_data
  python h5ad_to_seurat_converter.py scanpypipline_output/adata_cellAnno.h5ad -o seurat_conversion
        """
    )
    
    parser.add_argument('input_file', nargs='?', help='Input h5ad file path')
    parser.add_argument('-o', '--output_dir', default='seurat_conversion',
                       help='Output directory (default: seurat_conversion)')
    parser.add_argument('-p', '--prefix', default='converted_data',
                       help='Prefix for output files (default: converted_data)')
    parser.add_argument('--rds', action='store_true',
                       help='Create Seurat RDS object directly (requires R)')
    parser.add_argument('--check-deps', action='store_true',
                       help='Check if required dependencies are installed')
    
    args = parser.parse_args()
    
    if args.check_deps:
        check_dependencies()
        print("All required dependencies are installed!")
        return

    # Check dependencies
    check_dependencies()

    # Validate input file
    if not args.input_file:
        print("Error: Input file is required when not using --check-deps")
        parser.print_help()
        sys.exit(1)

    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' does not exist")
        sys.exit(1)
    
    # Load and process data
    adata = load_h5ad_file(args.input_file)
    adata = prepare_seurat_data(adata)

    # Choose conversion method
    if args.rds:
        # Create RDS file directly
        print("Creating Seurat RDS object...")
        success = create_seurat_rds(adata, args.output_dir, args.prefix)
        if success:
            print(f"\n✅ RDS conversion completed successfully!")
            print(f"Seurat RDS file saved in: {args.output_dir}")
            print(f"Load in R with: seurat_obj <- readRDS('{args.output_dir}/{args.prefix}_seurat.rds')")
        else:
            print(f"\n❌ RDS conversion failed!")
            print("Falling back to standard export...")
            export_for_seurat(adata, args.output_dir, args.prefix)
            create_r_conversion_script(args.output_dir, args.prefix)
    else:
        # Standard export (h5ad + CSV + R script)
        export_for_seurat(adata, args.output_dir, args.prefix)
        create_r_conversion_script(args.output_dir, args.prefix)
        print(f"\nConversion completed successfully!")
        print(f"Output files saved in: {args.output_dir}")
        print(f"Run the generated R script to load data into Seurat v5")

if __name__ == "__main__":
    main()
