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
    
    # Export data
    export_for_seurat(adata, args.output_dir, args.prefix)
    
    # Create R conversion script
    create_r_conversion_script(args.output_dir, args.prefix)
    
    print(f"\nConversion completed successfully!")
    print(f"Output files saved in: {args.output_dir}")
    print(f"Run the generated R script to load data into Seurat v5")

if __name__ == "__main__":
    main()
