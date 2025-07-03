# H5AD to Seurat Converter

[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-v5-green.svg)](https://satijalab.org/seurat/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A comprehensive toolkit for converting AnnData objects (h5ad files) to Seurat objects, designed for seamless integration between Python-based scanpy workflows and R-based Seurat analysis pipelines.

## 🚀 Features

- **Large Dataset Support**: Efficiently handles datasets with 90k+ cells
- **Complete Metadata Preservation**: Maintains all cell and gene metadata
- **Embedding Conservation**: Preserves UMAP, PCA, and t-SNE coordinates
- **Seurat v5 Compatible**: Optimized for the latest Seurat version
- **Dual Loading Methods**: Supports both SeuratDisk (h5ad) and CSV fallback
- **Automated Workflows**: Batch processing and integrated analysis scripts
- **Quality Assurance**: Comprehensive testing and validation tools

## 📋 Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage Examples](#usage-examples)
- [Output Files](#output-files)
- [Advanced Usage](#advanced-usage)
- [Testing](#testing)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [Citation](#citation)

## 🛠 Installation

### Python Dependencies

```bash
# Install required Python packages
pip install scanpy pandas numpy anndata

# Or using conda
conda install -c conda-forge scanpy pandas numpy anndata
```

### R Dependencies

```r
# Install Seurat v5
remotes::install_version("SeuratObject", version = "4.1.4", 
                        repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("Seurat", version = "4.4.0", 
                        repos = c("https://satijalab.r-universe.dev", getOption("repos")))

# Install SeuratDisk for h5ad support (recommended)
remotes::install_github("mojaveazure/seurat-disk")

# Additional packages for visualization
install.packages(c("dplyr", "ggplot2", "patchwork"))
```

## ⚡ Quick Start

### 1. Basic Conversion

```bash
# Convert a single h5ad file
python h5ad_to_seurat_converter.py your_data.h5ad

# With custom output directory and prefix
python h5ad_to_seurat_converter.py your_data.h5ad -o output_dir -p my_data
```

### 2. Load in Seurat

```r
# Method 1: Using SeuratDisk (recommended)
library(Seurat)
library(SeuratDisk)
seurat_obj <- LoadH5Seurat("seurat_conversion/converted_data.h5ad")

# Method 2: Using auto-generated script
source("seurat_conversion/converted_data_load_seurat.R")
```

### 3. Verify Conversion

```r
# Check the object
print(seurat_obj)
DimPlot(seurat_obj, reduction = "umap")
```

## 📖 Usage Examples

### Single File Conversion

```bash
# Basic conversion
python h5ad_to_seurat_converter.py data/adata_celltype.h5ad

# Advanced options
python h5ad_to_seurat_converter.py data/adata_celltype.h5ad \
  --output_dir seurat_objects \
  --prefix experiment_1
```

### Batch Processing

```bash
# Run automated workflow for multiple files
python example_workflow.py
```

### Testing Your Data

```bash
# Test conversion capabilities
python test_h5ad_conversion.py

# Check dependencies
python h5ad_to_seurat_converter.py --check-deps
```

## 📁 Output Files

Each conversion generates the following files:

### Core Files
- **`converted_data.h5ad`** - Optimized AnnData file for SeuratDisk
- **`converted_data_load_seurat.R`** - Auto-generated R loading script

### Metadata Files
- **`converted_data_metadata.csv`** - Cell metadata (barcodes, cell types, QC metrics)
- **`converted_data_gene_metadata.csv`** - Gene information and annotations

### Embeddings (when available)
- **`converted_data_umap.csv`** - UMAP coordinates
- **`converted_data_pca.csv`** - PCA coordinates

### Large Datasets
For datasets >10k cells or >5k genes, count matrices are kept in efficient h5ad format rather than CSV to optimize performance.

## 🔬 Advanced Usage

### Custom Workflow Integration

```python
from h5ad_to_seurat_converter import load_h5ad_file, prepare_seurat_data, export_for_seurat

# Load and process
adata = load_h5ad_file("your_data.h5ad")
adata = prepare_seurat_data(adata)
export_for_seurat(adata, "output_dir", "prefix")
```

### R Integration Workflow

```r
# Load multiple converted objects
source("testis_analysis_workflow.R")

# Access individual objects
cell_annotation_obj <- seurat_objects[["testis_cell_annotation"]]
integrated_obj <- seurat_objects[["testis_integrated_scvi"]]

# Comparative analysis
comparison_plot <- wrap_plots(
  DimPlot(cell_annotation_obj, group.by = "celltype"),
  DimPlot(integrated_obj, group.by = "celltype"),
  ncol = 2
)
```

## 🧪 Testing

### Validate Installation

```bash
# Check Python dependencies
python h5ad_to_seurat_converter.py --check-deps

# Test with sample data
python test_h5ad_conversion.py
```

### R Environment Testing

```r
# Test Seurat v5 compatibility
source("test_seurat_v5_conversion.R")
```

### Performance Benchmarks

The toolkit has been tested with:
- ✅ **90k+ cells** (human testis dataset)
- ✅ **37k+ genes** (full transcriptome)
- ✅ **Multiple embeddings** (UMAP, PCA, t-SNE)
- ✅ **Complex metadata** (20+ annotation columns)

## 🔧 Troubleshooting

### Common Issues

**1. SeuratDisk not available**
```r
# Install SeuratDisk
remotes::install_github("mojaveazure/seurat-disk")
```

**2. Memory issues with large datasets**
- The toolkit automatically uses sparse matrices
- For very large datasets (>100k cells), consider subsetting

**3. Gene name conflicts**
```python
# Automatically handled by the converter
adata.var_names_make_unique()
```

**4. Missing embeddings**
- Check if embeddings exist in original h5ad file
- Embeddings are preserved when available

### Error Messages

The toolkit provides detailed error reporting:
- File validation errors
- Memory usage warnings  
- Compatibility checks
- Conversion progress tracking

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/h5ad-to-seurat-converter.git
cd h5ad-to-seurat-converter

# Install development dependencies
pip install -r requirements-dev.txt

# Run tests
python -m pytest tests/
```

## 📊 Use Cases

This toolkit is particularly useful for:

- **Cross-platform Analysis**: Moving between Python (scanpy) and R (Seurat) workflows
- **Collaborative Research**: Sharing data between teams using different tools
- **Method Comparison**: Comparing results from different analysis pipelines
- **Publication Workflows**: Preparing data for different visualization requirements

## 🏆 Success Stories

The toolkit has been successfully used for:
- Human testis single-cell RNA-seq analysis (90k+ cells)
- Multi-sample integration studies
- Cell type annotation workflows
- Developmental trajectory analysis

## 📚 Documentation

- [Detailed Usage Guide](README_h5ad_conversion.md)
- [API Documentation](docs/api.md)
- [Troubleshooting Guide](docs/troubleshooting.md)
- [Best Practices](docs/best_practices.md)

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- [Scanpy](https://scanpy.readthedocs.io/) for Python single-cell analysis
- [Seurat](https://satijalab.org/seurat/) for R single-cell analysis
- [SeuratDisk](https://github.com/mojaveazure/seurat-disk) for h5ad support in R
- The single-cell genomics community for valuable feedback

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/yourusername/h5ad-to-seurat-converter/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/h5ad-to-seurat-converter/discussions)
- **Email**: your.email@domain.com

## 📈 Citation

If you use this toolkit in your research, please cite:

```bibtex
@software{h5ad_to_seurat_converter,
  title={H5AD to Seurat Converter: Seamless Integration Between Python and R Single-Cell Analysis},
  author={Your Name},
  year={2025},
  url={https://github.com/yourusername/h5ad-to-seurat-converter}
}
```

---

**Made with ❤️ for the single-cell genomics community**
