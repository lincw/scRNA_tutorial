# Hands-On Examples

Follow these examples to learn by doing with the PBMC 3K dataset. Notebooks and scripts are located in the `analysis/` directory, which contains detailed analysis workflows for both Python and R.

## Python Notebook (`analysis/pbmc_analysis.ipynb`)
This Jupyter notebook walks through:
1. Loading 10× PBMC data
2. Quality control and filtering
3. Normalization and scaling
4. Dimensionality reduction (PCA, UMAP)
5. Clustering and visualization

Key excerpts:
```python
import scanpy as sc
# Load the data from MTX format
adata = sc.read_10x_mtx('../data/filtered_gene_bc_matrices/hg19/', var_names='gene_symbols')
# QC
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=2500)  # Remove potential doublets
adata = adata[adata.obs.pct_counts_mt < 5, :]
# Normalize & log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# PCA & UMAP
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color='leiden')
```

## R Markdown (`analysis/pbmc_analysis.Rmd`)
This R Markdown file demonstrates an end-to-end Seurat pipeline:
```r
library(Seurat)
library(dplyr)
# Load the data from MTX format
pbmc.data <- Read10X(data.dir = '../data/filtered_gene_bc_matrices/hg19/')
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
# Calculate mitochondrial content
pbmc[['percent.mt']] <- PercentageFeatureSet(pbmc, pattern = '^MT-')
# QC filtering
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# Normalize & identify variable features
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = 'vst', nfeatures = 2000)
# Scale data and run PCA
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Cluster and visualize
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = 'umap', label = TRUE)
```

## Directory Structure
```
scRNA_tutorial/
├── data/                     # Example datasets (PBMC 3K)
│   └── README.md             # Instructions for downloading datasets
├── analysis/                 # Notebooks and Rmd for hands-on tutorials
│   ├── pbmc_analysis.ipynb   # Python Notebook
│   └── pbmc_analysis.Rmd     # R Markdown
```

Happy analyzing! 🚀
