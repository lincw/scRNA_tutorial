# Hands-On Examples

Follow these examples to learn by doing with the PBMC 3K dataset. Notebooks and scripts are located in the `analysis/` directory.

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
adata = sc.read_10x_h5('data/pbmc_3k_filtered_feature_bc_matrix.h5')
# QC
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs.pct_counts_mt < 5, :]
# Normalize & log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# PCA & UMAP
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color='leiden')
```

## R Markdown (`analysis/pbmc_analysis.Rmd`)
This R Markdown file demonstrates an end-to-end Seurat pipeline:
```r
library(Seurat)
pbmc.data <- Read10X_h5('data/pbmc_3k_filtered_feature_bc_matrix.h5')
pbmc <- CreateSeuratObject(pbmc.data)
# QC
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 5)
# Normalize & scale
pbmc <- NormalizeData(pbmc)
pbmc <- ScaleData(pbmc)
# PCA & UMAP
pbmc <- RunPCA(pbmc, npcs = 50)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = 'umap')
```

## Directory Structure
```
scRNAseq_guide/
├── data/                     # Example datasets (PBMC 3K)
├── analysis/                 # Notebooks and Rmd for hands-on tutorials
│   ├── pbmc_analysis.ipynb   # Python Notebook
│   └── pbmc_analysis.Rmd     # R Markdown
```

Happy analyzing! 🚀
