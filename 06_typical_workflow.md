# Typical Workflow & Analysis Pipeline

A typical single-cell RNA-seq analysis pipeline combines the experimental setup and computational workflow. Below is a step-by-step narrative of each phase, with example code snippets in R (Seurat) and Python (Scanpy) to illustrate the core analysis steps.

## 1. Experimental Design & Library Preparation

Plan cell isolation, barcoding, and sequencing depth according to your biological question. See [Library Preparation](./library_preparation.md) for protocol comparisons and QC metrics.

## 2. Sequencing

Generate raw sequencing reads (FASTQ files) using Illumina platforms (e.g., NovaSeq, NextSeq). These files contain the raw read sequences and quality scores.

## 3. Preprocessing & Quantification

Use Cell Ranger or STARsolo to demultiplex reads, align to the reference genome, and count UMIs per gene:

Python (Scanpy):
```python
import scanpy as sc
adata = sc.read_10x_h5('filtered_feature_bc_matrix.h5')
```

R (Seurat):
```r
library(Seurat)
pbmc.data <- Read10X_h5('filtered_feature_bc_matrix.h5')
pbmc <- CreateSeuratObject(pbmc.data)
```

## 4. Quality Control (QC)

Filter out low-quality cells and genes based on metrics (e.g., gene count, mitochondrial content):

Python:
```python
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs.pct_counts_mt < 5, :]
```

R:
```r
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 5)
```

## 5. Normalization & Scaling

Normalize counts to correct for sequencing depth and log-transform the data:

Python:
```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

R:
```r
pbmc <- NormalizeData(pbmc)
pbmc <- ScaleData(pbmc)
```

## 6. Batch Correction & Integration

Batch effects can obscure true biological signals when combining datasets from different experiments. Popular integration methods include:

- **Seurat Integration (Anchors in R)**
  ```r
  # Split by batch and preprocess
  pbmc.list <- SplitObject(pbmc, split.by = "batch")
  pbmc.list <- lapply(pbmc.list, NormalizeData)
  pbmc.list <- lapply(pbmc.list, FindVariableFeatures)
  
  # Find integration anchors and integrate
  anchors <- FindIntegrationAnchors(object.list = pbmc.list)
  integrated <- IntegrateData(anchors)
  # Use 'integrated' for downstream analysis
  ```
- **Harmony (R)**
  ```r
  library(harmony)
  pbmc <- RunHarmony(pbmc, group.by.vars = "batch")
  ```
- **BBKNN (Python)**
  ```python
  import bbknn
  sc.pp.pca(adata, n_comps=50)
  bbknn.bbknn(adata, batch_key='batch')
  sc.tl.umap(adata)
  sc.pl.umap(adata, color='batch')
  ```
- **Mutual Nearest Neighbors (MNN in R)**
  ```r
  library(batchelor)
  sce <- SingleCellExperiment(assays = list(counts = counts(pbmc)))
  corrected <- fastMNN(sce, batch = sce$batch)
  ```
- **Scanorama (Python)**
  ```python
  import scanorama
  adatas = [adata1, adata2]
  integrated, genes = scanorama.correct(adatas, return_dimred=False)
  ```

Choose an integration method based on dataset characteristics and computational resources.

## 7. Dimensionality Reduction

Identify major axes of variation and visualize in 2D:

Python:
```python
sc.pp.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_pcs=10)
sc.tl.umap(adata)
sc.pl.umap(adata)
```

R:
```r
pbmc <- RunPCA(pbmc, npcs = 50)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
```

## 8. Clustering & Cell Type Identification

Cluster cells and visualize cluster assignments:

Python:
```python
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color='leiden')
```

R:
```r
pbmc <- FindClusters(pbmc, resolution = 0.5)
DimPlot(pbmc, group.by = 'seurat_clusters')
```

## 9. Differential Expression Analysis

Identify marker genes for clusters:

Python:
```python
sc.tl.rank_genes_groups(adata, groupby='leiden', method='t-test')
sc.pl.rank_genes_groups(adata)
```

R:
```r
markers <- FindMarkers(pbmc, ident.1 = 0, ident.2 = 1)
head(markers)
```

## 10. Visualization & Interpretation

Generate publication-ready plots (feature plots, violin plots) and interpret results in the context of your biological question. For interactive exploration, consider tools like [Loupe Cell Browser](https://support.10xgenomics.com/) or [Cellxgene](https://chanzuckerberg.github.io/cellxgene/).

## References

- Stuart T, et al. (2019). Comprehensive Integration of Single-Cell Data. Cell, 177(7), 1888–1902.e21. https://www.cell.com/fulltext/S0092-8674(19)30559-8
- Wolf FA, et al. (2018). Scanpy: large-scale single-cell gene expression data analysis. Genome Biol, 19, 15. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1382-0
- Luecken MD & Theis FJ (2019). Current best practices in single-cell RNA-seq analysis: a tutorial. Mol Syst Biol, 15(6):e8746. https://www.embopress.org/doi/full/10.15252/msb.20188746
- Orchestrating Single-Cell Analysis with Bioconductor: https://osca.bioconductor.org/
- Korsunsky I, et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods, 16(12), 1289–1296. https://www.nature.com/articles/s41592-019-0619-0
- Haghverdi L, et al. (2018). Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. Nat Biotechnol, 36, 421–427. https://www.nature.com/articles/nbt.4091
- Polański K, et al. (2020). BBKNN: fast batch alignment of single cell transcriptomes. Bioinformatics, 36(5), 964–965. https://academic.oup.com/bioinformatics/article/36/5/964/5574211
- Hie B, et al. (2019). Efficient integration of heterogeneous single-cell transcriptomes using Scanorama. Nat Biotechnol, 37, 685–691. https://www.nature.com/articles/s41587-019-0113-5
