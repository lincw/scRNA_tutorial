# Typical Workflow & Analysis Pipeline

A typical single-cell RNA-seq analysis pipeline combines the experimental setup and computational workflow. Below is a step-by-step narrative of each phase, with example code snippets in R (Seurat) and Python (Scanpy) to illustrate the core analysis steps.

## 1. Experimental Design & Library Preparation

Plan cell isolation, barcoding, and sequencing depth according to your biological question. See [Library Preparation](./05_library_preparation.md) for protocol comparisons and QC metrics.

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

Filter out low-quality cells and genes based on metrics, including standard filtering and doublet detection:

### Standard QC Filtering

Python:
```python
# Calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # mitochondrial genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Visualize QC metrics before filtering
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4)

# Apply QC thresholds
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=5000)  # Remove potential doublets based on gene count
adata = adata[adata.obs.pct_counts_mt < 10, :]  # Filter cells with high mitochondrial content

# Filter genes
sc.pp.filter_genes(adata, min_cells=3)  # Keep genes expressed in at least 3 cells
```

R:
```r
# Calculate mitochondrial percentage
pbmc[['percent.mt']] <- PercentageFeatureSet(pbmc, pattern = '^MT-')

# Visualize QC metrics
VlnPlot(pbmc, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)

# Apply QC thresholds
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
```

### Doublet Detection (Highly Recommended)

Python (using Scrublet):
```python
import scrublet as scr

# Initialize Scrublet
scrub = scr.Scrublet(adata.X)

# Run doublet detection
doublet_scores, predicted_doublets = scrub.scrub_doublets()

# Add results to adata
adata.obs['doublet_scores'] = doublet_scores
adata.obs['predicted_doublets'] = predicted_doublets

# Visualize results
scrub.plot_histogram()

# Filter out predicted doublets
adata = adata[~adata.obs.predicted_doublets, :]
```

R (using DoubletFinder):
```r
library(DoubletFinder)

# Run PCA first (required for DoubletFinder)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)

# Find optimal pK value
sweetp <- paramSweep_v3(pbmc)
sweetscore <- summarizeSweep(sweetp)
bcmvn <- find.pK(sweetscore)
pk_value <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

# Run DoubletFinder
homotypic_prop <- modelHomotypic(pbmc$seurat_clusters) 
nExp <- round(0.08 * ncol(pbmc))  # Expect 8% doublets (adjust based on 10X loading density)
pbmc <- doubletFinder_v3(pbmc, pN=0.25, pK=pk_value, nExp=nExp, reuse.pANN=FALSE)

# Filter out predicted doublets (column name depends on parameters used)
doublet_col <- colnames(pbmc@meta.data)[grepl("DF.classifications", colnames(pbmc@meta.data))]
pbmc <- subset(pbmc, subset = get(doublet_col) == "Singlet")
```

## 5. Normalization & Feature Selection

Normalize counts to correct for sequencing depth, log-transform the data, and identify highly variable genes for downstream analysis:

Python:
```python
# Basic normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Find highly variable genes (critical step!)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print(f"Number of highly variable genes: {sum(adata.var.highly_variable)}")

# Keep only highly variable genes for dimensionality reduction
adata = adata[:, adata.var.highly_variable]

# Scale data
sc.pp.scale(adata, max_value=10)
```

R:
```r
# Standard normalization
pbmc <- NormalizeData(pbmc)

# Find variable features (standard workflow)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale data (using variable features by default)
pbmc <- ScaleData(pbmc)
```

### Alternative Normalization Methods

Modern workflows often use more sophisticated normalization approaches:

**SCTransform (Seurat):**
```r
# Alternative: SCTransform (performs normalization, variance stabilization, and feature selection)
pbmc <- SCTransform(pbmc, verbose = FALSE) 
```

**scran pooling-based normalization (Python):**
```python
# Alternative: scran-based normalization
sc.pp.normalize_per_cell(adata)  # Pre-normalize to equalize library size
sc.external.pp.scrna_normalize_scran(adata)  # Compute size factors and normalize
```

## 6. Batch Correction & Integration

Batch effects can obscure true biological signals when combining datasets from different experiments. Most integration methods should be applied **after feature selection** but **before dimensionality reduction**. The correct workflow order is critical for effective batch integration:

### Best Practices for Batch Integration
1. Perform QC and normalization on each batch separately
2. Identify variable features in each batch separately (or use a joint approach)
3. Apply integration methods to align batches in a shared feature space
4. Proceed with integrated data for dimensionality reduction and clustering

Popular integration methods include:

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

Cluster cells, visualize cluster assignments, and carefully annotate cell types using multiple validation approaches:

Python:
```python
# Try multiple clustering resolutions to find optimal granularity
resolutions = [0.2, 0.4, 0.6, 0.8, 1.0]
for res in resolutions:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_res{res}')
    
# Visualize clustering at resolution 0.5 (or another chosen value)
sc.tl.leiden(adata, resolution=0.5, key_added='leiden')
sc.pl.umap(adata, color='leiden')
```

R:
```r
# Test multiple resolutions
pbmc <- FindClusters(pbmc, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0))

# Choose a resolution for downstream analysis
DimPlot(pbmc, group.by = 'seurat_clusters')
```

### Cell Type Annotation Best Practices

1. **Use multiple markers for each cell type:** Don't rely on a single gene to define a cell population

   ```python
   # Example: Check expression of multiple markers per cell type
   sc.pl.umap(adata, color=['CD3D', 'CD3E', 'CD4', 'CD8A', 'MS4A1', 'CD79A', 'LYZ', 'FCGR3A'])
   ```

2. **Cross-validate with orthogonal sources:**
   - Published reference datasets (e.g., Azimuth, SingleR, scReference)
   - Flow cytometry data if available
   - Spatial data or proteomics if available

3. **Test automated annotation methods:**
   ```r
   # Example with SingleR
   library(SingleR)
   reference <- HumanPrimaryCellAtlasData()
   predictions <- SingleR(test=GetAssayData(pbmc), ref=reference, 
                         labels=reference$label.main)
   pbmc$SingleR.labels <- predictions$labels
   ```

4. **Validate with domain experts:** Collaborate with biologists who understand the expected cell types in your tissue

5. **Use consistent naming conventions:** Follow community standards like the Cell Ontology (https://www.ebi.ac.uk/ols/ontologies/cl)

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
