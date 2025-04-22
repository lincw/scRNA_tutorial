# Key Concepts

Understanding these foundational concepts is crucial for interpreting single-cell RNA-seq data:

## 1. Transcriptome
The complete set of RNA transcripts expressed in a cell, including mRNA, non-coding RNA, and other RNA species.

## 2. Cell Barcode & UMI
- **Cell Barcode:** A short, unique sequence attached to each cell’s cDNA during library preparation, enabling assignment of reads to individual cells.
- **Unique Molecular Identifier (UMI):** A random sequence tag added to each mRNA molecule before amplification to distinguish original transcripts from PCR duplicates.

## 3. Gene–Cell Count Matrix
A matrix where rows represent genes and columns represent cells; each entry is the number of UMIs (or reads) detected for that gene in that cell:

```
Gene \ Cell | Cell1 | Cell2 | Cell3
------------|-------|-------|------
GeneA       | 5     | 0     | 3    
GeneB       | 0     | 2     | 1    
GeneC       | 10    | 8     | 0    
```

## 4. Dropout (Zero Inflation)
Technical or biological events where a transcript present in a cell is not captured or sequenced, resulting in excess zeros. Dropouts complicate downstream analysis and require specialized handling.

## 5. Batch Effects
Systematic differences in data introduced by variations in experimental conditions, reagents, or processing dates. Batch correction methods (e.g., Seurat integration, Harmony, MNN) are used to mitigate these artifacts.

## 6. Quality Metrics
Key per-cell metrics include:
- **Library Size:** Total UMIs or reads per cell.
- **Gene Count:** Number of genes detected per cell.
- **Mitochondrial Content:** Fraction of reads mapping to mitochondrial genes, used to identify low-quality or dying cells.

## 7. Clustering
Grouping cells with similar expression profiles into clusters that often correspond to cell types or states. Common algorithms include Louvain, Leiden, and hierarchical clustering.

## 8. Dimensionality Reduction
Techniques to project high-dimensional gene expression data into 2D/3D for visualization and noise reduction:
- **PCA:** Principal Component Analysis
- **t-SNE:** t-distributed Stochastic Neighbor Embedding
- **UMAP:** Uniform Manifold Approximation and Projection

