# Glossary & FAQ

## Glossary

- **UMI (Unique Molecular Identifier):** Short random sequences added to each molecule during library preparation to distinguish PCR duplicates from true biological molecules.
- **Cell Barcode:** Unique sequence tags that label reads from the same cell.
- **Count Matrix:** Matrix of gene (rows) by cell (columns) containing UMI counts for each gene in each cell.
- **Normalization:** Process of adjusting for sequencing depth or cell-specific biases.
- **Scaling:** Centering and variance-standardizing data across genes.
- **PCA (Principal Component Analysis):** Linear dimensionality reduction capturing major axes of variation.
- **UMAP (Uniform Manifold Approximation and Projection):** Non-linear dimensionality reduction preserving local and global structure.
- **Clustering:** Grouping cells with similar expression profiles into clusters (cell types or states).
- **Pseudotime:** Relative ordering of cells along inferred trajectories reflecting progression through a biological process.
- **Batch Effect:** Technical variation between samples or runs that can confound biological signals.

## Frequently Asked Questions (FAQ)

**Q: How do I choose clustering resolution?**
A: Test multiple resolutions (e.g., 0.4–1.2) and validate clusters by marker gene expression and known biology. Use silhouette scores or stability metrics when possible.

**Q: PCA vs. UMAP vs. t-SNE—when to use each?**
A: Use PCA for initial variance analysis and noise reduction. UMAP offers fast computation and preserves global structure. t-SNE excels at local neighborhood structure but can distort global relationships.

**Q: What are good QC thresholds?**
A: Common filters: remove cells with <200 genes, >5–10% mitochondrial reads, and remove genes detected in <3 cells. Adjust thresholds based on dataset characteristics.

**Q: How do I correct batch effects?**
A: Use methods like Seurat Integration, Harmony, BBKNN, or MNN. Always inspect before/after correction with diagnostic plots to ensure you remove unwanted variation while preserving biology.

**Q: How is pseudotime different from actual time?**
A: Pseudotime is a relative measure ordering cells along a trajectory. It reflects progression through a process but not real chronological time.

**Q: Can I combine datasets from different platforms?**
A: Yes, but apply batch correction methods and ensure consistent preprocessing steps to minimize technical biases.

**Q: What file formats should I use for sharing data?**
A: Use HDF5-based formats like `.h5ad` (Scanpy) or `.h5seurat` for sharing processed AnnData or Seurat objects, and CSV/MTX for raw count matrices.
