# Visualization & Reporting

Effective visualization and interpretation are essential for drawing meaningful biological conclusions from single-cell RNA-seq data. This section covers:

1. Publication-quality plotting
2. Guidelines for interpreting key plot types
3. Interactive dashboards
4. Reproducible reporting

---

## 1. Publication-Quality Plots

### R (Seurat + ggplot2)
```r
# UMAP plot with labels and custom theme
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP of PBMC Clusters") +
  theme_minimal()

# Violin plot of marker gene expression
VlnPlot(pbmc, features = "CD3D", pt.size = 0) +
  labs(title = "CD3D Expression Across Clusters") +
  theme_classic()

# Dot plot showing expression and percent of cells
DotPlot(pbmc, features = c("CD3D", "MS4A1")) +
  scale_color_viridis_c() +
  labs(title = "Marker Expression per Cluster")
```

### Python (Scanpy + Matplotlib/Seaborn)
```python
# UMAP plot
sc.pl.umap(adata, color='leiden', title='UMAP of PBMC Clusters', size=20)

# Violin plot
sc.pl.violin(adata, keys='CD3D', groupby='leiden', stripplot=False, title='CD3D Expression Across Clusters')

# Dot plot
sc.pl.dotplot(adata, var_names=['CD3D','MS4A1'], groupby='leiden', title='Marker Expression per Cluster')
```

---

## 2. Interpreting Results

### Dimensionality Reduction Plots (UMAP/t-SNE)
- How to read:
  - **Axes**: UMAP1 and UMAP2 are abstract dimensions; only relative distances between points convey biological meaning.
  - **Clusters**: Colors or shapes identify clusters (cell types/states); distinct groupings imply discrete identities.
  - **Distance**: Proximity indicates similarity in gene expression profiles.
  - **Outliers**: Isolated points may represent rare cell types, doublets, or technical noise.
  - **Density**: High point density highlights abundant or common cell states.
  - **Legend**: Use the legend to match colors to cluster labels for accurate interpretation.
- Interpretation:
  - Well-separated clusters suggest distinct cellular identities.
  - Overlapping clusters may indicate transitional or closely related cell states.

### Violin & Box Plots
- Show distribution of gene expression within each cluster.
- **Median and quartiles** reveal heterogeneity within clusters.
- Wide violins indicate high variability; narrow suggest uniform expression.

### Dot & Heatmaps
- **Dot size** indicates the fraction of cells expressing the gene; color indicates average expression.
- Heatmaps of top markers help identify cluster-specific signature patterns.
- Use hierarchical clustering on heatmaps to reveal substructure.

### Pseudotime & Trajectory Plots
- Cells ordered by pseudotime reflect progression through a process.
- Branch points mark fate decisions; examine gene changes at these points.
- Overlay gene expression on trajectory to see dynamic regulation.

### Cell-Cell Communication Networks
- **Nodes** represent cell types/clusters; **edges** represent inferred ligand-receptor interactions.
- Edge thickness or color often encodes interaction strength or significance.
- Focus on top interacting pairs for biological validation.

---

## 3. Interactive Dashboards

- **Shiny (R)**: Launch with `runApp('shiny_app/')` to explore Seurat objects interactively.
- **Dash/Streamlit (Python)**: Use `dash` or `streamlit` to create lightweight web apps for custom visualization.
- **cellxgene**: `cellxgene launch data.h5ad` provides a ready-to-use explorer.
- **Loupe Cell Browser**: Import Cell Ranger output for annotation and interactive plots.

---

## 4. Reproducible Reporting

- Combine code and narrative using **R Markdown** (`.Rmd`) or **Jupyter Notebooks**.
- Export figures as **SVG** or **PDF** for high-resolution publication.
- Automate workflows with **Snakemake** or **Nextflow**, and version-control with **Git**.

---

## References

- Satija R, et al. (2015). Spatial reconstruction of single-cell gene expression data. _Nat Biotechnol_, 33(5):495–502.
- Wickham H. (2016). _ggplot2: Elegant Graphics for Data Analysis_. Springer.
- Bostock M, Ogievetsky V, Heer J. (2011). D3: Data-Driven Documents. _IEEE Trans Vis Comput Graph_, 17(12):2301–2309.
- Korotkevich G, et al. (2021). Fast, scalable, and accurate differential expression analysis for single-cell RNA-seq data. _Nat Methods_.
